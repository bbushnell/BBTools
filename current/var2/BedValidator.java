package var2;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.IntList;

/**
 * Stateless checks for BED core fields and UCSC custom-track files.
 * Supports unsigned64 coordinates, BED3-9/BED12, and explicit BEDn+ prefixes.
 * Chromosome names may contain printable punctuation; reference bounds, sorting,
 * custom-field semantics, and strict GA4GH chromosome naming are not checked.
 * Input scanning supports LF/CRLF through ByteFile1.
 * @author Shinobu
 * @date Sep 26, 2026
 */
public final class BedValidator {

	private BedValidator(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks an inferred BED core using fresh scratch; reuse scratch for repeated calls. */
	public static String validateFields(final byte[] line){return validateFields(line, new Scratch(), 0);}

	/**
	 * Checks one data line without modifying it. Scratch must be exclusive to the caller.
	 * A valid single-tab interpretation takes precedence; otherwise horizontal
	 * whitespace is tried. Explicit standardFields disambiguates BEDn+ custom columns.
	 * @param line Data bytes without a line terminator
	 * @param scratch Reusable private parsing state
	 * @param standardFields Zero to infer, otherwise3-9 or12; later columns are custom
	 * @return Null on success, otherwise the first core-field diagnostic
	 */
	public static String validateFields(final byte[] line, final Scratch scratch, final int standardFields){
		if(scratch==null || (standardFields!=0 && !standardCount(standardFields))){
			throw new IllegalArgumentException("BED scratch is required; standard fields must be auto, 3-9, or 12.");
		}
		if(line==null || line.length==0){return "Missing BED data line";}
		scratch.set(line, standardFields);
		final String error=validateParsed(scratch, standardFields);
		if(error!=null && scratch.tabDelimited){
			scratch.setWhitespace(line);
			if(validateParsed(scratch, standardFields)==null){return null;}
			scratch.tabDelimited=true;//Restore the interpretation associated with the first diagnostic.
		}
		return error;
	}

	/** Checks already selected field boundaries; allocation is confined to invalid diagnostics. */
	private static String validateParsed(final Scratch scratch, final int standardFields){
		final byte[] line=scratch.line;
		assert(line!=null) : "validateFields must install a nonnull BED line before field checks";
		final int count=scratch.fieldCount(), core=standardFields==0 ? Math.min(count, 12) : standardFields;
		if(!standardCount(core)){return "Expected BED3-9 or BED12; use bedfields=N to declare a custom BEDn+ prefix";}
		if(count<core){return "Fewer columns than the declared BED prefix: "+count+" < "+core;}
		for(byte b : line){if(b!='\t' && (b<32 || b>126)){return "BED fields must contain printable ASCII characters";}}
		for(int i=0; i<core; i++){if(scratch.end(i)==scratch.start(i)){return "Empty standard BED field "+(i+1);}}
		if(scratch.end(0)-scratch.start(0)>255 || scratch.spaceIn(0)){return "Chromosome name must be 1-255 non-whitespace characters";}
		if(core>=4 && scratch.end(3)-scratch.start(3)>255){return "BED name exceeds 255 characters";}
		try{
			final long start=number(scratch, 1), end=number(scratch, 2);
			if(Long.compareUnsigned(start, end)>0){return "chromStart exceeds chromEnd";}
			if(core>=5 && Long.compareUnsigned(number(scratch, 4), 1000)>0){return "Score must be an integer from 0 through 1000";}
			if(core>=6 && !(scratch.equals(5, "+") || scratch.equals(5, "-") || scratch.equals(5, "."))){
				return "Strand must be '+', '-', or '.'";
			}
			if(core>=7){
				final long thickStart=number(scratch, 6);
				if(Long.compareUnsigned(thickStart, start)<0 || Long.compareUnsigned(thickStart, end)>0){return "thickStart is outside the feature";}
				if(core>=8){
					final long thickEnd=number(scratch, 7);
					if(Long.compareUnsigned(thickEnd, thickStart)<0 || Long.compareUnsigned(thickEnd, end)>0){return "thickEnd is outside [thickStart, chromEnd]";}
				}
			}
			if(core>=9){
				final String error=rgb(scratch);
				if(error!=null){return error;}
			}
			if(core==12){return blocks(scratch, end-start);}
		}catch(NumberFormatException e){return e.getMessage();}
		return null;
	}

	/** Scans a BED/custom-track file with local state; custom column meanings remain unchecked. */
	public static Result validateFile(final FileFormat ff, final int standardFields){
		if(ff==null || !ff.read() || (standardFields!=0 && !standardCount(standardFields))){
			throw new IllegalArgumentException("Expected BED input and standard fields auto, 3-9, or 12.");
		}
		final ByteFile bf=new ByteFile1(ff);
		final Scratch scratch=new Scratch();
		long lines=0, records=0, invalid=0, comments=0, blanks=0, tracks=0, custom=0, firstLine=0;
		int expectedColumns=-1;
		boolean readError, sawWhitespace=false, needsTabs=false;
		String firstError=null;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lines++;
				if(blank(line)){blanks++; continue;}
				if(line[0]=='#'){comments++; continue;}
				final int directive=trackLine(line);
				if(directive>0){
					tracks++;
					if(directive==1){expectedColumns=-1; sawWhitespace=false; needsTabs=false;}
					else if(expectedColumns>=0){
						invalid++;
						if(firstError==null){firstError="Browser directives must precede track data"; firstLine=lines;}
					}
					continue;
				}
				records++;
				String error=validateFields(line, scratch, standardFields);
				final int count=scratch.fieldCount();
				if(expectedColumns<0){expectedColumns=count;}
				if(error==null && count!=expectedColumns){error="Inconsistent field count within a track: "+count+" versus "+expectedColumns;}
				if(!scratch.tabDelimited){sawWhitespace=true;}
				if(scratch.tabDelimited && scratch.spacedNameOrCustom(standardFields)){needsTabs=true;}
				if(error==null && sawWhitespace && needsTabs){error="Spaces inside fields require single-tab delimiters throughout a track";}
				if(count>(standardFields==0 ? 12 : standardFields)){custom++;}
				if(error!=null){
					invalid++;
					if(firstError==null){firstError=error; firstLine=lines;}
				}
			}
		}finally{readError=bf.close();}
		assert(lines==records+comments+blanks+tracks) : "BED scan must classify every physical line exactly once";
		return new Result(lines, records, invalid, comments, blanks, tracks, custom, firstLine, firstError, readError);
	}

	/** Prints explicit core coverage and bounded diagnostics; returns whether the scan passed. */
	public static boolean printResult(final Result r, final FileFormat ff, final PrintStream out, final PrintStream err){
		out.println("File\t\t"+ff.name());
		out.println("Format\t\tbed");
		out.println("Compression\t"+(ff.bgzip() ? "bgzip" : FileFormat.COMPRESSION_ARRAY[ff.compression()]));
		out.println("ValidationScope\tBED core fields and block structure; custom fields, reference bounds, sorting, and track attributes unchecked");
		out.println("Intervals\t"+r.records);
		out.println("InvalidLines\t"+r.invalid);
		out.println("CommentLines\t"+r.comments);
		out.println("BlankLines\t"+r.blanks);
		out.println("TrackLines\t"+r.tracks);
		out.println("CustomFieldRecords\t"+r.customRecords);
		if(r.firstError!=null){err.println(ff.name()+":"+r.firstErrorLine+": "+r.firstError);}
		if(r.readError){err.println("Read/decompression/close failure in "+ff.name());}
		final boolean writeError=out.checkError();
		if(writeError){err.println("Failed to write BED validation report for "+ff.name());}
		return r.passed() && !writeError && !err.checkError();
	}

	/** Standard BED10/11 lack complete block columns and require an explicit custom prefix. */
	private static boolean standardCount(final int n){return (n>=3 && n<=9) || n==12;}

	/** Parses an unsigned64 decimal field, preserving all bits rather than signed-long range. */
	private static long number(final Scratch s, final int field){return unsigned(s.line, s.start(field), s.end(field), field+1);}

	/** Checked uint64 parsing; failures allocate diagnostics only on the invalid-input path. */
	private static long unsigned(final byte[] line, final int a, final int b, final int field){
		assert(a>=0 && a<=b && b<=line.length) : "BED numeric bounds must remain within the parsed line";
		long value=0;
		if(a==b){throw new NumberFormatException("Empty number in BED field "+field);}
		for(int i=a; i<b; i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || Long.compareUnsigned(value, 1844674407370955161L)>0 ||
				(value==1844674407370955161L && digit>5)){
				throw new NumberFormatException("BED field "+field+" requires an unsigned 64-bit decimal integer");
			}
			value=value*10+digit;
		}
		return value;
	}

	/** Checks RGB triples, including the special scalar zero. */
	private static String rgb(final Scratch s){
		if(s.equals(8, "0")){return null;}
		int a=s.start(8), channels=0;
		final int end=s.end(8);
		while(a<end){
			int b=a; while(b<end && s.line[b]!=','){b++;}
			if(Long.compareUnsigned(unsigned(s.line, a, b, 9), 255)>0){return "RGB components must be from 0 through 255";}
			channels++;
			if(b==end){return channels==3 ? null : "RGB requires exactly three components or scalar 0";}
			a=b+1;
		}
		return "RGB requires exactly three components without a trailing comma";
	}

	/** Checks paired block lists without materializing arrays or overflowing unsigned extents. */
	private static String blocks(final Scratch s, final long span){
		final long n=number(s, 9);
		if(n==0 || Long.compareUnsigned(n, span)>0 || Long.compareUnsigned(n, s.line.length)>0){return "Block count must be positive and fit the feature and block lists";}
		int sizeA=s.start(10), startA=s.start(11);
		final int sizeEnd=s.end(10), startEnd=s.end(11);
		long previousEnd=0;
		for(int i=0; i<(int)n; i++){
			if(sizeA>=sizeEnd || startA>=startEnd){return "Block lists contain fewer entries than blockCount";}
			int sizeB=sizeA, startB=startA;
			while(sizeB<sizeEnd && s.line[sizeB]!=','){sizeB++;}
			while(startB<startEnd && s.line[startB]!=','){startB++;}
			final long size=unsigned(s.line, sizeA, sizeB, 11), start=unsigned(s.line, startA, startB, 12);
			if(i==0 && start!=0){return "The first block must begin at chromStart";}
			if(Long.compareUnsigned(start, previousEnd)<0){return "Blocks overlap or are out of order";}
			if(Long.compareUnsigned(start, span)>0 || Long.compareUnsigned(size, span-start)>0){return "Block extends beyond chromEnd";}
			previousEnd=start+size;
			sizeA=sizeB<sizeEnd ? sizeB+1 : sizeB;
			startA=startB<startEnd ? startB+1 : startB;
		}
		if(sizeA!=sizeEnd || startA!=startEnd){return "Block lists contain more entries than blockCount";}
		return previousEnd==span ? null : "The last block must end at chromEnd";
	}

	/** Blank BED lines may contain horizontal whitespace. */
	private static boolean blank(final byte[] line){for(byte b : line){if(b!=' ' && b!='\t'){return false;}} return true;}

	/** Returns1 for track,2 for browser,0 for data; no per-line token allocation. */
	private static int trackLine(final byte[] line){
		return directive(line, "track") ? 1 : directive(line, "browser") ? 2 : 0;
	}

	/** Recognizes a directive word while preserving numeric BED rows on identically named chromosomes. */
	private static boolean directive(final byte[] line, final String word){
		if(line.length<word.length()){return false;}
		int i=0; while(i<word.length() && line[i]==word.charAt(i)){i++;}
		if(i<word.length()){return false;}
		if(i==line.length){return false;}
		if(line[i]!=' ' && line[i]!='\t'){return false;}
		while(i<line.length && (line[i]==' ' || line[i]=='\t')){i++;}
		if(i<line.length && line[i]>='0' && line[i]<='9'){return false;}
		if(word.equals("browser")){
			for(String command : BROWSER_COMMANDS){
				int j=0;
				while(j<command.length() && i+j<line.length && line[i+j]==command.charAt(j)){j++;}
				if(j==command.length() && (i+j==line.length || line[i+j]==' ' || line[i+j]=='\t')){return true;}
			}
			return false;
		}
		for(; i<line.length; i++){if(line[i]=='='){return true;}}
		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Scratch State        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reusable, caller-owned parsing buffers; never share concurrently between calls. */
	public static final class Scratch {
		/** Number of fields in the most recently checked data line. */
		public int fieldCount(){return tabDelimited ? tabs.terms() : starts.size;}
		/** Selects unambiguous tab fields when possible, otherwise scans horizontal whitespace. */
		private void set(final byte[] value, final int core){
			line=value; tabs.set(line);
			tabDelimited=tabs.terms()>=Math.max(3, core);
			if(tabDelimited){
				for(int field=0; field<3; field++){
					tabs.setBounds(field);
					if(tabs.a()==tabs.b()){tabDelimited=false; break;}
					for(int i=tabs.a(); i<tabs.b(); i++){
						if(line[i]==' ' || (field>0 && (line[i]<'0' || line[i]>'9'))){tabDelimited=false; break;}
					}
				}
			}
			if(tabDelimited){return;}
			setWhitespace(line);
		}
		/** Selects reusable whitespace-token bounds, including repeated separators. */
		private void setWhitespace(final byte[] value){
			line=value; tabDelimited=false;
			starts.clear(); ends.clear();
			for(int i=0; i<line.length;){
				while(i<line.length && (line[i]==' ' || line[i]=='\t')){i++;}
				if(i==line.length){break;}
				starts.add(i);
				while(i<line.length && line[i]!=' ' && line[i]!='\t'){i++;}
				ends.add(i);
			}
		}
		/** Gets a field boundary from the selected parser. */
		private int start(final int field){if(tabDelimited){tabs.setBounds(field); return tabs.a();} return starts.get(field);}
		/** Gets the exclusive end of a parsed field. */
		private int end(final int field){if(tabDelimited){tabs.setBounds(field); return tabs.b();} return ends.get(field);}
		/** Compares a field without allocating a String. */
		private boolean equals(final int field, final String value){
			final int a=start(field), b=end(field);
			if(b-a!=value.length()){return false;}
			for(int i=a; i<b; i++){if(line[i]!=value.charAt(i-a)){return false;}}
			return true;
		}
		/** Checks whether a tab-delimited field contains spaces. */
		private boolean spaceIn(final int field){for(int i=start(field), b=end(field); i<b; i++){if(line[i]==' '){return true;}} return false;}
		/** Detects fields whose embedded spaces require a tab-only track. */
		private boolean spacedNameOrCustom(final int core){
			if(fieldCount()>=4 && spaceIn(3)){return true;}
			for(int i=core==0 ? 12 : core; i<fieldCount(); i++){
				if(start(i)==end(i) || spaceIn(i)){return true;}
			}
			return false;
		}
		private byte[] line;
		private boolean tabDelimited;
		private final LineParser1 tabs=new LineParser1('\t');
		private final IntList starts=new IntList(), ends=new IntList();
	}

	/** Immutable per-file counts and first diagnostic; no global validation state. */
	public static final class Result {
		private Result(final long l, final long r, final long i, final long c, final long b, final long t,
			final long custom, final long first, final String error, final boolean io){
			lines=l; records=r; invalid=i; comments=c; blanks=b; tracks=t; customRecords=custom;
			firstErrorLine=first; firstError=error; readError=io;
		}
		/** True if all documented checks and input reading succeeded. */
		public boolean passed(){return invalid==0 && !readError;}
		/** Scan counters; customRecords counts rows with unchecked extension columns. */
		public final long lines, records, invalid, comments, blanks, tracks, customRecords;
		/** One-based first failing line, or zero when no line failed. */
		public final long firstErrorLine;
		/** First bounded diagnostic, or null. */
		public final String firstError;
		/** Reader/decompression/close failure flag. */
		public final boolean readError;
	}

	/** Documented UCSC browser commands; attributes/values themselves remain unchecked. */
	private static final String[] BROWSER_COMMANDS={"position", "hide", "dense", "pack", "squish", "full"};
}
