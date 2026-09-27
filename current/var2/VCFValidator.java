package var2;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Tools;
import structures.IntList;

/**
 * Reusable VCF4.1-4.5 header and record-structure checks.
 * Checks field counts, POS syntax, REF bases, QUAL syntax and FORMAT/sample shape.
 * Detailed ALT/GT syntax, INFO/FILTER/ID semantics, metadata declarations, reference
 * agreement, ordering, numeric ranges and final-newline conformance are unchecked.
 * File scanning tolerates blank lines after the required first-line version.
 * @author Shinobu
 * @date Sep 26, 2026
 */
public final class VCFValidator {

	private VCFValidator(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks a record against a header's column count, allocating scratch for this call. */
	public static String validateFields(final byte[] line, final int expectedColumns){
		return validateFields(line, new Scratch(), expectedColumns);
	}

	/**
	 * Checks a data record without changing its bytes or creating field Strings.
	 * POS is checked as nonnegative decimal text, including telomere zero, without
	 * imposing a machine integer range. ALT tokens are nonempty and whitespace-free;
	 * their detailed symbolic/breakend grammar is deliberately not interpreted.
	 * @param line Record bytes without CR/LF
	 * @param scratch Caller-owned scratch, not shared concurrently
	 * @param expectedColumns Column count established by a validated header
	 * @return Null on success, otherwise the first structural diagnostic
	 */
	public static String validateFields(final byte[] line, final Scratch scratch, final int expectedColumns){
		if(scratch==null || expectedColumns<8){throw new IllegalArgumentException("VCF validation needs scratch and a header with at least eight columns.");}
		if(line==null || line.length==0){return "Missing VCF record";}
		final LineParser1 p=scratch.fields;
		p.set(line);
		if(p.terms()!=expectedColumns){return "Record has "+p.terms()+" columns; header declares "+expectedColumns;}
		for(int i=0; i<p.terms(); i++){if(p.length(i)==0){return "Empty VCF field "+(i+1)+"; represent missing values with '.'";}}
		if(controlByte(line)){return "Unescaped control character in VCF record";}
		p.setBounds(0);
		if(p.termEquals('.', 0) || whitespace(line, p.a(), p.b())){return "CHROM must be a nonmissing identifier without whitespace";}
		p.setBounds(1);
		for(int i=p.a(), b=p.b(); i<b; i++){if(line[i]<'0' || line[i]>'9'){return "POS must be a nonnegative decimal integer (zero is allowed for telomeres)";}}
		p.setBounds(3);
		for(int i=p.a(), b=p.b(); i<b; i++){
			final int c=line[i]|32;
			if(c!='a' && c!='c' && c!='g' && c!='t' && c!='n'){return "REF must contain only A, C, G, T, or N bases";}
		}
		p.setBounds(4);
		if(whitespace(line, p.a(), p.b()) || line[p.a()]==',' || line[p.b()-1]==','){return "ALT contains whitespace or an empty alternative";}
		for(int i=p.a()+1, b=p.b(); i<b; i++){if(line[i]==',' && line[i-1]==','){return "ALT contains an empty alternative";}}
		p.setBounds(5);
		if(!p.currentTermEquals((byte)'.') && !floatSyntax(line, p.a(), p.b())){return "QUAL must be '.', a decimal/exponent number, infinity, or NaN";}
		return p.terms()>8 ? sampleFields(scratch) : null;
	}

	/**
	 * Checks the exact fixed header names, optional FORMAT, and unique nonempty
	 * sample IDs. Sample Strings are allocated only while checking the header.
	 * @param line The #CHROM line without a terminator
	 * @param scratch Caller-owned parser; fieldCount is available after this call
	 * @return Null on success, otherwise the first header diagnostic
	 */
	public static String validateColumnHeader(final byte[] line, final Scratch scratch){
		if(scratch==null){throw new IllegalArgumentException("VCF header scratch is required.");}
		if(line==null || line.length==0){return "Missing #CHROM header";}
		final LineParser1 p=scratch.fields;
		p.set(line);
		if(p.terms()<8){return "VCF column header requires at least eight columns";}
		if(controlByte(line)){return "Unescaped control character in VCF header";}
		for(int i=0; i<HEADER_FIELDS.length; i++){
			if(!p.termEquals(HEADER_FIELDS[i], i)){return "Expected "+HEADER_FIELDS[i]+" in header column "+(i+1);}
		}
		if(p.terms()>8 && !p.termEquals("FORMAT", 8)){return "Ninth header column must be FORMAT";}
		final HashSet<String> names=new HashSet<String>();
		for(int i=9; i<p.terms(); i++){
			p.setBounds(i);
			if(p.a()==p.b()){return "Empty sample name in column "+(i+1);}
			final String name=new String(line, p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
			if(!names.add(name)){return "Duplicate sample name in column "+(i+1);}
		}
		return null;
	}

	/** Scans to EOF with independent header state; preserves BBTools summary metadata. */
	public static Result validateFile(final FileFormat ff){
		if(ff==null || !ff.read()){throw new IllegalArgumentException("VCF validation requires an input FileFormat.");}
		final ByteFile bf=new ByteFile1(ff);
		final Scratch scratch=new Scratch();
		final Scan scan=new Scan();
		boolean versionSeen=false, headerSeen=false, dataSeen=false;
		int columns=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				scan.lines++;
				if(scan.lines==1 && !Tools.startsWith(line, FILEFORMAT)){scan.fail("First line must declare ##fileformat=VCFv4.x", scan.lines);}
				if(line.length==0){scan.blanks++; continue;}
				if(line[0]=='#'){
					scan.headers++;
					if(controlByte(line)){scan.fail("Unescaped control character in VCF header", scan.lines);}
					if(Tools.startsWith(line, FILEFORMAT)){
						if(versionSeen || scan.lines!=1){scan.fail("Duplicate or misplaced VCF version declaration", scan.lines);}
						versionSeen=true;
						if(line.length==VERSION_PREFIX.length()+1 && Tools.startsWith(line, VERSION_PREFIX) &&
							line[line.length-1]>='1' && line[line.length-1]<='5'){
							scan.version=new String(line, FILEFORMAT.length(), line.length-FILEFORMAT.length(), StandardCharsets.US_ASCII);
						}else{scan.fail("Unsupported VCF version; this structural checker covers VCF4.1 through VCF4.5", scan.lines);}
					}else if(Tools.startsWith(line, "##")){
						if(headerSeen || dataSeen){scan.fail("Metadata must precede the column header and data", scan.lines);}
						scan.metadata(line);
					}else{
						if(headerSeen || dataSeen){scan.fail("Duplicate or misplaced VCF column header", scan.lines);}
						final String error=validateColumnHeader(line, scratch);
						if(error!=null){scan.fail(error, scan.lines);}
						else if(!headerSeen && !dataSeen){headerSeen=true; columns=scratch.fieldCount();}
					}
					continue;
				}
				scan.records++; dataSeen=true;
				if(!headerSeen){scan.fail("Data record precedes a valid #CHROM header", scan.lines);}
				else{
					final String error=validateFields(line, scratch, columns);
					if(error!=null){scan.fail(error, scan.lines);}
				}
			}
		}finally{scan.readError=bf.close();}
		if(!versionSeen){scan.fail("Missing VCF version declaration at EOF", scan.lines+1);}
		if(!headerSeen){scan.fail("Missing valid #CHROM header at EOF", scan.lines+1);}
		assert(scan.lines==scan.records+scan.headers+scan.blanks) : "VCF scan must classify every physical line exactly once";
		return new Result(scan, columns);
	}

	/** Prints per-file structural coverage, existing metadata summaries, and bounded diagnostics. */
	public static boolean printResult(final Result r, final FileFormat ff, final PrintStream out, final PrintStream err){
		out.println("File\t\t"+ff.name());
		out.println("Format\t\tvcf");
		out.println("Compression\t"+(ff.bgzip() ? "bgzip" : FileFormat.COMPRESSION_ARRAY[ff.compression()]));
		out.println("HeaderLines\t"+r.headers);
		out.println("Variants\t"+r.records);
		out.println("ValidationScope\tVCF header/record structure, POS/REF/QUAL syntax, FORMAT/sample shape; ALT/GT details, INFO/FILTER/ID semantics, declarations, numeric ranges, reference, and sorting unchecked");
		out.println("Version\t\t"+r.version);
		out.println("Samples\t\t"+Math.max(0, r.columns-9));
		out.println("ValidationErrors\t"+r.errors);
		out.println("BlankLinesSkipped\t"+r.blanks);
		if(r.ploidy>0){out.println("Ploidy\t\t"+r.ploidy);}
		if(r.pairingRate>0){out.println("PairingRate\t"+Tools.format("%.4f", r.pairingRate));}
		if(r.mapqAvg>0){out.println("MapqAvg\t\t"+Tools.format("%.2f", r.mapqAvg));}
		if(r.qualityAvg>0){out.println("QualityAvg\t"+Tools.format("%.2f", r.qualityAvg));}
		if(r.readLengthAvg>0){out.println("ReadLengthAvg\t"+Tools.format("%.2f", r.readLengthAvg));}
		if(r.firstError!=null){err.println(ff.name()+":"+r.firstErrorLine+": "+r.firstError);}
		if(r.readError){err.println("Read/decompression/close failure in "+ff.name());}
		final boolean writeError=out.checkError();
		if(writeError){err.println("Failed to write VCF validation report for "+ff.name());}
		return r.passed() && !writeError && !err.checkError();
	}

	/** Checks key syntax/uniqueness and sample shape; trailing omitted subfields are valid. */
	private static String sampleFields(final Scratch s){
		final LineParser1 p=s.fields;
		final byte[] line=p.line();
		if(p.termEquals('.', 8)){
			for(int i=9; i<p.terms(); i++){if(!p.termEquals('.', i)){return "Samples must be missing when FORMAT is missing";}}
			return null;
		}
		s.starts.clear(); s.ends.clear();
		p.setBounds(8);
		final int end=p.b();
		for(int a=p.a(); a<=end;){
			int b=a; while(b<end && line[b]!=':'){b++;}
			if(a==b){return "Empty FORMAT key";}
			if(!alpha(line[a]) && line[a]!='_'){return "FORMAT keys must start with an ASCII letter or underscore";}
			for(int i=a+1; i<b; i++){
				final byte c=line[i];
				if(!alpha(c) && (c<'0' || c>'9') && c!='_' && c!='.'){return "Invalid character in FORMAT key";}
			}
			if(b-a==2 && line[a]=='G' && line[a+1]=='T' && s.starts.size>0){return "GT must be the first FORMAT key";}
			for(int i=0; i<s.starts.size; i++){
				if(equal(line, a, b, s.starts.get(i), s.ends.get(i))){return "Duplicate FORMAT key";}
			}
			s.starts.add(a); s.ends.add(b);
			if(b==end){break;}
			a=b+1;
		}
		for(int field=9; field<p.terms(); field++){
			p.setBounds(field);
			int count=1;
			for(int i=p.a(); i<p.b(); i++){
				if(line[i]==':'){
					if(i==p.a() || i==p.b()-1 || line[i-1]==':'){return "Empty sample subfield; use '.' or omit trailing values";}
					count++;
				}
			}
			if(count>s.starts.size){return "Sample column "+(field+1)+" has more values than FORMAT keys";}
		}
		return null;
	}

	/** Tests VCF Float lexical syntax, including case-insensitive signed infinity/NaN. */
	private static boolean floatSyntax(final byte[] line, int a, final int b){
		assert(a>=0 && a<b && b<=line.length) : "QUAL syntax requires nonempty parsed field bounds";
		if(line[a]=='+' || line[a]=='-'){a++;}
		if(equalIgnoreCase(line, a, b, "inf") || equalIgnoreCase(line, a, b, "infinity") || equalIgnoreCase(line, a, b, "nan")){return true;}
		int digits=0;
		while(a<b && line[a]>='0' && line[a]<='9'){a++; digits++;}
		if(a<b && line[a]=='.'){
			a++;
			final int first=a;
			while(a<b && line[a]>='0' && line[a]<='9'){a++; digits++;}
			if(a==first){return false;}
		}
		if(digits==0){return false;}
		if(a<b && (line[a]=='e' || line[a]=='E')){
			a++;
			if(a<b && (line[a]=='+' || line[a]=='-')){a++;}
			final int first=a;
			while(a<b && line[a]>='0' && line[a]<='9'){a++;}
			if(a==first){return false;}
		}
		return a==b;
	}

	/** Rejects raw control bytes while permitting UTF-8 non-ASCII bytes and tab separators. */
	private static boolean controlByte(final byte[] line){for(byte b : line){if((b>=0 && b<32 && b!='\t') || b==127){return true;}} return false;}
	/** Tests ASCII whitespace inside a field; tabs have already been separated. */
	private static boolean whitespace(final byte[] line, final int a, final int b){for(int i=a; i<b; i++){if(line[i]>=0 && line[i]<=32){return true;}} return false;}
	/** ASCII alphabetic test for FORMAT identifiers. */
	private static boolean alpha(final int c){return (c>='A' && c<='Z') || (c>='a' && c<='z');}
	/** Compares two ranges in one record without allocating keys. */
	private static boolean equal(final byte[] line, final int a, final int b, final int c, final int d){
		if(b-a!=d-c){return false;}
		for(int i=0; i<b-a; i++){if(line[a+i]!=line[c+i]){return false;}}
		return true;
	}
	/** Case-insensitive ASCII comparison for numeric literals and BBTools metadata keys. */
	private static boolean equalIgnoreCase(final byte[] line, final int a, final int b, final String value){
		if(b-a!=value.length()){return false;}
		for(int i=a; i<b; i++){if((line[i]|32)!=(value.charAt(i-a)|32)){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Local State          ----------------*/
	/*--------------------------------------------------------------*/

	/** Caller-owned buffers; use one per concurrent validation invocation. */
	public static final class Scratch {
		/** Header or record column count established by the most recent call. */
		public int fieldCount(){return fields.terms();}
		private final LineParser1 fields=new LineParser1('\t');
		private final IntList starts=new IntList(), ends=new IntList();
	}

	/** Private per-file counters and retained summary metadata. */
	private static final class Scan {
		/** Counts each failing physical line/EOF position once and retains one bounded diagnostic. */
		private void fail(final String error, final long line){
			if(lastErrorLine!=line){errors++; lastErrorLine=line;}
			if(firstError==null){firstError=error; firstErrorLine=line;}
		}
		/** Parses only the legacy BBTools summary keys; structured VCF metadata stays opaque. */
		private void metadata(final byte[] line){
			metadataParser.set(line);
			if(metadataParser.terms()!=2){return;}
			metadataParser.setBounds(0);
			final int a=metadataParser.a(), b=metadataParser.b();
			final int key=equalIgnoreCase(line, a, b, "##ploidy") ? 1 :
				equalIgnoreCase(line, a, b, "##properPairRate") ? 2 :
				equalIgnoreCase(line, a, b, "##totalQualityAvg") ? 3 :
				equalIgnoreCase(line, a, b, "##mapqAvg") ? 4 :
				equalIgnoreCase(line, a, b, "##readLengthAvg") ? 5 : 0;
			if(key==0){return;}
			final String value=metadataParser.parseString(1);
			try{
				if(key==1){ploidy=Integer.parseInt(value);}
				else if(key==2){pairingRate=Double.parseDouble(value);}
				else if(key==3){qualityAvg=Double.parseDouble(value);}
				else if(key==4){mapqAvg=Double.parseDouble(value);}
				else{readLengthAvg=Double.parseDouble(value);}
			}catch(NumberFormatException e){fail("Malformed BBTools numeric summary metadata", lines);}
		}
		private final LineParser1 metadataParser=new LineParser1('=');
		private long lines, records, headers, blanks, errors, firstErrorLine, lastErrorLine=-1;
		private String firstError, version="unknown";
		private boolean readError;
		private int ploidy=-1;
		private double pairingRate=-1, qualityAvg=-1, mapqAvg=-1, readLengthAvg=-1;
	}

	/** Immutable scan result; passing does not imply complete VCF conformance. */
	public static final class Result {
		private Result(final Scan s, final int count){
			lines=s.lines; records=s.records; headers=s.headers; blanks=s.blanks; errors=s.errors;
			firstErrorLine=s.firstErrorLine; firstError=s.firstError; version=s.version; readError=s.readError;
			columns=count; ploidy=s.ploidy; pairingRate=s.pairingRate; qualityAvg=s.qualityAvg;
			mapqAvg=s.mapqAvg; readLengthAvg=s.readLengthAvg;
		}
		/** Whether all documented structural checks and reading succeeded. */
		public boolean passed(){return errors==0 && !readError;}
		/** Physical line classifications and distinct failing line/EOF positions. */
		public final long lines, records, headers, blanks, errors, firstErrorLine;
		/** First bounded diagnostic and declared supported version, if found. */
		public final String firstError, version;
		/** Read/decompression/close error flag. */
		public final boolean readError;
		/** Declared column count and optional BBTools summary metadata. */
		public final int columns, ploidy;
		/** Optional metadata values, negative when absent. */
		public final double pairingRate, qualityAvg, mapqAvg, readLengthAvg;
	}

	private static final String FILEFORMAT="##fileformat=";
	private static final String VERSION_PREFIX="##fileformat=VCFv4.";
	private static final String[] HEADER_FIELDS={"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
}
