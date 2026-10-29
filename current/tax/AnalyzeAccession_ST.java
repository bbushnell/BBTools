package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Single-threaded predecessor to {@link AnalyzeAccession}.  Counts accession
 * pattern frequencies, calculates combinatorial possibilities, and keeps a
 * minimal static pattern table for legacy disk-key digitization.
 * @author Brian Bushnell
 * @date May 9, 2018
 */
public class AnalyzeAccession_ST{
	
	/** Runs single-threaded accession-pattern analysis from the command line.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AnalyzeAccession_ST x=new AnalyzeAccession_ST(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses command-line arguments and prepares input and output formats.
	 * @param args Command-line arguments */
	public AnalyzeAccession_ST(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				if(b==null){in.clear();}
				else{
					String[] split2=b.split(",");
					for(String s2 : split2){
						in.add(s2);
					}
				}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffina=new FileFormat[in.size()];
		for(int i=0; i<in.size(); i++){
			ffina[i]=FileFormat.testInput(in.get(i), FileFormat.TXT, null, true, false);
		}
	}
	
	/** Counts accession patterns from all inputs, optionally writes a sorted summary,
	 * and throws if the reader or output writer reported an error.
	 * @param t Runtime timer started by the caller */
	void process(Timer t){

		for(FileFormat ffin : ffina){
			process_inner(ffin);
		}
		
		if(ffout!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout);
			/* FIXED [tax/AnalyzeAccession_ST#001]: start() must precede any
			 * print; addJob requires a started writer once enough pattern
			 * lines fill its internal buffer. */
			bsw.start();
			bsw.println("#Pattern\tCount\tCombos\tBits");
			ArrayList<StringNum> list=new ArrayList<StringNum>();
			list.addAll(countMap.values());
			Collections.sort(list);
			Collections.reverse(list);
			for(StringNum sn : list){
				double combos=1;
				for(int i=0; i<sn.s.length(); i++){
					char c=sn.s.charAt(i);
					if(c=='D'){combos*=10;}
					else if(c=='L'){combos*=26;}
				}
				bsw.print(sn.toString().getBytes());
				bsw.println("\t"+(long)combos+"\t"+Tools.format("%.2f", Tools.log2(combos)));
			}
			errorState|=bsw.poisonAndWait();
		}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesOut);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Counts one input file sequentially, skipping the first accession header.
	 * Always attempts to close the reader, suppressing close failures behind a
	 * primary parse failure.
	 * @param ffin Input file format descriptor */
	void process_inner(FileFormat ffin){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		Throwable failure=null;
		try{
			byte[] line=bf.nextLine();
			StringBuilder buffer=new StringBuilder(32);

			for(int lineNum=0; line!=null; lineNum++){
				if(line.length>0){
					if(maxLines>0 && linesProcessed>=maxLines){break;}
					linesProcessed++;
					bytesProcessed+=(line.length+1);

					assert((lineNum==0)==(Tools.startsWith(line, "accession"))) :
						"Line "+lineNum+": "+new String(line);
//				final boolean valid=(line[0]!='#');

					if(lineNum>0 || !Tools.startsWith(line, "accession")){
						linesOut++;
						bytesOut+=(line.length+1);
						increment(line, buffer);
					}
				}
				line=bf.nextLine();
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			try{
				errorState|=bf.close();
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
	}
	
	/** Remaps one accession token into a pattern and increments its count.
	 * @param line Input line beginning with the accession token
	 * @param buffer Reused scratch buffer */
	void increment(byte[] line, StringBuilder buffer){
		buffer.setLength(0);
		for(int i=0; i<line.length; i++){
			final byte b=line[i];
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		if(buffer.length()>=3 && line.length>=3
				&& (line[0]=='N' || line[0]=='n') && (line[1]=='Z' || line[1]=='z')
				&& (line[2]=='_' || line[2]=='-')
				&& buffer.charAt(0)=='L' && buffer.charAt(1)=='L' && buffer.charAt(2)=='-'){
			buffer.setCharAt(0, 'N');
			buffer.setCharAt(1, 'Z');
		}
		String key=buffer.toString();
		StringNum value=countMap.get(key);
		if(value!=null){value.increment();}
		else{countMap.put(key, new StringNum(key, 1));}
	}
	
	/** Counts possible concrete accessions represented by a pattern.
	 * @param s Pattern using {@code D} for digits and {@code L} for letters
	 * @return Product of all variable positions, clamped to Long.MAX_VALUE */
	public static long combos(String s){
		double combos=1;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? Long.MAX_VALUE : (long)Math.ceil(combos));
	}
	
	/** Counts possible concrete accessions represented by a byte pattern.
	 * @param s Pattern using {@code D} for digits and {@code L} for letters
	 * @return Product of all variable positions, or -1 if it overflows a long */
	public static long combos(byte[] s){
		double combos=1;
		for(int i=0; i<s.length; i++){
			byte c=s[i];
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? -1 : (long)Math.ceil(combos));
	}
	
	/*--------------------------------------------------------------*/
	
	/** Loads tab-delimited accession patterns and prepares legacy digitization
	 * tables for this class' direct map lookup.  The pattern file is closed
	 * before global table construction or publication; close failures identify
	 * the file and are suppressed behind a primary read failure.
	 * @param fname Pattern table path
	 * @return Pattern to dense-index map, with -1 for overflow patterns */
	public static HashMap<String, Integer> loadCodeMap(String fname){
		assert(codeMap==null);
		TextFile tf=new TextFile(fname);
		ArrayList<String> list=new ArrayList<String>();
		Throwable failure=null;
		try{
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(!line.startsWith("#")){
					String[] split=line.split("\t");
					list.add(split[0]);
				}
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			try{
				if(tf.close()){throw new RuntimeException("Error closing "+fname);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}

		final int numPatterns=list.size();
		long[] combosArr=new long[numPatterns];
		for(int i=0; i<numPatterns; i++){
			String s=list.get(i);
			longestPattern=Tools.max(longestPattern, s.length());
			combosArr[i]=combos(s);
		}

		boolean[] isTable1=new boolean[numPatterns];
		boolean changed=true;
		while(changed){
			changed=false;
			int count0=0;
			for(int i=0; i<numPatterns; i++){
				if(!isTable1[i]){count0++;}
			}
			if(count0==0){break;}
			int cb=Math.max(1, (int)Math.ceil(Tools.log2(count0)));
			long maxCombos0=(1L<<(62-cb));
			for(int i=0; i<numPatterns; i++){
				if(!isTable1[i] && (combosArr[i]<0 || combosArr[i]>=maxCombos0)){
					isTable1[i]=true;
					changed=true;
				}
			}
		}

		int count0=0, count1=0;
		for(int i=0; i<numPatterns; i++){
			if(isTable1[i]){count1++;}
			else{count0++;}
		}

		codeBits0=(count0>0) ? Math.max(1, (int)Math.ceil(Tools.log2(count0))) : 1;
		codeBits1=(count1>0) ? Math.max(1, (int)Math.ceil(Tools.log2(count1))) : 1;
		long maxCombos1=(1L<<(62-codeBits1));

		HashMap<String, Integer> map=new HashMap<String, Integer>(numPatterns*3);
		patternShift=new int[numPatterns];
		patternLowBits=new int[numPatterns];

		int code0=0, code1=0;
		for(int i=0; i<numPatterns; i++){
			String s=list.get(i);
			if(isTable1[i]){
				if(combosArr[i]<0 || combosArr[i]>=maxCombos1){
					map.put(s, -1);
					patternShift[i]=-1;
				}else{
					int pc=code1++;
					patternShift[i]=codeBits1+1;
					patternLowBits[i]=(pc<<1)|1;
					map.put(s, i);
				}
			}else{
				int pc=code0++;
				patternShift[i]=codeBits0+1;
				patternLowBits[i]=(pc<<1)|0;
				map.put(s, i);
			}
		}

		codeMap=map;
		return map;
	}
	
	/** Converts an ASCII accession string to this class' compact positive key
	 * after {@link #loadCodeMap(String)} has initialized the pattern table.
	 * @param s Accession string; parsing stops at space, tab, or '.'
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern */
	public static long digitize(String s){
		String pattern=remap(s);
		Integer idx=codeMap.get(pattern);
		if(idx==null){return -2;}
		if(idx.intValue()<0){return -1;}

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			char c=s.charAt(i);
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}

	/** Converts ASCII accession bytes to this class' compact positive key after
	 * {@link #loadCodeMap(String)} has initialized the pattern table.
	 * @param s Accession bytes; parsing stops at space, tab, or '.'
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern */
	public static long digitize(byte[] s){
		String pattern=remap(s);
		Integer idx=codeMap.get(pattern);
		if(idx==null){return -2;}
		if(idx.intValue()<0){return -1;}

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			byte c=s[i];
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : new String(s);
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}
	
	/** Converts an accession prefix into its D/L/-/?/NZ pattern.
	 * @param s ASCII accession string; parsing stops at space, tab, or '.'
	 * @return Remapped pattern */
	public static String remap(String s){
		ByteBuilder buffer=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		String result=buffer.toString();
		if(result.length()>=3 && s.length()>=3
				&& (s.charAt(0)=='N' || s.charAt(0)=='n')
				&& (s.charAt(1)=='Z' || s.charAt(1)=='z')
				&& (s.charAt(2)=='_' || s.charAt(2)=='-')
				&& result.startsWith("LL-")){
			result="NZ"+result.substring(2);
		}
		return result;
	}

	/** Converts accession bytes into their D/L/-/?/NZ pattern.
	 * @param s ASCII accession bytes; parsing stops at space, tab, or '.'
	 * @return Remapped pattern */
	public static String remap(byte[] s){
		ByteBuilder buffer=new ByteBuilder(s.length);
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		String result=buffer.toString();
		if(result.length()>=3 && s.length>=3
				&& (s[0]=='N' || s[0]=='n')
				&& (s[1]=='Z' || s[1]=='z')
				&& (s[2]=='_' || s[2]=='-')
				&& result.startsWith("LL-")){
			result="NZ"+result.substring(2);
		}
		return result;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Input filenames from {@code in=} or bare existing-file arguments. */
	private ArrayList<String> in=new ArrayList<String>();
	/** Optional summary output filename. */
	private String out=null;
	
	/*--------------------------------------------------------------*/

	/** Aggregated pattern counts. */
	private HashMap<String, StringNum> countMap=new HashMap<String, StringNum>();
	/** Pattern-to-index table populated by {@link #loadCodeMap(String)}. */
	public static HashMap<String, Integer> codeMap;
	/** Bits assigned to table-0 pattern codes. */
	static int codeBits0=-1;
	/** Bits assigned to table-1 pattern codes. */
	static int codeBits1=-1;
	/** Per-pattern bit shift used to leave room for the pattern code. */
	static int[] patternShift;
	/** Per-pattern low bits containing the pattern code and table flag. */
	static int[] patternLowBits;
	/** Longest loaded pattern length. */
	private static int longestPattern=-1;
	
	/** Total nonempty input lines seen. */
	private long linesProcessed=0;
	/** Total nonempty rows counted as patterns. */
	private long linesOut=0;
	/** Approximate input bytes including line terminators. */
	private long bytesProcessed=0;
	/** Approximate counted bytes including line terminators. */
	private long bytesOut=0;
	
	/** Maximum lines to process; nonpositive command-line values mean unlimited. */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	/** Input file descriptors. */
	private final FileFormat[] ffina;
	/** Optional summary output descriptor. */
	private final FileFormat ffout;
	
	/** ASCII byte-to-pattern remap table. */
	private static final byte[] remap=makeRemap();
	
	/** Builds the ASCII byte-to-pattern lookup table used by the remap helpers.
	 * @return Remap table for byte values below 128 */
	private static byte[] makeRemap(){
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)'?');
		for(int i='A'; i<='Z'; i++){array[i]='L';}
		for(int i='a'; i<='z'; i++){array[i]='L';}
		for(int i='0'; i<='9'; i++){array[i]='D';}
		array['_']=array['-']='-';
		return array;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Status and progress output stream. */
	private PrintStream outstream=System.err;
	/** Enables verbose logging in this class and the delegated I/O classes. */
	public static boolean verbose=false;
	/** True once the reader or writer reports a failure. */
	public boolean errorState=false;
	/** Allow output files to be overwritten. */
	private boolean overwrite=true;
	/** Append to existing output files instead of replacing them. */
	private boolean append=false;
	
}
