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
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Counts patterns in Accessions.
 * Handles hashing for Accession to TaxID lookups.
 * @author Brian Bushnell
 * @date May 9, 2018
 *
 */
public class AnalyzeAccession_ST {
	
	/**
	 * Program entry point for accession pattern analysis.
	 * Initializes timer, processes command-line arguments, and executes analysis.
	 * @param args Command-line arguments including input files and options
	 */
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
	
	/**
	 * Constructor that parses command-line arguments and configures analysis.
	 * Sets up input/output files, processing limits, and validation options.
	 * @param args Command-line arguments array containing file paths and parameters
	 */
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
	
	/**
	 * Main processing method that analyzes all input files.
	 * Processes each input file to count accession patterns, then outputs
	 * results with pattern statistics and combinatorial analysis.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){

		for(FileFormat ffin : ffina){
			process_inner(ffin);
		}
		
		if(ffout!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout);
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
			bsw.start();
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
	
	/**
	 * Processes a single input file for accession pattern analysis.
	 * Reads lines from the input file and extracts patterns using the increment method.
	 * @param ffin FileFormat object representing the input file to process
	 */
	void process_inner(FileFormat ffin){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		byte[] line=bf.nextLine();
		StringBuilder buffer=new StringBuilder(32);
		
		for(int lineNum=0; line!=null; lineNum++){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				assert((lineNum==0)==(Tools.startsWith(line, "accession"))) : "Line "+lineNum+": "+new String(line);
//				final boolean valid=(line[0]!='#');
				
				if(true){
					linesOut++;
					bytesOut+=(line.length+1);
					increment(line, buffer);
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
	}
	
	/**
	 * Extracts and counts accession patterns from a single line.
	 * Remaps characters to pattern symbols (D=digit, L=letter, -=separator)
	 * and increments the count for the resulting pattern.
	 *
	 * @param line Raw line bytes containing the accession number
	 * @param buffer StringBuilder for temporary pattern construction
	 */
	void increment(byte[] line, StringBuilder buffer){
		buffer.setLength(0);
		for(int i=0; i<line.length; i++){
			final byte b=line[i];
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		String key=buffer.toString();
		StringNum value=countMap.get(key);
		if(value!=null){value.increment();}
		else{countMap.put(key, new StringNum(key, 1));}
	}
	
	/**
	 * Calculates the number of possible combinations for a given pattern string.
	 * Each 'D' contributes 10 possibilities (0-9), each 'L' contributes 26 (A-Z).
	 * @param s Pattern string with 'D' for digits and 'L' for letters
	 * @return Number of possible combinations, or Long.MAX_VALUE if overflow
	 */
	public static long combos(String s){
		double combos=1;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? Long.MAX_VALUE : (long)Math.ceil(combos));
	}
	
	/**
	 * Calculates the number of possible combinations for a given pattern byte array.
	 * Each 'D' contributes 10 possibilities (0-9), each 'L' contributes 26 (A-Z).
	 * @param s Pattern byte array with 'D' for digits and 'L' for letters
	 * @return Number of possible combinations, or -1 if overflow
	 */
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
	
	/**
	 * Loads pattern code mapping from a file for efficient digitization.
	 * Creates a map from pattern strings to integer codes, calculating bit
	 * requirements and maximum pattern combinations.
	 *
	 * @param fname Filename containing pattern definitions
	 * @return HashMap mapping pattern strings to their integer codes
	 */
	public static HashMap<String, Integer> loadCodeMap(String fname){
		assert(codeMap==null);
		TextFile tf=new TextFile(fname);
		ArrayList<String> list=new ArrayList<String>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(!line.startsWith("#")){
				String[] split=line.split("\t");
				list.add(split[0]);
			}
		}
		HashMap<String, Integer> map=new HashMap<String, Integer>(list.size()*3);
		codeBits=(int)Math.ceil(Tools.log2(list.size()));
		final int patternBits=63-codeBits;
		final long maxCombos=((1L<<(patternBits-1))-1);
		for(int i=0; i<list.size(); i++){
			String s=list.get(i);
			longestPattern=Tools.max(longestPattern, s.length());
			long combos=combos(s);
			if(combos<0 || combos>=maxCombos){map.put(s, -1);}
			else{map.put(s, i);}
		}
		codeMap=map;
		return map;
	}
	
	/**
	 * Converts an accession string to a compact long integer representation.
	 * Combines the pattern code with the actual values using bit shifting
	 * for efficient storage and lookup operations.
	 *
	 * @param s Accession string to convert
	 * @return Long integer representation, or negative values for errors
	 */
	public static long digitize(String s){
		String pattern=remap(s);
		Integer code=codeMap.get(pattern);
		if(code==null){return -2;}
		if(code.intValue()<0){return -1;}
		
		long number=0;
		for(int i=0; i<pattern.length(); i++){
			char c=s.charAt(i);
			char p=pattern.charAt(i);
			if(p=='-'){
				//do nothing
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=(number<<codeBits)+code;
		return number;
	}
	
	/**
	 * Converts an accession byte array to a compact long integer representation.
	 * Combines the pattern code with the actual values using bit shifting
	 * for efficient storage and lookup operations.
	 *
	 * @param s Accession byte array to convert
	 * @return Long integer representation, or negative values for errors
	 */
	public static long digitize(byte[] s){
		String pattern=remap(s);
		Integer code=codeMap.get(pattern);
		if(code==null){return -2;}
		if(code.intValue()<0){return -1;}
		
		long number=0;
		for(int i=0; i<pattern.length(); i++){
			byte c=s[i];
			char p=pattern.charAt(i);
			if(p=='-'){
				//do nothing
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=(number<<codeBits)+code;
		return number;
	}
	
	/**
	 * Remaps an accession string to its pattern representation.
	 * Converts letters to 'L', digits to 'D', and separators to '-'.
	 * @param s Input accession string
	 * @return Pattern string using D/L/- notation
	 */
	public static String remap(String s){
		ByteBuilder buffer=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		return buffer.toString();
	}
	
	/**
	 * Remaps an accession byte array to its pattern representation.
	 * Converts letters to 'L', digits to 'D', and separators to '-'.
	 * @param s Input accession byte array
	 * @return Pattern string using D/L/- notation
	 */
	public static String remap(byte[] s){
		ByteBuilder buffer=new ByteBuilder(s.length);
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		return buffer.toString();
	}
	
	/*--------------------------------------------------------------*/
	
	/** List of input file paths to process */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path for pattern analysis results */
	private String out=null;
	
	/*--------------------------------------------------------------*/

	/** Map storing pattern counts for frequency analysis */
	private HashMap<String, StringNum> countMap=new HashMap<String, StringNum>();
	/** Global map from pattern strings to integer codes for digitization */
	public static HashMap<String, Integer> codeMap;
	/** Number of bits required to store pattern codes */
	private static int codeBits=-1;
	/** Length of the longest pattern found in the code map */
	private static int longestPattern=-1;
	
	/** Total number of lines processed from input files */
	private long linesProcessed=0;
	/** Number of valid lines written to output */
	private long linesOut=0;
	/** Total bytes read from input files */
	private long bytesProcessed=0;
	/** Total bytes written to output files */
	private long bytesOut=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	/** Array of input file format objects */
	private final FileFormat[] ffina;
	/** Output file format object */
	private final FileFormat ffout;
	
	/** Character remapping array for converting accessions to patterns */
	private static final byte[] remap=makeRemap();
	
	/**
	 * Creates the character remapping array for pattern extraction.
	 * Maps A-Z and a-z to 'L', 0-9 to 'D', underscore and hyphen to '-',
	 * and all other characters to '?'.
	 * @return Byte array for character-to-pattern mapping
	 */
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
	
	/** Output stream for status messages and results */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
