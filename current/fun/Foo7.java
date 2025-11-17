package fun;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser;
import shared.LineParser1;
import shared.LineParser2;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.LongList;

/**
 * Reads a text file.
 * Prints it to another text file.
 * Filters out invalid lines and prints them to an optional third file.
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class Foo7 {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Foo7 x=new Foo7(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Foo7(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		lp=(useLP2 ? new LineParser2(delimiter) : new LineParser1(delimiter));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
//		parser.out1="stdout";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("lp2")){
				useLP2=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(bf, bsw, bswInvalid);
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		
//		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
//		
//		outstream.println();
//		outstream.println("Valid Lines:       \t"+linesOut);
//		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core processing logic that reads file records, extracts size/timestamp
	 * data, and generates percentile statistics. Outputs various file size
	 * distributions and access time statistics.
	 * @param bf Input file reader
	 * @param bsw Output writer (currently unused)
	 * @param bswInvalid Invalid record writer (currently unused)
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){
		
		Timer t=new Timer();
		
		byte[] line=bf.nextLine();
		LongList list=new LongList(100000000);
		
		long sum=0;
		while(line!=null){
			linesProcessed++;
			bytesProcessed+=(line.length+1);
			long size=processLine(line, list);
			if(size>=0){sum+=size;}
			if(linesProcessed>=maxLines) {break;}
			line=bf.nextLine();
		}
		
		t.stop("parsing: ");
		t.start();
		
		final long tebi=1024L*1024L*1024L*1024L;
		final long tera=1000L*1000L*1000L*1000L;
		
		{
			list.sort();
			t.stop("sorting1: ");
			t.start();
			final int[] idxArray=new int[] {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
			final int[] pairArray=idxArray.clone();
			for(int i=0; i<idxArray.length; i++) {pairArray[i]=(int)(idxArray[i]*.01*list.size);}
			long tsum=0;
			
			for(int i=0, nextIdx=0, nextPair=pairArray[0]; i<list.size && nextIdx<idxArray.length; i++) {
				final long v=list.get(i);
				final long size=getSize(v), time=getTime(v);
				tsum+=size;
				if(i>=nextPair) {
					String s=idxArray[nextIdx]+" percent of files have not been accessed since: "+
							timeString(time*1000)+" ("+(tsum/tebi)+" tebibytes)";
					System.out.println(s);
					nextIdx++;
					if(nextIdx<idxArray.length) {nextPair=pairArray[nextIdx];}
				}
			}
			t.stop("stats1: ");
			t.start();
		}
		
		for(int i=0; i<list.size; i++) {list.array[i]=getSize(list.array[i]);}
		t.stop("recoding: ");
		t.start();
		list.sort();
		t.stop("sorting2: ");
		t.start();
		long mean=sum/list.size;
		long median=list.get((int)(list.size*0.5));
		System.out.println("total size: \t"+(sum/tera)+" TB \t("+sum+")"+"\t"+"("+((sum/tebi))+" tebibytes)");
		System.out.println("mean size:  \t"+mean+" bytes");
		System.out.println("P50 size:   \t"+median+" bytes");
		System.out.println("P80 size:   \t"+list.get((int)(list.size*0.8))+" bytes");
		System.out.println("P90 size:   \t"+list.get((int)(list.size*0.9))+" bytes");
		System.out.println("P95 size:   \t"+list.get((int)(list.size*0.95))+" bytes");
		t.stop("stats2: ");
		t.start();
	}
	
	/**
	 * Process a single line from input file. Parses file size, type, and
	 * timestamp. Only processes file records (type 'F'), ignoring others.
	 * @param line Raw line bytes from input
	 * @param list Collection to store processed timestamp/size combinations
	 * @return File size if valid file record, -1 if invalid or non-file
	 */
	long processLine(byte[] line, LongList list) {
		lp.set(line, 11);
		
		long size=lp.parseLong(3);
		byte type=lp.parseByte(6, 0);
		if(type!='F') {return -1;}
		
		long time=lp.parseLong(11);
		
		linesOut++;
		list.add(combine(time, size));
		return size;
	}
	
	/**
	 * Creates and starts a ByteStreamWriter for the given FileFormat.
	 * @param ff FileFormat specification, may be null
	 * @return Started ByteStreamWriter or null if ff is null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/**
	 * Formats Unix timestamp as human-readable PST date/time string.
	 * @param time Unix timestamp in milliseconds
	 * @return Formatted date string in PST timezone
	 */
	static String timeString(long time){
		SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
//		sdf.setTimeZone(TimeZone.getDefault());
		return sdf.format(new Date(time));
	}
	
	/** Number of bits used for compressed size storage */
	static final int LOWER_BITS=31;
	/** Number of mantissa bits for size compression */
	static final int MANTISSA_BITS=24;
	/** Number of exponent bits for size compression */
	static final int EXP_BITS=LOWER_BITS-MANTISSA_BITS;
	/** Number of bits available for timestamp storage */
	static final int UPPER_BITS=64-MANTISSA_BITS;
	/** Bit mask for extracting compressed size from combined value */
	static final long LOWER_MASK=~((-1L)<<LOWER_BITS);
	/** Bit mask for mantissa extraction during compression */
	static final long MANTISSA_MASK=~((-1L)<<MANTISSA_BITS);
	/**
	 * Compresses large file sizes using floating-point representation.
	 * Values <= 24-bit mantissa are stored directly; larger values use
	 * mantissa + exponent encoding to fit in lower 31 bits.
	 * @param raw Uncompressed file size
	 * @return Compressed size representation
	 */
	static final long compress(long raw) {
		if(raw<=MANTISSA_MASK){return raw;}
		int leading=Long.numberOfLeadingZeros(raw);
		int exp=UPPER_BITS-leading;
		assert(exp>=1);
		return (raw>>>exp)|(exp<<MANTISSA_BITS);
	}
	/**
	 * Decompresses file sizes encoded by compress(). Reverses the
	 * mantissa + exponent encoding for large values.
	 * @param f Compressed size representation
	 * @return Original file size (approximate for large values)
	 */
	static final long decompress(long f) {
		if(f<=MANTISSA_MASK){return f;}
		int exp=(int)(f>>>MANTISSA_BITS);
		assert(exp>=1);
		return (f&MANTISSA_MASK)<<exp;
	}
	/**
	 * Combines timestamp and compressed file size into single long value.
	 * Places timestamp in upper 33 bits, compressed size in lower 31 bits.
	 * @param time Unix timestamp
	 * @param size File size (will be compressed)
	 * @return Combined timestamp/size representation
	 */
	static final long combine(long time, long size) {
		return (time<<LOWER_BITS) | compress(size);
	}
	/**
	 * Extracts timestamp from combined timestamp/size value.
	 * @param combined Value created by combine()
	 * @return Original timestamp
	 */
	static final long getTime(long combined) {
		return combined>>>LOWER_BITS;
	}
	/**
	 * Extracts and decompresses file size from combined timestamp/size value.
	 * @param combined Value created by combine()
	 * @return Original file size (approximate for compressed values)
	 */
	static final long getSize(long combined) {
		return decompress(combined&LOWER_MASK);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Junk output file path */
	private String outInvalid=null;
	
	/** Whether to use LineParser2 instead of LineParser1 */
	private boolean useLP2=false;

	/** Field delimiter character for parsing input lines */
	private static final byte delimiter=(byte)'|';
	/** Parser for delimited input lines */
	private final LineParser lp;
	
	/*--------------------------------------------------------------*/
	
	/** Total number of input lines processed */
	private long linesProcessed=0;
	/** Number of valid output lines processed */
	private long linesOut=0;
	/** Total bytes read from input */
	private long bytesProcessed=0;
	/** Total bytes written to output */
	private long bytesOut=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input File */
	private final FileFormat ffin1;
	/** Output File */
	private final FileFormat ffout1;
	/** Optional Output File for Junk */
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
