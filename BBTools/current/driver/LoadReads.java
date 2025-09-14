package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * This class loads data to see how much memory is used.
 * 
 * @author Brian Bushnell
 * @date October 28, 2016
 *
 */
public class LoadReads {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		LoadReads x=new LoadReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public LoadReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBufferLen(100);
		Shared.capBuffers(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("earlyexit")){
				earlyExit=Parse.parseBoolean(b);
			}else if(a.equals("lowcomplexity")){
				lowComplexity=Parse.parseBoolean(b);
			}else if(a.equals("gc")){
				gc=Parse.parseBoolean(b);
			}else if(a.equals("overhead")){
				overhead=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			
			extin=parser.extin;
		}
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null) : "\nin1="+in1+"\nin2="+in2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		calcMem();
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		calcMem();
		
		double[] estimates=Tools.estimateFileMemory(in1, 200, overhead, earlyExit, lowComplexity);
		long memEst0=(long)estimates[0];
		long diskEst0=(long)estimates[1];
		double memRatio0=estimates[2];
		double diskRatio0=estimates[3];
		double readEst0=estimates[4];
		
		long size=new File(in1).length();
		long usedMem=maxMem-minMem;
		double memRatio1=memBytesProcessed/(double)size;
		double memRatio=usedMem/(double)size;
		double diskRatio=diskBytesProcessed/(double)size;
		double readRatio=readsProcessed/(double)size;
		
		long overhead=usedMem-diskBytesProcessed;
		
		double mult=1.0/readsProcessed;
		double memPerRead=usedMem*mult;
//		double mem0PerRead=usedMem*mult;
		double mem1PerRead=memBytesProcessed*mult;
		double basesPerRead=basesProcessed*mult;
		double qualsPerRead=qualitiesProcessed*mult;
		double headerPerRead=headersProcessed*mult;
		double overheadPerRead=overhead*mult;
		
		
		final long afterGC;
		if(gc){
			outstream.println("Final GC.");

			System.gc();
			afterGC=Shared.memUsed();
		}else{
			afterGC=0;
		}

		outstream.println("Initial Memory:     \t"+(initialMem/1000000)+" m");
		outstream.println("Final Memory:       \t"+(finalMem/1000000)+" m");
		if(gc){outstream.println("After GC:           \t"+(afterGC/1000000)+" m");}
		outstream.println("Min Memory:         \t"+(minMem/1000000)+" m");
		outstream.println("Max Memory:         \t"+(maxMem/1000000)+" m");
		outstream.println();
		outstream.println("Memory Estimate 0:  \t"+(memEst0/1000000)+" m");
		outstream.println("Memory Estimate 1:  \t"+(memBytesProcessed/1000000)+" m");
		outstream.println("Memory:             \t"+(usedMem/1000000)+" m");
		outstream.println();
		outstream.println("Disk Estimate 0:    \t"+(diskEst0/1000000+" m"));
		outstream.println("Disk Bytes:         \t"+((diskBytesProcessed)/1000000)+" m");
		outstream.println();
		outstream.println("Memory Ratio Est 0: \t"+(Tools.format("%.2f", memRatio0)));
		outstream.println("Memory Ratio Est 1: \t"+(Tools.format("%.2f", memRatio1)));
		outstream.println("Memory Ratio:       \t"+(Tools.format("%.2f", memRatio)));
		outstream.println();
		outstream.println("Disk Ratio Est 0:   \t"+(Tools.format("%.2f", diskRatio0)));
		outstream.println("Disk Ratio:         \t"+(Tools.format("%.2f", diskRatio)));
		outstream.println();
		outstream.println("Read Estimate 0:    \t"+(Tools.format("%d", (long)Math.ceil(readEst0))));
		outstream.println("Read Ratio 1:       \t"+(Tools.format("%.2f", readRatio)));
		outstream.println("Reads:              \t"+(Tools.format("%d", readsProcessed)));
		outstream.println();
//		outstream.println("Average Memory 0:   \t"+(Tools.format("%.2f", mem0PerRead)));
		outstream.println("Average Memory 1:   \t"+(Tools.format("%.2f", mem1PerRead)));
		outstream.println("Average Memory:     \t"+(Tools.format("%.2f", memPerRead)));
		outstream.println("Average Bases:      \t"+(Tools.format("%.2f", basesPerRead)));
		outstream.println("Average Q-Scores:   \t"+(Tools.format("%.2f", qualsPerRead)));
		outstream.println("Average Header Len: \t"+(Tools.format("%.2f", headerPerRead)));
		outstream.println("Average Overhead:   \t"+(Tools.format("%.2f", overheadPerRead)));
		outstream.println();
		
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double dpnano=diskBytesProcessed/(double)(t.elapsed);
			double mpnano=memBytesProcessed/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			String dpstring=(diskBytesProcessed<100000 ? ""+diskBytesProcessed : diskBytesProcessed<100000000 ? (diskBytesProcessed/1000)+"k" : (diskBytesProcessed/1000000)+"m");
			String mpstring=(memBytesProcessed<100000 ? ""+memBytesProcessed : memBytesProcessed<100000000 ? (memBytesProcessed/1000)+"k" : (memBytesProcessed/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(dpstring.length()<8){dpstring=" "+dpstring;}
			while(mpstring.length()<8){mpstring=" "+mpstring;}
			
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			outstream.println("Disk Bytes Processed: "+dpstring+" \t"+Tools.format("%.2fm bytes/sec", dpnano*1000));
			outstream.println("Mem Bytes Processed:  "+mpstring+" \t"+Tools.format("%.2fm bytes/sec", mpnano*1000));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				storage.add(reads);
				calcMem();
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					qualitiesProcessed+=(r1.qlength()+(r2==null ? 0 : r2.qlength()));
					headersProcessed+=(r1.id.length()+(r2==null ? 0 : r2.id.length()));
					
					diskBytesProcessed+=r1.countFastqBytes();
					if(r2!=null){diskBytesProcessed+=r2.countFastqBytes();}
					
					memBytesProcessed+=r1.countBytes();
					if(r2!=null){memBytesProcessed+=r2.countBytes();}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void calcMem(){
		final long used=Shared.memUsed();
		minMem=Tools.min(used, minMem);
		maxMem=Tools.max(used, maxMem);
		finalMem=used;
		if(initialMem<0){
//			System.err.println(minMem+", "+maxMem+", "+initialMem+", "+finalMem);
			initialMem=used;
		}
//		System.err.println(minMem+", "+maxMem+", "+initialMem+", "+finalMem);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;
	
	/** Override input file extension */
	private String extin=null;
	
	private ArrayList<ArrayList<Read>> storage=new ArrayList<ArrayList<Read>>();
	
	/*--------------------------------------------------------------*/
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of quality scores processed */
	protected long qualitiesProcessed=0;
	/** Number of header characters processed */
	protected long headersProcessed=0;
	
	/** Approximate number of fastq disk bytes processed */
	protected long diskBytesProcessed=0;
	/** Approximate number of read memory bytes processed */
	protected long memBytesProcessed=0;
	
	/** Minimal observed memory usage */
	protected long minMem=Long.MAX_VALUE;
	/** Maximal observed memory usage */
	protected long maxMem=0;
	
	protected long initialMem=-1;
	protected long finalMem=-1;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private int overhead=0;
	private boolean earlyExit=false;
	private boolean gc=false;
	private boolean lowComplexity=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
