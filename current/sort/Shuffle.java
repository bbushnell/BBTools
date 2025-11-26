package sort;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SamReadInputStream;
import structures.ListNum;
import tracker.ReadStats;


/**
 * Randomizes the order of reads.
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */
public class Shuffle {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program entry point that creates a Shuffle instance and executes processing.
	 * @param args Command-line arguments specifying input files, output files, and mode */
	public static void main(String[] args){
		Timer t=new Timer();
		Shuffle x=new Shuffle(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs a Shuffle instance and parses command-line arguments.
	 * Sets up input/output file formats, processing mode, and validation parameters.
	 * Configures threading options and determines paired vs unpaired processing.
	 * @param args Command-line arguments including file paths and mode settings
	 */
	public Shuffle(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, printClass ? getClass() : null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.SET_FROM_OK=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		
		int mode_=SHUFFLE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("shuffle")){
				mode_=SHUFFLE;
			}else if(a.equals("name")){
				mode_=SORT_NAME;
			}else if(a.equals("coordinate")){
				mode_=SORT_COORD;
			}else if(a.equals("sequence")){
				mode_=SORT_SEQ;
			}else if(a.equals("id")){
				mode_=SORT_ID;
			}else if(a.equals("mode")){
				if(b==null){
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}else if(b.equals("shuffle")){
					mode_=SHUFFLE;
				}else if(b.equals("name")){
					mode_=SORT_NAME;
				}else if(b.equals("coordinate")){
					mode_=SORT_COORD;
				}else if(b.equals("sequence")){
					mode_=SORT_SEQ;
				}else if(b.equals("id")){
					mode_=SORT_ID;
				}else{
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		assert(mode>=1 && mode<=5) : "mode must be shuffle, name, coordinate, sequence, or id.";
		
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

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		useSharedHeader=(ffout1.samOrBam() && ffout1!=null && ffout1.samOrBam());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that reads all sequences, applies sorting/shuffling, and writes output.
	 * Loads entire input into memory, performs the requested operation (shuffle or sort),
	 * then writes the reordered sequences to output files.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		ArrayList<Read> bigList=new ArrayList<Read>(65530);
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					bigList.add(r1);
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStream(cris);
		errorState|=ReadStats.writeAll();
		
		if(mode==SHUFFLE){
			Collections.shuffle(bigList);
		}else if(mode==SORT_NAME){
			Shared.sort(bigList, ReadComparatorName.comparator);
		}else if(mode==SORT_SEQ){
			Shared.sort(bigList, ReadComparatorTopological.comparator);
		}else if(mode==SORT_COORD){
			Shared.sort(bigList, new ReadComparatorMapping());
		}else if(mode==SORT_ID){
			Shared.sort(bigList, ReadComparatorID.comparator);
		}else{
			assert(false) : "No mode set.";
		}
		
		if(ffout1!=null){
			final ByteStreamWriter bsw1, bsw2;
			if(ffout1!=null){
				bsw1=new ByteStreamWriter(ffout1);
				bsw1.start();
				if(useSharedHeader){writeHeader(bsw1);}
			}else{bsw1=null;}
			if(ffout2!=null){
				bsw2=new ByteStreamWriter(ffout2);
				bsw2.start();
				if(useSharedHeader){writeHeader(bsw1);}
			}else{bsw2=null;}
			final boolean b=(bsw2==null);
			for(int i=0, lim=bigList.size(); i<lim; i++){
				final Read r1=bigList.set(i, null);
				final Read r2=r1.mate;
				bsw1.println(r1, b);
				if(r2!=null && !b){bsw2.println(r2);}
			}
			if(bsw1!=null){errorState|=bsw1.poisonAndWait();}
			if(bsw2!=null){errorState|=bsw2.poisonAndWait();}
		}
		
		t.stop();
		
		if(showSpeed){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Writes SAM/BAM header lines to the output stream if present.
	 * @param bsw Output stream writer for header data */
	private void writeHeader(ByteStreamWriter bsw){
		ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
		if(list==null){
			System.err.println("Header was null.");
		}else{
			for(byte[] line : list){
				bsw.print(line).print('\n');
			}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Threading wrapper for running shuffle operations in parallel.
	 * Encapsulates input/output paths and processing mode for execution
	 * in a separate thread with thread pool management.
	 */
	public static class ShuffleThread extends Thread{
		
		/**
		 * Creates a new shuffle thread with specified input/output files and processing mode.
		 *
		 * @param in1_ Primary input file path
		 * @param in2_ Secondary input file path for paired reads (may be null)
		 * @param out1_ Primary output file path
		 * @param out2_ Secondary output file path for paired reads (may be null)
		 * @param mode_ Processing mode (shuffle, sort by name, sequence, coordinates, or ID)
		 * @param ow_ Whether to overwrite existing output files
		 */
		public ShuffleThread(String in1_, String in2_, String out1_, String out2_, int mode_, boolean ow_){
			in1=in1_;
			in2=in2_;
			out1=out1_;
			out2=out2_;
			mode=mode_;
			ow=ow_;
		}
		
		@Override
		public synchronized void start(){
			addThread(1);
			super.start();
		}
		
		@Override
		public void run(){
			ArrayList<String> list=new ArrayList<String>();
			if(in1!=null){list.add("in1="+in1);}
			if(in2!=null){list.add("in1="+in2);} //Possible bug: should be "in2="+in2
			if(out1!=null){list.add("out1="+out1);}
			if(out2!=null){list.add("out2="+out2);}
			list.add("mode="+MODES[mode]);
			list.add("ow="+ow);
			try{
				Shuffle.main(list.toArray(new String[0]));
			}catch(Throwable e){
				System.err.println("Failed to shuffle "+in1+"\nException:"+e+"\n");
			}
			addThread(-1);
		}
		
		/** Secondary input file path for this thread */
		/** Primary input file path for this thread */
		final String in1, in2;
		/** Secondary output file path for this thread */
		/** Primary output file path for this thread */
		final String out1, out2;
		/** Processing mode for this thread */
		final int mode;
		/** Whether this thread should overwrite output files */
		final boolean ow;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path for paired reads */
	private String in2=null;
	
	/** Quality file path for primary input */
	private String qfin1=null;
	/** Quality file path for secondary input */
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path for paired reads */
	private String out2=null;

	/** Quality file path for primary output */
	private String qfout1=null;
	/** Quality file path for secondary output */
	private String qfout2=null;
	
	/** File extension for input files */
	private String extin=null;
	/** File extension for output files */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/** Processing mode: shuffle, sort by name, sequence, coordinates, or ID */
	private final int mode;
	
	/*--------------------------------------------------------------*/
	
	/** File format handler for primary input */
	private final FileFormat ffin1;
	/** File format handler for secondary input */
	private final FileFormat ffin2;

	/** File format handler for primary output */
	private final FileFormat ffout1;
	/** File format handler for secondary output */
	private final FileFormat ffout2;
	
	/** Whether to use shared SAM/BAM headers across outputs */
	private final boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum number of concurrent shuffle threads allowed */
	private static int maxShuffleThreads=1;
	/** Current number of active shuffle threads */
	private static int currentShuffleThreads=0;
	
	/** Sets the maximum number of concurrent shuffle threads allowed.
	 * @param x Maximum thread count (must be greater than 0) */
	public static void setMaxThreads(final int x){
		assert(x>0);
		synchronized(SHUFFLE_LOCK){
			maxShuffleThreads=x;
		}
	}
	
	/**
	 * Manages thread pool by adding or removing threads with blocking when at capacity.
	 * Waits if adding threads would exceed the maximum, decrements when removing threads.
	 * @param x Number of threads to add (positive) or remove (negative)
	 * @return Current number of active threads after the operation
	 */
	public static int addThread(final int x){
		synchronized(SHUFFLE_LOCK){
			while(x>0 && currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			currentShuffleThreads+=x;
			if(currentShuffleThreads<maxShuffleThreads){SHUFFLE_LOCK.notify();}
			return currentShuffleThreads;
		}
	}
	
	/** Blocks until all shuffle threads complete and thread count drops below maximum.
	 * Used for synchronization when coordinating multiple parallel operations. */
	public static void waitForFinish(){
		synchronized(SHUFFLE_LOCK){
			while(currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	/** Synchronization lock object for thread pool management */
	private static String SHUFFLE_LOCK=new String("SHUFFLE_LOCK");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose logging output */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	/** Whether to display processing speed statistics */
	public static boolean showSpeed=true;
	/** Whether to print class information in help messages */
	public static boolean printClass=true;

	/** Mode constant for sorting reads by numeric ID */
	/** Mode constant for sorting reads by mapping coordinates */
	/** Mode constant for sorting reads by sequence content */
	/** Mode constant for sorting reads by name */
	/** Mode constant for random shuffling of reads */
	public static final int SHUFFLE=1, SORT_NAME=2, SORT_SEQ=3, SORT_COORD=4, SORT_ID=5;
	/** String representations of available processing modes */
	public static final String[] MODES={"shuffle", "name", "sequence", "coordinate", "id"};
	
	
}
