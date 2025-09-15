package template;

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
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Jan 12, 2015
 *
 */
public abstract class BBTool_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		throw new RuntimeException("This method must be overridden.");
		/*
		//Example:
		Timer t=new Timer();
		BBTool_ST bbt=new BBTool_ST(args);
		bbt.process(t);
		*/
	}
	
	/**
	 * Constructor.  Parses argument list and sets relevant fields.
	 * Must be called by subclass.
	 * @param args Command line arguments
	 */
	public BBTool_ST(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		
		setDefaults();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(parseArgument(arg, a, b)){
				// do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
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
	}
	
	/** Must be overridden if variables are defined at bottom of class file. */
	protected abstract void setDefaults();
	
	protected void reparse(String[] args){
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parseArgument(arg, a, b)){
				// do nothing
			}
		}
	}
	
	/**
	 * Must be overridden; this body is just for example.
	 * Parses an argument from the command line.
	 * Assumed to be in "key=value" form, but this is not required.
	 * @param arg The full original argument.
	 * @param a Left hand side, to lower case.
	 * @param b Right hand side, unaltered.
	 * @return true if a matched some keyword.
	 */
	public abstract boolean parseArgument(String arg, String a, String b);
	
//	//Example
//	@Override
//	public boolean parseArgument(String arg, String a, String b){
//		if(a.equals("keepunmapped") | a.equals("ku")){
//			keepUnmapped=Parse.parseBoolean(b);
//			return true;
//		}else if(a.equals("ignorepairorder") | a.equals("ipo")){
//			usePairOrder=!Parse.parseBoolean(b);
//			return true;
//		}else if(a.equals("sorted")){
//			sorted=Parse.parseBoolean(b);
//			return true;
//		}
//		return false;
//	}
	
	/**
	 * And example of how to override parseArgument.
	 */
	private boolean parseArgument_EXAMPLE(String arg, String a, String b){
		if(true){throw new RuntimeException("parseArgument() must be overridden.");}
		
		//These are dummy values for demonstration purposes.
		//In real code they should be class fields.
		int value1;
		boolean value2;
		String value3;
		
		if(a.equals("key1")){
			value1=Parse.parseIntKMG(b);
			//do anything else necessary here
			return true;
		}else if(a.equals("key2")){
			value2=Parse.parseBoolean(b);
			//do anything else necessary here
			return true;
		}else if(a.equals("key3")){
			value3=b;
			//do anything else necessary here
			return true;
		}
		
		//There was no match to the argument
		return false;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){process(new Timer());}

	/** Create read streams and process all data */
	protected void process(Timer t){
		
		//Start the read streams
		startup();
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read streams
		processInner(cris_primary, ros_primary);
		
		//Close the read streams
		shutdown(cris_primary, ros_primary);
		
		showStats(t);
	}

	/** Create read streams */
	protected void startup(){
		startupSubclass();
		
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, out1, out2, qfout1, qfout2)){
			assert(false) : "Duplicate files.";
		}
		
		final ConcurrentReadInputStream cris;
		final ConcurrentReadOutputStream ros;
		
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		final boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader());
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		cris_primary=cris;
		ros_primary=ros;
	}
	
	/** Called before startup().
	 * Implement if necessary */
	protected abstract void startupSubclass();
	
	protected final void shutdown(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		shutdownSubclass();
		
		final boolean paired=cris.paired();
		if(verbose){outstream.println("Finished.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, ros);
	}
	
	/** Called before shutdown().
	 * Implement if necessary */
	protected abstract void shutdownSubclass();
	
	protected void showStats(final Timer t){
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		showStatsSubclass(t, readsProcessed, basesProcessed);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Called AFTER showStats().
	 * Implement if necessary */
	protected abstract void showStatsSubclass(final Timer t, long readsIn, long basesIn);
	
	/** Iterate through the reads.
	 * This may optionally be overridden. */
	protected void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
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
					
					boolean keep=processReadPair(r1, r2);
					if(!keep){reads.set(idx, null);}
					
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Process a single read pair.  Must be overidden.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	protected abstract boolean processReadPair(final Read r1, final Read r2);
	
	/** For sam format */
	protected abstract boolean useSharedHeader();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	protected String in1=null;
	protected String in2=null;
	
	protected String qfin1=null;
	protected String qfin2=null;

	protected String out1=null;
	protected String out2=null;

	protected String qfout1=null;
	protected String qfout2=null;
	
	protected String extin=null;
	protected String extout=null;
	
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadInputStream cris_primary;
	private ConcurrentReadOutputStream ros_primary;
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	protected final FileFormat ffin1;
	protected final FileFormat ffin2;

	protected FileFormat ffout1;
	protected FileFormat ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	protected PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	protected boolean overwrite=true;
	protected boolean append=false;
	
}
