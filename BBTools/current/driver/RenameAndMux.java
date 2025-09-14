package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

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

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class RenameAndMux {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameAndMux x=new RenameAndMux(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public RenameAndMux(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=true;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_UNBGZIP=false;
		ReadWrite.setZipThreads((Shared.threads()*3+1)/4);
		
		//Create a parser object
		Parser parser=new Parser();
		
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
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){readPaths.add(x2);}
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
			
			overwrite=parser.overwrite;
			setInterleaved=parser.setInterleaved;
			
			if(parser.in1!=null){
				for(String s : parser.in1.split(",")){
					readPaths.add(s);
				}
			}

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
//		System.err.println("out1="+out1+", out2="+out2);
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(readPaths.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){throw new RuntimeException("Error - output destination is required.");}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Reset counters
		readsProcessedA.set(0);
		basesProcessedA.set(0);
		nextListNumber.set(0);
		
		//Process the read stream
		renameAndMerge_MT();
		
		final long readsProcessed=readsProcessedA.get();
		final long basesProcessed=basesProcessedA.get();
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
//	private void renameAndMerge_ST(){
//
//		long readsProcessed=0;
//		long basesProcessed=0;
//
//		FileFormat ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, false, false);
//		FileFormat ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, false, false);
//
//		for(String path : readPaths){
//
//			String in1=null, in2=null;
//
//			//Do input file # replacement
//			if(path.indexOf('#')>-1 && !new File(path).exists()){
//				in2=path.replace("#", "2");
//				in1=path.replace("#", "1");
//			}
//
//			//Adjust interleaved settings based on number of output files
//			if(!setInterleaved){
//				assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
//				if(in2!=null){ //If there are 2 input streams.
//					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
//					if(!printedInterleavedMessage){
//						outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
//						printedInterleavedMessage=true;
//					}
//				}else{ //There is one input stream.
//					if(out2!=null){
//						FASTQ.FORCE_INTERLEAVED=true;
//						FASTQ.TEST_INTERLEAVED=false;
//						if(!printedInterleavedMessage){
//							outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
//							printedInterleavedMessage=true;
//						}
//					}
//				}
//			}
//
//			//Ensure input files can be read
//			if(!Tools.testInputFiles(false, true, in1, in2)){
//				throw new RuntimeException("\nCan't read some input files.\n");  
//			}
//
//			//Ensure that no file was specified multiple times
//			if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
//				throw new RuntimeException("\nSome file names were specified multiple times.\n");
//			}
//
//			final FileFormat ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
//			final FileFormat ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
//			final ConcurrentReadInputStream cris;
//			final String core=ReadWrite.stripToCore(in1);
//			{
//				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
//				if(verbose){outstream.println("Started cris");}
//				cris.start(); //4567
//			}
//			boolean paired=cris.paired();
//			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//			final ConcurrentReadOutputStream ros;
//
//			final int buff=4;
//
//			if(cris.paired()){
//				outstream.println("Writing interleaved.");
//			}
//
//			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, buff, null, false);
//			ros.start();
//
//			{
//
//				ListNum<Read> ln=cris.nextList();
//				ArrayList<Read> reads=(ln!=null ? ln.list : null);
//
//				if(reads!=null && !reads.isEmpty()){
//					Read r=reads.get(0);
//					assert((r.mate!=null)==cris.paired());
//				}
//
//				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//
//					for(int idx=0; idx<reads.size(); idx++){
//						final Read r1=reads.get(idx);
//						final Read r2=r1.mate;
//
//						final int initialLength1=r1.length();
//						final int initialLength2=(r1.mateLength());
//
//						{
//							readsProcessed++;
//							basesProcessed+=initialLength1;
//							r1.id=core+"_"+r1.numericID+" /1";
//						}
//						if(r2!=null){
//							readsProcessed++;
//							basesProcessed+=initialLength2;
//							r2.id=core+"_"+r1.numericID+" /2";
//						}
//
//
//					}
//
//					final ArrayList<Read> listOut=reads;
//
//					if(ros!=null){ros.add(listOut, ln.id);}
//
//					cris.returnList(ln);
//					ln=cris.nextList();
//					reads=(ln!=null ? ln.list : null);
//				}
//				if(ln!=null){
//					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
//				}
//			}
//
//			errorState|=ReadWrite.closeStreams(cris, ros);
//
//			ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, false, true, false);
//			ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, false, true, false);
//		}
//		readsProcessedA.addAndGet(readsProcessed);
//		basesProcessedA.addAndGet(basesProcessed);
//	}
	
	private void renameAndMerge_MT(){
		FileFormat ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, false, ordered);
		FileFormat ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, false, ordered);
		
		final ConcurrentReadOutputStream ros;
		
		final int buff=4;

		ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, buff, null, false);
		ros.start();
//		System.err.println("Started ros.");
		
		final int threads=Shared.threads();
		ArrayList<MuxThread> list=new ArrayList<MuxThread>(threads);
		for(int i=0; i<threads; i++){
			MuxThread mt=new MuxThread(ros);
			list.add(mt);
			mt.start();
		}
		
//		System.err.println("Started threads.");
		for(MuxThread mt : list){
			while(mt.getState()!=Thread.State.TERMINATED){
				try {
					mt.join();
//					System.err.println("Joined.");
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
//			System.err.println(mt.getState());
		}
//		System.err.println("All threads finished.");
		
		errorState|=ReadWrite.closeStream(ros);
	}
	
	void renameAndMergeOneFile(final String path, final ConcurrentReadOutputStream ros){

		long readsProcessed=0;
		long basesProcessed=0;

		String in1=path, in2=null;

		//Do input file # replacement
		if(path.indexOf('#')>-1 && !new File(path).exists()){
			in2=path.replace("#", "2");
			in1=path.replace("#", "1");
		}

		synchronized(getClass()){

			//Adjust interleaved settings based on number of output files
			if(!setInterleaved){
				assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
				if(in2!=null){ //If there are 2 input streams.
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
					if(!printedInterleavedMessage){
						outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
						printedInterleavedMessage=true;
					}
				}else{ //There is one input stream.
					if(out2!=null){
						FASTQ.FORCE_INTERLEAVED=true;
						FASTQ.TEST_INTERLEAVED=false;
						if(!printedInterleavedMessage){
							outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
							printedInterleavedMessage=true;
						}
					}
				}
			}

			//Ensure input files can be read
			if(!Tools.testInputFiles(false, true, in1, in2)){
				throw new RuntimeException("\nCan't read some input files.\n");  
			}

			//Ensure that no file was specified multiple times
			if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
				throw new RuntimeException("\nSome file names were specified multiple times.\n");
			}
		}

		final FileFormat ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		final FileFormat ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		final ConcurrentReadInputStream cris;
		final String core=ReadWrite.stripToCore(in1);
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}


		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((r.mate!=null)==cris.paired());
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
						r1.id=core+"_"+r1.numericID+" 1:";
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						r2.id=core+"_"+r1.numericID+" 2:";
					}


				}

				final ArrayList<Read> listOut=reads;
				
				if(ros!=null){ros.add(listOut, 0/*nextListNumber.getAndIncrement()*/);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		synchronized(getClass()){
			errorState|=ReadWrite.closeStream(cris);
		}
		readsProcessedA.addAndGet(readsProcessed);
		basesProcessedA.addAndGet(basesProcessed);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class MuxThread extends Thread{
		
		MuxThread(ConcurrentReadOutputStream ros_){
			ros=ros_;
		}
		
		@Override
		public void run(){
			for(int i=nextPathNumber.getAndIncrement(); i<readPaths.size(); i=nextPathNumber.getAndIncrement()){
				renameAndMergeOneFile(readPaths.get(i), ros);
			}
		}
		
		final ConcurrentReadOutputStream ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	protected ArrayList<String> readPaths=new ArrayList<String>();
	
	protected String out1, out2;
	
	protected String extin;
	protected String extout;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected AtomicLong readsProcessedA=new AtomicLong(0);
	/** Number of bases processed */
	protected AtomicLong basesProcessedA=new AtomicLong(0);
	
	protected AtomicLong nextListNumber=new AtomicLong(0);
	
	protected AtomicInteger nextPathNumber=new AtomicInteger(0);

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	private boolean printedInterleavedMessage=false;
}
