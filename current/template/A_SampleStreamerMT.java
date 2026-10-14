package template;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;

/**
 * Multithreaded sequence-processing template using shared Streamer and Writer APIs.
 * Each worker consumes numbered Read batches, processes roots and linked mates,
 * and submits the original batch ID even when filtering removes every read.
 * Copy this template and customize the private read-pair hook; the supplied hook
 * retains all reads. Configuration changes process-wide I/O settings.
 * <p>Created May 1, 2026.
 * @author Brian Bushnell, Eru, Shinobu
 */
public class A_SampleStreamerMT implements Accumulator<A_SampleStreamerMT.ProcessThread>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the example and closes redirected diagnostic output on normal return.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		A_SampleStreamerMT x=new A_SampleStreamerMT(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/** Parses options, configures global I/O settings and creates file descriptors.
	 * Requires a primary input; secondary input/output paths are optional.
	 * @param args Command-line arguments, including standard Parser options */
	public A_SampleStreamerMT(String[] args){

		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=parser.overwrite;
			append=parser.append;
			setInterleaved=parser.setInterleaved;
			workers=parser.workers();
			threadsIn=parser.threadsIn;
			threadsOut=parser.threadsOut;

			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
			out2=parser.out2;
		}

		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		//WriterFactory requires ordered descriptors for separate mate outputs.
		final boolean orderedOut=ordered || out2!=null;
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, orderedOut);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, orderedOut);
		//Retain input alignments when writing SAM/BAM without rebuilding a genome index.
		final boolean samIn=ffin1.samOrBam();
		SamLine.SET_FROM_OK=samIn;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=samIn && ffout1!=null && ffout1.samOrBam();
	}

	/** Parses local options and forwards recognized standard options to Parser.
	 * Positional input/output names are also accepted in the first two positions.
	 * @param args Preparsed arguments
	 * @return Parser containing the standard settings copied by the constructor */
	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			//Preserve the complete value after the first equals sign.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null && Tools.looksLikeOutputSequenceStream(arg)){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts shared I/O and workers, accumulates results and finalizes normal output.
	 * Workers always request Read objects, including for SAM input/output. Changes
	 * Read.VALIDATE_IN_CONSTRUCTOR globally; it is not restored here. Reader and
	 * writer completion errors join the accumulated worker status before reporting.
	 * @param t Timer started by main */
	private void process(Timer t){
		final boolean saveHeader=(ffin1!=null && ffin1.samOrBam() && ffout1!=null && ffout1.samOrBam());

		final int threads=Tools.max(1, workers>0 ? workers : Shared.threads());

		Read.VALIDATE_IN_CONSTRUCTOR=(threads<2 && threadsIn<2);

		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads,
			saveHeader, true, threadsIn);//Workers consume Reads even for SAM output or no output.
		st.start();

		Writer fw=WriterFactory.makeWriter(ffout1, ffout2, threadsOut, null, saveHeader);
		if(fw!=null){fw.start();}

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(st, fw, i));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;

		if(fw!=null){
			errorState|=fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}

		errorState|=ReadWrite.closeStream(st);

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		if(readsOut>0){
			outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds one completed worker under ThreadWaiter coordination; collect each worker once.
	 * Read totals count individual reads, including linked mates; base totals sum
	 * both mates. A worker that did not reach normal completion sets errorState.
	 * @param pt Worker whose final counters are ready to collect */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
		}
	}

	/** Reports accumulated worker status and any collected I/O completion errors. */
	@Override
	public final boolean success(){return !errorState;}

	/** Supplies the lock used by ThreadWaiter to coordinate accumulation. */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	/** Coordination lock; does not serialize workers' processing of independent batches. */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns each fetched batch until submission and keeps worker-local statistics.
	 * The reader and optional writer are shared across workers. */
	static class ProcessThread extends Thread{

		/** Binds shared I/O; lifecycle is managed by the enclosing driver.
		 * @param st_ Started reader
		 * @param fw_ Started writer, or null when output is disabled
		 * @param tid_ Worker index available for customization */
		ProcessThread(Streamer st_, Writer fw_, int tid_){
			st=st_;
			fw=fw_;
			tid=tid_;
		}

		/** Consumes until null, retaining empty batches and dense output IDs.
		 * Input counts precede the hook; retained counts follow it. Rejected roots
		 * become null list entries so writers omit them and their linked mates.
		 * Submitted lists and reads are not mutated afterward. */
		@Override
		public void run(){
			ListNum<Read> ln=st.nextList();
			while(ln!=null){//Empty batches retain their IDs and do not terminate input.
				for(int idx=0; idx<ln.size(); idx++){
					final Read r=ln.get(idx);
					final Read r2=r.mate;
					if(!r.validated()){r.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}

					readsProcessedT+=r.pairCount();
					basesProcessedT+=r.pairLength();

					boolean keep=processReadPair(r, r2);

					if(keep){
						readsOutT+=r.pairCount();
						basesOutT+=r.pairLength();
					}else{
						r.setDiscarded(true);
						ln.list.set(idx, null);//Writers skip null roots, not the discarded flag; retain the batch ID.
					}
				}
				if(fw!=null){fw.addReads(ln);}
				ln=st.nextList();
			}
			success=true;
		}

		/** Customization hook; the supplied implementation retains every root.
		 * @param r1 Nonnull input root, owned by this worker until submission
		 * @param r2 Linked mate, or null
		 * @return Whether the root and its linked mate should be retained */
		private boolean processReadPair(Read r1, Read r2){
			return true;
		}

		/** Individual reads encountered before customization, including mates. */
		long readsProcessedT=0;
		/** Bases encountered before customization, including mates. */
		long basesProcessedT=0;
		/** Individual reads retained after customization, including mates. */
		long readsOutT=0;
		/** Bases retained after customization, including mates. */
		long basesOutT=0;
		/** True only after the reader returns terminal null and run completes normally. */
		boolean success=false;

		/** Shared reader, consumed concurrently through nextList. */
		private final Streamer st;
		/** Shared writer; null disables output without disabling processing. */
		private final Writer fw;
		/** Worker index, currently unused by the default hook. */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input path. */
	private String in1=null;
	/** Optional second input path. */
	private String in2=null;
	/** Optional primary output path. */
	private String out1=null;
	/** Optional second output path. */
	private String out2=null;

	/** Primary input descriptor. */
	private FileFormat ffin1;
	/** Optional secondary input descriptor. */
	private FileFormat ffin2;
	/** Optional primary output descriptor; separate mate outputs remain ordered. */
	private FileFormat ffout1;
	/** Optional secondary output descriptor, ordered whenever present. */
	private FileFormat ffout2;

	/** Requested worker count; nonpositive uses Shared.threads, with a minimum of one. */
	private int workers=-1;
	/** Reader thread setting forwarded to StreamerFactory. */
	private int threadsIn=-1;
	/** Writer thread setting forwarded to WriterFactory. */
	private int threadsOut=-1;
	/** Requested input fragment limit; negative means unlimited. */
	private long maxReads=-1;

	/** Accumulated individual-read input count, including linked mates. */
	private long readsProcessed=0;
	/** Accumulated input bases before customization, including mates. */
	private long basesProcessed=0;
	/** Retained individual reads; replaced by writer counters when output is enabled. */
	private long readsOut=0;
	/** Retained bases; replaced by writer counters when output is enabled. */
	private long basesOut=0;

	/** Whether output descriptors permit replacement of existing files. */
	private boolean overwrite=false;
	/** Whether output descriptors request append mode. */
	private boolean append=false;
	/** Input and single-output order request; separate mate outputs remain ordered. */
	private boolean ordered=true;
	/** Parser marker retained here but not subsequently consulted. */
	private boolean setInterleaved=false;
	/** Accumulated worker and I/O completion error flag. */
	private boolean errorState=false;
	/** Diagnostic output, possibly redirected by PreParser. */
	private PrintStream outstream=System.err;

	/** Parsed verbosity setting; currently unused by this example. */
	public static boolean verbose=false;

}
