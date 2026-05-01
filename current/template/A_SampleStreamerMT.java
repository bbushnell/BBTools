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
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;

/**
 * This class does nothing.
 * It is designed to be easily modified into a program
 * that processes reads in multiple threads using a Streamer.
 *
 * @author Brian Bushnell, Eru
 * @date May 1, 2026
 */
public class A_SampleStreamerMT implements Accumulator<A_SampleStreamerMT.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) {
		Timer t=new Timer();
		A_SampleStreamerMT x=new A_SampleStreamerMT(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

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
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);
	}

	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
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
			}else if(i==1 && parser.in1!=null && Tools.looksLikeOutputSequenceStream(arg)) {
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

	private void process(Timer t) {
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean saveHeader=(ffin1!=null && ffin1.samOrBam() && ffout1!=null && ffout1.samOrBam());

		final int threads=Tools.max(1, workers>0 ? workers : Shared.threads());

		Read.VALIDATE_IN_CONSTRUCTOR=(threads<2 && threadsIn<2);

		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads,
			saveHeader, outputReads, threadsIn);
		st.start();

		Writer fw=WriterFactory.makeWriter(ffout1, ffout2, threadsOut, null, saveHeader);
		if(fw!=null){fw.start();}

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(st, fw, i));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;

		if(fw!=null){
			fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}

		ReadWrite.closeStream(st);

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

	@Override
	public final void accumulate(ProcessThread pt) {
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
		}
	}

	@Override
	public final boolean success() {return !errorState;}

	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	static class ProcessThread extends Thread {

		ProcessThread(Streamer st_, Writer fw_, int tid_) {
			st=st_;
			fw=fw_;
			tid=tid_;
		}

		@Override
		public void run() {
			ListNum<Read> ln=st.nextList();
			while(ln!=null && ln.size()>0) {
				for(Read r : ln) {
					final Read r2=r.mate;

					readsProcessedT+=r.pairCount();
					basesProcessedT+=r.pairLength();

					boolean keep=processReadPair(r, r2);

					if(keep){
						readsOutT+=r.pairCount();
						basesOutT+=r.pairLength();
					}else{
						r.setDiscarded(true);
					}
				}
				if(fw!=null){fw.addReads(ln);}
				ln=st.nextList();
			}
			success=true;
		}

		private boolean processReadPair(Read r1, Read r2) {
			return true;
		}

		long readsProcessedT=0;
		long basesProcessedT=0;
		long readsOutT=0;
		long basesOutT=0;
		boolean success=false;

		private final Streamer st;
		private final Writer fw;
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;
	private String out1=null;
	private String out2=null;

	private FileFormat ffin1;
	private FileFormat ffin2;
	private FileFormat ffout1;
	private FileFormat ffout2;

	private int workers=-1;
	private int threadsIn=-1;
	private int threadsOut=-1;
	private long maxReads=-1;

	private long readsProcessed=0;
	private long basesProcessed=0;
	private long readsOut=0;
	private long basesOut=0;

	private boolean overwrite=false;
	private boolean append=false;
	private boolean ordered=true;
	private boolean setInterleaved=false;
	private boolean errorState=false;
	private PrintStream outstream=System.err;

	public static boolean verbose=false;

}
