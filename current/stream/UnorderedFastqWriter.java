package stream;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Unordered single/interleaved FASTQ writer for multiple producer threads.
 * Each producer formats a complete batch before queue submission, retaining
 * mate adjacency while avoiding ordered-list head-of-line blocking.
 *
 * @author Collei
 * @date September 20, 2026
 */
public final class UnorderedFastqWriter implements Writer {

	public UnorderedFastqWriter(FileFormat ffout_, int producerThreads){
		if(ffout_==null || !ffout_.fastq() || ffout_.ordered()){
			throw new IllegalArgumentException("UnorderedFastqWriter requires unordered FASTQ output: "+ffout_);
		}
		ffout=ffout_;
		fname=ffout.name();
		final int capacity=Tools.max(4, (3*Tools.max(1, producerThreads))/2+4);
		queue=new ArrayBlockingQueue<OutputJob>(capacity);
		outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
	}

	@Override
	public void start(){
		synchronized(stateLock){
			if(started){throw new IllegalStateException("Writer already started: "+fname);}
			started=true;
			outputThread=new Thread(new Runnable(){
				@Override
				public void run(){writeLoop();}
			}, "UnorderedFastqWriter-Output");
			outputThread.start();
		}
	}

	@Override
	public void add(ArrayList<Read> reads, long id){addReads(new ListNum<Read>(reads, id));}

	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		beginAdd();
		try{
			put(format(reads.list));
		}catch(Throwable t){
			fail(t, false);
			throw asRuntime(t);
		}finally{
			endAdd();
		}
	}

	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		final ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			if(sl!=null){reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));}
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}

	private OutputJob format(ArrayList<Read> reads){
		final ByteBuilder bb=builders.get();
		assert(bb.length()==0) : "Thread-local FASTQ builder was not cleared after its prior batch";
		long bases=0;
		int count=0;
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(r1!=null){
				r1.toFastq(bb).nl();
				count++;
				bases+=r1.length();
			}
			if(r2!=null){
				r2.toFastq(bb).nl();
				count++;
				bases+=r2.length();
			}
		}
		final byte[] bytes=bb.toBytes();
		bb.clear();
		return new OutputJob(bytes, count, bases, false, false);
	}

	private void put(OutputJob job){
		try{
			while(!queue.offer(job, 100, TimeUnit.MILLISECONDS)){
				if(errorState){throw new RuntimeException("Writer failed while submitting output: "+fname);}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new RuntimeException("Interrupted while submitting output: "+fname, e);
		}
		if(errorState){throw new RuntimeException("Writer failed while submitting output: "+fname);}
	}

	private void beginAdd(){
		synchronized(stateLock){
			if(!started){throw new IllegalStateException("Writer has not started: "+fname);}
			if(!accepting || errorState){throw new IllegalStateException("Writer is closed: "+fname);}
			activeAdds++;
		}
	}

	private void endAdd(){
		synchronized(stateLock){
			activeAdds--;
			assert(activeAdds>=0) : "FASTQ producer count became negative; beginAdd/endAdd are unbalanced";
			if(activeAdds==0){stateLock.notifyAll();}
		}
	}

	@Override
	public void poison(){
		synchronized(stateLock){
			if(poisoned){return;}
			accepting=false;
			while(activeAdds>0 && !errorState){
				try{stateLock.wait();}
				catch(InterruptedException e){
					Thread.currentThread().interrupt();
					fail(e, false);
					return;
				}
			}
			poisoned=true;
		}
		if(!errorState){put(OutputJob.POISON);}
	}

	@Override
	public boolean waitForFinish(){
		synchronized(stateLock){
			while(!finished){
				try{stateLock.wait();}
				catch(InterruptedException e){
					Thread.currentThread().interrupt();
					fail(e, false);
					break;
				}
			}
		}
		return errorState;
	}

	@Override
	public boolean poisonAndWait(){poison();return waitForFinish();}

	@Override
	public void finishError(){fail(null, true);}

	private void fail(Throwable t, boolean interruptWriter){
		final Thread thread;
		synchronized(stateLock){
			errorState=true;
			accepting=false;
			thread=outputThread;
			stateLock.notifyAll();
		}
		queue.clear();
		queue.offer(OutputJob.ABORT);
		if(interruptWriter && thread!=null && thread!=Thread.currentThread()){thread.interrupt();}
		if(t!=null && thread==Thread.currentThread()){t.printStackTrace();}
	}

	private RuntimeException asRuntime(Throwable t){
		return t instanceof RuntimeException ? (RuntimeException)t : new RuntimeException(t);
	}

	private void writeLoop(){
		boolean normal=false;
		try{
			while(true){
				final OutputJob job=queue.take();
				if(job.abort){break;}
				if(job.poison){normal=true;break;}
				outstream.write(job.bytes);
				readsWritten+=job.reads;
				basesWritten+=job.bases;
			}
		}catch(Throwable t){
			if(!(t instanceof InterruptedException && errorState)){t.printStackTrace();}
			fail(t, false);
		}finally{
			final boolean closeError=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			synchronized(stateLock){
				errorState|=closeError || !normal;
				finished=true;
				stateLock.notifyAll();
			}
		}
	}

	@Override
	public long readsWritten(){return readsWritten;}
	@Override
	public long basesWritten(){return basesWritten;}
	@Override
	public String fname(){return fname;}
	@Override
	public boolean errorState(){return errorState;}
	@Override
	public boolean finishedSuccessfully(){return finished && !errorState;}

	private static final class OutputJob {
		OutputJob(byte[] bytes_, int reads_, long bases_, boolean poison_, boolean abort_){
			bytes=bytes_;reads=reads_;bases=bases_;poison=poison_;abort=abort_;
		}
		final byte[] bytes;
		final int reads;
		final long bases;
		final boolean poison;
		final boolean abort;
		static final OutputJob POISON=new OutputJob(null, 0, 0, true, false);
		static final OutputJob ABORT=new OutputJob(null, 0, 0, false, true);
	}

	private final ThreadLocal<ByteBuilder> builders=new ThreadLocal<ByteBuilder>(){
		@Override
		protected ByteBuilder initialValue(){return new ByteBuilder();}
	};
	private final Object stateLock=new Object();
	private final ArrayBlockingQueue<OutputJob> queue;
	private final FileFormat ffout;
	private final String fname;
	private final OutputStream outstream;
	private volatile boolean errorState=false;
	private volatile boolean finished=false;
	private boolean accepting=true;
	private boolean started=false;
	private boolean poisoned=false;
	private int activeAdds=0;
	private Thread outputThread;
	private volatile long readsWritten=0;
	private volatile long basesWritten=0;
}
