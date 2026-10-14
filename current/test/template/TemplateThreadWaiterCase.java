package test.template;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import template.Accumulator;
import template.ThreadWaiter;

/** Direct checks for ThreadWaiter accumulator lock release on exceptional exits.
 * Baseline mode observes the existing leaked locks; fixed mode requires both
 * read and write locks to be released before exceptions propagate.
 * @author Jean
 */
public final class TemplateThreadWaiterCase{

	/** Runs finite in-process ThreadWaiter exceptional-lock cases. */
	public static void main(String[] args) throws Exception{
		final String expect=parseExpect(args);

		startAndWaitStartThrow(expect);
		waitOnlyIteratorThrow(expect);
		startAndWaitAccumulateThrow(expect);
		waitOnlyAccumulateThrow(expect);
		startAndWaitSuccessThrow(expect);
		startAndWaitNormal(expect);
		waitOnlyNormal(expect);
		startAndWaitNormalNullRwlock(expect);
		waitOnlyNormalNullRwlock(expect);
	}

	/** Parses the advertised single expect=baseline|fixed argument. */
	private static String parseExpect(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		return expect;
	}

	/** Exercises startAndWait when the second start throws under the read lock. */
	private static void startAndWaitStartThrow(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		final Thread thread=new NoopThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			ThreadWaiter.startAndWait(Arrays.asList(thread, thread), acc);
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown instanceof IllegalThreadStateException,
				"Expected IllegalThreadStateException, got "+thrown);
		checkLocks("start-throw", expect, readHolds, writeHolds, LockPhase.READ);
	}

	/** Exercises the wait-only overload when iteration throws under the read lock. */
	private static void waitOnlyIteratorThrow(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		final Thread thread=terminatedThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			ThreadWaiter.waitForThreadsToFinish(new ThrowingIterable<Thread>(thread), acc);
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown instanceof MarkerException, "Expected MarkerException, got "+thrown);
		checkLocks("wait-iterator-throw", expect, readHolds, writeHolds, LockPhase.READ);
	}

	/** Exercises startAndWait when accumulate throws under the write lock. */
	private static void startAndWaitAccumulateThrow(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		acc.throwAccumulate=true;
		final Thread thread=new NoopThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			ThreadWaiter.startAndWait(Collections.singletonList(thread), acc);
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown instanceof MarkerException, "Expected MarkerException, got "+thrown);
		checkLocks("start-accumulate-throw", expect, readHolds, writeHolds,
				LockPhase.WRITE);
	}

	/** Exercises the wait-only overload when accumulate throws under the write lock. */
	private static void waitOnlyAccumulateThrow(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		acc.throwAccumulate=true;
		final Thread thread=terminatedThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			ThreadWaiter.waitForThreadsToFinish(Collections.singletonList(thread), acc);
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown instanceof MarkerException, "Expected MarkerException, got "+thrown);
		checkLocks("wait-accumulate-throw", expect, readHolds, writeHolds,
				LockPhase.WRITE);
	}

	/** Exercises acc.success throwing after all accumulate calls complete. */
	private static void startAndWaitSuccessThrow(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		acc.throwSuccess=true;
		final Thread thread=new NoopThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			ThreadWaiter.startAndWait(Collections.singletonList(thread), acc);
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown instanceof MarkerException, "Expected MarkerException, got "+thrown);
		checkLocks("start-success-throw", expect, readHolds, writeHolds,
				LockPhase.WRITE);
	}

	/** Confirms normal start-and-wait accumulation releases both locks. */
	private static void startAndWaitNormal(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		final Thread thread=new NoopThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			final boolean success=ThreadWaiter.startAndWait(
					Collections.singletonList(thread), acc);
			check(success, "Expected normal startAndWait success");
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown==null, "Unexpected normal startAndWait exception: "+thrown);
		check(acc.accumulated==1, "Expected one accumulated thread, got "+acc.accumulated);
		checkLocks("start-normal", expect, readHolds, writeHolds, LockPhase.NONE);
	}

	/** Confirms normal wait-only accumulation releases both locks. */
	private static void waitOnlyNormal(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator();
		final Thread thread=terminatedThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			final boolean success=ThreadWaiter.waitForThreadsToFinish(
					Collections.singletonList(thread), acc);
			check(success, "Expected normal waitForThreadsToFinish success");
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown==null, "Unexpected normal waitForThreadsToFinish exception: "+thrown);
		check(acc.accumulated==1, "Expected one accumulated thread, got "+acc.accumulated);
		checkLocks("wait-normal", expect, readHolds, writeHolds, LockPhase.NONE);
	}

	/** Confirms normal start-and-wait behavior when the accumulator has no lock. */
	private static void startAndWaitNormalNullRwlock(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator(false);
		final Thread thread=new NoopThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			final boolean success=ThreadWaiter.startAndWait(
					Collections.singletonList(thread), acc);
			check(success, "Expected normal null-rwlock startAndWait success");
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown==null, "Unexpected null-rwlock startAndWait exception: "+thrown);
		check(acc.accumulated==1, "Expected one accumulated thread, got "+acc.accumulated);
		checkLocks("start-normal-null-rwlock", expect, readHolds, writeHolds,
				LockPhase.NONE);
	}

	/** Confirms normal wait-only behavior when the accumulator has no lock. */
	private static void waitOnlyNormalNullRwlock(final String expect) throws Exception{
		final LockingAccumulator acc=new LockingAccumulator(false);
		final Thread thread=terminatedThread();
		Throwable thrown=null;
		int readHolds=0;
		int writeHolds=0;
		try{
			final boolean success=ThreadWaiter.waitForThreadsToFinish(
					Collections.singletonList(thread), acc);
			check(success, "Expected normal null-rwlock waitForThreadsToFinish success");
		}catch(Throwable t){
			thrown=t;
		}finally{
			try{
				join(thread);
			}finally{
				readHolds=acc.lock.getReadHoldCount();
				writeHolds=acc.lock.getWriteHoldCount();
				acc.drainLocks();
			}
		}
		check(thrown==null,
				"Unexpected null-rwlock waitForThreadsToFinish exception: "+thrown);
		check(acc.accumulated==1, "Expected one accumulated thread, got "+acc.accumulated);
		checkLocks("wait-normal-null-rwlock", expect, readHolds, writeHolds,
				LockPhase.NONE);
	}

	/** Starts a no-op thread and waits for it to terminate. */
	private static Thread terminatedThread() throws Exception{
		final Thread thread=new NoopThread();
		thread.start();
		join(thread);
		return thread;
	}

	/** Joins a thread if it was already started. */
	private static void join(final Thread thread) throws Exception{
		if(thread.getState()==Thread.State.NEW){return;}
		thread.join(30000);
		check(thread.getState()==Thread.State.TERMINATED,
				"Thread did not terminate: "+thread.getState());
	}

	/** Checks leaked lock counts before fixture cleanup. */
	private static void checkLocks(final String name, final String expect,
			final int readHolds, final int writeHolds, final LockPhase phase){
		final int expectedRead=("baseline".equals(expect) && phase==LockPhase.READ) ? 1 : 0;
		final int expectedWrite=("baseline".equals(expect) && phase==LockPhase.WRITE) ? 1 : 0;
		check(readHolds==expectedRead, name+" read holds "+readHolds
				+" != "+expectedRead);
		check(writeHolds==expectedWrite, name+" write holds "+writeHolds
				+" != "+expectedWrite);
		System.out.println("PASS\t"+expect+"\t"+name+"\tread="+readHolds
				+"\twrite="+writeHolds);
	}

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Optional accumulator lock phase expected to leak in the baseline. */
	private enum LockPhase{
		NONE, READ, WRITE
	}

	/** No-op worker thread for local ThreadWaiter cases. */
	private static final class NoopThread extends Thread{
		@Override
		public void run(){
			//No-op.
		}
	}

	/** Marker exception thrown from owned callbacks and iterators. */
	private static final class MarkerException extends RuntimeException{
		private static final long serialVersionUID=1L;
	}

	/** Accumulator with a concrete lock and optional throwing callbacks. */
	private static final class LockingAccumulator implements Accumulator<Thread>{
		LockingAccumulator(){
			this(true);
		}

		LockingAccumulator(boolean returnLock){
			this.returnLock=returnLock;
		}

		@Override
		public void accumulate(Thread t){
			if(throwAccumulate){throw new MarkerException();}
			accumulated++;
		}

		@Override
		public ReadWriteLock rwlock(){return returnLock ? lock : null;}

		@Override
		public boolean success(){
			if(throwSuccess){throw new MarkerException();}
			return true;
		}

		/** Releases any same-thread locks left by the baseline defect. */
		void drainLocks(){
			while(lock.getReadHoldCount()>0){lock.readLock().unlock();}
			while(lock.getWriteHoldCount()>0){lock.writeLock().unlock();}
		}

		final ReentrantReadWriteLock lock=new ReentrantReadWriteLock();
		final boolean returnLock;
		boolean throwAccumulate=false;
		boolean throwSuccess=false;
		int accumulated=0;
	}

	/** Iterable that throws while ThreadWaiter is still in its read-locked phase. */
	private static final class ThrowingIterable<T> implements Iterable<T>{
		ThrowingIterable(final T value){
			this.value=value;
		}

		@Override
		public Iterator<T> iterator(){
			return new Iterator<T>(){
				@Override
				public boolean hasNext(){
					if(returned){throw new MarkerException();}
					return true;
				}

				@Override
				public T next(){
					returned=true;
					return value;
				}

				boolean returned=false;
			};
		}

		private final T value;
	}
}
