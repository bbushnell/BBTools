package test.clump2;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLongArray;

import stream.Read;
import stream.SamLine;
import stream.Writer;
import structures.ListNum;

/**
 * Focused tests for KmerSort3's bounded parallel zero-threaded destination
 * submitter. Uses fake Writers through reflection so package-private helper
 * behavior is tested without touching real file outputs.
 * @author Sucrose
 */
public final class KmerSort3ParallelDestinationWriterTest{

	/** Runs all helper-contract cases. */
	public static void main(String[] args) throws Throwable{
		check(args.length==0, "Expected no arguments");
		System.out.println("case\tok");
		checkEachDestinationOnce();
		checkMultipleRoundsPreserveDestinationOrder();
		checkBoundedWorkerReuse();
		checkConcurrentDestinationAdds();
		checkEmptyRoundsSkipIds();
		checkWorkerInterruptedAfterPublicationCompletesRound();
		checkWorkerFailureWakesCoordinator();
	}

	/** One nonempty list per destination should submit once with ID zero. */
	private static void checkEachDestinationOnce() throws Throwable{
		final FakeWriter[] writers=fakes(4);
		final Harness harness=new Harness(writers, 4);
		try{
			final ArrayList<Read>[] round=fullRound(writers.length, "once");
			harness.submitRound(round);
			for(int i=0; i<writers.length; i++){
				checkIds(writers[i], 0);
				check(writers[i].call(0).list==round[i], "Destination "+i+" did not receive its list");
			}
			checkIds(harness.ids, 1, 1, 1, 1);
		}finally{harness.close();}
		System.out.println("each-destination-once\ttrue");
	}

	/** Dense IDs must advance independently and only for nonempty submissions. */
	private static void checkMultipleRoundsPreserveDestinationOrder() throws Throwable{
		final FakeWriter[] writers=fakes(3);
		final Harness harness=new Harness(writers, 2);
		try{
			ArrayList<Read>[] round=emptyRound(writers.length);
			round[0]=reads("first-0", 1);
			round[2]=reads("first-2", 1);
			harness.submitRound(round);
			round=emptyRound(writers.length);
			round[0]=reads("second-0", 1);
			round[1]=reads("second-1", 1);
			round[2]=reads("second-2", 1);
			harness.submitRound(round);
			checkIds(writers[0], 0, 1);
			checkIds(writers[1], 0);
			checkIds(writers[2], 0, 1);
			checkIds(harness.ids, 2, 1, 2);
		}finally{harness.close();}
		System.out.println("multi-round-order\ttrue");
	}

	/** More destinations than workers should reuse one bounded worker pool. */
	private static void checkBoundedWorkerReuse() throws Throwable{
		final FakeWriter[] writers=fakes(9);
		final CountDownLatch entered=new CountDownLatch(3);
		final CountDownLatch release=new CountDownLatch(1);
		for(FakeWriter writer : writers){writer.gate(entered, release);}
		final Harness harness=new Harness(writers, 3);
		try{
			final Submitter submitter=new Submitter(harness, fullRound(writers.length, "reuse"));
			submitter.start();
			try{await(entered, "Expected three concurrent helper workers");}
			finally{release.countDown();}
			submitter.joinSuccessfully();
			final HashSet<Thread> firstRound=threads(writers);
			check(totalCalls(writers)==writers.length, "Expected one call per destination");
			check(firstRound.size()==3, "Expected exactly three reused workers");

			final CountDownLatch enteredAgain=new CountDownLatch(3);
			final CountDownLatch releaseAgain=new CountDownLatch(1);
			for(FakeWriter writer : writers){writer.gate(enteredAgain, releaseAgain);}
			final Submitter second=new Submitter(harness, fullRound(writers.length, "reuse-again"));
			second.start();
			try{await(enteredAgain, "Expected same three helper workers in the second round");}
			finally{releaseAgain.countDown();}
			second.joinSuccessfully();
			check(totalCalls(writers)==2*writers.length, "Expected two calls per destination");
			final HashSet<Thread> secondRound=threads(writers, 1);
			check(secondRound.equals(firstRound),
					"Expected same helper Thread objects across rounds; first="+firstRound+", second="+secondRound);
		}finally{harness.close();}
		System.out.println("bounded-worker-reuse\ttrue");
	}

	/** Two slow destinations must be in Writer.add at the same time. */
	private static void checkConcurrentDestinationAdds() throws Throwable{
		final FakeWriter[] writers=fakes(2);
		final CountDownLatch entered=new CountDownLatch(2);
		final CountDownLatch release=new CountDownLatch(1);
		for(FakeWriter writer : writers){writer.gate(entered, release);}
		final Harness harness=new Harness(writers, 2);
		try{
			final Submitter submitter=new Submitter(harness, fullRound(writers.length, "parallel"));
			submitter.start();
			try{await(entered, "Expected both destination writes to overlap");}
			finally{release.countDown();}
			submitter.joinSuccessfully();
		}finally{harness.close();}
		System.out.println("concurrent-destinations\ttrue");
	}

	/** Null and empty destinations should not call add or consume output IDs. */
	private static void checkEmptyRoundsSkipIds() throws Throwable{
		final FakeWriter[] writers=fakes(4);
		final Harness harness=new Harness(writers, 4);
		try{
			final ArrayList<Read>[] round=emptyRound(writers.length);
			round[1]=new ArrayList<Read>();
			round[3]=new ArrayList<Read>();
			harness.submitRound(round);
			check(totalCalls(writers)==0, "Empty lists should not be submitted");
			checkIds(harness.ids, 0, 0, 0, 0);
		}finally{harness.close();}
		System.out.println("empty-round\ttrue");
	}

	/** Interrupting a worker after publication must not strand the round count. */
	private static void checkWorkerInterruptedAfterPublicationCompletesRound() throws Throwable{
		final FakeWriter[] writers=fakes(4);
		final Harness harness=new Harness(writers, 2);
		try{
			harness.awaitWorkersWaiting();
			final Submitter submitter=new Submitter(harness, fullRound(writers.length, "interrupt"), true);
			submitter.start();
			submitter.joinFinished();
			check(totalCalls(writers)==writers.length,
					"Interrupted worker stranded destinations before finishing the round");
			checkIds(harness.ids, 1, 1, 1, 1);
		}finally{harness.close();}
		System.out.println("post-publish-interrupt\ttrue");
	}

	/** A destination exception must be rethrown on submitRound and all workers must join. */
	private static void checkWorkerFailureWakesCoordinator() throws Throwable{
		final FakeWriter[] writers=fakes(4);
		writers[0].fail(new IllegalStateException("synthetic writer failure"));
		final Harness harness=new Harness(writers, 2);
		try{
			final Submitter submitter=new Submitter(harness, fullRound(writers.length, "fail"));
			submitter.start();
			submitter.joinWithFailure(IllegalStateException.class);
		}finally{harness.close();}
		System.out.println("worker-failure\ttrue");
	}

	/** Creates one nonempty list per destination. */
	private static ArrayList<Read>[] fullRound(final int destinations, final String prefix){
		final ArrayList<Read>[] round=emptyRound(destinations);
		for(int i=0; i<round.length; i++){round[i]=reads(prefix+"-"+i, 1);}
		return round;
	}

	/** Creates a null-filled round array. */
	@SuppressWarnings("unchecked")
	private static ArrayList<Read>[] emptyRound(final int destinations){
		return new ArrayList[destinations];
	}

	/** Creates synthetic reads whose payloads are irrelevant to the writer helper. */
	private static ArrayList<Read> reads(final String prefix, final int count){
		final ArrayList<Read> list=new ArrayList<Read>(count);
		for(int i=0; i<count; i++){
			list.add(new Read(BASES, QUALS, prefix+"-"+i, i));
		}
		return list;
	}

	/** Creates fake destinations. */
	private static FakeWriter[] fakes(final int count){
		final FakeWriter[] writers=new FakeWriter[count];
		for(int i=0; i<count; i++){writers[i]=new FakeWriter("fake-"+i);}
		return writers;
	}

	/** Checks the calls made to a fake writer. */
	private static void checkIds(final FakeWriter writer, final long... expected){
		final ArrayList<Call> calls=writer.calls();
		check(calls.size()==expected.length, writer.name+" calls: "+calls.size());
		for(int i=0; i<expected.length; i++){
			check(calls.get(i).id==expected[i],
					writer.name+" call "+i+" expected ID "+expected[i]+", got "+calls.get(i).id);
		}
	}

	/** Checks the post-round ID counters. */
	private static void checkIds(final AtomicLongArray ids, final long... expected){
		check(ids.length()==expected.length, "ID counter length changed");
		for(int i=0; i<expected.length; i++){
			check(ids.get(i)==expected[i],
					"Destination "+i+" expected next ID "+expected[i]+", got "+ids.get(i));
		}
	}

	/** Counts all fake-writer calls. */
	private static int totalCalls(final FakeWriter[] writers){
		int calls=0;
		for(FakeWriter writer : writers){calls+=writer.calls().size();}
		return calls;
	}

	/** Returns distinct helper threads observed by fake writers. */
	private static HashSet<Thread> threads(final FakeWriter[] writers){
		final HashSet<Thread> threads=new HashSet<Thread>();
		for(FakeWriter writer : writers){
			for(Call call : writer.calls()){threads.add(call.thread);}
		}
		return threads;
	}

	/** Returns distinct helper threads observed at one per-writer call index. */
	private static HashSet<Thread> threads(final FakeWriter[] writers, final int callIndex){
		final HashSet<Thread> threads=new HashSet<Thread>();
		for(FakeWriter writer : writers){threads.add(writer.call(callIndex).thread);}
		return threads;
	}

	/** Waits for a latch with a hard bound so failures cannot hang the test. */
	private static void await(final CountDownLatch latch, final String message){
		try{
			check(latch.await(5, TimeUnit.SECONDS), message);
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new RuntimeException(e);
		}
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reflects over the package-private helper. */
	private static final class Harness implements AutoCloseable{

		Harness(final Writer[] destinations, final int maxWorkers) throws Exception{
			ids=new AtomicLongArray(destinations.length);
			final Class<?> type=Class.forName("clump2.ParallelDestinationWriter");
			final Constructor<?> constructor=type.getDeclaredConstructor(Writer[].class,
					AtomicLongArray.class, Integer.TYPE);
			constructor.setAccessible(true);
			helper=constructor.newInstance(destinations, ids, maxWorkers);
			submitRound=type.getDeclaredMethod("submitRound", ArrayList[].class);
			close=type.getDeclaredMethod("close");
			lock=type.getDeclaredField("lock");
			workers=type.getDeclaredField("workers");
			submitRound.setAccessible(true);
			close.setAccessible(true);
			lock.setAccessible(true);
			workers.setAccessible(true);
		}

		void submitRound(final ArrayList<Read>[] round) throws Throwable{
			try{
				submitRound.invoke(helper, new Object[] {round});
			}catch(InvocationTargetException e){
				throw e.getCause();
			}
		}

		void submitRoundAfterInterruptingWaiters(final ArrayList<Read>[] round) throws Throwable{
			synchronized(helperLock()){
				for(Thread thread : workerThreads()){thread.interrupt();}
				submitRound(round);
			}
		}

		@Override
		public void close() throws Exception{
			final Throwable[] failure=new Throwable[1];
			final Thread closer=new Thread("ParallelDestinationWriter-close"){
				@Override
				public void run(){
					try{
						invoke(close);
						checkWorkersTerminated();
					}catch(Throwable t){
						failure[0]=t;
					}
				}
			};
			closer.setDaemon(true);
			closer.start();
			closer.join(5000);
			if(closer.isAlive()){
				System.err.println("Timed out closing ParallelDestinationWriter; dumping all threads");
				for(Map.Entry<Thread, StackTraceElement[]> entry : Thread.getAllStackTraces().entrySet()){
					System.err.println(entry.getKey());
					for(StackTraceElement element : entry.getValue()){System.err.println("\tat "+element);}
				}
				Runtime.getRuntime().halt(1);
			}
			rethrow(failure[0]);
		}

		private void invoke(final Method method, final Object... args) throws Exception{
			try{
				method.invoke(helper, args);
			}catch(InvocationTargetException e){
				final Throwable cause=e.getCause();
				if(cause instanceof Error){throw (Error)cause;}
				if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
				if(cause instanceof Exception){throw (Exception)cause;}
				throw new RuntimeException(cause);
			}
		}

		private void checkWorkersTerminated() throws Exception{
			for(Thread thread : workerThreads()){
				check(!thread.isAlive(), "Helper worker still alive: "+thread.getName());
			}
		}

		private void awaitWorkersWaiting() throws Exception{
			final Thread[] threads=workerThreads();
			for(long limit=System.nanoTime()+TimeUnit.SECONDS.toNanos(5); System.nanoTime()<limit; ){
				int waiting=0;
				for(Thread thread : threads){
					if(thread.getState()==Thread.State.WAITING){waiting++;}
				}
				if(waiting==threads.length){return;}
				Thread.sleep(10);
			}
			throw new AssertionError("Timed out waiting for idle helper workers");
		}

		private void rethrow(final Throwable cause) throws Exception{
			if(cause==null){return;}
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Exception){throw (Exception)cause;}
			throw new RuntimeException(cause);
		}

		private Object helperLock() throws IllegalAccessException{
			return lock.get(helper);
		}

		private Thread[] workerThreads() throws IllegalAccessException{
			return (Thread[])workers.get(helper);
		}

		final AtomicLongArray ids;
		final Object helper;
		final Method submitRound, close;
		final Field lock, workers;
	}

	/** Runs submitRound on a separate coordinator so concurrency cases can inspect it. */
	private static final class Submitter extends Thread{

		Submitter(final Harness harness_, final ArrayList<Read>[] round_){
			this(harness_, round_, false);
		}

		Submitter(final Harness harness_, final ArrayList<Read>[] round_,
				final boolean interruptWaiters_){
			harness=harness_;
			round=round_;
			interruptWaiters=interruptWaiters_;
		}

		@Override
		public void run(){
			try{
				if(interruptWaiters){
					harness.submitRoundAfterInterruptingWaiters(round);
				}else{
					harness.submitRound(round);
				}
			}catch(Throwable t){
				failure=t;
			}
		}

		void joinSuccessfully() throws Exception{
			joinBounded();
			if(failure!=null){throw new AssertionError("submitRound failed", failure);}
		}

		void joinWithFailure(final Class<? extends Throwable> expected) throws Exception{
			joinBounded();
			check(failure!=null, "submitRound should have failed");
			check(expected.isInstance(failure),
					"Expected "+expected.getName()+", got "+failure.getClass().getName());
		}

		void joinFinished() throws Exception{
			joinBounded();
		}

		private void joinBounded() throws Exception{
			join(5000);
			check(!isAlive(), "submitRound did not finish");
		}

		private final Harness harness;
		private final ArrayList<Read>[] round;
		private final boolean interruptWaiters;
		private volatile Throwable failure;
	}

	/** Minimal Writer fake that records add calls and can block or fail inside add. */
	private static final class FakeWriter implements Writer{

		FakeWriter(final String name_){
			name=name_;
		}

		void gate(final CountDownLatch entered_, final CountDownLatch release_){
			entered=entered_;
			release=release_;
		}

		void fail(final RuntimeException failure_){
			failure=failure_;
		}

		@Override
		public void add(final ArrayList<Read> reads, final long id){
			synchronized(this){calls.add(new Call(id, reads, Thread.currentThread()));}
			if(entered!=null){entered.countDown();}
			if(release!=null){await(release, "Timed out in "+name);}
			if(failure!=null){throw failure;}
		}

		synchronized Call call(final int index){
			return calls.get(index);
		}

		synchronized ArrayList<Call> calls(){
			return new ArrayList<Call>(calls);
		}

		@Override
		public void start(){}

		@Override
		public long readsWritten(){return 0;}

		@Override
		public long basesWritten(){return 0;}

		@Override
		public void addReads(final ListNum<Read> reads){
			add(reads.list, reads.id);
		}

		@Override
		public void addLines(final ListNum<SamLine> lines){
			throw new UnsupportedOperationException();
		}

		@Override
		public void poison(){}

		@Override
		public boolean waitForFinish(){return false;}

		@Override
		public boolean poisonAndWait(){return false;}

		@Override
		public void finishError(){}

		@Override
		public String fname(){return name;}

		@Override
		public boolean errorState(){return false;}

		@Override
		public boolean finishedSuccessfully(){return true;}

		final String name;
		private final ArrayList<Call> calls=new ArrayList<Call>();
		private CountDownLatch entered, release;
		private RuntimeException failure;
	}

	/** One fake add call. */
	private static final class Call{

		Call(final long id_, final ArrayList<Read> list_, final Thread thread_){
			id=id_;
			list=list_;
			thread=thread_;
		}

		final long id;
		final ArrayList<Read> list;
		final Thread thread;
	}

	private static final byte[] BASES="ACGT".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] QUALS=new byte[] {40, 40, 40, 40};

	private KmerSort3ParallelDestinationWriterTest(){}
}
