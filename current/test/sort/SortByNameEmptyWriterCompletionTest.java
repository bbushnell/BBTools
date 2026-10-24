package test.sort;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import sort.SortByName;
import stream.Read;

/** Focused completion test for SortByName's empty direct-output writer.
 * @author Jean
 */
public final class SortByNameEmptyWriterCompletionTest{

	/** Runs expect=baseline|fixed in a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		if(args.length!=2 || (!args[0].equals("expect=baseline") && !args[0].equals("expect=fixed"))){
			throw new IllegalArgumentException("Expected expect=baseline|expect=fixed and a fresh fixture directory");
		}
		final boolean fixed=args[0].equals("expect=fixed");
		final Path root=Paths.get(args[1]).toAbsolutePath();
		Files.createDirectory(root);
		checkEmptyWriterCompletion(root, fixed);
	}

	/** Blocks the real empty WriteThread and checks whether waitOnMemory waits for it. */
	private static void checkEmptyWriterCompletion(Path root, boolean fixed) throws Exception{
		final Path input=root.resolve("input.fq"), output=root.resolve("out.fq");
		Files.write(input, READ.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final SortByName tool=new SortByName(new String[]{"in="+input, "out="+output,
			"delete=f", "overwrite=t"});
		final CountDownLatch writerReady=new CountDownLatch(1), writerRelease=new CountDownLatch(1);
		final CountDownLatch waiterEntered=new CountDownLatch(1), waiterReturned=new CountDownLatch(1);
		final AtomicLong outstanding=new AtomicLong(0);
		final AtomicReference<Thread> writer=new AtomicReference<Thread>();
		final AtomicReference<Throwable> waiterFailure=new AtomicReference<Throwable>();
		final GatedEmptyList storage=new GatedEmptyList(writer, writerReady, writerRelease);
		final Method waitOnMemory=waitOnMemory();
		Thread waiter=null;
		boolean returnedEarly=false;
		try{
			sortAndDump().invoke(tool, storage, 0L, outstanding, output.toString(), null, false);
			await(writerReady, "SortByName WriteThread did not enter empty storage.size().");
			waiter=new Thread("SortByNameEmptyWriterCompletionTest-waiter"){
				@Override
				public void run(){
					try{
						waiterEntered.countDown();
						waitOnMemory.invoke(tool, outstanding, 0L);
					}
					catch(Throwable t){waiterFailure.set(t);}
					finally{waiterReturned.countDown();}
				}
			};
			waiter.setDaemon(true);
			waiter.start();
			await(waiterEntered, "waitOnMemory waiter did not start.");
			returnedEarly=waiterReturned.await(250, TimeUnit.MILLISECONDS);
			if(fixed){
				assert(!returnedEarly) : "waitOnMemory(..., 0) returned before the empty writer completed.";
			}else{
				assert(returnedEarly) : "Expected baseline waitOnMemory(..., 0) to miss the empty writer.";
				assert(writer.get()!=null && writer.get().isAlive()) : "Writer was not gated after early wait return.";
			}
			writerRelease.countDown();
		}finally{
			writerRelease.countDown();
			try{joinIfStarted(writer.get(), "SortByName WriteThread");}
			finally{joinIfStarted(waiter, "waitOnMemory waiter");}
		}
		if(waiterFailure.get()!=null){throw rethrow(waiterFailure.get());}
		assert(outstanding.get()==0) : "Writers joined with outstanding memory "+outstanding.get();
		System.out.println("empty-writer-completion\t"+returnedEarly+'\t'+outstanding.get());
	}

	/** Finds SortByName's private direct-output writer launcher. */
	private static Method sortAndDump() throws Exception{
		final Method method=SortByName.class.getDeclaredMethod("sortAndDump",
			ArrayList.class, Long.TYPE, AtomicLong.class, String.class, String.class, Boolean.TYPE);
		method.setAccessible(true);
		return method;
	}

	/** Finds SortByName's private memory gate. */
	private static Method waitOnMemory() throws Exception{
		final Method method=SortByName.class.getDeclaredMethod("waitOnMemory",
			AtomicLong.class, Long.TYPE);
		method.setAccessible(true);
		return method;
	}

	/** Waits for a latch without letting a deadlock become the test result. */
	private static void await(CountDownLatch latch, String message) throws InterruptedException{
		if(!latch.await(5, TimeUnit.SECONDS)){throw new AssertionError(message);}
	}

	/** Waits for a thread without letting a deadlock become the test result. */
	private static void join(Thread thread, String name) throws InterruptedException{
		if(thread==null){throw new AssertionError(name+" was never started.");}
		thread.join(5000);
		if(thread.isAlive()){throw new AssertionError(name+" did not finish after release.");}
	}

	/** Waits for a thread if its launch point was reached. */
	private static void joinIfStarted(Thread thread, String name) throws InterruptedException{
		if(thread!=null){join(thread, name);}
	}

	/** Unwraps reflection failures on a helper thread. */
	private static RuntimeException rethrow(Throwable t){
		final Throwable cause=(t instanceof InvocationTargetException ?
			((InvocationTargetException)t).getCause() : t);
		if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
		if(cause instanceof Error){throw (Error)cause;}
		throw new RuntimeException(cause);
	}

	private static final String READ="@read0\nACGT\n+\nIIII\n";

	/** Empty storage that blocks the first real writer thread size() call. */
	private static final class GatedEmptyList extends ArrayList<Read>{

		GatedEmptyList(AtomicReference<Thread> writer_, CountDownLatch ready_, CountDownLatch release_){
			owner=Thread.currentThread();
			writer=writer_;
			ready=ready_;
			release=release_;
		}

		@Override
		public int size(){
			if(Thread.currentThread()!=owner && writer.compareAndSet(null, Thread.currentThread())){
				ready.countDown();
				try{release.await();}
				catch(InterruptedException e){throw new RuntimeException(e);}
			}
			return super.size();
		}

		private final Thread owner;
		private final AtomicReference<Thread> writer;
		private final CountDownLatch ready;
		private final CountDownLatch release;

		private static final long serialVersionUID=1L;
	}
}
