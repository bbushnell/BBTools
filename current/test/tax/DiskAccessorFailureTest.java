package test.tax;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import tax.DiskAccessionTable;
import tax.NIODiskAccessor;
import tax.SerialDiskAccessor;

/** Finite failure-path checks for disk-backed accession accessors.
 * @author Jean */
public final class DiskAccessorFailureTest{

	/** Runs baseline or fixed expectations for selected disk-accessor failures.
	 * @param args Optional expect=baseline or expect=fixed and startup=baseline or
	 *             startup=fixed
	 * @throws Exception If fixture construction fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		Boolean startupBaseline=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("startup=baseline")){startupBaseline=Boolean.TRUE;}
			else if(arg.equals("startup=fixed")){startupBaseline=Boolean.FALSE;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final Path root=Files.createTempDirectory("disk-accessor-failure-");
		try{
			final Path full=root.resolve("accession.bin");
			final Path truncated=root.resolve("accession.truncated.bin");
			buildTinyTable(full);
			truncateToHeaderAndPartitions(full, truncated);
			final long[] multiWorkerKeys=multiWorkerKeys();
			final int[] expectedMultiWorker=expectedMultiWorkerResults(multiWorkerKeys);

			DiskAccessionTable table=null;
			try{
				table=new DiskAccessionTable(full.toString());
				final int[] serial=new SerialDiskAccessor(table).getBatch(new long[]{1, 2, 3, 4, 99});
				check(Arrays.equals(serial, new int[] {11, 22, 33, 44, -1}),
						"Serial mapped lookup control failed: "+Arrays.toString(serial));
			}finally{
				if(table!=null){table.close();}
			}

			NIODiskAccessor nio=null;
			try{
				nio=new NIODiskAccessor(full.toString(), 4);
				final int[] fullResult=nio.getBatch(new long[]{1, 2, 3, 4, 99});
				check(Arrays.equals(fullResult, new int[] {11, 22, 33, 44, -1}),
						"NIO positioned lookup control failed: "+Arrays.toString(fullResult));
				final int[] multiWorkerResult=nio.getBatch(multiWorkerKeys);
				check(Arrays.equals(multiWorkerResult, expectedMultiWorker),
						"NIO multiworker lookup control failed: "+
						Arrays.toString(multiWorkerResult));
			}finally{
				if(nio!=null){nio.close();}
			}

			checkTruncatedRead(truncated, baseline);
			checkInterruptedJoin(full, baseline);
			if(startupBaseline!=null){
				for(StartupFailureMode mode : StartupFailureMode.values()){
					checkStartupFailure(full, startupBaseline.booleanValue(), mode);
				}
			}
		}finally{
			deleteTempTree(root);
		}
		System.out.println("DISK_ACCESSOR_FAILURE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed")+
				(startupBaseline==null ? "" : " startup="+
				(startupBaseline.booleanValue() ? "baseline" : "fixed")));
	}

	/** Prevents construction of this standalone fixture. */
	private DiskAccessorFailureTest(){}

	/** Builds a small valid accession table with four present keys.
	 * @param path Output table path
	 * @throws Exception If table construction fails */
	private static void buildTinyTable(final Path path) throws Exception{
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(path.toString(), 4)){
			builder.add(1, 11);
			builder.add(2, 22);
			builder.add(3, 33);
			builder.add(4, 44);
		}
	}

	/** Keeps only the readable table header and partition table, removing all slot data.
	 * @param input Complete source table
	 * @param output Truncated destination table
	 * @throws Exception If the table cannot be read or written */
	private static void truncateToHeaderAndPartitions(final Path input, final Path output)
			throws Exception{
		final byte[] bytes=Files.readAllBytes(input);
		final int numPartitions=readInt(bytes, 16);
		// DiskAccessionTable writes 24 bytes of header values, then partition
		// metadata at byte 24, but it starts slot data after a 32-byte legacy
		// header pad plus 16 bytes per partition.
		final int headerBytes=32+numPartitions*16;
		Files.write(output, Arrays.copyOf(bytes, headerBytes));
	}

	/** Asserts EOF in a nonempty table is not silently converted to an absent key.
	 * @param truncated Table with metadata but no slot data
	 * @param baseline True to expect the current silent -1 behavior
	 * @throws Exception If the accessor cannot be opened or closed */
	private static void checkTruncatedRead(final Path truncated, final boolean baseline)
			throws Exception{
		NIODiskAccessor nio=null;
		try{
			nio=new NIODiskAccessor(truncated.toString(), 1);
			final int[] result=nio.getBatch(new long[]{1});
			if(baseline){
				check(result.length==1 && result[0]==-1,
						"Baseline must expose silent truncated-table failure: "+
						Arrays.toString(result));
			}else{
				throw new AssertionError("Truncated slot data must fail loudly, not return "+
						Arrays.toString(result));
			}
		}catch(RuntimeException e){
			if(baseline){throw e;}
			checkTruncatedFailure(e);
		}finally{
			if(nio!=null){nio.close();}
		}
		try{
			nio=new NIODiskAccessor(truncated.toString(), 4);
			final int[] result=nio.getBatch(multiWorkerKeys());
			if(baseline){
				check(allMissing(result),
						"Baseline multiworker truncated-table result should silently "+
						"report only absent keys: "+Arrays.toString(result));
			}else{
				throw new AssertionError("Multiworker truncated slot data must fail "+
						"loudly, not return "+Arrays.toString(result));
			}
		}catch(RuntimeException e){
			if(baseline){throw e;}
			checkTruncatedFailure(e);
		}finally{
			if(nio!=null){nio.close();}
		}
	}

	/** Asserts a join interruption is restored and propagated instead of ignored.
	 * @param full Complete table used by worker threads
	 * @param baseline True to expect current ignored-join behavior
	 * @throws Exception If the accessor cannot be opened or closed */
	private static void checkInterruptedJoin(final Path full, final boolean baseline)
			throws Exception{
		final CountDownLatch callerDone=new CountDownLatch(1);
		final AtomicReference<int[]> result=new AtomicReference<int[]>();
		final AtomicReference<Throwable> thrown=new AtomicReference<Throwable>();
		final AtomicBoolean callerInterrupted=new AtomicBoolean(false);
		DeterministicNIODiskAccessor nio=null;
		Thread caller=null;
		try{
			nio=new DeterministicNIODiskAccessor(full.toString(), 4);
			final NIODiskAccessor accessor=nio;
			caller=new Thread(()->{
				try{
					result.set(accessor.getBatch(missingKeys(100)));
				}catch(Throwable t){
					thrown.set(t);
				}finally{
					callerInterrupted.set(Thread.currentThread().isInterrupted());
					callerDone.countDown();
				}
			}, "disk-accessor-interrupted-join");
			caller.start();
			check(nio.awaitFirstWorkerStarted(),
					"Worker zero did not reach the deterministic test latch");
			caller.interrupt();
			if(baseline){
				check(callerDone.await(WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"Baseline getBatch must return before blocked worker zero is released");
				check(thrown.get()==null,
						"Baseline getBatch should ignore interruption, not throw "+thrown.get());
				check(containsUnfilledSlot(result.get()),
						"Baseline getBatch must expose a partially filled result: "+
						Arrays.toString(result.get()));
				check(!callerInterrupted.get(),
						"Baseline getBatch must expose that join cleared the interrupt");
				check(!nio.firstWorkerFinished(),
						"Baseline getBatch returned after worker zero was already joined");
			}else{
				check(!callerDone.await(200, TimeUnit.MILLISECONDS),
						"Fixed getBatch must not return while worker zero remains blocked");
				nio.releaseFirstWorker();
				check(callerDone.await(WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"Fixed getBatch did not propagate interruption after worker release");
				check(thrown.get() instanceof RuntimeException,
						"Fixed getBatch must fail loudly; observed "+thrown.get());
				check(callerInterrupted.get(),
						"Fixed getBatch must restore the interrupt before failure");
			}
		}finally{
			try{
				if(nio!=null){nio.releaseFirstWorker();}
				joinBounded(caller, "caller");
			}finally{
				try{
					if(nio!=null){joinBounded(nio.workers(), "NIO worker");}
				}finally{
					if(nio!=null){nio.close();}
				}
			}
		}
	}

	/** Asserts startup failures wait for already-started workers before returning.
	 * @param full Complete table used by worker threads
	 * @param baseline True to expect current unjoined startup-failure behavior
	 * @param mode Startup failure injection point to exercise
	 * @throws Exception If the accessor cannot be opened or closed */
	private static void checkStartupFailure(final Path full, final boolean baseline,
			final StartupFailureMode mode) throws Exception{
		final CountDownLatch callerDone=new CountDownLatch(1);
		final AtomicReference<int[]> result=new AtomicReference<int[]>();
		final AtomicReference<Throwable> thrown=new AtomicReference<Throwable>();
		final AtomicBoolean callerInterrupted=new AtomicBoolean(false);
		final AtomicBoolean callerReturnedBeforeRelease=new AtomicBoolean(false);
		final AtomicReference<Integer> liveWorkersAtCallerDone=
				new AtomicReference<Integer>();
		DeterministicNIODiskAccessor nio=null;
		Thread caller=null;
		try{
			nio=new DeterministicNIODiskAccessor(full.toString(), 4, mode);
			final DeterministicNIODiskAccessor accessor=nio;
			caller=new Thread(()->{
				try{
					result.set(accessor.getBatch(missingKeys(100)));
				}catch(Throwable t){
					thrown.set(t);
				}finally{
					liveWorkersAtCallerDone.set(accessor.liveStartedWorkerCount());
					callerInterrupted.set(Thread.currentThread().isInterrupted());
					callerDone.countDown();
				}
			}, "disk-accessor-startup-"+mode.name().toLowerCase(Locale.ROOT));
			caller.start();
			check(nio.awaitStartupFailureReached(),
					mode+" startup failure did not reach its injection point");
			if(baseline){
				check(callerDone.await(WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"Baseline "+mode+" startup failure did not reach the caller");
				callerReturnedBeforeRelease.set(true);
				check(thrown.get()==nio.startupFailure(),
						"Baseline "+mode+" must expose the exact startup sentinel, not "+
						thrown.get()+" with result "+Arrays.toString(result.get()));
				check(!nio.firstWorkerFinished(),
						"Baseline "+mode+" returned after worker zero was already joined");
				if(mode==StartupFailureMode.START_AFTER_SUPER){
					check(!nio.startupFailureWorkerFinished(),
							"Baseline START_AFTER_SUPER returned after worker one was joined");
				}
				checkFactoryInterruptPreserved(callerInterrupted, mode, "Baseline");
				nio.releaseFirstWorker();
				nio.releaseStartupFailureWorker();
				check(nio.awaitStartedWorkersFinished(),
						"Baseline "+mode+" workers did not finish after fixture release");
				joinBounded(nio.workers(), "baseline NIO worker");
			}else{
				check(!callerDone.await(200, TimeUnit.MILLISECONDS),
						"Fixed "+mode+" must not return while worker zero remains blocked");
				nio.releaseFirstWorker();
				if(mode==StartupFailureMode.START_AFTER_SUPER){
					check(!callerDone.await(200, TimeUnit.MILLISECONDS),
							"Fixed START_AFTER_SUPER must not return after only "+
							"worker zero is released");
					nio.releaseStartupFailureWorker();
				}
				check(callerDone.await(WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"Fixed "+mode+" did not rethrow startup failure after worker release");
				check(nio.firstWorkerFinished(),
						"Fixed "+mode+" did not join worker zero before returning");
				if(mode==StartupFailureMode.START_AFTER_SUPER){
					check(nio.startupFailureWorkerFinished(),
							"Fixed START_AFTER_SUPER did not join worker one before returning");
				}
				check(thrown.get()==nio.startupFailure(),
						"Fixed "+mode+" must preserve the exact startup sentinel, not "+
						thrown.get()+" with result "+Arrays.toString(result.get()));
				checkFactoryInterruptPreserved(callerInterrupted, mode, "Fixed");
			}
			checkLiveWorkersAtCallerDone(liveWorkersAtCallerDone.get(), baseline, mode);
			printStartupReceipt(mode, baseline, callerReturnedBeforeRelease.get(),
					thrown.get()==nio.startupFailure(), callerInterrupted.get(),
					liveWorkersAtCallerDone.get());
			if(mode==StartupFailureMode.START_BEFORE_SUPER){
				check(nio.startupFailureThreadState()==Thread.State.NEW,
						"Failed START_BEFORE_SUPER worker must remain NEW, not "+
						nio.startupFailureThreadState());
			}
			nio.checkNoStartedWorkerAlive();
		}finally{
			try{
				if(nio!=null){
					nio.releaseFirstWorker();
					nio.releaseStartupFailureWorker();
				}
				joinBounded(caller, "caller");
			}finally{
				try{
					if(nio!=null){joinBounded(nio.workers(), "NIO worker");}
				}finally{
					if(nio!=null){nio.close();}
				}
			}
		}
	}

	/** Checks the return-time live-worker receipt before printing it.
	 * @param liveWorkersAtCallerDone Retained live worker count when caller returned
	 * @param baseline True if the case expected current baseline behavior
	 * @param mode Startup failure injection point exercised */
	private static void checkLiveWorkersAtCallerDone(
			final Integer liveWorkersAtCallerDone, final boolean baseline,
			final StartupFailureMode mode){
		final int expected=baseline ?
				(mode==StartupFailureMode.START_AFTER_SUPER ? 2 : 1) : 0;
		check(liveWorkersAtCallerDone!=null &&
				liveWorkersAtCallerDone.intValue()==expected,
				"Expected "+expected+" live workers when "+mode+
				" caller returned under "+(baseline ? "baseline" : "fixed")+
				", not "+liveWorkersAtCallerDone);
	}

	/** Prints one startup-failure artifact line.
	 * @param mode Startup failure injection point exercised
	 * @param baseline True if the case expected current baseline behavior
	 * @param callerReturnedBeforeRelease True if caller returned before fixture release
	 * @param exactPrimary True if the observed primary was the retained sentinel
	 * @param callerInterrupted Final caller-thread interrupt state
	 * @param liveWorkersAtCallerDone Live started workers when the caller returned */
	private static void printStartupReceipt(final StartupFailureMode mode,
			final boolean baseline, final boolean callerReturnedBeforeRelease,
			final boolean exactPrimary, final boolean callerInterrupted,
			final Integer liveWorkersAtCallerDone){
		System.out.println("DISK_ACCESSOR_STARTUP mode="+mode+
				" expect="+(baseline ? "baseline" : "fixed")+
				" caller_before_release="+callerReturnedBeforeRelease+
				" exact_primary="+exactPrimary+
				" caller_interrupted="+callerInterrupted+
				" live_workers_at_caller="+liveWorkersAtCallerDone);
	}

	/** Checks that startup cleanup preserves a pending caller interrupt.
	 * @param callerInterrupted Final caller-thread interrupt state
	 * @param mode Startup failure injection point exercised
	 * @param expectationName Baseline or fixed diagnostic label */
	private static void checkFactoryInterruptPreserved(
			final AtomicBoolean callerInterrupted, final StartupFailureMode mode,
			final String expectationName){
		if(mode==StartupFailureMode.FACTORY_INTERRUPTED){
			check(callerInterrupted.get(),
					expectationName+" FACTORY_INTERRUPTED startup must preserve "+
					"the caller interrupt flag");
		}
	}

	/** Checks that truncated slot data raised the intended I/O failure.
	 * @param e Failure propagated by a fixed NIODiskAccessor */
	private static void checkTruncatedFailure(final RuntimeException e){
		final Throwable cause=e.getCause();
		check(cause instanceof IOException,
				"Truncated slot data must propagate an IOException cause, not "+cause);
		final String message=String.valueOf(e.getMessage()).toLowerCase(Locale.ROOT);
		final String causeMessage=String.valueOf(cause.getMessage()).toLowerCase(Locale.ROOT);
		check(containsReadFailureDetail(message),
				"Truncated slot RuntimeException message should describe the read failure: "+
				e.getMessage());
		check(containsReadFailureDetail(causeMessage),
				"Truncated slot IOException message should describe the read failure: "+
				cause.getMessage());
	}

	/** Reports whether an exception message describes the truncated-slot failure.
	 * @param message Lower-case message text
	 * @return True if the text names EOF, truncation, or slot reads */
	private static boolean containsReadFailureDetail(final String message){
		return message.contains("eof") || message.contains("end") ||
				message.contains("truncat") || message.contains("slot");
	}

	/** Reports whether a missing-key result returned before a worker filled its slice.
	 * @param result Observed missing-key lookup result
	 * @return True if at least one element remains at the int default */
	private static boolean containsUnfilledSlot(final int[] result){
		if(result==null){return false;}
		for(int value : result){
			if(value==0){return true;}
		}
		return false;
	}

	/** Reports whether every lookup was silently reported as the missing-key sentinel.
	 * @param result Observed truncated-table lookup result
	 * @return True if all results are exactly -1 */
	private static boolean allMissing(final int[] result){
		if(result==null || result.length!=MIN_MULTIWORKER_KEYS){return false;}
		for(int value : result){
			if(value!=-1){return false;}
		}
		return true;
	}

	/** Joins a thread with a bounded wait before fixture teardown.
	 * @param thread Thread to join, or null
	 * @param name Diagnostic thread role name
	 * @throws InterruptedException If the fixture thread is interrupted */
	private static void joinBounded(final Thread thread, final String name)
			throws InterruptedException{
		if(thread==null){return;}
		thread.join(WAIT_MILLIS);
		check(!thread.isAlive(), name+" did not finish within "+WAIT_MILLIS+" ms");
	}

	/** Joins all retained worker threads with bounded waits before closing the accessor.
	 * @param threads Worker snapshot
	 * @param name Diagnostic thread role name
	 * @throws InterruptedException If the fixture thread is interrupted */
	private static void joinBounded(final Thread[] threads, final String name)
			throws InterruptedException{
		for(Thread thread : threads){joinBounded(thread, name);}
	}

	/** Creates positive keys that are absent from the tiny table.
	 * @param length Number of keys to return
	 * @return Missing positive keys */
	private static long[] missingKeys(final int length){
		final long[] keys=new long[length];
		for(int i=0; i<keys.length; i++){keys[i]=1000+i;}
		return keys;
	}

	/** Creates a batch split across all four NIO test workers with one known ID per slice.
	 * @return Mixed present and missing keys */
	private static long[] multiWorkerKeys(){
		final long[] keys=missingKeys(MIN_MULTIWORKER_KEYS);
		keys[0]=1;
		keys[MIN_MULTIWORKER_KEYS/4]=2;
		keys[MIN_MULTIWORKER_KEYS/2]=3;
		keys[(MIN_MULTIWORKER_KEYS*3)/4]=4;
		return keys;
	}

	/** Builds expected IDs for the mixed multiworker lookup batch.
	 * @param keys Mixed present and missing keys
	 * @return Expected tax IDs by query position */
	private static int[] expectedMultiWorkerResults(final long[] keys){
		final int[] expected=new int[keys.length];
		Arrays.fill(expected, -1);
		expected[0]=11;
		expected[MIN_MULTIWORKER_KEYS/4]=22;
		expected[MIN_MULTIWORKER_KEYS/2]=33;
		expected[(MIN_MULTIWORKER_KEYS*3)/4]=44;
		return expected;
	}

	/** Reads a big-endian int from a byte array.
	 * @param bytes Source bytes
	 * @param offset First byte to read
	 * @return Big-endian integer value */
	private static int readInt(final byte[] bytes, final int offset){
		return ((bytes[offset]&0xFF)<<24) | ((bytes[offset+1]&0xFF)<<16) |
				((bytes[offset+2]&0xFF)<<8) | (bytes[offset+3]&0xFF);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Startup failure injection point. */
	private enum StartupFailureMode{
		/** Throw from the worker factory before returning worker one. */
		FACTORY,
		/** Interrupt the caller and throw from the worker factory. */
		FACTORY_INTERRUPTED,
		/** Throw from Thread.start() before worker one is actually started. */
		START_BEFORE_SUPER,
		/** Start worker one, then throw from the same Thread.start() call. */
		START_AFTER_SUPER
	}

	/** NIO accessor that can hold its first worker until the fixture releases it. */
	private static final class DeterministicNIODiskAccessor extends NIODiskAccessor{

		/** Opens a complete table and installs an instance-owned worker seam.
		 * @param path Table file path
		 * @param maxThreads Maximum worker threads per batch
		 * @throws IOException If the table cannot be opened */
		DeterministicNIODiskAccessor(final String path, final int maxThreads)
				throws IOException{
			this(path, maxThreads, null);
		}

		/** Opens a complete table and installs a startup-failure worker seam.
		 * @param path Table file path
		 * @param maxThreads Maximum worker threads per batch
		 * @param startupFailureMode_ Startup failure to inject, or null
		 * @throws IOException If the table cannot be opened */
		DeterministicNIODiskAccessor(final String path, final int maxThreads,
				final StartupFailureMode startupFailureMode_) throws IOException{
			super(path, maxThreads);
			startupFailureMode=startupFailureMode_;
		}

		/** Creates one retained worker, wrapping only the first with latches.
		 * @param task Slice reader to run
		 * @return New unstarted worker */
		@Override
		protected Thread makeWorker(final Runnable task){
			synchronized(workers){
				if(workers.size()==1 && startupFailureMode!=null){
					return makeStartupFailureWorker(task);
				}
				final Runnable wrappedTask=workers.isEmpty() ? wrapFirstWorker(task) : task;
				final Thread worker=super.makeWorker(wrappedTask);
				workers.add(worker);
				return worker;
			}
		}

		/** Releases the first worker if it is blocked at the fixture latch. */
		void releaseFirstWorker(){firstWorkerRelease.countDown();}

		/** Releases the failed START_AFTER_SUPER worker if it is blocked. */
		void releaseStartupFailureWorker(){startupFailureWorkerRelease.countDown();}

		/** Awaits first-worker startup with the fixture's bounded wait.
		 * @return True if the first worker reached the latch
		 * @throws InterruptedException If the fixture thread is interrupted */
		boolean awaitFirstWorkerStarted() throws InterruptedException{
			return firstWorkerStarted.await(WAIT_MILLIS, TimeUnit.MILLISECONDS);
		}

		/** Awaits the injected startup failure with the fixture's bounded wait.
		 * @return True if the failure reached its injection point
		 * @throws InterruptedException If the fixture thread is interrupted */
		boolean awaitStartupFailureReached() throws InterruptedException{
			return startupFailureReached.await(WAIT_MILLIS, TimeUnit.MILLISECONDS);
		}

		/** Awaits all workers that were deliberately started by this fixture.
		 * @return True if all started workers finished before the bounded timeout
		 * @throws InterruptedException If the fixture thread is interrupted */
		boolean awaitStartedWorkersFinished() throws InterruptedException{
			if(!firstWorkerFinished.await(WAIT_MILLIS, TimeUnit.MILLISECONDS)){
				return false;
			}
			return startupFailureMode!=StartupFailureMode.START_AFTER_SUPER ||
					startupFailureWorkerFinished.await(WAIT_MILLIS, TimeUnit.MILLISECONDS);
		}

		/** Reports whether the first worker has exited.
		 * @return True if the first worker exit latch has been counted down */
		boolean firstWorkerFinished(){return firstWorkerFinished.getCount()==0;}

		/** Reports whether the failed START_AFTER_SUPER worker has exited.
		 * @return True if worker one was never started or has exited */
		boolean startupFailureWorkerFinished(){
			return startupFailureMode!=StartupFailureMode.START_AFTER_SUPER ||
					startupFailureWorkerFinished.getCount()==0;
		}

		/** Returns the startup failure sentinel.
		 * @return Failure sentinel thrown by this fixture */
		RuntimeException startupFailure(){return startupFailure;}

		/** Returns the failed startup worker's state.
		 * @return Failed startup worker state, or null if factory creation failed */
		Thread.State startupFailureThreadState(){
			return startupFailureThread==null ? null : startupFailureThread.getState();
		}

		/** Returns a snapshot of workers created so far.
		 * @return Worker thread snapshot */
		Thread[] workers(){
			synchronized(workers){return workers.toArray(new Thread[0]);}
		}

		/** Counts retained workers that were started and have not yet exited.
		 * @return Number of retained started workers still alive */
		Integer liveStartedWorkerCount(){
			int live=0;
			synchronized(workers){
				for(Thread worker : workers){
					if(worker.getState()!=Thread.State.NEW && worker.isAlive()){live++;}
				}
			}
			return live;
		}

		/** Asserts no retained started worker is alive before fixture teardown joins. */
		void checkNoStartedWorkerAlive(){
			synchronized(workers){
				for(Thread worker : workers){
					if(worker.getState()!=Thread.State.NEW){
						check(!worker.isAlive(),
								"Fixture teardown would hide a live started worker: "+
								worker.getName()+" "+worker.getState());
					}
				}
			}
		}

		/** Creates or throws the second-worker startup failure.
		 * @param task Slice reader to run for start-after-super failures
		 * @return New retained worker that throws from start()
		 * @throws RuntimeException Always for factory failures */
		private Thread makeStartupFailureWorker(final Runnable task){
			awaitFirstWorkerBeforeStartupFailure();
			if(startupFailureMode==StartupFailureMode.FACTORY ||
					startupFailureMode==StartupFailureMode.FACTORY_INTERRUPTED){
				throwStartupFailure();
			}
			final Runnable wrappedTask=startupFailureMode==StartupFailureMode.START_AFTER_SUPER ?
					wrapStartupFailureWorker(task) : task;
			startupFailureThread=new Thread(wrappedTask){
				@Override
				public synchronized void start(){
					if(startupFailureMode==StartupFailureMode.START_AFTER_SUPER){
						super.start();
						awaitStartupFailureWorkerStarted();
					}
					throwStartupFailure();
				}
			};
			workers.add(startupFailureThread);
			return startupFailureThread;
		}

		/** Wraps the first worker so the fixture can hold it before slot reads.
		 * @param task Worker task to run after release
		 * @return Wrapped first-worker task */
		private Runnable wrapFirstWorker(final Runnable task){
			return ()->{
				try{
					firstWorkerStarted.countDown();
					firstWorkerRelease.await();
					task.run();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
				}finally{
					firstWorkerFinished.countDown();
				}
			};
		}

		/** Wraps worker one so START_AFTER_SUPER can prove cleanup joins it.
		 * @param task Worker task to run after release
		 * @return Wrapped second-worker task */
		private Runnable wrapStartupFailureWorker(final Runnable task){
			return ()->{
				try{
					startupFailureWorkerStarted.countDown();
					startupFailureWorkerRelease.await();
					task.run();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
				}finally{
					startupFailureWorkerFinished.countDown();
				}
			};
		}

		/** Waits for worker zero before injecting failure into worker one. */
		private void awaitFirstWorkerBeforeStartupFailure(){
			try{
				check(awaitFirstWorkerStarted(),
						"Worker zero did not start before "+startupFailureMode+
						" startup failure");
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new AssertionError("Interrupted before "+startupFailureMode+
						" startup failure", e);
			}
		}

		/** Waits for the failed START_AFTER_SUPER worker before throwing failure. */
		private void awaitStartupFailureWorkerStarted(){
			try{
				check(startupFailureWorkerStarted.await(WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"Worker one did not start before "+startupFailureMode+
						" startup failure");
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new AssertionError("Interrupted before worker one "+
						startupFailureMode+" startup failure", e);
			}
		}

		/** Marks and throws the retained startup failure sentinel. */
		private void throwStartupFailure(){
			if(startupFailureMode==StartupFailureMode.FACTORY_INTERRUPTED){
				Thread.currentThread().interrupt();
			}
			startupFailureReached.countDown();
			throw startupFailure;
		}

		/** First-worker startup latch. */
		private final CountDownLatch firstWorkerStarted=new CountDownLatch(1);
		/** First-worker release latch. */
		private final CountDownLatch firstWorkerRelease=new CountDownLatch(1);
		/** First-worker exit latch. */
		private final CountDownLatch firstWorkerFinished=new CountDownLatch(1);
		/** Startup failure injection-point latch. */
		private final CountDownLatch startupFailureReached=new CountDownLatch(1);
		/** START_AFTER_SUPER worker startup latch. */
		private final CountDownLatch startupFailureWorkerStarted=new CountDownLatch(1);
		/** START_AFTER_SUPER worker release latch. */
		private final CountDownLatch startupFailureWorkerRelease=new CountDownLatch(1);
		/** START_AFTER_SUPER worker exit latch. */
		private final CountDownLatch startupFailureWorkerFinished=new CountDownLatch(1);
		/** Workers created by the overridden factory. */
		private final ArrayList<Thread> workers=new ArrayList<Thread>();
		/** Startup failure to inject, or null for ordinary deterministic workers. */
		private final StartupFailureMode startupFailureMode;
		/** Startup failure sentinel. */
		private final RuntimeException startupFailure=
				new RuntimeException("Injected NIO worker startup failure.");
		/** Failed startup worker, or null for factory failures. */
		private Thread startupFailureThread;
	}

	/** Bounded wait for worker scheduling and cleanup. */
	private static final long WAIT_MILLIS=5000;
	/** Smallest batch that current NIODiskAccessor settings split over four workers. */
	private static final int MIN_MULTIWORKER_KEYS=200;
}
