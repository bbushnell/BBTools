package test.tax;

import java.io.File;
import java.io.RandomAccessFile;
import java.lang.management.LockInfo;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import shared.Shared;
import tax.AnalyzeAccession;
import tax.BuildDiskTables;
import tax.DiskAccessionTable;

/** Finite baseline for BuildDiskTables worker joins and input cleanup.
 * @author Jean, Shinobu */
public final class BuildDiskTablesWorkerFailureTest{

	/** Accession2taxid input with one header row and one valid accession row. */
	private static final String HEALTHY=
			"accession\taccession.version\ttaxid\tgi\n"+
			"ABC123\tABC123.1\t7\t0\n";
	/** Exact pattern table matching {@link #HEALTHY}. */
	private static final String PATTERNS=
			"#Pattern\tCount\tCombos\tBits\n"+
			"LLLDDD\t1\t17576000\t24.07\n";

	/** Runs the parent harness or one isolated child.
	 * @param args Optional expect=baseline|fixed, count=baseline|fixed,
	 * reader=baseline|fixed and startup=baseline|fixed; child mode is internal
	 * @throws Exception If fixture construction, child launch, or reflection fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false, countBaseline=false, countSet=false;
		boolean readerBaseline=false, readerSet=false, startupBaseline=false;
		boolean startupSet=false, child=false;
		Path childRoot=null;
		String childCase="failure";
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("count=baseline")){
				countBaseline=true;
				countSet=true;
			}else if(arg.equals("count=fixed")){
				countBaseline=false;
				countSet=true;
			}else if(arg.equals("reader=baseline")){
				readerBaseline=true;
				readerSet=true;
			}else if(arg.equals("reader=fixed")){
				readerBaseline=false;
				readerSet=true;
			}else if(arg.equals("startup=baseline")){
				startupBaseline=true;
				startupSet=true;
			}else if(arg.equals("startup=fixed")){
				startupBaseline=false;
				startupSet=true;
			}else if(arg.equals("child")){child=true;}
			else if(arg.startsWith("case=")){childCase=arg.substring(5);}
			else if(arg.startsWith("root=")){childRoot=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(!countSet){countBaseline=baseline;}
		if(!readerSet){readerBaseline=baseline;}
		if(!startupSet){startupBaseline=baseline;}

		if(child){
			check(childRoot!=null, "Child mode requires root=...");
			runChild(childCase, childRoot, baseline, countBaseline, readerBaseline,
					startupBaseline);
			return;
		}

		final Path root=Files.createTempDirectory("build-disk-tables-worker-");
		try{
			checkChild(root, "failure", baseline, countBaseline, readerBaseline, true);
			checkChild(root, "failure", baseline, countBaseline, readerBaseline, false);
			checkChild(root, "count-failure", baseline, countBaseline, readerBaseline, true);
			checkChild(root, "count-failure", baseline, countBaseline, readerBaseline, false);
			checkChild(root, "reader", baseline, countBaseline, readerBaseline, true);
			checkChild(root, "reader", baseline, countBaseline, readerBaseline, false);
			checkChild(root, "count-reader", baseline, countBaseline, readerBaseline, true);
			checkChild(root, "count-reader", baseline, countBaseline, readerBaseline, false);
			checkChild(root, "startup-count", baseline, countBaseline, readerBaseline,
					startupBaseline, true);
			checkChild(root, "startup-count", baseline, countBaseline, readerBaseline,
					startupBaseline, false);
			checkChild(root, "startup-build", baseline, countBaseline, readerBaseline,
					startupBaseline, true);
			checkChild(root, "startup-build", baseline, countBaseline, readerBaseline,
					startupBaseline, false);
			if(!startupBaseline){
				checkChild(root, "startup-count-interrupted", baseline, countBaseline,
						readerBaseline, startupBaseline, true);
				checkChild(root, "startup-count-interrupted", baseline, countBaseline,
						readerBaseline, startupBaseline, false);
				checkChild(root, "startup-build-interrupted", baseline, countBaseline,
						readerBaseline, startupBaseline, true);
				checkChild(root, "startup-build-interrupted", baseline, countBaseline,
						readerBaseline, startupBaseline, false);
			}
			if(!countBaseline){
				checkChild(root, "count-healthy", baseline, countBaseline, readerBaseline, true);
				checkChild(root, "count-healthy", baseline, countBaseline, readerBaseline, false);
				checkChild(root, "count-interrupted", baseline, countBaseline, readerBaseline, true);
				checkChild(root, "count-interrupted", baseline, countBaseline, readerBaseline, false);
			}
			if(!baseline){
				checkChild(root, "healthy", baseline, countBaseline, readerBaseline, true);
				checkChild(root, "healthy", baseline, countBaseline, readerBaseline, false);
				checkChild(root, "interrupted", baseline, countBaseline, readerBaseline, true);
				checkChild(root, "interrupted", baseline, countBaseline, readerBaseline, false);
			}
		}finally{
			deleteTempTree(root);
		}
		System.out.println("BUILD_DISK_TABLES_WORKER_FAILURE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed")+" count="+
				(countBaseline ? "baseline" : "fixed")+" reader="+
				(readerBaseline ? "baseline" : "fixed")+" startup="+
				(startupBaseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private BuildDiskTablesWorkerFailureTest(){}

	/** Launches one bounded child JVM for the requested assertion mode.
	 * @param root Temporary parent directory
	 * @param childCase Child case selector
	 * @param baseline True to expect the current early collector failure
	 * @param countBaseline True to expect the current early CountStats.sum failure
	 * @param readerBaseline True to expect the current worker ByteFile leak
	 * @param assertions True to run the child with -ea
	 * @throws Exception If the child fails or fixture cleanup fails */
	private static void checkChild(final Path root, final String childCase,
			final boolean baseline, final boolean countBaseline,
			final boolean readerBaseline, final boolean assertions) throws Exception{
		checkChild(root, childCase, baseline, countBaseline, readerBaseline,
				baseline, assertions);
	}

	/** Launches one bounded child JVM for the requested assertion mode.
	 * @param root Temporary parent directory
	 * @param childCase Child case selector
	 * @param baseline True to expect the current early collector failure
	 * @param countBaseline True to expect the current early CountStats.sum failure
	 * @param readerBaseline True to expect the current worker ByteFile leak
	 * @param startupBaseline True to expect partial startup to strand a worker
	 * @param assertions True to run the child with -ea
	 * @throws Exception If the child fails or fixture cleanup fails */
	private static void checkChild(final Path root, final String childCase,
			final boolean baseline, final boolean countBaseline,
			final boolean readerBaseline, final boolean startupBaseline,
			final boolean assertions) throws Exception{
		final String label=childCase+"-"+(assertions ? "ea" : "da");
		final Path childRoot=root.resolve(label);
		Files.createDirectories(childRoot);

		final List<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+
				File.separator+"java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx512m");
		command.add("-Xms512m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(BuildDiskTablesWorkerFailureTest.class.getName());
		command.add("child");
		command.add("case="+childCase);
		command.add(baseline ? "expect=baseline" : "expect=fixed");
		command.add(countBaseline ? "count=baseline" : "count=fixed");
		command.add(readerBaseline ? "reader=baseline" : "reader=fixed");
		command.add(startupBaseline ? "startup=baseline" : "startup=fixed");
		command.add("root="+childRoot.toString());

		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		Process child=null;
		final long started=System.nanoTime();
		try{
			child=new ProcessBuilder(command)
					.redirectOutput(stdout.toFile())
					.redirectError(stderr.toFile())
					.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				child.destroy();
				if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
					killChild(child, label, stdout, stderr);
				}
				throw new AssertionError(label+" child timed out\nstdout\n"+
						readAscii(stdout)+"\nstderr\n"+readAscii(stderr));
			}
			final String childOut=readAscii(stdout);
			final String childErr=readAscii(stderr);
			System.out.println(label+" child exit="+child.exitValue()+
					" elapsed_ms="+TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started)+" stdout='"+
					compact(childOut)+"' stderr='"+compact(childErr)+"'");
			check(child.exitValue()==0, label+" child exit "+child.exitValue()+
					"\nstdout\n"+childOut+"\nstderr\n"+childErr);
		}finally{
			if(child!=null && child.isAlive()){
				child.destroy();
				if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
					killChild(child, label, stdout, stderr);
				}
			}
		}
	}

	/** Runs one in-process reflective worker aggregation check.
	 * @param childCase Child case selector
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect waitForThreads to return before the FIFO reader
	 * @param countBaseline True to expect CountStats.sum to return before the FIFO reader
	 * @param readerBaseline True to expect the current worker ByteFile leak
	 * @param startupBaseline True to expect partial startup to strand a worker
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runChild(final String childCase, final Path root,
			final boolean baseline, final boolean countBaseline,
			final boolean readerBaseline, final boolean startupBaseline)
			throws Exception{
		if("failure".equals(childCase)){runFailureChild(root, baseline);}
		else if("reader".equals(childCase)){runReaderChild(root, readerBaseline);}
		else if("count-failure".equals(childCase)){
			runCountFailureChild(root, countBaseline);
		}else if("count-reader".equals(childCase)){
			runCountReaderChild(root, readerBaseline);
		}else if("startup-count".equals(childCase)){
			runStartupChild(root, startupBaseline, true, false);
		}else if("startup-build".equals(childCase)){
			runStartupChild(root, startupBaseline, false, false);
		}else if("startup-count-interrupted".equals(childCase)){
			runStartupChild(root, startupBaseline, true, true);
		}else if("startup-build-interrupted".equals(childCase)){
			runStartupChild(root, startupBaseline, false, true);
		}else if("count-healthy".equals(childCase)){runCountHealthyChild(root);}
		else if("count-interrupted".equals(childCase)){runCountInterruptedChild(root);}
		else if("healthy".equals(childCase)){runHealthyChild(root);}
		else if("interrupted".equals(childCase)){runInterruptedChild(root);}
		else{throw new IllegalArgumentException("Unknown child case: "+childCase);}
	}

	/** Runs one in-process reflective waitForThreads failure check.
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect waitForThreads to return before the FIFO reader
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runFailureChild(final Path root, final boolean baseline) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path poison=root.resolve("poison.tsv");
		final Path fifo=root.resolve("blocking.fifo");
		final Path table=root.resolve("accession-disk.bin");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(poison, poisonInput());
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		DiskAccessionTable.Builder builder=null;
		Thread fifoThread=null, failingThread=null;
		Collector collector=null;
		Throwable primary=null;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			builder=DiskAccessionTable.build(table.toString(), 1);
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(1);
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveProcessThread access=new ReflectiveProcessThread();

			fifoThread=access.create(fifo.toString(), error, buffers);
			fifoThread.start();
			waitForReader(fifo);
			check(countFds(fifo, O_RDONLY)==1, "Expected one worker FIFO reader");
			check(countFds(fifo, O_RDWR)==1, "Expected one O_RDWR FIFO holder");

			failingThread=access.create(poison.toString(), error, buffers);
			failingThread.start();
			final Object threads=Array.newInstance(access.threadClass, 2);
			Array.set(threads, 0, failingThread);
			Array.set(threads, 1, fifoThread);

			collector=new Collector(access.waitForThreads, threads, fifoThread);
			collector.start();
			waitTerminated(failingThread, "failing worker");

			check(countFds(fifo, O_RDONLY)==1,
					"FIFO reader should stay blocked until parent release");
			check(countFds(fifo, O_RDWR)==1,
					"O_RDWR FIFO holder should be counted separately");

			if(baseline){
				collector.join(COLLECTOR_TIMEOUT_MS);
				check(!collector.isAlive(),
						"Baseline waitForThreads did not return before FIFO release");
				check(Boolean.TRUE.equals(collector.fifoAliveAtReturn),
						"Baseline waitForThreads returned after FIFO worker ended: "+
						collector.fifoAliveAtReturn);
			}else{
				waitForCollectorJoin(collector, fifoThread);
			}

			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(), "Collector should finish after FIFO release");
			check(Boolean.valueOf(baseline).equals(collector.fifoAliveAtReturn),
					"waitForThreads returned at wrong FIFO-worker time; baseline="+
					baseline+" fifoAliveAtReturn="+collector.fifoAliveAtReturn);
			checkPrimaryFailure(collector.thrown);
			check(countFds(fifo, O_RDONLY)==0, "FIFO reader leaked after release");
			check(countFds(fifo, O_RDWR)==0, "FIFO holder leaked after release");
			builder.close();
			builder=null;
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(failingThread, "cleanup failing worker", cleanup);
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO worker", cleanup);
			cleanup=cleanupJoin(collector, "cleanup collector", cleanup);
			final boolean live=(failingThread!=null && failingThread.isAlive()) ||
					(fifoThread!=null && fifoThread.isAlive()) ||
					(collector!=null && collector.isAlive());
			if(builder!=null && !live){
				try{builder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}else if(builder!=null){
				cleanup=combine(cleanup, new AssertionError(
						"Cannot close builder while an owned thread is still active"));
			}
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_WORKER_FAILURE_CHILD_PASS expect="+
				(baseline ? "baseline" : "fixed")+
				" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Runs one second-pass ByteFile cleanup check without FIFO blocking.
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect the current worker read-descriptor leak
	 * @throws Exception If reflection or worker cleanup fails */
	private static void runReaderChild(final Path root, final boolean baseline)
			throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path poison=root.resolve("poison.tsv");
		final Path table=root.resolve("accession-disk.bin");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(poison, poisonInput());

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		DiskAccessionTable.Builder builder=null;
		Thread worker=null;
		Throwable primary=null;
		int readFds=-1;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			builder=DiskAccessionTable.build(table.toString(), 1);
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(1);
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveProcessThread access=new ReflectiveProcessThread();

			worker=access.create(poison.toString(), error, buffers);
			worker.start();
			waitTerminated(worker, "poison worker");

			final Object threads=Array.newInstance(access.threadClass, 1);
			Array.set(threads, 0, worker);
			final Throwable failure=aggregationFailure(access.waitForThreads, threads);
			checkPrimaryFailure(failure);
			readFds=countFds(poison, O_RDONLY);
			final int expected=baseline ? 1 : 0;
			check(readFds==expected, "Second-pass poison reader FDs should be "+
					expected+", not "+readFds);
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=cleanupJoin(worker, "cleanup poison worker", null);
			if(builder!=null && (worker==null || !worker.isAlive())){
				try{builder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}else if(builder!=null){
				cleanup=combine(cleanup, new AssertionError(
						"Cannot close builder while an owned thread is still active"));
			}
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_READER_CHILD_PASS mode=second expect="+
				(baseline ? "baseline" : "fixed")+" read_fds="+readFds+
				" error=ArrayIndexOutOfBoundsException");
	}

	/** Runs one in-process reflective CountStats.sum failure check.
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect CountStats.sum to return before the FIFO reader
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runCountFailureChild(final Path root, final boolean baseline)
			throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path poison=root.resolve("poison.tsv");
		final Path fifo=root.resolve("blocking.fifo");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(poison, poisonInput());
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		Thread fifoThread=null, failingThread=null;
		Collector collector=null;
		Throwable primary=null;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveCountThread access=new ReflectiveCountThread();

			fifoThread=access.create(fifo.toString(), error);
			fifoThread.start();
			waitForReader(fifo);
			check(countFds(fifo, O_RDONLY)==1,
					"Expected one first-pass FIFO reader");
			check(countFds(fifo, O_RDWR)==1,
					"Expected one O_RDWR first-pass FIFO holder");

			failingThread=access.create(poison.toString(), error);
			failingThread.start();
			final Object threads=Array.newInstance(access.threadClass, 2);
			Array.set(threads, 0, failingThread);
			Array.set(threads, 1, fifoThread);

			collector=new Collector(access.sum, threads, fifoThread);
			collector.start();
			waitTerminated(failingThread, "failing count worker");

			check(countFds(fifo, O_RDONLY)==1,
					"First-pass FIFO reader should stay blocked until parent release");
			check(countFds(fifo, O_RDWR)==1,
					"First-pass O_RDWR FIFO holder should be counted separately");

			if(baseline){
				collector.join(COLLECTOR_TIMEOUT_MS);
				check(!collector.isAlive(),
						"Baseline CountStats.sum did not return before FIFO release");
				check(Boolean.TRUE.equals(collector.fifoAliveAtReturn),
						"Baseline CountStats.sum returned after FIFO worker ended: "+
						collector.fifoAliveAtReturn);
			}else{
				waitForCollectorJoin(collector, fifoThread);
			}

			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO count worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(),
					"CountStats.sum collector should finish after FIFO release");
			check(Boolean.valueOf(baseline).equals(collector.fifoAliveAtReturn),
					"CountStats.sum returned at wrong FIFO-worker time; baseline="+
					baseline+" fifoAliveAtReturn="+collector.fifoAliveAtReturn);
			checkPrimaryFailure(collector.thrown);
			check(countFds(fifo, O_RDONLY)==0,
					"First-pass FIFO reader leaked after release");
			check(countFds(fifo, O_RDWR)==0,
					"First-pass FIFO holder leaked after release");
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(failingThread, "cleanup failing count worker", cleanup);
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO count worker", cleanup);
			cleanup=cleanupJoin(collector, "cleanup count collector", cleanup);
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_COUNT_FAILURE_CHILD_PASS expect="+
				(baseline ? "baseline" : "fixed")+
				" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Runs one first-pass ByteFile cleanup check without FIFO blocking.
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect the current worker read-descriptor leak
	 * @throws Exception If reflection or worker cleanup fails */
	private static void runCountReaderChild(final Path root, final boolean baseline)
			throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path poison=root.resolve("poison.tsv");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(poison, poisonInput());

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		Thread worker=null;
		Throwable primary=null;
		int readFds=-1;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveCountThread access=new ReflectiveCountThread();

			worker=access.create(poison.toString(), error);
			worker.start();
			waitTerminated(worker, "poison count worker");

			final Object threads=Array.newInstance(access.threadClass, 1);
			Array.set(threads, 0, worker);
			final Throwable failure=aggregationFailure(access.sum, threads);
			checkPrimaryFailure(failure);
			readFds=countFds(poison, O_RDONLY);
			final int expected=baseline ? 1 : 0;
			check(readFds==expected, "First-pass poison reader FDs should be "+
					expected+", not "+readFds);
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=cleanupJoin(worker, "cleanup poison count worker", null);
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_READER_CHILD_PASS mode=count expect="+
				(baseline ? "baseline" : "fixed")+" read_fds="+readFds+
				" error=ArrayIndexOutOfBoundsException");
	}

	/** Runs one actual startup-loop failure check.
	 * @param root Temporary child artifact directory
	 * @param baseline True to expect the startup failure to strand the first worker
	 * @param count True for first-pass startup, false for second-pass startup
	 * @param interrupted True to interrupt the startup caller during fixed cleanup
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runStartupChild(final Path root, final boolean baseline,
			final boolean count, final boolean interrupted) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path healthy=root.resolve("healthy.tsv");
		final Path fifo=root.resolve("blocking.fifo");
		final Path table=root.resolve("accession-disk.bin");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(healthy, HEALTHY.getBytes(StandardCharsets.US_ASCII));
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		final int oldThreads=Shared.threads();
		RandomAccessFile holder=null;
		DiskAccessionTable.Builder builder=null;
		StartupBuildDiskTables starter=null;
		StartupCaller caller=null;
		Throwable primary=null;
		try{
			check(!interrupted || !baseline,
					"Interrupted startup cleanup is fixed-only evidence");
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(2);
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			starter=new StartupBuildDiskTables(fifo);
			final String[] files={fifo.toString(), healthy.toString()};
			final Method spawn;
			final Object[] args;
			if(count){
				spawn=BuildDiskTables.class.getDeclaredMethod("spawnCountThreads",
						String[].class);
				args=new Object[]{files};
			}else{
				builder=DiskAccessionTable.build(table.toString(), 1);
				final DiskAccessionTable.BufferSet buffers=builder.bufferSet(1);
				spawn=BuildDiskTables.class.getDeclaredMethod("spawnThreads",
						String[].class, DiskAccessionTable.BufferSet.class);
				args=new Object[]{files, buffers};
			}
			spawn.setAccessible(true);

			caller=new StartupCaller(spawn, starter, args);
			caller.start();
			waitForStartupInjection(starter, caller);
			check(countFds(fifo, O_RDONLY)==1,
					"Startup FIFO worker should hold one read-only descriptor");
			check(countFds(fifo, O_RDWR)==1,
					"Startup FIFO holder should be counted separately");

			if(baseline){
				caller.join(COLLECTOR_TIMEOUT_MS);
				check(!caller.isAlive(),
						"Baseline startup caller did not return before FIFO release");
				check(Boolean.TRUE.equals(caller.fifoAliveAtReturn),
						"Baseline startup caller returned after the FIFO worker ended: "+
						caller.fifoAliveAtReturn);
			}else{
				waitForStartupJoin(caller, starter.firstWorker);
				if(interrupted){
					caller.interrupt();
					waitForInterruptedStartupRejoin(caller, starter.firstWorker);
				}
			}

			holder.close();
			holder=null;
			waitTerminated(starter.firstWorker, "startup FIFO worker");
			caller.join(COLLECTOR_TIMEOUT_MS);
			check(!caller.isAlive(), "Startup caller should finish after FIFO release");
			check(caller.thrown==starter.sentinel,
					"Startup caller should preserve the injected sentinel, not "+
					caller.thrown);
			check(Boolean.valueOf(baseline).equals(caller.fifoAliveAtReturn),
					"Startup caller returned at wrong FIFO-worker time; baseline="+
					baseline+" fifoAliveAtReturn="+caller.fifoAliveAtReturn);
			if(interrupted){
				check(Boolean.TRUE.equals(caller.interruptedAtReturn),
						"Startup caller did not restore its interrupt flag: "+
						caller.interruptedAtReturn);
				check(hasSuppressedInterrupted(caller.thrown),
						"Startup sentinel did not retain the interrupted join");
			}
			check(countFds(fifo, O_RDONLY)==0,
					"Startup FIFO reader leaked after release");
			check(countFds(fifo, O_RDWR)==0,
					"Startup FIFO holder leaked after release");
			if(builder!=null){
				builder.close();
				builder=null;
			}
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(starter==null ? null : starter.firstWorker,
					"cleanup startup FIFO worker", cleanup);
			cleanup=cleanupJoin(caller, "cleanup startup caller", cleanup);
			final boolean live=(starter!=null && starter.firstWorker!=null &&
					starter.firstWorker.isAlive()) || (caller!=null && caller.isAlive());
			if(builder!=null && !live){
				try{builder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}else if(builder!=null){
				cleanup=combine(cleanup, new AssertionError(
						"Cannot close builder while an owned startup thread is still active"));
			}
			Shared.setThreads(oldThreads);
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		if(interrupted){
			System.out.println("BUILD_DISK_TABLES_STARTUP_INTERRUPTED_CHILD_PASS mode="+
					(count ? "count" : "build")+" fifo_alive_at_return="+
					caller.fifoAliveAtReturn+" interrupted_at_return="+
					caller.interruptedAtReturn);
		}else{
			System.out.println("BUILD_DISK_TABLES_STARTUP_CHILD_PASS mode="+
					(count ? "count" : "build")+" expect="+
					(baseline ? "baseline" : "fixed")+" fifo_alive_at_return="+
					caller.fifoAliveAtReturn);
		}
	}

	/** Runs one fixed-mode successful CountStats.sum check.
	 * @param root Temporary child artifact directory
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runCountHealthyChild(final Path root) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path healthy=root.resolve("healthy.tsv");
		final Path fifo=root.resolve("blocking.fifo");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(healthy, HEALTHY.getBytes(StandardCharsets.US_ASCII));
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		Thread healthyThread=null, fifoThread=null;
		Collector collector=null;
		Throwable primary=null;
		long total=-1, digitizable=-1, nonDigitizable=-1, estimatedUnique=-1;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveCountThread access=new ReflectiveCountThread();

			healthyThread=access.create(healthy.toString(), error);
			fifoThread=access.create(fifo.toString(), error);
			fifoThread.start();
			waitForReader(fifo);
			healthyThread.start();

			final Object threads=Array.newInstance(access.threadClass, 2);
			Array.set(threads, 0, healthyThread);
			Array.set(threads, 1, fifoThread);
			collector=new Collector(access.sum, threads, fifoThread);
			collector.start();

			waitTerminated(healthyThread, "healthy count worker");
			waitForCollectorJoin(collector, fifoThread);
			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO count worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(),
					"Healthy CountStats.sum collector should finish after FIFO release");
			check(collector.thrown==null,
					"Healthy CountStats.sum collector should not throw: "+
					collector.thrown);
			check(Boolean.FALSE.equals(collector.fifoAliveAtReturn),
					"Healthy CountStats.sum returned before FIFO worker ended: "+
					collector.fifoAliveAtReturn);
			total=access.total(collector.result);
			digitizable=access.digitizable(collector.result);
			nonDigitizable=access.nonDigitizable(collector.result);
			estimatedUnique=access.estimatedUnique(collector.result);
			check(total==1, "Healthy CountStats total should be 1, not "+total);
			check(digitizable==1,
					"Healthy CountStats digitizable should be 1, not "+digitizable);
			check(nonDigitizable==0,
					"Healthy CountStats nonDigitizable should be 0, not "+
					nonDigitizable);
			check(estimatedUnique==1,
					"Healthy CountStats estimatedUnique should be 1, not "+
					estimatedUnique);
			check(countFds(fifo, O_RDONLY)==0,
					"First-pass FIFO reader leaked after healthy release");
			check(countFds(fifo, O_RDWR)==0,
					"First-pass FIFO holder leaked after healthy release");
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(healthyThread,
					"cleanup healthy count worker", cleanup);
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO count worker", cleanup);
			cleanup=cleanupJoin(collector, "cleanup healthy count collector", cleanup);
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_COUNT_HEALTHY_CHILD_PASS total="+total+
				" digitizable="+digitizable+" non_digitizable="+nonDigitizable+
				" estimated_unique="+estimatedUnique+
				" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Runs one fixed-mode interrupted CountStats.sum check.
	 * @param root Temporary child artifact directory
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runCountInterruptedChild(final Path root) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path fifo=root.resolve("blocking.fifo");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		Thread fifoThread=null;
		Collector collector=null;
		Throwable primary=null;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveCountThread access=new ReflectiveCountThread();

			fifoThread=access.create(fifo.toString(), error);
			fifoThread.start();
			waitForReader(fifo);

			final Object threads=Array.newInstance(access.threadClass, 1);
			Array.set(threads, 0, fifoThread);
			collector=new Collector(access.sum, threads, fifoThread);
			collector.start();
			waitForCollectorJoin(collector, fifoThread);
			collector.interrupt();
			waitForInterruptedCollectorRejoin(collector, fifoThread);

			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO count worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(),
					"Interrupted CountStats.sum collector should finish after FIFO release");
			check(collector.thrown instanceof RuntimeException,
					"Interrupted CountStats.sum collector should throw RuntimeException, not "+
					collector.thrown);
			check(collector.thrown.getCause()!=null &&
					collector.thrown.getCause().getClass()==InterruptedException.class,
					"Interrupted CountStats.sum collector should preserve exact "+
					"InterruptedException cause, not "+collector.thrown.getCause());
			check(Boolean.FALSE.equals(collector.fifoAliveAtReturn),
					"Interrupted CountStats.sum returned before FIFO worker ended: "+
					collector.fifoAliveAtReturn);
			check(Boolean.TRUE.equals(collector.interruptedAtReturn),
					"Interrupted CountStats.sum did not restore its interrupt flag: "+
					collector.interruptedAtReturn);
			check(countFds(fifo, O_RDONLY)==0,
					"First-pass FIFO reader leaked after interrupt release");
			check(countFds(fifo, O_RDWR)==0,
					"First-pass FIFO holder leaked after interrupt release");
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO count worker", cleanup);
			cleanup=cleanupJoin(collector,
					"cleanup interrupted count collector", cleanup);
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_COUNT_INTERRUPTED_CHILD_PASS "+
				"interrupted_at_return="+collector.interruptedAtReturn+
				" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Runs one fixed-mode no-error waitForThreads check.
	 * @param root Temporary child artifact directory
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runHealthyChild(final Path root) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path healthy=root.resolve("healthy.tsv");
		final Path fifo=root.resolve("blocking.fifo");
		final Path table=root.resolve("accession-disk.bin");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		Files.write(healthy, HEALTHY.getBytes(StandardCharsets.US_ASCII));
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		DiskAccessionTable.Builder builder=null;
		Thread healthyThread=null, fifoThread=null;
		Collector collector=null;
		Throwable primary=null;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			builder=DiskAccessionTable.build(table.toString(), 1);
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(1);
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveProcessThread access=new ReflectiveProcessThread();

			healthyThread=access.create(healthy.toString(), error, buffers);
			fifoThread=access.create(fifo.toString(), error, buffers);
			fifoThread.start();
			waitForReader(fifo);
			healthyThread.start();

			final Object threads=Array.newInstance(access.threadClass, 2);
			Array.set(threads, 0, healthyThread);
			Array.set(threads, 1, fifoThread);
			collector=new Collector(access.waitForThreads, threads, fifoThread);
			collector.start();

			waitTerminated(healthyThread, "healthy worker");
			waitForCollectorJoin(collector, fifoThread);
			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(), "Healthy collector should finish after FIFO release");
			check(collector.thrown==null, "Healthy collector should not throw: "+collector.thrown);
			check(Boolean.FALSE.equals(collector.fifoAliveAtReturn),
					"Healthy collector returned before FIFO worker ended: "+
					collector.fifoAliveAtReturn);
			check(Long.valueOf(1).equals(collector.collected),
					"Healthy collector should return one collected row, not "+
					collector.collected);
			buffers.flush();
			check(buffers.added()==1,
					"Healthy buffers should contain one row, not "+buffers.added());
			check(countFds(fifo, O_RDONLY)==0, "FIFO reader leaked after healthy release");
			check(countFds(fifo, O_RDWR)==0, "FIFO holder leaked after healthy release");
			builder.close();
			builder=null;
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(healthyThread, "cleanup healthy worker", cleanup);
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO worker", cleanup);
			cleanup=cleanupJoin(collector, "cleanup collector", cleanup);
			final boolean live=(healthyThread!=null && healthyThread.isAlive()) ||
					(fifoThread!=null && fifoThread.isAlive()) ||
					(collector!=null && collector.isAlive());
			if(builder!=null && !live){
				try{builder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}else if(builder!=null){
				cleanup=combine(cleanup, new AssertionError(
						"Cannot close builder while an owned thread is still active"));
			}
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_WORKER_HEALTHY_CHILD_PASS collected="+
				collector.collected+" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Runs one fixed-mode interrupted-collector waitForThreads check.
	 * @param root Temporary child artifact directory
	 * @throws Exception If reflection, FIFO setup or worker cleanup fails */
	private static void runInterruptedChild(final Path root) throws Exception{
		final Path pattern=root.resolve("patterns.tsv");
		final Path fifo=root.resolve("blocking.fifo");
		final Path table=root.resolve("accession-disk.bin");

		Files.write(pattern, PATTERNS.getBytes(StandardCharsets.US_ASCII));
		mkfifo(fifo);

		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		final boolean bf3=ByteFile.FORCE_MODE_BF3, bf4=ByteFile.FORCE_MODE_BF4;
		RandomAccessFile holder=null;
		DiskAccessionTable.Builder builder=null;
		Thread fifoThread=null;
		Collector collector=null;
		Throwable primary=null;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			AnalyzeAccession.loadCodeMap(pattern.toString());

			holder=new RandomAccessFile(fifo.toFile(), "rw");
			builder=DiskAccessionTable.build(table.toString(), 1);
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(1);
			final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
			final ReflectiveProcessThread access=new ReflectiveProcessThread();

			fifoThread=access.create(fifo.toString(), error, buffers);
			fifoThread.start();
			waitForReader(fifo);

			final Object threads=Array.newInstance(access.threadClass, 1);
			Array.set(threads, 0, fifoThread);
			collector=new Collector(access.waitForThreads, threads, fifoThread);
			collector.start();
			waitForCollectorJoin(collector, fifoThread);
			collector.interrupt();
			waitForInterruptedCollectorRejoin(collector, fifoThread);

			holder.close();
			holder=null;
			waitTerminated(fifoThread, "FIFO worker");
			collector.join(COLLECTOR_TIMEOUT_MS);
			check(!collector.isAlive(),
					"Interrupted collector should finish after FIFO release");
			check(collector.thrown instanceof RuntimeException,
					"Interrupted collector should throw RuntimeException, not "+
					collector.thrown);
			check(collector.thrown.getCause()!=null &&
					collector.thrown.getCause().getClass()==InterruptedException.class,
					"Interrupted collector should preserve exact InterruptedException "+
					"cause, not "+collector.thrown.getCause());
			check(Boolean.FALSE.equals(collector.fifoAliveAtReturn),
					"Interrupted collector returned before FIFO worker ended: "+
					collector.fifoAliveAtReturn);
			check(Boolean.TRUE.equals(collector.interruptedAtReturn),
					"Interrupted collector did not restore its interrupt flag: "+
					collector.interruptedAtReturn);
			check(countFds(fifo, O_RDONLY)==0, "FIFO reader leaked after interrupt release");
			check(countFds(fifo, O_RDWR)==0, "FIFO holder leaked after interrupt release");
			builder.close();
			builder=null;
		}catch(Throwable t){
			primary=t;
			rethrow(t);
		}finally{
			Throwable cleanup=null;
			if(holder!=null){
				try{holder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}
			cleanup=cleanupJoin(fifoThread, "cleanup FIFO worker", cleanup);
			cleanup=cleanupJoin(collector, "cleanup collector", cleanup);
			final boolean live=(fifoThread!=null && fifoThread.isAlive()) ||
					(collector!=null && collector.isAlive());
			if(builder!=null && !live){
				try{builder.close();}
				catch(Throwable t){cleanup=combine(cleanup, t);}
			}else if(builder!=null){
				cleanup=combine(cleanup, new AssertionError(
						"Cannot close builder while an owned thread is still active"));
			}
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			if(cleanup!=null){
				if(primary==null){rethrow(cleanup);}
				else if(primary!=cleanup){primary.addSuppressed(cleanup);}
			}
		}
		System.out.println("BUILD_DISK_TABLES_WORKER_INTERRUPTED_CHILD_PASS "+
				"interrupted_at_return="+collector.interruptedAtReturn+
				" fifo_alive_at_return="+collector.fifoAliveAtReturn);
	}

	/** Creates one reflected BuildDiskTables.ProcessThread factory. */
	private static final class ReflectiveProcessThread{

		/** Locates private BuildDiskTables worker entry points.
		 * @throws Exception If reflection fails */
		ReflectiveProcessThread() throws Exception{
			threadClass=Class.forName("tax.BuildDiskTables$ProcessThread");
			constructor=threadClass.getDeclaredConstructor(String[].class,
					AtomicInteger.class, AtomicReference.class,
					DiskAccessionTable.BufferSet.class);
			constructor.setAccessible(true);
			waitForThreads=BuildDiskTables.class.getDeclaredMethod("waitForThreads",
					Array.newInstance(threadClass, 0).getClass());
			waitForThreads.setAccessible(true);
		}

		/** Builds one worker that owns a single input file.
		 * @param file Input path for the worker
		 * @param error Shared worker-failure slot
		 * @param buffers Shared accession table buffers
		 * @return New unstarted worker thread
		 * @throws Exception If reflection fails */
		Thread create(final String file, final AtomicReference<Throwable> error,
				final DiskAccessionTable.BufferSet buffers) throws Exception{
			return (Thread)constructor.newInstance(new String[]{file},
					new AtomicInteger(0), error, buffers);
		}

		/** Private BuildDiskTables.ProcessThread class. */
		final Class<?> threadClass;
		/** Private worker constructor. */
		final Constructor<?> constructor;
		/** Private waitForThreads(ProcessThread[]) method. */
		final Method waitForThreads;

	}

	/** Creates one reflected BuildDiskTables.CountThread factory. */
	private static final class ReflectiveCountThread{

		/** Locates private BuildDiskTables first-pass worker entry points.
		 * @throws Exception If reflection fails */
		ReflectiveCountThread() throws Exception{
			threadClass=Class.forName("tax.BuildDiskTables$CountThread");
			final Class<?> statsClass=Class.forName("tax.BuildDiskTables$CountStats");
			constructor=threadClass.getDeclaredConstructor(String[].class,
					AtomicInteger.class, AtomicReference.class);
			constructor.setAccessible(true);
			sum=statsClass.getDeclaredMethod("sum",
					Array.newInstance(threadClass, 0).getClass());
			sum.setAccessible(true);
			totalField=statsClass.getDeclaredField("total");
			digitizableField=statsClass.getDeclaredField("digitizable");
			nonDigitizableField=statsClass.getDeclaredField("nonDigitizable");
			estimatedUniqueField=statsClass.getDeclaredField("estimatedUnique");
			totalField.setAccessible(true);
			digitizableField.setAccessible(true);
			nonDigitizableField.setAccessible(true);
			estimatedUniqueField.setAccessible(true);
		}

		/** Builds one worker that owns a single input file.
		 * @param file Input path for the worker
		 * @param error Shared worker-failure slot
		 * @return New unstarted worker thread
		 * @throws Exception If reflection fails */
		Thread create(final String file, final AtomicReference<Throwable> error)
				throws Exception{
			return (Thread)constructor.newInstance(new String[]{file},
					new AtomicInteger(0), error);
		}

		/** Reads a CountStats total field.
		 * @param stats Reflected CountStats result
		 * @return Total source-row count
		 * @throws Exception If reflection fails */
		long total(final Object stats) throws Exception{return totalField.getLong(stats);}

		/** Reads a CountStats digitizable field.
		 * @param stats Reflected CountStats result
		 * @return Digitizable row count
		 * @throws Exception If reflection fails */
		long digitizable(final Object stats) throws Exception{
			return digitizableField.getLong(stats);
		}

		/** Reads a CountStats nonDigitizable field.
		 * @param stats Reflected CountStats result
		 * @return Non-digitizable row count
		 * @throws Exception If reflection fails */
		long nonDigitizable(final Object stats) throws Exception{
			return nonDigitizableField.getLong(stats);
		}

		/** Reads a CountStats estimatedUnique field.
		 * @param stats Reflected CountStats result
		 * @return Estimated unique-key count
		 * @throws Exception If reflection fails */
		long estimatedUnique(final Object stats) throws Exception{
			return estimatedUniqueField.getLong(stats);
		}

		/** Private BuildDiskTables.CountThread class. */
		final Class<?> threadClass;
		/** Private worker constructor. */
		final Constructor<?> constructor;
		/** Private CountStats.sum(CountThread[]) method. */
		final Method sum;
		/** Private CountStats.total field. */
		final Field totalField;
		/** Private CountStats.digitizable field. */
		final Field digitizableField;
		/** Private CountStats.nonDigitizable field. */
		final Field nonDigitizableField;
		/** Private CountStats.estimatedUnique field. */
		final Field estimatedUniqueField;

	}

	/** BuildDiskTables subclass that injects a failure on the second production start. */
	private static final class StartupBuildDiskTables extends BuildDiskTables{

		/** Creates one startup-failure harness.
		 * @param fifo_ FIFO the first worker should block on */
		StartupBuildDiskTables(final Path fifo_){
			super(new String[]{"tree=f"});
			fifo=fifo_;
		}

		@Override
		protected void startWorker(final Thread worker){
			final int attempt=attempts.incrementAndGet();
			if(attempt==1){
				firstWorker=worker;
				worker.start();
			}else if(attempt==2){
				try{
					waitForReader(fifo);
					check(countFds(fifo, O_RDONLY)==1,
							"Startup injection needs one worker FIFO reader");
					check(countFds(fifo, O_RDWR)==1,
							"Startup injection should ignore the O_RDWR FIFO holder");
				}catch(Exception e){
					throw new RuntimeException(e);
				}
				secondAttempted=true;
				throw sentinel;
			}else{
				throw new AssertionError("Unexpected worker start attempt "+attempt);
			}
		}

		/** FIFO path opened by the first worker. */
		private final Path fifo;
		/** Sentinel thrown from the second worker start attempt. */
		private final RuntimeException sentinel=
				new RuntimeException("Injected BuildDiskTables startup failure");
		/** Number of production startWorker calls. */
		private final AtomicInteger attempts=new AtomicInteger(0);
		/** First production worker thread, started on the FIFO path. */
		volatile Thread firstWorker=null;
		/** True after the second start attempt has observed the FIFO reader. */
		volatile boolean secondAttempted=false;

	}

	/** Invokes an actual private BuildDiskTables startup method on a separate thread. */
	private static final class StartupCaller extends Thread{

		/** Creates one startup caller.
		 * @param spawn Private production spawn method
		 * @param target BuildDiskTables instance to invoke on
		 * @param args Spawn method arguments */
		StartupCaller(final Method spawn, final StartupBuildDiskTables target,
				final Object[] args){
			method=spawn;
			receiver=target;
			arguments=args;
		}

		@Override
		public void run(){
			try{
				result=method.invoke(receiver, arguments);
			}catch(InvocationTargetException e){
				thrown=e.getCause();
			}catch(Throwable e){
				thrown=e;
			}finally{
				final Thread first=receiver.firstWorker;
				fifoAliveAtReturn=Boolean.valueOf(first!=null && first.isAlive());
				interruptedAtReturn=Boolean.valueOf(Thread.currentThread().isInterrupted());
			}
		}

		/** Reflected production startup method. */
		private final Method method;
		/** BuildDiskTables instance receiving the reflected call. */
		private final StartupBuildDiskTables receiver;
		/** Arguments for the reflected startup method. */
		private final Object[] arguments;
		/** Returned worker array, or null when startup throws. */
		volatile Object result=null;
		/** Startup failure thrown by the reflected call. */
		volatile Throwable thrown=null;
		/** Whether the first worker was alive as startup returned. */
		volatile Boolean fifoAliveAtReturn=null;
		/** Whether startup restored interruption before returning. */
		volatile Boolean interruptedAtReturn=null;

	}

	/** Runs a reflected worker aggregation on a separate thread. */
	private static final class Collector extends Thread{

		/** Creates one collector.
		 * @param method Reflected worker aggregation method
		 * @param threads Worker array instance
		 * @param fifoWorker Worker that should block in its FIFO read
		 */
		Collector(final Method method, final Object threads, final Thread fifoWorker){
			aggregate=method;
			workerThreads=threads;
			fifoThread=fifoWorker;
		}

		@Override
		public void run(){
			try{
				result=aggregate.invoke(null, workerThreads);
				if(result instanceof Long){collected=(Long)result;}
			}catch(InvocationTargetException e){
				thrown=e.getCause();
			}catch(Throwable e){
				thrown=e;
			}finally{
				fifoAliveAtReturn=Boolean.valueOf(fifoThread.isAlive());
				interruptedAtReturn=Boolean.valueOf(Thread.currentThread().isInterrupted());
			}
		}

		/** Reflected worker aggregation method. */
		private final Method aggregate;
		/** Worker array instance. */
		private final Object workerThreads;
		/** FIFO worker observed at aggregation return. */
		private final Thread fifoThread;
		/** Returned reflected result, or null when it throws. */
		volatile Object result=null;
		/** Returned waitForThreads row count, or null when it throws. */
		volatile Long collected=null;
		/** Primary failure thrown by the aggregation. */
		volatile Throwable thrown=null;
		/** Whether the FIFO worker was alive as the aggregation returned. */
		volatile Boolean fifoAliveAtReturn=null;
		/** Whether the aggregation restored interruption before returning. */
		volatile Boolean interruptedAtReturn=null;

	}

	/** Creates a POSIX FIFO.
	 * @param path FIFO path to create
	 * @throws Exception If mkfifo fails */
	private static void mkfifo(final Path path) throws Exception{
		final Process process=new ProcessBuilder("mkfifo", path.toString())
				.redirectErrorStream(true).start();
		final byte[] output=readAll(process);
		check(process.waitFor()==0, "mkfifo failed: "+
				new String(output, StandardCharsets.UTF_8));
	}

	/** Waits until the FIFO worker has opened its own read descriptor.
	 * @param fifo FIFO path
	 * @throws Exception If the reader does not appear */
	private static void waitForReader(final Path fifo) throws Exception{
		final long deadline=System.nanoTime()+
				TimeUnit.MILLISECONDS.toNanos(FD_TIMEOUT_MS);
		while(System.nanoTime()<deadline){
			if(countFds(fifo, O_RDONLY)==1){return;}
			Thread.sleep(POLL_MS);
		}
		throw new AssertionError("FIFO reader did not appear; readFds="+
				countFds(fifo, O_RDONLY)+" holders="+countFds(fifo, O_RDWR));
	}

	/** Observes the collector waiting on the live FIFO worker's Thread.join monitor.
	 * This gates release for the existing join-based collector so release cannot
	 * race ahead of its invocation and hide an early-return defect.
	 * @param collector Thread invoking the reflected worker aggregation
	 * @param worker FIFO worker that must remain live until release
	 * @throws Exception If the collector returns early or never reaches its join */
	private static void waitForCollectorJoin(final Thread collector, final Thread worker)
			throws Exception{
		waitForWorkerJoin(collector, worker, "Collector");
	}

	/** Observes the startup caller waiting on the live FIFO worker's join monitor.
	 * @param caller Thread invoking the reflected production startup method
	 * @param worker FIFO worker that must remain live until release
	 * @throws Exception If the caller returns or never reaches its join */
	private static void waitForStartupJoin(final Thread caller, final Thread worker)
			throws Exception{
		waitForWorkerJoin(caller, worker, "Startup caller");
	}

	/** Waits until an interrupted startup caller consumes the interrupt and rejoins.
	 * @param caller Thread invoking the reflected production startup method
	 * @param worker FIFO worker that must remain live until release
	 * @throws Exception If the caller exits or does not rejoin the worker */
	private static void waitForInterruptedStartupRejoin(final Thread caller,
			final Thread worker) throws Exception{
		final ThreadMXBean threads=ManagementFactory.getThreadMXBean();
		final int monitor=System.identityHashCode(worker);
		final long deadline=System.nanoTime()+
				TimeUnit.MILLISECONDS.toNanos(COLLECTOR_TIMEOUT_MS);
		while(System.nanoTime()<deadline){
			check(caller.isAlive(),
					"Interrupted startup caller returned before FIFO release");
			if(!caller.isInterrupted()){
				final ThreadInfo info=threads.getThreadInfo(caller.getId());
				if(isWaitingOn(info, monitor, worker)){return;}
			}
			Thread.sleep(POLL_MS);
		}
		throw new AssertionError(
				"Interrupted startup caller did not rejoin the FIFO worker");
	}

	/** Observes one thread waiting on a live FIFO worker's Thread.join monitor.
	 * @param waiter Thread expected to join the worker
	 * @param worker FIFO worker that must remain live until release
	 * @param label Diagnostic label for the waiter thread
	 * @throws Exception If the waiter exits or never reaches its join */
	private static void waitForWorkerJoin(final Thread waiter, final Thread worker,
			final String label) throws Exception{
		final ThreadMXBean threads=ManagementFactory.getThreadMXBean();
		final int monitor=System.identityHashCode(worker);
		final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(COLLECTOR_TIMEOUT_MS);
		while(System.nanoTime()<deadline){
			check(waiter.isAlive(), label+" returned before the FIFO worker was released");
			final ThreadInfo info=threads.getThreadInfo(waiter.getId());
			if(isWaitingOn(info, monitor, worker)){return;}
			Thread.sleep(POLL_MS);
		}
		throw new AssertionError(label+" did not join the blocked FIFO worker");
	}

	/** Waits until the injected second start attempt has reached the FIFO gate.
	 * @param starter Startup-injection BuildDiskTables instance
	 * @param caller Thread invoking the reflected production startup method
	 * @throws Exception If the second attempt does not appear */
	private static void waitForStartupInjection(final StartupBuildDiskTables starter,
			final StartupCaller caller) throws Exception{
		final long deadline=System.nanoTime()+
				TimeUnit.MILLISECONDS.toNanos(COLLECTOR_TIMEOUT_MS);
		while(System.nanoTime()<deadline){
			if(starter.secondAttempted){return;}
			if(!caller.isAlive()){
				if(starter.secondAttempted){return;}
				throw new AssertionError("Startup caller returned before the gated "+
						"second start attempt: "+caller.thrown);
			}
			Thread.sleep(POLL_MS);
		}
		throw new AssertionError("Startup failure injection did not reach the "+
				"second worker start; attempts="+starter.attempts.get());
	}

	/** Waits until an interrupted collector consumes the interrupt and rejoins.
	 * @param collector Thread invoking the reflected worker aggregation
	 * @param worker FIFO worker that must remain live until release
	 * @throws Exception If the collector exits or does not rejoin the worker */
	private static void waitForInterruptedCollectorRejoin(final Thread collector,
			final Thread worker) throws Exception{
		final ThreadMXBean threads=ManagementFactory.getThreadMXBean();
		final int monitor=System.identityHashCode(worker);
		final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(COLLECTOR_TIMEOUT_MS);
		while(System.nanoTime()<deadline){
			check(collector.isAlive(), "Interrupted collector returned before FIFO release");
			if(!collector.isInterrupted()){
				final ThreadInfo info=threads.getThreadInfo(collector.getId());
				if(isWaitingOn(info, monitor, worker)){return;}
			}
			Thread.sleep(POLL_MS);
		}
		throw new AssertionError("Interrupted collector did not rejoin the FIFO worker");
	}

	/** Reports whether a thread is WAITING on a worker's join monitor.
	 * @param info Thread information to inspect
	 * @param monitor Expected monitor identity hash
	 * @param worker Worker whose monitor should be held
	 * @return True only for a Thread.join wait on the worker */
	private static boolean isWaitingOn(final ThreadInfo info, final int monitor,
			final Thread worker){
		if(info==null || info.getThreadState()!=Thread.State.WAITING){return false;}
		final LockInfo lock=info.getLockInfo();
		return lock!=null && lock.getIdentityHashCode()==monitor &&
				lock.getClassName().equals(worker.getClass().getName());
	}

	/** Reports whether a throwable suppressed an InterruptedException.
	 * @param thrown Throwable to inspect
	 * @return True if an InterruptedException was suppressed */
	private static boolean hasSuppressedInterrupted(final Throwable thrown){
		for(Throwable suppressed : thrown.getSuppressed()){
			if(suppressed instanceof InterruptedException){return true;}
		}
		return false;
	}

	/** Waits for a worker to terminate.
	 * @param thread Worker thread
	 * @param label Diagnostic thread label
	 * @throws Exception If the thread stays live */
	private static void waitTerminated(final Thread thread, final String label)
			throws Exception{
		thread.join(WORKER_TIMEOUT_MS);
		check(!thread.isAlive(), label+" did not terminate");
	}

	/** Requires the preserved first worker failure.
	 * @param thrown Throwable emitted by the worker aggregation
	 */
	private static void checkPrimaryFailure(final Throwable thrown){
		check(thrown instanceof RuntimeException,
				"Worker aggregation should throw RuntimeException, not "+thrown);
		final Throwable cause=thrown.getCause();
		check(cause instanceof ArrayIndexOutOfBoundsException,
				"Primary cause should be the poison accession error, not "+cause);
		final StackTraceElement[] stack=cause.getStackTrace();
		check(stack.length>0 &&
				stack[0].getClassName().equals("tax.AnalyzeAccession") &&
				stack[0].getMethodName().equals("remapKey"),
				"Poison accession failed at the wrong site: "+
				Arrays.toString(stack));
	}

	/** Invokes a reflected worker aggregation that is expected to fail.
	 * @param aggregate Reflected worker aggregation method
	 * @param threads Worker array instance
	 * @return Primary failure thrown by the reflected method
	 * @throws Exception If reflection fails before method entry */
	private static Throwable aggregationFailure(final Method aggregate,
			final Object threads) throws Exception{
		try{
			aggregate.invoke(null, threads);
		}catch(InvocationTargetException e){
			return e.getCause();
		}
		throw new AssertionError("Worker aggregation should have failed");
	}

	/** Counts this process' descriptors pointing at a path with a requested access mode.
	 * @param path Input path or FIFO path
	 * @param mode O_ACCMODE bits to match
	 * @return Matching descriptor count
	 * @throws Exception If fdinfo cannot be read */
	private static int countFds(final Path path, final int mode) throws Exception{
		final Path real=path.toRealPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(
				Paths.get("/proc/self/fd"))){
			for(Path fd : fds){
				final Path target;
				try{
					target=Files.readSymbolicLink(fd);
				}catch(Exception e){
					continue;
				}
				if(target.equals(real) && (fdFlags(fd)&O_ACCMODE)==mode){
					count++;
				}
			}
		}
		return count;
	}

	/** Reads one /proc fd flags field.
	 * @param fd Descriptor path under /proc/self/fd
	 * @return Parsed octal flags
	 * @throws Exception If fdinfo cannot be read */
	private static int fdFlags(final Path fd) throws Exception{
		final Path info=Paths.get("/proc/self/fdinfo",
				fd.getFileName().toString());
		for(String line : Files.readAllLines(info, StandardCharsets.US_ASCII)){
			if(line.startsWith("flags:")){
				return Integer.parseInt(line.substring(6).trim(), 8);
			}
		}
		throw new IllegalStateException("No flags line in "+info);
	}

	/** Adds the Vector incubator module when this JDK exposes it.
	 * @param command Child JVM command to extend
	 */
	private static void addVectorModule(final List<String> command){
		try{
			final Class<?> moduleLayer=Class.forName("java.lang.ModuleLayer");
			final Object boot=moduleLayer.getMethod("boot").invoke(null);
			final Object optional=moduleLayer.getMethod("findModule",
					String.class).invoke(boot, "jdk.incubator.vector");
			final boolean present=((Boolean)optional.getClass()
					.getMethod("isPresent").invoke(optional)).booleanValue();
			if(present){
				command.add("--add-modules");
				command.add("jdk.incubator.vector");
			}
		}catch(ClassNotFoundException e){
			return;
		}catch(Exception e){
			throw new RuntimeException("Could not probe jdk.incubator.vector", e);
		}
	}

	/** Joins a possibly-started owned thread during cleanup.
	 * @param thread Thread to join, or null
	 * @param label Diagnostic thread label
	 * @param prior Earlier cleanup failure, or null
	 * @return Prior plus any new cleanup failure
	 */
	private static Throwable cleanupJoin(final Thread thread, final String label,
			final Throwable prior){
		if(thread==null){return prior;}
		try{
			thread.join(WORKER_TIMEOUT_MS);
			return thread.isAlive() ?
					combine(prior, new AssertionError(label+" did not terminate")) :
					prior;
		}catch(Throwable t){
			return combine(prior, t);
		}
	}

	/** Builds an accession2taxid input with a valid row before a poison byte.
	 * @return Input bytes whose poison row fails BuildDiskTables workers */
	private static byte[] poisonInput(){
		final byte[] prefix=HEALTHY.getBytes(StandardCharsets.US_ASCII);
		final byte[] suffix="\tBAD.1\t8\t0\n".getBytes(StandardCharsets.US_ASCII);
		final byte[] result=new byte[prefix.length+1+suffix.length];
		System.arraycopy(prefix, 0, result, 0, prefix.length);
		result[prefix.length]=(byte)0x80;
		System.arraycopy(suffix, 0, result, prefix.length+1, suffix.length);
		return result;
	}

	/** Reads one process' stdout fully.
	 * @param process Process to read
	 * @return Output bytes
	 * @throws Exception If the stream cannot be read */
	private static byte[] readAll(final Process process) throws Exception{
		final byte[] buffer=new byte[4096];
		final java.io.ByteArrayOutputStream baos=
				new java.io.ByteArrayOutputStream();
		for(int len=process.getInputStream().read(buffer); len>=0;
				len=process.getInputStream().read(buffer)){
			baos.write(buffer, 0, len);
		}
		return baos.toByteArray();
	}

	/** Force-kills one child JVM with a bounded wait.
	 * @param child Child process to terminate
	 * @param label Child diagnostic label
	 * @param stdout Child stdout path
	 * @param stderr Child stderr path
	 * @throws Exception If the child ignores forced termination */
	private static void killChild(final Process child, final String label,
			final Path stdout, final Path stderr) throws Exception{
		child.destroyForcibly();
		if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
			throw new AssertionError(label+" child ignored destroyForcibly\nstdout\n"+
					readAscii(stdout)+"\nstderr\n"+readAscii(stderr));
		}
	}

	/** Reads one US-ASCII file, returning empty text when absent.
	 * @param path File to read
	 * @return Decoded file contents
	 * @throws Exception If the file exists but cannot be read */
	private static String readAscii(final Path path) throws Exception{
		return Files.exists(path) ?
				new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}

	/** Collapses process output for a parent-visible one-line child report.
	 * @param text Text to compact
	 * @return Single-line compact text
	 */
	private static String compact(final String text){
		return text.trim().replaceAll("\\s+", " ");
	}

	/** Combines two cleanup failures.
	 * @param prior Earlier failure, or null
	 * @param next Later failure, or null
	 * @return Failure to rethrow, or null
	 */
	private static Throwable combine(final Throwable prior, final Throwable next){
		if(prior==null){return next;}
		if(next!=null && next!=prior){prior.addSuppressed(next);}
		return prior;
	}

	/** Rethrows a cleanup failure using its original type when possible.
	 * @param thrown Failure to rethrow
	 * @throws Exception If the failure is checked
	 */
	private static void rethrow(final Throwable thrown) throws Exception{
		if(thrown instanceof Error){throw (Error)thrown;}
		if(thrown instanceof Exception){throw (Exception)thrown;}
		throw new RuntimeException(thrown);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(
					Comparator.reverseOrder())::iterator){
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

	/** O_ACCMODE value for read-only descriptors. */
	private static final int O_RDONLY=0;
	/** O_ACCMODE value for read-write descriptors. */
	private static final int O_RDWR=2;
	/** Mask for Linux file access mode bits. */
	private static final int O_ACCMODE=3;
	/** Parent child-process wall-time bound. */
	private static final long CHILD_TIMEOUT_MS=60000;
	/** Time to wait for graceful child destruction before kill. */
	private static final long DESTROY_TIMEOUT_MS=1000;
	/** Time to wait for waitForThreads after its expected release point. */
	private static final long COLLECTOR_TIMEOUT_MS=10000;
	/** Time to wait for FIFO descriptor observation. */
	private static final long FD_TIMEOUT_MS=10000;
	/** Time to wait for each worker thread to end. */
	private static final long WORKER_TIMEOUT_MS=10000;
	/** Descriptor polling interval. */
	private static final long POLL_MS=10;

}
