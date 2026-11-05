package test.kmer;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;

import bloom.KCountArray;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableInterface;
import shared.Shared;
import shared.Timer;
import stream.FASTQ;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Finite AbstractKmerTableSet base-method contracts for failure propagation and cleanup.
 * Run each case in a fresh JVM under both assertion modes. Minimal test doubles isolate
 * the final base-class code from real file loading and real tables. Prefilter cases
 * test width-selection termination with no memory available for sketch allocation.
 * @author Jean, Shinobu
 */
public final class AbstractKmerTableSetContractTest{

	/** Accepts a case name and output directory through case= and dir= arguments. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(mode==null || dir==null){throw new IllegalArgumentException("Expected case=<case name> and dir=<new directory>");}
		switch(mode){
			case "gc-direct-control": gc(false, false); break;
			case "gc-direct-failure": gc(false, true); break;
			case "gc-parallel-control": gc(true, false); break;
			case "gc-parallel-failure": gc(true, true); break;
			case "load-control": load(false, false); break;
			case "load-primary-failure": load(true, false); break;
			case "load-extra-failure": load(false, true); break;
			case "hist-control": hist(dir, false); break;
			case "hist-full": hist(dir, true); break;
			case "prefilter-one-control": prefilter(true, 65535); break;
			case "prefilter-one-high": prefilter(true, 65536); break;
			case "prefilter-multi-control": prefilter(false, 65534); break;
			case "prefilter-multi-high": prefilter(false, 65535); break;
			case "prefilter-max-multi": maxIntPrefilter(); break;
			case "stats-instance-false": statsInstanceFalse(); break;
			case "stats-instance-true": statsInstanceTrue(); break;
			case "stats-static-false": statsStaticFalse(); break;
			case "stats-progress-control": statsProgressControl(); break;
			default: throw new IllegalArgumentException("Unknown case: "+mode);
		}
		System.out.println("PASS case="+mode);
	}

	/**
	 * Requires finite width selection followed by the zero-memory early return.
	 * A timed-out daemon cannot keep this fresh test JVM alive after the main fails.
	 */
	private static void prefilter(final boolean onePass, final int limit) throws InterruptedException{
		AbstractKmerTableSet.MASK_MIDDLE=false;
		final PrefilterCall call=new PrefilterCall(onePass, limit);
		call.start();
		call.join(2000);
		final String settings="onePass="+onePass+", threshold="+limit;
		check(!call.isAlive(), "PREFILTER_NONRETURN: width selection did not return within two seconds; "+settings);
		if(call.failure!=null){
			throw new AssertionError("PREFILTER_UNEXPECTED_FAILURE: "+settings, call.failure);
		}
		check(call.returned, "PREFILTER_NO_RESULT: worker ended without returning or throwing; "+settings);
		check(call.result==null && !call.tableSet.prefilter,
			"PREFILTER_WRONG_RETURN: zero memory must disable filtering without allocating a sketch; "+settings);
		System.out.println("PREFILTER_RETURN: "+settings+", result=null, enabled=false");
	}

	/** Checks allocated counter width and the saturated threshold without opening any inputs. */
	private static void maxIntPrefilter(){
		AbstractKmerTableSet.MASK_MIDDLE=false;
		AbstractKmerTableSet.DISPLAY_PROGRESS=false;
		AbstractKmerTableSet.DISPLAY_STATS=false;
		final MinimalTableSet tableSet=new MinimalTableSet();
		tableSet.configureMaxPrefilter();
		check(tableSet.in1.isEmpty() && tableSet.in2.isEmpty() && tableSet.extra.isEmpty(),
			"Maximum-threshold fixture must construct an empty sketch without file input");
		check(tableSet.processInput()==0 && tableSet.maxIntGateObserved,
			"Maximum-threshold fixture must reach allocation and load zero input keys");
		System.out.println("PREFILTER_MAX_GATE_PASS");
	}

	/** showstats=f should suppress base reporting while leaving processInput progress independent. */
	private static void statsInstanceFalse(){
		final String text=captureStats(false, true, true, 37);
		check(!text.contains("Estimated kmer capacity") && !text.contains("Input:")
			&& !text.contains("Unique Kmers:") && !text.contains("Load Time:"),
			"showStats=false must suppress capacity and summary reporting: "+visible(text));
		check(text.contains("After table allocation:"),
			"DISPLAY_PROGRESS must still allow allocation progress outside showStats: "+visible(text));
	}

	/** showstats=t should preserve literal read, base and unique-count fields. */
	private static void statsInstanceTrue(){
		final String text=captureStats(true, true, false, 37);
		check(text.contains("Estimated kmer capacity: \t1000"),
			"DISPLAY_STATS must keep capacity reporting: "+visible(text));
		check(text.contains("Input:                      \t12 reads \t\t345 bases."),
			"Summary must print literal read/base counters: "+visible(text));
		check(text.contains("Unique Kmers:               \t37"),
			"Summary must print the supplied unique count: "+visible(text));
		check(text.contains("Load Time:                  \t"),
			"Summary must print load timing: "+visible(text));
		check(!text.contains("After table allocation:") && !text.contains("After loading:"),
			"DISPLAY_PROGRESS=false must suppress memory progress only: "+visible(text));
	}

	/** Static DISPLAY_STATS=false should suppress capacity and summary regardless of showstats=t. */
	private static void statsStaticFalse(){
		final String text=captureStats(true, false, true, 37);
		check(!text.contains("Estimated kmer capacity") && !text.contains("Input:")
			&& !text.contains("Unique Kmers:") && !text.contains("Load Time:"),
			"DISPLAY_STATS=false must suppress capacity and summary reporting: "+visible(text));
		check(text.contains("After table allocation:"),
			"DISPLAY_PROGRESS must still allow allocation progress outside DISPLAY_STATS: "+visible(text));
	}

	/** DISPLAY_PROGRESS should add memory blocks without controlling summary visibility. */
	private static void statsProgressControl(){
		final String text=captureStats(true, true, true, 37);
		check(text.contains("Estimated kmer capacity: \t1000")
			&& text.contains("Input:                      \t12 reads \t\t345 bases.")
			&& text.contains("Unique Kmers:               \t37"),
			"DISPLAY_PROGRESS=true must not suppress stats output: "+visible(text));
		check(text.contains("After table allocation:") && text.contains("After loading:"),
			"DISPLAY_PROGRESS=true must add both allocation and post-load memory blocks: "+visible(text));
	}

	/** Runs the base reporting code with a minimal no-input table and captures stderr/outstream bytes. */
	private static String captureStats(final boolean instanceStats, final boolean staticStats,
			final boolean displayProgress, final long added){
		final PrintStream oldOutstream=AbstractKmerTableSet.outstream;
		final PrintStream oldErr=System.err;
		final boolean oldDisplayStats=AbstractKmerTableSet.DISPLAY_STATS;
		final boolean oldDisplayProgress=AbstractKmerTableSet.DISPLAY_PROGRESS;
		final boolean oldShowSpeed=AbstractKmerTableSet.showSpeed;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(bytes);
		try{
			System.setErr(capture);
			AbstractKmerTableSet.outstream=capture;
			AbstractKmerTableSet.DISPLAY_STATS=staticStats;
			AbstractKmerTableSet.DISPLAY_PROGRESS=displayProgress;
			AbstractKmerTableSet.showSpeed=false;

			final MinimalTableSet tableSet=new MinimalTableSet();
			tableSet.showStats=instanceStats;
			tableSet.readsIn=12;
			tableSet.basesIn=345;
			check(tableSet.processInput()==0 && tableSet.in1.isEmpty() && tableSet.extra.isEmpty(),
				"Stats fixture must allocate an empty table without loading files");
			tableSet.showStats(new Timer((PrintStream)null), added);
			capture.flush();
			return new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		}finally{
			AbstractKmerTableSet.showSpeed=oldShowSpeed;
			AbstractKmerTableSet.DISPLAY_PROGRESS=oldDisplayProgress;
			AbstractKmerTableSet.DISPLAY_STATS=oldDisplayStats;
			AbstractKmerTableSet.outstream=oldOutstream;
			System.setErr(oldErr);
			capture.close();
		}
	}

	/** Makes captured newlines visible in failed assertions. */
	private static String visible(final String text){
		return text.replace("\n", "\\n");
	}

	/** Owns one isolated production call; captures failures independently of its normal result. */
	private static final class PrefilterCall extends Thread{

		/** Configures a memoryless fixture; the daemon boundary limits a defective infinite loop. */
		PrefilterCall(final boolean onePass, final int limit){
			super("prefilter-width-check");
			setDaemon(true);
			tableSet.configurePrefilter(onePass, limit);
		}

		/** Runs only width selection and the zero-memory gate, retaining any unexpected failure. */
		@Override
		public void run(){
			try{
				result=tableSet.makePrefilter(new KCountArray[1], null);
				returned=true;
			}catch(Throwable t){
				failure=t;
			}
		}

		/** Exclusive minimal state; no input files or table allocations. */
		final MinimalTableSet tableSet=new MinimalTableSet();
		/** Normal result, read only after the worker has terminated and join has returned. */
		KCountArray result;
		/** Distinguishes a normal null return from a call that never completed. */
		boolean returned;
		/** Unexpected failure; checked separately from timeout and normal return. */
		Throwable failure;
	}

	/** Checks direct and worker GC accumulation, including the KMR108 worker failure path. */
	private static void gc(final boolean parallel, final boolean fail){
		final int oldThreads=Shared.threads();
		try{
			Shared.setThreads(parallel ? 2 : 1);
			final RuntimeException injected=new RuntimeException("injected-gc");
			final MinimalTableSet tableSet=new MinimalTableSet();
			tableSet.tables=new KmerTableInterface[] {
				new CountingTable(new long[] {1, 2, 3}, null),
				new CountingTable(new long[] {4, 5, 6}, fail ? injected : null)
			};
			final long[] gcCounts={10, 20, 30};
			RuntimeException failure=null;
			try{
				tableSet.countGC(gcCounts, 2);
			}catch(RuntimeException re){
				failure=re;
			}
			if(fail){
				check(hasCause(failure, injected), "Injected GC failure must reach the caller: "+failure);
			}else{
				check(failure==null, "GC control must not fail: "+failure);
				check(Arrays.equals(gcCounts, new long[] {15, 27, 39}),
					"GC contributions must accumulate into existing bins: "+Arrays.toString(gcCounts));
			}
		}finally{
			Shared.setThreads(oldThreads);
		}
	}

	/** Checks processInput loading cleanup for normal, primary-failure and extra-failure paths. */
	private static void load(final boolean primaryFailure, final boolean extraFailure){
		final boolean oldValidate=Read.VALIDATE_IN_CONSTRUCTOR;
		final boolean oldForceInterleaved=FASTQ.FORCE_INTERLEAVED;
		final boolean oldDisplayProgress=AbstractKmerTableSet.DISPLAY_PROGRESS;
		final boolean oldDisplayStats=AbstractKmerTableSet.DISPLAY_STATS;
		try{
			Read.VALIDATE_IN_CONSTRUCTOR=true;
			FASTQ.FORCE_INTERLEAVED=true;
			AbstractKmerTableSet.DISPLAY_PROGRESS=false;
			AbstractKmerTableSet.DISPLAY_STATS=false;

			final RuntimeException injected=new RuntimeException("injected-load");
			final MinimalTableSet tableSet=new MinimalTableSet();
			tableSet.in1.add("primary.fa");
			tableSet.extra.add("extra.fa");
			tableSet.loadReturns=new long[] {5, 7};
			tableSet.loadFailureCall=primaryFailure ? 1 : extraFailure ? 2 : -1;
			tableSet.loadFailure=injected;

			RuntimeException failure=null;
			long loaded=-1;
			try{
				loaded=tableSet.processInput();
			}catch(RuntimeException re){
				failure=re;
			}
			if(primaryFailure || extraFailure){
				check(hasCause(failure, injected), "Injected load failure must reach the caller: "+failure);
				check(tableSet.observedValidate.equals(primaryFailure ?
					Arrays.asList(Boolean.FALSE) : Arrays.asList(Boolean.FALSE, Boolean.FALSE)),
					"Load failures must run with validation disabled: "+tableSet.observedValidate);
				check(tableSet.observedForceInterleaved.equals(primaryFailure ?
					Arrays.asList(Boolean.TRUE) : Arrays.asList(Boolean.TRUE, Boolean.FALSE)),
					"Primary loads must keep the prior interleaving flag; extra loads must disable it: "+
					tableSet.observedForceInterleaved);
				check(Read.VALIDATE_IN_CONSTRUCTOR, "processInput must restore Read.VALIDATE_IN_CONSTRUCTOR after a load failure");
				check(FASTQ.FORCE_INTERLEAVED, "processInput must restore FASTQ.FORCE_INTERLEAVED after an extra-load failure");
			}else{
				check(failure==null, "Load control must not fail: "+failure);
				check(loaded==12, "Load control must return the sum of primary and extra loader counts: "+loaded);
				check(tableSet.calls.equals(Arrays.asList("primary.fa|null", "extra.fa|null")),
					"Load control must call primary and extra inputs in order: "+tableSet.calls);
				check(tableSet.observedValidate.equals(Arrays.asList(Boolean.FALSE, Boolean.FALSE)),
					"Normal loads must run with validation disabled: "+tableSet.observedValidate);
				check(tableSet.observedForceInterleaved.equals(Arrays.asList(Boolean.TRUE, Boolean.FALSE)),
					"Primary loads must keep the prior interleaving flag; extra loads must disable it: "+
					tableSet.observedForceInterleaved);
				check(Read.VALIDATE_IN_CONSTRUCTOR && FASTQ.FORCE_INTERLEAVED,
					"Load control must restore both process-wide parser flags");
			}
		}finally{
			Read.VALIDATE_IN_CONSTRUCTOR=oldValidate;
			FASTQ.FORCE_INTERLEAVED=oldForceInterleaved;
			AbstractKmerTableSet.DISPLAY_PROGRESS=oldDisplayProgress;
			AbstractKmerTableSet.DISPLAY_STATS=oldDisplayStats;
		}
	}

	/** Checks makeKhist return identity, literal output bytes and /dev/full failure reporting. */
	private static void hist(final Path dir, final boolean full) throws Exception{
		if(full){requireDevFull();}
		final MinimalTableSet tableSet=new MinimalTableSet();
		final long[] histogram={0, 3, 0, 2};
		tableSet.histogram=histogram;
		final String fname;
		if(full){
			fname="/dev/full";
		}else{
			Files.createDirectories(dir);
			final Path output=dir.resolve("hist.txt");
			check(!Files.exists(output), "Histogram output must be new; stale bytes cannot satisfy the oracle");
			fname=output.toString();
		}

		RuntimeException failure=null;
		long[][] result=null;
		try{
			result=tableSet.makeKhist(fname, 2, 3, true, true, false, false, false, false, 2.0, 1, 1);
		}catch(RuntimeException re){
			if(full){re.printStackTrace(System.err);}
			failure=re;
		}
		if(full){
			check(failure!=null, "A /dev/full histogram write must throw a RuntimeException after writer completion");
		}else{
			check(failure==null, "Histogram control must not fail: "+failure);
			check(result!=null && result[0]==histogram && result[1]==null,
				"Histogram control must return the original count array and no GC array");
			final byte[] expected="#Depth\tCount\n1\t3\n2\t0\n3\t2\n".getBytes(StandardCharsets.US_ASCII);
			final byte[] actual=Files.readAllBytes(Paths.get(fname));
			check(Arrays.equals(actual, expected),
				"Histogram bytes differ: "+new String(actual, StandardCharsets.US_ASCII).replace("\n", "\\n"));
		}
	}

	/** Restricts /dev/full oracles to Linux hosts with the expected special file. */
	private static void requireDevFull(){
		check(System.getProperty("os.name").toLowerCase().contains("linux"),
			"/dev/full failure cases are Linux-only");
		check(Files.exists(Paths.get("/dev/full")) && !Files.isRegularFile(Paths.get("/dev/full")),
			"/dev/full must exist as a special write-failing device");
	}

	/** True when the exact injected exception appears in the causal chain. */
	private static boolean hasCause(Throwable t, final RuntimeException injected){
		for(; t!=null; t=t.getCause()){
			if(t==injected){return true;}
		}
		return false;
	}

	/** Keeps test oracles active under -da and supplies a concrete failed contract. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Minimal table-set test double for exercising final AbstractKmerTableSet methods. */
	private static final class MinimalTableSet extends AbstractKmerTableSet{

		/** Enables the requested threshold while retaining the zero-memory filter budget. */
		void configurePrefilter(final boolean onePass_, final int limit){
			prefilter=true;
			onePass=onePass_;
			filterMax=limit;
		}

		/** Sets protected pass configuration inside the subclass and enables the allocated oracle. */
		void configureMaxPrefilter(){
			configurePrefilter(false, Integer.MAX_VALUE);
			prepasses=1;
			filterMemoryBytes=1L<<20;
			maxIntGateCheck=true;
		}

		@Override
		public void clear(){tables=null;}

		@Override
		public long loadKmers(String fname1, String fname2){
			calls.add(fname1+"|"+fname2);
			observedValidate.add(Read.VALIDATE_IN_CONSTRUCTOR);
			observedForceInterleaved.add(FASTQ.FORCE_INTERLEAVED);
			loadCalls++;
			if(loadCalls==loadFailureCall){throw loadFailure;}
			return loadReturns[loadCalls-1];
		}

		@Override
		public long regenerate(int limit){throw new AssertionError("Unexpected regenerate call");}

		@Override
		public Object getTable(int tnum){return tables[tnum];}

		@Override
		public long[] fillHistogram(int histMax){return histogram;}

		@Override
		public void initializeOwnership(){throw new AssertionError("Unexpected initializeOwnership call");}

		@Override
		public void clearOwnership(){throw new AssertionError("Unexpected clearOwnership call");}

		@Override
		public int ways(){return tables.length;}

		@Override
		public int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, Kmer kmer){
			throw new AssertionError("Unexpected fillSpecificCounts call");
		}

		@Override
		public int regenerateCounts(byte[] bases, IntList counts, Kmer kmer, BitSet changed){
			throw new AssertionError("Unexpected regenerateCounts call");
		}

		@Override
		public boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining){
			throw new AssertionError("Unexpected dumpKmersAsBytes call");
		}

		@Override
		public boolean dumpKmersAsBytes_MT(String fname, int mincount, int maxcount,
				boolean printTime, AtomicLong remaining){
			throw new AssertionError("Unexpected dumpKmersAsBytes_MT call");
		}

		@SuppressWarnings("unchecked")
		@Override
		public <T extends KmerTableInterface> T[] tables(){return (T[])tables;}

		@Override
		public int kbig(){return 3;}

		@Override
		public long filterMemory(int pass){return filterMemoryBytes;}

		@Override
		public long tableMemory(){return 0;}

		@Override
		public long estimatedKmerCapacity(){return 1000;}

		@Override
		public boolean ecco(){return false;}

		@Override
		public boolean qtrimLeft(){return false;}

		@Override
		public boolean qtrimRight(){return false;}

		@Override
		public float minAvgQuality(){return 0;}

		@Override
		public boolean rcomp(){return false;}

		@Override
		protected void allocateTables(){
			if(maxIntGateCheck){
				check(prefilterArray!=null, "PREFILTER_MAX_WIDTH: expected an allocated sketch");
				final KCountArray filter=prefilterArray;
				final String observed="bits="+filter.cellBits+", maximum="+filter.maxValue+", gate="+filterMax2;
				System.out.println("PREFILTER_MAX_OBSERVED: "+observed);
				check(filter.cellBits==32 && filter.maxValue==Integer.MAX_VALUE && filterMax2==Integer.MAX_VALUE-1,
					"PREFILTER_MAX_WIDTH: maximum int requires 32-bit counters and a saturated gate; "+observed);
				check(filter.incrementAndReturnUnincremented(42L, Integer.MAX_VALUE-1)==0,
					"Empty sketch must return zero before its first increment");
				check(filter.read(42L)==Integer.MAX_VALUE-1 && filter.read(42L)<=filterMax2,
					"A count at the gate must not pass the strict post-filter comparison");
				final int before=filter.incrementAndReturnUnincremented(42L, 1);
				check(before==Integer.MAX_VALUE-1 && before>=filterMax2 && filter.read(42L)>filterMax2,
					"Reaching saturation must pass both old-count one-pass and new-count post-filter gates");
				check(filter.incrementAndReturnUnincremented(42L, 1)==Integer.MAX_VALUE &&
					filter.read(42L)==Integer.MAX_VALUE, "A saturated counter must remain saturated");
				maxIntGateObserved=true;
			}
			if(tables==null){tables=new KmerTableInterface[0];}
		}

		/** Filter budget; zero preserves the original early-return fixtures. */
		long filterMemoryBytes;
		/** Enables the allocated maximum-count oracle before processInput clears the sketch. */
		boolean maxIntGateCheck;
		/** Records that the complete width and saturation oracle ran successfully. */
		boolean maxIntGateObserved;
		KmerTableInterface[] tables=new KmerTableInterface[0];
		long[] histogram={0};
		long[] loadReturns={0};
		int loadFailureCall=-1, loadCalls=0;
		RuntimeException loadFailure=null;
		final ArrayList<String> calls=new ArrayList<String>();
		final ArrayList<Boolean> observedValidate=new ArrayList<Boolean>();
		final ArrayList<Boolean> observedForceInterleaved=new ArrayList<Boolean>();
	}

	/** Deterministic countGC delegate with optional distinct injected failure. */
	private static final class CountingTable implements KmerTableInterface{

		CountingTable(final long[] contribution_, final RuntimeException failure_){
			contribution=contribution_;
			failure=failure_;
		}

		@Override
		public void countGC(long[] gcCounts, int max){
			if(failure!=null){throw failure;}
			for(int i=0; i<contribution.length; i++){
				gcCounts[i]+=contribution[i];
			}
		}

		private final long[] contribution;
		private final RuntimeException failure;
	}
}
