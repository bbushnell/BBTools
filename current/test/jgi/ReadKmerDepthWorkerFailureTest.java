package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLongArray;

import bloom.KCountArray;
import dna.AminoAcid;
import jgi.ReadKmerDepthDistribution;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Child-JVM fixture for ReadKmerDepthDistribution worker and histogram failures.
 * Each invocation must run in a fresh JVM so deliberate exit2 and static histogram
 * fields stay local. The fixture also assumes the production default of 8 workers.
 * Histogram cases require {@code hist=<path>}; the capture controller validates the path
 * and independently checks the expected I/O diagnostics for hist-fail.
 * Count cases require a fresh {@code in=<path>} that this fixture creates before
 * invoking the real count entry point. Retained-output count cases also require
 * {@code out=<path>}.
 *
 * @author Jean
 */
public final class ReadKmerDepthWorkerFailureTest{

	/**
	 * Runs one requested child-JVM scenario.
	 *
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, ReadKmerDepthWorkerFailureTest.class, false);
		final String[] parsedArgs=(pp.args==null ? new String[0] : pp.args);

		try{
			final TestCase testCase=parseCase(parsedArgs);
			if(testCase==TestCase.CONTROL){
				runControl();
			}else if(testCase==TestCase.WORKER_NULL){
				runWorkerNull();
			}else if(testCase==TestCase.HIST_CONTROL){
				runHistogramControl();
			}else if(testCase==TestCase.HIST_FAIL){
				runHistogramFailure();
			}else if(testCase==TestCase.COUNT_CONTROL){
				runCountControl();
			}else if(testCase==TestCase.COUNT_PARENT){
				runCountParent();
			}else if(testCase==TestCase.COUNT_OUTPUT_CONTROL){
				runCountOutputControl();
			}else if(testCase==TestCase.COUNT_OUTPUT_FAIL){
				runCountOutputFail();
			}else{
				throw new AssertionError("Unhandled case "+testCase);
			}
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the successful control case. */
	private static void runControl(){
		final InMemoryCris cris=makeCris(new Read(new byte[] {'A'}, null, "control", 0));
		prepareReadKmerDepth(false);
		System.out.println("RKDD_WORKER_FAILURE_READY case=control");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.downsample(cris, null, 31, -1, null, null, true, 0);
		verify(bases==1, "Expected one returned base but saw "+bases);
		cris.verifyState(9, 1, 8, 0);
		System.out.println("RKDD_WORKER_CONTROL_RETURNED bases="+bases);
	}

	/**
	 * Runs the baseline worker-failure case and deliberately exits 2 after
	 * proving that the current owner thread returns despite a worker failure.
	 */
	private static void runWorkerNull(){
		final InMemoryCris cris=makeCris(null);
		prepareReadKmerDepth(false);
		System.out.println("RKDD_WORKER_FAILURE_READY case=worker-null");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.downsample(cris, null, 31, -1, null, null, true, 0);
		verify(bases==0, "Expected zero returned bases but saw "+bases);
		cris.verifyState(8, 0, 7, 1);
		System.out.println("RKDD_FATAL_WORKER_RETURNED bases="+bases);
		System.out.flush();
		System.exit(2);
	}

	/** Runs the successful real-count control case. */
	private static void runCountControl(){
		prepareReadKmerDepth(false);
		writeOneBaseFasta();
		System.out.println("RKDD_COUNT_CONTROL_READY case=count-control");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.count(inputFile, null, null,
				31, -1, null, true, null, 0);
		verify(bases==1, "Expected one returned base but saw "+bases);
		System.out.println("RKDD_COUNT_CONTROL_RETURNED bases="+bases);
	}

	/**
	 * Runs the baseline parent-failure case and deliberately exits 2 after
	 * proving that the parent exception reached the fixture.
	 */
	private static void runCountParent(){
		prepareReadKmerDepth(false);
		ReadKmerDepthDistribution.THREAD_HIST_LEN=-1;
		writeOneBaseFasta();
		System.out.println("RKDD_COUNT_PARENT_READY case=count-parent");
		System.out.flush();
		try{
			ReadKmerDepthDistribution.count(inputFile, null, null, 31, -1,
					null, true, null, 0);
		}catch(NegativeArraySizeException e){
			System.out.println("RKDD_COUNT_PARENT_CAUGHT java.lang.NegativeArraySizeException");
			e.printStackTrace(System.err);
			System.err.flush();
			System.out.flush();
			System.exit(2);
		}
		throw new AssertionError("count-parent unexpectedly returned");
	}

	/** Runs the successful retained-read attachment-output control case. */
	private static void runCountOutputControl(){
		prepareReadKmerDepth(false);
		writeRetainedReadFasta();
		final KCountArray kca=makeRetainedReadCounts();
		System.out.println("RKDD_COUNT_OUTPUT_READY case=count-output-control");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.count(inputFile, null, kca,
				RETAINED_K, -1, outputFile, true, null, 0);
		verify(bases==RETAINED_BASES.length(),
				"Expected "+RETAINED_BASES.length()+" returned bases but saw "+bases);
		verifyRetainedOutput();
		System.out.println("RKDD_COUNT_OUTPUT_RETURNED case=count-output-control bases="+bases);
	}

	/**
	 * Runs the baseline retained-read output-failure case and deliberately exits 2
	 * if count returns despite the close error. That exit2 is harness-controlled,
	 * not evidence that the production code naturally terminated.
	 */
	private static void runCountOutputFail(){
		prepareReadKmerDepth(false);
		writeRetainedReadFasta();
		final KCountArray kca=makeRetainedReadCounts();
		System.out.println("RKDD_COUNT_OUTPUT_READY case=count-output-fail");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.count(inputFile, null, kca,
				RETAINED_K, -1, outputFile, true, null, 0);
		verify(bases==RETAINED_BASES.length(),
				"Expected "+RETAINED_BASES.length()+" returned bases but saw "+bases);
		System.out.println("RKDD_COUNT_OUTPUT_RETURNED case=count-output-fail bases="+bases);
		System.out.flush();
		System.exit(2);
	}

	/** Runs the successful histogram-control case. */
	private static void runHistogramControl(){
		final InMemoryCris cris=makeCris(new Read(new byte[] {'A'}, null, "hist-control", 0));
		prepareReadKmerDepth(true);
		System.out.println("RKDD_WORKER_FAILURE_READY case=hist-control");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.downsample(cris, null, 31,
				-1, null, histFile, true, 0);
		verify(bases==1, "Expected one returned base but saw "+bases);
		cris.verifyState(9, 1, 8, 0);
		verifyHistogramFile();
		System.out.println("RKDD_HISTOGRAM_CONTROL_RETURNED bases="+bases);
	}

	/**
	 * Runs the histogram-failure capture and deliberately exits 2 on normal return.
	 * The controller independently verifies that the supplied output path failed.
	 */
	private static void runHistogramFailure(){
		final InMemoryCris cris=makeCris(new Read(new byte[] {'A'}, null, "hist-fail", 0));
		prepareReadKmerDepth(true);
		System.out.println("RKDD_WORKER_FAILURE_READY case=hist-fail");
		System.out.flush();
		final long bases=ReadKmerDepthDistribution.downsample(cris, null, 31,
				-1, null, histFile, true, 0);
		verify(bases==1, "Expected one returned base but saw "+bases);
		cris.verifyState(9, 1, 8, 0);
		System.out.println("RKDD_HISTOGRAM_FAILURE_WRONG_RETURN bases="+bases);
		System.out.flush();
		System.exit(2);
	}

	/** Writes the literal one-base input for real-count cases. */
	private static void writeOneBaseFasta(){
		try{
			Files.write(Paths.get(inputFile), ">count\nA\n".getBytes(StandardCharsets.US_ASCII),
					StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/** Writes the literal A31C read used by retained-read output cases. */
	private static void writeRetainedReadFasta(){
		try{
			final String text=">kept\n"+RETAINED_BASES+"\n";
			Files.write(Paths.get(inputFile), text.getBytes(StandardCharsets.US_ASCII),
					StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/**
	 * Builds a public small KCountArray for the retained read and verifies exact
	 * public 100/1 readback before the caller prints its READY line.
	 */
	private static KCountArray makeRetainedReadCounts(){
		Shared.setThreads(1);
		final KCountArray kca=KCountArray.makeNew(1024, 8, 1);
		final long[] keys=retainedCanonicalKeys();
		verify(keys[0]!=keys[1], "Retained k-mers unexpectedly canonicalized together");
		kca.increment(keys[0], 100);
		kca.increment(keys[1], 1);
		verify(kca.read(keys[0], RETAINED_K, true)==100,
				"High-depth retained k-mer did not read back exactly");
		verify(kca.read(keys[1], RETAINED_K, true)==1,
				"Low-depth retained k-mer did not read back exactly");
		return kca;
	}

	/** Returns the two canonical k=31 keys in the literal retained read. */
	private static long[] retainedCanonicalKeys(){
		final long[] keys=new long[2];
		long kmer=0;
		int len=0;
		int keyCount=0;
		final byte[] bases=RETAINED_BASES.getBytes(StandardCharsets.US_ASCII);
		for(final byte base : bases){
			final int x=AminoAcid.baseToNumber[base];
			verify(x>=0, "Unexpected non-ACGT base: "+(char)base);
			kmer=((kmer<<2)|x)&RETAINED_MASK;
			len++;
			if(len>=RETAINED_K){
				verify(keyCount<keys.length, "Unexpected retained k-mer count");
				keys[keyCount]=KCountArray.makeCanonical2(kmer, RETAINED_K);
				keyCount++;
			}
		}
		verify(keyCount==keys.length, "Expected "+keys.length+" retained k-mers but saw "+keyCount);
		return keys;
	}

	/** Verifies the exact six attachment bytes for a retained A31C read. */
	private static void verifyRetainedOutput(){
		try{
			final Path path=Paths.get(outputFile);
			final byte[] expected="100\t1\n".getBytes(StandardCharsets.US_ASCII);
			verify(Files.isRegularFile(path),
					"Retained output requires a regular file: "+outputFile);
			final long length=Files.size(path);
			verify(length==expected.length,
					"Expected "+expected.length+" retained output bytes before reading, saw "
					+length);
			final byte[] actual=Files.readAllBytes(path);
			verify(Arrays.equals(expected, actual),
					"Unexpected retained output: "+new String(actual, StandardCharsets.US_ASCII));
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/**
	 * Parses the required case plus the histogram or input path required by
	 * cases that create external files.
	 *
	 * @param args Command line arguments
	 * @return Selected test case
	 */
	private static TestCase parseCase(final String[] args){
		TestCase testCase=null;
		histFile=null;
		inputFile=null;
		outputFile=null;
		for(final String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase(Locale.ROOT);
			final String value=(split.length>1 ? split[1] : null);
			if("case".equals(key)){
				require(value!=null && value.length()>0, "case requires a value");
				final TestCase parsed=TestCase.parse(value.toLowerCase(Locale.ROOT));
				require(testCase==null || testCase==parsed,
						"Conflicting case values: "+testCase+" and "+parsed);
				testCase=parsed;
			}else if("hist".equals(key)){
				require(value!=null && value.length()>0, "hist requires a value");
				require(histFile==null || histFile.equals(value),
						"Conflicting hist values: "+histFile+" and "+value);
				histFile=value;
			}else if("in".equals(key)){
				require(value!=null && value.length()>0, "in requires a value");
				require(inputFile==null || inputFile.equals(value),
						"Conflicting in values: "+inputFile+" and "+value);
				inputFile=value;
			}else if("out".equals(key)){
				require(value!=null && value.length()>0, "out requires a value");
				require(outputFile==null || outputFile.equals(value),
						"Conflicting out values: "+outputFile+" and "+value);
				outputFile=value;
			}else{
				require(false, "Unknown argument: "+arg);
			}
		}
		require(testCase!=null, "Missing required case argument");
		if(testCase==TestCase.HIST_CONTROL || testCase==TestCase.HIST_FAIL){
			require(histFile!=null, "Histogram cases require hist=<path>");
		}else{
			require(histFile==null, "hist is only valid for histogram cases");
		}
		if(testCase==TestCase.COUNT_CONTROL || testCase==TestCase.COUNT_PARENT){
			require(inputFile!=null, "Count cases require in=<path>");
		}else if(testCase==TestCase.COUNT_OUTPUT_CONTROL || testCase==TestCase.COUNT_OUTPUT_FAIL){
			require(inputFile!=null, "Count output cases require in=<path>");
		}else{
			require(inputFile==null, "in is only valid for count cases");
		}
		if(testCase==TestCase.COUNT_OUTPUT_CONTROL || testCase==TestCase.COUNT_OUTPUT_FAIL){
			require(outputFile!=null, "Count output cases require out=<path>");
		}else{
			require(outputFile==null, "out is only valid for count output cases");
		}
		return testCase;
	}

	/** Sets bounded public histogram sizes and whether downsample writes a histogram. */
	private static void prepareReadKmerDepth(final boolean useHistogram){
		ReadKmerDepthDistribution.THREAD_HIST_LEN=4;
		ReadKmerDepthDistribution.HIST_LEN=4;
		ReadKmerDepthDistribution.HIST_LEN_PRINT=4;
		ReadKmerDepthDistribution.PRINT_ZERO_COVERAGE=false;
		ReadKmerDepthDistribution.USE_HISTOGRAM=useHistogram;
		ReadKmerDepthDistribution.histogram_total=new AtomicLongArray(4);
	}

	/** Verifies that the short-read control wrote only the histogram header. */
	private static void verifyHistogramFile(){
		try{
			final Path path=Paths.get(histFile);
			verify(Files.isRegularFile(path), "Histogram control requires a regular file: "+histFile);
			final long length=Files.size(path);
			final int expectedLength=HISTOGRAM_HEADER.getBytes(StandardCharsets.UTF_8).length;
			verify(length==expectedLength,
					"Expected "+expectedLength+" histogram bytes before reading, saw "+length);
			final byte[] bytes=Files.readAllBytes(path);
			final String text=new String(bytes, StandardCharsets.UTF_8);
			verify(HISTOGRAM_HEADER.equals(text), "Unexpected histogram text: "+text);
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/**
	 * Makes one data batch followed by eight terminal batches, matching the
	 * production default worker count exactly.
	 *
	 * @param read Data read, which may be null for a worker-failure batch
	 * @return Synthetic input stream
	 */
	private static InMemoryCris makeCris(final Read read){
		@SuppressWarnings("unchecked")
		final ListNum<Read>[] lists=new ListNum[LIST_COUNT];
		final ArrayList<Read> data=new ArrayList<Read>(1);
		data.add(read);
		lists[DATA_ID]=new ListNum<Read>(data, DATA_ID);
		for(int i=1; i<lists.length; i++){
			lists[i]=new ListNum<Read>(new ArrayList<Read>(0), i);
		}
		return new InMemoryCris(lists);
	}

	/** Fails for malformed test input. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Fails for a violated oracle. */
	private static void verify(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Supported child-JVM cases. */
	private enum TestCase{
		CONTROL,
		WORKER_NULL,
		HIST_CONTROL,
		HIST_FAIL,
		COUNT_CONTROL,
		COUNT_PARENT,
		COUNT_OUTPUT_CONTROL,
		COUNT_OUTPUT_FAIL;

		/**
		 * Parses a command-line case value.
		 *
		 * @param value Value to parse
		 * @return Matching case
		 */
		private static TestCase parse(final String value){
			if("control".equals(value)){return CONTROL;}
			if("worker-null".equals(value)){return WORKER_NULL;}
			if("hist-control".equals(value)){return HIST_CONTROL;}
			if("hist-fail".equals(value)){return HIST_FAIL;}
			if("count-control".equals(value)){return COUNT_CONTROL;}
			if("count-parent".equals(value)){return COUNT_PARENT;}
			if("count-output-control".equals(value)){return COUNT_OUTPUT_CONTROL;}
			if("count-output-fail".equals(value)){return COUNT_OUTPUT_FAIL;}
			throw new IllegalArgumentException("Unknown case: "+value);
		}
	}

	/** ConcurrentReadInputStream that delivers fixed in-memory lists once. */
	private static final class InMemoryCris extends ConcurrentReadInputStream{

		/** @param lists_ Lists to deliver */
		private InMemoryCris(final ListNum<Read>[] lists_){
			super("rkdd-fatal-in-memory");
			lists=lists_;
			delivered=new boolean[lists.length];
			returned=new boolean[lists.length];
		}

		@Override
		public synchronized ListNum<Read> nextList(){
			deliveryAttempts++;
			verify(cursor<lists.length, "Unexpected delivery request "+deliveryAttempts);
			final ListNum<Read> list=lists[cursor];
			cursor++;
			final int id=toIndex(list.id);
			verify(!delivered[id], "Duplicate delivery id "+list.id);
			delivered[id]=true;
			return list;
		}

		@Override
		public synchronized void returnList(final long id, final boolean poison){
			returnAttempts++;
			final int index=toIndex(id);
			verify(delivered[index], "Returned undelivered id "+id);
			verify(!returned[index], "Returned id "+id+" more than once");
			returned[index]=true;
			if(index==DATA_ID){
				verify(!poison, "Data batch was returned as terminal");
				dataReturns++;
			}else{
				verify(poison, "Terminal batch was not returned as terminal");
				terminalReturns++;
			}
		}

		/**
		 * Verifies delivery and return counts.
		 *
		 * @param expectedDeliveries Expected delivered lists
		 * @param expectedDataReturns Expected data-list returns
		 * @param expectedTerminalReturns Expected terminal-list returns
		 * @param expectedSpareTerminals Expected undelivered terminal lists
		 */
		private synchronized void verifyState(final int expectedDeliveries,
				final int expectedDataReturns, final int expectedTerminalReturns,
				final int expectedSpareTerminals){
			final int expectedReturnAttempts=expectedDataReturns+expectedTerminalReturns;
			verify(deliveryAttempts==expectedDeliveries,
					"Expected "+expectedDeliveries+" delivery attempts but saw "+deliveryAttempts);
			verify(cursor==expectedDeliveries,
					"Expected "+expectedDeliveries+" deliveries but saw "+cursor);
			verify(returnAttempts==expectedReturnAttempts,
					"Expected "+expectedReturnAttempts+" return attempts but saw "+returnAttempts);
			verify(dataReturns==expectedDataReturns,
					"Expected "+expectedDataReturns+" data returns but saw "+dataReturns);
			verify(terminalReturns==expectedTerminalReturns,
					"Expected "+expectedTerminalReturns+" terminal returns but saw "
					+terminalReturns);
			verify(spareTerminalCount()==expectedSpareTerminals,
					"Expected "+expectedSpareTerminals+" spare terminals but saw "
					+spareTerminalCount());
		}

		/** @return Number of undelivered terminal batches */
		private int spareTerminalCount(){
			int count=0;
			for(int i=1; i<delivered.length; i++){
				if(!delivered[i]){count++;}
			}
			return count;
		}

		/** Converts a long list ID to a local array index. */
		private int toIndex(final long id){
			verify(id>=0 && id<lists.length, "Unexpected list id "+id);
			return (int)id;
		}

		@Override
		public void run(){
			// Nothing to start; lists are delivered directly from nextList().
		}

		@Override
		public void shutdown(){
			// Nothing to stop.
		}

		@Override
		public void restart(){
			throw new UnsupportedOperationException();
		}

		@Override
		public void close(){
			// Nothing to close.
		}

		@Override
		public boolean paired(){return false;}

		@Override
		public Object[] producers(){return null;}

		@Override
		public boolean errorState(){return false;}

		@Override
		public void setSampleRate(final float rate, final long seed){
			// The synthetic stream is deterministic and not sampled.
		}

		@Override
		public long basesIn(){return 0;}

		@Override
		public long readsIn(){return 0;}

		@Override
		public boolean verbose(){return false;}

		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/

		/** Lists to deliver. */
		private final ListNum<Read>[] lists;

		/** Delivered list IDs. */
		private final boolean[] delivered;

		/** Returned list IDs. */
		private final boolean[] returned;

		/** Next delivery index. */
		private int cursor=0;

		/** Delivery attempts, including invalid extra requests. */
		private int deliveryAttempts=0;

		/** Return attempts, including invalid duplicate attempts. */
		private int returnAttempts=0;

		/** Non-terminal returns. */
		private int dataReturns=0;

		/** Terminal returns. */
		private int terminalReturns=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Production default worker count and terminal-list count. */
	private static final int DEFAULT_THREADS=8;

	/** Data-list ID. */
	private static final int DATA_ID=0;

	/** Total delivered lists for the control path. */
	private static final int LIST_COUNT=DEFAULT_THREADS+1;

	/** The only histogram bytes produced by a one-base read with zero coverage bins. */
	private static final String HISTOGRAM_HEADER="#Depth\tRaw_Count\tUnique_Kmers\n";

	/** k value for the A31C retained-read output cases. */
	private static final int RETAINED_K=31;

	/** Bitmask for two-bit k=31 primitive kmers. */
	private static final long RETAINED_MASK=~((-1L)<<(2*RETAINED_K));

	/** Literal 32-base read with two k=31 windows and no spike-fix interior. */
	private static final String RETAINED_BASES="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC";

	/** Histogram output path for histogram cases. */
	private static String histFile=null;

	/** Literal FASTA input path for real-count cases. */
	private static String inputFile=null;

	/** Retained-read attachment output path for count output cases. */
	private static String outputFile=null;
}
