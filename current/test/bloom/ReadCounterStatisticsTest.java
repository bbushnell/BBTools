package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import bloom.KCountArray;
import bloom.KCountArray3;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;

/** Checks finite ReadCounter statistics and long histogram aggregation.
 * The overflow adapter is a formatter-boundary sentinel, not a giant real-table
 * reproduction: it returns a fixed 100000-bin histogram through KCountArray.
 * Raw reports are retained for every case so baseline failures stay inspectable.
 * @author Brian Bushnell, Jean */
public final class ReadCounterStatisticsTest{

	/** Takes one new scratch directory and records all three report cases. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Files.createDirectory(Paths.get(args[0]));
		int failures=0;
		failures+=emptyActualCounter(root);
		failures+=nonemptyActualCounter(root);
		failures+=overflowFormatterBoundary(root);
		if(failures!=0){throw new AssertionError("ReadCounter statistics failed "+failures+" of 3 cases.");}
		System.out.println("READ_COUNTER_STATISTICS_TEST_PASS cases=3");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Empty occupied cells must preserve shared keysCounted and print finite zero ratios. */
	private static int emptyActualCounter(final Path root) throws Exception{
		final KCountArray counts=KCountArray.makeNew(2048, 16, 1);
		final String report=captureStatistics(counts, 11);
		Files.write(root.resolve("empty-actual-counter.txt"), report.getBytes(StandardCharsets.UTF_8));
		return check("empty-actual-counter", report, new String[]{
				"Keys Counted:  \t         \t11",
				"Unique:        \t         \t0",
				"Avg Sites/Key: \t         \t0.000",
				"Singleton:     \t0.000%   \t0",
				"Useful:        \t0.000%   \t0"
		});
	}

	/** One occupied real cell should keep the existing nonempty ratios unchanged. */
	private static int nonemptyActualCounter(final Path root) throws Exception{
		final KCountArray counts=KCountArray.makeNew(2048, 16, 1);
		counts.increment(42, 3);
		final String report=captureStatistics(counts, 3);
		Files.write(root.resolve("nonempty-actual-counter.txt"), report.getBytes(StandardCharsets.UTF_8));
		return check("nonempty-actual-counter", report, new String[]{
				"Keys Counted:  \t         \t3",
				"Unique:        \t         \t1",
				"Avg Sites/Key: \t         \t3.000",
				"Singleton:     \t0.000%   \t0",
				"Useful:        \t100.000%   \t1"
		});
	}

	/** A long-valued histogram bucket must not overflow while formatting the 8-15 row. */
	private static int overflowFormatterBoundary(final Path root) throws Exception{
		final long[] freq=new long[100000];
		freq[0]=1;
		freq[8]=2147483648L;
		freq[9]=7;
		final KCountArray counts=new FormatterBoundaryKCountArray(freq);
		final String report=captureStatistics(counts, 2147483655L);
		Files.write(root.resolve("overflow-formatter-boundary.txt"), report.getBytes(StandardCharsets.UTF_8));
		return check("overflow-formatter-boundary", report, new String[]{
				"8-15    \t100.000%   \t2147483655",
				"Keys Counted:  \t         \t2147483655",
				"Unique:        \t         \t2147483655"
		});
	}

	/*--------------------------------------------------------------*/
	/*----------------      Capture and Checks      ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures printStatistics while restoring global output and shared keys. */
	private static String captureStatistics(final KCountArray counts, final long keysCounted) throws Exception{
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream prior=System.out;
		final long oldKeys=KmerCountAbstract.keysCounted;
		try{
			KmerCountAbstract.keysCounted=keysCounted;
			try(PrintStream capture=new PrintStream(buffer, true, "UTF-8")){
				System.setOut(capture);
				new ReadCounter(1).printStatistics(counts);
			}
		}finally{
			System.setOut(prior);
			KmerCountAbstract.keysCounted=oldKeys;
		}
		return buffer.toString("UTF-8");
	}

	/** Reports all missing literal fixed-oracle lines instead of stopping on the first one. */
	private static int check(final String name, final String report, final String[] expectedLines){
		int missing=0;
		for(String line : expectedLines){
			if(!hasLine(report, line)){
				missing++;
				System.out.println("case="+name+" missing="+line.replace('\t', '|'));
			}
		}
		final boolean pass=missing==0 && finite(report);
		System.out.println("case="+name+" finite="+finite(report)+" missing="+missing+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Matches complete output lines so a longer incorrect count cannot pass by prefix. */
	private static boolean hasLine(final String report, final String expected){
		for(String line : report.split("\\r?\\n")){
			if(trimRight(line).equals(trimRight(expected))){return true;}
		}
		return false;
	}

	/** Removes only the formatter's trailing padding. */
	private static String trimRight(final String line){
		int end=line.length();
		while(end>0 && Character.isWhitespace(line.charAt(end-1))){end--;}
		return line.substring(0, end);
	}

	/** The fixed report must not contain an undefined or infinite percentage. */
	private static boolean finite(final String report){
		return !report.contains("NaN") && !report.contains("Infinity");
	}

	/*--------------------------------------------------------------*/
	/*----------------     Histogram Test Adapter   ----------------*/
	/*--------------------------------------------------------------*/

	/** Supplies a fixed histogram at the public printStatistics boundary. */
	private static final class FormatterBoundaryKCountArray extends KCountArray3{

		FormatterBoundaryKCountArray(final long[] freq_){
			super(2048, 16);
			freq=freq_;
		}

		@Override
		public long[] transformToFrequency(){
			return freq;
		}

		private final long[] freq;
	}
}
