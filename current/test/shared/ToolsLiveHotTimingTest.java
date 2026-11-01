package test.shared;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import shared.Shared;
import shared.Tools;
import stream.SiteScore;

/** Matched timing for per-read/per-line Tools helpers with known boundaries.
 * Literal wrong-result cases live in ToolsLocationBoundaryTest; these workloads
 * use current-correct same-chromosome mapping cases and delimited searches.
 * @author Jean
 */
public final class ToolsLiveHotTimingTest{

	/** Runs mode=scores|delimited with positive rounds and warmup counts. */
	public static void main(String[] args){
		long rounds=10000000, warmup=1000000;
		String name="scores";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("rounds")){rounds=Long.parseLong(split[1]);}
			else if(split[0].equals("warmup")){warmup=Long.parseLong(split[1]);}
			else if(split[0].equals("mode")){name=split[1];}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}

		final int mode=name.equals("scores") ? 0 : name.equals("delimited") ? 1 : -1;
		if(mode<0 || rounds<1 || warmup<1){
			throw new IllegalArgumentException("Positive rounds/warmup and mode=scores|delimited required");
		}

		final ToolsLiveHotTimingTest test=new ToolsLiveHotTimingTest();
		final int operations=mode==0 ? test.scoreCases.length : test.delimitedCases.length;
		final long expected=mode==0 ? test.expectedScoreRound : test.expectedDelimitedRound;
		if(rounds>Long.MAX_VALUE/Math.max(1, operations)
				|| warmup>Long.MAX_VALUE/Math.max(1, operations)
				|| rounds>Long.MAX_VALUE/expected || warmup>Long.MAX_VALUE/expected){
			throw new IllegalArgumentException("Repetition count would overflow reported counters");
		}

		require(test.loop(warmup, mode)==expected*warmup, "warmup checksum");
		final long start=System.nanoTime();
		final long sum=test.loop(rounds, mode);
		final long elapsed=System.nanoTime()-start;
		require(sum==expected*rounds, "measured checksum");
		System.out.println("TOOLS_LIVE_HOT_TIMING mode="+name+" rounds="+rounds+
				" warmup="+warmup+" operations="+(rounds*operations)+
				" nanos="+elapsed+" checksum="+sum);
	}

	/** Dispatches once per timing mode, outside the measured inner loops. */
	private long loop(long rounds, int mode){
		return mode==0 ? loopScores(rounds) : loopDelimited(rounds);
	}

	/** Times only countTopScores calls on preallocated same-chromosome lists. */
	private long loopScores(long rounds){
		long sum=0;
		for(long round=0; round<rounds; round++){
			for(ScoreCase test : scoreCases){
				sum+=Tools.countTopScores(test.list, test.threshold);
			}
		}
		return sum;
	}

	/** Times only indexOfDelimited calls on preallocated ASCII byte arrays. */
	private long loopDelimited(long rounds){
		long sum=0;
		for(long round=0; round<rounds; round++){
			for(DelimitedCase test : delimitedCases){
				sum+=Tools.indexOfDelimited(test.ref, test.query, test.start, test.delimiter);
			}
		}
		return sum;
	}

	/** Prebuilds every input and independent per-round expected sum. */
	private ToolsLiveHotTimingTest(){
		scoreCases=new ScoreCase[] {
				scoreCase(2, 3, site(100, 109, 100), site(200, 209, 99), site(300, 309, 98)),
				scoreCase(5, 2, site(110, 119, 100), site(210, 219, 95), site(310, 319, 94)),
				scoreCase(3, 2, site(120, 129, 100), site(220, 229, 97), site(320, 329, 93)),
				scoreCase(5, 3, site(130, 139, 120), site(230, 239, 120),
						site(330, 339, 115), site(430, 439, 113)),
				scoreCase(5, 1, site(140, 149, 88), site(240, 249, 82), site(340, 349, 81)),
				scoreCase(3, 2, site(150, 159, 75), site(150, 159, 74), site(250, 259, 73))
		};
		delimitedCases=new DelimitedCase[] {
				delimitedCase("RNAME=chr1;POS=100;NM=0", "RNAME=", 0, 0),
				delimitedCase("AA;CHROM=1;TYP=B;END=Z", "TYP=", 0, 11),
				delimitedCase("XXTYP=A;TAG=B;TYP=C", "TYP=", 0, 14),
				delimitedCase("junk;TYP=X;END=Q", "TYP=", 5, 5),
				delimitedCase("AA;ZZ=1;TYP=2", "TYP=", 4, 8),
				delimitedCase("AS=1;XS=2;NM=3", "YT=", 0, -1),
				delimitedCase("XXTYP=A;TAG=B", "TYP=", 0, -1)
		};
		verifyScores(scoreCases);
		verifyDelimited(delimitedCases);
		expectedScoreRound=expectedScores(scoreCases);
		expectedDelimitedRound=expectedDelimited(delimitedCases);
	}

	/** Creates one preallocated countTopScores case. */
	private static ScoreCase scoreCase(int threshold, int expected, SiteScore... sites){
		final ArrayList<SiteScore> list=new ArrayList<SiteScore>(sites.length);
		for(SiteScore site : sites){list.add(site);}
		return new ScoreCase(list, threshold, expected);
	}

	/** Creates one same-chromosome plus-strand site with an initial score. */
	private static SiteScore site(int start, int stop, int score){
		return new SiteScore(1, Shared.PLUS, start, stop, 1, score);
	}

	/** Creates one preallocated indexOfDelimited case using semicolon delimiters. */
	private static DelimitedCase delimitedCase(String ref, String query, int start, int expected){
		return new DelimitedCase(ref.getBytes(StandardCharsets.US_ASCII), query, start, (byte)';', expected);
	}

	/** Sums independent expected score counts for one full corpus pass. */
	private static long expectedScores(ScoreCase[] tests){
		long sum=0;
		for(ScoreCase test : tests){sum+=test.expected;}
		return sum;
	}

	/** Sums independent expected delimiter indexes for one full corpus pass. */
	private static long expectedDelimited(DelimitedCase[] tests){
		long sum=0;
		for(DelimitedCase test : tests){sum+=test.expected;}
		return sum;
	}

	/** Checks every countTopScores case before aggregate checksum timing. */
	private static void verifyScores(ScoreCase[] tests){
		for(ScoreCase test : tests){
			final int count=Tools.countTopScores(test.list, test.threshold);
			require(count==test.expected, "countTopScores case: "+count+" != "+test.expected);
		}
	}

	/** Checks every indexOfDelimited case before aggregate checksum timing. */
	private static void verifyDelimited(DelimitedCase[] tests){
		for(DelimitedCase test : tests){
			final int index=Tools.indexOfDelimited(test.ref, test.query, test.start, test.delimiter);
			require(index==test.expected, "indexOfDelimited case: "+index+" != "+test.expected);
		}
	}

	/** Checks expected output independently of JVM assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One immutable mapping-count timing case. */
	private static final class ScoreCase{
		/** Stores a descending-score site list. */
		ScoreCase(ArrayList<SiteScore> list_, int threshold_, int expected_){
			list=list_;
			threshold=threshold_;
			expected=expected_;
		}
		/** Descending-score sites. */
		final ArrayList<SiteScore> list;
		/** Inclusive score distance from the best site. */
		final int threshold;
		/** Independent expected count. */
		final int expected;
	}

	/** One immutable delimited-search timing case. */
	private static final class DelimitedCase{
		/** Stores an ASCII byte array and its expected matching position. */
		DelimitedCase(byte[] ref_, String query_, int start_, byte delimiter_, int expected_){
			ref=ref_;
			query=query_;
			start=start_;
			delimiter=delimiter_;
			expected=expected_;
		}
		/** Bytes to search. */
		final byte[] ref;
		/** Token to find. */
		final String query;
		/** Search start offset. */
		final int start;
		/** Required preceding delimiter. */
		final byte delimiter;
		/** Independent expected matching position. */
		final int expected;
	}

	/** Prebuilt score-count workloads. */
	private final ScoreCase[] scoreCases;
	/** Prebuilt delimited-search workloads. */
	private final DelimitedCase[] delimitedCases;
	/** Expected sum for one score-count corpus pass. */
	private final long expectedScoreRound;
	/** Expected sum for one delimited-search corpus pass. */
	private final long expectedDelimitedRound;
}
