package test.consensus;

import java.lang.reflect.Method;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLongArray;

import shared.Tools;

/** Exact insert-distribution ranks used by both consensus scaffolding tools. @author Shinobu */
public final class ConsensusHistogramTest{

	/** Baseline mode describes the old helper; fixed mode checks the intended quantiles. */
	public static void main(final String[] args) throws Exception{
		boolean baseline=false;
		for(final String arg : args){
			if(arg.equals("phase=baseline")){baseline=true;}
			else if(!arg.equals("phase=fixed")){throw new IllegalArgumentException(arg);}
		}
		final AtomicLongArray point=new AtomicLongArray(101);
		point.set(100, 11);
		final int[] expected=new int[1001];
		Arrays.fill(expected, baseline ? 0 : 100);
		check(Tools.makeHistogram(point, 1000), expected, "Eleven identical inserts are all length100.");
		point.set(100, 1001);
		Arrays.fill(expected, baseline ? 101 : 100);
		if(baseline){expected[0]=0;}
		check(Tools.makeHistogram(point, 1000), expected, "Returned lengths must be observed indices, not the next index.");
		final AtomicLongArray mixed=new AtomicLongArray(21);
		mixed.set(4, 3);
		mixed.set(10, 2);
		mixed.set(20, 2);
		check(Tools.makeHistogram(mixed, 3), baseline ? new int[]{0, 5, 5, 5} : new int[]{4, 4, 10, 20},
			"Four bin centers select ranks1,3,5,7 from seven expanded observations.");
		check(Tools.makeHistogram(mixed, 0), new int[]{baseline ? 5 : 10},
			"One bin selects rank4, the middle observation of seven.");
		final long large=1L<<53;
		final AtomicLongArray precise=new AtomicLongArray(new long[]{0, large, 0, large+2});
		check(Tools.makeHistogram(precise, 2), baseline ? new int[]{2, 4, 4} : new int[]{1, 3, 3},
			"The exact midpoint rank is 2^53+1, beyond the first occupied bin.");
		check(Tools.makeHistogram(new AtomicLongArray(0), 3), new int[4],
			"An empty distribution retains the existing zero-filled convention.");
		boolean rejected=false;
		try{Tools.makeHistogram(new AtomicLongArray(new long[]{-1}), 1);}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		if(rejected==baseline){throw new AssertionError("Negative observation counts must be rejected.");}
		checks++;
		rejected=false;
		try{Tools.makeHistogram(new AtomicLongArray(new long[]{Long.MAX_VALUE, 1}), 0);}
		catch(ArithmeticException expectedFailure){rejected=true;}
		if(rejected==baseline){throw new AssertionError("A total exceeding long range must be rejected.");}
		checks++;
		if(!baseline){extremeRanks();}
		System.out.println("CONSENSUS_HISTOGRAM_PASS phase="+(baseline ? "baseline" : "fixed")+" checks="+checks);
	}

	/** Cross-checks bounded integer arithmetic with an independent arbitrary-precision oracle. */
	private static void extremeRanks() throws Exception{
		final long half=Long.MAX_VALUE/2;
		check(Tools.makeHistogram(new AtomicLongArray(new long[]{0, half, 0, Long.MAX_VALUE-half}), 3),
			new int[]{1, 1, 3, 3}, "A representable Long.MAX_VALUE population retains its occupied values.");
		final Method rank=Tools.class.getDeclaredMethod("histogramRank", long.class, int.class, int.class);
		rank.setAccessible(true);
		for(final long total : new long[]{1, 2, 7, (1L<<53)+2, Long.MAX_VALUE}){
			for(final int slots : new int[]{1, 2, 1001, Integer.MAX_VALUE}){
				final BigInteger denominator=BigInteger.valueOf(2L*slots);
				for(final int index : new int[]{0, slots/2, slots-1}){
					final BigInteger numerator=BigInteger.valueOf(total).multiply(BigInteger.valueOf(2L*index+1));
					final long expected=numerator.add(denominator).subtract(BigInteger.ONE).divide(denominator).longValueExact();
					final long actual=(Long)rank.invoke(null, total, slots, index);
					if(actual!=expected){throw new AssertionError("Exact rank mismatch: total="+total+" slots="+slots+
						" index="+index+" expected="+expected+" actual="+actual);}
					checks++;
				}
			}
		}
	}

	/** Reports the first differing quantile rather than dumping a thousand bins. */
	private static void check(final int[] actual, final int[] expected, final String message){
		if(actual.length!=expected.length){throw new AssertionError(message+" Length="+actual.length);}
		for(int i=0; i<actual.length; i++){
			if(actual[i]!=expected[i]){throw new AssertionError(message+" Bin="+i+" expected="+expected[i]+" actual="+actual[i]);}
		}
		checks++;
	}

	private static int checks;
}
