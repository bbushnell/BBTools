package test.structures;

import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.SuperLongList;

/**
 * Literal percentile oracles, including exactly representable totals that round
 * upward or downward when converted to double. Huge multiplicities use zero bins;
 * no test expands them or enters a per-sample loop for a huge selected bin.
 * @author Shinobu
 */
public final class SuperLongListPercentileTest{

	/** Runs qualified old failures with mode=baseline, or repaired results with mode=fixed. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SuperLongListPercentileTest.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses standard flags and runs only finite literal checks. */
	private static void run(final String[] args){
		boolean baseline=false, configured=false;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(arg.equals("mode=baseline")){baseline=true; configured=true;}
			else if(arg.equals("mode=fixed")){baseline=false; configured=true;}
			else if(!Parser.parseStatic(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		check(configured, "Specify mode=baseline or mode=fixed.");
		boolean assertions=false;
		assert(assertions=true);
		check(assertions, "-ea is required to qualify the old endpoint failures.");
		final SuperLongList single=list(Long.MAX_VALUE-1);
		endpoint("single-large-value", single, 1, Long.MAX_VALUE-1, baseline);
		endpoint("single-large-count", single, 3, 1, baseline);
		final SuperLongList mixed=list(1, Long.MAX_VALUE-2);
		endpoint("mixed-large-value", mixed, 1, Long.MAX_VALUE-2, baseline);
		endpoint("mixed-large-count", mixed, 3, 2, baseline);
		final SuperLongList hugeCount=new SuperLongList(2);
		hugeCount.increment(0, Long.MAX_VALUE-2);
		hugeCount.add(1);
		//The final selected bin contains one value: even the fixed sum walk needs only one iteration.
		endpoint("large-count-value", hugeCount, 0, 1, baseline);
		endpoint("large-count-sum", hugeCount, 2, 1, baseline);
		final SuperLongList roundedDown=new SuperLongList(2);
		roundedDown.increment(0, 9007199254740992L);
		roundedDown.add(1);
		check(roundedDown.percentileValueByCount(1)==(baseline ? 0 : 1), "Rounded-down count target selected the wrong bin.");
		System.out.println("PASS rounded-down-count value="+(baseline ? 0 : 1));
		final SuperLongList ordinary=list(1, 1, 8, 90);
		values("zero", ordinary, 0, 0, 0, 0, 0);
		values("tiny", ordinary, Double.MIN_VALUE, 0, 0, 0, 0);
		values("lower-tail", ordinary, 0.1, 0, 8, 0, 3);
		values("half", ordinary, 0.5, 1, 90, 2, 4);
		values("full", ordinary, 1, 90, 90, 100, 4);
		values("dense-half", list(1, 2, 3), 0.5, 1, 2, 1, 2);
		values("dense-full", list(1, 2, 3), 1, 3, 3, 6, 3);
		for(final double p : new double[]{0, 0.5, 1}){
			values("empty-"+p, list(), p, 0, 0, 0, 0);
			values("all-zero-"+p, list(0, 0, 0), p, 0, 0, 0, 0);
		}
		check(ordinary.count()==4 && ordinary.sum()==100, "Percentile queries changed histogram totals.");
		System.out.println("SUPER_LONG_LIST_PERCENTILE mode="+(baseline ? "baseline" : "fixed")+" failed=0");
	}

	/** Qualifies the exact terminal assertion in old code, not a timeout or unrelated exception. */
	private static void endpoint(final String label, final SuperLongList list,
			final int method, final long expected, final boolean baseline){
		AssertionError failure=null;
		long actual=-1;
		try{actual=query(list, method, 1);}
		catch(final AssertionError e){failure=e;}
		if(baseline){
			check(failure!=null && ("1.0, "+list.count()+", "+list.sum()).equals(failure.getMessage()),
					label+": expected the rounded target to exceed its exact total, got "+failure);
			check(failure.getStackTrace()[0].getClassName().equals("structures.SuperLongList"),
					label+": failure did not originate in SuperLongList.");
		}else{check(failure==null && actual==expected, label+": expected "+expected+", got "+actual+", "+failure);}
		System.out.println("PASS "+label+" expected="+expected+" old_assertion="+baseline);
	}

	/** Checks all four APIs against independent literal values. */
	private static void values(final String label, final SuperLongList list, final double p,
			final long valueCount, final long valueSum, final long sumCount, final long countSum){
		final long[] expected={valueCount, valueSum, sumCount, countSum};
		for(int i=0; i<expected.length; i++){
			final long actual=query(list, i, p);
			check(actual==expected[i], label+": method="+i+", expected="+expected[i]+", actual="+actual);
		}
		System.out.println("PASS "+label+" four_methods");
	}

	/** Dispatches only the four public percentile methods under test. */
	private static long query(final SuperLongList list, final int method, final double p){
		switch(method){
			case 0: return list.percentileValueByCount(p);
			case 1: return list.percentileValueBySum(p);
			case 2: return list.percentileSumByCount(p);
			case 3: return list.percentileCountBySum(p);
			default: throw new IllegalArgumentException("Unknown method: "+method);
		}
	}

	/** Uses a small histogram plus sorted overflow values. */
	private static SuperLongList list(final long... values){
		final SuperLongList list=new SuperLongList(4);
		for(final long value : values){list.add(value);}
		list.sort();
		return list;
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private SuperLongListPercentileTest(){}
}
