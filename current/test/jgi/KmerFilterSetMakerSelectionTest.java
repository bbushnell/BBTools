package test.jgi;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import jgi.KmerFilterSetMaker;

/**
 * Literal finite checks for KmerFilterSetMaker's histogram cutoff and quota selector.
 * Cases use small histograms with synthetic boundary counts rather than giant inputs.
 * @author Jean
 */
public final class KmerFilterSetMakerSelectionTest{

	/** Runs every literal selection fixture and reports observed threshold/quota triples. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final KmerFilterSetMakerSelectionTest test=
				new KmerFilterSetMakerSelectionTest(selectKmersMethod());
		for(final Case c : CASES){test.check(c);}
		System.out.println("KmerFilterSetMakerSelectionTest cases="+test.cases+
				" failures="+test.failures);
		if(test.failures!=0){
			throw new AssertionError("KmerFilterSetMaker selection failures="+test.failures);
		}
	}

	/** Stores the reflected production helper and its immutable result fields. */
	private KmerFilterSetMakerSelectionTest(final Method selectKmers_) throws Exception{
		selectKmers=selectKmers_;
		final Class<?> selection=Class.forName("jgi.KmerFilterSetMaker$Selection");
		minCountToKeep=selection.getDeclaredField("minCountToKeep");
		numGood=selection.getDeclaredField("numGood");
		numKept=selection.getDeclaredField("numKept");
		minCountToKeep.setAccessible(true);
		numGood.setAccessible(true);
		numKept.setAccessible(true);
	}

	/** Invokes the production helper once and compares the literal oracle. */
	private void check(final Case c) throws Exception{
		cases++;
		final Object observed=select(c);
		final int observedMinCountToKeep=minCountToKeep.getInt(observed);
		final long observedNumGood=numGood.getLong(observed);
		final int observedNumKept=numKept.getInt(observed);
		System.out.println(c.name+"\tminCountToKeep="+observedMinCountToKeep+
				"\tnumGood="+observedNumGood+"\tnumKept="+observedNumKept);
		if(observedMinCountToKeep!=c.minCountToKeep || observedNumGood!=c.numGood ||
				observedNumKept!=c.numKept){
			failures++;
			System.err.println("FAIL "+c.name+": expected minCountToKeep="+
					c.minCountToKeep+", numGood="+c.numGood+", numKept="+c.numKept);
		}
	}

	/** Calls the private static helper, preserving its underlying failure cause. */
	private Object select(final Case c) throws Exception{
		try{
			return selectKmers.invoke(null, c.counts, c.max, c.minKmersPerIteration,
					c.maxKmersPerIteration, c.readsIn);
		}catch(final InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Finds the private production helper without changing production accessibility. */
	private static Method selectKmersMethod() throws NoSuchMethodException{
		final Method method=KmerFilterSetMaker.class.getDeclaredMethod("selectKmers",
				long[].class, int.class, int.class, int.class, long.class);
		method.setAccessible(true);
		return method;
	}

	/** Literal histogram and expected selection tuple. */
	private static final class Case{

		Case(final String name_, final long[] counts_, final int max_,
				final int minKmersPerIteration_, final int maxKmersPerIteration_,
				final long readsIn_, final int minCountToKeep_, final long numGood_,
				final int numKept_){
			name=name_;
			counts=counts_;
			max=max_;
			minKmersPerIteration=minKmersPerIteration_;
			maxKmersPerIteration=maxKmersPerIteration_;
			readsIn=readsIn_;
			minCountToKeep=minCountToKeep_;
			numGood=numGood_;
			numKept=numKept_;
		}

		final String name;
		final long[] counts;
		final int max;
		final int minKmersPerIteration;
		final int maxKmersPerIteration;
		final long readsIn;
		final int minCountToKeep;
		final long numGood;
		final int numKept;
	}

	/** Private production helper under test. */
	private final Method selectKmers;
	/** Immutable result field carrying the selected minimum count. */
	private final Field minCountToKeep;
	/** Immutable result field carrying the selected high-kmer count. */
	private final Field numGood;
	/** Immutable result field carrying the selected dump quota. */
	private final Field numKept;
	/** Number of finite fixtures executed. */
	private int cases=0;
	/** Number of fixture results that missed their literal oracle. */
	private int failures=0;

	private static final Case[] CASES={
		new Case("all-zero", new long[]{0, 0, 0, 0}, 0,
				1, 2, 10, -1, 0, 0),
		new Case("ordinary-high-bin", new long[]{0, 0, 0, 0, 3}, 4,
				2, 2, 10, 4, 3, 2),
		new Case("descending-cutoff", new long[]{0, 0, 0, 0, 0, 2, 1}, 6,
				3, 9, 9, 5, 3, 3),
		new Case("exclude-singletons", new long[]{0, 999, 1}, 2,
				10, 10, 10, 2, 1, 1),
		new Case("cap-to-reads", new long[]{0, 0, 0, 0, 0, 10}, 5,
				1, 3, 2, 5, 10, 2),
		new Case("single-bin-int-min", new long[]{0, 0, 0, 0, 0, 1L<<31}, 5,
				1, 7, 7, 5, 1L<<31, 7),
		new Case("single-bin-wrap-one", new long[]{0, 0, 0, 0, 0, (1L<<32)+1}, 5,
				1, 7, 7, 5, (1L<<32)+1, 7),
		new Case("cumulative-overflow", new long[]{0, 0, 0, 0, 0, 2, Integer.MAX_VALUE-1L}, 6,
				Integer.MAX_VALUE, 9, 9, 5, Integer.MAX_VALUE+1L, 9)
	};
}
