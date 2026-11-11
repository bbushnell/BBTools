package test.assemble;

import java.lang.reflect.Constructor;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Checks finite signs and ordering laws on TadpoleWrapper's private assembly snapshots.
 * Reflection constructs Record instances only; comparisons use Comparable directly.
 * Uniform-contig snapshots have L50=L90=max=length and total=length*count.
 * Fixed mode (the default) uses literal expectations for 64-base bins and a small
 * deterministic corpus for ordering laws; it does not implement a reference comparator.
 * Explicit baseline mode retains the original observations, including known defects,
 * and is intended only for the original production comparator.
 * No assembly, input files, statistics calculation or sorting is performed.
 * Run in an isolated JVM with assertions enabled through the development launcher.
 * @author Brian Bushnell, Sucrose
 */
public final class TadpoleWrapperComparatorTest{

	/** Accepts mode=fixed (default) or mode=baseline, plus normal launcher JVM flags. */
	public static void main(final String[] args) throws ReflectiveOperationException{
		final PreParser pp=new PreParser(args, TadpoleWrapperComparatorTest.class, false);
		try{
			boolean baseline=false;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2 || !pair[0].equalsIgnoreCase("mode") ||
						(!pair[1].equalsIgnoreCase("baseline") && !pair[1].equalsIgnoreCase("fixed"))){
					throw new IllegalArgumentException("Expected mode=fixed or mode=baseline: "+arg);
				}
				baseline=pair[1].equalsIgnoreCase("baseline");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Comparator observations require -ea, as supplied by the launcher.");}
			if(baseline){runBaseline();}else{runFixed();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private TadpoleWrapperComparatorTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercises count/tolerance witnesses, a cycle, equal snapshots and a clear length preference. */
	private static void runBaseline() throws ReflectiveOperationException{
		final Class<?> record=Class.forName("assemble.TadpoleWrapper$Record");
		final Constructor<?> constructor=record.getDeclaredConstructor(int.class, long.class, long.class,
				long.class, long.class, long.class);
		constructor.setAccessible(true);
		final Comparable<Object> two=uniform(constructor, 31, 100, 2);
		final Comparable<Object> one=uniform(constructor, 41, 100, 1);
		final Comparable<Object> low=uniform(constructor, 31, 250000, 1);
		final Comparable<Object> high=uniform(constructor, 41, 250501, 1);
		final Comparable<Object> a=uniform(constructor, 31, 10180, 3);
		final Comparable<Object> b=uniform(constructor, 41, 10090, 2);
		final Comparable<Object> c=uniform(constructor, 51, 10000, 1);
		final Comparable<Object> equalA=uniform(constructor, 31, 10180, 3);
		final Comparable<Object> longer=uniform(constructor, 61, 10200, 4);
		final ByteBuilder out=new ByteBuilder(512);
		out.append("case\texpected_sign\tobserved_sign\tpassed\n");
		int checks=0;
		checks+=observe(out, "two100_vs_one100", two, one, -1);
		checks+=observe(out, "one100_vs_two100", one, two, -1);
		checks+=observe(out, "tolerance_high_vs_low", high, low, 1);
		checks+=observe(out, "tolerance_low_vs_high", low, high, 1);
		checks+=observe(out, "cycle_a_vs_b", a, b, -1);
		checks+=observe(out, "cycle_b_vs_c", b, c, -1);
		checks+=observe(out, "cycle_c_vs_a", c, a, -1);
		checks+=observe(out, "self_a", a, a, 0);
		checks+=observe(out, "equal_metrics_same_k", a, equalA, 0);
		checks+=observe(out, "longer_by_two_percent", longer, c, 1);
		checks+=observe(out, "shorter_by_two_percent", c, longer, -1);
		System.out.print(out);
		System.out.println("TADPOLE_WRAPPER_COMPARATOR_BASELINE checks="+checks+" failed=0");
	}

	/** Literal policy cases; raw values and total bases never compute the expected signs. */
	private static void runFixed() throws ReflectiveOperationException{
		final Class<?> record=Class.forName("assemble.TadpoleWrapper$Record");
		final Constructor<?> constructor=record.getDeclaredConstructor(int.class, long.class, long.class,
				long.class, long.class, long.class);
		constructor.setAccessible(true);
		final Comparable<Object> two=uniform(constructor, 31, 100, 2), one=uniform(constructor, 41, 100, 1);
		final Comparable<Object> low=uniform(constructor, 31, 250000, 1), high=uniform(constructor, 41, 250501, 1);
		final Comparable<Object> a=uniform(constructor, 31, 10180, 3), b=uniform(constructor, 41, 10090, 2);
		final Comparable<Object> c=uniform(constructor, 51, 10000, 1), equalA=uniform(constructor, 31, 10180, 3);
		final Comparable<Object> longer=uniform(constructor, 61, 10200, 4);
		final ByteBuilder out=new ByteBuilder(2048);
		out.append("case\texpected_sign\tobserved_sign\tpassed\n");
		int checks=0;
		checks+=fixedSign(out, "two100_vs_one100", two, one, -1);
		checks+=fixedSign(out, "one100_vs_two100", one, two, 1);
		checks+=fixedSign(out, "tolerance_high_vs_low", high, low, 1);
		checks+=fixedSign(out, "tolerance_low_vs_high", low, high, -1);
		checks+=fixedSign(out, "cycle_a_vs_b", a, b, 1);
		checks+=fixedSign(out, "cycle_b_vs_c", b, c, 1);
		checks+=fixedSign(out, "cycle_c_vs_a", c, a, -1);
		checks+=fixedSign(out, "self_a", a, a, 0);
		checks+=fixedSign(out, "equal_metrics_same_k", a, equalA, 0);
		checks+=fixedSign(out, "longer_by_two_percent", longer, c, 1);
		checks+=fixedSign(out, "shorter_by_two_percent", c, longer, -1);

		final Comparable<Object> bin63=uniform(constructor, 31, 63, 1), bin64=uniform(constructor, 31, 64, 1);
		final Comparable<Object> bin127=uniform(constructor, 31, 127, 1), bin128=uniform(constructor, 31, 128, 1);
		checks+=fixedSign(out, "boundary_63_64", bin63, bin64, -1);
		checks+=fixedSign(out, "boundary_127_128", bin127, bin128, -1);
		checks+=fixedSign(out, "within_bin_all_lengths", bin64, bin127, 0);
		checks+=fixedSign(out, "within_bin_l50", snapshot(constructor, 31, 128, 64, 512, 4, 256),
				snapshot(constructor, 31, 191, 64, 512, 4, 256), 0);
		checks+=fixedSign(out, "within_bin_l90", snapshot(constructor, 31, 192, 64, 768, 4, 256),
				snapshot(constructor, 31, 192, 127, 768, 4, 256), 0);
		checks+=fixedSign(out, "within_bin_max", snapshot(constructor, 31, 128, 64, 640, 4, 256),
				snapshot(constructor, 31, 128, 64, 640, 4, 319), 0);
		// Each winning primary key faces worse lower-priority keys, count and K.
		final Comparable<Object> l50First=snapshot(constructor, 61, 192, 64, 700, 5, 256);
		final Comparable<Object> l50Other=snapshot(constructor, 31, 191, 128, 700, 4, 320);
		final Comparable<Object> l90First=snapshot(constructor, 61, 256, 128, 760, 5, 320);
		final Comparable<Object> l90Other=snapshot(constructor, 31, 256, 127, 780, 4, 384);
		final Comparable<Object> maxFirst=snapshot(constructor, 61, 128, 64, 650, 5, 320);
		final Comparable<Object> maxOther=snapshot(constructor, 31, 128, 64, 639, 4, 319);
		checks+=fixedSign(out, "l50_priority", l50First, l50Other, 1);
		checks+=fixedSign(out, "l90_priority", l90First, l90Other, 1);
		checks+=fixedSign(out, "max_priority", maxFirst, maxOther, 1);
		final Comparable<Object> countOther=uniform(constructor, 31, 127, 2);
		final Comparable<Object> kOther=uniform(constructor, 41, 127, 1);
		checks+=fixedSign(out, "count_after_equal_bins", countOther, uniform(constructor, 61, 64, 1), -1);
		checks+=fixedSign(out, "k_after_equal_bins", kOther, bin64, -1);
		checks+=fixedSign(out, "total_bases_ignored", l50First, snapshot(constructor, 61, 192, 64, 701, 5, 256), 0);
		final Comparable<Object> top=uniform(constructor, 31, Long.MAX_VALUE, 1);
		final Comparable<Object> topFloor=uniform(constructor, 31, 9223372036854775744L, 1);
		checks+=fixedSign(out, "long_max_same_bin", top, topFloor, 0);
		checks+=fixedSign(out, "long_max_previous_bin", top, uniform(constructor, 31, 9223372036854775743L, 1), 1);
		final Comparable<Object> empty31=snapshot(constructor, 31, 0, 0, 0, 0, 0);
		final Comparable<Object> empty41=snapshot(constructor, 41, 0, 0, 0, 0, 0);
		checks+=fixedSign(out, "empty_vs_short", empty31, bin63, -1);
		checks+=fixedSign(out, "short_vs_empty", bin63, empty31, 1);
		checks+=fixedSign(out, "emptyK31_vs_emptyK41", empty31, empty41, 1);
		System.out.print(out);
		checkLaws(bin63, bin64, bin127, bin128, l50First, l50Other, l90First, l90Other,
				maxFirst, maxOther, countOther, kOther, top, topFloor, empty31, empty41);
		System.out.println("TADPOLE_WRAPPER_COMPARATOR_FIXED checks="+checks+" failed=0");
	}

	/**
	 * Constructs synthetic scalar statistics for key-priority checks. These records
	 * satisfy positive scalar/length-order bounds or have all-zero empty metrics;
	 * they are not generated assemblies.
	 * @param constructor Accessible private Record constructor
	 * @param k Positive trial identifier; no normalization is performed
	 * @param l50 Length at the 50-percent threshold
	 * @param l90 Length at the 90-percent threshold
	 * @param total Total bases, deliberately varied in one ignored-field control
	 * @param count Positive contig count, or zero with all other statistics zero
	 * @param maximum Maximum contig length
	 * @return Snapshot compared through Comparable, never reflective method invocation
	 */
	@SuppressWarnings("unchecked")
	private static Comparable<Object> snapshot(final Constructor<?> constructor, final int k, final long l50,
			final long l90, final long total, final long count, final long maximum) throws ReflectiveOperationException{
		assert(constructor!=null && k>0 && ((count==0 && l50==0 && l90==0 && total==0 && maximum==0) ||
				(count>0 && l90>0 && l90<=l50 && l50<=maximum && maximum<=total))) :
			"Synthetic ordering cases require positive ordered lengths or consistently empty statistics; k="+k+" L50="+l50+" L90="+l90+" max="+maximum+" total="+total+" count="+count;
		return (Comparable<Object>)constructor.newInstance(k, l50, l90, total, count, maximum);
	}

	/** Compares actual signs with caller-supplied literals, without computing quantized keys. */
	private static int fixedSign(final ByteBuilder out, final String name, final Comparable<Object> left,
			final Comparable<Object> right, final int expected){
		assert(left!=null && right!=null && expected>=-1 && expected<=1) :
			"Fixed sign checks require constructed snapshots and a literal sign: "+name;
		final int observed=Integer.signum(left.compareTo(right));
		out.append(name).tab().append(expected).tab().append(observed).tab().append(observed==expected ? 1 : 0).nl();
		if(observed!=expected){
			System.out.print(out);
			throw new AssertionError("Comparator fixed expectation failed: "+name+" expected="+expected+" observed="+observed);
		}
		return 1;
	}

	/**
	 * Checks every pair/triple in a fixed small corpus using observed signs only.
	 * Zero-equivalent records must compare identically against every third record
	 * in both directions. The corpus is borrowed and neither stored nor modified.
	 * @param records Fixed snapshots containing boundary, priority, tie and extreme cases
	 */
	@SafeVarargs
	private static void checkLaws(final Comparable<Object>... records){
		assert(records.length>1) : "Comparator laws require distinct snapshots, including zero-equivalent keys.";
		final int n=records.length;
		final int[][] signs=new int[n][n];
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){signs[i][j]=Integer.signum(records[i].compareTo(records[j]));}
		}
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				final int ij=signs[i][j];
				if(ij!=-signs[j][i]){throw new AssertionError("Sign antisymmetry failed: i="+i+" j="+j);}
				for(int k=0; k<n; k++){
					final int jk=signs[j][k], ik=signs[i][k];
					if((ij>0 && jk>0 && ik<=0) || (ij<0 && jk<0 && ik>=0)){
						throw new AssertionError("Sign transitivity failed: i="+i+" j="+j+" k="+k);
					}
					if(ij==0 && (ik!=jk || signs[k][i]!=signs[k][j])){
						throw new AssertionError("Zero-equivalence consistency failed: i="+i+" j="+j+" k="+k);
					}
				}
			}
		}
		System.out.println("TADPOLE_WRAPPER_COMPARATOR_LAWS records="+n+" pairs="+(n*n)+" triples="+(n*n*n)+" failed=0");
	}

	/**
	 * Constructs a private Record for positive, equal-length contigs without invoking assembly.
	 * @param constructor Actual six-scalar Record constructor, made accessible by the caller
	 * @param k Positive trial identifier; no K normalization is performed
	 * @param length Length of every contig, in bases; also L50, L90 and maximum
	 * @param count Number of contigs; total bases are computed with checked multiplication
	 * @return Snapshot viewed only through its public Comparable interface
	 */
	@SuppressWarnings("unchecked")
	private static Comparable<Object> uniform(final Constructor<?> constructor, final int k,
			final long length, final long count) throws ReflectiveOperationException{
		assert(constructor!=null && k>0 && length>0 && count>0) :
			"Uniform-contig witnesses need a Record constructor and positive K/length/count; k="+k+" length="+length+" count="+count;
		return (Comparable<Object>)constructor.newInstance(k, length, length, Math.multiplyExact(length, count), count, length);
	}

	/** Appends the actual sign and rejects any departure from the literal original-code expectation. */
	private static int observe(final ByteBuilder out, final String name, final Comparable<Object> left,
			final Comparable<Object> right, final int expected){
		assert(left!=null && right!=null && expected>=-1 && expected<=1) :
			"Baseline sign observations compare constructed snapshots against a literal sign: "+name+" expected="+expected;
		final int observed=Integer.signum(left.compareTo(right));
		out.append(name).tab().append(expected).tab().append(observed).tab().append(observed==expected ? 1 : 0).nl();
		if(observed!=expected){
			System.out.print(out);
			throw new AssertionError("Comparator baseline changed: "+name+" expected sign="+expected+" observed="+observed);
		}
		return 1;
	}
}
