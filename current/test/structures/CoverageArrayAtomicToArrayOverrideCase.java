package test.structures;

import java.util.Arrays;

import structures.CoverageArray2A;
import structures.CoverageArray3A;

/** Regression fixture for atomic median snapshots with overridden public toArray().
 * @author Sucrose
 */
public final class CoverageArrayAtomicToArrayOverrideCase{

	/** Runs descendant override checks for the atomic coverage backends. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=unchanged");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]) && "unchanged".equals(split[1]),
				"Expected expect=unchanged");
		check2A();
		check3A();
	}

	/** Checks CoverageArray2A descendants that expose retained char storage. */
	private static void check2A(){
		final ExposedCoverageArray2A ca=new ExposedCoverageArray2A();
		fill(ca, CA2A_DEPTHS);
		final char[] originalExposed=ca.exposed.clone();
		checkEquals("ca2a first median", CA2A_MEDIAN, ca.median());
		checkArray("ca2a depths after first median", ca, CA2A_DEPTHS);
		checkArray("ca2a exposed after first median", originalExposed, ca.exposed);
		checkBounds("ca2a after first median", ca);
		ca.setMedian(-1);
		checkEquals("ca2a second median", CA2A_MEDIAN, ca.median());
		checkArray("ca2a depths after second median", ca, CA2A_DEPTHS);
		checkArray("ca2a exposed after second median", originalExposed, ca.exposed);
		checkBounds("ca2a after second median", ca);
		pass("ca2a", CA2A_MEDIAN, ca.toString(), toString(ca.exposed));
	}

	/** Checks CoverageArray3A descendants that expose retained int storage. */
	private static void check3A(){
		final ExposedCoverageArray3A ca=new ExposedCoverageArray3A();
		fill(ca, CA3A_DEPTHS);
		final int[] originalExposed=ca.exposed.clone();
		checkEquals("ca3a first median", CA3A_MEDIAN, ca.median());
		checkArray("ca3a depths after first median", ca, CA3A_DEPTHS);
		checkArray("ca3a exposed after first median", originalExposed, ca.exposed);
		checkBounds("ca3a after first median", ca);
		ca.setMedian(-1);
		checkEquals("ca3a second median", CA3A_MEDIAN, ca.median());
		checkArray("ca3a depths after second median", ca, CA3A_DEPTHS);
		checkArray("ca3a exposed after second median", originalExposed, ca.exposed);
		checkBounds("ca3a after second median", ca);
		pass("ca3a", CA3A_MEDIAN, ca.toString(), Arrays.toString(ca.exposed));
	}

	/** Fills a 16-bit atomic coverage array with exact values. */
	private static void fill(final ExposedCoverageArray2A ca, final int[] values){
		for(int i=0; i<values.length; i++){
			ca.set(i, values[i]);
		}
	}

	/** Fills a 32-bit atomic coverage array with exact values. */
	private static void fill(final ExposedCoverageArray3A ca, final int[] values){
		for(int i=0; i<values.length; i++){
			ca.set(i, values[i]);
		}
	}

	/** Checks exact logical coverage values. */
	private static void checkArray(final String name, final CoverageArray2A ca,
			final int[] expected){
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], ca.get(i));
		}
	}

	/** Checks exact logical coverage values. */
	private static void checkArray(final String name, final CoverageArray3A ca,
			final int[] expected){
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], ca.get(i));
		}
	}

	/** Checks exact retained char array contents. */
	private static void checkArray(final String name, final char[] expected,
			final char[] observed){
		check(expected.length==observed.length,
				name+" expected length "+expected.length+", got "+observed.length);
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], observed[i]);
		}
	}

	/** Checks exact retained int array contents. */
	private static void checkArray(final String name, final int[] expected,
			final int[] observed){
		check(expected.length==observed.length,
				name+" expected length "+expected.length+", got "+observed.length);
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], observed[i]);
		}
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks exact touched bounds. */
	private static void checkBounds(final String name, final CoverageArray2A ca){
		checkEquals(name+" minIndex", 0, ca.minIndex);
		checkEquals(name+" maxIndex", CA2A_DEPTHS.length-1, ca.maxIndex);
	}

	/** Checks exact touched bounds. */
	private static void checkBounds(final String name, final CoverageArray3A ca){
		checkEquals(name+" minIndex", 0, ca.minIndex);
		checkEquals(name+" maxIndex", CA3A_DEPTHS.length-1, ca.maxIndex);
	}

	/** Formats char values numerically without embedding control characters. */
	private static String toString(final char[] array){
		final StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<array.length; i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array[i]);
		}
		sb.append(']');
		return sb.toString();
	}

	/** Prints one successful backend row. */
	private static void pass(final String name, final int median,
			final String depths, final String exposed){
		System.out.println("PASS\t"+name+'\t'+median+'\t'+depths+'\t'+exposed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** CoverageArray2A descendant with public toArray() exposing retained storage. */
	private static final class ExposedCoverageArray2A extends CoverageArray2A{

		/** Constructs a bounded atomic coverage array. */
		ExposedCoverageArray2A(){
			super(0, CA2A_DEPTHS.length);
		}

		/** Exposes the retained array median() must not sort. */
		@Override
		public char[] toArray(){return exposed;}

		final char[] exposed=CA2A_EXPOSED.clone();
	}

	/** CoverageArray3A descendant with public toArray() exposing retained storage. */
	private static final class ExposedCoverageArray3A extends CoverageArray3A{

		/** Constructs a bounded atomic coverage array. */
		ExposedCoverageArray3A(){
			super(0, CA3A_DEPTHS.length);
		}

		/** Exposes the retained array median() must not sort. */
		@Override
		public int[] toArray(){return exposed;}

		final int[] exposed=CA3A_EXPOSED.clone();
	}

	private static final int[] CA2A_DEPTHS={8, 1, 5, 2, 7};
	private static final char[] CA2A_EXPOSED={9, 3, 5, 1, 7};
	private static final int CA2A_MEDIAN=5;
	private static final int[] CA3A_DEPTHS={80, 10, 50, 20, 70};
	private static final int[] CA3A_EXPOSED={90, 30, 50, 10, 70};
	private static final int CA3A_MEDIAN=50;

}
