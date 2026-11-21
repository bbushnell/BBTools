package test.jgi;

import jgi.CalcUniqueness;

/** Observes CalcUniqueness construction without running k-mer processing.
 * Forwards arguments unchanged in a fresh JVM; the production constructor may
 * inspect input metadata and allocate k-mer tables. Successful construction prints
 * CONSTRUCTED; construction failures propagate unchanged.
 * @author Shinobu
 */
public final class CalcUniquenessConstructorTest{

	/** Invokes only the constructor, never CalcUniqueness.main or process.
	 * @param args Literal CalcUniqueness arguments supplied by the dev launcher
	 */
	public static void main(final String[] args){
		assert(args!=null) : "The JVM supplies a nonnull argument array; forward it unchanged to CalcUniqueness.";
		new CalcUniqueness(args);
		System.out.println("CONSTRUCTED");
	}

	/** Utility entry point; no fixture instances are needed. */
	private CalcUniquenessConstructorTest(){}
}
