package test.jgi;

import jgi.RenameReads;

/** Observes RenameReads construction without invoking its processing entry points.
 * Forwards all arguments unchanged; use a fresh JVM and disposable ordinary files.
 * The production constructor still parses options and may inspect input metadata.
 * A normal return prints CONSTRUCTED; construction failures propagate unchanged.
 * @author Sucrose
 */
public final class RenameReadsConstructorTest{

	/** Constructs only; never calls RenameReads.main/process or starts read/output streams.
	 * @param args Literal RenameReads arguments supplied by the dev launcher
	 */
	public static void main(final String[] args){
		assert(args!=null) : "The JVM supplies a nonnull argument array; forward it unchanged to RenameReads.";
		new RenameReads(args);
		System.out.println("CONSTRUCTED");
	}

	/** Utility entry point; no fixture instances are needed. */
	private RenameReadsConstructorTest(){}
}
