package test.jgi;

import jgi.ReplaceHeaders;
import shared.KillSwitch;

/** Observes ReplaceHeaders.main with bounded failure containment in a fresh JVM.
 * An escaping cause is printed before immediate whole-process failure; this does
 * not demonstrate production cleanup or flushing. Successful returns stay natural,
 * so surviving production threads still trigger the armed daemon watchdog.
 * @author Sucrose
 */
public final class ReplaceHeadersMissingHeaderTest{

	/** Forwards unchanged CLI arguments after arming the failure-only watchdog.
	 * @param args Literal ReplaceHeaders arguments supplied by the dev launcher
	 */
	public static void main(final String[] args){
		assert(args!=null) : "The JVM supplies a nonnull argument array; forward it unchanged to ReplaceHeaders.";
		final Thread watchdog=new Thread(new Runnable(){
			@Override
			public void run(){watchDeadline();}
		}, "replaceheaders-missing-header-watchdog");
		watchdog.setDaemon(true);
		watchdog.start();
		try{
			ReplaceHeaders.main(args);
			System.out.println("TEST_MAIN_RETURNED");
		}catch(Throwable cause){
			cause.printStackTrace(System.err);
			System.err.println("TEST_GUARD_CAUGHT: ReplaceHeaders.main returned exceptionally; this guard does not establish cleanup.");
			System.err.flush();
			KillSwitch.kill("Test containment halts this JVM after the original cause; no cleanup or flush guarantee.");
		}
	}

	/** Fails this test JVM if it survives the deadline, including after main returns.
	 * An interrupted watchdog also invalidates the test; it never permits success.
	 */
	private static void watchDeadline(){
		try{
			Thread.sleep(WATCHDOG_MILLIS);
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
		}
		System.err.println("TEST_GUARD_TIMEOUT: deadline reached or watchdog interrupted; discard this test result.");
		System.err.flush();
		Runtime.getRuntime().halt(124);
	}

	/** Utility entry point; no fixture instances are needed. */
	private ReplaceHeadersMissingHeaderTest(){}

	/** Failure-containment deadline, not a performance or natural-teardown claim. */
	private static final long WATCHDOG_MILLIS=30000;
}
