package test.sort;

import java.io.File;

import sort.Shuffle;
import sort.Shuffle.ShuffleThread;

/** Checks worker failure visibility and slot release; extends Jean's SYN024 regression.
 * @author Jean, Shinobu
 */
public final class ShuffleThreadTest{

	/** Takes one valid FASTQ input and a fresh output path for the successful control. */
	public static void main(String[] args) throws InterruptedException{
		if(args.length!=2){throw new IllegalArgumentException("Expected input and fresh output paths");}
		Shuffle.setMaxThreads(1);
		check(new ShuffleThread(null, null, null, null, Shuffle.SHUFFLE, true), true, "missing-input");
		check(new ShuffleThread(null, null, null, null, 0, true), true, "invalid-mode");
		check(new ShuffleThread(args[0], null, args[1], null, Shuffle.SHUFFLE, false), false, "valid-input");
		if(!new File(args[1]).isFile()){throw new AssertionError("Successful worker produced no output");}
		Shuffle.waitForFinish();
	}

	/** Joins each worker and checks reported failure and capacity accounting independently. */
	private static void check(ShuffleThread thread, boolean expectedError, String name) throws InterruptedException{
		thread.start();
		thread.join();
		if(thread.errorState()!=expectedError){throw new AssertionError(name+": incorrect worker status");}
		if(Shuffle.addThread(0)!=0){throw new AssertionError(name+": terminated worker retained a slot");}
		System.out.println("PASS\t"+name+"\terror="+thread.errorState());
	}
}
