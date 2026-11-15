package test.pacbio;

import java.util.Arrays;

import pacbio.RemoveAdapters2;
import pacbio.RemoveAdapters3;
import stream.FASTQ;

/** Invokes a legacy adapter main with forced interleaving configured before input construction.
 * The production command lines do not accept an interleaved flag. This separate-process
 * test entry changes only the input pairing configuration and forwards every remaining argument.
 * @author Brian Bushnell, Sucrose
 */
public final class RemoveAdaptersForcedMainTest{

	/** Selects the legacy main and lets all production failures propagate to the process exit.
	 * @param args First program=2 or program=3, followed by unchanged production arguments
	 */
	public static void main(String[] args){
		if(args.length<2 || !(args[0].equals("program=2") || args[0].equals("program=3"))){
			throw new IllegalArgumentException("Use program=2|3 followed by legacy adapter arguments");
		}
		if(!RemoveAdapters2.class.desiredAssertionStatus() || !RemoveAdapters3.class.desiredAssertionStatus()){
			throw new AssertionError("The configured-main fixture requires -ea to retain production invariants");
		}
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=true;
		final String[] forwarded=Arrays.copyOfRange(args, 1, args.length);
		if(args[0].equals("program=2")){RemoveAdapters2.main(forwarded);}
		else{RemoveAdapters3.main(forwarded);}
	}

	private RemoveAdaptersForcedMainTest(){}
}
