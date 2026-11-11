package test.assemble;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import assemble.Tadpole2;
import assemble.TargetKmerBridge;
import kmer.AbstractKmerTableSet;
import ukmer.Kmer;

/**
 * Checks that TargetKmerBridge rejects unsupported maxDistance values before sizing its visited table.
 * Run each invocation in a fresh JVM: the test seeds shared kmer state, constructs one
 * Tadpole2 owner and never loads reads, starts workers or emits graph output.
 * @author Brian Bushnell, Jean
 */
public final class TargetKmerBridgeCapacityTest{

	/** Exercises small safe capacities and invalid constructor bounds. */
	public static void main(final String[] args) throws Exception{
		boolean fixed=false;
		for(String arg : args){
			if("fixed=t".equalsIgnoreCase(arg)){fixed=true;}
			else if("fixed=f".equalsIgnoreCase(arg)){fixed=false;}
			else{throw new IllegalArgumentException("Optional mode is fixed=t or fixed=f.");}
		}
		final Path input=Files.createTempFile("target-kmer-bridge-capacity-", ".fa");
		Tadpole2 tadpole=null;
		int failures=0;
		try{
			Files.write(input, Arrays.asList(">tiny", "ACGTACGTACGTACGTACGT"), StandardCharsets.UTF_8);
			System.err.println("TargetKmerBridgeCapacityTest input="+input);
			Kmer.PACKED=false;
			AbstractKmerTableSet.MASK_CORE=false;
			tadpole=new Tadpole2(tadpoleArgs(input), false);
			check(tadpole.k()==63, "No-read Tadpole2 owner must preserve the requested k.");
			failures+=checkCapacity("zero", tadpole, 0, 4);
			failures+=checkCapacity("one", tadpole, 1, 8);
			failures+=checkCapacity("thirty", tadpole, 30, 64);
			failures+=checkRejected("negative", tadpole, -1);
			failures+=checkRejected("int-max", tadpole, Integer.MAX_VALUE);
			if(fixed && failures==0){failures+=checkRejected("overflow-next", tadpole, 536870911);}
		}finally{
			if(tadpole!=null){tadpole.clearData();}
			Files.deleteIfExists(input);
		}
		if(failures>0){throw new AssertionError("TargetKmerBridge capacity failures="+failures);}
		System.out.println("TARGET_KMER_BRIDGE_CAPACITY_PASS checks="+checks);
	}

	/** Constructor arguments matching the validated no-read Tadpole2 setup. */
	private static String[] tadpoleArgs(final Path input){
		return new String[]{
				"in="+input, "out=null", "k=63", "mode=extend", "t=1",
				"prefilter=0", "prealloc=f", "ways=1", "initialsize=128",
				"ecc=f", "ecco=f", "merge=f", "markerrors=f", "hpindel=f",
				"fixindels=f", "hollow=f", "hollowbranchstop=f",
				"hollowbranchrollback=f", "hollowrescue=f", "el=0", "er=1",
				"showstats=f"};
	}

	/** Verifies a supported bound's rounded power-of-two visited table length. */
	private static int checkCapacity(final String name, final Tadpole2 tadpole,
			final int maxDistance, final int expected) throws Exception{
		final TargetKmerBridge bridge=new TargetKmerBridge(tadpole, maxDistance, false);
		final int actual=visitedCapacity(bridge);
		final boolean ok=actual==expected;
		System.out.println("case="+name+" maxDistance="+maxDistance+
				" expectedRejected=false rejected=false expectedCapacity="+expected+
				" actualCapacity="+actual+" capacityOk="+ok);
		checks++;
		return ok ? 0 : 1;
	}

	/** Verifies that an unsupported bound fails before allocating a visited table. */
	private static int checkRejected(final String name, final Tadpole2 tadpole,
			final int maxDistance) throws Exception{
		int actual=-1;
		boolean rejected=false;
		try{
			actual=visitedCapacity(new TargetKmerBridge(tadpole, maxDistance, false));
		}catch(IllegalArgumentException expected){
			rejected=true;
		}
		System.out.println("case="+name+" maxDistance="+maxDistance+
				" expectedRejected=true rejected="+rejected+" actualCapacity="+actual+
				" rejectionOk="+rejected);
		checks++;
		return rejected ? 0 : 1;
	}

	/** Reads the private visited hash table length without changing production visibility. */
	private static int visitedCapacity(final TargetKmerBridge bridge) throws Exception{
		final Field field=TargetKmerBridge.class.getDeclaredField("visitedHashes");
		field.setAccessible(true);
		return ((long[])field.get(bridge)).length;
	}

	/** Unit failures stay fatal even if assertions are accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Number of constructor cases checked. */
	private static int checks;
}
