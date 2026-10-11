package synth;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import shared.Random;
import shared.Shared;
import stream.Read;

/** Finite candidate-generation observations without entering the chimera retry loop.
 * @author Brian Bushnell, Jean
 */
public final class MakeChimerasTest{

	/** Emits six observations for comparison with literal external oracles.
	 * @param args One in= argument naming an existing fixture for constructor validation
	 * @throws ReflectiveOperationException If the production candidate method cannot be called
	 */
	public static void main(String[] args) throws ReflectiveOperationException{
		require(args.length==1 && args[0].startsWith("in=") && args[0].length()>3,
				"Specify in=<fixture> for MakeChimeras constructor validation");
		require(MakeChimeras.class.desiredAssertionStatus(), "Production assertions must be enabled");
		final Method candidate=MakeChimeras.class.getDeclaredMethod("makeChimera",
				Read.class, Read.class, Random.class, long.class);
		candidate.setAccessible(true);
		System.out.println("case\tfirst_length\tsecond_length\tforce_length\tactual");
		runCase(candidate, args[0], "one_equal", 1, 1, 1);
		runCase(candidate, args[0], "one_beyond", 1, 1, 2);
		runCase(candidate, args[0], "two_equal", 2, 2, 2);
		runCase(candidate, args[0], "two_beyond", 2, 2, 3);
		runCase(candidate, args[0], "two_feasible", 2, 2, 1);
		runCase(candidate, args[0], "clipped_first", 1, 3, 2);
	}

	/** Calls only the candidate helper; construction does not start a reader or writer. */
	private static void runCase(Method candidate, String input, String label,
			int firstLength, int secondLength, int forced) throws ReflectiveOperationException{
		require(firstLength>0 && secondLength>0 && forced>0,
				"Fixtures model retained nonempty sources with positive forced lengths");
		final MakeChimeras maker=new MakeChimeras(new String[]{input, "out=null", "chimeras=1",
				"forcelength="+forced, "t=1", "simd=f", "bf1=t"});
		final byte[] first=new byte[firstLength], second=new byte[secondLength];
		Arrays.fill(first, (byte)'N');
		Arrays.fill(second, (byte)'N');
		final Read a=new Read(first, null, "first", 0);
		final Read b=new Read(second, null, "second", 1);
		final Read result=(Read)candidate.invoke(maker, a, b, Shared.threadLocalRandom(), 0L);
		final String actual=result==null ? "NULL" : new String(result.bases, StandardCharsets.US_ASCII);
		System.out.println(label+"\t"+firstLength+"\t"+secondLength+"\t"+forced+"\t"+actual);
	}

	/** Keeps fixture preconditions active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Static-only test entry point. */
	private MakeChimerasTest(){}
}
