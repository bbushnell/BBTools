package test.icecream;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import icecream.IceCreamFinder;
import icecream.ReformatPacBio;
import icecream.ZMWStreamer;
import shared.Shared;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Tests the two PacBio trimming implementations against literal end/quality/header oracles.
 * @author Brian Bushnell, Sucrose
 */
public final class PacBioTrimTest{

	/** Constructs workers without starting them, then invokes their trimming helpers.
	 * @param args One in= argument naming an existing tiny input for constructor validation
	 * @throws ReflectiveOperationException If the production worker/helper signatures change
	 */
	public static void main(String[] args) throws ReflectiveOperationException{
		require(args.length==1 && args[0].startsWith("in=") && args[0].length()>3,
				"Specify in=<existing fixture> for constructor validation; no input processing occurs");
		require(IceCreamFinder.class.desiredAssertionStatus() && ReformatPacBio.class.desiredAssertionStatus(),
				"Production assertions must be enabled");
		final boolean originalJNI=Shared.USE_JNI;
		Shared.USE_JNI=false;
		int failures=0;
		try{
			final String[] options={args[0], "out=null", "t=1", "simd=f", "pigz=f", "unpigz=f"};
			final Worker finder=new Worker(new IceCreamFinder(options));
			final Worker reformat=new Worker(new ReformatPacBio(options));
			System.out.println("owner\tcheck\tcase\texpected\tobserved\tresult");
			failures+=runCases(finder);
			failures+=runCases(reformat);
		}finally{
			Shared.USE_JNI=originalJNI;
		}
		require(failures==0, failures+" PacBio trimming observations differ from the literal oracle; see all case rows");
	}

	/** Gives each side and the actual trimming method independent finite expectations. */
	private static int runCases(Worker worker) throws ReflectiveOperationException{
		int failures=0;
		failures+=helperCase(worker, "no_N", "ACGTAC", 0, 0);
		failures+=helperCase(worker, "one_right_N", "ACGTAN", 0, 1);
		failures+=helperCase(worker, "two_right_N", "ACGTNN", 0, 2);
		failures+=helperCase(worker, "one_left_N", "NACGTA", 1, 0);
		failures+=helperCase(worker, "two_left_N", "NNACGT", 2, 0);
		failures+=helperCase(worker, "both_ends", "NACGTN", 1, 1);
		failures+=helperCase(worker, "empty", "", 0, 0);
		failures+=trimCase(worker, "no_N", "ACGTAC", 0, "ACGTAC", "movie/11/100_106",
				new byte[]{10, 11, 12, 13, 14, 15});
		failures+=trimCase(worker, "one_right_N", "ACGTAN", 1, "ACGTA", "movie/11/100_105",
				new byte[]{10, 11, 12, 13, 14});
		failures+=trimCase(worker, "two_right_N", "ACGTNN", 2, "ACGT", "movie/11/100_104",
				new byte[]{10, 11, 12, 13});
		failures+=trimCase(worker, "one_left_N", "NACGTA", 1, "ACGTA", "movie/11/101_106",
				new byte[]{11, 12, 13, 14, 15});
		failures+=trimCase(worker, "two_left_N", "NNACGT", 2, "ACGT", "movie/11/102_106",
				new byte[]{12, 13, 14, 15});
		failures+=trimCase(worker, "both_ends", "NACGTN", 2, "ACGT", "movie/11/101_105",
				new byte[]{11, 12, 13, 14});
		return failures;
	}

	/** Empty input is checked only at the byte-array helpers, without requiring an empty Read. */
	private static int helperCase(Worker worker, String label, String bases, int expectedLeft,
			int expectedRight) throws ReflectiveOperationException{
		final byte[] bytes=bases.getBytes(StandardCharsets.US_ASCII);
		final int left=(Integer)worker.left.invoke(worker.instance, bytes, LOOKAHEAD);
		final int right=(Integer)worker.right.invoke(worker.instance, bytes, LOOKAHEAD);
		return report(worker.name, "left_count", label, ""+expectedLeft, ""+left)
				+report(worker.name, "right_count", label, ""+expectedRight, ""+right);
	}

	/** Invokes the real read mutator and checks the returned count plus all three mutable fields. */
	private static int trimCase(Worker worker, String label, String bases, int expectedTrimmed,
			String expectedBases, String expectedHeader, byte[] expectedQuality) throws ReflectiveOperationException{
		final byte[] bytes=bases.getBytes(StandardCharsets.US_ASCII);
		final byte[] quality=new byte[bytes.length];
		for(int i=0; i<quality.length; i++){quality[i]=(byte)(10+i);}
		final Read read=new Read(bytes, quality, "movie/11/100_106", 0);
		final int trimmed=(Integer)worker.trim.invoke(worker.instance, read, LOOKAHEAD);
		final String expected=expectedTrimmed+"|"+expectedBases+"|"+expectedHeader+"|"+Arrays.toString(expectedQuality);
		final String observed=trimmed+"|"+new String(read.bases, StandardCharsets.US_ASCII)+"|"+read.id+"|"+Arrays.toString(read.quality);
		return report(worker.name, "trim_read", label, expected, observed);
	}

	/** Reports every finite observation before main's final assertion. */
	private static int report(String owner, String check, String label, String expected, String observed){
		final boolean pass=expected.equals(observed);
		System.out.println(owner+"\t"+check+"\t"+label+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Retains private production visibility while giving both implementations the same oracle. */
	private static final class Worker{
		Worker(Object owner) throws ReflectiveOperationException{
			final Class<?> ownerClass=owner.getClass();
			name=ownerClass.getSimpleName();
			final Class<?> workerClass=Class.forName(ownerClass.getName()+"$ProcessThread");
			final Constructor<?> constructor;
			if(owner instanceof IceCreamFinder){
				constructor=workerClass.getDeclaredConstructor(IceCreamFinder.class, ZMWStreamer.class,
						ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
						ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, int.class);
				constructor.setAccessible(true);
				instance=constructor.newInstance(owner, null, null, null, null, null, 0);
			}else{
				require(owner instanceof ReformatPacBio, "Only the two assigned PacBio trimmers belong in this fixture");
				constructor=workerClass.getDeclaredConstructor(ReformatPacBio.class, ZMWStreamer.class,
						ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, int.class);
				constructor.setAccessible(true);
				instance=constructor.newInstance(owner, null, null, null, 0);
			}
			left=workerClass.getDeclaredMethod("calcLeftTrim", byte[].class, int.class);
			right=workerClass.getDeclaredMethod("calcRightTrim", byte[].class, int.class);
			trim=workerClass.getDeclaredMethod("trimRead", Read.class, int.class);
			left.setAccessible(true);
			right.setAccessible(true);
			trim.setAccessible(true);
		}

		final String name;
		final Object instance;
		final Method left, right, trim;
	}

	/** Keeps fixture preconditions and oracle failures active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Shorter than the defined middle, so left and right end scans remain independent. */
	private static final int LOOKAHEAD=3;

	/** Static-only fixture. */
	private PacBioTrimTest(){}
}
