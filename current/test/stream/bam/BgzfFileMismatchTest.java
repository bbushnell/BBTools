package test.stream.bam;

import java.io.IOException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import stream.bam.TestBgzfA80k;

/** Literal byte/offset cases for the Java8 diagnostic comparison.
 * On newer Java, the original Files.mismatch implementation is a second oracle.
 * No compression stress loop or biological input is run.
 * @author Keqing
 */
public final class BgzfFileMismatchTest{

	/** Runs exact mismatch and prefix cases in one explicitly fresh directory. */
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh absolute fixture directory");}
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root), "Fixture directory must be fresh and absolute");
		Files.createDirectory(root);
		Method original=null;
		try{original=Files.class.getMethod("mismatch", Path.class, Path.class);}
		catch(NoSuchMethodException expected){/* Java8 has no original comparison API. */}
		final Path left=root.resolve("left.bin"), right=root.resolve("right.bin");
		int cases=0;
		compare(left, right, new byte[0], new byte[0], -1, original);cases++;
		compare(left, right, new byte[0], new byte[]{0}, 0, original);cases++;
		compare(left, right, new byte[]{0}, new byte[0], 0, original);cases++;
		compare(left, right, new byte[]{0, (byte)255}, new byte[]{0, (byte)255}, -1, original);cases++;
		compare(left, right, new byte[]{0, (byte)255}, new byte[]{0, 0}, 1, original);cases++;
		final byte[] full=new byte[65537];
		for(int i=0; i<full.length; i++){full[i]=(byte)i;}
		compare(left, right, full, full, -1, original);cases++;
		for(int offset:new int[]{0, 32767, 32768, 65536}){
			final byte[] changed=full.clone();changed[offset]^=(byte)255;
			compare(left, right, full, changed, offset, original);cases++;
		}
		final byte[] prefix=Arrays.copyOf(full, 65536);
		compare(left, right, prefix, full, 65536, original);cases++;
		compare(left, right, full, prefix, 65536, original);cases++;
		boolean rejected=false;
		try{TestBgzfA80k.firstMismatch(root.resolve("absent.bin"), right);}
		catch(IOException expected){rejected=true;}
		check(rejected, "A missing input must fail rather than count as an equal empty file");
		System.out.println("BGZF_FILE_MISMATCH_TEST_PASS cases="+cases+" jdk_oracle="+(original!=null)+" missing_input_rejected=true");
	}

	/** Writes literal inputs and checks both the stated offset and the available JDK oracle. */
	private static void compare(final Path left, final Path right, final byte[] a, final byte[] b,
			final long expected, final Method original)throws Exception{
		Files.write(left, a);
		Files.write(right, b);
		final long actual=TestBgzfA80k.firstMismatch(left, right);
		check(actual==expected, "First differing byte or shorter prefix: expected="+expected+", actual="+actual+
			", left="+a.length+", right="+b.length);
		if(original!=null){
			final long reference=((Long)original.invoke(null, left, right)).longValue();
			check(actual==reference, "Java8 comparison disagrees with original Files.mismatch: "+actual+" vs "+reference);
		}
	}

	/** A failed expectation must terminate the fixture with its diagnostic reason. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
