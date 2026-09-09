package ml;

/**
 * One-line runtime probe for gate receipts: prints whether this JVM's BBTools will
 * take the SIMD feed-forward path ({@code Shared.SIMD} is true only when
 * {@code simd.Vector.simd256} reports AVX2 AND the JVM was started with
 * {@code --add-modules jdk.incubator.vector}; {@code Shared.SIMD_FEED_FORWARD} is the
 * per-run switch) plus the Java version, so a recorded parity number names the code
 * path that produced it.
 *
 * <p>Packaged copy of bbtools-dev {@code magqc_gpu_sweep/SimdProbe.java} (default package,
 * SHA-256 {@code 9b0b160b...} at bbtools-dev {@code c26965b}); identical output.</p>
 *
 * <p>Usage: {@code java [--add-modules jdk.incubator.vector] -cp <BBTools/current> ml.SimdProbe}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class SimdProbe {

	public static void main(String[] args){
		System.out.println("Shared.SIMD="+shared.Shared.SIMD+" SIMD_FEED_FORWARD="+shared.Shared.SIMD_FEED_FORWARD
			+" simd256="+simd.Vector.simd256+" java="+System.getProperty("java.version"));
	}
}
