package ml;

/**
 * Prints the current SIMD master flag, feed-forward flag, detected vector-width
 * capability, and Java version. These values describe configuration, not an
 * observed inference path: array length, per-job flags, and the calling kernel
 * also determine whether vector code runs. The capability flag reports a detected
 * maximum vector width of at least 256 bits, not a particular instruction set.
 *
 * <p>Packaged copy of bbtools-dev {@code magqc_gpu_sweep/SimdProbe.java} (default package,
 * SHA-256 {@code 9b0b160b...} at bbtools-dev {@code c26965b}); identical output.</p>
 *
 * <p>Usage: {@code java [--add-modules jdk.incubator.vector] -cp <BBTools/current> ml.SimdProbe}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class SimdProbe{

	/** Prints one line of configuration values; arguments are ignored. */
	public static void main(String[] args){
		System.out.println("Shared.SIMD="+shared.Shared.SIMD+" SIMD_FEED_FORWARD="+shared.Shared.SIMD_FEED_FORWARD
			+" simd256="+simd.Vector.simd256+" java="+System.getProperty("java.version"));
	}
}
