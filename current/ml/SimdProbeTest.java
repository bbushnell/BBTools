package ml;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

/**
 * Fixture for {@link SimdProbe}: one line, the four expected keys, values that are
 * booleans / a version string, and {@code Shared.SIMD} consistent with
 * {@code simd.Vector.simd256} (SIMD can only be on when AVX2 is reported).
 *
 * <p>Run with assertions on: {@code java -ea [--add-modules jdk.incubator.vector] ml.SimdProbeTest}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class SimdProbeTest {

	public static void main(String[] args) throws Exception{
		final PrintStream old=System.out;
		final ByteArrayOutputStream bo=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(bo, true, "UTF-8");
		System.setOut(ps);
		try{SimdProbe.main(new String[0]);}
		finally{System.setOut(old); ps.flush();}
		final String out=new String(bo.toByteArray(), StandardCharsets.UTF_8);
		final String[] lines=out.split("\n", -1);
		check(lines.length==2 && lines[1].length()==0, "expected exactly one line: "+out);
		final String[] kv=lines[0].split(" ");
		check(kv.length==4, "expected 4 key=value fields: "+lines[0]);
		check(kv[0].startsWith("Shared.SIMD=") && kv[1].startsWith("SIMD_FEED_FORWARD=") && kv[2].startsWith("simd256=")
			&& kv[3].startsWith("java="), "unexpected keys: "+lines[0]);
		final boolean simd=Boolean.parseBoolean(kv[0].substring(12)), simd256=Boolean.parseBoolean(kv[2].substring(8));
		check(kv[0].endsWith("=true") || kv[0].endsWith("=false"), "Shared.SIMD not boolean: "+kv[0]);
		check(kv[1].endsWith("=true") || kv[1].endsWith("=false"), "SIMD_FEED_FORWARD not boolean: "+kv[1]);
		check(kv[2].endsWith("=true") || kv[2].endsWith("=false"), "simd256 not boolean: "+kv[2]);
		check(!simd || simd256, "Shared.SIMD=true while simd256=false: "+lines[0]);
		check(kv[3].length()>5, "empty java version: "+kv[3]);
		System.err.println("SimdProbeTest: PASS ("+lines[0]+")");
	}

	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
