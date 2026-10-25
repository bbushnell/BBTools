package test.cardinality;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import cardinality.CardinalityParser;

/** Checks that only the DDL worker owning index zero truncates and writes out4.
 * @author Sucrose
 */
public final class DDLCalibrationDriver2Out4OwnershipTest{

	/** Runs the captured baseline or the corrected ownership contract. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation: "+args[0]);
		final Path root=Files.createTempDirectory("bbtools-ddlcal2-out4-");
		final Path out4=root.resolve("out4.tsv");
		Files.write(out4, SENTINEL.getBytes(StandardCharsets.UTF_8));
		final Thread nonOwner=newCalThread(out4);
		nonOwner.run();
		final boolean success=success(nonOwner);
		final String observed=new String(Files.readAllBytes(out4), StandardCharsets.UTF_8);
		if(fixed){
			check(success, "Non-owner worker hit a swallowed exception");
			check(SENTINEL.equals(observed),
				"Non-owner worker should not truncate out4; observed length "+observed.length());
			checkOwnerRowSurvives(root.resolve("owner.tsv"));
		}else{
			check(success, "Baseline non-owner worker hit a swallowed exception");
			check(observed.length()==0,
				"Baseline unexpectedly preserved out4 ownership sentinel");
		}
		System.out.println("PASS\tfixed="+fixed+"\tsuccess="+success+
			"\tout4_length="+observed.length());
	}

	/** Instantiates a worker that has an out4 path but cannot claim any DDL index. */
	private static Thread newCalThread(final Path out4) throws Exception{
		return newCalThread(out4, new AtomicLong(1));
	}

	/** Instantiates a worker that shares the caller's next-DDL index. */
	private static Thread newCalThread(final Path out4, final AtomicLong nextDDL) throws Exception{
		final Class<?> clazz=Class.forName("cardinality.DDLCalibrationDriver2$CalThread");
		final Constructor<?> cons=clazz.getDeclaredConstructor(long.class, int.class,
			AtomicLong.class, long[].class, int.class, int.class, long.class, String.class,
			String.class, int.class, boolean.class);
		cons.setAccessible(true);
		return (Thread)cons.newInstance(1L, 1, nextDDL, new long[]{1},
			64, 1, 1L, "dll4", out4.toString(), 1, false);
	}

	/** Checks a real out4 row from worker zero is not lost when a later worker exits. */
	private static void checkOwnerRowSurvives(final Path out4) throws Exception{
		CardinalityParser.initializeAll("dll4", 64, 1, null, null);
		final AtomicLong nextDDL=new AtomicLong(0);
		final Thread owner=newCalThread(out4, nextDDL);
		owner.run();
		check(success(owner), "Owner worker hit a swallowed exception");
		final byte[] written=Files.readAllBytes(out4);
		final String observed=new String(written, StandardCharsets.UTF_8);
		final String[] lines=observed.split("\n", -1);
		check(lines.length==2 && lines[1].isEmpty(),
			"Expected one out4 owner row; observed "+observed.replace('\n', '|'));
		final String[] fields=lines[0].split("\t", -1);
		check(fields.length==42, "Expected 42 out4 fields; observed "+fields.length);
		check("1".equals(fields[0]), "Expected true-cardinality 1; observed "+fields[0]);
		final Thread nonOwner=newCalThread(out4, nextDDL);
		nonOwner.run();
		check(success(nonOwner), "Late non-owner worker hit a swallowed exception");
		check(Arrays.equals(written, Files.readAllBytes(out4)),
			"Late non-owner worker should not truncate the owner out4 row");
	}

	/** Returns the worker completion flag that would be false after a swallowed exception. */
	private static boolean success(final Thread thread) throws Exception{
		final Field f=thread.getClass().getDeclaredField("success");
		f.setAccessible(true);
		return f.getBoolean(thread);
	}

	/** Keeps oracle checks active under both assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String SENTINEL="owner sentinel\n";
}
