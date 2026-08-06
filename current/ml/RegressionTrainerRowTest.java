package ml;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Stream;

import shared.Timer;

/** Focused end-to-end checks for row-oriented data and external validation. */
public class RegressionTrainerRowTest {

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("regressiontrainer-row-test-");
		try{
			final Path train=dir.resolve("train.tsv");
			final Path valid=dir.resolve("valid.tsv");
			Files.write(train, TRAIN.getBytes(StandardCharsets.US_ASCII));
			Files.write(valid, VALID.getBytes(StandardCharsets.US_ASCII));

			testExternalValidation(dir, train, valid);
			testLegacyFraction(dir, train);
			testSimdExternalValidation(dir, train, valid);
			testConflictingValidationArguments(dir, train, valid);
			System.err.println("RegressionTrainerRowTest: PASS");
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	private static void testExternalValidation(Path dir, Path train, Path valid) throws Exception{
		final Path out=dir.resolve("external.bbnet");
		final RegressionTrainer trainer=new RegressionTrainer(new String[]{
			"in="+train, "valin="+valid, "out="+out, "dims=3,4,2",
			"epochs=2", "batch=2", "seed=3", "final=sigmoid"
		});
		trainer.process(new Timer());

		assertEquals(6, getInt(trainer, "numTrain"), "external training rows");
		assertEquals(2, getInt(trainer, "numValid"), "external validation rows");
		final double[] mean=(double[])getField(trainer, "mean");
		assertClose(2.5, mean[0], "training-only mean[0]");
		assertClose(1.0, mean[1], "training-only mean[1]");
		assertClose(0.5, mean[2], "training-only mean[2]");
		if(CellNetParser.load(out.toString())==null){
			throw new AssertionError("external-validation output net did not load");
		}

		final Path continued=dir.resolve("continued.bbnet");
		final RegressionTrainer netinTrainer=new RegressionTrainer(new String[]{
			"in="+train, "valin="+valid, "out="+continued, "netin="+out,
			"epochs=1", "batch=2", "seed=3", "final=sigmoid"
		});
		netinTrainer.process(new Timer());
		if(CellNetParser.load(continued.toString())==null){
			throw new AssertionError("netin output net did not load");
		}
	}

	private static void testLegacyFraction(Path dir, Path train) throws Exception{
		final Path out=dir.resolve("legacy.bbnet");
		final RegressionTrainer trainer=new RegressionTrainer(new String[]{
			"in="+train, "out="+out, "dims=3,4,2", "vfraction=0.5",
			"epochs=1", "batch=2", "seed=3", "final=sigmoid"
		});
		trainer.process(new Timer());
		assertEquals(3, getInt(trainer, "numTrain"), "legacy training rows");
		assertEquals(3, getInt(trainer, "numValid"), "legacy validation rows");
	}

	private static void testSimdExternalValidation(Path dir, Path train, Path valid) throws Exception{
		final Path out=dir.resolve("simd.bbnet");
		final RegressionTrainer trainer=new RegressionTrainer(new String[]{
			"in="+train, "valin="+valid, "out="+out, "dims=3,8,2",
			"epochs=1", "batch=2", "seed=3", "final=sigmoid", "simd=t"
		});
		trainer.process(new Timer());
		if(CellNetParser.load(out.toString())==null){
			throw new AssertionError("SIMD external-validation output net did not load");
		}
	}

	private static void testConflictingValidationArguments(Path dir, Path train, Path valid){
		try{
			new RegressionTrainer(new String[]{
				"in="+train, "valin="+valid, "out="+dir.resolve("conflict.bbnet"),
				"dims=3,4,2", "vfraction=0.2"
			});
			throw new AssertionError("valin plus vfraction>0 should be rejected");
		}catch(IllegalArgumentException expected){
			if(!expected.getMessage().contains("either valin=")){throw expected;}
		}
	}

	private static Object getField(Object owner, String name) throws Exception{
		final Field field=RegressionTrainer.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(owner);
	}

	private static int getInt(Object owner, String name) throws Exception{
		return ((Integer)getField(owner, name)).intValue();
	}

	private static void assertEquals(int expected, int actual, String label){
		if(expected!=actual){
			throw new AssertionError(label+": expected "+expected+", observed "+actual);
		}
	}

	private static void assertClose(double expected, double actual, String label){
		if(Math.abs(expected-actual)>1e-12){
			throw new AssertionError(label+": expected "+expected+", observed "+actual);
		}
	}

	private static final String TRAIN=
		"#dims\t3\t2\n"
		+"0\t0\t0\t0\t0\n"
		+"1\t0\t1\t0.2\t0.1\n"
		+"2\t1\t0\t0.4\t0.2\n"
		+"3\t1\t1\t0.6\t0.3\n"
		+"4\t2\t0\t0.8\t0.4\n"
		+"5\t2\t1\t1\t0.5\n";

	private static final String VALID=
		"#dims\t3\t2\n"
		+"100\t100\t100\t0\t1\n"
		+"200\t200\t200\t1\t0\n";
}
