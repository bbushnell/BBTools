package test.scalar;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import shared.Timer;

/** Auto-partition SAM header preservation check for PartitionReads2/3. */
public class PartitionReadsAutoSamHeaderTest{

	private static final String SAM_TEXT="@HD\tVN:1.6\tSO:unsorted\n"
			+"@SQ\tSN:ref\tLN:100\n"
			+"low0\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low1\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low2\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low3\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"high0\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high1\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high2\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high3\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n";

	/** Exercises one SAM output creation path for one auto-capable variant. */
	public static void main(String[] args) throws Exception{
		final String variant=args[0], mode=args[1];
		if(variant.equals("1")){
			throw new IllegalArgumentException("PartitionReads v1 has no auto-partition SAM path");
		}

		final Path dir=Paths.get(args[2]);
		Files.createDirectories(dir);
		final Path input=dir.resolve("reads.sam");
		write(input, SAM_TEXT);

		if(mode.equals("explicit")){
			runExplicit(variant, input, dir);
		}else if(mode.equals("auto")){
			runAuto(variant, input, dir);
		}else{
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		System.out.println("PASS "+variant+" "+mode+"-sam-header");
	}

	/** Runs the construction-time FileFormat path for the metric SAM control. */
	private static void runExplicit(String variant, Path input, Path dir) throws Exception{
		final Object explicit=construct(variant, new String[]{"in="+input,
				dir.resolve("explicit%.sam").toString(), "mode=gc", "ways=2", "ow=t"});
		process(explicit);
		requireSamHeader(dir.resolve("explicit0.sam"));
		requireSamHeader(dir.resolve("explicit1.sam"));
	}

	/** Runs the auto-partition path that defers FileFormat creation until ways is known. */
	private static void runAuto(String variant, Path input, Path dir) throws Exception{
		final Object auto=construct(variant, new String[]{"in="+input,
				dir.resolve("auto%.sam").toString(), "mode=gc", "auto=t", "maxp=2", "ow=t"});
		process(auto);
		requireSamHeader(dir.resolve("auto0.sam"));
		requireSamHeader(dir.resolve("auto1.sam"));
	}

	/** Reflectively constructs the selected auto-capable partitioner. */
	private static Object construct(String variant, String[] args) throws Exception{
		final Class<?> clazz=Class.forName(className(variant));
		final Constructor<?> c=clazz.getConstructor(String[].class);
		try{
			return c.newInstance((Object)args);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Maps supported fixture variant numbers to auto-capable production classes. */
	private static String className(String variant){
		if(variant.equals("2")){return "scalar.PartitionReads2";}
		if(variant.equals("3")){return "scalar.PartitionReads3";}
		throw new IllegalArgumentException("Unknown auto-capable PartitionReads variant: "+variant);
	}

	/** Invokes package-private processing without changing production visibility. */
	private static void process(Object target) throws Exception{
		final Method m=target.getClass().getDeclaredMethod("process", Timer.class);
		m.setAccessible(true);
		m.invoke(target, new Timer());
	}

	/** Requires the SAM output to preserve the input @SQ header and mapped ref body. */
	private static void requireSamHeader(Path path) throws Exception{
		final String text=read(path);
		require(text.contains("@SQ\tSN:ref\tLN:100\n"), path+" lost @SQ ref header");
		require(text.contains("	ref	"), path+" lost mapped ref alignments");
	}

	/** Reads an expected output shard as UTF-8 text. */
	private static String read(Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	/** Writes deterministic UTF-8 fixture input. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}

	/** These checks remain active under both -ea and -da. */
	private static void require(boolean valid, String message){
		if(!valid){throw new AssertionError(message);}
	}
}
