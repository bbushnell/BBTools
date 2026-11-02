package test.scalar;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import shared.Timer;

/** Finite constructor-level checks for scalar PartitionReads variants. */
public class PartitionReadsContractTest{

	/** Runs one independent PartitionReads argument-parsing or setup case. */
	public static void main(String[] args) throws Exception{
		final String variant=args[0], mode=args[1];
		final Path dir=Paths.get(args[2]);
		Files.createDirectories(dir);
		final Path input=dir.resolve("reads.fa");
		final Path equalsInput=dir.resolve("reads=1.fa");
		final Path hashInput=dir.resolve("reads#tag.fa");
		write(input, ">read1\nACGTACGTACGT\n");
		write(equalsInput, ">read2\nACGTACGTACGT\n");
		write(hashInput, ">read3\nACGTACGTACGT\n");

		if(mode.equals("ordinary-output")){
			final String pattern=dir.resolve("part%.fq").toString();
			final Object partition=construct(variant, input, pattern);
			require(pattern.equals(field(partition, "out1")), "ordinary output pattern changed");
		}else if(mode.equals("equals-output")){
			final String pattern=dir.resolve("part=%#.fq").toString();
			final Object partition=construct(variant, new String[]{"in="+input,
					"out="+pattern, "ways=2", "ow=t"});
			require(pattern.replace("#", "1").equals(field(partition, "out1")),
					"out1 lost text after the out= key");
			require(pattern.replace("#", "2").equals(field(partition, "out2")),
					"out2 lost text after the out= key");
		}else if(mode.equals("null-output")){
			final Object partition=construct(variant, new String[]{"in="+input,
					"out=", "ways=2", "ow=t"});
			require(field(partition, "out1")==null, "out= should be normalized to null");
		}else if(mode.equals("equals-input")){
			construct(variant, equalsInput, dir.resolve("equals%.fq").toString());
		}else if(mode.equals("hash-input")){
			construct(variant, new String[]{"in="+hashInput, dir.resolve("hash%.fq").toString(),
					"ways=2", "ow=t"});
		}else if(mode.equals("valid-coverage")){
			final Path cov=dir.resolve("coverage.tsv");
			write(cov, "#ID\tAvg_fold\ncontig\t2.5\n");
			construct(variant, new String[]{"in="+input, dir.resolve("cov%.fq").toString(),
					"mode=depth", "cov="+cov, "ways=2", "ow=t"});
		}else if(mode.equals("custom-cutoff")){
			final Path pattern=dir.resolve("cutoff%.fq");
			final Object partition=construct(variant, new String[]{"in="+input,
					pattern.toString(), "mode=gc", "cutoff=0.5", "ow=t"});
			process(partition);
			final String low=read(dir.resolve("cutoff0.fq"));
			final String high=read(dir.resolve("cutoff1.fq"));
			require(!low.contains("read1"), "read at cutoff boundary was written to partition0");
			require(high.contains("read1"), "read at cutoff boundary was not written to partition1");
		}else if(mode.equals("custom-cutoff-count")){
			final Object partition=construct(variant, new String[]{"in="+input,
					dir.resolve("count%.fq").toString(), "mode=gc", "cutoff=0.5", "ow=t"});
			process(partition);
			requireCount(partition, 1, 12, "custom cutoff");
		}else if(mode.equals("metric-count")){
			final Object partition=construct(variant, new String[]{"in="+input,
					dir.resolve("metric%.fq").toString(), "mode=gc", "ways=2", "ow=t"});
			process(partition);
			requireCount(partition, 1, 12, "normal metric");
		}else if(mode.equals("empty-coverage")){
			final Path cov=dir.resolve("empty.cov");
			write(cov, "");
			try{
				construct(variant, new String[]{"in="+input, dir.resolve("empty%.fq").toString(),
						"mode=depth", "depth="+cov, "ways=2", "ow=t"});
			}catch(IllegalArgumentException expected){
				System.out.println("PASS "+variant+" "+mode);
				return;
			}
			throw new AssertionError("Empty coverage must be rejected with a clear argument error");
		}else{
			throw new IllegalArgumentException(mode);
		}
		System.out.println("PASS "+variant+" "+mode);
	}

	/** Constructs the requested PartitionReads implementation. */
	private static Object construct(String variant, Path input, String output) throws Exception{
		return construct(variant, new String[]{"in="+input, output, "ways=2", "ow=t"});
	}

	/** Reflectively constructs the selected historical or current partitioner. */
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

	/** Maps fixture variant numbers to production classes. */
	private static String className(String variant){
		if(variant.equals("1")){return "scalar.PartitionReads";}
		if(variant.equals("2")){return "scalar.PartitionReads2";}
		if(variant.equals("3")){return "scalar.PartitionReads3";}
		throw new IllegalArgumentException("Unknown PartitionReads variant: "+variant);
	}

	/** Reads a private string field without changing production visibility. */
	private static Object field(Object target, String name) throws Exception{
		final Field f=target.getClass().getDeclaredField(name);
		f.setAccessible(true);
		return f.get(target);
	}

	/** Invokes the package-private process method without changing production visibility. */
	private static void process(Object target) throws Exception{
		final Method m=target.getClass().getDeclaredMethod("process", Timer.class);
		m.setAccessible(true);
		m.invoke(target, new Timer());
	}

	/** Checks one read and its 12 bases were counted exactly once. */
	private static void requireCount(Object partition, long reads, long bases, String label)
			throws Exception{
		require(((Long)field(partition, "readsProcessed")).longValue()==reads,
				label+" run counted reads incorrectly");
		require(((Long)field(partition, "basesProcessed")).longValue()==bases,
				label+" run counted bases incorrectly");
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
