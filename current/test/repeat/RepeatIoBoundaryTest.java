package test.repeat;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import fileIO.FileFormat;
import repeat.Repeat;
import repeat.RepeatFinder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Finite constructor, short-region statistics and owned-writer shutdown checks.
 * No input workers start. Writer cases invoke the actual TSV writer and wait for
 * its shutdown. The Linux /dev/full case isolates a buffered flush/close failure;
 * it does not inject failures into unrelated processes or test borrowed stdout.
 * @author Sucrose
 */
public final class RepeatIoBoundaryTest{

	/** Requires work=fresh-directory with an existing parent; preserves all small fixture files. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("work=")){
			throw new IllegalArgumentException("Expected work=fresh-directory");
		}
		final Path work=Paths.get(args[0].substring(5));
		Files.createDirectory(work);
		final Path ordinary=work.resolve("input.fa"), equals=work.resolve("input=literal.fa");
		final byte[] fasta=">literal\nAAAAAA\n".getBytes(StandardCharsets.US_ASCII);
		Files.write(ordinary, fasta);
		Files.write(equals, fasta);
		final RepeatIoBoundaryTest test=new RepeatIoBoundaryTest();
		test.check("ordinary filename control", () -> constructor(ordinary));
		test.check("equals filename", () -> constructor(equals));
		test.check("short region below entropy k", () -> shortRegion(ordinary, "AAA"));
		test.check("region equal to entropy k", () -> shortRegion(ordinary, "AAAAA"));
		test.check("longer homogeneous region control", () -> shortRegion(ordinary, "AAAAAA"));
		test.check("disabled preview columns", () -> preview(ordinary, false));
		test.check("enabled preview columns control", () -> preview(ordinary, true));
		test.check("regular writer control", () -> writer(ordinary, work.resolve("report.tsv"), false));
		test.check("owned writer flush failure", () -> writer(ordinary, Paths.get("/dev/full"), true));
		System.out.println("RepeatIoBoundaryTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Repeat I/O boundary failures="+test.failures);}
	}

	/** Supplies a real legal filename; the constructor must retain it through path validation. */
	private static RepeatFinder constructor(Path input){
		final RepeatFinder finder=new RepeatFinder(new String[] {"in="+input, "out=null", "depth=f",
				"tandem=t", "entropy=f", "k=3", "minrepeat=0", "mask=f"});
		require(!finder.errorState, "Fresh detector unexpectedly starts in an error state");
		return finder;
	}

	/** Retains a monomer tandem and checks finite zero entropy for homogeneous sequence.
	 * The shorter input has no complete default 5mer; it must not require an invalid
	 * tracker window merely to retain/report an otherwise valid three-base region. */
	private static void shortRegion(Path input, String sequence) throws Exception{
		final RepeatFinder finder=constructor(input);
		final Class<?> type=Class.forName("repeat.RepeatFinder$ProcessThread");
		final Constructor<?> ctor=type.getDeclaredConstructor(RepeatFinder.class,
				ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, int.class);
		ctor.setAccessible(true);
		final Object worker=ctor.newInstance(finder, null, null, 0);
		final Method method=type.getDeclaredMethod("processReadShortTandem", Read.class,
				int.class, int.class, int.class, int.class);
		method.setAccessible(true);
		final Read read=new Read(sequence.getBytes(StandardCharsets.US_ASCII), null, "literal", 0, false);
		final int count=((Integer)invoke(method, worker, read, 1, 1, 2, 2)).intValue();
		final Field setField=type.getDeclaredField("repeatSet");
		setField.setAccessible(true);
		final Object set=setField.get(worker);
		final Field closedField=set.getClass().getDeclaredField("closedRepeats");
		closedField.setAccessible(true);
		final List<?> closed=(List<?>)closedField.get(set);
		require(count==sequence.length() && closed.size()==1, "Expected one full homogeneous region");
		final Repeat repeat=(Repeat)closed.get(0);
		require(repeat.start==0 && repeat.stop==sequence.length()-1 && repeat.entropy==0f,
				"Expected full inclusive region and finite zero entropy, actual="+repeat.entropy);
	}

	/** Exercises the parsed preview option and literal column expectations, restoring its global setting. */
	private static void preview(Path input, boolean enabled){
		final int prior=Repeat.SEQ_AFFIX_LEN;
		try{
			new RepeatFinder(new String[] {"in="+input, "out=null", "depth=f", "k=3",
					"preview="+(enabled ? 27 : 0)});
			final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null, "literal", 0, false);
			final Repeat repeat=new Repeat(read, 0, 2, 3, 0, 0, 'R');
			repeat.stop=3;
			final String[] header=Repeat.tsvHeader().split("\t", -1), row=repeat.toString().split("\t", -1);
			final int columns=(enabled ? 14 : 13);
			require(header.length==columns && row.length==columns,
					"Expected "+columns+" columns; header="+header.length+" row="+row.length);
			require(header[columns-1].equals(enabled ? "seq" : "cName")
					&& row[columns-1].equals(enabled ? "ACGT" : "literal"), "Unexpected final report column");
		}finally{Repeat.SEQ_AFFIX_LEN=prior;}
	}

	/** A tiny header remains buffered until finalization; checks the caller-visible error flag. */
	private static void writer(Path input, Path output, boolean fail) throws Exception{
		if(fail){require(Files.exists(output) && Files.isWritable(output), "Linux writable /dev/full is required");}
		final RepeatFinder finder=constructor(input);
		final FileFormat ff=FileFormat.testOutput(output.toString(), FileFormat.TXT, null, false, true, false, false);
		final Method method=RepeatFinder.class.getDeclaredMethod("printRepeats", FileFormat.class, ArrayList.class);
		method.setAccessible(true);
		invoke(method, finder, ff, new ArrayList<Repeat>());
		require(finder.errorState==fail, "Writer failure flag must reach parent; expected="+fail+" actual="+finder.errorState);
		if(!fail){
			final String actual=new String(Files.readAllBytes(output), StandardCharsets.US_ASCII);
			require(actual.equals(Repeat.tsvHeader()+"\n"), "Empty report must contain exactly its header and newline");
		}
	}

	/** Reflects actual package-private implementations, preserving their original failure types. */
	private static Object invoke(Method method, Object target, Object... args) throws Exception{
		try{return method.invoke(target, args);}
		catch(InvocationTargetException failure){
			final Throwable cause=failure.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new AssertionError("Unexpected reflected failure", cause);
		}
	}

	/** Runs independent groups even if a preceding baseline case fails. */
	private void check(String name, Action action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(Exception | AssertionError failure){failures++;System.out.println(name+" FAIL: "+failure);}
	}

	/** Keeps literal acceptance checks enabled under both assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One independent case, possibly doing bounded file I/O. */
	private interface Action{void run() throws Exception;}
	/** Attempted and failed case counts. */
	private int cases, failures;
}
