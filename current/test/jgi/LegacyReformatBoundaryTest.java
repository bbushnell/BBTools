package test.jgi;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;
import parse.Parse;
import stream.Read;

/** Direct legacy Reformat edit-filter controls.
 * @author Jean
 */
public final class LegacyReformatBoundaryTest{

	/** Requires phase and a fresh output directory for constructor-owned inputs. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline|fixed");
		}
		if(directory==null){throw new IllegalArgumentException("Set outdir to a new directory.");}
		final boolean baseline=phase.equals("baseline");
		final Path root=Files.createDirectory(Paths.get(directory)).toAbsolutePath();
		final Path input=root.resolve("input.fq");
		Files.write(input, "@r0\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII));
		minimumSubs(ReformatReads.class, input, baseline);
		minimumSubs(ReformatReadsOld.class, input, baseline);
		System.out.println("LEGACY_REFORMAT_BOUNDARY_PASS phase="+phase+" checks="+checks);
	}

	/** A stored-perfect legacy read has zero substitutions and must obey a positive minimum too. */
	private static void minimumSubs(final Class<?> type, final Path input,
			final boolean baseline) throws Exception{
		final Object min1=driver(type, input, 1);
		final Object min0=driver(type, input, 0);
		check(passes(min1, read("perfect", "mmmm", true))==baseline,
			type.getSimpleName()+": baseline accepts a perfect zero-substitution read, the repair must reject minsubs=1.");
		check(passes(min0, read("perfect_min0", "mmmm", true)),
			type.getSimpleName()+": minsubs=0 must still accept a perfect alignment.");
		check(!passes(min1, read("zero_subs", "mmmm", false)),
			type.getSimpleName()+": minsubs=1 must reject a non-perfect zero-substitution alignment.");
		check(passes(min1, read("one_sub", "mSmm", false)),
			type.getSimpleName()+": one literal substitution meets minsubs=1.");
	}

	/** Uses the real public constructor with a tiny valid input and disabled output. */
	private static Object driver(final Class<?> type, final Path input, final int minsubs) throws Exception{
		final String[] args=new String[]{"in="+input, "out=null", "minsubs="+minsubs,
			"qin=33", "qout=33", "int=f", "minlength=0", "simd=f"};
		final Constructor<?> constructor=type.getConstructor(String[].class);
		try{return constructor.newInstance((Object)args);}
		catch(InvocationTargetException e){throw new AssertionError("Legacy constructor failed", e.getCause());}
	}

	/** Invokes the actual helper while keeping this test independent of either concrete class. */
	private static boolean passes(final Object owner, final Read read) throws Exception{
		final Method method=owner.getClass().getMethod("passesEditFilter", Read.class, boolean.class);
		try{return (Boolean)method.invoke(owner, read, false);}
		catch(InvocationTargetException e){throw new AssertionError("passesEditFilter failed", e.getCause());}
	}

	/** Makes the same mapped match-string read oracle as the shared processor control. */
	private static Read read(final String name, final String match, final boolean perfect){
		final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null, name, 0);
		read.match=match.getBytes(StandardCharsets.US_ASCII);
		read.setMapped(true);
		read.setPerfect(perfect);
		return read;
	}

	/** Independent literal controls must fail loudly even if assertions were accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
