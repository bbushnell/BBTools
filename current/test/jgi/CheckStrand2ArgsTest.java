package test.jgi;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.FileFormat;
import jgi.CheckStrand2;

/**
 * Constructor-only filename argument cases for {@link CheckStrand2}.
 * @author Jean
 */
public final class CheckStrand2ArgsTest{

	/** Writes literal SAM inputs, constructs CheckStrand2, and verifies parsed constructor fields. */
	public static void main(String[] args) throws Exception{
		Path out=null;
		String phase=null;
		Path root=Paths.get(".");
		for(String arg : args){
			if(arg.startsWith("outdir=")){out=Paths.get(arg.substring("outdir=".length()));}
			else if(arg.startsWith("phase=")){phase=arg.substring("phase=".length());}
			else if(arg.startsWith("root=")){root=Paths.get(arg.substring("root=".length()));}
			else{throw new IllegalArgumentException("Unknown argument "+arg);}
		}
		root=root.toRealPath();
		require(out!=null, "Set outdir=<new evidence directory>.");
		require(phase!=null && (phase.equals("baseline") || phase.equals("fixed")),
				"Set phase=baseline or phase=fixed.");
		Files.createDirectory(out);
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run CheckStrand2ArgsTest with assertions enabled.");
		write(out.resolve("cases.tsv"), "case\tphase\texpected\tactual\tpass\n");
		write(out.resolve("root.txt"), root+"\n");

		final CheckStrand2ArgsTest test=new CheckStrand2ArgsTest(out, phase);
		test.run("ordinary input out=null", "ordinary_null", "input.sam", "out=null", null, true);
		test.run("ordinary input out=", "ordinary_empty", "input.sam", "out=", null, true);
		test.run("equals input out=null", "equals_null", "reads=one.sam", "out=null", "reads",
				phase.equals("fixed"));

		final String summary="CheckStrand2ArgsTest phase="+phase+" cases="+test.cases+
				" checks="+test.checks+" failures="+test.failures;
		System.out.println(summary);
		if(test.failures!=0){throw new AssertionError(summary);}
	}

	private CheckStrand2ArgsTest(final Path out_, final String phase_) throws Exception{
		out=out_;
		phase=phase_;
		in1Field=field("in1");
		out1Field=field("out1");
		ffin1Field=field("ffin1");
		ffout1Field=field("ffout1");
	}

	/** Runs one independent constructor case and compares exact retained paths or exception details. */
	private void run(final String name, final String slug, final String inputName, final String outArg,
			final String absentPrefix, final boolean shouldConstruct) throws Exception{
		cases++;
		final Path caseDir=out.resolve(slug);
		Files.createDirectory(caseDir);
		final Path input=caseDir.resolve(inputName);
		write(input, SAM);

		final Path prefix=absentPrefix==null ? null : caseDir.resolve(absentPrefix);
		require(Files.exists(input) && Files.isReadable(input),
				"Full literal input is not readable: "+input);
		if(prefix!=null){
			require(!Files.exists(prefix), "Absent prefix exists before constructor: "+prefix);
		}

		final String[] constructorArgs=new String[] {
				"in="+input, outArg, "callreads=f", "orf=f", "polya=f", "t=1", "bf1=t"};
		write(caseDir.resolve("constructor.argv.tsv"), argv(constructorArgs));
		write(caseDir.resolve("input.tsv"),
				"full\texists\treadable\tabsentPrefix\texists\n"+
				input+"\t"+Files.exists(input)+"\t"+Files.isReadable(input)+"\t"+
				(prefix==null ? "n/a" : prefix.toString())+"\t"+
				(prefix==null ? "n/a" : Boolean.toString(Files.exists(prefix)))+"\n");

		final String expected=shouldConstruct ? success(input) : failure(prefix);
		String actual="not invoked";
		boolean pass=false;
		try{
			final CheckStrand2 x=new CheckStrand2(constructorArgs);
			actual=actual(x);
			pass=shouldConstruct && actual.equals(expected);
			check(pass, name+" expected "+expected+" but observed "+actual);
		}catch(RuntimeException re){
			actual=actual(re);
			pass=!shouldConstruct && actual.equals(expected);
			check(pass, name+" expected "+expected+" but observed "+actual);
		}catch(Throwable t){
			actual=t.getClass().getSimpleName()+": "+t.getMessage();
			fail(name+" unexpected "+actual);
			t.printStackTrace(System.out);
		}finally{
			append(out.resolve("cases.tsv"), name+"\t"+phase+"\t"+expected+"\t"+actual+"\t"+pass+"\n");
		}
	}

	private String actual(final CheckStrand2 x) throws Exception{
		final FileFormat ffin1=(FileFormat)ffin1Field.get(x);
		final Object out1=out1Field.get(x);
		final Object ffout1=ffout1Field.get(x);
		return "constructed in1="+in1Field.get(x)+" ffin1="+
				(ffin1==null ? null : ffin1.name())+
				" out1="+out1+" out1IsNull="+(out1==null)+
				" ffout1="+ffout1+" ffout1IsNull="+(ffout1==null);
	}

	private static String success(final Path input){
		return "constructed in1="+input+" ffin1="+input+
				" out1=null out1IsNull=true ffout1=null ffout1IsNull=true";
	}

	private static String actual(final RuntimeException re){
		final StackTraceElement[] stack=re.getStackTrace();
		return re.getClass().getSimpleName()+": "+re.getMessage()+" at "+
				(stack.length<1 ? "empty stack" :
					stack[0].getClassName()+"."+stack[0].getMethodName());
	}

	private static String failure(final Path prefix){
		return "RuntimeException: Can't read file '"+prefix+"' at shared.Tools.testInputFiles";
	}

	private static Field field(final String name) throws Exception{
		final Field f=CheckStrand2.class.getDeclaredField(name);
		f.setAccessible(true);
		return f;
	}

	private void check(final boolean pass, final String message){
		checks++;
		if(!pass){fail(message);}
	}

	private void fail(final String message){
		failures++;
		System.out.println("FAIL "+message);
	}

	private static void require(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}

	private static String argv(final String[] args){
		final StringBuilder sb=new StringBuilder("index\targ\n");
		for(int i=0; i<args.length; i++){sb.append(i).append('\t').append(args[i]).append('\n');}
		return sb.toString();
	}

	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}

	private static void append(final Path path, final String text){
		try{
			Files.write(path, text.getBytes(StandardCharsets.UTF_8),
					java.nio.file.StandardOpenOption.APPEND);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	private static final String SAM="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:2000\n";

	private final Path out;
	private final String phase;
	private final Field in1Field;
	private final Field out1Field;
	private final Field ffin1Field;
	private final Field ffout1Field;
	private int cases;
	private int checks;
	private int failures;
}
