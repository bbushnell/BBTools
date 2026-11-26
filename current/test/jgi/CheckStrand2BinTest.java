package test.jgi;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.security.MessageDigest;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;

import dna.Data;
import fileIO.FileFormat;
import jgi.CheckStrand2;

/** Baseline-only native bin-path helper cases for {@link CheckStrand2}.
 * @author Jean
 */
public final class CheckStrand2BinTest{

	/** Writes one literal SAM header, constructs CheckStrand2, and calls only native helpers. */
	public static void main(String[] args) throws Exception{
		Path out=null;
		Path root=Paths.get(".");
		for(String arg : args){
			if(arg.startsWith("root=")){root=Paths.get(arg.substring("root=".length()));}
			else if(arg.startsWith("outdir=")){out=Paths.get(arg.substring("outdir=".length()));}
			else{throw new IllegalArgumentException("Expected root=<BBTools root> or outdir=<new evidence directory>, not "+arg);}
		}
		root=root.toRealPath();
		require(out!=null, "Set outdir=<new evidence directory>.");
		Files.createDirectory(out);

		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run CheckStrand2BinTest with assertions enabled.");
		verifyMinorAlleleResource(root, out);

		final Path sam=out.resolve("literal.sam");
		Files.write(sam, SAM.getBytes(StandardCharsets.US_ASCII));

		final String[] constructorArgs=new String[] {
				"in="+sam, "out=null", "callreads=f", "orf=f", "polya=f",
				"t=1", "bf1=t", "bin=1000"};
		write(out.resolve("constructor.argv.tsv"), argv(constructorArgs));
		write(out.resolve("cases.tsv"), "case\texpected\tactual\tpass\n");
		final CheckStrand2 x=new CheckStrand2(constructorArgs);
		final CheckStrand2BinTest test=new CheckStrand2BinTest(x, out);
		test.run("null map -> null", null, 2, new Increment[0], null);
		test.run("empty map -> null", test.emptyMap(), 2, new Increment[0], null);
		test.run("header empty bins -> assert", test.headerMap(), 2, new Increment[0], ASSERT);
		test.run("one read default minreads2 -> assert", test.headerMap(), 2, new Increment[] {
				new Increment("ref", 0, 0)}, ASSERT);
		test.run("singleton minreads1 -> helper assert", test.headerMap(), 1, new Increment[] {
				new Increment("ref", 0, 0)}, ASSERT);
		test.run("two plus threshold2 -> ones", test.headerMap(), 2, new Increment[] {
				new Increment("ref", 0, 0), new Increment("ref", 0, 0)}, ONES);
		test.run("two plus and separate singleton -> ones", test.headerMap(), 2, new Increment[] {
				new Increment("ref", 0, 0), new Increment("ref", 0, 0),
				new Increment("ref", 1000, 0)}, ONES);

		final String summary="CheckStrand2BinTest cases="+test.cases+
				" checks="+test.checks+" failures="+test.failures;
		System.out.println(summary);
		if(test.failures!=0){throw new AssertionError(summary);}
	}

	/** Resolves the production expected-minor table before CheckStrand can fall back to simulation. */
	private static void verifyMinorAlleleResource(final Path root, final Path out) throws Exception{
		final Path expected=root.resolve("resources").resolve("minorAlleleCount.txt").toRealPath();
		final String resolved=Data.findPath("?minorAlleleCount.txt", true);
		require(resolved!=null, "Data.findPath did not resolve minorAlleleCount.txt.");
		final Path actual=Paths.get(resolved).toRealPath();
		require(actual.equals(expected),
				"Expected minorAlleleCount.txt under "+expected+" but resolved "+actual);
		final long bytes=Files.size(actual);
		final String sha80=sha80(actual);
		require(bytes==90088 && sha80.equals("6274e3df103430c1e408"),
				"Unexpected minorAlleleCount.txt pin: bytes="+bytes+", sha80="+sha80);
		write(out.resolve("resource.tsv"), "path\tbytes\tsha80\n"+actual+"\t"+bytes+"\t"+sha80+"\n");
	}

	private CheckStrand2BinTest(final CheckStrand2 x_, final Path out_) throws Exception{
		x=x_;
		out=out_;
		ffin1=(FileFormat)field("ffin1").get(x);
		binMap=field("binMap");
		minReads=field("minReads");
		binlen=field("binlen");
		incrementBinMap=CheckStrand2.class.getDeclaredMethod("incrementBinMap",
				String.class, ConcurrentHashMap.class, int.class, int.class, int.class);
		incrementBinMap.setAccessible(true);
		binStrandedness=CheckStrand2.class.getDeclaredMethod("binStrandedness");
		binStrandedness.setAccessible(true);
		final Class<?> seq=Class.forName("jgi.CheckStrand2$Seq");
		loadSamHeader=seq.getDeclaredMethod("loadSamHeader", FileFormat.class);
		loadSamHeader.setAccessible(true);
		seqLen=seq.getDeclaredField("len");
		seqLen.setAccessible(true);
	}

	/** Runs one independent literal native state and compares to a literal expected outcome. */
	private void run(final String name, final ConcurrentHashMap<?, ?> map,
			final int minReads_, final Increment[] increments, final double[] expected){
		cases++;
		String actual="not invoked";
		boolean pass=false;
		try{
			minReads.setInt(x, minReads_);
			binlen.setInt(x, 1000);
			binMap.set(x, map);
			if(map!=null){
				for(Increment inc : increments){
					incrementBinMap.invoke(x, inc.gene, map, inc.pos, inc.strand, 1);
				}
			}
			if(expected==ASSERT){
				final int failuresBefore=failures;
				actual=expectBareAssertion(name);
				pass=failures==failuresBefore;
			}else{
				final double[] observed=invokeBinStrandedness();
				actual=string(observed);
				pass=expected==null ? observed==null : Arrays.equals(expected, observed);
				check(pass, name+" expected "+string(expected)+" but observed "+actual);
			}
		}catch(Throwable t){
			actual=t.getClass().getSimpleName()+": "+t.getMessage();
			fail(name+" unexpected "+t.getClass().getSimpleName()+": "+t.getMessage());
			t.printStackTrace(System.out);
		}finally{
			append(out.resolve("cases.tsv"), name+"\t"+expected(expected)+"\t"+actual+"\t"+pass+"\n");
		}
	}

	/** Verifies only the expected InvocationTargetException wrapping a bare AssertionError. */
	private String expectBareAssertion(final String name) throws Exception{
		checks++;
		try{
			binStrandedness.invoke(x);
			fail(name+" expected the native CheckStrand.strandedness assertion.");
		}catch(InvocationTargetException ite){
			final Throwable target=ite.getTargetException();
			final StackTraceElement[] stack=target.getStackTrace();
			if(!(target instanceof AssertionError) || target.getMessage()!=null ||
					stack.length<1 || !stack[0].getClassName().equals("jgi.CheckStrand") ||
					!stack[0].getMethodName().equals("strandedness")){
				final String actual=target.getClass().getSimpleName()+": "+target.getMessage()+" at "+
						(stack.length<1 ? "empty stack" : stack[0].toString());
				fail(name+" expected a bare AssertionError at jgi.CheckStrand.strandedness but observed "+
						actual);
				return actual;
			}
			return "bare AssertionError at jgi.CheckStrand.strandedness";
		}
		return "no exception";
	}

	private double[] invokeBinStrandedness() throws Exception{
		return (double[])binStrandedness.invoke(x);
	}

	private ConcurrentHashMap<?, ?> emptyMap(){
		return new ConcurrentHashMap<Object, Object>();
	}

	private ConcurrentHashMap<?, ?> headerMap() throws Exception{
		final ConcurrentHashMap<?, ?> map=(ConcurrentHashMap<?, ?>)loadSamHeader.invoke(null, ffin1);
		require(map.size()==1 && map.containsKey("ref"),
				"Literal SAM header should load exactly one ref sequence: "+map.keySet());
		require(seqLen.getInt(map.get("ref"))==2000, "Literal SAM ref length should be 2000.");
		return map;
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

	private static String string(final double[] array){
		return array==null ? "null" : Arrays.toString(array);
	}

	private static String expected(final double[] array){
		return array==ASSERT ? "bare AssertionError at jgi.CheckStrand.strandedness" : string(array);
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
			Files.write(path, text.getBytes(StandardCharsets.UTF_8), java.nio.file.StandardOpenOption.APPEND);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	private static String sha80(final Path path) throws Exception{
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final byte[] hash=digest.digest(Files.readAllBytes(path));
		final StringBuilder sb=new StringBuilder(hash.length*2);
		for(byte b : hash){sb.append(String.format("%02x", b&0xff));}
		return sb.substring(sb.length()-20);
	}

	private static final class Increment{
		Increment(final String gene_, final int pos_, final int strand_){
			gene=gene_;
			pos=pos_;
			strand=strand_;
		}
		final String gene;
		final int pos;
		final int strand;
	}

	private static final String SAM="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:2000\n";
	private static final double[] ASSERT=new double[0];
	private static final double[] ONES=new double[] {1.0, 1.0, 1.0, 1.0};

	private final CheckStrand2 x;
	private final Path out;
	private final FileFormat ffin1;
	private final Field binMap;
	private final Field minReads;
	private final Field binlen;
	private final Field seqLen;
	private final Method incrementBinMap;
	private final Method binStrandedness;
	private final Method loadSamHeader;
	private int cases;
	private int checks;
	private int failures;
}
