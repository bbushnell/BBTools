package test.jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteStreamWriter;
import jgi.AssemblyStats3;
import parse.Parse;
import shared.Shared;

/**
 * Exercises AssemblyStats3 cleanup after unchecked Assembly construction failure
 * and after a latched ByteStreamWriter error, plus their public-main controls.
 * Each row runs in a fresh child JVM. A failure child records only the matching
 * ByteStreamWriter survivor, then poisons that survivor before asserting.
 * @author Jean
 */
public final class AssemblyStats3CleanupTest{

	/** Parent requires phase=baseline|fixed plus a fresh outdir; child is private to this fixture. */
	public static void main(final String[] args) throws Exception{
		String mode="parent", phase=null, directory=null, kind=null, output=null, work=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("kind")){kind=pair[1];}
			else if(pair[0].equals("output")){output=pair[1];}
			else if(pair[0].equals("work")){work=pair[1];}
			else{throw new IllegalArgumentException("Unknown AssemblyStats3 cleanup option: "+arg);}
		}
		if("child".equals(mode)){
			child(kind, output, work);
			return;
		}
		require("parent".equals(mode), "Set mode=parent|child.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/assemblystats3cleanuptest.sh with assertions enabled.");
		require(Shared.threads()==1, "The cleanup fixture requires one tool thread.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyStats3CleanupTest test=new AssemblyStats3CleanupTest(root, phase.equals("baseline"));
		test.scenario("clean-main", false, false);
		test.scenario("writer-error", false, true);
		test.scenario("unchecked-main", true, false);
		test.scenario("unchecked-writer-error", true, true);
		test.check(test.cases==4 && test.children==4,
			"Expected four cases and four children: "+test.cases+"/"+test.children);
		test.check(test.uncheckedSurvivors==(test.baseline ? 2 : 0),
			"Unexpected unchecked survivor total: "+test.uncheckedSurvivors);
		test.check(test.writerErrors==2, "Expected two writer-error rows: "+test.writerErrors);
		final String summary="ASSEMBLY_STATS3_CLEANUP_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" uncheckedSurvivors="+test.uncheckedSurvivors+" writerErrors="+test.writerErrors+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	private AssemblyStats3CleanupTest(final Path root_, final boolean baseline_){
		root=root_; baseline=baseline_;
	}

	/** Runs one fresh child JVM and checks its retained evidence. */
	private void scenario(final String name, final boolean inject, final boolean writerError) throws Exception{
		final Path work=Files.createDirectory(root.resolve(name));
		final Path output=work.resolve("stats3.tsv");
		final int exit=run(work, kind(inject, writerError));
		final String observed=read(work.resolve("observed.tsv"));
		final String[] values=observed.split("\n", -1)[1].split("\t", -1);
		final int bulk=Integer.parseInt(values[1]), single=Integer.parseInt(values[2]);
		final int delivered=Integer.parseInt(values[3]), errors=Integer.parseInt(values[4]);
		final int writerMatches=Integer.parseInt(values[5]), before=Integer.parseInt(values[6]), after=Integer.parseInt(values[7]);
		final boolean restored=Boolean.parseBoolean(values[8]);
		final boolean unchanged=Boolean.parseBoolean(values[9]);
		final boolean sameInjected=Boolean.parseBoolean(values[10]);
		final int rescueErrorFlags=Integer.parseInt(values[11]);
		final String errorClass=values[12], errorMessage=values[13];
		final int suppressedCount=Integer.parseInt(values[14]);
		final String suppressedClass=values[15], suppressedMessage=values[16], cleanup=values[17];
		if(writerError){writerErrors++;}
		if(inject){uncheckedSurvivors+=before;}
		check(exit==0, name+": child must exit naturally after rescuing writers: "+exit);
		check(bulk==2 && single==0 && delivered==PREFIX.length,
			name+": expected two bulk reads and one literal prefix: "+bulk+"/"+single+"/"+delivered);
		check(errors==(inject ? 1 : 0), name+": unexpected injected error count: "+errors);
		check(writerMatches==(writerError ? 1 : 0), name+": unexpected BSW error-latch targets: "+writerMatches);
		check(before==(inject && baseline ? 1 : 0),
			name+": unexpected matching writer survivors before rescue: "+before);
		check(after==0, name+": matching writer survived harness cleanup: "+after);
		check(rescueErrorFlags==(inject && writerError && baseline ? 1 : 0),
			name+": unexpected harness rescue error flags: "+rescueErrorFlags);
		check(restored && unchanged, name+": fixture state or source bytes changed.");
		check("none".equals(cleanup), name+": writer rescue failed: "+cleanup);
		check(Files.isRegularFile(output), name+": child did not create raw Stats3 output.");
		check(Arrays.equals(Files.readAllBytes(output), (inject ? HEADER_ONLY : CLEAN_OUTPUT).getBytes(StandardCharsets.US_ASCII)),
			name+": unexpected Stats3 output bytes.");
		if(inject){
			check(sameInjected && errorClass.equals(InjectedReadException.class.getName()) && ERROR.equals(errorMessage),
				name+": unchecked read failure did not escape intact: "+sameInjected+" / "+errorClass+" / "+errorMessage);
			checkSuppressed(name, writerError && !baseline, suppressedCount, suppressedClass, suppressedMessage, output);
		}else if(writerError && !baseline){
			check(errorClass.equals(RuntimeException.class.getName()) && errorMessage.contains(output.toString()),
				name+": fixed writer failure did not name its output: "+errorClass+" / "+errorMessage);
			checkSuppressed(name, false, suppressedCount, suppressedClass, suppressedMessage, output);
		}else{
			check("none".equals(errorClass) && errorMessage.isEmpty(), name+": clean case threw "+errorClass+" / "+errorMessage);
			checkSuppressed(name, false, suppressedCount, suppressedClass, suppressedMessage, output);
		}
		cases++;
	}

	/** Launches this same fixture through its dev launcher, with all child evidence retained. */
	private int run(final Path work, final String kind) throws Exception{
		final Path launcher=launcher();
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", launcher.toString(),
			"mode=child", "kind="+kind, "output="+work.resolve("stats3.tsv"), "work="+work));
		Files.write(work.resolve("command.argv"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(work.toFile())
			.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile());
		final Process process=builder.start();
		process.getOutputStream().close();
		final int exit=process.waitFor();
		Files.write(work.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		children++;
		return exit;
	}

	private static String kind(final boolean inject, final boolean writerError){
		return inject ? (writerError ? "unchecked-writer" : "unchecked") : (writerError ? "writer" : "clean");
	}

	/** Child mode calls the public driver under forced BF1 and always rescues matching writers. */
	private static void child(final String kind, final String output, final String work) throws Exception{
		require(kind!=null && output!=null && work!=null, "Child requires kind/output/work.");
		final boolean inject=kind.startsWith("unchecked"), writerError=kind.endsWith("writer");
		require(kind.equals(kind(inject, writerError)), "Unknown child kind: "+kind);
		final Path root=Paths.get(work).toAbsolutePath().normalize();
		require(Shared.threads()==1, "Child cleanup fixture requires one tool thread.");
		require(!new File(INPUT).exists(), "A real stdin.fa would override borrowed standard input.");
		require(ByteFile1.bufferlen>=PREFIX.length, "The literal prefix must fit in BF1's first buffer.");
		final InputStream originalIn=System.in;
		final boolean originalForce=ByteFile.FORCE_MODE_BF1, originalSIMD=Shared.SIMD, originalBuffered=ByteFile1.BUFFERED;
		final TrackedInput input=new TrackedInput(inject, writerError, output);
		Throwable error=null;
		Rescue rescue=Rescue.NONE;
		String matchingBefore="fname\tstate\tdaemon\talive\terrorState\n";
		int before=-1, after=-1;
		boolean restored=false;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			Shared.SIMD=false;
			ByteFile1.BUFFERED=false;
			System.setIn(input);
			try{AssemblyStats3.main(new String[]{"in="+INPUT, "out="+output});}
			catch(final Throwable caught){error=caught;}
		}finally{
			try{
				final List<ByteStreamWriter> writers=matchingWriters(output);
				before=writers.size();
				matchingBefore=describe(writers);
				rescue=poisonAll(writers);
			}finally{
				System.setIn(originalIn);
				ByteFile.FORCE_MODE_BF1=originalForce;
				Shared.SIMD=originalSIMD;
				ByteFile1.BUFFERED=originalBuffered;
				restored=System.in==originalIn && ByteFile.FORCE_MODE_BF1==originalForce &&
					Shared.SIMD==originalSIMD && ByteFile1.BUFFERED==originalBuffered;
				after=matchingWriters(output).size();
			}
		}
		write(root.resolve("matching-writers-before.tsv"), matchingBefore);
		final Throwable[] suppressed=error==null ? new Throwable[0] : error.getSuppressed();
		write(root.resolve("observed.tsv"), "kind\tbulk\tsingle\tdelivered\terrors\twriterMatches"+
			"\tmatchedBefore\tmatchedAfter\tglobalsRestored\tprefixUnchanged\tsameInjected\trescueErrorFlags"+
			"\terrorClass\terrorMessage\tsuppressedCount\tsuppressedClass\tsuppressedMessage\tcleanupFailure\n"+
			kind+"\t"+input.bulkCalls+"\t"+input.singleCalls+"\t"+input.delivered+"\t"+input.errors+
			"\t"+input.writerMatches+"\t"+before+"\t"+after+"\t"+restored+"\t"+Arrays.equals(PREFIX, input.bytes)+
			"\t"+(error==input.injected)+"\t"+rescue.errorFlags+"\t"+(error==null ? "none" : error.getClass().getName())+"\t"+
			clean(error==null ? "" : String.valueOf(error.getMessage()))+"\t"+suppressed.length+"\t"+
			(suppressed.length==0 ? "none" : suppressed[0].getClass().getName())+"\t"+
			clean(suppressed.length==0 ? "" : String.valueOf(suppressed[0].getMessage()))+"\t"+
			(rescue.failure==null ? "none" : rescue.failure.getClass().getName()+":"+
				clean(String.valueOf(rescue.failure.getMessage())))+"\n");
		if(error!=null){write(root.resolve("primary-stack.txt"), stack(error));}
		if(rescue.failure!=null){write(root.resolve("cleanup-stack.txt"), stack(rescue.failure));}
		if(rescue.failure!=null){
			final AssertionError failure=new AssertionError("Harness writer rescue failed");
			failure.initCause(rescue.failure);
			throw failure;
		}
	}

	/** Returns only still-live ByteStreamWriter threads for this case's unique output path. */
	private static ArrayList<ByteStreamWriter> matchingWriters(final String output){
		final ArrayList<ByteStreamWriter> writers=new ArrayList<ByteStreamWriter>();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter){
				final ByteStreamWriter writer=(ByteStreamWriter)thread;
				if(thread.isAlive() && output.equals(writer.fname)){writers.add(writer);}
			}
		}
		return writers;
	}

	private void checkSuppressed(final String name, final boolean expected, final int count,
			final String className, final String message, final Path output){
		check(count==(expected ? 1 : 0), name+": unexpected suppressed cleanup failure count: "+count);
		if(expected){
			check(className.equals(RuntimeException.class.getName()) && message.contains(output.toString()),
				name+": suppressed cleanup failure did not name its output: "+className+" / "+message);
		}else{
			check("none".equals(className) && message.isEmpty(),
				name+": unexpected suppressed cleanup failure: "+className+" / "+message);
		}
	}

	/** Attempts every matched writer rescue and retains any cleanup failures. */
	private static Rescue poisonAll(final List<ByteStreamWriter> writers){
		Throwable failure=null;
		int errorFlags=0;
		for(final ByteStreamWriter writer : writers){
			try{
				if(writer.poisonAndWait()){errorFlags++;}
			}
			catch(final Throwable caught){
				if(failure==null){failure=caught;}
				else if(caught!=failure){failure.addSuppressed(caught);}
			}
		}
		return failure==null && errorFlags==0 ? Rescue.NONE : new Rescue(failure, errorFlags);
	}

	private static final class Rescue{
		private static final Rescue NONE=new Rescue(null, 0);

		private Rescue(final Throwable failure_, final int errorFlags_){
			failure=failure_;
			errorFlags=errorFlags_;
		}

		final Throwable failure;
		final int errorFlags;
	}

	/** Locates the dev launcher from this fixture's current/ classpath. */
	private static Path launcher(){
		String cp=System.getProperty("java.class.path");
		final int split=cp.indexOf(File.pathSeparatorChar);
		if(split>=0){cp=cp.substring(0, split);}
		return Paths.get(cp).toAbsolutePath().normalize().getParent().resolve("dev/assemblystats3cleanuptest.sh");
	}

	private static String describe(final List<ByteStreamWriter> writers){
		final StringBuilder sb=new StringBuilder("fname\tstate\tdaemon\talive\terrorState\n");
		for(final ByteStreamWriter writer : writers){
			sb.append(writer.fname).append('\t').append(writer.getState()).append('\t')
				.append(writer.isDaemon()).append('\t').append(writer.isAlive()).append('\t')
				.append(writer.errorState).append('\n');
		}
		return sb.toString();
	}

	private static String stack(final Throwable t){
		final java.io.StringWriter out=new java.io.StringWriter();
		t.printStackTrace(new java.io.PrintWriter(out));
		return out.toString();
	}
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	private static String clean(final String value){return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');}
	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	/** Supplies one complete contig, then EOF or an unchecked failure on the next bulk read. */
	private static final class TrackedInput extends InputStream{
		TrackedInput(final boolean inject_, final boolean writerError_, final String output_){
			inject=inject_; writerError=writerError_; output=output_;
		}
		@Override
		public int read() throws IOException{
			singleCalls++;
			throw new IOException("Unexpected single-byte Stats3 fixture access");
		}
		@Override
		public int read(final byte[] target, final int offset, final int length) throws IOException{
			if(target==null){throw new NullPointerException();}
			if(offset<0 || length<0 || offset>target.length-length){throw new IndexOutOfBoundsException();}
			if(length==0){return 0;}
			bulkCalls++;
			if(delivered<bytes.length){
				final int count=Math.min(length, bytes.length-delivered);
				System.arraycopy(bytes, delivered, target, offset, count);
				delivered+=count;
				return count;
			}
			if(writerError && !latched){
				final List<ByteStreamWriter> writers=matchingWriters(output);
				writerMatches=writers.size();
				for(final ByteStreamWriter writer : writers){writer.errorState=true;}
				latched=true;
			}
			if(inject){errors++; throw injected;}
			return -1;
		}
		final byte[] bytes=PREFIX.clone();
		final boolean inject, writerError;
		final String output;
		final InjectedReadException injected=new InjectedReadException(ERROR);
		boolean latched;
		int bulkCalls, singleCalls, delivered, errors, writerMatches;
	}

	private static final class InjectedReadException extends RuntimeException{
		InjectedReadException(final String message){super(message);}
		private static final long serialVersionUID=1L;
	}

	private final Path root;
	private final boolean baseline;
	private int cases, children, uncheckedSurvivors, writerErrors, checks, failures;
	private static final String INPUT="stdin.fa";
	private static final byte[] PREFIX=">ok\nAC\n".getBytes(StandardCharsets.US_ASCII);
	private static final String ERROR="Injected unchecked Stats3 Assembly failure";
	private static final String HEADER_ONLY="fname\tsize\tcontigs\tgc\tmaxContig\t5kplus\t10kplus\t25kplus\t50kplus\n";
	private static final String CLEAN_OUTPUT=HEADER_ONLY+"stdin.fa\t2\t1\t0.500\t2\t0\t0\t0\t0\n";
}
