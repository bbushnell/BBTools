package test.pacbio;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import pacbio.ProcessStackedSitesNormalized;
import pacbio.ProcessStackedSitesNormalized.Buffer;
import shared.Shared;
import stream.SiteScoreR;

/** Bounded no-NaN and short-buffer controls for ProcessStackedSitesNormalized.
 * @author Brian Bushnell, Jean
 */
public final class ProcessStackedSitesNormalizedBoundaryTest{

	/** Runs finite score, empty-buffer, one-group and three-group controls.
	 * @param args Phase and fresh absolute retained-evidence directory
	 * @throws Exception if reflection or retained evidence files fail
	 */
	public static void main(String[] args) throws Exception{
		String phase="baseline";
		String child=null;
		String outdir=null;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length==2 && "phase".equalsIgnoreCase(split[0])){phase=split[1].toLowerCase();}
			else if(split.length==2 && "child".equalsIgnoreCase(split[0])){child=split[1].toLowerCase();}
			else if(split.length==2 && "outdir".equalsIgnoreCase(split[0])){outdir=split[1];}
			else{throw new IllegalArgumentException(
					"Expected phase=baseline|fixed, child=name or outdir=path, not "+arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Unknown phase: "+phase);
		}
		if(outdir==null){throw new IllegalArgumentException("Missing outdir=path");}
		OUTDIR=new File(outdir);
		File parent=OUTDIR.getParentFile();
		if(!OUTDIR.isAbsolute() || OUTDIR.exists() || parent==null ||
				!parent.isDirectory() || !OUTDIR.mkdir()){
			throw new IllegalArgumentException(
					"outdir must be a fresh absolute directory: "+OUTDIR);
		}
		if(child!=null){
			runChild(child);
			return;
		}

		int failures=0;
		failures+=checkScores(phase);
		failures+=checkIntervalZero(phase);
		failures+=checkEmptyBuffer();
		failures+=checkOneGroup(phase);
		failures+=checkFixedCliChildren(phase);
		failures+=checkThreeGroupCli();
		if(failures!=0){throw new AssertionError(failures+" PSSN boundary failure(s)");}
	}

	/** Captures current finite, NaN, infinity and readlen-zero score behavior. */
	private static int checkScores(String phase){
		reset();
		int failures=0;
		failures+=checkScore("normal", site(0, 23, 24), 6, "101.5", phase);
		failures+=checkScore("equal-span", site(0, 11, 12), 0,
				"baseline".equals(phase) ? "nan" : "throw", phase);
		failures+=checkScore("short-span", site(0, 10, 11), 1,
				"baseline".equals(phase) ? "infinite" : "throw", phase);
		failures+=checkScore("readlen0", site(0, 23, 0), 6,
				"baseline".equals(phase) ? "finite" : "throw", phase);
		// Old int multiplication overflows and clamps the read-length term low;
		// the fixed float multiply clamps it high and still keeps the score finite.
		failures+=checkScore("large-readlen", site(0, 23, Integer.MAX_VALUE), 6,
				"baseline".equals(phase) ? "100.7" : "102.3", phase);
		return failures;
	}

	/** Checks one normalizedScore classification, literal, or expected exception. */
	private static int checkScore(String name, SiteScoreR ssr, int endDist,
			String expected, String phase){
		try{
			float observed=ProcessStackedSitesNormalized.normalizedScore(ssr, endDist);
			String type=Float.isNaN(observed) ? "nan" :
				(Float.isInfinite(observed) ? "infinite" : "finite");
			String literal=String.valueOf(observed);
			boolean pass=expected.equals(type) || expected.equals(literal);
			System.out.println("score-"+name+'\t'+phase+'\t'+expected+'\t'+type+
					'\t'+literal+'\t'+result(pass));
			return pass ? 0 : 1;
		}catch(IllegalArgumentException e){
			boolean pass="throw".equals(expected);
			System.out.println("score-"+name+'\t'+phase+'\t'+expected+
					"\tthrow\t"+e.getClass().getSimpleName()+'\t'+result(pass));
			return pass ? 0 : 1;
		}
	}

	/** Records the old interval=0 liveness hazard without running it. */
	private static int checkIntervalZero(String phase) throws Exception{
		File in=file("interval-zero-input.txt");
		File out=file("interval-zero-output.txt");
		write(in, site(0, 23, 24).toText().toString(),
				site(1, 24, 24).toText().toString(),
				site(2, 25, 24).toText().toString());
		if("baseline".equals(phase)){
			System.out.println("process-interval-zero	baseline	source-supported	not-run	PASS");
			return 0;
		}
		reset();
		ProcessStackedSitesNormalized.INTERVAL=0;
		try{
			ProcessStackedSitesNormalized.process(in.getPath(), out.getPath());
			System.out.println("process-interval-zero	fixed	throw	none	FAIL");
			return 1;
		}catch(IllegalArgumentException e){
			System.out.println("process-interval-zero	fixed	throw	"+
					e.getClass().getSimpleName()+"\tPASS");
			return 0;
		}
	}

	/** Empty input may construct a Buffer, but must be closed without entering process(). */
	private static int checkEmptyBuffer() throws Exception{
		reset();
		File in=file("empty-input.txt");
		File out=file("empty-output.txt");
		write(in);
		Buffer buffer=new Buffer(3, in.getPath(), out.getPath());
		try{
			boolean pass=buffer.array[0]==null && buffer.chrom==-1 &&
					buffer.min==-1 && buffer.max==-1;
			System.out.println("buffer-empty\tbaseline\tnull-array0\tchrom="+
					buffer.chrom+",min="+buffer.min+",max="+buffer.max+'\t'+
					result(pass));
			return pass ? 0 : 1;
		}finally{
			closeAndJoin(buffer);
		}
	}

	/** One stacked group leaves null slots; baseline processInterval dereferences one. */
	private static int checkOneGroup(String phase) throws Exception{
		reset();
		File in=file("one-group-input.txt");
		File out=file("one-group-output.txt");
		write(in, site(0, 23, 24).toText().toString());
		Buffer buffer=new Buffer(3, in.getPath(), out.getPath());
		try{
			Throwable thrown=invokeProcessInterval(buffer);
			String observed=thrown==null ? "none" : thrown.getClass().getSimpleName();
			String expected="baseline".equals(phase) ? "NullPointerException" : "none";
			boolean pass=expected.equals(observed);
			System.out.println("buffer-one-group\t"+phase+'\t'+expected+'\t'+
					observed+'\t'+result(pass));
			return pass ? 0 : 1;
		}finally{
			closeAndJoin(buffer);
		}
	}

	/** Three stacked groups fill the old fixed window and write the retained group. */
	private static int checkThreeGroupCli() throws Exception{
		reset();
		File in=file("three-group-input.txt");
		File out=file("three-group-output.txt");
		String expected="0,0,2,25,24,3,0,00,0,100\n";
		write(in, site(0, 23, 24).toText().toString(),
				site(1, 24, 24).toText().toString(),
				site(2, 25, 24).toText().toString());
		try{
			ProcessStackedSitesNormalized.main(new String[] {in.getPath(), out.getPath(),
					"retainall=t", "minvotestoretain=1"});
			String observed=read(out);
			boolean pass=expected.equals(observed);
			System.out.println("cli-three-group\tbaseline\t"+expected.trim()+
					'\t'+observed.trim()+'\t'+result(pass));
			return pass ? 0 : 1;
		}catch(Throwable t){
			System.out.println("cli-three-group\tbaseline\t"+expected.trim()+'\t'+
					t.getClass().getSimpleName()+"\tFAIL");
			return 1;
		}
	}


	/** Runs fixed CLI child checks without touching the old hanging baseline paths. */
	private static int checkFixedCliChildren(String phase) throws Exception{
		if("baseline".equals(phase)){
			System.out.println("child-empty-main\tbaseline\tnot-run\tnot-run\tPASS");
			System.out.println("child-invalid-score\tbaseline\tnot-run\tnot-run\tPASS");
			System.out.println("child-invalid-setting\tbaseline\tnot-run\tnot-run\tPASS");
			return 0;
		}
		int failures=0;
		failures+=checkChild("empty-main", 0, "Sites In:\t0    \t0.000% correct",
				"Sites Out:\t0    \t0.000% correct", true, true);
		failures+=checkChild("invalid-score", 1, "IllegalArgumentException",
				"ProcessStackedSitesNormalized failed", false, false);
		failures+=checkChild("invalid-setting", 1, "IllegalArgumentException",
				"centerweight must be finite", false, false);
		return failures;
	}

	/** Checks one child JVM for expected exit and diagnostic text. */
	private static int checkChild(String child, int expectedExit, String required,
			String secondRequired, boolean checkEmptyOutput,
			boolean rejectNonfiniteOutput) throws Exception{
		File dir=file(child);
		ChildResult result=runChildJvm(child, dir);
		String combined=result.stdout+'\n'+result.stderr;
		File out=new File(new File(dir, "work"), child+"-output.txt");
		boolean outputOk=!checkEmptyOutput || (out.exists() && out.length()==0);
		boolean finiteOk=!rejectNonfiniteOutput ||
				(!combined.contains("NaN") && !combined.contains("Infinity"));
		boolean pass=result.exit==expectedExit && combined.contains(required) &&
				combined.contains(secondRequired) && outputOk && finiteOk;
		System.out.println("child-"+child+"\tfixed\texit"+expectedExit+"\texit"+
				result.exit+'\t'+result.status+'\t'+result(pass));
		if(!pass){
			System.out.println("child-"+child+"-stdout\t"+safe(result.stdout));
			System.out.println("child-"+child+"-stderr\t"+safe(result.stderr));
		}
		return pass ? 0 : 1;
	}

	/** Starts one child JVM and captures its bounded output. */
	private static ChildResult runChildJvm(String child, File dir) throws Exception{
		if(!dir.mkdir()){throw new AssertionError("Could not create child evidence dir "+dir);}
		File work=new File(dir, "work");
		String wrapper=System.getProperty("pssn.wrapper");
		if(wrapper==null || wrapper.length()<1){
			throw new IllegalStateException("Missing -Dpssn.wrapper path");
		}
		ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(wrapper);
		command.add("phase=fixed");
		command.add("child="+child);
		command.add("outdir="+work.getAbsolutePath());
		File stdout=new File(dir, "stdout.txt");
		File stderr=new File(dir, "stderr.txt");
		File exit=new File(dir, "exit.txt");
		write(new File(dir, "command.txt"), command.toArray(new String[0]));
		Process process=new ProcessBuilder(command).
				redirectOutput(stdout).redirectError(stderr).start();
		boolean finished=process.waitFor(10, TimeUnit.SECONDS);
		if(!finished){
			process.destroyForcibly();
			process.waitFor();
		}
		int exitValue=finished ? process.exitValue() : -1;
		String status=finished ? "finished" : "timeout";
		write(exit, String.valueOf(exitValue), status);
		return new ChildResult(exitValue, status, read(stdout), read(stderr));
	}

	/** Executes one child case inside a fresh JVM. */
	private static void runChild(String child) throws Exception{
		reset();
		File in=file(child+"-input.txt");
		File out=file(child+"-output.txt");
		if("empty-main".equals(child)){
			write(in);
			ProcessStackedSitesNormalized.main(new String[] {in.getPath(), out.getPath()});
		}else if("invalid-score".equals(child)){
			write(in, site(0, 11, 12).toText().toString());
			ProcessStackedSitesNormalized.main(new String[] {in.getPath(), out.getPath(),
					"retainall=t", "minvotestoretain=1"});
		}else if("invalid-setting".equals(child)){
			write(in, site(0, 23, 24).toText().toString());
			ProcessStackedSitesNormalized.main(new String[] {in.getPath(), out.getPath(),
					"centerweight=1e999"});
		}else{
			throw new IllegalArgumentException("Unknown child: "+child);
		}
	}

	/** Flattens captured output for a single diagnostic row. */
	private static String safe(String s){
		return s.replace('\t', ' ').replace('\n', '|');
	}

	/** Invokes the private interval worker and unwraps the target exception. */
	private static Throwable invokeProcessInterval(Buffer buffer) throws Exception{
		Method method=ProcessStackedSitesNormalized.class.getDeclaredMethod(
				"processInterval", Buffer.class, int.class, int.class, int.class);
		method.setAccessible(true);
		try{
			method.invoke(null, buffer, buffer.chrom, buffer.min,
					buffer.min+ProcessStackedSitesNormalized.INTERVAL-1);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Closes a Buffer and joins its reflected TextStreamWriter. */
	private static void closeAndJoin(Buffer buffer) throws Exception{
		Field field=Buffer.class.getDeclaredField("tsw");
		field.setAccessible(true);
		TextStreamWriter tsw=(TextStreamWriter)field.get(buffer);
		buffer.close();
		tsw.join(10000);
		if(tsw.isAlive()){throw new AssertionError("TextStreamWriter did not finish");}
	}

	/** Restores static counters and knobs used by this fixture. */
	private static void reset(){
		ProcessStackedSitesNormalized.INTERVAL=12;
		ProcessStackedSitesNormalized.MIN_LENGTH_TO_RETAIN=0;
		ProcessStackedSitesNormalized.RETAIN_ALL=false;
		ProcessStackedSitesNormalized.MIN_SITES_TO_DISCARD=0;
		ProcessStackedSitesNormalized.SITES_TO_RETAIN1=8;
		ProcessStackedSitesNormalized.SITES_TO_RETAIN2=16;
		ProcessStackedSitesNormalized.MIN_VOTES_TO_RETAIN=5;
		ProcessStackedSitesNormalized.FRACTION_TO_RETAIN1=0.75f;
		ProcessStackedSitesNormalized.FRACTION_TO_RETAIN2=0.3f;
		ProcessStackedSitesNormalized.MIN_FRACTION_FROM_READ_ENDS=0.35f;
		ProcessStackedSitesNormalized.SCORE_THRESH=0.034f;
		ProcessStackedSitesNormalized.CENTER_WEIGHT=0.015f;
		ProcessStackedSitesNormalized.sitesIn=0;
		ProcessStackedSitesNormalized.correctIn=0;
		ProcessStackedSitesNormalized.sitesOut=0;
		ProcessStackedSitesNormalized.correctOut=0;
	}

	/** Builds a plus-strand scored site with an inclusive reference span. */
	private static SiteScoreR site(int start, int stop, int readlen){
		return new SiteScoreR(0, Shared.PLUS, start, stop, readlen, start+1, (byte)0,
				100, 0, false, false);
	}

	/** Writes one retained stacked-site input with BBTools text output. */
	private static void write(File file, String... lines){
		if(file.exists()){throw new AssertionError("Refusing to overwrite "+file);}
		TextStreamWriter tsw=new TextStreamWriter(file.getPath(), true, false, true);
		tsw.start();
		for(String line : lines){tsw.println(line);}
		if(tsw.poisonAndWait()){throw new RuntimeException("Error writing "+file);}
	}

	/** Reads an entire retained evidence file with BBTools text input. */
	private static String read(File file){
		StringBuilder sb=new StringBuilder();
		for(String s : TextFile.toStringLines(file.getPath())){sb.append(s).append('\n');}
		return sb.toString();
	}

	/** Returns one named evidence file in the retained output directory. */
	private static File file(String name){return new File(OUTDIR, name);}


	/** Formats one PASS/FAIL token. */
	private static String result(boolean pass){return pass ? "PASS" : "FAIL";}

	/** Retained input/output evidence directory. */
	private static File OUTDIR;

	/** One child JVM result. */
	private static final class ChildResult{
		ChildResult(int exit_, String status_, String stdout_, String stderr_){
			exit=exit_;
			status=status_;
			stdout=stdout_;
			stderr=stderr_;
		}
		final int exit;
		final String status;
		final String stdout;
		final String stderr;
	}

	/** Static-only fixture. */
	private ProcessStackedSitesNormalizedBoundaryTest(){}
}
