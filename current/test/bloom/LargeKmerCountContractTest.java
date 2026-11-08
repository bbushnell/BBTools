package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import bloom.KCountArray2;
import bloom.LargeKmerCount;
import bloom.LargeKmerCount2;
import shared.Shared;
import stream.FASTQ;

/** Finite real-FASTQ contracts for legacy reports, estimator finiteness and per-call state.
 * Uses tiny eight-bit tables and isolated children, not mocks or an exact-cardinality
 * replacement for the two heuristic variants. Preserves documented readsProcessed
 * accumulation. A hidden local NaN cast to a long is not observable through report text.
 * @author Brian Bushnell, Sucrose */
public final class LargeKmerCountContractTest{

	/** Runs all eleven children before reporting failure; child arguments are internal. */
	public static void main(final String[] args) throws Exception{
		if(args.length==2 && args[0].startsWith("child=") && args[1].startsWith("outdir=")){
			child(findCase(args[0].substring(6)), Paths.get(args[1].substring(7)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "A new result directory prevents old child evidence from satisfying this run.");
		Files.createDirectories(root);
		int failures=0;
		for(CaseSpec spec : CASES){failures+=runChild(root, spec);}
		if(failures>0){throw new AssertionError("Legacy large-counter contracts failed "+failures+" of "+CASES.length+" cases.");}
		System.out.println("LARGE_KMER_COUNT_CONTRACT_TEST_PASS cases="+CASES.length);
	}

	/** Bounds only the child this fixture starts; a timeout or nonzero exit is failed evidence. */
	private static int runChild(final Path root, final CaseSpec spec) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(spec.name));
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(LargeKmerCountContractTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		if(!System.getProperty("java.specification.version").startsWith("1.")){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(LargeKmerCountContractTest.class.getName());
		command.add("child="+spec.name);
		command.add("outdir="+dir);
		Files.write(dir.resolve("command.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).redirectInput(Paths.get("/dev/null").toFile())
				.redirectOutput(dir.resolve("stdout.txt").toFile()).redirectError(dir.resolve("stderr.txt").toFile()).start();
		boolean completed=false;
		try{completed=process.waitFor(15, TimeUnit.SECONDS);}
		finally{
			if(process.isAlive()){
				process.destroy();
				if(!process.waitFor(2, TimeUnit.SECONDS)){
					process.destroyForcibly();
					check(process.waitFor(5, TimeUnit.SECONDS), "Fixture-owned child must terminate after forced destruction.");
				}
			}
		}
		final Properties result=new Properties();
		final Path resultFile=dir.resolve("result.properties");
		if(Files.exists(resultFile)){
			try(InputStream input=Files.newInputStream(resultFile)){result.load(input);}
		}
		final int exit=process.exitValue();
		final boolean pass=completed && exit==0 && "true".equals(result.getProperty("pass")) &&
				spec.name.equals(result.getProperty("case")) &&
				Boolean.toString(LargeKmerCountContractTest.class.desiredAssertionStatus()).equals(result.getProperty("assertions"));
		final String summary="case="+spec.name+" completed="+completed+" exit="+exit+" result="+result+" pass="+pass;
		Files.write(dir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Calls real counters, retains their state after failure, and never closes their streams for them. */
	private static void child(final CaseSpec spec, final Path dir) throws Exception{
		Locale.setDefault(Locale.ROOT);
		Shared.setThreads(1);
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=FASTQ.DETECT_QUALITY=false;
		FASTQ.ASCII_OFFSET=33;
		LargeKmerCount.verbose=LargeKmerCount2.verbose=false;
		LargeKmerCount.maxReads=LargeKmerCount2.maxReads=-1;
		check(LargeKmerCount.readsProcessed==0 && LargeKmerCount2.readsProcessed==0,
				"Each child starts with fresh static counters; do not reset them between the retained-state calls.");
		final ArrayList<LiteralInput> inputs=new ArrayList<LiteralInput>();
		final Properties result=new Properties();
		Throwable failure=null;
		try{
			final Path first=writeInput(inputs, dir.resolve("first.fq"), spec.first, "pair/1");
			final Path second=(spec.second==null ? null : writeInput(inputs, dir.resolve("second.fq"), spec.second, "pair/2"));
			if(spec.kind==REPORT){reportCase(spec, dir, first);}
			else if(spec.kind==REPEAT){repeatCase(spec, first, second, result);}
			else{stateCase(spec, dir, inputs, first, result);}
		}catch(Throwable caught){failure=caught;}
		for(LiteralInput input : inputs){
			try{check(Arrays.equals(input.bytes, Files.readAllBytes(input.path)), "Counter changed literal input "+input.path);}
			catch(Throwable changed){if(failure==null){failure=changed;}else{failure.addSuppressed(changed);}}
		}
		recordState(result, "final.", spec.variant);
		result.setProperty("case", spec.name);
		result.setProperty("assertions", Boolean.toString(LargeKmerCountContractTest.class.desiredAssertionStatus()));
		result.setProperty("pass", Boolean.toString(failure==null));
		result.setProperty("failure", failure==null ? "none" : failure.toString());
		try(OutputStream output=Files.newOutputStream(dir.resolve("result.properties"))){result.store(output, "Observed state after real counter calls");}
		if(failure!=null){throw new AssertionError("Child contract failed: "+spec.name, failure);}
		//Return naturally: a leaked non-daemon reader remains visible to the parent deadline.
	}

	/** Checks observable report rows, including an unaffected nonempty denominator control. */
	private static void reportCase(final CaseSpec spec, final Path dir, final Path input) throws Exception{
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream original=System.out;
		try(PrintStream capture=new PrintStream(buffer, true, "UTF-8")){
			System.setOut(capture);
			final String[] flags={input.toString(), "6", "8", "3"};
			if(spec.variant==1){LargeKmerCount.main(flags);}else{LargeKmerCount2.main(flags);}
		}finally{
			System.setOut(original);
			Files.write(dir.resolve("report.txt"), buffer.toByteArray());
		}
		final String report=buffer.toString("UTF-8");
		check(!report.contains("NaN") && !report.contains("Infinity"), "Valid empty/finite input must not print nonfinite statistics.");
		final boolean occupied=spec.first.equals("CCC");
		check(row(report, "Unique:").equals(occupied ? "1" : "0"), "Report occupancy must match the literal zero/one-key input.");
		check(row(report, "Singleton:").equals(occupied ? "100.000% 1" : "0.000% 0"),
				"Singleton report must be finite and use the literal occupied-cell count.");
		check(row(report, "Useful:").equals("0.000% 0"), "These report controls have no cell with depth above one.");
		final long roots=spec.first.length()==0 ? 0 : 1;
		check(readsProcessed(spec.variant)==roots, "One fresh invocation must count input roots, including undefined-only reads.");
		if(spec.variant==2){
			finiteState();
			check(LargeKmerCount2.keysCounted==(occupied ? 1 : 0) && LargeKmerCount2.load==(occupied ? 1 : 0),
					"A single CCC window occupies one cell; empty/undefined inputs contribute none.");
		}
	}

	/** Requires loud rejection at the singular estimator, with a finite neighboring-table control. */
	private static void repeatCase(final CaseSpec spec, final Path first, final Path second, final Properties result){
		KCountArray2 table=null;
		Throwable rejected=null;
		try{table=count(spec.variant, first, second, spec.indexBits);}
		catch(RuntimeException | AssertionError failure){rejected=failure; failure.printStackTrace(System.err);}
		result.setProperty("estimatorException", rejected==null ? "none" : rejected.toString());
		if(spec.indexBits==6){
			check(rejected!=null && explainsEstimatorFailure(rejected),
					"A zero/undefined collision-estimator denominator must throw an explanatory exception, not return NaN or invent a fallback: "+rejected);
			return;
		}
		check(rejected==null && table!=null, "The neighboring cells32 control must count normally: "+rejected);
		checkTable(table, 1L<<spec.indexBits, 21, 3);
		check(LargeKmerCount2.keysCounted==3 && LargeKmerCount2.load==1 && LargeKmerCount2.actualCollisions==2,
				"Three identical CCC windows imply three increments, one occupied cell and two old-value hits.");
		check(LargeKmerCount2.readsProcessed==1, "The mate-only repeat case still contains one input root.");
		check(LargeKmerCount2.invCells==1d/32, "The smaller-table control must use the existing nonzero inverse-cell probability.");
		finiteState();
		check(LargeKmerCount2.expectedCollisions>0, "The smaller-table control must exercise a positive finite denominator.");
	}

	/** Rejects unrelated reader/setup failures; the diagnostic must identify the invalid estimator. */
	private static boolean explainsEstimatorFailure(final Throwable failure){
		for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
			if(cause.getMessage()==null){continue;}
			final String message=cause.getMessage().toLowerCase(Locale.ROOT);
			final boolean context=message.contains("estimat") || message.contains("collision") || message.contains("denominator");
			final boolean reason=message.contains("zero") || message.contains("undefined") || message.contains("non-finite") ||
					message.contains("nan") || message.contains("infin") || message.contains("singular") || message.contains("positive");
			if(context && reason){return true;}
		}
		return false;
	}

	/** Separates new-table statistics from intentionally accumulated input-root reporting. */
	private static void stateCase(final CaseSpec spec, final Path dir, final ArrayList<LiteralInput> inputs,
			final Path first, final Properties result) throws IOException{
		final KCountArray2 before=count(spec.variant, first, null, spec.indexBits);
		checkTable(before, 1L<<spec.indexBits, spec.variant==1 ? 8 : 21, 3);
		if(spec.variant==1){
			check(LargeKmerCount.collisionsA==1 && LargeKmerCount.collisionsB==2,
					"At sixteen cells, literal CCC maps both probes to cell8: first hit A, then two B hits.");
		}else{
			finiteState();
			check(LargeKmerCount2.keysCounted==3 && LargeKmerCount2.actualCollisions==2 && LargeKmerCount2.expectedCollisions>0,
					"The first call must populate every retained-state counter under test.");
		}
		recordState(result, "first.", spec.variant);
		final Path empty=writeInput(inputs, dir.resolve("empty.fq"), "", "empty");
		final int nextBits=spec.variant==1 ? 4 : 6;
		final KCountArray2 after=count(spec.variant, empty, null, nextBits);
		check(before!=after, "Each countFastq invocation promises a fresh table.");
		checkTable(after, 1L<<nextBits, -1, 0);
		check(readsProcessed(spec.variant)==1, "Preserve the documented accumulated root count rather than changing its policy here.");
		if(spec.variant==1){
			check(LargeKmerCount.collisionsA==0 && LargeKmerCount.collisionsB==0, "The first variant already resets both per-call collision counters.");
		}else{
			finiteState();
			check(LargeKmerCount2.keysCounted==0 && LargeKmerCount2.actualCollisions==0 && LargeKmerCount2.expectedCollisions==0 &&
					LargeKmerCount2.load==0 && LargeKmerCount2.probNewKeyCollisions==0,
					"An empty fresh table must not inherit key/collision estimates from the prior invocation.");
		}
	}

	/** Uses real file-owning public APIs; no table or stream implementation is replaced. */
	private static KCountArray2 count(final int variant, final Path first, final Path second, final int indexBits){
		final String mate=second==null ? null : second.toString();
		return variant==1 ? LargeKmerCount.countFastq(first.toString(), mate, indexBits, 8, 3) :
				LargeKmerCount2.countFastq(first.toString(), mate, indexBits, 8, 3);
	}

	/** Reads every allocated cell, avoiding the separately documented cached-occupancy field. */
	private static void checkTable(final KCountArray2 table, final long cells, final int occupied, final int depth){
		check(table.cells==cells && table.cellBits==8 && table.maxValue==255, "Fixture requires the literal tiny eight-bit table geometry.");
		for(long cell=0; cell<cells; cell++){
			final int expected=cell==occupied ? depth : 0;
			check(table.read(cell)==expected, "Literal table cell "+cell+" must have depth "+expected+", observed "+table.read(cell));
		}
	}

	/** Only finiteness is asserted for heuristic doubles, not exact distinct-key estimation. */
	private static void finiteState(){
		for(double value : new double[]{LargeKmerCount2.expectedCollisions, LargeKmerCount2.probNewKeyCollisions,
				LargeKmerCount2.invCells, LargeKmerCount2.invKmerSpace}){
			check(!Double.isNaN(value) && !Double.isInfinite(value), "Public estimator state must remain finite for valid literal inputs: "+value);
		}
	}

	/** Saves raw shared state even when a finite-contract check fails. */
	private static void recordState(final Properties result, final String prefix, final int variant){
		result.setProperty(prefix+"readsProcessed", Long.toString(readsProcessed(variant)));
		if(variant==1){
			result.setProperty(prefix+"collisionsA", Long.toString(LargeKmerCount.collisionsA));
			result.setProperty(prefix+"collisionsB", Long.toString(LargeKmerCount.collisionsB));
		}else{
			result.setProperty(prefix+"keysCounted", Long.toString(LargeKmerCount2.keysCounted));
			result.setProperty(prefix+"load", Long.toString(LargeKmerCount2.load));
			result.setProperty(prefix+"actualCollisions", Long.toString(LargeKmerCount2.actualCollisions));
			result.setProperty(prefix+"expectedCollisions", Double.toString(LargeKmerCount2.expectedCollisions));
			result.setProperty(prefix+"probNewKeyCollisions", Double.toString(LargeKmerCount2.probNewKeyCollisions));
			result.setProperty(prefix+"invCells", Double.toString(LargeKmerCount2.invCells));
			result.setProperty(prefix+"invKmerSpace", Double.toString(LargeKmerCount2.invKmerSpace));
		}
	}

	/** Returns the documented accumulated root count for the selected implementation. */
	private static long readsProcessed(final int variant){
		return variant==1 ? LargeKmerCount.readsProcessed : LargeKmerCount2.readsProcessed;
	}

	/** Extracts exactly one labeled report row, ignoring spacing but no numeric content. */
	private static String row(final String report, final String label){
		String value=null;
		for(String line : report.split("\\r?\\n")){
			if(line.startsWith(label)){
				check(value==null, "Duplicate report row: "+label);
				value=line.substring(label.length()).trim().replaceAll("\\s+", " ");
			}
		}
		check(value!=null, "Missing report row: "+label);
		return value;
	}

	/** Writes either a genuine empty file or one ASCII-33 FASTQ record, retaining its exact bytes. */
	private static Path writeInput(final ArrayList<LiteralInput> inputs, final Path path,
			final String sequence, final String id) throws IOException{
		final char[] qualities=new char[sequence.length()];
		Arrays.fill(qualities, 'I');
		final byte[] bytes=sequence.length()==0 ? new byte[0] :
				("@"+id+"\n"+sequence+"\n+\n"+new String(qualities)+"\n").getBytes(StandardCharsets.US_ASCII);
		Files.write(path, bytes);
		inputs.add(new LiteralInput(path, bytes));
		return path;
	}

	/** Rejects unknown child selectors rather than silently running a different case. */
	private static CaseSpec findCase(final String name){
		for(CaseSpec spec : CASES){if(spec.name.equals(name)){return spec;}}
		throw new IllegalArgumentException("Unknown child case: "+name);
	}

	/** Keeps independent fixture checks active under both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Literal input and its expected bytes; not a copy of any user dataset. */
	private static final class LiteralInput{
		LiteralInput(final Path path_, final byte[] bytes_){path=path_; bytes=bytes_;}
		final Path path;
		final byte[] bytes;
	}

	/** Immutable case specification; all expected key counts remain explicit above. */
	private static final class CaseSpec{
		CaseSpec(final String name_, final int variant_, final int kind_, final String first_, final String second_, final int indexBits_){
			name=name_; variant=variant_; kind=kind_; first=first_; second=second_; indexBits=indexBits_;
		}
		final String name, first, second;
		final int variant, kind, indexBits;
	}

	private static final int REPORT=0, REPEAT=1, STATE=2;
	private static final CaseSpec[] CASES={
		new CaseSpec("large1-empty", 1, REPORT, "", null, 6),
		new CaseSpec("large1-undefined", 1, REPORT, "NNNNN", null, 6),
		new CaseSpec("large1-single", 1, REPORT, "CCC", null, 6),
		new CaseSpec("large2-empty", 2, REPORT, "", null, 6),
		new CaseSpec("large2-undefined", 2, REPORT, "NNNNN", null, 6),
		new CaseSpec("large2-single", 2, REPORT, "CCC", null, 6),
		new CaseSpec("large2-repeat", 2, REPEAT, "CCCCC", null, 6),
		new CaseSpec("large2-mate-repeat", 2, REPEAT, "CC", "CCCCC", 6),
		new CaseSpec("large2-smaller-table", 2, REPEAT, "CCCCC", null, 5),
		new CaseSpec("large2-retained-state", 2, STATE, "CCCCC", null, 5),
		new CaseSpec("large1-reset-control", 1, STATE, "CCCCC", null, 4)
	};
}
