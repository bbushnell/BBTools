package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import dna.FastaToChromArrays2;
import jgi.AssemblyStats2;
import parse.Parse;
import shared.Shared;
import structures.LongList;

/**
 * Small direct AssemblyStats2 scanner cases for the record-length overflow guard.
 * FASTA and FASTQ literals cross the small-histogram/large-list cutoff without
 * huge inputs. Read-failure rows capture and rescue only the optional GC writer
 * whose fname matches the case's unique path. Guard-unit rows are reflection-only
 * evidence for the private increment helper and never claim a full parser reached
 * Integer.MAX_VALUE. Separate parser modes retain correct literal expectations
 * beside baseline observations; capture success never means a buggy parser passed.
 * Refill-guards mode calls only private helpers with tiny buffers and near-MAX
 * carried lengths; it never claims full-parser MAX coverage.
 * GC-header names-only modes retain accounting observations without requiring
 * the separate JGI117 header-accounting repair.
 * @author Jean, Sucrose
 */
public final class AssemblyStats2RecordLengthTest{

	/** Requires phase=baseline|fixed and a fresh output directory; mode defaults to record-length. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, mode="record-length";
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else{throw new IllegalArgumentException("Unknown record-length fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		require("record-length".equals(mode) || "parser-boundaries".equals(mode) ||
			"gc-headers".equals(mode) || "gc-header-controls".equals(mode) ||
			"gc-header-names".equals(mode) || "gc-header-name-controls".equals(mode) ||
			"advisory-memory".equals(mode) || "refill-guards".equals(mode), "Unknown fixture mode: "+mode);
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats2recordlengthtest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		final boolean gcNamesOnly="gc-header-names".equals(mode) || "gc-header-name-controls".equals(mode);
		require(!"gc-header-controls".equals(mode) || "fixed".equals(phase), "GC-header controls require phase=fixed.");
		require(!gcNamesOnly || "fixed".equals(phase), "GC-header names-only modes require phase=fixed.");
		require(!"refill-guards".equals(mode) || "fixed".equals(phase), "Refill-guard helper checks require phase=fixed.");
		Files.createDirectory(root);
		if("refill-guards".equals(mode)){
			new RefillGuardRows(root).run();
			return;
		}
		if("advisory-memory".equals(mode)){
			new AdvisoryMemoryRows(root, phase).run();
			return;
		}
		if(gcNamesOnly){
			final GcHeaderRows rows=new GcHeaderRows(root, phase, false);
			if("gc-header-name-controls".equals(mode)){rows.runControls();}
			else{rows.run();}
			return;
		}
		if("gc-header-controls".equals(mode)){
			new GcHeaderRows(root, phase).runControls();
			return;
		}
		if("gc-headers".equals(mode)){
			new GcHeaderRows(root, phase).run();
			return;
		}
		if("parser-boundaries".equals(mode)){
			new ParserRows(root, phase).run();
			return;
		}
		final AssemblyStats2RecordLengthTest test=new AssemblyStats2RecordLengthTest(root, phase.equals("baseline"));
		test.runScanner("fasta-storage-no-gc", true, null, false);
		test.runScanner("fastq-storage-no-gc", false, null, false);
		test.runScanner("fasta-storage-file-gc", true, root.resolve("fasta-storage.gc").toString(), false);
		test.runScanner("fastq-storage-file-gc", false, root.resolve("fastq-storage.gc").toString(), false);
		test.runScanner("fasta-stdout-gc", true, "stdout", true);
		test.runFailure("fasta-read-failure-file-gc", true, root.resolve("fasta-failure.gc").toString());
		test.runFailure("fastq-read-failure-file-gc", false, root.resolve("fastq-failure.gc").toString());
		test.runGuard(new GuardRow("guard-zero", 0, "1"));
		test.runGuard(new GuardRow("guard-max-minus-one", Integer.MAX_VALUE-1, Integer.toString(Integer.MAX_VALUE)));
		test.runGuard(new GuardRow("guard-max", Integer.MAX_VALUE, "reject"));
		test.check(test.scannerCases==5, "Expected five normal scanner cases: "+test.scannerCases);
		test.check(test.failureCases==2, "Expected two read-failure cases: "+test.failureCases);
		test.check(test.guardCases==3, "Expected three guard-unit cases: "+test.guardCases);
		test.check(test.guardUnavailable==(test.baseline ? 3 : 0),
			"Unexpected unavailable guard rows: "+test.guardUnavailable);
		test.check(test.writerLeaks==(test.baseline ? 2 : 0),
			"Unexpected live owned GC writers: "+test.writerLeaks);
		test.check(test.writerRescues==(test.baseline ? 2 : 0),
			"Unexpected fixture writer rescues: "+test.writerRescues);
		final String summary="ASSEMBLY_STATS2_RECORD_LENGTH_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" scannerCases="+test.scannerCases+
			" failureCases="+test.failureCases+" guardCases="+test.guardCases+
			" guardUnavailable="+test.guardUnavailable+" writerLeaks="+test.writerLeaks+
			" writerRescues="+test.writerRescues+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Binds private state observations and the optional future guard helper. */
	private AssemblyStats2RecordLengthTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_;
		baseline=baseline_;
		clist=field("clist"); slist=field("slist"); sclist1=field("sclist1");
		sclist2=field("sclist2"); llist=field("llist"); tlist=field("tlist");
		guard=guardMethod();
		Files.write(root.resolve("cases.tsv"), ("case\tphase\tformat\tgcTarget\texpectedOutcome"+
			"\tactualOutcome\tcounts\tclist2\tclist3\tslist3\tsclist1_3\tsclist2_3"+
			"\tlargeContigs\tlargeTriples\tgcLines\tstdoutBytes\tstdoutMarker\tstdoutArtifact"+
			"\tpreRescueLiveWriters\tpostRescueLiveWriters\trescuedWriters\trescueErrorFlags"+
			"\trescueFailure\tinputBulkReads\tinputSingleReads\tinputDelivered\tinputErrors"+
			"\tprimaryClass\tprimaryMessage\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("guard.tsv"), ("case\tphase\tcoverage\tinput\texpected"+
			"\tactualOutcome\tactualValue\texceptionClass\texceptionMessage\tfullParserMaxReached\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Runs one normal finite FASTA or FASTQ scanner and checks exact storage branches. */
	private void runScanner(final String name, final boolean fasta, final String gc,
			final boolean captureStdout) throws Exception{
		scannerCases++;
		final PrintStream original=System.out;
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream();
		StaticState state=null;
		PrintStream capture=null;
		AssemblyStats2 stats=null;
		long[] counts=null;
		Throwable error=null;
		ArrayList<TextStreamWriter> preRescue=new ArrayList<TextStreamWriter>();
		Rescue rescue=Rescue.NONE;
		int postRescue=0;
		String gcText="";
		String stdoutText="";
		try{
			state=new StaticState();
			try{
				stats=newStats();
				final InputStream input=new ByteArrayInputStream(fasta ? FASTA : FASTQ);
				if(captureStdout){
					capture=new PrintStream(stdout, true, "UTF-8");
					System.setOut(capture);
				}
				counts=fasta ? stats.countFasta(input, gc) : stats.countFastq(input, gc);
				if(captureStdout){
					capture.print(STDOUT_REUSE_MARKER);
					capture.flush();
				}
			}catch(final Throwable caught){
				error=caught;
			}finally{
				if(gc!=null){
					preRescue=matchingWriters(gc);
					rescue=poisonAll(preRescue);
					postRescue=matchingWriters(gc).size();
				}
				stdoutText=captureStdout ? stdout.toString("UTF-8") : "";
				gcText=captureStdout ? gcBeforeMarker(stdoutText) : readGc(gc);
			}
		}finally{
			System.setOut(original);
			if(state!=null){state.restore();}
			if(capture!=null){capture.close();}
		}
		writerLeaks+=preRescue.size();
		writerRescues+=rescue.rescued;
		final boolean stdoutBlankLine=!captureStdout || gcText.endsWith("\n\n");
		final boolean stdoutMarker=captureStdout &&
			stdoutText.substring(gcText.length()).equals(STDOUT_REUSE_MARKER);
		final String stdoutArtifact=captureStdout ? name+".stdout.txt" : "NA";
		if(captureStdout){
			Files.write(root.resolve(stdoutArtifact), stdoutText.getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE_NEW);
		}
		final String expectedGc=gc==null ? "" : expectedGc(fasta)+(captureStdout ? "\n" : "");
		check(error==null, name+": normal scanner threw "+error);
		check(Arrays.equals(counts, EXPECTED_COUNTS), name+": counts expected "+
			Arrays.toString(EXPECTED_COUNTS)+", saw "+Arrays.toString(counts));
		check(stats!=null && storageMatches(stats), name+": cutoff branch storage differed");
		check(gc==null || expectedGc.equals(gcText), name+": GC text differed: "+gcText);
		check(stdoutBlankLine, name+": stdout GC spacer was not printed after rows");
		check(!captureStdout || stdoutMarker, name+": borrowed stdout was not reusable after GC writer finalization");
		check(preRescue.isEmpty(), name+": normal scanner needed fixture writer rescue");
		check(postRescue==0, name+": completed normal scanner left a live owned writer");
		recordScanner(name, fasta, gc, "return", error==null ? "return" : "throw", counts, stats,
			gcText, stdoutText.getBytes(StandardCharsets.UTF_8).length,
			captureStdout ? Boolean.toString(stdoutMarker) : "NA", stdoutArtifact,
			preRescue.size(), postRescue, rescue, 0, 0, 0, 0, error);
	}

	/** Runs one owner-side read failure and rescues only the optional writer for the unique GC path. */
	private void runFailure(final String name, final boolean fasta, final String gc) throws Exception{
		failureCases++;
		final IOException primary=new IOException("Injected AssemblyStats2 "+name);
		StaticState state=null;
		AssemblyStats2 stats=null;
		long[] counts=null;
		Throwable error=null;
		ThrowingInputStream input=null;
		ArrayList<TextStreamWriter> preRescue=new ArrayList<TextStreamWriter>();
		Rescue rescue=Rescue.NONE;
		int postRescue=0;
		String gcText="";
		try{
			state=new StaticState();
			try{
				stats=newStats();
				final byte[] prefix=fasta ? FASTA_PREFIX : FASTQ_PREFIX;
				input=new ThrowingInputStream(prefix, primary);
				counts=fasta ? stats.countFasta(input, gc) : stats.countFastq(input, gc);
			}catch(final Throwable caught){
				error=caught;
			}finally{
				preRescue=matchingWriters(gc);
				rescue=poisonAll(preRescue);
				postRescue=matchingWriters(gc).size();
				gcText=readGc(gc);
			}
		}finally{
			if(state!=null){state.restore();}
		}
		writerLeaks+=preRescue.size();
		writerRescues+=rescue.rescued;
		check(error==primary, name+": expected retained IOException primary, saw "+error);
		check(counts==null, name+": read failure returned counts "+Arrays.toString(counts));
		check(preRescue.size()==(baseline ? 1 : 0),
			name+": expected "+(baseline ? 1 : 0)+" live owned writers, saw "+preRescue.size());
		check(postRescue==0, name+": rescue left live owned writers");
		check(rescue.failure==null, name+": fixture writer rescue failed "+rescue.failure);
		check(rescue.errorFlags==0, name+": fixture writer rescue saw errorState flags "+rescue.errorFlags);
		check(expectedFailureGc(fasta).equals(gcText), name+": post-cleanup GC text differed: "+gcText);
		check(input!=null && input.bulkReads==2 && input.singleReads==0 &&
			input.delivered==(fasta ? FASTA_PREFIX.length : FASTQ_PREFIX.length) && input.errors==1,
			name+": unexpected bounded input metrics");
		recordScanner(name, fasta, gc, "throw", error==null ? "return" : "throw", counts, stats,
			gcText, 0, "NA", "NA", preRescue.size(), postRescue, rescue,
			input==null ? 0 : input.bulkReads, input==null ? 0 : input.singleReads,
			input==null ? 0 : input.delivered, input==null ? 0 : input.errors, error);
	}

	/** Checks the future tiny helper directly; baseline source is expected not to have it. */
	private void runGuard(final GuardRow row) throws Exception{
		guardCases++;
		Object result=null;
		Throwable error=null;
		if(guard==null){
			guardUnavailable++;
		}else{
			try{result=guard.invoke(null, row.input);}
			catch(final InvocationTargetException caught){error=caught.getCause();}
		}
		final String outcome=guard==null ? "unavailable" : error==null ? "return" : "throw";
		if(baseline){
			check(guard==null, row.name+": baseline must not fabricate a guard helper");
		}else if("reject".equals(row.expected)){
			check(error instanceof IllegalArgumentException && GUARD_MESSAGE.equals(error.getMessage()),
				row.name+": expected exact MAX_VALUE IllegalArgumentException, saw "+error);
		}else{
			check(error==null && result instanceof Integer && row.expected.equals(result.toString()),
				row.name+": expected "+row.expected+", saw "+result+" / "+error);
		}
		final String observed=row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\tprivate-incrementRecordLength-reflection\t"+row.input+"\t"+row.expected+
			"\t"+outcome+"\t"+(result==null ? "unavailable" : result.toString())+
			"\t"+type(error)+"\t"+clean(error==null ? "" : String.valueOf(error.getMessage()))+
			"\tfalse\n";
		Files.write(root.resolve("guard.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Writes one scanner/failure observation row. */
	private void recordScanner(final String name, final boolean fasta, final String gc, final String expected,
			final String outcome, final long[] counts, final AssemblyStats2 stats, final String gcText,
			final int stdoutBytes, final String stdoutMarker, final String stdoutArtifact,
			final int preRescue, final int postRescue, final Rescue rescue, final int bulkReads,
			final int singleReads, final int delivered, final int inputErrors,
			final Throwable primary) throws Exception{
		final String observed=name+"\t"+(baseline ? "baseline" : "fixed")+
				"\t"+(fasta ? "fasta" : "fastq")+"\t"+(gc==null ? "none" : clean(gc))+
			"\t"+expected+"\t"+outcome+"\t"+Arrays.toString(counts)+
			"\t"+(stats==null ? "unavailable" : Long.toString(value(stats, clist, 2)))+
			"\t"+(stats==null ? "unavailable" : Long.toString(value(stats, clist, 3)))+
			"\t"+(stats==null ? "unavailable" : Long.toString(value(stats, slist, 3)))+
			"\t"+(stats==null ? "unavailable" : Long.toString(value(stats, sclist1, 3)))+
			"\t"+(stats==null ? "unavailable" : Long.toString(value(stats, sclist2, 3)))+
			"\t"+(stats==null ? "unavailable" : Arrays.toString(array(stats, llist)))+
			"\t"+(stats==null ? "unavailable" : triples(stats))+
			"\t"+gcRows(gcText)+"\t"+stdoutBytes+"\t"+stdoutMarker+"\t"+stdoutArtifact+
			"\t"+preRescue+"\t"+postRescue+
			"\t"+rescue.rescued+"\t"+rescue.errorFlags+
			"\t"+(rescue.failure==null ? "none" : type(rescue.failure)+":"+clean(String.valueOf(rescue.failure.getMessage())))+
			"\t"+bulkReads+"\t"+singleReads+"\t"+delivered+"\t"+inputErrors+
			"\t"+type(primary)+"\t"+clean(primary==null ? "" : String.valueOf(primary.getMessage()))+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Creates a scanner instance with the literal cutoff and N-split options. */
	private static AssemblyStats2 newStats(){
		return new AssemblyStats2(new String[]{"in=stdin", "cutoff=4", "n=2",
			"format=1", "gcformat=1", "gcfraction=f", "showbbmap=f",
			"showspeed=f", "n_=f"});
	}

	/** Checks the exact LongList and Triple storage that proves both cutoff branches. */
	private boolean storageMatches(final AssemblyStats2 stats) throws Exception{
		return value(stats, clist, 0)==0 && value(stats, clist, 1)==0 &&
			value(stats, clist, 2)==1 && value(stats, clist, 3)==1 &&
			value(stats, clist, 4)==0 && value(stats, slist, 3)==1 &&
			value(stats, slist, 4)==0 && value(stats, sclist1, 3)==1 &&
			value(stats, sclist2, 3)==3 && Arrays.equals(array(stats, llist), new long[]{4, 4}) &&
			"4/1/4,8/2/6".equals(triples(stats));
	}

	/** Returns a LongList value from a private AssemblyStats2 list field. */
	private static long value(final AssemblyStats2 stats, final Field field, final int loc) throws Exception{
		return ((LongList)field.get(stats)).get(loc);
	}

	/** Returns the logical contents of a private AssemblyStats2 LongList field. */
	private static long[] array(final AssemblyStats2 stats, final Field field) throws Exception{
		final LongList list=(LongList)field.get(stats);
		return Arrays.copyOf(list.array, list.size);
	}

	/** Returns compact Triple length/contig-count/contig-length observations. */
	private static String triples(final AssemblyStats2 stats) throws Exception{
		final StringBuilder sb=new StringBuilder();
		final ArrayList<?> list=(ArrayList<?>)tlist.get(stats);
		for(int i=0; i<list.size(); i++){
			if(i>0){sb.append(',');}
			final Object triple=list.get(i);
			final Class<?> clazz=triple.getClass();
			final Field length=clazz.getDeclaredField("length");
			final Field contigs=clazz.getDeclaredField("contigs");
			final Field contiglen=clazz.getDeclaredField("contiglen");
			length.setAccessible(true);
			contigs.setAccessible(true);
			contiglen.setAccessible(true);
			sb.append(length.getLong(triple)).append('/').append(contigs.getLong(triple)).
				append('/').append(contiglen.getLong(triple));
		}
		return sb.toString();
	}

	/** Reads a GC file or returns an empty string for absent/non-file targets. */
	private static String readGc(final String gc) throws IOException{
		if(gc==null || "stdout".equalsIgnoreCase(gc)){return "";}
		final Path path=Paths.get(gc);
		return Files.exists(path) ? new String(Files.readAllBytes(path), StandardCharsets.UTF_8) : "";
	}

	/** Returns only the AssemblyStats2 stdout GC output before the fixture reuse marker. */
	private static String gcBeforeMarker(final String stdoutText){
		final int marker=stdoutText.indexOf(STDOUT_REUSE_MARKER);
		return marker<0 ? stdoutText : stdoutText.substring(0, marker);
	}

	/** Counts data rows in a GC output string. */
	private static int gcRows(final String gcText){
		if(gcText==null || gcText.isEmpty()){return 0;}
		int rows=0;
		final String[] lines=gcText.split("\n");
		for(final String line : lines){
			if(line.length()>0 && line.charAt(0)!='#'){rows++;}
		}
		return rows;
	}

	/** Returns the exact GC output for the three literal scanner records. */
	private static String expectedGc(final boolean fasta){
		final String prefix=fasta ? "" : "@";
		return "#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC\n"+
			prefix+"small\t3\t1.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\n"+
			prefix+"edge\t4\t1.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\n"+
			prefix+"split\t8\t0.3333\t0.0000\t0.6667\t0.0000\t0.2500\t0.0000\t0.0000\t0.6667\n";
	}

	/** Returns the exact optional GC bytes that remain after an owner-side read failure. */
	private static String expectedFailureGc(final boolean fasta){
		return "#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC\n"+
			(fasta ? "" : "@ok\t2\t0.5000\t0.5000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.5000\n");
	}

	/** Resolves a private instance field without widening production visibility. */
	private static Field field(final String name) throws Exception{
		final Field field=AssemblyStats2.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Resolves the future fixed-only increment helper if this source has it. */
	private static Method guardMethod() throws Exception{
		try{
			final Method method=AssemblyStats2.class.getDeclaredMethod("incrementRecordLength", Integer.TYPE);
			method.setAccessible(true);
			require(Modifier.isStatic(method.getModifiers()) && method.getReturnType()==Integer.TYPE,
				"The record-length helper must be static int incrementRecordLength(int).");
			return method;
		}catch(final NoSuchMethodException nsme){
			return null;
		}
	}

	/** Returns only still-live TextStreamWriter threads for a case's unique GC output. */
	private static ArrayList<TextStreamWriter> matchingWriters(final String output){
		final ArrayList<TextStreamWriter> writers=new ArrayList<TextStreamWriter>();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof TextStreamWriter){
				final TextStreamWriter writer=(TextStreamWriter)thread;
				if(thread.isAlive() && output.equals(writer.fname)){writers.add(writer);}
			}
		}
		return writers;
	}

	/** Attempts fixture rescue for every captured owned writer. */
	private static Rescue poisonAll(final ArrayList<TextStreamWriter> writers){
		Throwable failure=null;
		int errorFlags=0;
		for(final TextStreamWriter writer : writers){
			try{
				if(writer.poisonAndWait()){errorFlags++;}
			}catch(final Throwable caught){
				if(failure==null){failure=caught;}
				else if(caught!=failure){failure.addSuppressed(caught);}
			}
		}
		return new Rescue(writers.size(), errorFlags, failure);
	}

	/** Makes one exception message fit in a TSV cell. */
	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}

	/** Records absent exceptions explicitly. */
	private static String type(final Throwable error){
		return error==null ? "none" : error.getClass().getName();
	}

	/** Accumulates mismatches so all finite rows leave evidence. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Saves and restores every mutable static AssemblyStats2 field used by direct scanner rows. */
	private static final class StaticState{
		StaticState() throws Exception{
			for(final StaticRef ref : STATIC_FIELDS){
				final Field field=ref.clazz.getDeclaredField(ref.name);
				field.setAccessible(true);
				values.put(field, field.get(null));
			}
		}
		void restore() throws Exception{
			for(final Map.Entry<Field, Object> entry : values.entrySet()){
				entry.getKey().set(null, entry.getValue());
			}
		}
		private final LinkedHashMap<Field, Object> values=new LinkedHashMap<Field, Object>();
	}

	/** A mutable static field whose value must be restored after direct constructor calls. */
	private static final class StaticRef{
		StaticRef(final Class<?> clazz_, final String name_){
			clazz=clazz_;
			name=name_;
		}
		final Class<?> clazz;
		final String name;
	}

	/** Finite JGI115/JGI114 rows, isolated from the accepted record-length/writer path. */
	private static final class ParserRows{
		ParserRows(final Path root_, final String phase_) throws Exception{
			root=root_; phase=phase_;
			tlist=field("tlist");
			Files.write(root.resolve("parser-cases.tsv"), ("case\tphase\tformat\tchunkCap\tinputBytes"+
				"\tcorrectCounts\tcorrectStorage\tbaselineCounts\tbaselineStorage\tactualCounts\tactualStorage"+
				"\texceptionClass\texceptionMessage\tbulkReads\tsingleReads\tdelivered\tparserCorrect\texpectedObservationMatched\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		}

		/** Retains the original eighteen rows, adds bounded repair controls, and prepares the same three CLI inputs. */
		void run() throws Exception{
			for(final int cap : new int[]{Integer.MAX_VALUE, 1}){
				for(final String ending : new String[]{"lf", "cr", "crlf"}){
					final boolean defect="crlf".equals(ending);
					row("fastq-"+ending+"-"+cap, false, fastq(ending), cap,
						FASTQ_COUNTS, FASTQ_STORAGE, defect ? ZERO_COUNTS : FASTQ_COUNTS,
						defect ? EMPTY_STORAGE : FASTQ_STORAGE);
				}
				for(final int length : new int[]{3, 4}){
					for(final String tail : new String[]{"none", "terminated", "partial"}){
						final String prefix=length==3 ? ">a\nAAA\n" : ">a\nAAAA\n";
						final String text=prefix+("none".equals(tail) ? "" : "terminated".equals(tail) ? ">tail\n" : ">tail");
						final long[] counts=length==3 ? THREE_A : FOUR_A;
						final String correct=length==3 ? THREE_STORAGE : FOUR_STORAGE;
						final String observed="partial".equals(tail) ?
							(length==3 ? THREE_DUPLICATE : FOUR_DUPLICATE) : correct;
						row("fasta-"+length+"-"+tail+"-"+cap, true,
							text.getBytes(StandardCharsets.US_ASCII), cap, counts, correct, counts, observed);
					}
				}
			}
			for(final int cap : new int[]{Integer.MAX_VALUE, 1}){
				for(final String ending : new String[]{"lf", "cr", "crlf"}){
					final boolean defect="crlf".equals(ending);
					row("fastq-empty-"+ending+"-"+cap, false, emptyFastq(ending), cap,
						EMPTY_RECORD_COUNTS, TWO_STORAGE, defect ? ZERO_COUNTS : EMPTY_RECORD_COUNTS,
						defect ? EMPTY_STORAGE : TWO_STORAGE);
				}
			}
			for(final int cap : new int[]{3, 16}){
				// CR ends the header read at cap3 and the first quality read at cap16.
				row("fastq-crlf-"+cap, false, fastq("crlf"), cap,
					FASTQ_COUNTS, FASTQ_STORAGE, ZERO_COUNTS, EMPTY_STORAGE);
			}
			// Header CRLF ends at byte8; the following LF is an empty logical sequence line.
			row("fastq-empty-mixed-crlf-header-lf-body-8", false,
				"@empty\r\n\n+\n\n@r\nAC\n+\nII\n".getBytes(StandardCharsets.US_ASCII), 8,
				EMPTY_RECORD_COUNTS, TWO_STORAGE, MIXED_BASELINE_COUNTS, TWO_STORAGE);
			Files.write(root.resolve("cli-expected.tsv"), ("case\tinput\texpectedExit\tscaffolds\tcontigs\tscaf_bp\tcontig_bp\tscaf_max\tctg_max\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			for(final String ending : new String[]{"lf", "cr", "crlf"}){
				final String name="fastq-"+ending;
				Files.write(root.resolve(name+".fq"), fastq(ending), StandardOpenOption.CREATE_NEW);
				Files.write(root.resolve("cli-expected.tsv"), (name+"\t"+name+".fq\t0\t2\t2\t7\t7\t4\t4\n")
					.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
			}
			finish("direct", 27);
		}

		/** Records exact correct and source-predicted observations without deriving either from output. */
		private void row(final String name, final boolean fasta, final byte[] inputBytes, final int cap,
				final long[] correctCounts, final String correctStorage,
				final long[] baselineCounts, final String baselineStorage) throws Exception{
			assert(correctCounts.length==8 && baselineCounts.length==8) : "AS2 scanner contract has eight count categories (makeCharToNum).";
			final CappedInputStream input=new CappedInputStream(inputBytes, cap);
			final StaticState saved=new StaticState();
			AssemblyStats2 stats=null;
			long[] actual=null;
			String storage="unavailable";
			Throwable failure=null;
			try{
				stats=new AssemblyStats2(new String[]{"in=stdin", "cutoff=4", "minscaf=1", "n=2",
					"format=1", "gcformat=1", "gcfraction=f", "showbbmap=f", "showspeed=f", "n_=f"});
				try{actual=fasta ? stats.countFasta(input, null) : stats.countFastq(input, null);}
				catch(final Throwable caught){failure=caught;}
				storage=storage(stats);
			}finally{saved.restore();}
			final boolean correct=failure==null && Arrays.equals(correctCounts, actual) && correctStorage.equals(storage);
			final boolean match=failure==null && Arrays.equals("baseline".equals(phase) ? baselineCounts : correctCounts, actual) &&
				("baseline".equals(phase) ? baselineStorage : correctStorage).equals(storage) &&
				input.singleReads==0 && input.delivered==inputBytes.length;
			rows++; if(correct){correctRows++;} if(!match){unexpected++;}
			Files.write(root.resolve(name+".input"), inputBytes, StandardOpenOption.CREATE_NEW);
			final String text=name+"\t"+phase+"\t"+(fasta ? "fasta" : "fastq")+"\t"+cap+"\t"+inputBytes.length+
				"\t"+Arrays.toString(correctCounts)+"\t"+correctStorage+"\t"+Arrays.toString(baselineCounts)+"\t"+baselineStorage+
				"\t"+Arrays.toString(actual)+"\t"+storage+"\t"+type(failure)+"\t"+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+
				"\t"+input.bulkReads+"\t"+input.singleReads+"\t"+input.delivered+"\t"+correct+"\t"+match+"\n";
			Files.write(root.resolve("parser-cases.tsv"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		}

		/** A capture pass asserts its observation contract; parserCorrectRows reports scientific correctness separately. */
		private void finish(final String group, final int expectedRows) throws IOException{
			final boolean match=rows==expectedRows && unexpected==0;
			final String summary="ASSEMBLY_STATS2_PARSER_CAPTURE_"+(match ? "PASS" : "FAIL")+" group="+group+
				" phase="+phase+" rows="+rows+" parserCorrectRows="+correctRows+" unexpectedObservations="+unexpected+"\n";
			Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			System.out.print(summary);
			if(!match){throw new AssertionError(summary+"Inspect saved correct/observed columns in "+root);}
		}

		/** Reads all nonzero bins and the entire logical large-list family, not just expected populated indexes. */
		private static String storage(final AssemblyStats2 stats) throws Exception{
			return "clist="+bins(stats, "clist")+";slist="+bins(stats, "slist")+
				";sclist1="+bins(stats, "sclist1")+";sclist2="+bins(stats, "sclist2")+
				";llist="+Arrays.toString(array(stats, field("llist")))+";triples="+triples(stats);
		}

		/** Sparse literal-friendly representation includes every nonzero logical bin. */
		private static String bins(final AssemblyStats2 stats, final String name) throws Exception{
			final long[] values=array(stats, field(name));
			final StringBuilder result=new StringBuilder();
			for(int i=0; i<values.length; i++){
				if(values[i]!=0){if(result.length()>0){result.append(',');} result.append(i).append(':').append(values[i]);}
			}
			return result.toString();
		}

		/** Only line-ending bytes vary; bases, record boundaries and quality lengths stay literal. */
		private static byte[] fastq(final String ending){
			assert("lf".equals(ending) || "cr".equals(ending) || "crlf".equals(ending)) : "Finite FASTQ encoding table: "+ending;
			return PARSER_FASTQ.replace("\n", "crlf".equals(ending) ? "\r\n" : "cr".equals(ending) ? "\r" : "\n")
				.getBytes(StandardCharsets.US_ASCII);
		}

		/** Empty logical sequence and quality lines must not consume the following AC record. */
		private static byte[] emptyFastq(final String ending){
			assert("lf".equals(ending) || "cr".equals(ending) || "crlf".equals(ending)) : "Finite empty-record encoding table: "+ending;
			return PARSER_EMPTY_FASTQ.replace("\n", "crlf".equals(ending) ? "\r\n" : "cr".equals(ending) ? "\r" : "\n")
				.getBytes(StandardCharsets.US_ASCII);
		}

		private final Path root;
		private final String phase;
		private int rows, correctRows, unexpected;
		private static final String PARSER_FASTQ="@r\nACG\n+\nIII\n@s\nTTAA\n+\nIIII\n";
		private static final String PARSER_EMPTY_FASTQ="@empty\n\n+\n\n@r\nAC\n+\nII\n";
		private static final long[] FASTQ_COUNTS={3, 1, 1, 2, 0, 0, 0, 0}, ZERO_COUNTS=new long[8];
		private static final long[] EMPTY_RECORD_COUNTS={1, 1, 0, 0, 0, 0, 0, 0};
		private static final long[] MIXED_BASELINE_COUNTS={0, 0, 0, 0, 1, 0, 0, 1};
		private static final long[] THREE_A={3, 0, 0, 0, 0, 0, 0, 0}, FOUR_A={4, 0, 0, 0, 0, 0, 0, 0};
		private static final String EMPTY_STORAGE="clist=;slist=;sclist1=;sclist2=;llist=[];triples=";
		private static final String TWO_STORAGE="clist=2:1;slist=2:1;sclist1=2:1;sclist2=2:2;llist=[];triples=";
		private static final String THREE_STORAGE="clist=3:1;slist=3:1;sclist1=3:1;sclist2=3:3;llist=[];triples=";
		private static final String FOUR_STORAGE="clist=;slist=;sclist1=;sclist2=;llist=[4];triples=4/1/4";
		private static final String FASTQ_STORAGE="clist=3:1;slist=3:1;sclist1=3:1;sclist2=3:3;llist=[4];triples=4/1/4";
		private static final String THREE_DUPLICATE="clist=3:3;slist=3:2;sclist1=3:3;sclist2=3:9;llist=[];triples=";
		private static final String FOUR_DUPLICATE="clist=;slist=;sclist1=;sclist2=;llist=[4, 4, 4];triples=4/1/4,4/2/8";
	}

	/** GC-header name checks with optional literal header accounting (JGI116/JGI117). */
	private static final class GcHeaderRows{
		/** Retains literal accounting checks for the established GC-header modes. */
		GcHeaderRows(final Path root_, final String phase_) throws Exception{
			this(root_, phase_, true);
		}

		/** Names-only modes record the observed header total without asserting its correctness. */
		GcHeaderRows(final Path root_, final String phase_, final boolean checkHeaderAccounting_) throws Exception{
			root=root_; phase=phase_;
			checkHeaderAccounting=checkHeaderAccounting_;
			headerlen=field("HEADERLENSUM");
			tlist=field("tlist");
			Files.write(root.resolve("gc-headers.tsv"), ("case\tphase\tformat\tpattern\tgcTarget\tcorrectCounts\tcorrectStorage"+
				"\tactualCounts\tactualStorage\tintendedGcName\tobservedGcNames\tgcRows\tintendedNameMatches\tintendedGcBytesMatch"+
				"\tliteralHeaderIds\tliteralHeaderBytes\theaderLenSumObserved\texceptionClass\texceptionMessage\tcountsStorageCorrect\texpectedObservationMatched\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		}

		/** Two formats, ordinary versus empty-surrounded records, and null versus file GC output. */
		void run() throws Exception{
			for(final boolean fasta : new boolean[]{true, false}){
				for(final boolean surrounded : new boolean[]{false, true}){
					for(final boolean fileGc : new boolean[]{false, true}){row(fasta, surrounded, fileGc);}
				}
			}
			writeSummary(8, checkHeaderAccounting ? "ASSEMBLY_STATS2_GC_HEADERS_CAPTURE_" : "ASSEMBLY_STATS2_GC_HEADER_NAMES_");
		}

		/** Four GC-file controls isolate filtered records, empty IDs, partial EOF and one-byte CRLF reads. */
		void runControls() throws Exception{
			require("fixed".equals(phase), "GC-header controls require repaired source.");
			row("fasta-filtered-min2", true, false, true,
				">short\nA\n>kept\nAC\n>tail\nT\n".getBytes(StandardCharsets.US_ASCII),
				"short|kept|tail", 13, 2, Integer.MAX_VALUE, "filtered-before+kept+filtered-after");
			row("fastq-filtered-min2", false, false, true,
				"@short\nA\n+\nI\n@kept\nAC\n+\nII\n@tail\nT\n+\nI\n".getBytes(StandardCharsets.US_ASCII),
				"@short|@kept|@tail", 16, 2, Integer.MAX_VALUE, "filtered-before+kept+filtered-after");
			row("fasta-empty-partial-crlf-cap1", true, false, true,
				">\r\n\r\n>kept\r\nAC\r\n>tail".getBytes(StandardCharsets.US_ASCII),
				"|kept|tail", 8, 1, 1, "empty-id+kept+unterminated-tail");
			row("fastq-empty-crlf-cap1", false, false, true,
				"@\r\n\r\n+\r\n\r\n@kept\r\nAC\r\n+\r\nII\r\n@tail\r\n\r\n+\r\n\r\n".getBytes(StandardCharsets.US_ASCII),
				"@|@kept|@tail", 11, 1, 1, "empty-id+kept+empty-tail");
			writeSummary(4, checkHeaderAccounting ? "ASSEMBLY_STATS2_GC_HEADER_CONTROLS_" : "ASSEMBLY_STATS2_GC_HEADER_NAME_CONTROLS_");
		}

		/** Strict fixed modes require literal totals; baseline and names-only modes record observations. */
		private void writeSummary(final int expectedRows, final String prefix) throws Exception{
			final boolean match=rows==expectedRows && correctRows==expectedRows && unexpected==0;
			final String summary=prefix+(match ? "PASS" : "FAIL")+" phase="+phase+
				" rows="+rows+" countsStorageCorrect="+correctRows+" gcNameChecks="+nameChecks+" gcNamesCorrect="+correctNames+
				" unexpectedObservations="+unexpected+" headerAccounting="+("baseline".equals(phase) || !checkHeaderAccounting ? "OBSERVED_ONLY" : "LITERAL_MATCH")+"\n";
			Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			System.out.print(summary);
			if(!match){throw new AssertionError(summary+"Inspect intended names, raw GC files and observed accounting in "+root);}
		}

		/** Keeps the correct name oracle separate from the predicted old prefix and the observed accounting scalar. */
		private void row(final boolean fasta, final boolean surrounded, final boolean fileGc) throws Exception{
			final String name=(fasta ? "fasta" : "fastq")+"-"+(surrounded ? "empty-kept-empty" : "kept")+"-"+(fileGc ? "file" : "null");
			final byte[] input=input(fasta, surrounded).getBytes(StandardCharsets.US_ASCII);
			final String intendedName=fasta ? "kept" : "@kept";
			final String literalIds=surrounded ? (fasta ? "empty-before|kept|empty-after" : "@empty-before|@kept|@empty-after") : intendedName;
			// AS2 skips FASTA > but counts FASTQ @ in its null-GC header branch; all IDs here are ASCII.
			final int literalHeaderBytes=surrounded ? (fasta ? 27 : 30) : (fasta ? 4 : 5);
			row(name, fasta, surrounded, fileGc, input, literalIds, literalHeaderBytes, 1, Integer.MAX_VALUE,
				surrounded ? "empty-before+kept+empty-after" : "kept");
		}

		/** Keeps exact GC/count/storage checks independent of the selected header-accounting policy. */
		private void row(final String name, final boolean fasta, final boolean surrounded, final boolean fileGc,
				final byte[] input, final String literalIds, final int literalHeaderBytes,
				final int minScaffold, final int chunkCap, final String pattern) throws Exception{
			assert(GC_COUNTS.length==8) : "AssemblyStats2's makeCharToNum defines eight count categories.";
			final Path gc=root.resolve(name+".gc");
			final String intendedName=fasta ? "kept" : "@kept";
			final byte[] intendedGc=fileGc ? gcBytes(intendedName) : new byte[0];
			final byte[] baselineGc=fileGc && surrounded ? gcBytes(fasta ? "empty-beforekept" : "@empty-before@kept") : intendedGc;
			final StaticState saved=new StaticState();
			long[] counts=null;
			String storage="unavailable", observedHeaderlen="unavailable";
			Throwable failure=null;
			try{
				final AssemblyStats2 stats=new AssemblyStats2(new String[]{"in=stdin", "cutoff=4", "minscaf="+minScaffold, "n=2",
					"format=1", "gcformat=1", "gcfraction=f", "showbbmap=f", "showspeed=f", "n_=f"});
				try{
					final InputStream stream=chunkCap==Integer.MAX_VALUE ? new ByteArrayInputStream(input) : new CappedInputStream(input, chunkCap);
					counts=fasta ? stats.countFasta(stream, fileGc ? gc.toString() : null) : stats.countFastq(stream, fileGc ? gc.toString() : null);
				}catch(final Throwable caught){failure=caught;}
				storage=ParserRows.storage(stats);
				observedHeaderlen=Long.toString(headerlen.getLong(null));
			}finally{saved.restore();}
			final byte[] actualGc=fileGc && Files.exists(gc) ? Files.readAllBytes(gc) : new byte[0];
			final String observedNames=gcNames(actualGc);
			final boolean correct=failure==null && Arrays.equals(GC_COUNTS, counts) && GC_STORAGE.equals(storage);
			final boolean namesMatch=fileGc && intendedName.equals(observedNames);
			final boolean bytesMatch=Arrays.equals(intendedGc, actualGc);
			final boolean accountingMatch=!checkHeaderAccounting || "baseline".equals(phase) || Integer.toString(literalHeaderBytes).equals(observedHeaderlen);
			final boolean match=correct && accountingMatch && Arrays.equals("baseline".equals(phase) ? baselineGc : intendedGc, actualGc);
			rows++; if(correct){correctRows++;} if(fileGc){nameChecks++; if(namesMatch){correctNames++;}} if(!match){unexpected++;}
			Files.write(root.resolve(name+".input"), input, StandardOpenOption.CREATE_NEW);
			Files.write(root.resolve(name+".gc-intended"), intendedGc, StandardOpenOption.CREATE_NEW);
			Files.write(root.resolve(name+".gc-actual"), actualGc, StandardOpenOption.CREATE_NEW);
			final String text=name+"\t"+phase+"\t"+(fasta ? "fasta" : "fastq")+"\t"+pattern+
				"\t"+(fileGc ? gc.toString() : "null")+"\t"+Arrays.toString(GC_COUNTS)+"\t"+GC_STORAGE+"\t"+Arrays.toString(counts)+"\t"+storage+
				"\t"+intendedName+"\t"+observedNames+"\t"+gcRows(new String(actualGc, StandardCharsets.UTF_8))+
				"\t"+(fileGc ? Boolean.toString(namesMatch) : "NA")+"\t"+bytesMatch+"\t"+literalIds+"\t"+literalHeaderBytes+
				"\t"+observedHeaderlen+"\t"+type(failure)+
				"\t"+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+"\t"+correct+"\t"+match+"\n";
			Files.write(root.resolve("gc-headers.tsv"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		}

		/** Literal LF records isolate header carry from the already-repaired CRLF state machine. */
		private static String input(final boolean fasta, final boolean surrounded){
			final String kept=fasta ? ">kept\nAC\n" : "@kept\nAC\n+\nII\n";
			if(!surrounded){return kept;}
			return fasta ? ">empty-before\n\n"+kept+">empty-after\n\n" :
				"@empty-before\n\n+\n\n"+kept+"@empty-after\n\n+\n\n";
		}

		/** Exact GCFORMAT1 bytes for a single kept AC record; no accounting value is derived from this oracle. */
		private static byte[] gcBytes(final String name){
			return ("#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC\n"+
				name+"\t2\t0.5000\t0.5000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.5000\n")
				.getBytes(StandardCharsets.UTF_8);
		}

		/** All emitted data-row names are retained; extra rows cannot hide behind a first-name check. */
		private static String gcNames(final byte[] bytes){
			final StringBuilder names=new StringBuilder();
			for(final String line : new String(bytes, StandardCharsets.UTF_8).split("\n")){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final int tab=line.indexOf('\t');
				if(names.length()>0){names.append('|');}
				names.append(tab<0 ? line : line.substring(0, tab));
			}
			return names.toString();
		}

		private final Path root;
		private final String phase;
		private final boolean checkHeaderAccounting;
		private final Field headerlen;
		private int rows, correctRows, nameChecks, correctNames, unexpected;
		private static final long[] GC_COUNTS={1, 1, 0, 0, 0, 0, 0, 0};
		private static final String GC_STORAGE="clist=2:1;slist=2:1;sclist1=2:1;sclist2=2:2;llist=[];triples=";
	}

	/** Fixed-only private-helper boundaries; no scanner, large input, writer or constructor is invoked. */
	private static final class RefillGuardRows{
		RefillGuardRows(final Path root_) throws Exception{
			root=root_;
			fasta=method("failIfFastaRecordExceedsMax");
			fastq=method("failIfFastqRecordExceedsMax");
			Files.write(root.resolve("refill-guards.tsv"), ("case\tphase\tcoverage\tformat\tinputBefore\tlim\tcarriedLength"+
				"\texpected\tactualOutcome\texceptionClass\texceptionMessage\tinputAfter\tbufferUnchanged\tmatched\tfullParserMaxReached\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		}

		/** Literal return/throw expectations do not scan or count the test bytes to derive an oracle. */
		void run() throws Exception{
			for(final boolean isFasta : new boolean[]{true, false}){
				row("exact-max-fit", isFasta, new byte[]{'A', 'A'}, 2, Integer.MAX_VALUE-2, "return");
				row("one-counted-byte-over", isFasta, new byte[]{'A', 'A', 'A'}, 3, Integer.MAX_VALUE-2, "overflow");
				row("prefix-excludes-overflow-byte", isFasta, new byte[]{'A', 'A', 'A'}, 2, Integer.MAX_VALUE-2, "return");
				row("already-max-one-byte", isFasta, new byte[]{'A'}, 1, Integer.MAX_VALUE, "overflow");
				row("byte14-exact-fit", isFasta, new byte[]{14}, 1, Integer.MAX_VALUE-1, "return");
				row("byte14-counts-overflow", isFasta, new byte[]{14, 'A'}, 2, Integer.MAX_VALUE-1, "overflow");
				row("space-exact-fit", isFasta, new byte[]{' '}, 1, Integer.MAX_VALUE-1, "return");
				row("space-counts-overflow", isFasta, new byte[]{' ', 'A'}, 2, Integer.MAX_VALUE-1, "overflow");
				row("entry-null-buffer", isFasta, null, 1, Integer.MAX_VALUE-1, "assertion");
				row("entry-zero-limit", isFasta, new byte[]{'A'}, 0, Integer.MAX_VALUE-1, "assertion");
				row("entry-negative-limit", isFasta, new byte[]{'A'}, -1, Integer.MAX_VALUE-1, "assertion");
				row("entry-limit-exceeds-buffer", isFasta, new byte[]{'A'}, 2, Integer.MAX_VALUE-1, "assertion");
				row("entry-negative-carried-length", isFasta, new byte[]{'A'}, 1, -1, "assertion");
			}
			row("header-defers-suffix", true, new byte[]{'A', 'A', '>', 'A', 'A', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("overflow-before-header", true, new byte[]{'A', 'A', 'A', '>'}, 4, Integer.MAX_VALUE-2, "overflow");
			row("negative-byte-defers-suffix", true, new byte[]{'A', 'A', -128, 'A', 'A', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("overflow-before-negative-byte", true, new byte[]{'A', 'A', 'A', -128}, 4, Integer.MAX_VALUE-2, "overflow");
			row("low-controls-not-counted", true, new byte[]{'A', 0, '\t', '\n', '\r', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("overflow-after-low-controls", true, new byte[]{'A', 0, '\t', '\n', '\r', 'A', 'A'}, 7, Integer.MAX_VALUE-2, "overflow");
			row("lf-stops-before-separator-quality", false, new byte[]{'A', 'A', '\n', '+', '\n', 'I', 'I', 'I'}, 8, Integer.MAX_VALUE-2, "return");
			row("cr-stops-before-separator-quality", false, new byte[]{'A', 'A', '\r', '+', '\n', 'I', 'I', 'I'}, 8, Integer.MAX_VALUE-2, "return");
			row("crlf-stops-before-separator-quality", false, new byte[]{'A', 'A', '\r', '\n', '+', '\n', 'I', 'I', 'I'}, 9, Integer.MAX_VALUE-2, "return");
			row("nul-defers-suffix", false, new byte[]{'A', 'A', 0, 'A', 'A', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("negative-byte-defers-suffix", false, new byte[]{'A', 'A', -128, 'A', 'A', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("tab-defers-suffix", false, new byte[]{'A', 'A', '\t', 'A', 'A', 'A'}, 6, Integer.MAX_VALUE-2, "return");
			row("overflow-before-line-end", false, new byte[]{'A', 'A', 'A', '\n', '+', '\n', 'I', 'I', 'I'}, 9, Integer.MAX_VALUE-2, "overflow");
			final boolean match=rows==39 && returned==17 && overflowed==12 && asserted==10 && unchanged==39 && unexpected==0;
			final String summary="ASSEMBLY_STATS2_REFILL_GUARDS_"+(match ? "PASS" : "FAIL")+" phase=fixed"+
				" coverage=private-refill-helper-reflection rows="+rows+" returned="+returned+" overflowed="+overflowed+
				" asserted="+asserted+" buffersUnchanged="+unchanged+" unexpected="+unexpected+" fullParserMaxReached=false\n";
			Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			System.out.print(summary);
			if(!match){throw new AssertionError(summary+"Inspect literal helper cases in "+root);}
		}

		/** Unwraps only the reflected target failure; reflection/setup errors remain fixture failures. */
		private void row(final String name, final boolean isFasta, final byte[] input, final int lim,
				final int carriedLength, final String expected) throws Exception{
			assert("return".equals(expected) || "overflow".equals(expected) || "assertion".equals(expected)) :
				"The refill-helper literal table declares one of three outcomes: "+name+" / "+expected;
			final byte[] before=input==null ? null : input.clone();
			Throwable failure=null;
			try{(isFasta ? fasta : fastq).invoke(null, input, lim, carriedLength);}
			catch(final InvocationTargetException caught){failure=caught.getCause();}
			final boolean overflow=failure!=null && failure.getClass()==IllegalArgumentException.class && GUARD_MESSAGE.equals(failure.getMessage());
			final boolean assertion=failure!=null && failure.getClass()==AssertionError.class;
			final boolean preserved=Arrays.equals(before, input);
			final boolean match=preserved && ("return".equals(expected) ? failure==null : "overflow".equals(expected) ? overflow : assertion);
			rows++; if(failure==null){returned++;} if(overflow){overflowed++;} if(assertion){asserted++;}
			if(preserved){unchanged++;} if(!match){unexpected++;}
			final String text=(isFasta ? "fasta-" : "fastq-")+name+"\tfixed\tprivate-refill-helper-reflection\t"+
				(isFasta ? "fasta" : "fastq")+"\t"+Arrays.toString(before)+"\t"+lim+"\t"+carriedLength+"\t"+expected+
				"\t"+(failure==null ? "return" : "throw")+"\t"+type(failure)+"\t"+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+
				"\t"+Arrays.toString(input)+"\t"+preserved+"\t"+match+"\tfalse\n";
			Files.write(root.resolve("refill-guards.tsv"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		}

		/** Requires the private static void signature so this mode cannot silently widen or replace the product API. */
		private static Method method(final String name) throws Exception{
			final Method method=AssemblyStats2.class.getDeclaredMethod(name, byte[].class, Integer.TYPE, Integer.TYPE);
			method.setAccessible(true);
			require(Modifier.isPrivate(method.getModifiers()) && Modifier.isStatic(method.getModifiers()) && method.getReturnType()==Void.TYPE,
				"Refill checks require a private static void helper(byte[], int, int): "+name);
			return method;
		}

		private final Path root;
		private final Method fasta, fastq;
		private int rows, returned, overflowed, asserted, unchanged, unexpected;
	}

	/** Exact advisory-byte observations, kept separate from minimum-heap or formatted-MB claims. */
	private static final class AdvisoryMemoryRows{
		AdvisoryMemoryRows(final Path root_, final String phase_) throws Exception{
			root=root_; phase=phase_;
			memory=AssemblyStats2.class.getDeclaredMethod("bbmapMemoryBytes", long[].class, Long.TYPE, Long.TYPE, Integer.TYPE);
			memory.setAccessible(true);
			require(Modifier.isStatic(memory.getModifiers()) && memory.getReturnType()==Long.TYPE,
				"The advisory oracle calls private static long bbmapMemoryBytes directly.");
			Files.write(root.resolve("bbmap-memory.tsv"), ("case\tphase\tinput\tcorrectCounts\tactualCounts\tscaffolds\tscaffoldsObserved"+
				"\theaderBytes\theaderBytesObserved\tk\tbaselineBytes\tintendedBytes\tactualBytes\tcountsPreserved\tparserCorrect"+
				"\tadvisoryCorrect\texpectedObservationMatched\tpaddingBefore\tpaddingAfter\tpaddingRestored\texceptionClass\texceptionMessage\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		}

		/** Literal returned-byte expectations were selected independently in the canonical advisory plan. */
		void run() throws Exception{
			row("defined-only", "", new long[]{20, 0, 0, 0, 0, 0, 0, 0}, 7386, 7386);
			row("iupac-only", "RRRRRRRRRR", new long[]{20, 0, 0, 0, 10, 0, 0, 0}, 7401, 7401);
			row("n-only", "NNNNNNNNNN", new long[]{20, 0, 0, 0, 0, 10, 0, 0}, 7386, 7401);
			row("iupac3-n7", "RRRNNNNNNN", new long[]{20, 0, 0, 0, 3, 7, 0, 0}, 7390, 7401);
			row("other-only", "??????????", new long[]{20, 0, 0, 0, 0, 0, 10, 0}, 7386, 7401);
			row("consumed-control-only", "++++++++++", new long[]{20, 0, 0, 0, 0, 0, 0, 10}, 7386, 7401);
			row("all-undefined", "RRRNNN??++", new long[]{20, 0, 0, 0, 3, 3, 2, 2}, 7390, 7401);
			final boolean match=rows==7 && parserRows==7 && advisoryRows==("baseline".equals(phase) ? 2 : 7) && unexpected==0 && restored;
			final String summary="ASSEMBLY_STATS2_ADVISORY_MEMORY_CAPTURE_"+(match ? "PASS" : "FAIL")+" phase="+phase+
				" rows="+rows+" parserCorrectRows="+parserRows+" advisoryCorrectRows="+advisoryRows+
				" unexpectedObservations="+unexpected+" dependencyStateRestored="+restored+"\n";
			Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			System.out.print(summary);
			if(!match){throw new AssertionError(summary+"Inspect literal byte oracles in "+root);}
		}

		/** The suffix is consumed FASTA sequence, not a header, line ending, separator or quality string. */
		private void row(final String name, final String suffix, final long[] correctCounts,
				final long baselineBytes, final long intendedBytes) throws Exception{
			assert(correctCounts.length==8) : "AssemblyStats2.makeCharToNum exposes eight count categories.";
			final String inputName=name+".input";
			final byte[] input=(">one\nAAAAAAAAAA"+suffix+"\n>four\nAAAAAAAAAA\n").getBytes(StandardCharsets.US_ASCII);
			Files.write(root.resolve(inputName), input, StandardOpenOption.CREATE_NEW);
			final int[] before=padding();
			final StaticState saved=new StaticState();
			long[] counts=null;
			long observedBytes=-1, observedHeader=-1;
			int observedScaffolds=-1;
			boolean countsPreserved=false, parserCorrect=false, paddingRestored=false;
			int[] after=before;
			Throwable failure=null;
			try{
				require(Arrays.equals(before, PADDING), "Advisory literals pin MID/START/END/MAX_LENGTH: "+Arrays.toString(before));
				final AssemblyStats2 stats=newStats();
				counts=stats.countFasta(new ByteArrayInputStream(input), null);
				observedHeader=field("HEADERLENSUM").getLong(null);
				observedScaffolds=((ArrayList<?>)field("tlist").get(stats)).size();
				parserCorrect=Arrays.equals(correctCounts, counts) && observedHeader==7 && observedScaffolds==2;
				if(parserCorrect){
					final long[] snapshot=counts.clone();
					observedBytes=((Long)memory.invoke(null, counts, 2L, 7L, 5)).longValue();
					countsPreserved=Arrays.equals(snapshot, counts);
				}
			}catch(final InvocationTargetException caught){failure=caught.getCause();}
			catch(final Throwable caught){failure=caught;}
			finally{
				after=padding();
				try{saved.restore();}
				finally{
					FastaToChromArrays2.MID_PADDING=before[0]; FastaToChromArrays2.START_PADDING=before[1];
					FastaToChromArrays2.END_PADDING=before[2]; FastaToChromArrays2.MAX_LENGTH=before[3];
					paddingRestored=Arrays.equals(before, padding());
				}
			}
			final boolean correct=failure==null && observedBytes==intendedBytes;
			final long expected="baseline".equals(phase) ? baselineBytes : intendedBytes;
			final boolean match=failure==null && parserCorrect && countsPreserved && Arrays.equals(before, after) && paddingRestored && observedBytes==expected;
			rows++; if(parserCorrect){parserRows++;} if(correct){advisoryRows++;} if(!match){unexpected++;} restored&=paddingRestored;
			final String text=name+"\t"+phase+"\t"+inputName+"\t"+Arrays.toString(correctCounts)+"\t"+Arrays.toString(counts)+
				"\t2\t"+observedScaffolds+"\t7\t"+observedHeader+"\t5\t"+baselineBytes+"\t"+intendedBytes+"\t"+observedBytes+
				"\t"+countsPreserved+"\t"+parserCorrect+"\t"+correct+"\t"+match+"\t"+Arrays.toString(before)+"\t"+Arrays.toString(after)+
				"\t"+paddingRestored+"\t"+type(failure)+"\t"+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+"\n";
			Files.write(root.resolve("bbmap-memory.tsv"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		}

		/** Snapshot order is explicit in the TSV and independent literal pin. */
		private static int[] padding(){
			return new int[]{FastaToChromArrays2.MID_PADDING, FastaToChromArrays2.START_PADDING,
				FastaToChromArrays2.END_PADDING, FastaToChromArrays2.MAX_LENGTH};
		}

		private final Path root;
		private final String phase;
		private final Method memory;
		private int rows, parserRows, advisoryRows, unexpected;
		private boolean restored=true;
		private static final int[] PADDING={300, 8000, 8000, 536670912};
	}

	/** Same positive-bulk-read cap contract as CountGCBaselineTest's accepted fragmented-FASTQ fixture. */
	private static final class CappedInputStream extends InputStream{
		CappedInputStream(final byte[] data_, final int cap_){
			require(data_!=null && cap_>0, "Bounded fixture stream requires data and a positive read cap.");
			data=data_; cap=cap_;
		}
		@Override
		public int read(){singleReads++; if(pos>=data.length){return -1;} delivered++; return data[pos++]&0xFF;}
		@Override
		public int read(final byte[] target, final int offset, final int length){
			if(target==null){throw new NullPointerException();}
			if(offset<0 || length<0 || offset>target.length-length){throw new IndexOutOfBoundsException();}
			bulkReads++;
			if(length==0){return 0;}
			if(pos>=data.length){return -1;}
			final int count=Math.min(Math.min(length, cap), data.length-pos);
			System.arraycopy(data, pos, target, offset, count);
			pos+=count; delivered+=count;
			return count;
		}
		private final byte[] data;
		private final int cap;
		private int pos, bulkReads, singleReads, delivered;
	}

	/** Emits a finite prefix once, then throws the retained owner read failure. */
	private static final class ThrowingInputStream extends InputStream{
		ThrowingInputStream(final byte[] data_, final IOException failure_){
			data=data_;
			failure=failure_;
		}
		@Override
		public int read() throws IOException{
			singleReads++;
			if(pos<data.length){
				delivered++;
				return data[pos++]&0xFF;
			}
			errors++;
			throw failure;
		}
		@Override
		public int read(final byte[] target, final int offset, final int length) throws IOException{
			if(target==null){throw new NullPointerException();}
			if(offset<0 || length<0 || offset>target.length-length){throw new IndexOutOfBoundsException();}
			bulkReads++;
			if(length==0){return 0;}
			if(pos>=data.length){
				errors++;
				throw failure;
			}
			final int count=Math.min(length, data.length-pos);
			System.arraycopy(data, pos, target, offset, count);
			pos+=count;
			delivered+=count;
			return count;
		}
		private final byte[] data;
		private final IOException failure;
		private int pos=0, bulkReads=0, singleReads=0, delivered=0, errors=0;
	}

	/** One fixed-only helper boundary row. */
	private static final class GuardRow{
		GuardRow(final String name_, final int input_, final String expected_){
			name=name_;
			input=input_;
			expected=expected_;
		}
		final String name, expected;
		final int input;
	}

	/** Writer-rescue result for the fixture's own captured writers. */
	private static final class Rescue{
		static final Rescue NONE=new Rescue(0, 0, null);
		Rescue(final int rescued_, final int errorFlags_, final Throwable failure_){
			rescued=rescued_;
			errorFlags=errorFlags_;
			failure=failure_;
		}
		final int rescued, errorFlags;
		final Throwable failure;
	}

	private final Path root;
	private final boolean baseline;
	private final Method guard;
	private static Field tlist;
	private final Field clist, slist, sclist1, sclist2, llist;
	private int scannerCases, failureCases, guardCases, guardUnavailable, writerLeaks, writerRescues;
	private int checks, failures;

	private static final byte[] FASTA=(">small\nAAA\n>edge\nAAAA\n>split\nAANNGGGG\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ=("@small\nAAA\n+\nIII\n@edge\nAAAA\n+\nIIII\n@split\nAANNGGGG\n+\nIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTA_PREFIX=">ok\nAC\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ_PREFIX="@ok\nAC\n+\nII\n".getBytes(StandardCharsets.US_ASCII);
	private static final String STDOUT_REUSE_MARKER="#fixture stdout marker after AssemblyStats2\n";
	private static final String GUARD_MESSAGE=
		"AssemblyStats2 cannot count a record longer than 2147483647 sequence characters.";
	private static final long[] EXPECTED_COUNTS={9, 0, 4, 0, 0, 2, 0, 0};
	private static final StaticRef[] STATIC_FIELDS={
		new StaticRef(AssemblyStats2.class, "GCFORMAT"),
		new StaticRef(AssemblyStats2.class, "GC_FRACTION"),
		new StaticRef(AssemblyStats2.class, "FORMAT"),
		new StaticRef(AssemblyStats2.class, "cutoff"),
		new StaticRef(AssemblyStats2.class, "LIMSUM"),
		new StaticRef(AssemblyStats2.class, "HEADERLENSUM"),
		new StaticRef(AssemblyStats2.class, "bbmapkmer"),
		new StaticRef(AssemblyStats2.class, "overwrite"),
		new StaticRef(AssemblyStats2.class, "append"),
		new StaticRef(AssemblyStats2.class, "useheader"),
		new StaticRef(AssemblyStats2.class, "addfilename"),
		new StaticRef(AssemblyStats2.class, "printL90"),
		new StaticRef(AssemblyStats2.class, "swapNL"),
		new StaticRef(AssemblyStats2.class, "printExtended"),
		new StaticRef(AssemblyStats2.class, "logSumOffset"),
		new StaticRef(AssemblyStats2.class, "powSumPower"),
		new StaticRef(AssemblyStats2.class, "logSumBase"),
		new StaticRef(AssemblyStats2.class, "squareLog"),
		new StaticRef(AssemblyStats2.class, "logPower"),
		new StaticRef(AssemblyStats2.class, "showspeed"),
		new StaticRef(AssemblyStats2.class, "printheadersize"),
		new StaticRef(AssemblyStats2.class, "skipDuplicateLines"),
		new StaticRef(AssemblyStats2.class, "N_UNDERSCORE"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMinLen"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMaxLen"),
		new StaticRef(AssemblyStats2.class, "alignedFraction"),
		new StaticRef(AssemblyStats2.class, "printAssemblyScore"),
		new StaticRef(AssemblyStats2.class, "lastLogSum"),
		new StaticRef(AssemblyStats2.class, "lastL90"),
		new StaticRef(AssemblyStats2.class, "lastL50"),
		new StaticRef(AssemblyStats2.class, "lastSize"),
		new StaticRef(AssemblyStats2.class, "lastContigs"),
		new StaticRef(AssemblyStats2.class, "lastMaxContig"),
		new StaticRef(AssemblyStats2.class, "outstream"),
		new StaticRef(ReadWrite.class, "USE_PIGZ"),
		new StaticRef(ReadWrite.class, "USE_UNPIGZ"),
		new StaticRef(Shared.class, "mainClass"),
		new StaticRef(Shared.class, "COMMAND_LINE")
	};
}
