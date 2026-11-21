package test.jgi;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.FileFormat;
import fileIO.ByteStreamWriter;
import jgi.BBMerge;
import parse.Parse;
import stream.Read;
import stream.Streamer;
import stream.SamLine;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Paired baseline/fixed evidence for the normal five-column insert report repair.
 * Baseline isolates the unsupported attachment factory before constructing workers.
 * Five finite cases retain literal inputs and paired no-report good/bad oracles.
 * Fixed phase compares actual report rows, checks normal CLI policies and runs one
 * owned worker over an empty middle batch. Other probes invoke NEW workers directly;
 * the fixed report path is also exercised through separate public commands.
 * Retry and output-correction witnesses distinguish selected-pair mismatches from
 * original/restored or subsequently corrected reads. They are predictions until run.
 * No model, graph, training export, throughput or generic shutdown claim is made.
 * @author Sucrose
 */
public final class BBMergeInsertReportTest{

	/** Requires phase=baseline|fixed, root and new outdir; fixed also needs the saved baseline fixture reference. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "case", "reference").contains(pair[0]), "Expected documented fixture flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false;
		assert(assertions=true);
		fixed="fixed".equals(options.get("phase"));
		require(assertions && (fixed || "baseline".equals(options.get("phase"))), "Use the launcher with -ea and phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite";
		if(mode.equals("probe")){
			require(Files.isDirectory(out), "Probe requires an existing fixture-owned case directory.");
			probe(root, out, spec(options.get("case")));return;
		}
		if(mode.equals("empty-batches")){require(fixed && Files.isDirectory(out), "Owned empty-batch probe is fixed-only.");emptyBatches(root, out);return;}
		require(mode.equals("suite") && !options.containsKey("case"), "Use suite or internal probe mode.");
		Files.createDirectory(out);
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(fixed ? reference!=null && text(reference.resolve("summary.txt")).startsWith("BBMERGE_INSERT_REPORT_BASELINE_PASS cases=5 ") : reference==null, "Fixed phase requires the accepted saved baseline fixture; baseline phase takes no reference.");
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(sourcePin()), "Fixture requires its reviewed production source.");
		require(pin(root.resolve("networks/bbmerge.bbnet.gz")).equals(MODEL_PIN), "Do not replace or retrain the model for report verification.");
		write(out.resolve("source.sha80"), sourcePin()+"\n");write(out.resolve("model.sha80"), MODEL_PIN+"\n");
		write(out.resolve("header.txt"), BBMerge.header());
		require(text(out.resolve("header.txt")).equals(HEADER), "The actual normal header must retain all five columns.");
		if(!fixed){isolatedDispatch(out);}else{formatterRows(out);write(out.resolve("reports.tsv"), "case\texit\tliteral_match\n");}
		write(out.resolve("cases.tsv"), "case\tprobe_exit\tcli_exit\tselected_insert\tselected_mismatches\texpected_status\tgood_pairs\tbad_pairs\n");
		for(final String name : NAMES){run(root, out, reference, spec(name));}
		if(fixed){controls(root, out);}
		require(children==(fixed ? 27 : 10), "Finite child inventory changed: "+children);
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(sourcePin()) && pin(root.resolve("networks/bbmerge.bbnet.gz")).equals(MODEL_PIN), "Source/model changed during collection.");
		final String summary=fixed ? "BBMERGE_INSERT_REPORT_FIXED_PASS cases=5 children=27 formatterRows=11 reportProcesses=5 cliControls=11 emptyBatchProbe=1 noReportControls=5 baselineBytesPreserved=true\n" : "BBMERGE_INSERT_REPORT_BASELINE_PASS cases=5 children=10 isolatedUnsupportedDispatch=1 headerColumns=5 noReportControls=5 retryWitness=1 eccoWitness=1 productionChanges=0\n";
		write(out.resolve("summary.txt"), summary);System.out.print(summary);
	}

	/** Calls only the unsupported factory; no input, good/bad writer or worker exists here. */
	private static void isolatedDispatch(final Path out) throws IOException{
		final Path target=out.resolve("unsupported.info");
		final FileFormat ff=FileFormat.testOutput(target.toString(), FileFormat.ATTACHMENT, ".info", true, false, false, true);
		require(ff.attachment(), "The isolated descriptor must reach the legacy ATTACHMENT route.");
		RuntimeException failure=null;
		try{WriterFactory.getStreamS(ff, null, null, null, 16, HEADER, false, 1);}
		catch(final RuntimeException e){failure=e;}
		require(failure!=null && failure.getMessage()!=null && failure.getMessage().startsWith("Unsupported file format: "), "Expected the specific unsupported factory rejection.");
		require(!Files.exists(target), "Isolated factory rejection must not create output.");
		write(out.resolve("dispatch.tsv"), "descriptor\texception\tmessage\toutput_exists\tprestarted_workers\nATTACHMENT\t"+failure.getClass().getName()+"\t"+failure.getMessage()+"\tfalse\t0_BY_CONSTRUCTION\n");
	}

	/** Retains exact inputs and literal no-report output oracles before starting children. */
	private static void run(final Path root, final Path out, final Path reference, final Spec spec) throws Exception{
		final Path work=Files.createDirectory(out.resolve(spec.name));
		final Read a=spec.first(), b=spec.second();
		write(work.resolve("A.fq"), fastq(a));write(work.resolve("B.fq"), fastq(b));
		final byte[] ab=Files.readAllBytes(work.resolve("A.fq")), bb=Files.readAllBytes(work.resolve("B.fq"));
		write(work.resolve("expected-report.tsv"), HEADER+spec.name+"/1\t0\t"+(spec.bad ? -1 : 100)+"\t"+(spec.bad ? "F" : (spec.selectedErrors==0 ? "P" : "I"))+"\t"+(spec.bad ? -1 : spec.selectedErrors)+"\n");
		write(work.resolve("expected-good1.fq"), spec.bad ? "" : fastq(spec.ecco ? spec.correctedFirst() : a));
		write(work.resolve("expected-good2.fq"), spec.bad ? "" : fastq(spec.ecco ? spec.correctedSecond() : b));
		write(work.resolve("expected-bad1.fq"), spec.bad ? fastq(a) : "");write(work.resolve("expected-bad2.fq"), spec.bad ? fastq(b) : "");
		if(reference!=null){for(final String file : new String[]{"A.fq", "B.fq", "expected-report.tsv", "expected-good1.fq", "expected-good2.fq", "expected-bad1.fq", "expected-bad2.fq"}){same(work.resolve(file), reference.resolve(spec.name).resolve(file), "Baseline literal/input changed: "+spec.name+"/"+file);}}
		final int probe=child(root, work, "probe", new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeinsertreporttest.sh").toString(), "mode=probe", "phase="+(fixed ? "fixed" : "baseline"), "root="+root, "outdir="+work, "case="+spec.name)));
		require(probe==0, "The selected-pair witness failed; inspect its saved observations before changing any expectation: "+spec.name);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbmerge.sh").toString(), "-Xmx128m", "-Xms32m"));command.addAll(arguments(work, spec));
		final int cli=child(root, work, "cli", command);
		require(cli==0, "Finite no-report public process failed: "+spec.name);
		for(final String channel : new String[]{"good1", "good2", "bad1", "bad2"}){
			require(Arrays.equals(Files.readAllBytes(work.resolve("expected-"+channel+".fq")), Files.readAllBytes(work.resolve(channel+".fq"))), "Literal no-report output mismatch: "+spec.name+"/"+channel);
			if(reference!=null){same(work.resolve(channel+".fq"), reference.resolve(spec.name).resolve(channel+".fq"), "No-report regression: "+spec.name+"/"+channel);}
		}
		if(fixed){
			final ArrayList<String> report=command(root, reportArguments(work, spec, work.resolve("report.tsv").toString()));
			final int code=child(root, work, "report", report);require(code==0, "Normal report CLI failed: "+spec.name);
			same(work.resolve("report.tsv"), work.resolve("expected-report.tsv"), "Wrong selected-pair report: "+spec.name);
			for(final String channel : new String[]{"good1", "good2", "bad1", "bad2"}){same(work.resolve("report-"+channel+".fq"), work.resolve("expected-"+channel+".fq"), "Enabling report changed read output: "+spec.name+"/"+channel);}
			append(out.resolve("reports.tsv"), spec.name+"\t0\ttrue\n");
		}
		require(Arrays.equals(ab, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(bb, Files.readAllBytes(work.resolve("B.fq"))), "Children changed fixture inputs.");
		append(out.resolve("cases.tsv"), spec.name+"\t"+probe+"\t"+cli+"\t"+(spec.bad ? BBMerge.RET_BAD : 100)+"\t"+spec.selectedErrors+"\t"+(spec.bad ? "F" : (spec.selectedErrors==0 ? "P" : "I"))+"\t"+(spec.bad ? 0 : 1)+"\t"+(spec.bad ? 1 : 0)+"\n");
	}

	/** Observes actual candidate and output boundaries on independent fresh literal pairs. */
	private static void probe(final Path root, final Path work, final Spec spec) throws Exception{
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(sourcePin()), "Probe production source differs from reviewed phase.");
		final BBMerge merge=new BBMerge(arguments(work, spec).toArray(new String[0]));
		require(!(Boolean)value(merge, "MAKE_VECTOR") && !(Boolean)value(merge, "TRAINING") && !(Boolean)value(merge, "TAG_CUSTOM"), "This baseline isolates the normal non-training path.");
		require((Boolean)value(merge, "qtrim2")==spec.retry && (Boolean)value(merge, "ecco")==spec.ecco && !(Boolean)value(merge, "join"), "Parsed retry/output mode differs from its witness.");
		require((Boolean)value(merge, "useRatioMode") && !(Boolean)value(merge, "useFlatMode") && !(Boolean)value(merge, "useEntropy") && !(Boolean)value(merge, "useEfilter") &&
			(Boolean)value(merge, "useQuality") && !(Boolean)value(merge, "overlapUsingQuality") && (Boolean)value(merge, "overlapWithoutQuality"), "The actual worker must use the planned unweighted, unfiltered ratio search.");
		require(value(merge, "tadpole")==null && value(merge, "bloomFilter")==null && !(Boolean)value(merge, "kmerBridge") && !(Boolean)value(merge, "qtrim1") && !(Boolean)value(merge, "changeQuality"), "Graph/correction/initial-trim state would hide the intended finite witness.");
		final String[] fields={"MAKE_VECTOR", "TRAINING", "TAG_CUSTOM", "qtrim1", "qtrim2", "qtrimLeft", "qtrimRight", "trimq", "ecco", "join", "changeQuality", "useRatioMode", "useFlatMode", "useEntropy", "useEfilter", "useQuality", "overlapUsingQuality", "overlapWithoutQuality", "pfilterRatio", "minReadLength", "minInsert", "minInsert0", "MIN_OVERLAPPING_BASES", "MIN_OVERLAPPING_BASES_0", "MIN_OVERLAPPING_BASES_RATIO_REDUCTION", "MAX_MISMATCHES_R", "MAX_RATIO", "THREADS"};
		final StringBuilder settings=new StringBuilder("field\tvalue\n");
		for(final String field : fields){final Object v=value(merge, field);settings.append(field).append('\t').append(v instanceof float[] ? Arrays.toString((float[])v) : v).append('\n');}
		write(work.resolve("parsed-settings.tsv"), settings.toString());
		final Read a=spec.first(), b=spec.second();link(a, b);
		final int inputErrors=BBMerge.countErrors(a, b, 100);
		final Object worker=worker(merge);
		final int result=(Integer)invoke(worker, "processReadPair", new Class<?>[]{Read.class, Read.class}, a, b);
		final int errors=(result>0 ? BBMerge.countErrors(a, b, result) : -1);
		write(work.resolve("selected.tsv"), "result\tinput_errors_at_insert100\tselected_errors\tselected_length1\tselected_length2\tscratch_errors\n"+result+"\t"+inputErrors+"\t"+errors+"\t"+a.length()+"\t"+b.length()+"\t"+a.errors+"\n");
		write(work.resolve("selected1.fq"), fastq(a));write(work.resolve("selected2.fq"), fastq(b));
		require(result==(spec.bad ? BBMerge.RET_BAD : 100) && errors==spec.selectedErrors, "Selected result/count did not reach its literal oracle: "+spec.name);
		require(inputErrors==spec.inputErrors, "Literal input mismatch count differs from the planned geometry.");
		require(Arrays.equals(a.bases, ascii(spec.retry ? F80 : F100+ADAPTER_A)) && Arrays.equals(b.bases, spec.second().bases), "Selected arrays or R2 raw orientation differ from the independent oracle.");
		if(spec.retry){
			require(a.length()==80 && b.length()==100 && inputErrors==20 && errors==0, "Retry must remove20low-quality mismatches before selecting insert100.");
			final Read x=spec.first(), y=spec.second();link(x, y);y.reverseComplementFast();
			final int initial=(Integer)invoke(worker(merge), "processReadPair_inner", new Class<?>[]{Read.class, Read.class}, x, y);
			write(work.resolve("initial.tsv"), "before_retry_result\n"+initial+"\n");
			require(initial<0, "A positive first search would not prove the qtrim2 retry path.");
		}
		final Read x=spec.first(), y=spec.second();link(x, y);
		final ArrayList<Read> good=new ArrayList<Read>(), bad=new ArrayList<Read>();
		final byte[][] originals=spec.retry ? new byte[4][] : null;
		final ByteBuilder captured=fixed ? new ByteBuilder() : null;
		final int returned=(Integer)(fixed ? invoke(worker(merge), "findOverlapInThread", new Class<?>[]{Read.class, byte[][].class, ArrayList.class, ArrayList.class, ByteBuilder.class}, x, originals, good, bad, captured) : invoke(worker(merge), "findOverlapInThread", new Class<?>[]{Read.class, byte[][].class, ArrayList.class, ArrayList.class}, x, originals, good, bad));
		if(fixed){write(work.resolve("captured-report.tsv"), HEADER+captured.toString());same(work.resolve("captured-report.tsv"), work.resolve("expected-report.tsv"), "Worker captured count after output mutation: "+spec.name);}
		final int after=returned>0 ? BBMerge.countErrors(x, y, returned) : -1;
		write(work.resolve("after-output.tsv"), "result\tmismatches\tlength1\tlength2\tgood_pairs\tbad_pairs\n"+returned+"\t"+after+"\t"+x.length()+"\t"+y.length()+"\t"+good.size()+"\t"+bad.size()+"\n");
		require(returned==result && good.size()==(spec.bad ? 0 : 1) && bad.size()==(spec.bad ? 1 : 0), "Output boundary changed candidate/routing unexpectedly.");
		require(after==(spec.bad ? -1 : (spec.ecco ? 0 : spec.inputErrors)), "Output restoration/correction witness differs from the literal mismatch count.");
		if(spec.ecco){require(errors==1 && after==0 && Arrays.equals(x.bases, ascii(F100)) && Arrays.equals(y.bases, ascii(RAW_F100)), "Output ecco must change the selected1mismatch into corrected0mismatches.");}
		if(spec.retry){require(errors==0 && after==20 && x.length()==100, "Restored no-ecco output must differ from selected80-base mismatch semantics.");}
		require(((Thread)worker).getState()==Thread.State.NEW, "A probe must not start a worker thread.");
		write(work.resolve("probe-summary.txt"), "BBMERGE_INSERT_REPORT_PROBE_PASS case="+spec.name+" selectedErrors="+errors+" afterOutputErrors="+after+" threadStarts=0\n");
	}

	/** Same explicit non-neural search settings for constructor probes and public commands. */
	private static ArrayList<String> arguments(final Path work, final Spec spec){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+work.resolve("A.fq"), "in2="+work.resolve("B.fq"), "out1="+work.resolve("good1.fq"), "out2="+work.resolve("good2.fq"), "outu1="+work.resolve("bad1.fq"), "outu2="+work.resolve("bad2.fq"),
			"t=1", "simd=f", "qin=33", "qout=33", "ow=f", "app=f", "merge=f", "ecco="+spec.ecco, "changequality=f", "qtrim=f", "maq=0", "minlength="+(spec.bad ? 200 : 1), "usequality=t", "ouq=f", "owq=t", "ratiomode=t", "flatmode=f", "entropy=f", "minoverlap="+(spec.retry ? 70 : 20), "minoverlap0="+(spec.retry ? 70 : 20), "ratiominoverlapreduction=0", "mininsert="+(spec.retry ? 100 : 1), "mininsert0="+(spec.retry ? 100 : 1), "efilter=f", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "nn=f", "maxratio=0.05", "maxmismatches=5"));
		if(spec.retry){args.add("qtrim2=r");args.add("trimq=20");}
		return args;
	}

	/** Owns only the returned process handle; no timeout, broad lookup or signaling. */
	private static int child(final Path root, final Path work, final String name, final ArrayList<String> command) throws Exception{
		children++;
		final StringBuilder text=new StringBuilder();for(int i=0; i<command.size(); i++){text.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve(name+".command.tsv"), text.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile()).redirectOutput(work.resolve(name+".stdout").toFile()).redirectError(work.resolve(name+".stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS");builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int code=builder.start().waitFor();write(work.resolve(name+".exit"), code+"\n");return code;
	}

	/** Eleven literal row contracts, including signed rejections and ignored stale scratch. */
	private static void formatterRows(final Path out) throws Exception{
		final int[] inserts={100, 100, -1, -2, -3, -4, -5, 0, 100, 200};
		final String[] expected={"format-0\t7\t100\tP\t0\n", "format-1\t7\t100\tI\t1\n", "format-2\t7\t-1\tF\t-1\n", "format-3\t7\t-1\tA\t-1\n", "format-4\t7\t-1\tF\t-1\n", "format-5\t7\t-1\tS\t-1\n", "format-6\t7\t-1\tF\t-1\n", "format-7\t7\t-1\tF\t-1\n", "7\t7\t100\tP\t0\n", "format-9\t7\t200\tP\t0\n"};
		final Method method=BBMerge.class.getDeclaredMethod("appendInsertReport", ByteBuilder.class, Read.class, int.class);method.setAccessible(true);
		final ByteBuilder all=new ByteBuilder();
		for(int i=0; i<inserts.length; i++){
			final Read a=new Read(ascii(F100), quality(100, 20), i==8 ? null : "format-"+i, 7, false);
			final char[] raw=RAW_F100.toCharArray();if(i==1){raw[49]='A';}
			final Read b=new Read(ascii(new String(raw)), quality(100, 20), "mate", 7, false);link(a, b);a.errors=42;b.errors=17;
			final byte[] ab=a.bases.clone(), bb=b.bases.clone(), aq=a.quality.clone(), bq=b.quality.clone();
			final ByteBuilder row=new ByteBuilder();method.invoke(null, row, a, inserts[i]);
			require(row.toString().equals(expected[i]), "Wrong literal formatter row: "+i+"/"+row);
			require(a.errors==42 && b.errors==17 && a.obj==null && b.obj==null && Arrays.equals(ab, a.bases) && Arrays.equals(bb, b.bases) && Arrays.equals(aq, a.quality) && Arrays.equals(bq, b.quality), "Report formatting changed reads or scratch.");
			all.append(row);
		}
		//Literal UTF-8 bytes cover Latin-1 and a supplementary code point without a shared encoder oracle.
		final Read unicode=new Read(ascii(F100), quality(100, 20), "caf\u00e9-\ud83e\uddea", 7, false);
		link(unicode, new Read(ascii(RAW_F100), quality(100, 20), "mate", 7, false));
		final ByteBuilder row=new ByteBuilder();method.invoke(null, row, unicode, 100);
		final byte[] expectedUtf8={'c', 'a', 'f', (byte)0xc3, (byte)0xa9, '-', (byte)0xf0, (byte)0x9f, (byte)0xa7, (byte)0xaa, '\t', '7', '\t', '1', '0', '0', '\t', 'P', '\t', '0', '\n'};
		require(Arrays.equals(row.toBytes(), expectedUtf8), "Report identifier must retain literal UTF-8 including the complete surrogate pair.");
		Files.write(out.resolve("formatter-utf8.tsv"), row.toBytes(), StandardOpenOption.CREATE_NEW);all.append(row);
		Files.write(out.resolve("formatter-rows.tsv"), all.toBytes(), StandardOpenOption.CREATE_NEW);
	}

	/** Exercises the changed worker's submission of an actual empty middle batch. */
	private static void emptyBatches(final Path root, final Path work) throws Exception{
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(FIXED_SOURCE_PIN), "Owned worker source differs from the reviewed fixed version.");
		final Spec perfect=spec("perfect"), imperfect=spec("imperfect");
		write(work.resolve("A.fq"), fastq(perfect.first()));write(work.resolve("B.fq"), fastq(perfect.second()));
		final BBMerge merge=new BBMerge(arguments(work, perfect).toArray(new String[0]));
		final FiniteInput input=new FiniteInput(perfect, imperfect);
		final ByteStreamWriter writer=new ByteStreamWriter(FileFormat.testOutput(work.resolve("worker.tsv").toString(), FileFormat.TEXT, ".txt", true, false, false, true));
		writer.forcePrint(HEADER);writer.start();input.start();
		final Object worker=worker(merge, input, writer);final Thread thread=(Thread)worker;
		thread.start();thread.join();input.close();
		final boolean error=writer.poisonAndWait();
		require(thread.getState()==Thread.State.TERMINATED && writer.getState()==Thread.State.TERMINATED && value(worker, "failure")==null && !error, "Owned finite worker/report writer failed.");
		require(((Number)value(writer, "nextJobID")).longValue()==3 && input.index==3, "Empty middle batch must be submitted with its dense ID.");
		require(text(work.resolve("worker.tsv")).equals(HEADER+"perfect/1\t0\t100\tP\t0\nimperfect/1\t0\t100\tI\t1\n"), "Owned batch report order/count changed.");
		write(work.resolve("summary.txt"), "BBMERGE_INSERT_EMPTY_BATCH_PASS dataBatches=2 emptyBatches=1 nextJobID=3 workersJoined=1 writerJoined=true\n");
	}

	/** Finite fixture-owned input; no parser, background producer, queue or injected failure. */
	private static final class FiniteInput implements Streamer{
		FiniteInput(final Spec first, final Spec last){
			final Read a=first.first(), b=first.second(), c=last.first(), d=last.second();link(a, b);link(c, d);
			batches.add(new ListNum<Read>(new ArrayList<Read>(Arrays.asList(a)), 0));batches.add(new ListNum<Read>(new ArrayList<Read>(), 1));batches.add(new ListNum<Read>(new ArrayList<Read>(Arrays.asList(c)), 2));
		}
		@Override public String fname(){return "fixture-owned-finite-batches";}
		@Override public void start(){require(!started, "Finite input starts once.");started=true;}
		@Override public void close(){closed=true;}
		@Override public boolean paired(){return true;}
		@Override public int pairnum(){return 0;}
		@Override public long readsProcessed(){return 0;}
		@Override public long basesProcessed(){return 0;}
		@Override public void setSampleRate(final float rate, final long seed){require(rate==1, "Finite input has literal batches, not a sampler.");}
		@Override public synchronized ListNum<Read> nextList(){require(started && !closed, "Consume owned finite input only between start and close.");return index<batches.size() ? batches.get(index++) : null;}
		@Override public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException("Fixture supplies reads only.");}
		@Override public boolean hasMore(){return index<batches.size();}
		@Override public boolean errorState(){return false;}
		private final ArrayList<ListNum<Read>> batches=new ArrayList<ListNum<Read>>();
		private int index;private boolean started, closed;
	}

	/** Public-process controls tied to report path, ordering, header and append changes. */
	private static void controls(final Path root, final Path out) throws Exception{
		final Path work=Files.createDirectory(out.resolve("controls"));final Spec spec=spec("perfect");
		write(work.resolve("A.fq"), fastq(spec.first()));write(work.resolve("B.fq"), fastq(spec.second()));
		final byte[] inputA=Files.readAllBytes(work.resolve("A.fq")), inputB=Files.readAllBytes(work.resolve("B.fq"));
		final Path requested=work.resolve("report#.tsv"), report=work.resolve("report1.tsv");
		final String row="perfect/1\t0\t100\tP\t0\n";
		control(root, work, "report-only", reportOnly(work, spec, requested.toString()), 0, null);
		require(text(report).equals(HEADER+row) && !Files.exists(requested), "Resolved report-only path/header differs.");
		write(work.resolve("before-append.tsv"), text(report));
		final ArrayList<String> append=reportOnly(work, spec, requested.toString());append.addAll(Arrays.asList("ow=t", "app=t"));
		control(root, work, "append", append, 0, null);require(text(report).equals(HEADER+row+row), "Append must preserve prior bytes and suppress the existing header.");
		final String saved=text(report);
		control(root, work, "overwrite-rejected", reportOnly(work, spec, requested.toString()), 1, "exists and overwrite=false");
		require(text(report).equals(saved), "Rejected overwrite changed the existing report.");
		final ArrayList<String> noheader=reportOnly(work, spec, work.resolve("noheader.tsv").toString());noheader.add("noheader=t");
		control(root, work, "noheader", noheader, 0, null);require(text(work.resolve("noheader.tsv")).equals(row), "noheader must suppress only the header.");
		write(work.resolve("empty-A.fq"), "");write(work.resolve("empty-B.fq"), "");
		final ArrayList<String> empty=reportOnly(work, spec, work.resolve("empty.tsv").toString());empty.add("in1="+work.resolve("empty-A.fq"));empty.add("in2="+work.resolve("empty-B.fq"));
		control(root, work, "empty-input", empty, 0, null);require(text(work.resolve("empty.tsv")).equals(HEADER), "Empty input must leave exactly the normal header.");
		final StringBuilder interleaved=new StringBuilder(), expected=new StringBuilder(HEADER);
		for(int i=0; i<2048; i++){
			final Spec pair=spec((i&1)==0 ? "perfect" : "imperfect");final Read a=pair.first(), b=pair.second();
			interleaved.append(fastq(new Read(a.bases, a.quality, "order-"+i+"/1", i, false))).append(fastq(new Read(b.bases, b.quality, "order-"+i+"/2", i, false)));
			expected.append("order-").append(i).append("/1\t").append(i).append("\t100\t").append((i&1)==0 ? "P\t0\n" : "I\t1\n");
		}
		write(work.resolve("interleaved.fq"), interleaved.toString());write(work.resolve("ordered-expected.tsv"), expected.toString());
		require(Files.size(work.resolve("interleaved.fq"))>2*262144, "Ordering input must span multiple pinned FastqScanStreamer byte buffers.");
		final ArrayList<String> ordered=reportOnly(work, spec, work.resolve("ordered.tsv").toString());require(ordered.remove("in2="+work.resolve("B.fq")), "Interleaved control must remove its separate R2 argument.");ordered.add("in1="+work.resolve("interleaved.fq"));ordered.add("t=2");
		control(root, work, "ordered-two-workers", ordered, 0, null);same(work.resolve("ordered.tsv"), work.resolve("ordered-expected.tsv"), "Ordered two-worker rows differ from literal IDs/counts.");
		final ArrayList<String> sampled=reportOnly(work, spec, work.resolve("sampled.tsv").toString());require(sampled.remove("in2="+work.resolve("B.fq")), "Sampled interleaved control must remove its separate R2 argument.");sampled.add("in1="+work.resolve("interleaved.fq"));sampled.add("samplerate=0");
		control(root, work, "sampled-empty-batches", sampled, 0, null);require(text(work.resolve("sampled.tsv")).equals(HEADER), "All-sampled-out input must finish without losing a required batch ID.");
		final Path unexpected=work.resolve("unexpected-good.fq");
		final ArrayList<String> training=reportOnly(work, spec, work.resolve("training.tsv").toString());training.add("training=t");training.add("out1="+unexpected);training.add("k=31");
		control(root, work, "training-rejected", training, 1, "training export through outinsert is not implemented");require(!Files.exists(unexpected) && !Files.exists(work.resolve("training.tsv")), "Training rejection started requested output.");
		write(work.resolve("alias1.fq"), fastq(spec.first()));
		final ArrayList<String> alias=reportOnly(work, spec, work.resolve("alias#.fq").toString());alias.add("in1="+work.resolve("alias1.fq"));alias.add("ow=t");alias.add("out1="+unexpected);
		control(root, work, "input-path-rejected", alias, 1, "was specified multiple times");require(!Files.exists(unexpected) && Arrays.equals(inputA, Files.readAllBytes(work.resolve("alias1.fq"))), "Resolved input/report collision changed files.");
		final Path collision=work.resolve("same1.tsv");
		final ArrayList<String> duplicate=reportOnly(work, spec, work.resolve("same#.tsv").toString());duplicate.add("out1="+collision);duplicate.add("ow=t");
		control(root, work, "output-path-rejected", duplicate, 1, "was specified multiple times");require(!Files.exists(collision), "Resolved output/report collision created output.");
		final ArrayList<String> neural=reportOnly(work, spec, work.resolve("neural.tsv").toString());neural.addAll(Arrays.asList("nn=t", "efilter=6", "pfilter=0.00004", "cutoff=0", "net="+root.resolve("networks/bbmerge.bbnet.gz")));
		control(root, work, "neural-report", neural, 0, null);require(text(work.resolve("neural.tsv")).equals(HEADER+row), "Ordinary vector mode was rejected or changed the selected normal report.");
		require(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), "Report controls changed input mates.");
		final Path batches=Files.createDirectory(out.resolve("empty-batches"));
		final int code=child(root, batches, "worker", new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeinsertreporttest.sh").toString(), "mode=empty-batches", "phase=fixed", "root="+root, "outdir="+batches)));
		require(code==0, "Actual worker empty-middle-batch contract failed.");
		write(work.resolve("summary.txt"), "BBMERGE_INSERT_REPORT_CONTROLS_PASS cliCases=11 rejected=4 orderedPairs=2048 ownedWorkerEmptyMiddleBatch=true neuralCutoff=0 defaultNeuralAccuracy=NOT_TESTED\n");
	}

	private static ArrayList<String> command(final Path root, final ArrayList<String> arguments){final ArrayList<String> result=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbmerge.sh").toString(), "-Xmx128m", "-Xms32m"));result.addAll(arguments);return result;}
	private static ArrayList<String> reportArguments(final Path work, final Spec spec, final String path){final ArrayList<String> args=arguments(work, spec);args.addAll(Arrays.asList("out1="+work.resolve("report-good1.fq"), "out2="+work.resolve("report-good2.fq"), "outu1="+work.resolve("report-bad1.fq"), "outu2="+work.resolve("report-bad2.fq"), "outinsert="+path));return args;}
	private static ArrayList<String> reportOnly(final Path work, final Spec spec, final String path){final ArrayList<String> args=arguments(work, spec);args.addAll(Arrays.asList("out1=null", "out2=null", "outu1=null", "outu2=null", "outinsert="+path, "ordered=t"));return args;}
	/** Expected rejection is bounded to the named constructor failure, before output workers. */
	private static void control(final Path root, final Path work, final String name, final ArrayList<String> args, final int expectedExit, final String diagnostic) throws Exception{
		final int code=child(root, work, name, command(root, args));
		require(expectedExit==0 ? code==0 : code!=0, "Wrong normal report control exit: "+name+"/"+code);
		if(diagnostic!=null){final String error=text(work.resolve(name+".stderr"));require(error.contains(diagnostic) && !error.contains("Started output threads."), "Expected early constructor rejection: "+name);}
	}

	/** Pure literal geometry; no production reverse complement or join builds an oracle. */
	private static final class Spec{
		Spec(final String name_){
			name=name_;retry=name.equals("retry");ecco=name.equals("ecco");bad=name.equals("reject");inputErrors=retry ? 20 : (name.equals("perfect") ? 0 : 1);selectedErrors=bad ? -1 : (retry ? 0 : inputErrors);
			require(F100.length()==100 && F80.length()==80 && F100.startsWith(F80) && RAW_F100.length()==100 && BAD_TAIL20.length()==20 && ADAPTER_A.length()==60 && ADAPTER_B_RAW.length()==60,
				"Literal geometry must preserve insert100, selected80-base retry and160-base adapter controls.");
			require(F100.charAt(50)=='G' && RAW_F100.charAt(49)=='C', "The selected mismatch is fragment50 G/T, raw mate2 site49 C/A.");
		}
		Read first(){
			final byte[] q=quality(retry ? 100 : 160, 20);
			if(retry){Arrays.fill(q, 0, 80, (byte)40);Arrays.fill(q, 80, 100, (byte)0);}
			else{Arrays.fill(q, 100, 160, (byte)40);if(inputErrors>0){q[50]=30;}}
			return read(retry ? F80+BAD_TAIL20 : F100+ADAPTER_A, q, name+"/1");
		}
		Read second(){
			if(retry){return read(RAW_F100, quality(100, 40), name+"/2");}
			final char[] bases=(RAW_F100+ADAPTER_B_RAW).toCharArray();if(inputErrors>0){bases[49]='A';}
			final byte[] q=quality(160, 20);Arrays.fill(q, 100, 160, (byte)40);return read(new String(bases), q, name+"/2");
		}
		Read correctedFirst(){final byte[] q=quality(100, 20);q[50]=30;return read(F100, q, name+"/1");}
		Read correctedSecond(){return read(RAW_F100, quality(100, 20), name+"/2");}
		final String name;final boolean retry, ecco, bad;final int inputErrors, selectedErrors;
	}
	private static Spec spec(final String name){require(Arrays.asList(NAMES).contains(name), "Unknown finite case: "+name);return new Spec(name);}
	private static Read read(final String bases, final byte[] quality, final String name){return new Read(ascii(bases), quality, name, 0, false);}
	private static byte[] ascii(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static byte[] quality(final int length, final int q){final byte[] result=new byte[length];Arrays.fill(result, (byte)q);return result;}
	private static void link(final Read a, final Read b){a.mate=b;b.mate=a;b.setPairnum(1);}
	/** Independent raw ASCII33 serialization used only for literal input/output files. */
	private static String fastq(final Read r){assert(r!=null && r.bases!=null && r.quality!=null && r.length()==r.quality.length) : "Literal ASCII33 output requires a complete length-matched fixture read.";final StringBuilder text=new StringBuilder("@").append(r.id).append('\n').append(new String(r.bases, StandardCharsets.US_ASCII)).append("\n+\n");for(final byte q : r.quality){text.append((char)(q+33));}return text.append('\n').toString();}
	private static Object worker(final BBMerge merge) throws Exception{return worker(merge, null, null);}
	private static Object worker(final BBMerge merge, final Streamer input, final ByteStreamWriter report) throws Exception{final Class<?> type=Class.forName("jgi.BBMerge$MateThread");final Constructor<?> c=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, fixed ? ByteStreamWriter.class : Writer.class, boolean.class, boolean.class);c.setAccessible(true);return c.newInstance(merge, input, null, null, report, false, false);}
	private static Object value(final Object owner, final String name) throws Exception{final Field f=owner.getClass().getDeclaredField(name);f.setAccessible(true);return f.get(owner);}
	private static Object invoke(final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{final Method m=owner.getClass().getDeclaredMethod(name, types);m.setAccessible(true);return m.invoke(owner, args);}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required fixture path.");return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}
	private static String text(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	private static void same(final Path first, final Path second, final String message) throws IOException{require(Arrays.equals(Files.readAllBytes(first), Files.readAllBytes(second)), message);}
	private static String sourcePin(){return fixed ? FIXED_SOURCE_PIN : SOURCE_PIN;}
	/** Prints only the final80bits of a small source/model file digest. */
	private static String pin(final Path path) throws Exception{final byte[] bytes=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(path));final StringBuilder text=new StringBuilder();for(int i=bytes.length-10; i<bytes.length; i++){text.append(Character.forDigit((bytes[i]>>>4)&15, 16)).append(Character.forDigit(bytes[i]&15, 16));}return text.toString();}
	private static final String[] NAMES={"perfect", "imperfect", "retry", "ecco", "reject"};
	private static final String HEADER="#id\tnumericID\tinsert\tstatus\tmismatches\n";
	private static final String SOURCE_PIN="989c3ef47aebab252185", MODEL_PIN="24998a8aa257552d4a22";
	private static final String FIXED_SOURCE_PIN="5a4527a67db65afb8ab1";
	private static boolean fixed;
	private static int children;
	private static final String F100="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCAGTACCGTATCATGCTAGGTCCGATCTACGATACGGTTCAGCAGTACGATG";
	private static final String F80="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCAGTACCGTATCATGCTAGGTCCGATCTACGA";
	private static final String BAD_TAIL20="ATGCCAAGTCGTCATGCTAC";
	private static final String RAW_F100="CATCGTACTGCTGAACCGTATCGTAGATCGGACCTAGCATGATACGGTACTGCTATCGGATCAGGTCATCGATCGAAGCTAGCATGTCGACATGCAACGT";
	private static final String ADAPTER_A="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTG";
	private static final String ADAPTER_B_RAW="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC";
}
