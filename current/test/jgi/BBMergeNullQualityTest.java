package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.FileFormat;
import parse.Parse;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;

/** Three gated, fresh-JVM public BBMerge cases isolate null-quality restoration.
 * Literal 40/56 inputs encode insert 40. A default-quality FASTA control precedes
 * the quality-disabled FASTA witness and a nonnull FASTQ restoration control.
 * Baseline accepts only the wrapped restoreQualities NullPointerException witness;
 * fixed requires a successful literal merge. FASTQ consensus quality arithmetic has
 * no invented oracle: its native output must be byte-identical to the saved baseline.
 * Extension filling, mixed-quality pairs and general worker cleanup are outside scope.
 * @author Sucrose
 */
public final class BBMergeNullQualityTest{

	/** phase=baseline|fixed, root=RUNTIME, outdir=NEW; fixed requires reference=BASELINE. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("phase", "root", "outdir", "reference").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergenullqualitytest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed requires a successful baseline reference.");
		require(Files.isRegularFile(root.resolve("bbmerge.sh")), "Set root to the single BBTools runtime.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBMERGE_NULL_QUALITY_PASS phase=baseline version=1 cases=3 children=3 controls=2 discrepancies=1 checks=[0-9]+ failures=0\\n"), "Reference must be the successful three-case null-quality baseline.");
		}
		final BBMergeNullQualityTest test=new BBMergeNullQualityTest(root, out, reference, phase);
		try{test.exercise();}
		finally{test.saveDispatch();}
		test.finish();
	}

	private BBMergeNullQualityTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws IOException{
		root=root_; out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("inputs.tsv"), "mate\tbases\tinsert\tFASTQ_input_quality\nR1\t"+FRAGMENT+"\t40\t30\nR2\t"+RAW_INPUT+"\t40\t30\n");
		write(out.resolve("oracles.tsv"), "case\tbaseline\tfixed\tmerged_bases\tquality_policy\n"+
			"fasta-default\tmerge40\tmerge40\t"+FRAGMENT+"\tnative_FASTA_bytes\n"+
			"fasta-retained\twrapped_restoreQualities_NPE\tmerge40\t"+FRAGMENT+"\tnative_FASTA_bytes\n"+
			"fastq-retained\tmerge40\tmerge40\t"+FRAGMENT+"\texact_saved_baseline_FASTQ_bytes\n");
		final StringBuilder settings=new StringBuilder("case\targument\n");
		for(final Spec spec : SPECS){
			for(final String arg : COMMON){settings.append(spec.name).append('\t').append(arg).append('\n');}
			settings.append(spec.name).append("\tusequality=").append(spec.quality ? 't' : 'f').append('\n');
		}
		write(out.resolve("settings.tsv"), settings.toString());
		final String header="case\texit\tresult\n";
		write(out.resolve("cases.tsv"), header); write(out.resolve("controls.tsv"), header);
	}

	/** A failed prerequisite throws before any later child is launched. */
	private void exercise() throws Exception{
		for(int i=0; i<SPECS.length; i++){run(i, SPECS[i]);}
		check(cases==3 && children==3 && controls==2, "Keep exactly three fresh public children and two unaffected controls.");
		check(discrepancies==(baseline ? 1 : 0), "Expected only the baseline null-restoration discrepancy.");
		if(reference!=null){compareReference();}
	}

	private void run(final int index, final Spec spec) throws Exception{
		final Path work=out.resolve(spec.name); Files.createDirectory(work);
		final String suffix=spec.fastq ? ".fq" : ".fa";
		final String input1=spec.fastq ? fastq(FRAGMENT) : fasta(FRAGMENT);
		final String input2=spec.fastq ? fastq(RAW_INPUT) : fasta(RAW_INPUT);
		write(work.resolve("input1"+suffix), input1); write(work.resolve("input2"+suffix), input2);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbmerge.sh").toString(), "-Xmx128m", "-Xms32m", "-ea"));
		command.addAll(Arrays.asList(COMMON));
		command.addAll(Arrays.asList("usequality="+(spec.quality ? "t" : "f"), "in1="+work.resolve("input1"+suffix),
			"in2="+work.resolve("input2"+suffix), "out="+work.resolve("merged"+suffix)));
		final StringBuilder invocation=new StringBuilder("index\targument\n");
		for(int i=0; i<command.size(); i++){invocation.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve("command.tsv"), invocation.toString());
		final Process process=new ProcessBuilder(command).directory(root.toFile())
			.redirectOutput(work.resolve("child.stdout").toFile()).redirectError(work.resolve("child.stderr").toFile()).start();
		children++; statuses[index]="STARTED";
		process.getOutputStream().close();
		final int exit=process.waitFor();
		write(work.resolve("child.exit"), exit+"\n"); statuses[index]="RETURNED";
		check(read(work.resolve("input1"+suffix)).equals(input1) && read(work.resolve("input2"+suffix)).equals(input2), spec.name+": child changed literal inputs.");
		final String stderr=read(work.resolve("child.stderr"));
		final String result;
		if(spec.witness && baseline){
			check(exit==1, "REACHABILITY_FAIL: only the ordinary Java failure exit is expected for the baseline witness, got "+exit);
			check(expectedNullCause(stderr), "REACHABILITY_FAIL: expected only the wrapped worker restoreQualities NPE; inspect actual child.stderr.");
			check(Files.isRegularFile(work.resolve("merged.fa")) && Files.size(work.resolve("merged.fa"))==0, "The failed one-pair worker must not emit a merged record.");
			discrepancies++; result="EXPECTED_RESTORE_NULL_POINTER";
		}else{
			check(exit==0, spec.name+": unexpected child exit "+exit+"; later work is gated.");
			check(!stderr.contains("Exception in thread") && !stderr.contains("Caused by:"), spec.name+": unexpected exception stack despite a zero exit.");
			if(spec.fastq){
				final ArrayList<Read> reads=readFastq(work.resolve("merged.fq"));
				check(reads.size()==1, "REACHABILITY_FAIL: native FASTQ reader must find exactly one merged record.");
				final Read merged=reads.get(0);
				check(merged.id.equals("insert=40") && merged.mate==null && Arrays.equals(merged.bases, bytes(FRAGMENT)), "REACHABILITY_FAIL: native FASTQ merged identity/bases differ from the literal 40-base fragment.");
				check(merged.quality!=null && merged.quality.length==40, "Nonnull-quality control must retain a 40-base quality output; consensus values have no invented arithmetic oracle.");
			}else{
				check(read(work.resolve("merged.fa")).equals(fasta(FRAGMENT)), "REACHABILITY_FAIL: native FASTA output differs from the literal merged fragment.");
				if(spec.witness){same(work.resolve("merged.fa"), out.resolve("fasta-default/merged.fa"));}
			}
			result="MERGED_LITERAL_40";
		}
		final String row=spec.name+"\t"+exit+"\t"+result+"\n";
		append(out.resolve("cases.tsv"), row);
		if(!spec.witness){controls++; append(out.resolve("controls.tsv"), row);}
		cases++; statuses[index]="PASSED";
	}

	/** Requires the public main wrapper and the first and only cause's actual top frame. */
	private static boolean expectedNullCause(final String text){
		if(!text.contains("java.lang.RuntimeException: BBMerge worker thread failed.") || !text.contains("at jgi.BBMerge.main(") || text.contains("Suppressed:")){return false;}
		int causes=0; boolean top=false, frame=false;
		for(final String raw : text.split("\\r?\\n")){
			final String line=raw.trim();
			if(line.startsWith("Caused by:")){
				causes++;
				if(!line.equals("Caused by: java.lang.NullPointerException") && !line.startsWith("Caused by: java.lang.NullPointerException:")){return false;}
				top=true;
			}else if(top && line.startsWith("at ")){
				frame=line.startsWith("at jgi.BBMerge$MateThread.restoreQualities("); top=false;
			}
		}
		return causes==1 && frame;
	}

	/** Uses the native reader; borrowed decoding/interleaving settings are restored. */
	private static ArrayList<Read> readFastq(final Path file){
		final byte offset=FASTQ.ASCII_OFFSET;
		final boolean detect=FASTQ.DETECT_QUALITY, custom=FASTQ.PARSE_CUSTOM;
		final boolean test=FASTQ.TEST_INTERLEAVED, force=FASTQ.FORCE_INTERLEAVED;
		try{
			FASTQ.ASCII_OFFSET=33; FASTQ.DETECT_QUALITY=false; FASTQ.PARSE_CUSTOM=false;
			FASTQ.TEST_INTERLEAVED=false; FASTQ.FORCE_INTERLEAVED=false;
			return ReadInputStream.toReads(file.toString(), FileFormat.FASTQ, -1);
		}finally{
			FASTQ.ASCII_OFFSET=offset; FASTQ.DETECT_QUALITY=detect; FASTQ.PARSE_CUSTOM=custom;
			FASTQ.TEST_INTERLEAVED=test; FASTQ.FORCE_INTERLEAVED=force;
		}
	}

	private void compareReference() throws IOException{
		for(final String file : new String[]{"inputs.tsv", "oracles.tsv", "settings.tsv", "controls.tsv"}){same(out.resolve(file), reference.resolve(file));}
		for(final Spec spec : SPECS){
			final String suffix=spec.fastq ? ".fq" : ".fa";
			for(final String mate : new String[]{"input1", "input2"}){same(out.resolve(spec.name).resolve(mate+suffix), reference.resolve(spec.name).resolve(mate+suffix));}
			if(!spec.witness){same(out.resolve(spec.name).resolve("merged"+suffix), reference.resolve(spec.name).resolve("merged"+suffix));}
		}
	}

	private void saveDispatch() throws IOException{
		final StringBuilder b=new StringBuilder("case\tstatus\n");
		for(int i=0; i<SPECS.length; i++){b.append(SPECS[i].name).append('\t').append(statuses[i]).append('\n');}
		write(out.resolve("dispatch.tsv"), b.toString());
	}
	private void finish() throws IOException{
		final String summary="BBMERGE_NULL_QUALITY_PASS phase="+phase+" version=1 cases="+cases+" children="+children+
			" controls="+controls+" discrepancies="+discrepancies+" checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}
	private void same(final Path a, final Path b) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), "Saved baseline input/control output changed: "+a);}
	private void check(final boolean condition, final String why){checks++; if(!condition){throw new IllegalStateException(why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static String fasta(final String bases){return ">insert=40\n"+bases+"\n";}
	private static String fastq(final String bases){final char[] q=new char[bases.length()]; Arrays.fill(q, '?'); return "@insert=40\n"+bases+"\n+\n"+new String(q)+"\n";}

	private static final class Spec{
		Spec(final String name_, final boolean fastq_, final boolean quality_, final boolean witness_){name=name_; fastq=fastq_; quality=quality_; witness=witness_;}
		final String name;
		final boolean fastq, quality, witness;
	}
	private final Path root, out, reference;
	private final String phase;
	private final boolean baseline;
	private int cases, children, controls, discrepancies, checks;
	private final String[] statuses={"SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE"};
	private static final Spec[] SPECS={new Spec("fasta-default", false, true, false), new Spec("fasta-retained", false, false, true), new Spec("fastq-retained", true, false, false)};
	private static final String[] COMMON={"t=1", "simd=f", "bf1=t", "nn=f", "qtrim=f", "maq=0", "minlength=1", "ordered=t", "ow=f", "app=f", "qin=33", "qout=33",
		"ratiomode=t", "flatmode=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0", "mininsert=1", "mininsert0=1",
		"ouq=f", "owq=t", "efilter=f", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "ecc=f", "join=t", "mix=f",
		"parsecustom=t", "samplerate=1", "sampleseed=1", "reads=-1", "recalibrate=f", "quantize=f", "train=f"};
	private static final String FRAGMENT="ACGTTGCATGTCAGTACGATCGTACCTGAGCTAGTCGATC";
	private static final String RAW_INPUT="GATCGACTAGCTCAGGTACGATCGTACTGACATGCAACGTAGATCGGAAGAGCACA";
}
