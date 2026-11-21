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
import java.util.List;
import java.util.Map;

import fileIO.FileFormat;
import jgi.BBMerge;
import parse.Parse;
import stream.FASTQ;
import stream.Read;
import stream.StreamerFactory;

/** Finite retained-quality coordinate witnesses for JGI095.
 * Five direct cases and two real public overlap calls use fresh linked reads.
 * One BBDuk CLI is dispatched only after the direct and public-wrapper checks pass.
 * Literal source and corrected arrays distinguish the known defect from a control;
 * no overlap search or quality-slicing expression generates an expected result.
 * Fixed runs require the saved successful baseline and preserve its control rows.
 * This does not exercise a full Clump pipeline or measure merger accuracy/speed.
 * @author Sucrose
 */
public final class BBMergeRetainedQualityTest{

	/** phase=baseline|fixed, root=RUNTIME, outdir=NEW, fixed-only reference=BASELINE. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("phase", "root", "outdir", "reference").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergeretainedqualitytest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed requires reference=BASELINE_FIXTURE_DIRECTORY.");
		require(Files.isRegularFile(root.resolve("bbduk.sh")), "Set root to the single runtime.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBMERGE_RETAINED_QUALITY_PASS phase=baseline version=1 cases=8 direct=5 wrappers=2 cli=1 children=1 discrepancies=3 checks=81 failures=0\\n"), "Reference must be a successful eight-case baseline.");
		}
		final BBMergeRetainedQualityTest test=new BBMergeRetainedQualityTest(root, out, reference, phase);
		test.exercise();
		test.check(test.cases==8 && test.direct==5 && test.wrappers==2 && test.cli==1 && test.children==1,
			"Expected five direct cases, two public wrappers and one actual BBDuk child.");
		test.check(test.discrepancies==(test.baseline ? 3 : 0), "Expected three baseline quality-coordinate witnesses and no fixed discrepancies.");
		if(reference!=null){
			test.same(out.resolve("inputs.tsv"), reference.resolve("inputs.tsv"), "Literal inputs changed.");
			test.same(out.resolve("oracles.tsv"), reference.resolve("oracles.tsv"), "Literal expected outputs changed.");
			test.same(out.resolve("controls.tsv"), reference.resolve("controls.tsv"), "Saved baseline controls changed.");
			for(final String file : new String[]{"input1.fq", "input2.fq"}){test.same(out.resolve(file), reference.resolve(file), "CLI input bytes changed: "+file);}
			if(test.failures==0){test.same(out.resolve("output1.fq"), reference.resolve("output1.fq"), "Unchanged CLI R1 output bytes changed.");}
		}
		final String summary="BBMERGE_RETAINED_QUALITY_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" version=1 cases="+test.cases+" direct="+test.direct+" wrappers="+test.wrappers+" cli="+test.cli+
			" children="+test.children+" discrepancies="+test.discrepancies+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect retained prerequisites and literal arrays; do not fit inputs to output.");}
	}

	private BBMergeRetainedQualityTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws IOException{
		root=root_; out=out_; reference=reference_; baseline=phase_.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("inputs.tsv"), "case\tapi\tinsert\tchangeQuality\tbefore\n");
		write(out.resolve("oracles.tsv"), "case\tsource\tcorrected\n");
		final String header="case\treturned\tactual\tsourceMatches\tcorrectedMatches\tcompleted\n";
		write(out.resolve("cases.tsv"), header); write(out.resolve("controls.tsv"), header);
	}

	/** Changes only owned public settings and always restores the borrowed adapter list. */
	private void exercise() throws Exception{
		final boolean oldQuality=BBMerge.changeQuality;
		final ArrayList<byte[]> oldAdapters=BBMerge.staticAdapterList;
		try{
			check(Read.MAX_MERGE_QUALITY>=37 && Read.MAX_CALLED_QUALITY>=39, "Literal qualities must not be capped by the current runtime.");
			BBMerge.staticAdapterList=BBMerge.getAdapterList("default");
			run("short-r2-retained", 0, small(true), 3, false,
				expect("ACG", q(12, 16, 20), "ACG", q(3, 5, 24), 0, 0, 0, -1),
				expect("ACG", q(12, 16, 20), "ACG", q(24, 28, 32), 0, 0, 0, -1), false, false);
			final Expected full=expect("ACG", q(12, 16, 20), "ACG", q(24, 28, 32), 0, 0, 0, -1);
			run("full-r2-control", 0, pair("ACG", q(12, 16, 20), "ACG", q(24, 28, 32)), 3, false, full, full, true, false);
			final Expected consensus=expect("ACG", q(27, 32, 37), "ACG", q(27, 32, 37), 0, 0, 0, -1);
			run("consensus-quality-control", 0, small(true), 3, true, consensus, consensus, true, false);
			final Expected absent=expect("ACG", null, "ACG", null, 0, 0, 0, -1);
			run("absent-quality-control", 0, small(false), 3, false, absent, absent, true, false);
			final Expected nonoverlap=expect("ACG", q(12, 16, 20), "TTACG", q(3, 5, 24, 28, 32), 0, 7, 9, -1);
			run("nonoverlap-control", 0, small(true), 8, false, nonoverlap, nonoverlap, true, true);
			final boolean strict=run("strict-public-ecc", 1, raw(), 40, false,
				expect(FRAGMENT, repeat(40, 30), RAW_FRAGMENT, RAW_WRONG, 40, 0, 0, 40),
				expect(FRAGMENT, repeat(40, 30), RAW_FRAGMENT, RAW_CORRECT, 40, 0, 0, 40), false, false);
			final Expected loose=expect(FRAGMENT, repeat(40, 30), RAW_INPUT, RAW_QUALITY, 40, 7, 9, 40);
			final boolean permissive=run("loose-public-control", 2, raw(), 40, false, loose, loose, true, true);
			cli(strict && permissive && failures==0);
		}finally{
			BBMerge.changeQuality=oldQuality;
			BBMerge.staticAdapterList=oldAdapters;
		}
		check(BBMerge.changeQuality==oldQuality && BBMerge.staticAdapterList==oldAdapters, "Owned public settings were not restored.");
	}

	/** Invokes actual public code; failed wrapper prerequisites remain visible and gate the CLI. */
	private boolean run(final String name, final int api, final Pair pair, final int insert, final boolean change,
			final Expected source, final Expected corrected, final boolean control, final boolean keepArrays) throws IOException{
		cases++; if(api==0){direct++;}else{wrappers++;}
		BBMerge.changeQuality=change;
		final Snapshot before=new Snapshot(pair);
		append(out.resolve("inputs.tsv"), name+"\t"+api+"\t"+insert+"\t"+change+"\t"+before.text+"\n");
		append(out.resolve("oracles.tsv"), name+"\t"+source.text()+"\t"+corrected.text()+"\n");
		int returned=Integer.MIN_VALUE; boolean completed=false;
		final int priorFailures=failures;
		try{
			returned=api==0 ? BBMerge.errorCorrectWithInsert(pair.a, pair.b, insert) :
				api==1 ? BBMerge.findOverlapStrict(pair.a, pair.b, true) : BBMerge.findOverlapLoose(pair.a, pair.b, false);
			completed=true;
		}catch(RuntimeException|AssertionError e){e.printStackTrace(System.err);}
		check(completed, name+": public call failed; dependent CLI must not run.");
		final Expected wanted=baseline ? source : corrected;
		check(returned==wanted.returned, name+": return must match the literal candidate/correction count.");
		check(Arrays.equals(pair.a.bases, wanted.ab), name+": R1 bases differ.");
		check(Arrays.equals(pair.b.bases, wanted.bb), name+": R2 bases differ.");
		check(Arrays.equals(pair.a.quality, wanted.aq), name+": R1 qualities differ.");
		check(Arrays.equals(pair.b.quality, wanted.bq), name+": R2 qualities differ.");
		check(pair.a.errors==wanted.ae && pair.b.errors==wanted.be, name+": error counters differ.");
		check(pair.a.insert()==wanted.insert && pair.b.insert()==wanted.insert, name+": stored inserts differ.");
		check(pair.a.strand()==before.as && pair.b.strand()==before.bs, name+": original strand was not restored.");
		check(pair.a.mate==pair.b && pair.b.mate==pair.a && pair.a.id.equals(ID) && pair.b.id.equals(ID) && pair.a.numericID==0 && pair.b.numericID==0,
			name+": read/mate identity changed.");
		if(keepArrays){check(before.sameArrays(pair), name+": no-correction control replaced or changed arrays.");}
		final boolean sourceMatch=completed && source.matches(pair, returned), correctedMatch=completed && corrected.matches(pair, returned);
		if(!correctedMatch){discrepancies++;}
		final String row=name+"\t"+returned+"\t"+packet(pair)+"\t"+sourceMatch+"\t"+correctedMatch+"\t"+completed+"\n";
		append(out.resolve("cases.tsv"), row); if(control){append(out.resolve("controls.tsv"), row);}
		return priorFailures==failures;
	}

	/** Exactly one real current BBDuk CLI, after both actual public-wrapper prerequisites. */
	private void cli(final boolean ready) throws Exception{
		final Pair input=raw();
		writeFastq(out.resolve("input1.fq"), input.a); writeFastq(out.resolve("input2.fq"), input.b);
		final byte[] first=Files.readAllBytes(out.resolve("input1.fq")), second=Files.readAllBytes(out.resolve("input2.fq"));
		append(out.resolve("inputs.tsv"), "bbduk-cli\tCLI\t40\tfalse\t"+packet(input)+"\n");
		final Expected source=expect(FRAGMENT, repeat(40, 30), RAW_FRAGMENT, RAW_WRONG, 0, 0, 0, -1);
		final Expected corrected=expect(FRAGMENT, repeat(40, 30), RAW_FRAGMENT, RAW_CORRECT, 0, 0, 0, -1);
		append(out.resolve("oracles.tsv"), "bbduk-cli\t"+source.text()+"\t"+corrected.text()+"\n");
		if(!ready){
			write(out.resolve("cli.status"), "SKIPPED: public-wrapper prerequisites failed\n");
			append(out.resolve("cases.tsv"), "bbduk-cli\tSKIPPED\tNA\tfalse\tfalse\tfalse\n");
			return;
		}
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbduk.sh").toString(), "-Xmx128m", "-Xms32m",
			"in1="+out.resolve("input1.fq"), "in2="+out.resolve("input2.fq"), "out1="+out.resolve("output1.fq"), "out2="+out.resolve("output2.fq"),
			"ecc=t", "cq=f", "qin=33", "qout=33", "ordered=t", "t=1", "simd=f", "ow=f", "qtrim=f", "ktrim=f", "tbo=f", "tpe=f",
			"ftl=0", "ftr=0", "ftr2=0", "ftm=0", "maq=0", "minlength=0", "minlenfraction=0", "maxlength=2147483647"));
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(out.resolve("cli.command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(out.resolve("cli.stdout").toFile()); builder.redirectError(out.resolve("cli.stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process child=builder.start(); children++; cli++; cases++;
		final int code=child.waitFor(); write(out.resolve("cli.exit"), code+"\n"); write(out.resolve("cli.status"), "DISPATCHED\n");
		check(code==0, "Actual BBDuk child failed; inspect saved stderr.");
		check(Arrays.equals(first, Files.readAllBytes(out.resolve("input1.fq"))) && Arrays.equals(second, Files.readAllBytes(out.resolve("input2.fq"))), "CLI changed an input FASTQ.");
		if(code!=0){append(out.resolve("cases.tsv"), "bbduk-cli\t"+code+"\tNA\tfalse\tfalse\tfalse\n"); return;}
		final List<Read> a=nativeReads(out.resolve("output1.fq")), b=nativeReads(out.resolve("output2.fq"));
		check(a.size()==1 && b.size()==1, "Actual CLI must emit exactly one read per mate file.");
		if(a.size()!=1 || b.size()!=1){return;}
		final Pair observed=new Pair(a.get(0), b.get(0));
		final Expected wanted=baseline ? source : corrected;
		check(a.get(0).id.equals(ID) && b.get(0).id.equals(ID), "CLI changed literal read IDs.");
		check(wanted.arrays(observed), "Native FASTQ records differ from phase-local literal base/quality oracles.");
		final boolean sourceMatch=source.arrays(observed), correctedMatch=corrected.arrays(observed);
		if(!correctedMatch){discrepancies++;}
		write(out.resolve("cli.outputs.tsv"), packet(observed)+"\n");
		append(out.resolve("cases.tsv"), "bbduk-cli\t"+code+"\t"+packet(observed)+"\t"+sourceMatch+"\t"+correctedMatch+"\ttrue\n");
	}

	/** Uses the native FASTQ reader; its process-local decoding settings are restored. */
	private static List<Read> nativeReads(final Path file){
		require(Files.isRegularFile(file), "Missing CLI output: "+file);
		final byte offset=FASTQ.ASCII_OFFSET;
		final boolean detect=FASTQ.DETECT_QUALITY, test=FASTQ.TEST_INTERLEAVED, force=FASTQ.FORCE_INTERLEAVED;
		try{
			FASTQ.ASCII_OFFSET=33; FASTQ.DETECT_QUALITY=false; FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
			return StreamerFactory.getReads(-1, false, FileFormat.testInput(file.toString(), FileFormat.FASTQ, null, true, false), null, null, null);
		}finally{FASTQ.ASCII_OFFSET=offset; FASTQ.DETECT_QUALITY=detect; FASTQ.TEST_INTERLEAVED=test; FASTQ.FORCE_INTERLEAVED=force;}
	}

	private static Pair small(final boolean quality){return pair("ACG", quality ? q(12, 16, 20) : null, "TTACG", quality ? q(3, 5, 24, 28, 32) : null);}
	private static Pair raw(){return pair(FRAGMENT, repeat(40, 30), RAW_INPUT, RAW_QUALITY);}
	private static Pair pair(final String a, final byte[] aq, final String b, final byte[] bq){
		final Read first=new Read(bytes(a), copy(aq), ID, 0, false), second=new Read(bytes(b), copy(bq), ID, 0, false);
		first.mate=second; second.mate=first; second.setPairnum(1); first.errors=7; second.errors=9;
		return new Pair(first, second);
	}
	private static Expected expect(final String a, final byte[] aq, final String b, final byte[] bq, final int ret, final int ae, final int be, final int insert){return new Expected(bytes(a), aq, bytes(b), bq, ret, ae, be, insert);}
	private static String packet(final Pair p){return p.a.id+":"+p.a.length()+":"+new String(p.a.bases, StandardCharsets.US_ASCII)+":"+Arrays.toString(p.a.quality)+":"+p.a.errors+":"+p.a.insert()+"|"+p.b.id+":"+p.b.length()+":"+new String(p.b.bases, StandardCharsets.US_ASCII)+":"+Arrays.toString(p.b.quality)+":"+p.b.errors+":"+p.b.insert();}
	private static void writeFastq(final Path file, final Read r) throws IOException{
		final StringBuilder text=new StringBuilder("@").append(r.id).append('\n').append(new String(r.bases, StandardCharsets.US_ASCII)).append("\n+\n");
		for(final byte q : r.quality){text.append((char)(q+33));} write(file, text.append('\n').toString());
	}
	private static byte[] repeat(final int n, final int value){final byte[] a=new byte[n]; Arrays.fill(a, (byte)value); return a;}
	private static byte[] q(final int... scores){final byte[] a=new byte[scores.length]; for(int i=0; i<a.length; i++){a[i]=(byte)scores[i];} return a;}
	private static byte[] copy(final byte[] a){return a==null ? null : a.clone();}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String s) throws IOException{Files.write(file, bytes(s), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String s) throws IOException{Files.write(file, bytes(s), StandardOpenOption.APPEND);}
	private void same(final Path a, final Path b, final String why) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), why);}
	private void check(final boolean condition, final String why){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Pair{
		Pair(final Read a_, final Read b_){a=a_; b=b_;}
		final Read a, b;
	}
	private static final class Snapshot{
		Snapshot(final Pair p){ab=p.a.bases; bb=p.b.bases; aq=p.a.quality; bq=p.b.quality; ac=copy(ab); bc=copy(bb); aqa=copy(aq); bqa=copy(bq); as=p.a.strand(); bs=p.b.strand(); text=packet(p);}
		boolean sameArrays(final Pair p){return p.a.bases==ab && p.b.bases==bb && p.a.quality==aq && p.b.quality==bq && Arrays.equals(ab, ac) && Arrays.equals(bb, bc) && Arrays.equals(aq, aqa) && Arrays.equals(bq, bqa);}
		final byte[] ab, bb, aq, bq, ac, bc, aqa, bqa;
		final int as, bs;
		final String text;
	}
	private static final class Expected{
		Expected(final byte[] ab_, final byte[] aq_, final byte[] bb_, final byte[] bq_, final int ret, final int ae_, final int be_, final int insert_){ab=ab_; aq=aq_; bb=bb_; bq=bq_; returned=ret; ae=ae_; be=be_; insert=insert_;}
		boolean arrays(final Pair p){return Arrays.equals(p.a.bases, ab) && Arrays.equals(p.b.bases, bb) && Arrays.equals(p.a.quality, aq) && Arrays.equals(p.b.quality, bq);}
		boolean matches(final Pair p, final int ret){return ret==returned && arrays(p) && p.a.errors==ae && p.b.errors==be && p.a.insert()==insert && p.b.insert()==insert;}
		String text(){return returned+":"+new String(ab, StandardCharsets.US_ASCII)+":"+Arrays.toString(aq)+":"+new String(bb, StandardCharsets.US_ASCII)+":"+Arrays.toString(bq)+":"+ae+":"+be+":"+insert;}
		final byte[] ab, aq, bb, bq;
		final int returned, ae, be, insert;
	}

	private final Path root, out, reference;
	private final boolean baseline;
	private int cases, direct, wrappers, cli, children, discrepancies, checks, failures;
	private static final String ID="retained-quality-pair";
	private static final String FRAGMENT="ACGTTGCATGTCAGTACGATCGTACCTGAGCTAGTCGATC";
	private static final String RAW_FRAGMENT="GATCGACTAGCTCAGGTACGATCGTACTGACATGCAACGT";
	private static final String RAW_INPUT="GATCGACTAGCTCAGGTACGATCGTACTGACATGCAACGTAGATCGGAAGAGCACA";
	private static final byte[] RAW_QUALITY=q(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
		20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35);
	private static final byte[] RAW_CORRECT=q(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
		20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39);
	private static final byte[] RAW_WRONG=q(36, 37, 38, 39, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
		36, 37, 38, 39, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35);
}
