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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fileIO.FileFormat;
import parse.Parse;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;

/** Two fresh legacy BBDuk children isolate the mate selected by filterpolyc.
 * Four literal pairs run under either-mate-bad and both-mates-bad policies.
 * Baseline requires the two wrong-mate disposition witnesses and six controls;
 * fixed requires the independent per-mate oracle, with exact unchanged payloads.
 * No reflection, private worker calls, new production API or timing is involved.
 * @author Sucrose
 */
public final class BBDukPolyCFilterTest{

	/** phase=baseline|fixed, root=RUNTIME, outdir=NEW; fixed requires reference=BASELINE. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("phase", "root", "outdir", "reference").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukpolycfiltertest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed requires the accepted baseline reference.");
		require(Files.isRegularFile(root.resolve("bbdukOld.sh")), "Set root to the single BBTools runtime containing bbdukOld.sh.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBDUK_POLYC_FILTER_PASS phase=baseline version=1 policies=2 children=2 pairs=8 controls=6 discrepancies=2 sourceMatches=8 correctedMatches=6 checks=[0-9]+ failures=0\\n"), "Reference must be the successful two-policy poly-C baseline.");
		}
		final BBDukPolyCFilterTest test=new BBDukPolyCFilterTest(root, out, reference, phase);
		try{test.exercise();}
		finally{test.saveDispatch();}
		test.finish();
	}

	private BBDukPolyCFilterTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws IOException{
		root=root_; out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("input1.fq"), input(1)); write(out.resolve("input2.fq"), input(2));
		final StringBuilder inputs=new StringBuilder("pair\tR1\tR2\tquality\n");
		final StringBuilder oracles=new StringBuilder("policy\tpair\tbaseline_kept\tfixed_kept\twitness\n");
		final StringBuilder settings=new StringBuilder("policy\targument\n");
		for(final Pair pair : PAIRS){inputs.append(pair.name).append('\t').append(pair.first).append('\t').append(pair.second).append("\tQ40\n");}
		for(int policy=0; policy<2; policy++){
			for(int p=0; p<PAIRS.length; p++){
				oracles.append(POLICIES[policy]).append('\t').append(PAIRS[p].name).append('\t').append(SOURCE_KEEP[policy][p]).append('\t').append(FIXED_KEEP[policy][p]).append('\t').append(witness(policy, p)).append('\n');
			}
			for(final String arg : COMMON){settings.append(POLICIES[policy]).append('\t').append(arg).append('\n');}
			settings.append(POLICIES[policy]).append("\trequirebothbad=").append(policy==1 ? 't' : 'f').append('\n');
		}
		write(out.resolve("inputs.tsv"), inputs.toString()); write(out.resolve("oracles.tsv"), oracles.toString()); write(out.resolve("settings.tsv"), settings.toString());
		write(out.resolve("cases.tsv"), "policy\tpair\tkept\tsource_kept\tfixed_kept\twitness\n");
		write(out.resolve("controls.tsv"), "policy\tpair\tkept\tR1_id\tR2_id\tR1_bases\tR2_bases\tR1_quality\tR2_quality\n");
	}

	/** Unexpected child behavior throws before the next public child can start. */
	private void exercise() throws Exception{
		for(int policy=0; policy<2; policy++){run(policy);}
		check(children==2 && pairs==8 && controls==6, "Keep two policies, eight pair dispositions and six unchanged controls.");
		check(discrepancies==(baseline ? 2 : 0), "Only the two baseline wrong-mate dispositions may disagree with the fixed oracle.");
		check(sourceMatches==(baseline ? 8 : 6) && correctedMatches==(baseline ? 6 : 8), "Source/fixed oracle counts must distinguish both mate orientations.");
		if(reference!=null){
			for(final String file : new String[]{"input1.fq", "input2.fq", "inputs.tsv", "oracles.tsv", "settings.tsv", "controls.tsv"}){same(out.resolve(file), reference.resolve(file));}
			for(final String policy : POLICIES){
				for(final String oracle : new String[]{"source", "fixed"}){
					for(final String stream : STREAMS){same(out.resolve(policy).resolve(oracle+"-"+stream), reference.resolve(policy).resolve(oracle+"-"+stream));}
				}
			}
		}
	}

	/** Checks all four actual output files against literal whole-file and native-read oracles. */
	private void run(final int policy) throws Exception{
		final Path work=out.resolve(POLICIES[policy]); Files.createDirectory(work);
		for(int stream=0; stream<STREAMS.length; stream++){
			write(work.resolve("source-"+STREAMS[stream]), expected(SOURCE_KEEP[policy], stream));
			write(work.resolve("fixed-"+STREAMS[stream]), expected(FIXED_KEEP[policy], stream));
		}
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbdukOld.sh").toString(), "-Xmx128m", "-Xms32m", "-ea"));
		command.addAll(Arrays.asList(COMMON));
		command.addAll(Arrays.asList("requirebothbad="+(policy==1 ? "t" : "f"), "in1="+out.resolve("input1.fq"), "in2="+out.resolve("input2.fq"),
			"out="+work.resolve("kept1.fq"), "out2="+work.resolve("kept2.fq"), "outm="+work.resolve("rejected1.fq"), "outm2="+work.resolve("rejected2.fq")));
		final StringBuilder argv=new StringBuilder("index\targument\n");
		for(int i=0; i<command.size(); i++){argv.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve("command.tsv"), argv.toString());
		write(work.resolve("environment.tsv"), "variable\tvalue\nJAVA_TOOL_OPTIONS\t-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8\n_JAVA_OPTIONS\tUNSET\nJDK_JAVA_OPTIONS\tUNSET\n");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile())
			.redirectOutput(work.resolve("child.stdout").toFile()).redirectError(work.resolve("child.stderr").toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		final Process process=builder.start(); children++; statuses[policy]="STARTED";
		process.getOutputStream().close();
		final int exit=process.waitFor();
		write(work.resolve("child.exit"), exit+"\n"); statuses[policy]="RETURNED";
		check(read(out.resolve("input1.fq")).equals(input(1)) && read(out.resolve("input2.fq")).equals(input(2)), "Child changed the literal paired inputs.");
		check(exit==0, POLICIES[policy]+": unexpected child exit "+exit+"; later child gated.");
		final String stderr=read(work.resolve("child.stderr"));
		check(Files.size(work.resolve("child.stdout"))==0, "Legacy tool must write reads only to the selected files.");
		check(!stderr.contains("Exception") && !stderr.contains("AssertionError") && !stderr.contains("Caused by:"), "Unexpected child failure diagnostic: "+stderr);
		final boolean[] keep=baseline ? SOURCE_KEEP[policy] : FIXED_KEEP[policy];
		final ArrayList<Read> kept1=verifyFile(work, 0, keep), kept2=verifyFile(work, 1, keep);
		final ArrayList<Read> rejected1=verifyFile(work, 2, keep), rejected2=verifyFile(work, 3, keep);
		check(kept1.size()+rejected1.size()==4 && kept2.size()+rejected2.size()==4, "Each mate must occur exactly once in kept or rejected output.");
		counter(stderr, INPUT, 8, 96, "input");
		counter(stderr, RESULT, kept1.size()*2L, kept1.size()*24L, "kept result");
		int good=0, bad=0;
		for(int p=0; p<PAIRS.length; p++){
			final Read a=keep[p] ? kept1.get(good) : rejected1.get(bad);
			final Read b=keep[p] ? kept2.get(good++) : rejected2.get(bad++);
			final boolean observed=contains(kept1, PAIRS[p].name+"/1");
			check(observed==keep[p], "Observed pair disposition differs from phase oracle: "+PAIRS[p].name);
			if(observed==SOURCE_KEEP[policy][p]){sourceMatches++;}
			if(observed==FIXED_KEEP[policy][p]){correctedMatches++;}else{discrepancies++;}
			final String row=POLICIES[policy]+"\t"+PAIRS[p].name+"\t"+observed;
			append(out.resolve("cases.tsv"), row+"\t"+SOURCE_KEEP[policy][p]+"\t"+FIXED_KEEP[policy][p]+"\t"+witness(policy, p)+"\n");
			if(!witness(policy, p)){
				controls++;
				append(out.resolve("controls.tsv"), row+"\t"+a.id+"\t"+b.id+"\t"+ascii(a.bases)+"\t"+ascii(b.bases)+"\t"+quality(a.quality)+"\t"+quality(b.quality)+"\n");
			}
			pairs++;
		}
		statuses[policy]="PASSED";
	}

	private ArrayList<Read> verifyFile(final Path work, final int stream, final boolean[] keep) throws IOException{
		final Path file=work.resolve(STREAMS[stream]);
		check(Files.isRegularFile(file), "Missing paired output: "+file);
		check(read(file).equals(expected(keep, stream)), "Whole-file ID/order/bases/quality oracle mismatch: "+file);
		final ArrayList<Read> reads=readFastq(file);
		int index=0;
		for(int p=0; p<PAIRS.length; p++){
			if(keep[p]!=(stream<2)){continue;}
			check(index<reads.size(), "Native reader lost expected record: "+file);
			final Read r=reads.get(index++);
			final int mate=(stream&1)+1;
			final String bases=mate==1 ? PAIRS[p].first : PAIRS[p].second;
			check(r.id.equals(PAIRS[p].name+"/"+mate) && r.mate==null, "Native reader ID/pairing mismatch: "+file);
			check(Arrays.equals(r.bases, bytes(bases)), "Native reader bases mismatch: "+file);
			check(r.quality!=null && Arrays.equals(r.quality, QUALITIES), "Native reader must decode exactly twelve Q40 qualities: "+file);
		}
		check(index==reads.size(), "Unexpected duplicate or extra native FASTQ records: "+file);
		return reads;
	}

	/** Native decoding borrows only these globals and restores every borrowed value. */
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

	private void counter(final String stderr, final Pattern pattern, final long reads, final long bases, final String name){
		final Matcher match=pattern.matcher(stderr);
		check(match.find(), "Missing final "+name+" counter.");
		check(Long.parseLong(match.group(1))==reads && Long.parseLong(match.group(2))==bases, "Unexpected "+name+" read/base counts.");
		check(!match.find(), "Duplicate "+name+" counter.");
	}
	private static boolean contains(final ArrayList<Read> reads, final String id){for(final Read r : reads){if(r.id.equals(id)){return true;}}return false;}
	private static String expected(final boolean[] keep, final int stream){
		final StringBuilder text=new StringBuilder();
		for(int p=0; p<PAIRS.length; p++){if(keep[p]==(stream<2)){text.append(record(p, (stream&1)+1));}}
		return text.toString();
	}
	private static String input(final int mate){final StringBuilder text=new StringBuilder();for(int p=0; p<PAIRS.length; p++){text.append(record(p, mate));}return text.toString();}
	private static String record(final int p, final int mate){return "@"+PAIRS[p].name+"/"+mate+"\n"+(mate==1 ? PAIRS[p].first : PAIRS[p].second)+"\n+\nIIIIIIIIIIII\n";}
	private static String quality(final byte[] values){final StringBuilder text=new StringBuilder();for(final byte q : values){text.append((char)(q+33));}return text.toString();}
	private static boolean witness(final int policy, final int pair){return SOURCE_KEEP[policy][pair]!=FIXED_KEEP[policy][pair];}
	private void saveDispatch() throws IOException{final StringBuilder text=new StringBuilder("policy\tstatus\n");for(int i=0; i<POLICIES.length; i++){text.append(POLICIES[i]).append('\t').append(statuses[i]).append('\n');}write(out.resolve("dispatch.tsv"), text.toString());}
	private void finish() throws IOException{
		final String summary="BBDUK_POLYC_FILTER_PASS phase="+phase+" version=1 policies=2 children="+children+" pairs="+pairs+" controls="+controls+" discrepancies="+discrepancies+
			" sourceMatches="+sourceMatches+" correctedMatches="+correctedMatches+" checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}
	private void same(final Path a, final Path b) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), "Baseline literal/control changed: "+a);}
	private void check(final boolean condition, final String why){checks++;if(!condition){throw new IllegalStateException(why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static Path path(final String text){require(text!=null && !text.isEmpty(), "Missing required path.");return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String ascii(final byte[] values){return new String(values, StandardCharsets.US_ASCII);}
	private static String read(final Path file) throws IOException{return ascii(Files.readAllBytes(file));}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}

	private static final class Pair{
		Pair(final String name_, final String first_, final String second_){name=name_;first=first_;second=second_;require(first.length()==12 && second.length()==12, "Literal pairs must remain twelve bases per mate.");}
		final String name, first, second;
	}
	private final Path root, out, reference;
	private final String phase;
	private final boolean baseline;
	private int children, pairs, controls, discrepancies, sourceMatches, correctedMatches, checks;
	private final String[] statuses={"SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE"};
	private static final String[] POLICIES={"either-bad", "both-bad"};
	private static final Pair[] PAIRS={new Pair("none_none", "ATGATGATGATG", "ATGATGATGATG"), new Pair("poly_none", "CCCCCCATGATG", "ATGATGATGATG"),
		new Pair("none_poly", "ATGATGATGATG", "CCCCCCATGATG"), new Pair("poly_poly", "CCCCCCATGATG", "CCCCCCATGATG")};
	private static final boolean[][] SOURCE_KEEP={{true, false, true, false}, {true, false, true, false}};
	private static final boolean[][] FIXED_KEEP={{true, false, false, false}, {true, true, true, false}};
	private static final String[] STREAMS={"kept1.fq", "kept2.fq", "rejected1.fq", "rejected2.fq"};
	private static final byte[] QUALITIES={40,40,40,40,40,40,40,40,40,40,40,40};
	private static final String[] COMMON={"t=1", "simd=f", "bf1=t", "ordered=t", "ow=f", "app=f", "qin=33", "qout=33", "interleaved=f",
		"filterpolyc=6", "maxnonpoly=0", "qtrim=f", "minlen=1", "maq=0", "mbq=0", "maxns=-1", "maxnrate=1", "minconsecutivebases=0", "minbasefrequency=0",
		"trimbyoverlap=f", "ecc=f", "entropy=-1", "entropymask=f", "entropytrim=f", "entropymark=f", "trimpolya=0",
		"trimpolygleft=0", "trimpolygright=0", "filterpolyg=0", "trimpolycleft=0", "trimpolycright=0", "trimfailures=f", "pairedtosingle=f", "quantize=f",
		"samplerate=1", "sampleseed=1", "reads=-1", "recalibrate=f", "rename=f", "json=f"};
	private static final Pattern INPUT=Pattern.compile("(?m)^Input:\\s+([0-9]+) reads\\s+([0-9]+) bases[.]\\s*$");
	private static final Pattern RESULT=Pattern.compile("(?m)^Result:\\s+([0-9]+) reads [^\\r\\n]*?\\s+([0-9]+) bases\\b[^\\r\\n]*$");
}
