package test.jgi;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import parse.Parse;

/** Equal-mate legacy BBDuk CLI timing support, separate from the correctness fixture.
 * Writes compact paired gzip inputs once, or validates ONE fresh public-launcher child.
 * Each 64-pair block has 32 non-poly-C pairs followed by 32 poly-C pairs, twelve bases
 * per mate and Q40 throughout. Block-local IDs deliberately repeat between blocks.
 * Either source operand gives the same dispositions on these equal-mate controls.
 * Bash owns calibration, the frozen input count, repetitions and phase comparisons.
 * GNU time measures only the ordinary launcher and its descendants, excluding this
 * fixture's generation and validation. Results include CLI startup/JIT/decompression;
 * they do not isolate the changed operand or establish biological accuracy.
 * @author Sucrose
 */
public final class BBDukPolyCTimingTest{

	/** mode=write outdir=NEW blocks=N; mode=run root=RUNTIME inputdir=INPUT outdir=NEW phase=pilot|baseline|fixed. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && !pair[1].isEmpty(), "Use documented nonempty flag=value arguments: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukpolyctimingtest.sh with assertions enabled.");
		final String mode=options.get("mode");
		if("write".equals(mode)){
			keys(options, "mode", "outdir", "blocks");
			writeInputs(path(options.get("outdir")), positive(options.get("blocks")));
		}else{
			require("run".equals(mode), "Set mode=write|run.");
			keys(options, "mode", "root", "inputdir", "outdir", "phase");
			run(path(options.get("root")), path(options.get("inputdir")), path(options.get("outdir")), options.get("phase"));
		}
	}

	/** One immutable 64-pair block per mate is reused; no per-record allocation in the large write. */
	private static void writeInputs(final Path out, final long blocks) throws Exception{
		final long pairs=Math.multiplyExact(blocks, 64L);
		Math.multiplyExact(pairs, 24L);// Require representable read/base counter oracles before creating anything.
		final byte[] first=block(1), second=block(2);
		final long firstRaw=Math.multiplyExact(blocks, first.length), secondRaw=Math.multiplyExact(blocks, second.length);
		Files.createDirectory(out);
		final String firstHash=gzip(out.resolve("input1.fq.gz"), first, blocks);
		final String secondHash=gzip(out.resolve("input2.fq.gz"), second, blocks);
		final String spec="version\t1\nblock_pairs\t64\nblocks\t"+blocks+"\npairs\t"+pairs+
			"\ninput1_raw_bytes\t"+firstRaw+"\ninput2_raw_bytes\t"+secondRaw+
			"\ninput1_gzip_bytes\t"+Files.size(out.resolve("input1.fq.gz"))+"\ninput2_gzip_bytes\t"+Files.size(out.resolve("input2.fq.gz"))+
			"\ninput1_sha80\t"+firstHash+"\ninput2_sha80\t"+secondHash+"\n";
		write(out.resolve("input-spec.tsv"), spec);
		System.out.println("BBDUK_POLYC_TIMING_INPUT version=1 blocks="+blocks+" pairs="+pairs+" failures=0");
	}

	/** Checks metadata each sample; Bash must hash both complete gzip files before/after each phase.
	 * No tool output files are requested: counters validate equal-mate workload completion,
	 * while the separate four-pair correctness fixture owns exact output-byte validation.
	 */
	private static void run(final Path root, final Path input, final Path out, final String phase) throws Exception{
		require(Arrays.asList("pilot", "baseline", "fixed").contains(phase), "Set phase=pilot|baseline|fixed.");
		require(Files.isRegularFile(root.resolve("bbdukOld.sh")), "Missing ordinary legacy launcher.");
		require(Files.isExecutable(Paths.get("/usr/bin/time")), "GNU /usr/bin/time is required.");
		final String spec=read(input.resolve("input-spec.tsv"));
		final Map<String, String> fields=fields(spec);
		keys(fields, "version", "block_pairs", "blocks", "pairs", "input1_raw_bytes", "input2_raw_bytes", "input1_gzip_bytes", "input2_gzip_bytes", "input1_sha80", "input2_sha80");
		require("1".equals(fields.get("version")) && "64".equals(fields.get("block_pairs")), "Wrong timing-input version or block composition.");
		final long blocks=positive(fields.get("blocks")), pairs=positive(fields.get("pairs"));
		require(pairs==Math.multiplyExact(blocks, 64L), "Pair count must be a positive multiple of the literal 64-pair block.");
		final long bases=Math.multiplyExact(pairs, 24L), reads=Math.multiplyExact(pairs, 2L);
		for(int mate=1; mate<=2; mate++){
			final String prefix="input"+mate;
			require(positive(fields.get(prefix+"_raw_bytes"))==Math.multiplyExact(blocks, block(mate).length), "Unexpected uncompressed byte count: "+prefix);
			require(fields.get(prefix+"_sha80").matches("[0-9a-f]{20}"), "Missing compressed input fingerprint: "+prefix);
			require(Files.size(input.resolve(prefix+".fq.gz"))==positive(fields.get(prefix+"_gzip_bytes")), "Compressed input size changed: "+prefix);
		}
		Files.createDirectory(out); write(out.resolve("input-spec.tsv"), spec);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/time", "-f", "wall\t%e\nuser\t%U\nsystem\t%S", "-o", out.resolve("child.time").toString(),
			"bash", root.resolve("bbdukOld.sh").toString(), "-Xmx128m", "-Xms32m", "-ea"));
		command.addAll(Arrays.asList(COMMON));
		command.add("in1="+input.resolve("input1.fq.gz")); command.add("in2="+input.resolve("input2.fq.gz"));
		final StringBuilder argv=new StringBuilder("index\targument\n");
		for(int i=0; i<command.size(); i++){argv.append(i).append('\t').append(command.get(i).replace("\t", "\\t").replace("\n", "\\n")).append('\n');}
		write(out.resolve("command.tsv"), argv.toString());
		write(out.resolve("environment.tsv"), "variable\tvalue\nLC_ALL\tC\nJAVA_TOOL_OPTIONS\t-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8\n_JAVA_OPTIONS\tUNSET\nJDK_JAVA_OPTIONS\tUNSET\n");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile())
			.redirectOutput(out.resolve("child.stdout").toFile()).redirectError(out.resolve("child.stderr").toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		builder.environment().put("LC_ALL", "C");
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		final Process process=builder.start(); process.getOutputStream().close();
		final int exit=process.waitFor(); write(out.resolve("child.exit"), exit+"\n");
		require(exit==0, "Legacy timing child failed: "+exit+"; inspect child.stderr before another sample.");
		require(Files.size(out.resolve("child.stdout"))==0, "Timing child must not emit sequence output.");
		final String stderr=read(out.resolve("child.stderr"));
		require(!stderr.contains("Exception") && !stderr.contains("AssertionError") && !stderr.contains("Caused by:"), "Unexpected timing-child diagnostic.");
		counter(stderr, "Input", reads, bases, true);
		counter(stderr, "Polymer-trimmed", pairs, bases/2, false);
		counter(stderr, "Total Removed", pairs, bases/2, false);
		counter(stderr, "Result", pairs, bases/2, false);
		require(spec.equals(read(input.resolve("input-spec.tsv"))), "Input metadata changed during the sample.");
		final Map<String, String> times=fields(read(out.resolve("child.time")));
		keys(times, "wall", "user", "system");
		final double wall=seconds(times.get("wall")), user=seconds(times.get("user")), system=seconds(times.get("system")), cpu=user+system;
		require("pilot".equals(phase) || cpu>=10, "Measured sample below ten user+system CPU seconds; enlarge the common input and redo BOTH phases.");
		final String summary="BBDUK_POLYC_TIMING_PASS version=1 phase="+phase+" blocks="+blocks+" pairs="+pairs+" inputReads="+reads+" inputBases="+bases+
			" keptReads="+pairs+" keptBases="+(bases/2)+" wall="+wall+" user="+user+" system="+system+" cpu="+cpu+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Java GZIP output is hashed as compressed bytes are written; only its final 80 bits are emitted. */
	private static String gzip(final Path file, final byte[] block, final long blocks) throws Exception{
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		try(final OutputStream raw=Files.newOutputStream(file, StandardOpenOption.CREATE_NEW);
			final GZIPOutputStream zip=new GZIPOutputStream(new DigestOutputStream(new BufferedOutputStream(raw, 65536), digest), 65536)){
			for(long i=0; i<blocks; i++){zip.write(block);}
		}
		final byte[] hash=digest.digest(); final StringBuilder result=new StringBuilder(20);
		for(int i=hash.length-10; i<hash.length; i++){result.append(Character.forDigit((hash[i]>>>4)&15, 16)).append(Character.forDigit(hash[i]&15, 16));}
		return result.toString();
	}

	private static byte[] block(final int mate){
		require(mate==1 || mate==2, "The literal workload has exactly two mates.");
		final StringBuilder text=new StringBuilder();
		for(int i=0; i<64; i++){
			text.append('@').append(i<32 ? "nn" : "cc").append(i).append('/').append(mate).append('\n')
				.append(i<32 ? "ATGATGATGATG" : "CCCCCCATGATG").append("\n+\nIIIIIIIIIIII\n");
		}
		return text.toString().getBytes(StandardCharsets.US_ASCII);
	}
	private static void counter(final String text, final String label, final long reads, final long bases, final boolean input){
		final String suffix=input ? " reads\\s+([0-9]+) bases[.]\\s*$" : " reads [^\\r\\n]*?\\s+([0-9]+) bases\\b[^\\r\\n]*$";
		final Matcher match=Pattern.compile("(?m)^"+Pattern.quote(label)+":\\s+([0-9]+)"+suffix).matcher(text);
		require(match.find(), "Missing "+label+" counters.");
		require(Long.parseLong(match.group(1))==reads && Long.parseLong(match.group(2))==bases, "Unexpected "+label+" read/base counts.");
		require(!match.find(), "Duplicate "+label+" counters.");
	}
	private static Map<String, String> fields(final String text){
		final Map<String, String> result=new LinkedHashMap<String, String>();
		for(final String line : text.split("\n")){
			final String[] pair=line.split("\t", -1);
			require(pair.length==2 && !pair[1].isEmpty() && result.put(pair[0], pair[1])==null, "Malformed or duplicate metadata/time row: "+line);
		}
		return result;
	}
	private static void keys(final Map<String, String> options, final String... names){require(options.size()==names.length && options.keySet().containsAll(Arrays.asList(names)), "Expected exactly these options/fields: "+Arrays.toString(names));}
	private static long positive(final String text){require(text!=null && text.matches("[0-9]+"), "Expected positive integer: "+text); final long value=Long.parseLong(text); require(value>0, "Count must be positive."); return value;}
	private static double seconds(final String text){require(text!=null && text.matches("[0-9]+[.][0-9]+"), "Expected GNU time decimal seconds: "+text); final double value=Double.parseDouble(text); require(!Double.isInfinite(value) && value>=0, "Invalid time: "+text); return value;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\n')<0 && text.indexOf('\t')<0, "Legacy paths must be nonempty without equals, tab or newline."); return Paths.get(text).toAbsolutePath().normalize();}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final String[] COMMON={"t=1", "simd=f", "bf1=t", "ordered=t", "ow=f", "app=f", "qin=33", "qout=33", "interleaved=f",
		"requirebothbad=f", "filterpolyc=6", "maxnonpoly=0", "qtrim=f", "minlen=1", "maq=0", "mbq=0", "maxns=-1", "maxnrate=1", "minconsecutivebases=0", "minbasefrequency=0",
		"trimbyoverlap=f", "ecc=f", "entropy=-1", "entropymask=f", "entropytrim=f", "entropymark=f", "trimpolya=0",
		"trimpolygleft=0", "trimpolygright=0", "filterpolyg=0", "trimpolycleft=0", "trimpolycright=0", "trimfailures=f", "pairedtosingle=f", "quantize=f",
		"samplerate=1", "sampleseed=1", "reads=-1", "recalibrate=f", "rename=f", "json=f", "unpigz=f", "gunzip=f", "unbgzip=f", "nativebgzip=f"};
}
