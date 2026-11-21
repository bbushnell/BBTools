package test.processor;

import java.io.BufferedOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parse;
import parse.Parser;
import processor.ReformatProcessor;
import stream.Read;
import structures.ByteBuilder;

/** Literal name-allocation, scope and input-description-trimming controls.
 * Baseline uses only existing options; new selector cases run only after repair.
 * Driver children use fresh JVMs and real input files, with no timeout or signal.
 * A separate write mode prepares bounded performance input and independent baseline
 * output hashes; it never invokes the processor or times work in this fixture.
 * @author Sucrose
 */
public final class ReformatNameCollisionTest{

	/** Requires phase=baseline|fixed and a new outdir; driver mode delegates an actual CLI child. */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && args[0].startsWith("driver=")){
			driver(args);
			return;
		}
		for(final String arg : args){
			if(arg.startsWith("mode=")){writeBenchmark(args); return;}
		}
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		if(directory==null){throw new IllegalArgumentException("Set a new outdir.");}
		final boolean baseline=phase.equals("baseline");
		final Path root=Files.createDirectory(Paths.get(directory)).toAbsolutePath();
		try(PrintWriter observations=new PrintWriter(Files.newBufferedWriter(root.resolve("observations.tsv"), StandardCharsets.UTF_8))){
			observations.println("case\texpected\tobserved");
			directCases(baseline, observations);
			paired(baseline, observations);
			cloneAndSelection(baseline, observations);
			if(!baseline){selectorControls(observations);}
			cliCases(root, baseline);
			if(observations.checkError()){throw new AssertionError("Could not retain naming observations.");}
		}
		final int expected=baseline ? BASELINE_CHECKS : FIXED_CHECKS;
		check(checks==expected-1, "Frozen check count before final sentry: "+checks+" != "+(expected-1));
		System.out.println("REFORMAT_NAME_COLLISION_PASS phase="+phase+" checks="+checks+" children="+children);
	}

	/** Writes repeated literal sixteen-read blocks and hashes independent current-behavior output.
	 * Requires mode=write phase=baseline blocks=1..1048576 outdir=NEW. The runner times
	 * reformat.sh directly with trd=f; this mode generates no enabled-output file and
	 * does not run the functional matrix, any driver, or the production name allocator.
	 */
	private static void writeBenchmark(final String[] args) throws Exception{
		String mode=null, phase=null, directory=null;
		long blocks=-1;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("blocks")){blocks=Long.parseLong(pair[1]);}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"write".equals(mode) || !"baseline".equals(phase) || directory==null || blocks<1 || blocks>1048576){
			throw new IllegalArgumentException("Set mode=write phase=baseline blocks=1..1048576 outdir=NEW.");
		}
		//Literal key equivalence is independent of ObjectIntMap and the production allocation code.
		for(int i=0; i<BENCH_HEADERS.length; i++){
			for(int j=0; j<i; j++){
				if(BENCH_HEADERS[i].equals(BENCH_HEADERS[j])!=(BENCH_KEYS[i]==BENCH_KEYS[j])){
					throw new AssertionError("Baseline suffix histories must group exactly equal original headers: "+i+", "+j);
				}
			}
			for(int j=0; j<BENCH_HEADERS[i].length(); j++){
				if(BENCH_HEADERS[i].charAt(j)>127){throw new AssertionError("The literal benchmark oracle requires ASCII headers.");}
			}
		}
		final byte[] block=fastq(BENCH_HEADERS).getBytes(StandardCharsets.US_ASCII);
		if(block.length!=271){throw new AssertionError("The frozen sixteen-read block must contain271 bytes, found "+block.length);}
		final long reads=Math.multiplyExact(16L, blocks), bases=Math.multiplyExact(4L, reads);
		final long inputBytes=Math.multiplyExact((long)block.length, blocks);
		final Path root=Files.createDirectory(Paths.get(directory)).toAbsolutePath();
		final MessageDigest inputDigest=MessageDigest.getInstance("SHA-256");
		final MessageDigest enabledDigest=MessageDigest.getInstance("SHA-256");
		final int[] counts=new int[11];
		final ByteBuilder expected=new ByteBuilder(64);
		long enabledBytes=0;
		try(OutputStream input=new BufferedOutputStream(new DigestOutputStream(Files.newOutputStream(root.resolve("input.fq"),
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE), inputDigest))){
			for(long b=0; b<blocks; b++){
				input.write(block);
				for(int i=0; i<BENCH_HEADERS.length; i++){
					//Current primitive-map absence yields suffix0; subsequent occurrences increment per full header.
					expected.clear().append('@').append(BENCH_HEADERS[i]).append('_').append(counts[BENCH_KEYS[i]]++)
						.append("\nACGT\n+\nIIII\n");
					enabledDigest.update(expected.array, 0, expected.length());
					enabledBytes+=expected.length();
				}
			}
		}
		if(Files.size(root.resolve("input.fq"))!=inputBytes){throw new AssertionError("Generated FASTQ byte count differs from the literal block count.");}
		final String inputHash=sha80(inputDigest.digest()), enabledHash=sha80(enabledDigest.digest());
		final String spec="phase\tbaseline\nblocks\t"+blocks+"\nreads\t"+reads+"\nbases\t"+bases+
			"\ninput_bytes\t"+inputBytes+"\ninput_sha80\t"+inputHash+
			"\nuniquenames_false_bytes\t"+inputBytes+"\nuniquenames_false_sha80\t"+inputHash+
			"\nuniquenames_true_bytes\t"+enabledBytes+"\nuniquenames_true_sha80\t"+enabledHash+"\ntrd\tfalse\n";
		Files.write(root.resolve("input-spec.tsv"), spec.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
		System.out.println("REFORMAT_NAME_INPUT_PASS phase=baseline blocks="+blocks+" reads="+reads+
			" bases="+bases+" bytes="+inputBytes);
	}

	/** Formats only the final80 bits of SHA-256; never emits a full SHA-256 fingerprint. */
	private static String sha80(final byte[] digest){
		final StringBuilder text=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			text.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));
		}
		return text.toString();
	}

	/** Expected names are literal and independent of the proposed allocation mechanism. */
	private static void directCases(final boolean baseline, final PrintWriter out){
		names(out, "normalize-forward", baseline,
			new String[]{"a b", "a_b"}, new String[]{"a_b_0", "a_b_0"}, new String[]{"a_b", "a_b_2"}, "underscore=t");
		names(out, "normalize-reverse", baseline,
			new String[]{"a_b", "a b"}, new String[]{"a_b_0", "a_b_0"}, new String[]{"a_b", "a_b_2"}, "underscore=t");
		names(out, "normalize-natural", baseline,
			new String[]{"a b", "a_b", "a_b_2"}, new String[]{"a_b_0", "a_b_0", "a_b_2_0"},
			new String[]{"a_b", "a_b_2", "a_b_2_2"}, "underscore=t");
		names(out, "normalize-generated", baseline,
			new String[]{"a b", "a b", "a_b"}, new String[]{"a_b_0", "a_b_1", "a_b_0"},
			new String[]{"a_b", "a_b_2", "a_b_3"}, "underscore=t");
		names(out, "generated-first", baseline,
			new String[]{"x", "x", "x_2"}, new String[]{"x_0", "x_1", "x_2_0"}, new String[]{"x", "x_2", "x_2_2"});
		names(out, "natural-first", baseline,
			new String[]{"x_2", "x", "x"}, new String[]{"x_2_0", "x_0", "x_1"}, new String[]{"x_2", "x", "x_3"});
		names(out, "reserved-chain", baseline,
			new String[]{"x_3", "x", "x", "x", "x_2"}, new String[]{"x_3_0", "x_0", "x_1", "x_2", "x_2_0"},
			new String[]{"x_3", "x", "x_2", "x_4", "x_2_2"});
		names(out, "noncanonical-chain", baseline,
			new String[]{"x", "x", "x_02", "x_1", "x_2_2", "x_2"},
			new String[]{"x_0", "x_1", "x_02_0", "x_1_0", "x_2_2_0", "x_2_0"},
			new String[]{"x", "x_2", "x_02", "x_1", "x_2_2", "x_2_3"});
		names(out, "nested-generated", baseline,
			new String[]{"x", "x", "x_2", "x_2", "x_2_2", "x"},
			new String[]{"x_0", "x_1", "x_2_0", "x_2_1", "x_2_2_0", "x_2"},
			new String[]{"x", "x_2", "x_2_2", "x_2_3", "x_2_2_2", "x_3"});
		names(out, "description-default", baseline,
			new String[]{"x a", "x\tb", "x a"}, new String[]{"x a_0", "x\tb_0", "x a_1"},
			new String[]{"x a", "x_2\tb", "x_3 a"});
		names(out, "description-natural", baseline,
			new String[]{"x a", "x_2 b", "x c"}, new String[]{"x a_0", "x_2 b_0", "x c_0"},
			new String[]{"x a", "x_2 b", "x_3 c"});
		names(out, "existing-slash", baseline,
			new String[]{"p/1", "p/1", "p_2/1"}, new String[]{"p/1_0", "p/1_1", "p_2/1_0"},
			new String[]{"p/1", "p_2/1", "p_2_2/1"});
		names(out, "spaced-slash", baseline,
			new String[]{"p", "p", "p_2"}, new String[]{"p_0 /1", "p_1 /1", "p_2_0 /1"},
			new String[]{"p /1", "p_2 /1", "p_2_2 /1"}, "addslash=t", "slashspace=t");
		names(out, "bare-slash", baseline,
			new String[]{"p", "p", "p_2"}, new String[]{"p_0/1", "p_1/1", "p_2_0/1"},
			new String[]{"p/1", "p_2/1", "p_2_2/1"}, "addslash=t", "slashspace=f");
		names(out, "colon-precedence", baseline,
			new String[]{"p", "p", "p_2"}, new String[]{"p_0 1:", "p_1 1:", "p_2_0 1:"},
			new String[]{"p 1:", "p_2 1:", "p_2_2 1:"}, "addslash=t", "addcolon=t");
		names(out, "nonunique-normalization", baseline,
			new String[]{"a b", "a_b"}, new String[]{"a_b", "a_b"}, new String[]{"a_b", "a_b"}, "uniquenames=f", "underscore=t");
	}

	/** Each direct read has explicit untouched sequence and quality oracles. */
	private static void names(final PrintWriter out, final String label, final boolean baseline,
		final String[] input, final String[] before, final String[] after, final String... options){
		final ReformatProcessor proc=processor(options);
		final String[] expected=baseline ? before : after;
		for(int i=0; i<input.length; i++){
			final Read read=read(input[i], i);
			final int selected=proc.processReadPair(read, null);
			observe(out, label+"-"+i, expected[i], read.id);
			check(selected==1, label+": valid single read must retain bit one.");
			check(expected[i].equals(read.id), label+": unexpected name "+read.id);
			check(Arrays.equals(read.bases, BASES) && Arrays.equals(read.quality, QUALITIES), label+": naming changed bases/qualities.");
		}
	}

	/** The same final key may occur once in each independent mate history. */
	private static void paired(final boolean baseline, final PrintWriter out){
		final String[] left={"a b", "a_b", "a_b_2"}, right={"a_b_2", "a_b", "a b"};
		final String[] beforeLeft={"a_b_0", "a_b_0", "a_b_2_0"}, beforeRight={"a_b_2_0", "a_b_0", "a_b_0"};
		final String[] afterLeft={"a_b", "a_b_2", "a_b_2_2"}, afterRight={"a_b_2", "a_b", "a_b_3"};
		final ReformatProcessor proc=processor("underscore=t");
		for(int i=0; i<left.length; i++){
			final Read r1=read(left[i], i), r2=read(right[i], i);
			pair(r1, r2);
			final int selected=proc.processReadPair(r1, r2);
			final String expectedLeft=(baseline ? beforeLeft : afterLeft)[i], expectedRight=(baseline ? beforeRight : afterRight)[i];
			observe(out, "mate-left-"+i, expectedLeft, r1.id);
			observe(out, "mate-right-"+i, expectedRight, r2.id);
			check(selected==3, "Both valid mates must survive.");
			check(expectedLeft.equals(r1.id) && expectedRight.equals(r2.id), "Independent mate histories differ.");
			check(r1.mate==r2 && r2.mate==r1 && r1.pairnum()==0 && r2.pairnum()==1, "Naming changed mate links/numbers.");
			check(Arrays.equals(r1.bases, BASES) && Arrays.equals(r2.bases, BASES) &&
				Arrays.equals(r1.quality, QUALITIES) && Arrays.equals(r2.quality, QUALITIES), "Naming changed mate payloads.");
		}
		for(final String option : new String[]{"slashspace=t", "slashspace=f", "addcolon=t"}){
			final ReformatProcessor tags=processor("addslash=t", option);
			for(int i=0; i<2; i++){
				final Read r1=read("tag", i), r2=read("tag", i);
				pair(r1, r2);
				final int selected=tags.processReadPair(r1, r2);
				final String stem=baseline ? (i==0 ? "tag_0" : "tag_1") : (i==0 ? "tag" : "tag_2");
				final String one=option.equals("addcolon=t") ? " 1:" : option.equals("slashspace=f") ? "/1" : " /1";
				final String two=option.equals("addcolon=t") ? " 2:" : option.equals("slashspace=f") ? "/2" : " /2";
				observe(out, "tag-left-"+option+"-"+i, stem+one, r1.id);
				observe(out, "tag-right-"+option+"-"+i, stem+two, r2.id);
				check(selected==3 && (stem+one).equals(r1.id) && (stem+two).equals(r2.id), "Mate tags must follow unique stems.");
			}
		}
	}

	/** Serial clones share histories; removed pairs do not reserve names, surviving asymmetric pairs do. */
	private static void cloneAndSelection(final boolean baseline, final PrintWriter out){
		final ReformatProcessor parent=processor(), clone=parent.clone();
		final Read first=read("shared", 0), second=read("shared", 1);
		final int a=clone.processReadPair(first, null), b=parent.processReadPair(second, null);
		final String expectedFirst=baseline ? "shared_0" : "shared", expectedSecond=baseline ? "shared_1" : "shared_2";
		observe(out, "clone-first", expectedFirst, first.id);
		observe(out, "clone-second", expectedSecond, second.id);
		check(a==1 && b==1 && expectedFirst.equals(first.id) && expectedSecond.equals(second.id), "Serial clones must share one naming history.");
		check(parent.recommendedWorkers()==1 && clone.recommendedWorkers()==1, "Shared histories require serial workers.");
		final ReformatProcessor missing=processor();
		final Read r1=read(null, 42), r2=read(null, 42);
		pair(r1, r2);
		check(missing.processReadPair(r1, r2)==3, "Null-ID pair must survive.");
		final String nullExpected=baseline ? "42_0" : "42";
		observe(out, "null-left", nullExpected, r1.id);
		observe(out, "null-right", nullExpected, r2.id);
		check(nullExpected.equals(r1.id) && nullExpected.equals(r2.id), "Null mate IDs must use the same original numeric fallback.");

		final ReformatProcessor filtered=processor();
		filtered.minReadLength=5;
		final Read dropped=read("filter", 0);
		check(filtered.processReadPair(dropped, null)==0 && dropped.id.equals("filter"), "Removed reads return before name allocation.");
		filtered.minReadLength=0;
		final Read kept=read("filter", 1);
		final String filteredExpected=baseline ? "filter_0" : "filter";
		check(filtered.processReadPair(kept, null)==1 && filteredExpected.equals(kept.id), "Removed reads must not reserve keys.");

		final ReformatProcessor asymmetric=processor();
		asymmetric.minReadLength=3;
		asymmetric.requireBothBad=true;
		final Read longRead=read("asym", 0), shortRead=new Read(new byte[]{'A', 'C'}, new byte[]{40, 40}, "asym", 0);
		pair(longRead, shortRead);
		check(asymmetric.processReadPair(longRead, shortRead)==1, "Only the long mate should be selected.");
		final Read laterLeft=read("asym", 1), laterRight=read("asym", 1);
		pair(laterLeft, laterRight);
		final String laterExpected=baseline ? "asym_1" : "asym_2";
		final int selected=asymmetric.processReadPair(laterLeft, laterRight);
		observe(out, "asymmetric-later-left", laterExpected, laterLeft.id);
		observe(out, "asymmetric-later-right", laterExpected, laterRight.id);
		check(selected==3 && laterExpected.equals(laterLeft.id) && laterExpected.equals(laterRight.id), "A surviving pair reserves both mate histories.");
	}

	/** New selector calls are deliberately absent from baseline execution. */
	private static void selectorControls(final PrintWriter out){
		final String[] input={"x a", "x\tb", "x a"};
		names(out, "explicit-identifier", false, input, null, new String[]{"x a", "x_2\tb", "x_3 a"}, "uniquenamescope=identifier");
		names(out, "complete-header", false, input, null, new String[]{"x a", "x\tb", "x a_2"}, "uniquenamescope=header");
		names(out, "header-mate-chain", false,
			new String[]{"p desc /1", "p desc /1", "p desc_2 /1"}, null,
			new String[]{"p desc /1", "p desc_2 /1", "p desc_2_2 /1"}, "uniquenamescope=header");
		names(out, "header-normalized", false,
			new String[]{"a b", "a_b", "a_b_2"}, null,
			new String[]{"a_b", "a_b_2", "a_b_2_2"}, "uniquenamescope=header", "underscore=t");
		names(out, "scope-without-enable", false, new String[]{"x a", "x a"}, null,
			new String[]{"x a", "x a"}, "uniquenames=f", "uniquenamescope=header");
		for(final String value : new String[]{null, "", "other"}){
			final ReformatProcessor proc=new ReformatProcessor();
			boolean rejected=false;
			try{proc.parse(value==null ? "uniquenamescope" : "uniquenamescope="+value, "uniquenamescope", value);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Missing/unknown uniqueness scope must fail clearly: "+value);
		}
	}

	/** Four real drivers read mixed space/tab headers; each child starts with fresh global parser state. */
	private static void cliCases(final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("input.fq");
		Files.write(input, fastq(new String[]{"x a", "x\tb", "x a"}).getBytes(StandardCharsets.US_ASCII));
		final Path launcher=Paths.get("dev/reformatnamecollisiontest.sh").toAbsolutePath();
		for(final String variant : new String[]{"streamer", "cris", "legacy", "legacyold"}){
			for(final boolean trim : new boolean[]{false, true}){
				for(final boolean trimFirst : new boolean[]{false, true}){
					runCli(root, input, launcher, variant, trim, trimFirst, false, baseline);
					if(!baseline){runCli(root, input, launcher, variant, trim, trimFirst, true, false);}
				}
			}
		}
	}

	/** Retains literal expected bytes before launch and actual child artifacts before assertions. */
	private static void runCli(final Path root, final Path input, final Path launcher, final String variant,
		final boolean trim, final boolean trimFirst, final boolean header, final boolean baseline) throws Exception{
		final String label=variant+"-trd"+trim+"-first"+trimFirst+(header ? "-header" : "-default");
		final Path output=root.resolve(label+".fq");
		final boolean legacy=variant.equals("legacy") || variant.equals("legacyold");
		//The saved tin=1 streamer baseline bypasses description trimming; fixed output must still honor trd.
		final boolean baselineTrim=trim && !variant.equals("streamer");
		//Legacy HashMap absence already preserves the first full header; modern primitive absence does not.
		final String[] expected=baseline && legacy ? (trim ? new String[]{"x", "x_2", "x_3"} : new String[]{"x a", "x\tb", "x a_2"}) :
			baseline ? (baselineTrim ? new String[]{"x_0", "x_1", "x_2"} : new String[]{"x a_0", "x\tb_0", "x a_1"}) :
			trim ? new String[]{"x", "x_2", "x_3"} : header ? new String[]{"x a", "x\tb", "x a_2"} : new String[]{"x a", "x_2\tb", "x_3 a"};
		final byte[] expectedBytes=fastq(expected).getBytes(StandardCharsets.US_ASCII);
		Files.write(root.resolve(label+".expected.fq"), expectedBytes);
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("bash", launcher.toString(), "driver="+variant,
			"in="+input, "out="+output, "t=1", "wt=1", "tin=1", "tout=1", "int=f",
			"qin=33", "qout=33", "ow=f"));
		//Legacy drivers process lists serially and do not recognize the modern ordered option.
		if(!legacy){args.add("ordered=t");}
		if(trimFirst){args.add("trd="+trim);}
		args.add("uniquenames=t");
		if(header){args.add("uniquenamescope=header");}
		if(!trimFirst){args.add("trd="+trim);}
		Files.write(root.resolve(label+".args"), args, StandardCharsets.UTF_8);
		final Process child=new ProcessBuilder(args).redirectOutput(root.resolve(label+".stdout").toFile())
			.redirectError(root.resolve(label+".stderr").toFile()).start();
		final int code=child.waitFor();
		children++;
		Files.write(root.resolve(label+".exit"), (code+"\n").getBytes(StandardCharsets.US_ASCII));
		check(code==0, label+": driver must complete successfully; see retained stderr.");
		check(Arrays.equals(expectedBytes, Files.readAllBytes(output)), label+": exact FASTQ output differs.");
	}

	/** Delegates to the requested actual driver in the child's fresh JVM; preserves its original failure. */
	private static void driver(final String[] args) throws Exception{
		final String variant=args[0].substring("driver=".length());
		final String className;
		if(variant.equals("streamer")){className="jgi.ReformatStreamer";}
		else if(variant.equals("cris")){className="jgi.Reformat2";}
		else if(variant.equals("legacy")){className="jgi.ReformatReads";}
		else if(variant.equals("legacyold")){className="jgi.ReformatReadsOld";}
		else{throw new IllegalArgumentException(args[0]);}
		final Class<?> type=Class.forName(className);
		try{type.getMethod("main", String[].class).invoke(null, (Object)Arrays.copyOfRange(args, 1, args.length));}
		catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof Exception){throw (Exception)cause;}
			throw new AssertionError("Driver failed", cause);
		}
	}

	/** Actual parser/postParse setup, with explicit reset of the process-wide slash style for each case. */
	private static ReformatProcessor processor(final String... options){
		final ReformatProcessor proc=new ReformatProcessor();
		proc.parse("slashspace=t", "slashspace", "t");
		proc.parse("uniquenames=t", "uniquenames", "t");
		for(final String option : options){
			final String[] pair=Parse.splitOnFirst(option, '=');
			if(!proc.parse(option, pair[0], pair[1])){throw new AssertionError("Unconsumed fixture option: "+option);}
		}
		proc.setFromParser(new Parser());
		proc.minReadLength=0;
		proc.maxReadLength=0;
		proc.postParse();
		return proc;
	}

	/** Fresh mutable payloads keep naming and payload assertions independent. */
	private static Read read(final String name, final long id){
		return new Read(BASES.clone(), QUALITIES.clone(), name, id);
	}

	/** Establishes reciprocal mates with explicit pair numbers. */
	private static void pair(final Read first, final Read second){
		first.mate=second; second.mate=first;
		first.setPairnum(0); second.setPairnum(1);
	}

	/** Literal four-base FASTQ records; headers alone vary across these cases. */
	private static String fastq(final String[] names){
		final StringBuilder text=new StringBuilder();
		for(final String name : names){text.append('@').append(name).append("\nACGT\n+\nIIII\n");}
		return text.toString();
	}

	/** Flushes evidence before a failed oracle can stop the run. */
	private static void observe(final PrintWriter out, final String label, final String expected, final String actual){
		out.println(label+'\t'+escape(expected)+'\t'+escape(actual));
		out.flush();
		if(out.checkError()){throw new AssertionError("Could not save observation: "+label);}
	}

	/** Escapes delimiters without changing the actual tested ID. */
	private static String escape(final String text){
		return text==null ? "<null>" : text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
	}

	/** Literal oracles remain active even with assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static final byte[] BASES={'A', 'C', 'G', 'T'}, QUALITIES={40, 40, 40, 40};
	/** Literal sixteen-record block; duplicates and natural/mate suffixes retain their original spelling. */
	private static final String[] BENCH_HEADERS={"x", "x", "x_2", "x", "x_02", "x_2", "x_2_2", "x_2",
		"n a", "n_2 b", "n c", "p/1", "p/1", "p_2/1", "q", "q_2"};
	/** Eleven independent full-header histories; maximum counter3*1048576 fits an int. */
	private static final int[] BENCH_KEYS={0, 0, 1, 0, 2, 1, 3, 1, 4, 5, 6, 7, 7, 8, 9, 10};
	//53 direct records*3 +18 paired +8 lifecycle +32 baseline CLI +1 sentry=218.
	//Fixed adds45 selector checks and32 additional CLI checks:295 total.
	private static final int BASELINE_CHECKS=218, FIXED_CHECKS=295;
	private static int checks, children;
}
