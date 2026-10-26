package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Exercises final-polish routing through real isolated fake child processes.
 * @author Brian Bushnell, Fischl */
public final class TadPipePolishTest {

	/** Accepts the installation path; test data and logs remain in a unique directory. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify the BBTools installation directory.");
		final TadPipePolishTest test=new TadPipePolishTest(Paths.get(args[0]).toAbsolutePath());
		test.run();
		test.runScaffolds();
		System.out.println("TADPIPE_POLISH_TEST_PASS checks="+checks+" dir="+test.root);
	}

	/** Tiny launchers expose dataflow without needing an assembler or mapper. */
	private TadPipePolishTest(final Path code_) throws Exception{
		code=code_;
		root=Files.createTempDirectory("tadpipe-polish-test-");
		tools=Files.createDirectory(root.resolve("tools with spaces"));
		Files.createDirectory(tools.resolve("networks"));
		Files.write(tools.resolve("networks/tadpole_fusion.bbnet"), Arrays.asList("fixture"), StandardCharsets.UTF_8);
		Files.write(tools.resolve("networks/tadpole_contig_rescue.bbnet"), Arrays.asList("fixture"), StandardCharsets.UTF_8);
		reads=Files.write(root.resolve("input.fq"), Arrays.asList("untouched reads"), StandardCharsets.UTF_8);
		final String common="#!/bin/bash\nset -eu\nprintf 'ARG:<%s>\\n' \"$@\"\n"+
				"for arg in \"$@\"; do case $arg in out=*|outu=*) printf '>unpolished\\nACGT\\n' > \"${arg#*=}\";; esac; done\n";
		for(String tool : Arrays.asList("clumpify.sh", "bbmerge.sh", "tadpole.sh", "reformat.sh", "bbduk.sh", "cat.sh")){
			Files.write(tools.resolve(tool), common.getBytes(StandardCharsets.UTF_8));
		}
		final String polish="#!/bin/bash\nset -eu\nprintf 'ARG:<%s>\\n' \"$@\"\n"+
				"in= ref= out= report=\nfor arg in \"$@\"; do case $arg in in=*) in=${arg#*=};; assembly=*) ref=${arg#*=};; out=*) out=${arg#*=};; report=*) report=${arg#*=};; zl=*) exit 31;; esac; done\n"+
				"[[ $in == */paired.fq.gz ]]\ntest -s \"$in\"\ntest -s \"$ref\"\n"+
				"test ! -e \"$(dirname \"$0\")/fail-polish\" || exit 7\n"+
				"printf '>polished\\nACGA\\n' > \"$out\"\nprintf 'edits\\n' > \"$report\"\n";
		Files.write(tools.resolve("tippolish.sh"), polish.getBytes(StandardCharsets.UTF_8));
		Files.createDirectory(tools.resolve("dev"));
		final String scaffold="#!/bin/bash\nset -eu\nprintf 'ARG:<%s>\\n' \"$@\"\n"+
				"contigs= out= joins= candidates=\nfor arg in \"$@\"; do case $arg in contigs=*) contigs=${arg#*=};; outassembly=*) out=${arg#*=};; outjoins=*) joins=${arg#*=};; out=*) candidates=${arg#*=};; in=*|extra=*) test -s \"${arg#*=}\";; esac; done\n"+
				"test -s \"$contigs\"\ntest ! -e \"$(dirname \"$0\")/fail-scaffold\" || exit 9\n"+
				"printf '>scaffold\\nACGTNNNNNNNNNNACGA\\n' > \"$out\"\nprintf 'links\\n' > \"$joins\"\nprintf 'candidates\\n' > \"$candidates\"\n";
		Files.write(tools.resolve("dev/hollowbridge.sh"), scaffold.getBytes(StandardCharsets.UTF_8));
	}

	/** Scaffolding defaults on, runs last, retains contigs and never publishes a failed-stage fallback. */
	private void runScaffolds() throws Exception{
		Files.delete(tools.resolve("fail-polish"));
		final Path on=launchInternal("scaffold_default", true, false, "delete=f");
		check(text(root.resolve("scaffold_default.fa")).contains("NNNNNNNNNN"), "The default pipeline must publish scaffold output.");
		check(text(on.resolve("assembly.fa")).contains("unpolished"), "Scaffolding must retain original contigs.");
		final String log=text(on.resolve("12_scaffold.log"));
		for(String option : Arrays.asList("k=224", "fillgaps=f", "slide=0", "flankcheck=t", "join=t", "interleaved=f", "t=1", "-Xmx256m")){
			check(log.contains("ARG:<"+option+">"), "Scaffolding lost a default or resource limit: "+option);
		}
		check(log.contains("ARG:<extra="), "Merged and residual unmerged reads must both support scaffolding.");
		check(text(on.resolve("10_assemble.log")).contains("ARG:<scaffold=f>"), "Assembly must not scaffold before the final pipeline stage.");
		final Path polished=launchInternal("scaffold_polished", true, false, "polish=t", "scaffold_threads=3", "scaffold_t=1");
		check(text(polished.resolve("12_scaffold.log")).contains("ARG:<contigs="+polished.resolve("polished.fa")+">"),
				"Scaffolding must use the polished assembly, not stale prepolish contigs.");
		check(!text(polished.resolve("12_scaffold.log")).contains("ARG:<t=3>"), "Thread aliases must use the last value.");
		final Path unmerged=launchInternal("scaffold_unmerged", true, false, "merge=f", "extend=f", "scaffold_fillgaps=t");
		final String unmergedLog=text(unmerged.resolve("12_scaffold.log"));
		check(!unmergedLog.contains("ARG:<extra=") && unmergedLog.contains("ARG:<fillgaps=t>"),
				"An unmerged-only pipeline must not duplicate evidence and must honor explicit filling.");
		final Path off=launchInternal("scaffold_disabled", true, false, "scaffold=f");
		check(!Files.exists(off.resolve("12_scaffold.log")), "scaffold=f must restore the unscaffolded workflow.");
		for(String option : Arrays.asList("scaffold_in=x", "scaffold_outassembly=x", "scaffold_ow=t", "assemble_scaffold=t")){
			launchInternal("scaffold_invalid"+(++invalid), false, false, option);
		}
		Files.write(tools.resolve("dev/fail-scaffold"), Arrays.asList("fail"), StandardCharsets.UTF_8);
		final Path failed=launchInternal("scaffold_failed", false, false);
		check(!Files.exists(root.resolve("scaffold_failed.fa")), "Failed scaffolding must not publish an unscaffolded fallback.");
		check(Files.exists(failed.resolve("assembly.fa")) && Files.exists(failed.resolve("paired.fq.gz")),
				"Failed scaffolding must preserve contigs and read evidence.");
	}

	/** Default off, success, failure and bounded override cases have independent outputs. */
	private void run() throws Exception{
		final Path off=launch("off", true);
		check(text(root.resolve("off.fa")).contains("unpolished"), "Default-off publication must retain assembly output.");
		check(!Files.exists(off.resolve("11_polish.log")), "Default off must not launch the polisher.");
		final Path on=launch("on", true, "polish=t", "polish_minoverlap=40", "polish_maxoverlap=180",
				"polish_indexed=f", "polish_threads=3", "polish_t=1");
		check(text(root.resolve("on.fa")).contains("polished\nACGA"), "Only the final polished result may be published.");
		check(text(on.resolve("assembly.fa")).contains("unpolished"), "The original assembly must remain available.");
		check(Files.exists(on.resolve("polish.edits.tsv")), "The edit report must survive pipeline completion.");
		check(!Files.exists(on.resolve("paired.fq.gz")), "Managed reads may be deleted only after polishing succeeds.");
		final String log=text(on.resolve("11_polish.log"));
		for(String flag : Arrays.asList("interleaved=t", "minoverlap=40", "maxoverlap=180", "indexed=f", "t=1", "-Xmx256m")){
			check(log.contains("ARG:<"+flag+">"), "Polishing lost its routed flag: "+flag);
		}
		check(!log.contains("ARG:<t=3>"), "polish_t and polish_threads must share last-value-wins semantics.");
		for(String phase : Arrays.asList("01_dedupe", "02_ecco", "03_clump", "04_ecc", "10_assemble")){
			final String before=text(off.resolve(phase+".log")).replace(off.toString(), "WORK");
			final String after=text(on.resolve(phase+".log")).replace(on.toString(), "WORK");
			check(before.equals(after), "Enabling polishing must not change previous stage commands: "+phase);
		}
		final Path noDedupe=launch("no_dedupe", true, "polish=t", "dedupe=f", "delete=f");
		check(Files.exists(noDedupe.resolve("paired.fq.gz")), "delete=f retains paired evidence even when deduplication is disabled.");
		final Path dry=launch("dry", true, "polish=t", "dryrun=t");
		check(!Files.exists(root.resolve("dry.fa")) && !Files.exists(dry.resolve("11_polish.log")), "Dry run must execute no polishing.");
		check(text(dry.resolve("commands.sh")).contains("tippolish.sh"), "Dry run must describe the optional final stage.");
		for(String option : Arrays.asList("polish_in=x", "polish_assembly=x", "polish_ref=x", "polish_report=x",
				"polish_temp=x", "polish_bbtools=x", "polish_ow=t", "dedupe_ecc=t", "dedupe_ecco=t", "dedupe_mergefirst=t")){
			launch("invalid"+(++invalid), false, "polish=t", option);
		}
		Files.write(tools.resolve("fail-polish"), Arrays.asList("fail"), StandardCharsets.UTF_8);
		final Path failed=launch("failed", false, "polish=t");
		check(!Files.exists(root.resolve("failed.fa")), "A failed polish must not publish an unpolished fallback.");
		check(Files.exists(failed.resolve("paired.fq.gz")) && Files.exists(failed.resolve("assembly.fa")), "Failure must retain evidence and unpolished sequence.");
		check(text(reads).equals("untouched reads\n"), "Pipeline must not mutate user reads.");
	}

	/** Runs the real launcher and finds only its explicitly printed unique work path. */
	private Path launch(final String name, final boolean success, final String... extras) throws Exception{
		return launchInternal(name, success, true, extras);
	}

	/** Uses native defaults when requested, while legacy polishing fixtures explicitly disable scaffolding. */
	private Path launchInternal(final String name, final boolean success, final boolean disableScaffold, final String... extras) throws Exception{
		assert(name!=null) : "Every test execution needs a distinct output and diagnostic log.";
		final Path log=root.resolve(name+".log");
		final List<String> command=new ArrayList<String>(Arrays.asList("bash", code.resolve("tadpipe.sh").toString(),
				"in="+reads, "out="+root.resolve(name+".fa"), "bbtools="+tools, "temp="+root, "t=1", "-Xmx256m"));
		if(disableScaffold){command.add("scaffold=f");}
		command.addAll(Arrays.asList(extras));
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		child.getOutputStream().close();
		check((child.waitFor()==0)==success, "Unexpected result for "+name+"; see "+log);
		for(String line : Files.readAllLines(log, StandardCharsets.UTF_8)){
			if(line.startsWith("TadPipe work directory: ")){return Paths.get(line.substring("TadPipe work directory: ".length()));}
		}
		check(!success, "Successful pipeline must print its work directory.");
		return null;
	}

	/** Reads small test diagnostics, never sequencing data. */
	private static String text(final Path file) throws Exception{return new String(Files.readAllBytes(file), StandardCharsets.UTF_8);}

	/** Test expectations remain fatal without enabled assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private final Path code, root, tools, reads;
	private int invalid;
	private static int checks;
}
