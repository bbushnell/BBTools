package test.assemble;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.stream.Stream;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

/**
 * Finite Linux launcher tests for assembly-wrapper preservation and threshold errors.
 * Every child uses tiny generated inputs in a fresh directory, a bounded heap,
 * and GNU timeout to terminate its process group if a failure path hangs.
 * Expected results concern publication and input preservation, not assembly quality.
 * @author Brian Bushnell, Shinobu
 */
public final class AssembleWrapperContractTest{

	/** Runs wrappers by default; case=thresholds or multi-thresholds selects numeric-option cases.
	 * Retains commands, child logs and observations for the selected group or named case. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, AssembleWrapperContractTest.class, false);
		try{
			Path tree=null, output=null;
			String selected="all";
			boolean assertions=true;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String name=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(name.equals("bbtools")){tree=Paths.get(value).toRealPath();}
				else if(name.equals("outdir")){output=Paths.get(value).toAbsolutePath();}
				else if(name.equals("case")){selected=value;}
				else if(name.equals("assertions")){assertions=Parse.parseBoolean(value);}
				else{throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(tree==null || output==null){throw new IllegalArgumentException("Specify bbtools=DIR and a new outdir=DIR.");}
			final AssembleWrapperContractTest test=new AssembleWrapperContractTest(tree, output, assertions);
			final boolean multi=selected.equals("multi-thresholds") || selected.startsWith("multi-threshold-");
			final boolean thresholds=selected.equals("thresholds") || selected.startsWith("threshold-");
			final String[] cases=multi ? MULTI_THRESHOLD_CASES : thresholds ? THRESHOLD_CASES : CASES;
			boolean found=false;
			for(String name : cases){
				if(selected.equals("all") || selected.equals("thresholds") || selected.equals("multi-thresholds") || selected.equals(name)){
					test.runCase(name); found=true;
				}
			}
			if(!found){throw new IllegalArgumentException("Unknown case: "+selected);}
			System.out.println("ASSEMBLE_WRAPPER_CASES checks="+test.checks+" failed="+test.failed);
			if(test.failed>0){throw new AssertionError("Wrapper contract failures; inspect "+output);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Refuses existing result roots, so even a mistaken invocation cannot reuse user data. */
	private AssembleWrapperContractTest(final Path tree_, final Path root_, final boolean assertions_) throws IOException{
		tree=tree_;root=root_;assertions=assertions_;
		if(!Files.isExecutable(Paths.get("/usr/bin/timeout"))){throw new IOException("This launcher fixture requires GNU timeout.");}
		Files.createDirectory(root);
		results=root.resolve("checks.tsv");
		Files.write(results, "case\tchild_exit\tpassed\tdetail\n".getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Each case owns all files that the deliberately unfixed wrapper might overwrite or delete. */
	private void runCase(final String name) throws Exception{
		assert(Arrays.asList(CASES).contains(name) || Arrays.asList(THRESHOLD_CASES).contains(name)
				|| Arrays.asList(MULTI_THRESHOLD_CASES).contains(name)) :
				"Only the named finite cases may create test paths.";
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path input=dir.resolve(name.endsWith("-equals") ? "tid_123=input.fa" : "tid_123_input.fa");
		final byte[] original=fixture();
		Files.write(input, original, StandardOpenOption.CREATE_NEW);
		if(name.startsWith("threshold-") || name.startsWith("multi-threshold-")){
			runThreshold(name, dir, input, original); return;
		}
		final Path out=dir.resolve("output.fa"), temps=dir.resolve("temps");
		int code;
		if(name.equals("reassemble-alias")){
			code=launch(dir, "reassemble.sh", "in="+input, "out="+input);
			record(name, code, code!=0 && unchanged(input, original), "input must survive rejected output alias");
		}else if(name.equals("reassemble-overwrite")){
			final byte[] saved="existing output must survive\n".getBytes(StandardCharsets.US_ASCII);
			Files.write(out, saved, StandardOpenOption.CREATE_NEW);
			code=launch(dir, "reassemble.sh", "in="+input, "out="+out, "ow=f");
			record(name, code, code!=0 && unchanged(out, saved) && unchanged(input, original), "overwrite=f must preserve both existing files");
		}else if(name.equals("reassemble-concat-failure")){
			final Path missing=blockedOutput(dir);
			code=launch(dir, "reassemble.sh", "in="+input, "out="+missing, "tempdir="+temps);
			final boolean completed=log(dir).contains("Assembly succeeded:");
			record(name, code, code!=0 && (!completed || hasTemporaryAssembly(temps)) && unchanged(input, original),
					"failed final write must retain any completed assembly; assembly_completed="+completed);
		}else if(name.startsWith("tadwrapper-")){
			final Path pattern=dir.resolve("assembly-%.fa"), best=dir.resolve("assembly-31.fa");
			final boolean failure=name.equals("tadwrapper-copy-failure");
			final Path destination=failure ? blockedOutput(dir) : out;
			code=launch(dir, "tadwrapper.sh", "in="+input, "out="+pattern, "outfinal="+destination,
					"assemblydefaults=f", "branchrescue=f");
			final boolean selected=log(dir).contains("Recommended K:");
			final boolean passed=failure ? code!=0 && (!selected || hasSequence(best)) : code==0 && selected && hasSequence(destination);
			record(name, code, passed && unchanged(input, original), failure ?
					"failed copy must retain any selected assembly; selection_completed="+selected : "normal final publication control");
		}else{
			final ArrayList<String> options=new ArrayList<String>(Arrays.asList("in="+input, "out="+out));
			if(name.equals("reassemble-temp-control")){options.add("tempdir="+temps);}
			if(name.startsWith("reassemble-long-")){options.add("k=62"); options.add("packed=t");}
			code=launch(dir, "reassemble.sh", options.toArray(new String[0]));
			record(name, code, code==0 && hasSequence(out) && unchanged(input, original), "valid input must produce an assembly without changing input");
		}
	}

	/** Checks real CLI rejection before output, preserving valid finite and unbounded ceilings. */
	private void runThreshold(final String name, final Path dir, final Path input, final byte[] original) throws Exception{
		final boolean multi=name.startsWith("multi-threshold-");
		final String prefix=multi ? "multi-threshold-" : "threshold-";
		final String[] flags=multi ? MULTI_THRESHOLD_FLAGS : THRESHOLD_FLAGS;
		final String[] names=multi ? MULTI_THRESHOLD_NAMES : THRESHOLD_NAMES;
		assert(flags.length==names.length) : "Each tested CLI alias needs its matching diagnostic name.";
		final Path out=dir.resolve("output.fa");
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList("in="+input, "out="+out,
				"classifygraphcontigs=t", "lowdepthcontigdiag=t", "pop=f", "sweeplen=0"));
		if(multi){
			//Exercise Config and its direct final-graph forwarding without automatic fusion models.
			options.addAll(Arrays.asList("assemblydefaults=f", "branchrescue=f", "k=31,21",
					"fusek=none", "bridgek=none", "graphk=21", "fusenet=null", "showstats=f"));
		}
		int invalid=-1;
		for(int i=0; i<flags.length; i++){
			if(name.equals(prefix+flags[i]+"-nan")){invalid=i;}
		}
		if(invalid>=0){
			options.add(flags[invalid]+"=NaN");
		}else if(name.equals(prefix+"zero")){
			for(String flag : flags){options.add(flag+"=0");}
		}else if(name.equals(prefix+"unbounded")){
			options.addAll(Arrays.asList("ldcmaxcov=Infinity", "gclmc=Infinity", "gcmf=Infinity", "gchf=Infinity"));
			if(multi){options.add("ckmdr=Infinity");}
		}else if(!name.equals(prefix+"control")){
			throw new IllegalArgumentException("Unknown threshold fixture: "+name);
		}
		final int code=launch(dir, "tadpole.sh", options.toArray(new String[0]));
		if(invalid<0){
			record(name, code, code==0 && hasSequence(out) && unchanged(input, original),
					"valid threshold control must produce an assembly and preserve input");
			return;
		}
		final boolean rejected=code!=0 && log(dir).contains(names[invalid]+" must") && !Files.exists(out);
		if((code!=0 && !rejected) || (code==0 && !hasSequence(out))){
			throw new IOException("Unexpected child failure or empty assembly in "+name+"; inspect "+dir);
		}
		record(name, code, rejected && unchanged(input, original),
				"NaN must be rejected before output with a diagnostic naming "+names[invalid]);
	}

	/** A regular file cannot be the parent of an output, even when writers create missing directories. */
	private static Path blockedOutput(final Path directory) throws IOException{
		final Path blocker=directory.resolve("blocked-parent");
		Files.write(blocker, new byte[]{'x'}, StandardOpenOption.CREATE_NEW);
		return blocker.resolve("final.fa");
	}

	/** Saved child diagnostics also distinguish early rejection from loss after successful assembly. */
	private static String log(final Path directory) throws IOException{
		return new String(Files.readAllBytes(directory.resolve("child.log")), StandardCharsets.UTF_8);
	}

	/** Executes the real shell launcher; timeout exits cannot masquerade as a correct rejection. */
	private int launch(final Path dir, final String script, final String... options) throws Exception{
		assert(script.equals("reassemble.sh") || script.equals("tadwrapper.sh") || script.equals("tadpole.sh")) :
				"Only the reviewed assembly launchers belong to this fixture.";
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--kill-after=5s", "60s",
				"bash", tree.resolve(script).toString(), "-Xmx256m", "-Xms32m", "-eoom", assertions ? "-ea" : "-da",
				"t=1", "buildthreads=1", "k=31", "prealloc=f", "prefilter=0", "minprob=0", "mcs=1", "mce=1",
				"mincontig=1", "minextension=1", "trimends=0", "shave=f", "rinse=f", "ow=t", "interleaved=f"));
		command.addAll(Arrays.asList(options));
		Files.write(dir.resolve("command.args"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final Process child=new ProcessBuilder(command).directory(tree.toFile()).redirectErrorStream(true)
				.redirectOutput(dir.resolve("child.log").toFile()).start();
		child.getOutputStream().close();
		final int code=child.waitFor();
		Files.write(dir.resolve("child.exit"), (code+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(code==124 || code==137){throw new IOException("Wrapper exceeded the finite child deadline; inspect "+dir);}
		return code;
	}

	/** Builds a deterministic 512-base sequence repeated four times to supply depth without bulk input. */
	private static byte[] fixture(){
		final Random random=new Random(714031);
		final char[] alphabet={'A','C','G','T'}, sequence=new char[512];
		for(int i=0; i<sequence.length; i++){sequence[i]=alphabet[random.nextInt(4)];}
		final StringBuilder fasta=new StringBuilder();
		for(int i=0; i<4; i++){fasta.append(">sample").append(i).append('\n').append(sequence).append('\n');}
		return fasta.toString().getBytes(StandardCharsets.US_ASCII);
	}

	/** Exact preservation oracle; absent files are a failure, never an empty-file success. */
	private static boolean unchanged(final Path file, final byte[] expected) throws IOException{
		return Files.isRegularFile(file) && Arrays.equals(expected, Files.readAllBytes(file));
	}

	/** Requires at least one sequence line after a FASTA header, not merely a nonempty file. */
	private static boolean hasSequence(final Path file) throws IOException{
		if(!Files.isRegularFile(file)){return false;}
		boolean header=false;
		for(String line : Files.readAllLines(file, StandardCharsets.US_ASCII)){
			if(line.startsWith(">")){header=true;}
			else if(header && !line.trim().isEmpty()){return true;}
		}
		return false;
	}

	/** Inspects only this case's private temporary directory after the child terminates. */
	private static boolean hasTemporaryAssembly(final Path directory) throws IOException{
		if(!Files.isDirectory(directory)){return false;}
		try(Stream<Path> paths=Files.list(directory)){
			for(Path file : (Iterable<Path>)paths::iterator){
				if(file.getFileName().toString().startsWith("reassemble_temp_") && hasSequence(file)){return true;}
			}
		}
		return false;
	}

	/** Persists every independent oracle even when earlier expected-contract checks fail. */
	private void record(final String name, final int code, final boolean passed, final String detail) throws IOException{
		checks++;if(!passed){failed++;}
		Files.write(results, (name+"\t"+code+"\t"+passed+"\t"+detail+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	private final Path tree, root, results;
	private final boolean assertions;
	private int checks, failed;
	private static final String[] CASES={"reassemble-control", "reassemble-temp-control", "reassemble-alias",
			"reassemble-overwrite", "reassemble-equals", "reassemble-concat-failure", "tadwrapper-control", "tadwrapper-copy-failure",
			"reassemble-long-control", "reassemble-long-equals"};
	private static final String[] THRESHOLD_FLAGS={"ldcmaxcov", "ldcfrac", "gclmc", "gclf", "gcmf", "gchf"};
	private static final String[] THRESHOLD_NAMES={"lowDepthContigMaxCov", "lowDepthContigFraction", "graphClassLowMaxCov",
			"graphClassLowFraction", "graphClassMediumFraction", "graphClassHighFraction"};
	private static final String[] THRESHOLD_CASES={"threshold-control", "threshold-zero", "threshold-unbounded",
			"threshold-ldcmaxcov-nan", "threshold-ldcfrac-nan", "threshold-gclmc-nan",
			"threshold-gclf-nan", "threshold-gcmf-nan", "threshold-gchf-nan"};
	/** Separate names preserve the old single-K group and its saved baseline artifacts. */
	private static final String[] MULTI_THRESHOLD_FLAGS={"ldcmaxcov", "ldcfrac", "gclmc", "gclf", "gcmf", "gchf", "ckmdr"};
	private static final String[] MULTI_THRESHOLD_NAMES={"lowDepthContigMaxCov", "lowDepthContigFraction", "graphClassLowMaxCov",
			"graphClassLowFraction", "graphClassMediumFraction", "graphClassHighFraction", "crosskmaxdepthratio"};
	private static final String[] MULTI_THRESHOLD_CASES={"multi-threshold-control", "multi-threshold-zero", "multi-threshold-unbounded",
			"multi-threshold-ldcmaxcov-nan", "multi-threshold-ldcfrac-nan", "multi-threshold-gclmc-nan",
			"multi-threshold-gclf-nan", "multi-threshold-gcmf-nan", "multi-threshold-gchf-nan", "multi-threshold-ckmdr-nan"};
}
