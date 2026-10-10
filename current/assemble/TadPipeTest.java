package assemble;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteStreamWriter;
import stream.ConcurrentGenericReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ListNum;
/**
 * Deterministic orchestration tests with tiny isolated child processes.
 * These test dataflow and failure contracts, not biological assembly quality.
 * @author Brian Bushnell, Fischl
 */
public class TadPipeTest {

	public static void main(String[] args) throws Exception{
		final TadPipeTest test=new TadPipeTest();
		test.testReadQuota();
		test.testPlan();
		test.testAssemblyOverrides();
		test.testSuccess();
		test.testFailures();
		test.testValidation();
		test.testEmptyCompressed();
		System.out.println("TADPIPE_TEST_PASS checks="+test.checks+" scratch="+test.root);
	}

	/** Reproduces the extra sampling prefetch deterministically, without timing races. */
	private void testReadQuota() throws Exception{
		for(int quota : new int[] {1,3}){
			final QuotaInput input=new QuotaInput(quota);
			final ConcurrentGenericReadInputStream cris=new ConcurrentGenericReadInputStream(input, null, quota);
			cris.start();
			long count=0;
			for(ListNum<Read> list=cris.nextList(); list!=null; list=cris.nextList()){
				final boolean end=list.list.isEmpty();
				count+=list.list.size();
				cris.returnList(list);
				if(end){break;}
			}
			//Join the captured producer BEFORE close can set shutdown and hide its
			//post-quota hasMore call. The raw queue has room for its terminal list.
			input.reader.join();
			cris.close();
			check(count==quota && !cris.errorState(), "Limited reads must preserve the requested root count.");
			check(input.postQuotaProbes==0, "A completed sampling quota must not prefetch into a concurrently closed producer.");
		}
	}

	/** Tiny source whose hasMore call is observable even after the requested quota. */
	private static final class QuotaInput extends ReadInputStream {
		QuotaInput(final int quota_){quota=quota_;}
		@Override
		public ArrayList<Read> nextList(){
			reader=Thread.currentThread();
			final ArrayList<Read> reads=new ArrayList<Read>();
			for(int i=0; i<2 && cursor<10; i++, cursor++){
				reads.add(new Read(new byte[] {'A'}, null, "probe"+cursor, cursor));
			}
			return reads.isEmpty() ? null : reads;
		}
		@Override
		public boolean hasMore(){
			if(cursor>=quota){postQuotaProbes++;}
			return cursor<10;
		}
		@Override
		public boolean close(){return false;}
		@Override
		public boolean paired(){return false;}
		@Override
		public String fname(){return "quota-probe";}
		@Override
		public void restart(){throw new UnsupportedOperationException("Test source is single-use.");}
		final int quota;
		int cursor, postQuotaProbes;
		volatile Thread reader;
	}

	/** Builds a fake installation whose launchers each execute in their own process. */
	TadPipeTest() throws IOException{
		root=Files.createTempDirectory("tadpipe-test ");
		tools=Files.createDirectories(root.resolve("tools with spaces"));
		Files.createDirectories(tools.resolve("networks"));
		Files.write(tools.resolve("networks/tadpole_fusion.bbnet"), Arrays.asList("fake model"), StandardCharsets.UTF_8);
		in=Files.write(root.resolve("reads R1.fq"), Arrays.asList("input one"), StandardCharsets.UTF_8);
		in2=Files.write(root.resolve("reads R2.fq"), Arrays.asList("input two"), StandardCharsets.UTF_8);
		final String script="#!/bin/bash\nset -eu\nprintf '%s\\n' \"$$\" >> pids.txt\n"+
				"printf 'ARG:<%s>\\n' \"$@\"\n"+
				"for arg in \"$@\"; do [[ $arg != testfail=t ]] || exit 7; done\n"+
				"for arg in \"$@\"; do [[ $arg != omit=t ]] || exit 0; done\n"+
				"for arg in \"$@\"; do case $arg in out=*|outu=*) printf '>fixture\\nACGT\\n' > \"${arg#*=}\";; esac; done\n";
		for(String tool : Arrays.asList("clumpify.sh", "bbmerge.sh", "tadpole.sh", "reformat.sh", "bbduk.sh", "cat.sh")){
			Files.write(tools.resolve(tool), script.getBytes(StandardCharsets.UTF_8));
		}
		check(Files.isRegularFile(in2), "Fixture must contain a distinct mate file.");
	}

	/** Checks recipe order, mate routing, resources, aliases, and specialized overrides. */
	private void testPlan() throws Exception{
		final TadPipe.Config config=config("planned.fa", "dedupe_subs=0", "extend_el=20", "extendu_el=7");
		config.validate();
		final TadPipe pipe=new TadPipe(config, System.err);
		pipe.work=Files.createDirectory(root.resolve("plan"));
		final ArrayList<TadPipe.Stage> plan=pipe.plan();
		check(plan.size()==13, "Default pipeline must retain all 13 processing stages.");
		check(plan.get(0).args.get("in2").equals(in2.toString()), "Read two must reach deduplication.");
		check(plan.get(0).args.get("optical").equals("f"), "SRA deduplication cannot require optical coordinates.");
		check(plan.get(0).args.get("subs").equals("0"), "Duplicate criteria must be overridable.");
		check(plan.get(2).args.get("unpair").equals("t") && plan.get(2).args.get("repair").equals("t"),
				"Clumpify must correct both mates and restore their pairing.");
		check(plan.get(6).args.get("extra").split(",").length==2, "Second REM pass must retain accumulated merge evidence.");
		check(plan.get(10).args.get("el").equals("20"), "Generic extension control applies to merged reads.");
		check(plan.get(11).args.get("el").equals("7"), "Specific extension control wins for unmerged reads.");
		check(plan.get(10).args.get("extra").equals(plan.get(11).args.get("in")), "Merged extension uses unextended residuals.");
		check(plan.get(11).args.get("extra").equals(plan.get(10).args.get("in")), "Unmerged extension uses unextended merged reads.");
		check(plan.get(12).args.get("k").equals("124,300,96,64,32"), "Assembly uses the ordered A124 recipe, not TadpoleWrapper.");
		for(TadPipe.Stage stage : plan){
			check(stage.command.contains("-Xmx256m") && stage.command.contains("t=1"), "Resources must reach each child.");
			check(stage.command.get(1).startsWith(tools.toString()), "All phases must use the same BBTools installation.");
		}
		final TadPipe.Config disabled=config("disabled.fa", "dedupe=f", "ecco=f", "clump=f", "ecc=f", "merge=f", "extend=f", "qtrim=f", "nn=f");
		disabled.validate();
		final TadPipe simple=new TadPipe(disabled, System.err);
		simple.work=Files.createDirectory(root.resolve("simple"));
		final ArrayList<TadPipe.Stage> small=simple.plan();
		check(small.size()==2 && small.get(0).tool.equals("reformat.sh"), "Skipping dedupe must still normalize mate files.");
		check(!small.get(1).args.containsKey("fusenet"), "nn=f removes the neural gate, not merely its cutoff.");
	}

	/** Assembly-only aliases and arbitrary prefixed controls cannot alter preparation. */
	private void testAssemblyOverrides() throws Exception{
		final TadPipe baseline=new TadPipe(config("baseline.fa"), System.err);
		baseline.work=Files.createDirectory(root.resolve("override-plan"));
		final ArrayList<TadPipe.Stage> original=baseline.plan();
		for(String prefix : Arrays.asList("", "assemble_")){
			final TadPipe.Config config=config("override.fa", prefix+"k=96,124,64,32,96",
					prefix+"assemblek=96", prefix+"graphk=64", prefix+"bridgek=124,64",
					prefix+"fusek=64,32", "assemble_mincontig=750", "assemble_mcs=2", "assemble_mce=1");
			config.validate();
			final TadPipe pipe=new TadPipe(config, System.err);
			pipe.work=baseline.work;
			final ArrayList<TadPipe.Stage> changed=pipe.plan();
			check(changed.size()==original.size(), "Assembly overrides must not skip or add preparation stages.");
			for(int i=0; i<changed.size()-1; i++){
				check(changed.get(i).command.equals(original.get(i).command),
						"Final assembly overrides must leave every intermediate command unchanged: "+changed.get(i).name);
			}
			final TadPipe.Stage assembly=changed.get(changed.size()-1);
			for(String arg : Arrays.asList("k=96,124,64,32,96", "assemblek=96", "graphk=64",
					"bridgek=124,64", "fusek=64,32", "mincontig=750", "mcs=2", "mce=1")){
				check(assembly.command.contains(arg), "The final child must receive this override without its prefix: "+arg);
			}
		}
		final TadPipe.Config ordered=config("ordered.fa", "assemble_k=124,64", "k=96,32",
				"assemblek=124", "assemble_assemblek=96", "ecc_k=55");
		final TadPipe pipe=new TadPipe(ordered, System.err);
		pipe.work=baseline.work;
		final ArrayList<TadPipe.Stage> plan=pipe.plan();
		final TadPipe.Stage assembly=plan.get(plan.size()-1);
		check(assembly.command.contains("k=96,32") && assembly.command.contains("assemblek=96"),
				"Bare and prefixed aliases share CLI order; the last value wins in either spelling.");
		check(plan.get(3).command.contains("k=55"), "An explicit ecc_k override still configures correction independently.");
		check(plan.get(10).command.contains("k=145") && plan.get(11).command.contains("k=124"),
				"Neither assembly nor ECC kmer overrides may leak into read extension.");
	}

	/** Confirms child isolation, publication, quoting, and successful FASTQ-only cleanup. */
	private void testSuccess() throws Exception{
		final TadPipe.Config config=config("success.fa", "dedupe_testlabel=a b=$literal",
				"k=96,124,64,32", "assemblek=96", "assemble_mincontig=750");
		final TadPipe pipe=new TadPipe(config, System.err);
		pipe.process();
		check(Files.size(config.out)>0, "Successful assembly must be published.");
		final List<String> pids=Files.readAllLines(pipe.work.resolve("pids.txt"), StandardCharsets.UTF_8);
		check(pids.size()==13 && new HashSet<String>(pids).size()==13, "Every phase must run in a fresh process.");
		check(Files.exists(pipe.work.resolve("commands.sh")), "Commands survive intermediate cleanup.");
		check(!Files.exists(pipe.work.resolve("paired.fq.gz")), "Successful cleanup removes only managed intermediate reads.");
		check(new String(Files.readAllBytes(pipe.work.resolve("01_dedupe.log")), StandardCharsets.UTF_8)
				.contains("ARG:<testlabel=a b=$literal>"), "Subprocess arguments must not undergo shell expansion or word splitting.");
		check(new String(Files.readAllBytes(in2), StandardCharsets.UTF_8).equals("input two\n"), "Mate input must remain untouched.");
		final String assemblyLog=new String(Files.readAllBytes(pipe.work.resolve("10_assemble.log")), StandardCharsets.UTF_8);
		check(assemblyLog.contains("ARG:<k=96,124,64,32>") && assemblyLog.contains("ARG:<assemblek=96>")
				&& assemblyLog.contains("ARG:<mincontig=750>"), "The launched final process must receive assembly overrides.");
		final String eccLog=new String(Files.readAllBytes(pipe.work.resolve("04_ecc.log")), StandardCharsets.UTF_8);
		check(eccLog.contains("ARG:<k=62>") && !eccLog.contains("assemblek=") && !eccLog.contains("mincontig="),
				"The launched correction process must keep its own parameters.");
		final TadPipe.Config dry=config("dry.fa", "dryrun=t");
		final TadPipe dryPipe=new TadPipe(dry, System.err);
		dryPipe.process();
		check(!Files.exists(dry.out) && !Files.exists(dryPipe.work.resolve("pids.txt")), "Dry run must not launch children or publish output.");
	}

	/** A nonzero exit or missing output must stop downstream work and preserve inputs. */
	private void testFailures() throws Exception{
		for(String flag : Arrays.asList("ecco_testfail=t", "ecco_omit=t")){
			final TadPipe.Config config=config(flag.contains("testfail") ? "failed.fa" : "missing.fa", flag);
			final TadPipe pipe=new TadPipe(config, System.err);
			boolean failed=false;
			try{pipe.process();}catch(IOException expected){failed=true;}
			check(failed && !Files.exists(config.out), "Failed child output cannot be presented as a completed assembly.");
			check(Files.exists(pipe.work.resolve("paired.fq.gz")), "Failure must retain intermediate reads for diagnosis.");
			check(!Files.exists(pipe.work.resolve("03_clump.log")), "Downstream consumers must not run after failure.");
		}
	}

	/** Validates input/output aliases and rejects stage options that steal pipeline filenames. */
	private void testValidation() throws Exception{
		final TadPipe.Config split=config("split.fa", "interleaved=f");
		split.validate();
		check(Boolean.FALSE.equals(split.interleaved) && split.in2!=null, "Explicit interleaved=f is valid with two mate files.");
		for(String flag : Arrays.asList("ecco_in=other.fq", "assemble_out=x.fa", "merge_config=hidden.config", "filter_ref=phix")){
			boolean failed=false;
			try{config("bad.fa", flag);}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid managed-path or obsolete phase option must fail before launch: "+flag);
		}
		final TadPipe.Config alias=config("unused.fa", "overwrite=t");
		alias.out=in2;
		boolean failed=false;
		try{alias.validate();}catch(IOException expected){failed=true;}
		check(failed, "Output may not overwrite either input, including read two.");
		final TadPipe.Config custom=config("custom.fa", "assemble_fusenet="+tools.resolve("networks/tadpole_fusion.bbnet"));
		failed=false;
		try{custom.validate();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A custom model requires an explicitly paired cutoff.");
		final TadPipe.Config modelAlias=config("unused.fa", "overwrite=t");
		modelAlias.out=tools.resolve("networks/tadpole_fusion.bbnet");
		failed=false;
		try{new TadPipe(modelAlias, System.err).process();}catch(IOException expected){failed=true;}
		check(failed, "Publication must not replace the model needed by the assembly.");
	}

	/** A valid compressed empty file must not pass the final nonempty-assembly gate. */
	private void testEmptyCompressed() throws Exception{
		final Path empty=root.resolve("empty.fa.gz");
		final ByteStreamWriter writer=new ByteStreamWriter(empty.toString(), false, false, true);
		writer.start();
		check(!writer.poisonAndWait(), "Native compressed fixture writer must finish without an I/O error.");
		check(Files.size(empty)>0 && !TadPipe.hasSequence(empty), "Gzip/BGZF header bytes alone cannot count as assembled sequence.");
	}

	/** Supplies small explicit resources so these tests never launch computational tools. */
	private TadPipe.Config config(final String output, final String... more){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+in, "in2="+in2,
				"out="+root.resolve(output), "bbtools="+tools, "temp="+root, "t=1", "childheap=256m"));
		args.addAll(Arrays.asList(more));
		return new TadPipe.Config(args.toArray(new String[args.size()]));
	}

	/** Counts assertions while keeping a useful explanation at every call site. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private final Path root, tools, in, in2;
	private int checks;
}
