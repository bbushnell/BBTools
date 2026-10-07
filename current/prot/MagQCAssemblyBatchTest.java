package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;

import ml.CellNet;
import ml.CellNetParser;

/** Small public input/configuration fixtures, without network or biological model access.
 * @author Yoimiya */
public final class MagQCAssemblyBatchTest {

	/** Checks deterministic inputs, relocatable resource names and guarded public options. */
	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Input-contract fixture accepts no arguments");}
		checkCallerMode();
		checkSearchOptions();
		checkMissingAssets();
		MagQCAssemblyReportTest.main(new String[0]);
		check(!MagQCAssemblyBatch.lockedMode("worker") && MagQCAssemblyBatch.lockedMode("locked"), "D240 comparison modes");
		reject("compositemode=automatic"); reject("subnetmode=shared");
		check(!MagQCAssemblyBatch.parallelMode("serial") && MagQCAssemblyBatch.parallelMode("parallel"), "Explicit resource loading modes");
		reject("loadmode=automatic");
		final Path root=Files.createTempDirectory("prokcc-input-contract-");
		final Path first=root.resolve("a.fa"), second=root.resolve("b.fna"), config=root.resolve("release.config");
		try{
			Files.write(first, ">a\nACGT\n".getBytes(StandardCharsets.US_ASCII));
			Files.write(second, ">b\nACGT\n".getBytes(StandardCharsets.US_ASCII));
			Files.write(config, "bundle=resources/subnets.bbnets\ncomperrormultiplier=2.0\n".getBytes(StandardCharsets.US_ASCII));
			final ArrayList<MagQCAssemblyBatch.Job> ordered=MagQCAssemblyBatch.inputs(root.toString());
			check(ordered.size()==2 && ordered.get(0).path.equals(first.toString()) && ordered.get(1).path.equals(second.toString()),
				"Directory enumeration must retain only FASTA inputs in stable order");
			final ArrayList<MagQCAssemblyBatch.Job> explicit=MagQCAssemblyBatch.inputs(second+","+first);
			check(explicit.get(0).path.equals(second.toString()), "An explicit batch must preserve requested order");
			boolean rejected=false;
			try{MagQCAssemblyBatch.inputs(first+","+first);}catch(IllegalArgumentException e){rejected=true;}
			check(rejected, "A repeated assembly must not produce duplicate bin rows");
			final HashMap<String,String> parsed=MagQCAssemblyBatch.parseOptions(new String[]{"config="+config, "in="+first, "out=stdout"});
			check(parsed.get("bundle").equals(root.resolve("resources/subnets.bbnets").toString()), "Relocated release config lost its resource base");
			check(parsed.get("in").equals(first.toString()), "Release config must not relocate user inputs");
			final HashMap<String,String> positional=MagQCAssemblyBatch.parseOptions(new String[]{"config="+config, first.toString()});
			check(!positional.containsKey("out") && parsed.get("out").equals("stdout"),
				"Omitted out disables data output; explicit out=stdout retains pipeable TSV output");
			final HashMap<String,String> named=MagQCAssemblyBatch.parseOptions(new String[]{"config="+config, "in="+first});
			check(positional.equals(named), "Bare input after the wrapper's config must equal explicit in=");
			check(new MagQCAssemblyBatch(new String[]{first.toString()}).jobs.size()==1
				&& new MagQCAssemblyBatch(new String[]{"in="+first}).jobs.size()==1,
				"A single input must not require out= in either public syntax");
			for(String[] duplicate:new String[][]{{first.toString(), second.toString()},
				{first.toString(), "in="+second}, {"in="+first, second.toString()}}){
				rejected=false;
				try{MagQCAssemblyBatch.parseOptions(duplicate);}catch(IllegalArgumentException e){rejected=true;}
				check(rejected, "A positional argument must not silently replace or append another input");
			}
			reject(root.resolve("missing.fa").toString());
			check(MagQCAssemblyBatch.parseOptions(new String[]{"swapln=t"}).get("swapnl").equals("t")
				&& MagQCAssemblyBatch.parseOptions(new String[0]).get("swapnl").equals("f"),
				"The swapln alias and BBTools default must select the same N/L convention as stats.sh");
			rejected=false;
			try{MagQCAssemblyBatch.parseOptions(new String[]{"swapnl=t", "swapln=f"});}catch(IllegalArgumentException e){rejected=true;}
			check(rejected, "Conflicting N/L aliases must not silently override each other");
			check(parsed.get("compositemode").equals("locked") && parsed.get("subnetmode").equals("locked"), "Measured public defaults must share composite and subnet inference");
			check(parsed.get("policy").equals("BOUNDED_LOOKAHEAD") && parsed.get("lookahead").equals("4"), "Public release assignment defaults changed");
			check(parsed.get("errorfitset").equals("UNCALIBRATED") && parsed.get("contamerrormultiplier").equals("1.0"), "Error calibration was silently invented");
			check(MagQCAssemblyBatch.multiplier(parsed.get("comperrormultiplier"))==2.0, "Explicit error factor lost");
			for(String invalid:new String[]{"0", "-1", "NaN", "Infinity"}){
				rejected=false;
				try{MagQCAssemblyBatch.multiplier(invalid);}catch(IllegalArgumentException e){rejected=true;}
				check(rejected, "Invalid error factor accepted: "+invalid);
			}
			reject("dummy=t"); reject("fasta=hidden.fa"); reject("unknown=value"); reject("taxphylum=a\tb");
			rejected=false;
			try{MagQCAssemblyBatch.parseOptions(new String[]{"t=1", "threads=2"});}catch(IllegalArgumentException e){rejected=true;}
			check(rejected, "Aliases must not hide duplicate worker settings");
			checkPartialStart(first, root.resolve("unused.tsv"));
			final Path occupied=root.resolve("existing.tsv");
			Files.write(occupied, "preserve".getBytes(StandardCharsets.US_ASCII));
			try{
				reject("outstream="+occupied);
				check(new String(Files.readAllBytes(occupied), StandardCharsets.US_ASCII).equals("preserve"),
					"Unsupported outstream must not truncate a file before public option validation");
				Files.write(config, ("outstream="+occupied+"\n").getBytes(StandardCharsets.US_ASCII));
				reject("config="+config);
				check(new String(Files.readAllBytes(occupied), StandardCharsets.US_ASCII).equals("preserve"),
					"Config-expanded outstream must be rejected without opening its destination");
				reject("metadatafile="+occupied);
				reject("proxyhost=unsupported");
				reject("bufferbf=f");
				rejected=false;
				try{new MagQCAssemblyBatch(new String[]{"in="+first, "out="+occupied, "ow=f"});}catch(IllegalArgumentException e){rejected=true;}
				check(rejected && new String(Files.readAllBytes(occupied), StandardCharsets.US_ASCII).equals("preserve"),
					"ow=f must reject existing reports before resource loading and leave them untouched");
				new MagQCAssemblyBatch(new String[]{"in="+first, "out="+occupied});
				new MagQCAssemblyBatch(new String[]{"in="+first, "out=NuLl"});
				check(new String(Files.readAllBytes(occupied), StandardCharsets.US_ASCII).equals("preserve"),
					"Default overwrite is allowed but output opens only after every bin succeeds");
				check(MagQCAssemblyBatch.parseOptions(new String[]{"overwrite=f"}).get("ow").equals("f"),
					"The standard overwrite alias must retain explicit false");
				rejected=false;
				try{new MagQCAssemblyBatch(new String[]{"in="+first, "out="+first});}catch(RuntimeException e){rejected=true;}
				check(rejected, "Standard duplicate-file checks must prevent overwriting the input assembly");
			}finally{Files.deleteIfExists(occupied);}
			System.out.println("MagQCAssemblyBatchTest PASS: input order, duplicate rejection, relocatable config, raw-error defaults and standard output handling");
		}finally{
			Files.deleteIfExists(first); Files.deleteIfExists(second); Files.deleteIfExists(config); Files.deleteIfExists(root);
		}
	}

	/** Public search defaults and the explicit legacy switch survive option parsing. */
	private static void checkSearchOptions(){
		final HashMap<String,String> normal=MagQCAssemblyBatch.parseOptions(new String[0]);
		final HashMap<String,String> legacy=MagQCAssemblyBatch.parseOptions(new String[]{"normalsearch=f"});
		check(MagQCAssemblyInput.normalSearch(normal) && !MagQCAssemblyInput.normalSearch(legacy),
			"Public taxonomy must default to normal search and retain an explicit legacy switch");
		boolean rejected=false;
		try{MagQCAssemblyBatch.parseOptions(new String[]{"normalsearch=t", "normalsearch=f"});}
		catch(IllegalArgumentException e){rejected=true;}
		check(rejected, "Duplicate search flags must not silently change inference behavior");
		reject("normalsearch=automatic"); reject("normalsearch=treu");
		final String row="#Query1\nmagqc_bin\t1\t10\t1\thit\t123\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\tp__Bacillota\t0\t0\n";
		check(MagQCAssemblyInput.parseResponse(row, 10, 1, true).status.equals("unknown"),
			"Local whole-bin parsing must retain unknown when a display lineage lacks its domain");
	}

	/** Missing optional downloads are reported together, without opening outputs or contacting services. */
	private static void checkMissingAssets() throws Exception{
		final Path root=Files.createTempDirectory("prokcc-missing-assets-");
		final HashMap<String,String> options=new HashMap<String,String>();
		options.put("net", root.resolve("composite_v1.2.1_shrunk_0.8pct_18bit.bbnet.gz").toString());
		options.put("bundle", root.resolve("magqc_subnets_v1.2.1_shrunk_57.4pct_18bit.bbnets.gz").toString());
		options.put("hbmbundle", root.resolve("magqc_hbm_v1.rare01.hbmt.gz").toString());
		options.put("sidecar", root.resolve("magqc_sidecar_v1.tsv.gz").toString());
		boolean rejected=false;
		try{MagQCAssemblyBatch.requireOptionalAssets(options);}
		catch(IllegalArgumentException e){
			final String message=e.getMessage();
			rejected=message.contains("composite_v1.2.1_shrunk_0.8pct_18bit.bbnet.gz") && message.contains("magqc_subnets_v1.2.1_shrunk_57.4pct_18bit.bbnets.gz")
				&& message.contains("magqc_hbm_v1.rare01.hbmt.gz") && message.contains("magqc_sidecar_v1.tsv.gz")
				&& message.contains("https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.2.1.tar")
				&& message.contains("Extract its contents into BBTools/resources/ to create resources/prokcc/")
				&& !message.contains("PROKCC_V1_TAG_PENDING");
		}
		check(rejected, "Missing assets must identify all four files and their download instructions");
		final String archive="https://sourceforge.net/projects/bbmap/files/Resources/prokcc_v1.2.1.tar";
		check(archive.equals(shared.Resources.downloadURL("?prokcc/release.config"))
			&& archive.equals(shared.Resources.downloadURL("magqc_subnets_v1.2.1_shrunk_57.4pct_18bit.bbnets.gz")),
			"The release config and full subnet bundle must resolve to the same manual archive");
		check("https://sourceforge.net/projects/bbmap/files/Resources/".equals(shared.Resources.downloadURL("ssuSketchDDL.tsv.gz")),
			"Unrelated resources must retain their existing download location");
		for(String path:options.values()){Files.write(java.nio.file.Paths.get(path), new byte[]{1});}
		MagQCAssemblyBatch.requireOptionalAssets(options);
		for(String path:options.values()){Files.delete(java.nio.file.Paths.get(path));}
		Files.delete(root);
	}

	/** Model loading must not change the mode used by a previously loaded legacy caller net. */
	private static void checkCallerMode() throws Exception{
		final boolean original=CellNet.DENSE;
		final String header="##bbnet\n#version 1\n#concise\n#coding decimal\n#seed 23\n#layers 2\n#blocksize 1\n#dims 3 1\n";
		final String dense=header+"#dense\nC4 LINEAR 0.25 0.5 0 -0.25\n";
		final String sparse=header+"#sparse\nI4 0 2\nC4 LINEAR 0.25 0.5 -0.25\n";
		final float[] input={2, 1, 4};
		try{
			for(boolean callerDense:new boolean[]{true, false}){
				final ArrayList<byte[]> lines=new ArrayList<byte[]>();
				for(String line:(callerDense ? dense : sparse).split("\n")){lines.add(line.getBytes(StandardCharsets.US_ASCII));}
				final CellNet caller=CellNetParser.loadFromLines(lines);
				caller.applyInput(input); final float expected=caller.feedForward();
				final MagQCNetBundle.InferenceModel model=MagQCNetBundle.parseInferenceModel(
					(callerDense ? sparse : dense).getBytes(StandardCharsets.US_ASCII), "opposite-mode fixture");
				check(CellNet.DENSE==callerDense, "Composite parsing changed the legacy caller mode");
				final CellNet copy=caller.copy(false), worker=model.newWorker();
				copy.applyInput(input); worker.applyInput(input);
				check(copy.feedForward()==expected && worker.feedForward()==expected,
					"Legacy copy or instance-specific inference changed after an opposite-mode load");
				boolean rejected=false;
				try{
					MagQCNetBundle.parseInferenceModel((header+(callerDense ? "#sparse\nI99 0 2\n" : "#dense\n")+
						"C99 LINEAR 0.25 0.5 0 -0.25\n").getBytes(StandardCharsets.US_ASCII), "broken-mode fixture");
				}catch(java.io.IOException e){rejected=true;}
				check(rejected && CellNet.DENSE==callerDense, "Failed master parsing leaked its temporary network mode");
			}
		}finally{CellNet.DENSE=original;}
	}

	/** Injects failure on the second start while the first worker is still using its session. */
	private static void checkPartialStart(Path input, Path output) throws Exception{
		final MagQCAssemblyBatch batch=new MagQCAssemblyBatch(new String[]{"in="+input, "out="+output});
		final CountDownLatch entered=new CountDownLatch(1), release=new CountDownLatch(1);
		final AtomicBoolean finished=new AtomicBoolean();
		final Thread first=new Thread(){
			@Override public void run(){
				entered.countDown();
				try{release.await(); Thread.sleep(50); finished.set(true);}
				catch(InterruptedException e){throw new AssertionError(e);}
			}
		};
		final Thread second=new Thread(){
			@Override public synchronized void start(){
				try{entered.await();}catch(InterruptedException e){throw new AssertionError(e);}
				release.countDown();
				throw new OutOfMemoryError("Injected native-thread creation failure");
			}
		};
		final ArrayList<Thread> workers=new ArrayList<Thread>();
		workers.add(first); workers.add(second);
		boolean rejected=false;
		try{batch.startAndJoin(workers);}catch(RuntimeException e){rejected=e.getCause() instanceof OutOfMemoryError;}
		check(rejected && finished.get() && first.getState()==Thread.State.TERMINATED && second.getState()==Thread.State.NEW,
			"Partial thread startup must join every started worker before propagating failure");
	}

	/** Unknown/developer-only public arguments fail before resource or service access. */
	private static void reject(String arg){
		boolean rejected=false;
		try{MagQCAssemblyBatch.parseOptions(new String[]{arg});}catch(IllegalArgumentException e){rejected=true;}
		check(rejected, "Invalid public option accepted: "+arg);
	}

	/** Fixture checks remain active when testing caller rejection without JVM assertions. */
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
}
