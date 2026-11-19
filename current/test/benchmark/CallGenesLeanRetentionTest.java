package test.benchmark;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.CountDownLatch;
import fileIO.ReadWrite;

import benchmark.CallGenesLeanCalls;
import benchmark.CallGenesLeanReplayBenchmark;
import benchmark.CallGenesRetainedExport;
import benchmark.CallGenesBenchmarkSupport;
import benchmark.CallGenesLengthScorecard;
import benchmark.CallGenesStopAudit;
import gff.GffLine;
import fileIO.ByteFile;
import prot.DigestSuffix;
import shared.Shared;
import structures.ByteBuilder;
import prok.Orf;
import prok.ProkObject;

/** Actual post-worker grading with retention off/on; no caller or classifier run.
 * Reflection supplies completed compact results without widening production APIs.
 * @author Sayu
 */
public final class CallGenesLeanRetentionTest{
	/** Checks synthetic retained-output parity, or runs the explicitly supplied saved-cohort export comparison. */
	public static void main(String[] args)throws Exception{
		check(args.length>=3 && args.length<=6, "Expected fresh fixture directory, bgzip path, pigz path and optional full983 inventory; legacy pin arguments are ignored");
		final boolean legacy=args.length>=5;
		compressorOptions=new String[]{"retentionbgzip="+args[1], "retentionpigz="+args[legacy ? 3 : 2]};
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute(), "Fixture directory must be absolute");
		Files.createDirectory(root);
		if(args.length==4 || args.length==6){
			fullCohort(root, Paths.get(args[legacy ? 5 : 3]));
			return;
		}
		new CallGenesRetainedExport.Tools(Paths.get(args[1]), "ignored", Paths.get(args[legacy ? 3 : 2]), "ignored");
		final Path reference=root.resolve("reference.gff"), loci=root.resolve("loci.tsv");
		Files.write(reference, ("c\tfixture\tCDS\t10\t132\t.\t+\t0\t.\n"+
			"c\tfixture\tCDS\t400\t597\t.\t-\t0\t.\n"+
			"c\tfixture\tCDS\t3901\t3999\t.\t+\t0\tpartial=01\n").getBytes(StandardCharsets.US_ASCII));
		loci(loci);
		final boolean callCds=ProkObject.callCDS;
		ProkObject.callCDS=true;
		try{
			final Path off=Files.createDirectory(root.resolve("off")), on=Files.createDirectory(root.resolve("on"));
			final String flags=compressionFlags();
			grade(root, reference, loci, off, null);
			grade(root, reference, loci, on, "true");
			check(flags.equals(compressionFlags()), "Small retention changed global compression settings");
			for(String name:new String[]{"scorecard.tsv", "by_genome.tsv", "features.tsv"}){
				check(Arrays.equals(Files.readAllBytes(off.resolve(name)), Files.readAllBytes(on.resolve(name))),
						"Retention altered actual grading/count output: "+name);
			}
			check(!Files.exists(off.resolve("retained_gffs")) && !Files.exists(off.resolve("retained_gffs.tsv")),
					"Default retention must create no GFF output");
			final Path nativeGff=on.resolve("retained_gffs/g0001.native.gff.gz"), formulaGff=on.resolve("retained_gffs/g0001.formula.gff.gz");
			final String inventory="genome\tnative_gff\tformula_gff\ng0001\t"+nativeGff+"\t"+formulaGff+"\n";
			check(inventory.equals(new String(Files.readAllBytes(on.resolve("retained_gffs.tsv")), StandardCharsets.US_ASCII)),
					"Retained inventory must preserve the accepted genome/native/formula schema");
			checkCosts(on, false);
			final Path parallel=Files.createDirectory(root.resolve("parallel"));
			grade(root, reference, loci, parallel, "true", "retentionworkers=2", "retainrna=true");
			for(String name:new String[]{"scorecard.tsv", "by_genome.tsv", "features.tsv"}){
				check(Arrays.equals(Files.readAllBytes(off.resolve(name)), Files.readAllBytes(parallel.resolve(name))), "Parallel/RNA retention changed all16 grading: "+name);
			}
			for(String arm:new String[]{"native", "formula"}){
				check(Arrays.equals(decoded(on.resolve("retained_gffs/g0001."+arm+".gff.gz")), decoded(parallel.resolve("retained_gffs/g0001."+arm+".gff.gz"))), "Parallel CDS bytes differ: "+arm);
				check(Arrays.equals(rnaBytes(calls()), decoded(parallel.resolve("retained_rna_gffs/g0001."+arm+".gff.gz"))), "RNA capture changed canonical attributes/order: "+arm);
			}
			checkCosts(parallel, false);
			checkRna(parallel, 2);
			final Path serialRna=Files.createDirectory(root.resolve("serial_rna"));
			grade(root, reference, loci, serialRna, "true", "retainrna=true");
			checkRna(serialRna, 2);
			for(String name:new String[]{"scorecard.tsv", "by_genome.tsv", "features.tsv"}){
				check(Arrays.equals(Files.readAllBytes(off.resolve(name)), Files.readAllBytes(serialRna.resolve(name))), "Serial RNA retention changed all16 grading: "+name);
			}
			for(String arm:new String[]{"native", "formula"}){
				check(Arrays.equals(rnaBytes(calls()), decoded(serialRna.resolve("retained_rna_gffs/g0001."+arm+".gff.gz"))), "Serial RNA bytes differ");
			}
			check(!Files.exists(on.resolve("retained_rna_gffs")), "RNA retention must remain opt-in");
			check(flags.equals(compressionFlags()), "Parallel workers changed global compression flags");
			workerFailures(root);
			final ArrayList<GffLine> exported=new ArrayList<GffLine>();
			final ByteFile reader=ByteFile.makeByteFile1(nativeGff.toString(), false);
			try{for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				if(line.length>0 && line[0]!='#'){exported.add(new GffLine(line));}
			}}finally{check(!reader.close(), "Retained gzip reader failed");}
			check(exported.size()==5 && exported.get(0).start==exported.get(1).start
					&& exported.get(0).stop==exported.get(1).stop, "Retained GFF dropped/reordered a duplicate or CDS");
			check(exported.get(4).attributes.equals("partial=01"), "Retained GFF lost physical partial provenance");
			final CallGenesLengthScorecard.ReferenceBands bands=CallGenesLengthScorecard.readSavedBands(loci, Collections.singleton("g0001"))
					.bind("g0001", CallGenesStopAudit.loadReference(reference.toString()));
			final CallGenesLeanCalls original=new CallGenesLeanCalls();
			original.append("c", calls());
			check(Arrays.equals(bands.grade(exported).rows("native", "g0001").toBytes(),
					bands.grade(original.lines()).rows("native", "g0001").toBytes()), "Export/reload changed any of the16 score fields");
			final Path collision=Files.createDirectory(root.resolve("collision"));
			Files.createDirectory(collision.resolve("retained_gffs"));
			boolean refused=false;
			try{grade(root, reference, loci, collision, "true");}
			catch(InvocationTargetException failure){
				check(failure.getCause() instanceof java.io.UncheckedIOException, "Wrong output failure");
				refused=true;
			}
			check(refused && !Files.exists(collision.resolve("retained_gffs.tsv")), "Output collision must refuse without a complete inventory");
			boolean invalid=false;
			try{instance(root, root.resolve("invalid"), "maybe");}
			catch(InvocationTargetException failure){
				check(failure.getCause() instanceof IllegalArgumentException, "Wrong option refusal");
				invalid=true;
			}
			check(invalid, "Invalid retention selector was accepted");
			for(String[] options:new String[][]{{"false", "retainrna=true"}, {"true", "retainrna=maybe"},
				{"true", "retentionworkers=0"}, {"true", "retentionworkers=3"}, {"false", "retentionworkers=2"}}){
				invalid=false;
				try{instance(root, root.resolve("invalid_extra"), options[0], options[1]);}
				catch(InvocationTargetException failure){
					check(failure.getCause() instanceof IllegalArgumentException, "Wrong RNA/worker option refusal");
					invalid=true;
				}
				check(invalid, "Invalid RNA/worker combination accepted: "+Arrays.toString(options));
			}
			//The same actual grading path exercises the large BGZF branch without biological work.
			repetitions=2000;
			final Path largeOff=Files.createDirectory(root.resolve("large_off")), largeOn=Files.createDirectory(root.resolve("large_on"));
			grade(root, reference, loci, largeOff, null);
			grade(root, reference, loci, largeOn, "true");
			for(String name:new String[]{"scorecard.tsv", "by_genome.tsv", "features.tsv"}){
				check(Arrays.equals(Files.readAllBytes(largeOff.resolve(name)), Files.readAllBytes(largeOn.resolve(name))),
						"BGZF retention changed any count output: "+name);
			}
			checkCosts(largeOn, true);
			final Path largeParallel=Files.createDirectory(root.resolve("large_parallel"));
			grade(root, reference, loci, largeParallel, "true", "retentionworkers=2", "retainrna=true");
			checkCosts(largeParallel, true);
			checkRna(largeParallel, 2L*repetitions);
			for(String name:new String[]{"scorecard.tsv", "by_genome.tsv", "features.tsv"}){
				check(Arrays.equals(Files.readAllBytes(largeOff.resolve(name)), Files.readAllBytes(largeParallel.resolve(name))), "Large parallel/RNA retention changed all16 grading: "+name);
			}
			for(String arm:new String[]{"native", "formula"}){
				check(Arrays.equals(decoded(largeOn.resolve("retained_gffs/g0001."+arm+".gff.gz")), decoded(largeParallel.resolve("retained_gffs/g0001."+arm+".gff.gz"))), "Large parallel CDS bytes differ");
				check(Arrays.equals(rnaBytes(calls()), decoded(largeParallel.resolve("retained_rna_gffs/g0001."+arm+".gff.gz"))), "Large RNA bytes/duplicates differ");
			}
			check(flags.equals(compressionFlags()), "BGZF retention changed global compression settings");
			repetitions=1;
			Files.write(root.resolve("PASS"), "LEAN_RETAINED_GFF_FIXTURE_PASS output_on_off_all16=true caller_runs=0\n".getBytes(StandardCharsets.US_ASCII));
			System.out.println("LEAN_RETAINED_GFF_FIXTURE_PASS output_on_off_all16=true caller_runs=0");
		}finally{
			ProkObject.callCDS=callCds;
			repetitions=1;
		}
	}

	/** Inject completed worker records, then exercise the actual unchanged grading barrier. */
	private static void grade(Path root, Path reference, Path loci, Path output, String retention, String... extra)throws Exception{
		final Object benchmark=instance(root, output, retention, extra);
		final Class<?> genomeType=Class.forName("benchmark.CallGenesLeanReplayBenchmark$Genome");
		final Constructor<?> genomeCtor=genomeType.getDeclaredConstructor(String.class, Path.class, String.class, Path.class, String.class, Path.class);
		genomeCtor.setAccessible(true);
		final ArrayList<Object> genomes=new ArrayList<Object>();
		genomes.add(genomeCtor.newInstance("g0001", root.resolve("unused.fna"), "unmeasured", reference, "unmeasured", root.resolve("unused.annotations.tsv")));
		final Class<?> resultType=Class.forName("benchmark.CallGenesLeanReplayBenchmark$Result");
		final Constructor<?> resultCtor=resultType.getDeclaredConstructor(String.class);
		resultCtor.setAccessible(true);
		final Object result=resultCtor.newInstance("g0001");
		final Method capture=CallGenesLeanReplayBenchmark.class.getDeclaredMethod("captureCalls", resultType, int.class, String.class, ArrayList.class);
		capture.setAccessible(true);
		for(int arm=0; arm<4; arm++){capture.invoke(benchmark, result, arm, "c", calls());}
		final Field enabled=CallGenesLeanReplayBenchmark.class.getDeclaredField("retainRna"), rna=resultType.getDeclaredField("rna");
		enabled.setAccessible(true);
		rna.setAccessible(true);
		check((rna.get(result)!=null)==enabled.getBoolean(benchmark), "RNA capture must allocate only when explicitly enabled");
		final ArrayList<Object> results=new ArrayList<Object>();
		results.add(result);
		final Method grade=CallGenesLeanReplayBenchmark.class.getDeclaredMethod("grade", ArrayList.class, ArrayList.class, CallGenesLengthScorecard.SavedBands.class);
		grade.setAccessible(true);
		grade.invoke(benchmark, genomes, results, CallGenesLengthScorecard.readSavedBands(loci, Collections.singleton("g0001")));
		check(results.get(0)==null, "Grading failed to release completed compact results");
	}

	/** Construction validates options but performs no initialization or biological input read. */
	private static Object instance(Path root, Path output, String retention, String... extra)throws Exception{
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList("in="+root.resolve("unused.in"),
			"loci="+root.resolve("loci.tsv"), "locisha80=unmeasured", "pgm="+root.resolve("unused.pgm"), "pgmpolicy=fixed",
			"tailreference="+root.resolve("unused.resources"), "tailreferencesha80=unmeasured", "fmnroot="+root.resolve("unused.fmn"),
			"out="+output, "genomes=1", "workers=1", "searchthreads=1", "cpus=2", "form=BLOSUMratio",
			"outer=0.6", "scale=0.525", "a=-30", "b=0.1125", "cap=400", "i0=4", "signal=blosum_raw"));
		if(retention!=null){options.add("retaingffs="+retention);}
		if("true".equals(retention)){options.addAll(Arrays.asList(compressorOptions));}
		options.addAll(Arrays.asList(extra));
		final Constructor<?> ctor=CallGenesLeanReplayBenchmark.class.getDeclaredConstructor(String[].class);
		ctor.setAccessible(true);
		return ctor.newInstance((Object)options.toArray(new String[0]));
	}

	/** Ordered duplicate, plus/minus matches, internal extra, physical partial and RNA. */
	private static ArrayList<Orf> calls(){
		final byte[] bases=new byte[4000];
		Arrays.fill(bases, (byte)'A');
		final ArrayList<Orf> result=new ArrayList<Orf>();
		result.add(new Orf("c", 9, 131, 0, 0, bases, false, ProkObject.CDS, false, false));
		result.add(new Orf("c", 9, 131, 0, 0, bases, false, ProkObject.CDS, false, false));
		result.add(new Orf("c", 399, 596, 1, 0, bases, false, ProkObject.CDS, false, false));
		result.add(new Orf("c", 3299, 3397, 0, 0, bases, false, ProkObject.CDS, false, false));
		result.add(new Orf("c", 3900, 3998, 0, 0, bases, false, ProkObject.CDS, false, true));
		result.add(new Orf("c", 2000, 2200, 0, 0, bases, false, ProkObject.r16S, false, false));
		final Orf trna=new Orf("c", 2400, 2474, 1, 0, bases, false, ProkObject.tRNA, true, false);
		trna.trnaModel="fixture_model";
		trna.trnaAnticodon="CAT";
		trna.ncrnaFamily="tmrna";
		result.add(trna);
		if(repetitions>1){
			final ArrayList<Orf> one=new ArrayList<Orf>(result);
			for(int i=1; i<repetitions; i++){result.addAll(one);}
		}
		return result;
	}

	/** Independent canonical serialization of original emitted non-CDS features. */
	private static byte[] rnaBytes(ArrayList<Orf> calls){
		final ByteBuilder expected=new ByteBuilder("##gff-version 3\n");
		for(Orf call:calls){
			if(call.type!=ProkObject.CDS){
				call.appendGff(expected);
				expected.nl();
			}
		}
		return expected.toBytes();
	}
	/** Decode concatenated gzip/BGZF members without normalizing any original byte. */
	private static byte[] decoded(Path file)throws java.io.IOException{
		final ByteBuilder bytes=new ByteBuilder();
		try(java.io.InputStream input=new java.util.zip.GZIPInputStream(Files.newInputStream(file))){
			final byte[] buffer=new byte[65536];
			for(int count=input.read(buffer); count>=0; count=input.read(buffer)){if(count>0){bytes.append(buffer, 0, count);}}
		}
		return bytes.toBytes();
	}
	/** Reconcile the separate two-arm inventory and row/byte provenance with actual RNA files. */
	private static void checkRna(Path output, long rows)throws Exception{
		final java.util.List<String> inventory=Files.readAllLines(output.resolve("retained_rna_gffs.tsv"), StandardCharsets.US_ASCII);
		check(inventory.size()==2 && inventory.get(0).equals("genome\tnative_rna_gff\tformula_rna_gff"), "RNA inventory changed CDS schema or lost coverage");
		final java.util.List<String> provenance=Files.readAllLines(output.resolve("retained_rna_provenance.tsv"), StandardCharsets.US_ASCII);
		check(provenance.size()==3, "RNA provenance does not cover both emitted arms");
		for(int i=1; i<provenance.size(); i++){
			final String[] cells=provenance.get(i).split("\t", -1);
			check(cells.length==7 && cells[0].equals("g0001") && Long.parseLong(cells[2])==rows
				&& cells[5].equals("Orf.appendGff_before_LeanCalls") && cells[6].equals("all_emitted_non_CDS"), "RNA row/formatter provenance differs");
			final Path file=output.resolve("retained_rna_gffs/g0001."+cells[1]+".gff.gz");
			check(Long.parseLong(cells[3])==decoded(file).length && Long.parseLong(cells[4])==Files.size(file), "RNA byte provenance differs from actual payload");
		}
	}
	/** Exercise real process failure and duplicate claims; the dispatcher must join before throwing. */
	private static void workerFailures(Path root)throws Exception{
		final Path falseTool=Paths.get("/usr/bin/false").toRealPath();
		final CallGenesRetainedExport.Tools failing=new CallGenesRetainedExport.Tools(falseTool, falseTool);
		final Path failed=Files.createDirectory(root.resolve("failed_compressor"));
		boolean refused=false;
		try{CallGenesRetainedExport.write(2, 2, i->new CallGenesRetainedExport.Payload(failed, "g0001", i==0 ? "native" : "formula", new ByteBuilder("##gff-version 3\n"), 0), failing);}
		catch(RuntimeException expected){refused=true;}
		check(refused, "A nonzero compressor exit published successful receipts");
		final CallGenesRetainedExport.Tools tools=tools();
		final Path duplicate=Files.createDirectory(root.resolve("duplicate_output"));
		refused=false;
		try{CallGenesRetainedExport.write(2, 2, i->new CallGenesRetainedExport.Payload(duplicate, "g0001", "native", new ByteBuilder("##gff-version 3\n"), 0), tools);}
		catch(RuntimeException expected){refused=true;}
		check(refused, "Duplicate parallel output paths were silently overwritten");
		final Path delayed=root.resolve("delayed_compressor.sh");
		Files.write(delayed, "#!/bin/sh\nsleep 0.2\ncat\n".getBytes(StandardCharsets.US_ASCII));
		check(delayed.toFile().setExecutable(true), "Cannot prepare native join-on-failure fixture");
		final CallGenesRetainedExport.Tools delayTools=new CallGenesRetainedExport.Tools(delayed, delayed);
		final Path joined=Files.createDirectory(root.resolve("join_on_failure"));
		final CountDownLatch started=new CountDownLatch(2);
		refused=false;
		try{
			CallGenesRetainedExport.write(2, 2, i->{
				started.countDown();
				try{started.await();}catch(InterruptedException problem){
					Thread.currentThread().interrupt();
					throw new IllegalStateException(problem);
				}
				if(i==0){throw new IllegalArgumentException("deliberate payload failure");}
				return new CallGenesRetainedExport.Payload(joined, "g0001", "formula", new ByteBuilder("joined-worker\n"), 0);
			}, delayTools);
		}catch(RuntimeException expected){refused=true;}
		check(refused && Arrays.equals(Files.readAllBytes(joined.resolve("g0001.formula.gff.gz")), "joined-worker\n".getBytes(StandardCharsets.US_ASCII)),
			"A failed payload returned before the other active writer/process joined");
	}
	/** Rebind the fixture's explicit native compressor paths without executable stamps. */
	private static CallGenesRetainedExport.Tools tools(){
		return new CallGenesRetainedExport.Tools(Paths.get(compressorOptions[0].substring("retentionbgzip=".length())),
			Paths.get(compressorOptions[1].substring("retentionpigz=".length())));
	}

	/** Optional full983 I/O comparison: original final GFFs are loaded only once.
	 * Serial and parallel phases are sequential to avoid timing contention. No caller runs. */
	private static void fullCohort(Path root, Path inventory)throws Exception{
		check(inventory.isAbsolute() && Files.isRegularFile(inventory) && !Files.isSymbolicLink(inventory), "Full-cohort inventory must be an existing canonical file");
		Shared.setThreads(1);
		final String flags=compressionFlags();
		final String inventoryPin=DigestSuffix.file(inventory.toString());
		final ArrayList<Loaded> loaded=new ArrayList<Loaded>();
		final HashMap<String, Boolean> ids=new HashMap<String, Boolean>();
		final HashMap<Path, String> pins=new HashMap<Path, String>();
		final ByteFile input=ByteFile.makeByteFile1(inventory.toString(), false);
		try{
			check(Arrays.equals(input.nextLine(), "genome\tnative_gff\tformula_gff".getBytes(StandardCharsets.US_ASCII)), "Wrong original retained inventory schema");
			final parse.LineParser1 parser=new parse.LineParser1('\t');
			for(byte[] row=input.nextLine(); row!=null; row=input.nextLine()){
				parser.set(row);
				final String id=parser.parseString(0);
				check(parser.terms()==3 && id.matches("g[0-9]{4}") && ids.put(id, Boolean.TRUE)==null, "Invalid repeated retained genome");
				for(int arm=0; arm<2; arm++){
					final Path file=Paths.get(parser.parseString(arm+1));
					check(file.isAbsolute() && Files.isRegularFile(file) && !Files.isSymbolicLink(file) && file.equals(file.toRealPath()), "Original retained file must be canonical");
					check(pins.put(file, DigestSuffix.file(file.toString()))==null, "Repeated original retained file");
					final ArrayList<GffLine> lines=GffLine.loadGffFile(file.toString(), null, false);
					for(GffLine line:lines){check("CDS".equals(line.type), "Full983 CDS validation received an RNA row");}
					loaded.add(new Loaded(id, arm==0 ? "native" : "formula", file, lines));
				}
			}
		}finally{check(!input.close(), "Original inventory reader failed");}
		check(ids.size()==983 && loaded.size()==1966, "Full parity requires exactly983 genomes and1966 CDS files");
		final Path serial=Files.createDirectory(root.resolve("serial")), parallel=Files.createDirectory(root.resolve("parallel"));
		final Object benchmark=instance(root, root.resolve("unused"), "true");
		final Method validate=CallGenesLeanReplayBenchmark.class.getDeclaredMethod("validateRetentionCompression");
		validate.setAccessible(true);
		validate.invoke(benchmark);
		final Method legacy=CallGenesLeanReplayBenchmark.class.getDeclaredMethod("retainCalls", Path.class, String.class, String.class, ArrayList.class, ByteBuilder.class, long[].class);
		legacy.setAccessible(true);
		final ByteBuilder costs=new ByteBuilder("genome\tpolicy\traw_payload_bytes\tstored_compressed_bytes\texport_wall_seconds\tcompression\tcommand\n");
		final long[] totals=new long[5];
		final ByteBuilder phases=new ByteBuilder("phase\tworkers\twall_seconds\tjvm_cpu_seconds\n");
		long start=System.nanoTime(), cpu=CallGenesBenchmarkSupport.processCpuNanos();
		for(Loaded item:loaded){legacy.invoke(null, serial, item.genome, item.arm, item.lines, costs, totals);}
		phases.append("legacy_serial\t1\t").appendSlow((System.nanoTime()-start)*1e-9).tab().appendSlow((CallGenesBenchmarkSupport.processCpuNanos()-cpu)*1e-9).nl();
		check(flags.equals(compressionFlags()), "Legacy serial comparison failed to restore compression flags");
		start=System.nanoTime();
		cpu=CallGenesBenchmarkSupport.processCpuNanos();
		final ArrayList<CallGenesRetainedExport.Receipt> receipts=CallGenesRetainedExport.write(loaded.size(), 8, i->{
			final Loaded item=loaded.get(i);
			return new CallGenesRetainedExport.Payload(parallel, item.genome, item.arm, CallGenesRetainedExport.cdsText(item.lines), item.lines.size());
		}, tools());
		phases.append("bounded_parallel\t8\t").appendSlow((System.nanoTime()-start)*1e-9).tab().appendSlow((CallGenesBenchmarkSupport.processCpuNanos()-cpu)*1e-9).nl();
		check(flags.equals(compressionFlags()), "Parallel comparison changed global compression flags");
		final ByteBuilder comparisons=new ByteBuilder("genome\tpolicy\tdecoded_bytes\toriginal_sha80_before\toriginal_sha80_after\n");
		final ByteBuilder parallelCosts=new ByteBuilder("genome\tpolicy\traw_payload_bytes\tstored_compressed_bytes\texport_wall_seconds\tcompression\tcommand\n");
		final ByteBuilder processes=new ByteBuilder("genome\tpolicy\tpopulation\toutput\texit\tcommand\n");
		for(int i=0; i<loaded.size(); i++){
			final Loaded item=loaded.get(i);
			final byte[] original=decoded(item.original);
			check(Arrays.equals(original, decoded(serial.resolve(item.genome+"."+item.arm+".gff.gz"))), "Legacy serialized bytes differ from accepted original");
			check(Arrays.equals(original, decoded(receipts.get(i).file)), "Parallel serialized bytes differ from accepted original");
			final String after=DigestSuffix.file(item.original.toString());
			check(after.equals(pins.get(item.original)), "Original GFF changed during parity validation");
			comparisons.append(item.genome).tab().append(item.arm).tab().append(original.length).tab().append(pins.get(item.original)).tab().append(after).nl();
			receipts.get(i).appendCost(parallelCosts);
			receipts.get(i).appendProcess(processes, "CDS");
		}
		check(inventoryPin.equals(DigestSuffix.file(inventory.toString())), "Original retained inventory changed during validation");
		validate.invoke(benchmark);
		Files.write(root.resolve("serial.costs.tsv"), costs.toBytes());
		Files.write(root.resolve("parallel.costs.tsv"), parallelCosts.toBytes());
		Files.write(root.resolve("phases.tsv"), phases.toBytes());
		Files.write(root.resolve("parity.tsv"), comparisons.toBytes());
		Files.write(root.resolve("parallel.processes.tsv"), processes.toBytes());
		Files.write(root.resolve("PASS"), "LEAN_RETAINED_EXPORT_FULL983_PASS genomes=983 files=1966 decoded_byte_identity=true source_loaded_once=true caller_runs=0 classifier_runs=0\n".getBytes(StandardCharsets.US_ASCII));
		System.out.println("LEAN_RETAINED_EXPORT_FULL983_PASS genomes=983 files=1966 decoded_byte_identity=true source_loaded_once=true caller_runs=0 classifier_runs=0");
	}
	/** One original file and the records loaded once for both sequential exporter phases. */
	private static final class Loaded{
		/** Retains one original payload and its parsed records for both exporter phases. */
		Loaded(String genome_, String arm_, Path original_, ArrayList<GffLine> lines_){
			genome=genome_;
			arm=arm_;
			original=original_;
			lines=lines_;
		}
		/** Genome identity and native/formula arm from the retained inventory. */
		final String genome, arm;
		/** Canonical original file, checked for preservation after both phases. */
		final Path original;
		/** Parsed original CDS records, reused by both sequential exporter phases. */
		final ArrayList<GffLine> lines;
	}

	/** Verify actual stored bytes, decoded payload size and observed backend in the cost record. */
	private static void checkCosts(Path output, boolean bgzf)throws Exception{
		final HashMap<String, String> values=new HashMap<String, String>();
		for(String row:Files.readAllLines(output.resolve("retained_gff_cost.tsv"), StandardCharsets.US_ASCII)){
			final String[] fields=row.split("\t", -1);
			check(fields.length==2 && values.put(fields[0], fields[1])==null, "Repeated or malformed cost key");
		}
		long raw=0, stored=0;
		for(String arm:new String[]{"native", "formula"}){
			final Path file=output.resolve("retained_gffs/g0001."+arm+".gff.gz");
			stored+=Files.size(file);
			final byte[] compressed=Files.readAllBytes(file);
			check(compressed.length>18 && (compressed[0]&255)==31 && (compressed[1]&255)==139, "Retained artifact is not gzip");
			if(bgzf){check((compressed[3]&4)!=0 && compressed[12]=='B' && compressed[13]=='C', "Large payload did not use BGZF");}
			final ByteFile reader=ByteFile.makeByteFile1(file.toString(), false);
			try{for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){raw+=line.length+1L;}}
			finally{check(!reader.close(), "Compressed retained record failed to decode");}
		}
		check(values.get("retained_files").equals("2") && Long.parseLong(values.get("raw_payload_bytes"))==raw
				&& Long.parseLong(values.get("stored_compressed_bytes"))==stored, "Export byte costs do not match actual artifacts");
		check(values.get("bgzf_files").equals(bgzf ? "2" : "0") && values.get("pigz_files").equals(bgzf ? "0" : "2"), "Wrong size-based compression backend");
		check(Double.parseDouble(values.get("export_wall_seconds"))>=0 && values.get("raw_files_retained").equals("false"), "Missing export cost or unexpected raw-file duplicate");
	}

	/** Existing42-column saved-locus format with three exact effective reference winners. */
	private static void loci(Path file)throws Exception{
		final String[] header=new String[42];
		final String[] names={"genome", "contig", "true_stop_1based", "strand", "reference_start_1based", "reference_left_1based", "reference_right_1based", "phase", "annotation_band"};
		System.arraycopy(names, 0, header, 0, names.length);
		for(int i=names.length; i<42; i++){header[i]="unused"+i;}
		final StringBuilder text=new StringBuilder(String.join("\t", header)).append('\n');
		final String[][] known={{"g0001", "c", "132", "0", "10", "10", "132", "0", "40_64"},
			{"g0001", "c", "400", "1", "597", "400", "597", "0", "65_99"},
			{"g0001", "c", "3999", "0", "3901", "3901", "3999", "0", "unknown"}};
		for(String[] row:known){
			final String[] cells=new String[42];
			Arrays.fill(cells, "NA");
			System.arraycopy(row, 0, cells, 0, row.length);
			text.append(String.join("\t", cells)).append('\n');
		}
		Files.write(file, text.toString().getBytes(StandardCharsets.US_ASCII));
	}
	/** Captures global compression settings as a literal restoration comparison, not a timing measurement. */
	private static String compressionFlags(){
		return ReadWrite.ZIPLEVEL+":"+ReadWrite.MAX_ZIP_THREADS()+":"+ReadWrite.PIGZ_BLOCKSIZE+":"+ReadWrite.PIGZ_ITERATIONS+":"
			+ReadWrite.ALLOW_ZIPLEVEL_CHANGE+":"+ReadWrite.USE_BGZF+":"+ReadWrite.USE_BGZIP+":"+ReadWrite.FORCE_BGZIP+":"
			+ReadWrite.PREFER_BGZIP+":"+ReadWrite.ALLOW_NATIVE_BGZF+":"+ReadWrite.USE_PIGZ+":"+ReadWrite.FORCE_PIGZ+":"
			+ReadWrite.RAWMODE+":"+ReadWrite.SET_ZIP_THREADS;
	}

	/** Explicit compressor path arguments selected by main; shared only within this sequential fixture. */
	private static String[] compressorOptions;
	/** Number of repeated synthetic call blocks; raised for the BGZF case and restored in main's finally. */
	private static int repetitions=1;
	/** Fails immediately with the fixture-specific explanation. */
	private static void check(boolean condition, String why){if(!condition){throw new AssertionError(why);}}
}
