package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import ddl.FamilySketchBuilder;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import prot.DdlClustCalibrationAdapter;
import prot.DigestSuffix;
import prot.ProteinSequence;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.Family;
import shared.Shared;

/**
 * Finite input, routing, crop, rejection and compression checks for the native
 * DDLClust-to-FamilyCalibration adapter. All modes retain their generated files
 * beneath a fresh evidence root. Selected-resource mode also exercises real
 * blacklist geometries and compressor integration; compression-cost is a separate
 * synthetic-byte cost probe, not a biological fitting or evaluation panel.
 * Run in an isolated JVM: evidence state is static, compression tests temporarily
 * change process-wide settings, and the optional DORI_RESOURCES override is not restored.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class DdlClustCalibrationAdapterTest{

	/** Dispatches ordinary validation, input-only generation or the explicit compression probe.
	 * @param args ROOT [RESOURCE_DIR], input-only ROOT, or compression-cost ROOT MIB WORKERS;
	 * roots must not exist and their parents must exist. The cost probe accepts 16..1024 MiB
	 * and one or four workers. Only completed ordinary/cost runs write the top-level PASS. */
	public static void main(String[] args) throws Exception{
		if(args.length==4 && args[0].equals("compression-cost")){
			compressionCost(Paths.get(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]));
			return;
		}
		if(args.length>0 && args[0].equals("input-only")){
			if(args.length!=2){throw new IllegalArgumentException("Usage: DdlClustCalibrationAdapterTest input-only FIXTURE_DIR");}
			inputOnly(Paths.get(args[1]));
			return;
		}
		if(args.length!=1 && args.length!=2){
			throw new IllegalArgumentException(
				"Usage: DdlClustCalibrationAdapterTest FIXTURE_DIR [SELECTED_RESOURCE_DIR] | input-only FIXTURE_DIR");
		}
		final Path root=Paths.get(args[0]).toAbsolutePath().normalize();
		require(!Files.exists(root), "Fixture root must be fresh");
		Files.createDirectory(root);
		evidenceRoot=root;
		evidence.clear();
		links.clear();
		passCase(root.resolve("pass"));
		defaultStructuralResidentCase(root.resolve("default_structural_resident"));
		tailCropCase(root.resolve("tailcrop"), false);
		multifileCase(root.resolve("multifile"), false);
		failCase(root.resolve("bad_id"), false, true, false);
		failCase(root.resolve("bad_numeric_map"), true, false, false);
		failCase(root.resolve("bad_stable_id"), true, true, true);
		failResidueCase(root.resolve("bad_stop"), line(100, '*'));
		failResidueCase(root.resolve("bad_illegal"), line(100, 'U'));
		if(args.length==2){
			Data.DORI_RESOURCES=Paths.get(args[1]).toAbsolutePath().normalize().toString()+"/";
			multifileCase(root.resolve("multifile_structural"), true);
			tailCropCase(root.resolve("tailcrop_structural"), true);
			final ProteinKmerGeometry[] geometries={ProteinBlacklistPolicy.GEOMETRY,
				ProteinKmerGeometry.make(5, "all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T"),
				ProteinKmerGeometry.make(6, ProteinBlacklistPolicy.DIMER120_ALPHABET)};
			for(int i=0; i<geometries.length; i++){
				structuralCase(root.resolve("selected"+i+"_short"), geometries[i], 40, 99);
				structuralCase(root.resolve("selected"+i+"_long"), geometries[i], 100, 100);
			}
			compressionBoundary(root.resolve("compression_at_boundary"), 200000);
			compressionBoundary(root.resolve("compression_above_boundary"), 200001);
			compressionRestoration(root.resolve("compression_restoration"));
			System.err.println("SELECTED_STRUCTURAL_ADAPTER_PASS geometries=3 views=2");
		}
		finish();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Compression cost       ----------------*/
	/*--------------------------------------------------------------*/

	/** Measures construction/write/close and settings checks for a repeated printable-ASCII block.
	 * Decoding and exact byte checks occur outside the measured interval. Requires
	 * threaded bgzip with level support; records whether the one-worker baseline
	 * reaches ten seconds rather than treating every run as a valid timing baseline. */
	private static void compressionCost(Path root, int mib, int workers)throws Exception{
		require(mib>=16 && mib<=1024, "Compression cost size must be16 through1024MiB");
		require(workers==1 || workers==4, "Compression cost compares1 versus4 workers only");
		root=root.toAbsolutePath().normalize();
		require(!Files.exists(root), "Compression cost root must be fresh");
		require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag,
			"Threaded bgzip with level support required");
		Files.createDirectory(root);
		evidenceRoot=root;
		evidence.clear();
		links.clear();
		final byte[] block=new byte[1<<20];
		final java.util.Random random=new java.util.Random(20261113L);
		for(int i=0; i<block.length; i++){block[i]=(byte)(32+random.nextInt(95));}
		final long bytes=(long)mib*block.length;
		final Path output=root.resolve("payload.tsv.gz");
		final CompressionSettings settings=new CompressionSettings();
		final boolean verbose=ReadWrite.verbose;
		ReadWrite.verbose=true;
		final long start=System.nanoTime();
		try{
			final ByteStreamWriter writer=DdlClustCalibrationAdapter.compressedWriter(output, bytes, workers);
			try{
				for(int i=0; i<mib; i++){writer.print(block);}
			}finally{require(!writer.poisonAndWait(), "Compression cost output failed");}
		}finally{
			ReadWrite.verbose=verbose;
			settings.check();
		}
		final long elapsedNanos=System.nanoTime()-start;
		retainIfPresent(output);
		// Check every decompressed byte outside the measured construction/write/close interval.
		try(java.io.InputStream input=new java.util.zip.GZIPInputStream(Files.newInputStream(output), 1<<16)){
			final byte[] decoded=new byte[1<<16];
			long offset=0;
			int length;
			while((length=input.read(decoded))!=-1){
				for(int i=0; i<length; i++){
					require(offset<bytes && decoded[i]==block[(int)(offset&(block.length-1))],
						"Compression cost byte mismatch");
					offset++;
				}
			}
			require(offset==bytes, "Compression cost output length mismatch");
		}
		writeBytes(root.resolve("cost.tsv"),
			("workers\tinput_bytes\tcompressed_bytes\twrite_close_elapsed_nanoseconds\tbaseline_elapsed_at_least_10_seconds\n"+
			workers+"\t"+bytes+"\t"+Files.size(output)+"\t"+elapsedNanos+"\t"+
				(workers==1 ? Boolean.toString(elapsedNanos>=10000000000L) : "not_baseline")+
				"\n").getBytes(StandardCharsets.US_ASCII));
		writeBytes(root.resolve("recipe.tsv"),
			("seed\t20261113\nblock_bytes\t1048576\nalphabet\tASCII32_through126\nrepetitions\t"+mib+
			"\ncompression\tbgzip9\nworkers\t"+workers+
			"\nscientific_panel\tfalse\n").getBytes(StandardCharsets.US_ASCII));
		finish();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Input generation       ----------------*/
	/*--------------------------------------------------------------*/

	/** Generates the shared three-family/five-assigned inputs and literal roster/ordinal tables.
	 * Does not run the adapter or call finish; success is the INPUTS_PASS diagnostic. */
	private static void inputOnly(Path root)throws Exception{
		root=root.toAbsolutePath().normalize();
		require(!Files.exists(root), "Input-only fixture root must be fresh");
		Files.createDirectory(root);
		evidenceRoot=root;
		evidence.clear();
		links.clear();
		writeInputs(root, true, true);
		writeBytes(root.resolve("roster.tsv"),
			("stable_family_id\troster_rank\tfamily_id\tassigned_members\tassigned_residues\n"+
			"20\t0\tf20\t2\t200\n"+
			"7\t1\tf7\t2\t200\n"+
			"100\t2\tf100\t1\t100\n").getBytes(StandardCharsets.US_ASCII));
		writeBytes(root.resolve("ordinals.tsv"),
			("input_index\tinput_id\tfamily_id\tstatus\troster_rank\tfamily_id_out\tbucket_at_buckets2\tlength\n"+
			"0\tg2|q0\t20\tassigned\t0\tf20\t0\t100\n"+
			"1\tg1|q1\t7\tassigned\t1\tf7\t1\t100\n"+
			"2\tg1|q2\t-1\tunassigned\t-1\t.\t-1\t100\n"+
			"3\tg3|q3\t-1\tfiltered_length\t-1\t.\t-1\t90\n"+
			"4\tg2|q4\t7\tassigned\t1\tf7\t1\t100\n"+
			"5\tg1|q5\t20\tassigned\t0\tf20\t0\t100\n"+
			"6\tg4|q6\t100\tassigned\t2\tf100\t0\t100\n").getBytes(StandardCharsets.US_ASCII));
		System.err.println("DDLCLUST_CALIBRATION_ADAPTER_INPUTS_PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------    Input and crop checks     ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal suffix oracle, including duplicate IDs and invalid skipped rows. */
	private static void tailCropCase(Path dir, boolean structural)throws Exception{
		Files.createDirectory(dir);
		final String first=repeat("ACDEFGHIKLMNPQRSTVWY", 5), second=repeat("YWVTSRQPNMLKIHGFEDCA", 5);
		final Path fasta=dir.resolve("input.faa"), members=dir.resolve("members.tsv"),
			codec=dir.resolve("families.pddl");
		final String raw=">g|dup first\nQQQ"+first+"*\n>g|dup second\nU?*"+second+"\n>u|skip\nRR"+first+"*\n";
		writeBytes(fasta, raw.getBytes(StandardCharsets.US_ASCII));
		writeBytes(members, ("input_index\tinput_id\tfamily_id\tstatus\n0\tg|dup\t42\tassigned\n"+
			"1\tg|dup\t42\tassigned\n2\tu|skip\t-1\tunassigned\n").getBytes(StandardCharsets.US_ASCII));
		final ProteinSequence[] sequences={new ProteinSequence("a", first), new ProteinSequence("b", second)};
		final Family family=new Family(42, "f42", 0, sequences[0].enc);
		family.addMember(1);
		family.refreshStats(sequences, new FamilySketchBuilder(16, 2, 1, false, 0, 0, null));
		ClusterSketchCodec.write(codec.toString(), false, new Family[]{family}, 2, 16, 1, false, 0, 0);
		String[] args={"input="+fasta, "members="+members, "codec="+codec, "out="+dir.resolve("out"),
			"buckets=2", "expectedfamilies=1", "expectedmembers=2", "minmemberaa=100", "maxmemberaa=100",
			"tailcrop=100",
			"verification=structural", "maxresidentrows=2", "maxresidentbytes=1000000"};
		if(!structural){
			args=Arrays.copyOf(args, 11);
			args[10]="verification=sha80";
		}
		runAdapter(args);
		final ByteFile bucket=ByteFile.makeByteFile1(dir.resolve("out/bucket_000.tsv.gz").toString(), false);
		try{
			require(new String(bucket.nextLine(), StandardCharsets.US_ASCII).equals(
				"family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence"), "Crop bucket header");
			require(new String(bucket.nextLine(), StandardCharsets.US_ASCII).equals("0\tf42\tg|dup\tg\t0\t"+first),
				"Literal first cropped sequence/ordinal/ID");
			require(new String(bucket.nextLine(), StandardCharsets.US_ASCII).equals("0\tf42\tg|dup\tg\t1\t"+second),
				"Duplicate ID must retain distinct second source ordinal");
			require(bucket.nextLine()==null, "Unassigned source must not become an assigned member");
		}finally{require(!bucket.close(), "Crop bucket read failed");}
		require(new String(Files.readAllBytes(dir.resolve("out/recipe.tsv")), StandardCharsets.US_ASCII)
			.contains("tailcrop\t100\nterminal_stop\tstrip_one\n"), "Cropping policy must be recorded");
		final prot.cluster.ClusteringInput loaded=prot.cluster.ClusteringInput.load(new String[]{fasta.toString()}, -1,
			false, 100);
		require(loaded.sequences.length==3 && Arrays.equals(loaded.sequences[0].enc, sequences[0].enc) &&
			Arrays.equals(loaded.sequences[1].enc, sequences[1].enc),
				"Literal adapter oracle must also match construction suffixes");
		final String[] invalid={repeat(99, 'A')+"*", first.substring(0, 99)+"U",
			first.substring(0, 99)+"*A", first+"**"};
		for(int i=0; i<invalid.length; i++){
			writeBytes(fasta, raw.replace("RR"+first+"*", invalid[i]).getBytes(StandardCharsets.US_ASCII));
			args[3]="out="+dir.resolve("invalid_skipped_"+i);
			rejectStructural(args, i==0 ? "below tailcrop" : "Unscorable retained residue");
		}
		writeBytes(fasta, raw.replace("QQQ"+first+"*", repeat(99, 'A')).getBytes(StandardCharsets.US_ASCII));
		args[3]="out="+dir.resolve("invalid_assigned");
		rejectStructural(args, "below tailcrop");
		writeBytes(fasta, raw.getBytes(StandardCharsets.US_ASCII));
		args[3]="out="+dir.resolve("disabled");
		args[9]="tailcrop=0";
		rejectStructural(args, "Assigned sequence length outside");
		args[3]="out="+dir.resolve("invalid_setting");
		args[9]="tailcrop=1";
		rejectStructural(args, "tailcrop must be0 or100");
		System.err.println("TAILCROP_ADAPTER_PASS structural="+structural+
			" literal=true duplicate_ids=true invalid_unassigned=true");
	}

	/** Preserve global ordinals across assigned, skipped and empty input-file boundaries. */
	private static void multifileCase(Path dir, boolean structural)throws Exception{
		Files.createDirectory(dir);
		writeInputs(dir, true, true);
		final String original=new String(Files.readAllBytes(dir.resolve("input.faa")), StandardCharsets.US_ASCII);
		final int firstEnd=original.indexOf(">g1|q2"), secondEnd=original.indexOf(">g2|q4");
		require(firstEnd>0 && secondEnd>firstEnd,
			"Multifile fixture must split before skipped and later assigned records");
		final Path[] parts={dir.resolve("first.faa"), dir.resolve("empty.faa"), dir.resolve("middle.faa"),
			dir.resolve("last.faa"), dir.resolve("empty_tail.faa")};
		final String[] text={original.substring(0, firstEnd), "", original.substring(firstEnd, secondEnd),
			original.substring(secondEnd), ""};
		final StringBuilder specification=new StringBuilder();
		for(int i=0; i<parts.length; i++){
			writeBytes(parts[i], text[i].getBytes(StandardCharsets.US_ASCII));
			if(i>0){specification.append(',');}
			specification.append(parts[i]);
		}
		String[] args={"input="+dir.resolve("input.faa"), "members="+dir.resolve("members.tsv"),
			"codec="+dir.resolve("families.pddl"), "out="+dir.resolve("single"), "buckets=2",
			"expectedfamilies=3", "expectedmembers=5", "verification=structural", "maxresidentrows=5",
			"maxresidentbytes=1000000"};
		if(!structural){
			args=Arrays.copyOf(args, 8);
			args[7]="verification=sha80";
		}
		runAdapter(args);
		final Path multiple=dir.resolve("multiple");
		args[0]="input="+specification;
		args[3]="out="+multiple;
		runAdapter(args);
		require(Files.isRegularFile(multiple.resolve("PASS")), "Ordered multifile conversion must complete");
		require(!Files.exists(multiple.resolve("work")),
			"Completed multifile conversion must not retain staging files");
		for(String name:new String[]{"families.tsv", "family_map.tsv", "consensus.faa", "summary.tsv"}){
			require(Arrays.equals(Files.readAllBytes(dir.resolve("single").resolve(name)), Files.readAllBytes(multiple.resolve(name))),
				"Ordered input boundaries changed final product: "+name);
		}
		for(int bucket=0; bucket<2; bucket++){
			final String name=String.format(java.util.Locale.ROOT, "bucket_%03d.tsv.gz", bucket);
			compareDecoded(dir.resolve("single").resolve(name), multiple.resolve(name));
		}
		final java.util.List<String> manifest=Files.readAllLines(multiple.resolve("inputs.sha80.tsv"),
			StandardCharsets.US_ASCII);
		require(manifest.size()==parts.length+3 && manifest.get(0).equals("input\tpath\tsha80"),
			"Every ordered source must appear in input provenance");
		for(int i=0; i<parts.length; i++){
			final String pin=structural ? "unmeasured" : DigestSuffix.file(parts[i].toString());
			require(manifest.get(i+2).equals("fasta"+i+"\t"+parts[i]+"\t"+pin),
				"Multifile provenance changed source order or identity");
		}
		final Path literal=dir.resolve("literal,comma.faa");
		writeBytes(literal, original.getBytes(StandardCharsets.US_ASCII));
		args[0]="input="+literal;
		args[3]="out="+dir.resolve("literal");
		runAdapter(args);
		final java.util.List<String> singleManifest=Files.readAllLines(dir.resolve("literal/inputs.sha80.tsv"),
			StandardCharsets.US_ASCII);
		require(singleManifest.size()==4 && singleManifest.get(2).equals("fasta\t"+literal+"\t"+
			(structural ? "unmeasured" : DigestSuffix.file(literal.toString()))),
				"Existing comma-containing filename must retain single-file provenance");
		compareDecoded(dir.resolve("single/bucket_000.tsv.gz"), dir.resolve("literal/bucket_000.tsv.gz"));
		args[0]="input="+parts[2]+","+parts[0]+","+parts[3];
		args[3]="out="+dir.resolve("wrong_order");
		rejectStructural(args, "input_id differs");
		args[0]="input="+parts[0]+","+dir.resolve("missing.faa");
		args[3]="out="+dir.resolve("missing_source");
		rejectStructural(args, "Missing regular input file");
		args[0]="input="+specification+",";
		args[3]="out="+dir.resolve("empty_name");
		rejectStructural(args, "empty filename");
		args[0]="input="+parts[0]+","+parts[2];
		args[3]="out="+dir.resolve("short_input");
		rejectStructural(args, "more rows than the input FASTA");
		writeBytes(dir.resolve("members.tsv"),
			new String(Files.readAllBytes(dir.resolve("members.tsv")), StandardCharsets.US_ASCII)
			.replace("2\tg1|q2", "0\tg1|q2").getBytes(StandardCharsets.US_ASCII));
		args[0]="input="+specification;
		args[3]="out="+dir.resolve("reset_ordinal");
		rejectStructural(args, "out of input order");
		System.err.println("ORDERED_MULTIFILE_ADAPTER_PASS structural="+structural+" sources="+parts.length);
	}

	/** Exercise each actual selected codec geometry through both original length bands.
	 * Compare decoded results with the unchanged measured path, including empty buckets. */
	private static void structuralCase(Path dir, ProteinKmerGeometry geometry, int minimum,
		int maximum)throws Exception{
		Files.createDirectory(dir);
		final String first=repeat("ACDEFGHIKLMNPQRSTVWY", 5).substring(0, minimum);
		final String second=repeat("YWVTSRQPNMLKIHGFEDCA", 5).substring(0, maximum);
		writeBytes(dir.resolve("input.faa"), (">g2|q0 description\n"+first+"\n>g1|q1\n"+second+
			"\n>g3|q2\n"+first+"\n").getBytes(StandardCharsets.US_ASCII));
		final String members="input_index\tinput_id\tfamily_id\tstatus\n0\tg2|q0\t42\tassigned\n"+
			"1\tg1|q1\t42\tassigned\n2\tg3|q2\t-1\tunassigned\n";
		writeBytes(dir.resolve("members.tsv"), members.getBytes(StandardCharsets.US_ASCII));
		final ProteinSequence[] sequences={new ProteinSequence("a", first), new ProteinSequence("b", second)};
		final Family family=new Family(42, "f42", 0, sequences[0].enc);
		family.addMember(1);
		final FamilySketchBuilder builder=new FamilySketchBuilder(512, geometry, 0, ProteinPositionMode.FIXED, 4, 2,
			null);
		family.refreshStats(sequences, builder);
		ClusterSketchCodec.write(dir.resolve("families.pddl").toString(), false, new Family[]{family}, geometry, 512, 0,
			ProteinPositionMode.FIXED, 4, 2, ProteinBlacklistPolicy.policy(geometry).sha80);
		final Path out=dir.resolve("structural"), measured=dir.resolve("measured");
		final String[] arguments={"input="+dir.resolve("input.faa"), "members="+dir.resolve("members.tsv"),
			"codec="+dir.resolve("families.pddl"), "out="+out, "buckets=2", "expectedfamilies=1", "expectedmembers=2",
			"minmemberaa="+minimum, "maxmemberaa="+maximum, "verification=structural",
			"maxresidentrows=2", "maxresidentbytes=1000000"};
		runAdapter(arguments);
		require(Files.isRegularFile(out.resolve("PASS")), "Structural conversion must finish");
		require(!Files.exists(out.resolve("work")), "Resident conversion must never create staging files");
		final String[] baseline=Arrays.copyOf(arguments, 11);
		baseline[3]="out="+measured;
		baseline[9]="verification=sha80";
		baseline[10]="tailcrop=0";
		runAdapter(baseline);
		if(minimum==40 && maximum==99){
			final Path disabled=dir.resolve("crop_disabled");
			final String[] explicit=Arrays.copyOf(arguments, 13);
			explicit[3]="out="+disabled;
			explicit[12]="tailcrop=0";
			runAdapter(explicit);
			for(String name:new String[]{"families.tsv", "family_map.tsv", "consensus.faa", "summary.tsv",
				"recipe.tsv"}){
				require(Arrays.equals(Files.readAllBytes(out.resolve(name)), Files.readAllBytes(disabled.resolve(name))),
					"Explicit tailcrop=0 changed SHORT product: "+name);
			}
			for(int bucket=0; bucket<2; bucket++){
				final String name=String.format(java.util.Locale.ROOT, "bucket_%03d.tsv.gz", bucket);
				compareDecoded(out.resolve(name), disabled.resolve(name));
			}
		}
		for(String name:new String[]{"families.tsv", "family_map.tsv", "consensus.faa", "summary.tsv"}){
			require(Arrays.equals(Files.readAllBytes(out.resolve(name)), Files.readAllBytes(measured.resolve(name))),
				"Structural conversion changed decoded product: "+name);
		}
		for(int bucket=0; bucket<2; bucket++){
			final String name=String.format(java.util.Locale.ROOT, "bucket_%03d.tsv.gz", bucket);
			final long[] counts=countBucket(out.resolve(name), bucket);
			require(counts[0]==(bucket==0 ? 2 : 0) && counts[1]==(bucket==0 ? minimum+maximum : 0),
				"Selected structural bucket census differs");
			compareDecoded(out.resolve(name), measured.resolve(name));
		}
		for(String row:Files.readAllLines(out.resolve("inputs.sha80.tsv"), StandardCharsets.US_ASCII)){
			require(row.startsWith("input\t") || row.endsWith("\tunmeasured"),
				"Structural input identity must be honest");
		}
		final java.util.List<String> chunks=Files.readAllLines(out.resolve("chunks.tsv"), StandardCharsets.US_ASCII);
		require(chunks.size()==3 && chunks.get(1).split("\t", -1)[3].equals("unmeasured") &&
			chunks.get(2).split("\t", -1)[3].equals("unmeasured"), "Every final chunk must retain unmeasured identity");
		require(Files.readAllLines(out.resolve("models.tsv"), StandardCharsets.US_ASCII).get(1).split("\t", -1)[6].equals("unmeasured"),
			"Consensus model identity must be unmeasured");
		require(new String(Files.readAllBytes(out.resolve("recipe.tsv")), StandardCharsets.US_ASCII).contains("verification\tstructural\n"),
			"Structural provenance must be explicit");
		final String[] invalid=arguments.clone();
		invalid[3]="out="+dir.resolve("bad_mode");
		invalid[9]="verification=unknown";
		rejectStructural(invalid, "verification must be");
		invalid[3]="out="+dir.resolve("wrong_count");
		invalid[9]="verification=structural";
		invalid[6]="expectedmembers=3";
		rejectStructural(invalid, "Expected member count differs");
		invalid[6]=arguments[6];
		invalid[3]="out="+dir.resolve("keepwork_resident");
		invalid[6]="keepwork=t";
		rejectStructural(invalid, "Resident storage has no work files");
		invalid[6]=arguments[6];
		final String[] unpaired=Arrays.copyOf(arguments, 11);
		unpaired[3]="out="+dir.resolve("unpaired_resident_limit");
		rejectStructural(unpaired, "requires maxresidentrows and maxresidentbytes together");
		final String[] sha80Unpaired=unpaired.clone();
		sha80Unpaired[3]="out="+dir.resolve("unpaired_sha80_resident_limit");
		sha80Unpaired[9]="verification=sha80";
		rejectStructural(sha80Unpaired, "requires maxresidentrows and maxresidentbytes together");
		invalid[3]="out="+dir.resolve("row_budget");
		invalid[10]="maxresidentrows=1";
		rejectStructural(invalid, "exceeds maxresidentrows");
		invalid[10]=arguments[10];
		invalid[3]="out="+dir.resolve("array_byte_budget");
		invalid[11]="maxresidentbytes=1";
		rejectStructural(invalid, "exceeds maxresidentbytes");
		invalid[3]="out="+dir.resolve("row_byte_budget");
		invalid[11]="maxresidentbytes=1400";
		rejectStructural(invalid, "exceeds maxresidentbytes");
		invalid[11]=arguments[11];
		invalid[3]="out="+dir.resolve("missing_limit");
		invalid[10]="keepwork=f";
		rejectStructural(invalid, "requires maxresidentrows and maxresidentbytes together");
		invalid[10]=arguments[10];
		final String[] retain=Arrays.copyOf(arguments, 13);
		retain[3]="out="+dir.resolve("keepwork");
		retain[12]="keepwork=t";
		rejectStructural(retain, "has no work files");
		writeBytes(dir.resolve("members.tsv"),
			members.replace("0\tg2|q0", "0\tg2|wrong").getBytes(StandardCharsets.US_ASCII));
		invalid[3]="out="+dir.resolve("wrong_id");
		rejectStructural(invalid, "input_id differs");
		writeBytes(dir.resolve("members.tsv"), members.getBytes(StandardCharsets.US_ASCII));
		symbolicLink(dir.resolve("input_link.faa"), dir.resolve("input.faa"));
		invalid[0]="input="+dir.resolve("input_link.faa");
		invalid[3]="out="+dir.resolve("symlink");
		rejectStructural(invalid, "non-symlink");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Compression checks      ----------------*/
	/*--------------------------------------------------------------*/

	/** Verify actual gzip format and decoded size on both sides of the exact raw-byte policy boundary. */
	private static void compressionBoundary(Path dir, int targetBytes)throws Exception{
		Files.createDirectory(dir);
		final int count=1000;
		final String sequence=repeat("ACDEFGHIKLMNPQRSTVWY", 5);
		final String[] ids=new String[count];
		long bytes="family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence\n".length();
		for(int i=0; i<count; i++){
			ids[i]=(i%2==0 ? "g2" : "g1")+"|q"+i;
			bytes+=("0\tf0\t"+ids[i]+"\t"+ids[i].substring(0, 2)+"\t"+i+"\t"+sequence+"\n").length();
		}
		require(bytes<targetBytes, "Compression fixture base must fit below boundary");
		ids[count-1]+=repeat(Math.toIntExact(targetBytes-bytes), 'x');
		final StringBuilder fasta=new StringBuilder(),
			members=new StringBuilder("input_index\tinput_id\tfamily_id\tstatus\n");
		final ProteinSequence[] sequences=new ProteinSequence[count];
		final Family family=new Family(0, "f0", 0, new ProteinSequence("seed", sequence).enc);
		for(int i=0; i<count; i++){
			fasta.append('>').append(ids[i]).append('\n').append(sequence).append('\n');
			members.append(i).append('\t').append(ids[i]).append("\t0\tassigned\n");
			sequences[i]=new ProteinSequence(ids[i], sequence);
			if(i>0){family.addMember(i);}
		}
		family.refreshStats(sequences, new FamilySketchBuilder(16, 2, 1, false, 0, 0, null));
		writeBytes(dir.resolve("input.faa"), fasta.toString().getBytes(StandardCharsets.US_ASCII));
		writeBytes(dir.resolve("members.tsv"), members.toString().getBytes(StandardCharsets.US_ASCII));
		ClusterSketchCodec.write(dir.resolve("families.pddl").toString(), false, new Family[]{family}, 2, 16, 1, false,
			0, 0);
		final String[] args={"input="+dir.resolve("input.faa"), "members="+dir.resolve("members.tsv"),
			"codec="+dir.resolve("families.pddl"), "out="+dir.resolve("out"), "buckets=2", "expectedfamilies=1",
			"expectedmembers="+count, "verification=structural", "maxresidentrows="+count, "maxresidentbytes=1000000"};
		runCompressionAdapter(args);
		final Path bucket=dir.resolve("out/bucket_000.tsv.gz");
		final byte[] compressed=Files.readAllBytes(bucket);
		final boolean bgzf=(compressed[3]&4)!=0 && compressed[12]=='B' && compressed[13]=='C';
		require(bgzf==(targetBytes>200000), "Wrong compression format at raw-byte boundary");
		final ByteFile in=ByteFile.makeByteFile1(bucket.toString(), false);
		long decoded=0, lines=0;
		try{
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				decoded+=line.length+1;
				lines++;
			}
		}finally{require(!in.close(), "Compressed boundary fixture read failed");}
		require(decoded==targetBytes && lines==count+1, "Compression boundary fixture lost decoded bytes or rows");
		final long[] census=countBucket(bucket, 0);
		require(census[0]==count && census[1]==100L*count, "Compression boundary fixture changed members or residues");
		for(int workers:new int[]{1, 4}){
			final Path changed=dir.resolve("workers"+workers);
			final String[] configured=Arrays.copyOf(args, args.length+1);
			configured[3]="out="+changed;
			configured[args.length]="compressionthreads="+workers;
			runCompressionAdapter(configured);
			for(String name:new String[]{"families.tsv", "family_map.tsv", "consensus.faa", "summary.tsv"}){
				require(Arrays.equals(Files.readAllBytes(dir.resolve("out").resolve(name)), Files.readAllBytes(changed.resolve(name))),
					"Compression workers changed collection metadata: "+name);
			}
			for(int index=0; index<2; index++){
				final String name=String.format(java.util.Locale.ROOT, "bucket_%03d.tsv.gz", index);
				compareDecoded(dir.resolve("out").resolve(name), changed.resolve(name));
				if(workers==1){compareGzipBytesExceptMtime(dir.resolve("out").resolve(name), changed.resolve(name));}
			}
			final byte[] header=Files.readAllBytes(changed.resolve("bucket_000.tsv.gz"));
			final boolean changedBgzf=(header[3]&4)!=0 && header[12]=='B' && header[13]=='C';
			require(changedBgzf==bgzf, "Compression workers changed codec at boundary");
			final String recipe=new String(Files.readAllBytes(changed.resolve("recipe.tsv")),
				StandardCharsets.US_ASCII);
			final String baselineRecipe=new String(Files.readAllBytes(dir.resolve("out/recipe.tsv")),
				StandardCharsets.US_ASCII);
			require(recipe.equals(baselineRecipe.replace("compression_threads\t1\n", "compression_threads\t"+workers+"\n")),
				"Only the recorded compressor worker count may change");
			require(Shared.threads()==1 && Shared.numBuffers()==4 && ReadWrite.MAX_ZIP_THREADS()==1,
				"Compressor workers leaked into collection settings");
		}
		for(int invalid:new int[]{-1, 0, 97}){
			final String[] rejected=Arrays.copyOf(args, args.length+1);
			rejected[3]="out="+dir.resolve("invalid_workers_"+invalid);
			rejected[args.length]="compressionthreads="+invalid;
			rejectStructural(rejected, "compressionthreads must be1 through96");
			require(!Files.exists(Paths.get(rejected[3].substring(4))), "Invalid compression count created outputs");
		}
		final String[] sha80Resident=Arrays.copyOf(args, args.length+1);
		sha80Resident[3]="out="+dir.resolve("sha80_resident_workers");
		sha80Resident[7]="verification=sha80";
		sha80Resident[args.length]="compressionthreads=4";
		runCompressionAdapter(sha80Resident);
		compareDecoded(dir.resolve("out/bucket_000.tsv.gz"), dir.resolve("sha80_resident_workers/bucket_000.tsv.gz"));
		final String sha80Recipe=new String(Files.readAllBytes(
			dir.resolve("sha80_resident_workers/recipe.tsv")), StandardCharsets.US_ASCII);
		require(sha80Recipe.contains("verification\tsha80\n") &&
			sha80Recipe.contains("storage\tbounded_resident\n"),
			"SHA80 resident conversion must record independent verification and storage policies");
		final String[] defaultResidentWorkers=Arrays.copyOf(args, 9);
		defaultResidentWorkers[3]="out="+dir.resolve("default_resident_workers");
		defaultResidentWorkers[defaultResidentWorkers.length-1]="compressionthreads=4";
		runCompressionAdapter(defaultResidentWorkers);
		compareDecoded(dir.resolve("out/bucket_000.tsv.gz"), dir.resolve("default_resident_workers/bucket_000.tsv.gz"));
		final String defaultRecipe=new String(Files.readAllBytes(
			dir.resolve("default_resident_workers/recipe.tsv")), StandardCharsets.US_ASCII);
		require(defaultRecipe.contains("storage\tbounded_resident\n") &&
			defaultRecipe.contains("compression_threads\t4\n"),
			"Default resident conversion must accept explicit compression workers");
	}

	/** Retain actual child command diagnostics for independent codec/worker review. */
	private static void runCompressionAdapter(String[] args)throws Exception{
		final boolean verbose=ReadWrite.verbose;
		ReadWrite.verbose=true;
		try{
			runAdapter(args);
		}finally{
			ReadWrite.verbose=verbose;
		}
	}

	/** Exercise exact global restoration, including Shared.setThreads' buffer side effect and failed stream creation. */
	private static void compressionRestoration(Path dir)throws Exception{
		Files.createDirectory(dir);
		final CompressionSettings original=new CompressionSettings();
		try{
			Shared.setThreads(3);
			Shared.setBuffers(17);
			ReadWrite.setZipThreads(7);
			ReadWrite.setZipThreadMult(0.5f);
			ReadWrite.SET_ZIP_THREADS=false;
			ReadWrite.SET_ZIP_THREAD_MULT=false;
			ReadWrite.ZIPLEVEL=2;
			ReadWrite.USE_BGZF=ReadWrite.USE_BGZIP=ReadWrite.USE_PIGZ=false;
			ReadWrite.FORCE_BGZIP=ReadWrite.FORCE_PIGZ=ReadWrite.PREFER_BGZIP=false;
			ReadWrite.ALLOW_NATIVE_BGZF=true;
			final CompressionSettings sentinel=new CompressionSettings();
			for(int raw:new int[]{200000, 200001}){
				final Path path=dir.resolve("restored_"+raw+".gz");
				final ByteStreamWriter writer=DdlClustCalibrationAdapter.compressedWriter(path, raw, 4);
				sentinel.check();
				writer.println("restoration");
				require(!writer.poisonAndWait(), "Restoration fixture compression failed");
				retainIfPresent(path);
				sentinel.check();
				final Path parent=dir.resolve("not_a_directory_"+raw);
				writeBytes(parent, new byte[]{1});
				boolean failed=false;
				try{
					DdlClustCalibrationAdapter.compressedWriter(parent.resolve("failure.gz"), raw, 4);
				}catch(AssertionError expected){failed=true;}catch(RuntimeException expected){failed=true;}
				require(failed, "Stream construction through a regular file must fail");
				sentinel.check();
			}
		}finally{original.restore();}
	}

	/** Snapshot the constructor's complete mutable compression configuration, including explicit-setting markers. */
	private static final class CompressionSettings{
		/** Caller thread/buffer and compressor limits, captured in declaration order. */
		final int threads=Shared.threads(), buffers=Shared.numBuffers(), zip=ReadWrite.MAX_ZIP_THREADS(),
			level=ReadWrite.ZIPLEVEL;
		/** Compressor thread multiplier. */
		final float multiplier=ReadWrite.ZIP_THREAD_MULT();
		/** Whether the caller explicitly selected each compressor concurrency setting. */
		final boolean setZip=ReadWrite.SET_ZIP_THREADS, setMultiplier=ReadWrite.SET_ZIP_THREAD_MULT;
		/** Compression implementation, forced/preferred tool and native fallback switches. */
		final boolean bgzf=ReadWrite.USE_BGZF, bgzip=ReadWrite.USE_BGZIP, pigz=ReadWrite.USE_PIGZ,
			forceBgzip=ReadWrite.FORCE_BGZIP, forcePigz=ReadWrite.FORCE_PIGZ, prefer=ReadWrite.PREFER_BGZIP,
			nativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
		/** Restores threads before buffers and setters before explicit-setting markers. */
		void restore(){
			Shared.setThreads(threads);
			Shared.setBuffers(buffers);
			ReadWrite.setZipThreads(zip);
			ReadWrite.setZipThreadMult(multiplier);
			ReadWrite.SET_ZIP_THREADS=setZip;
			ReadWrite.SET_ZIP_THREAD_MULT=setMultiplier;
			ReadWrite.ZIPLEVEL=level;
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.USE_BGZIP=bgzip;
			ReadWrite.USE_PIGZ=pigz;
			ReadWrite.FORCE_BGZIP=forceBgzip;
			ReadWrite.FORCE_PIGZ=forcePigz;
			ReadWrite.PREFER_BGZIP=prefer;
			ReadWrite.ALLOW_NATIVE_BGZF=nativeBgzf;
		}
		/** Requires every captured setting to match; this does not restore changed state. */
		void check(){
			require(Shared.threads()==threads && Shared.numBuffers()==buffers && ReadWrite.MAX_ZIP_THREADS()==zip &&
				ReadWrite.ZIP_THREAD_MULT()==multiplier && ReadWrite.SET_ZIP_THREADS==setZip &&
				ReadWrite.SET_ZIP_THREAD_MULT==setMultiplier &&
				ReadWrite.ZIPLEVEL==level && ReadWrite.USE_BGZF==bgzf && ReadWrite.USE_BGZIP==bgzip &&
				ReadWrite.USE_PIGZ==pigz &&
				ReadWrite.FORCE_BGZIP==forceBgzip && ReadWrite.FORCE_PIGZ==forcePigz &&
				ReadWrite.PREFER_BGZIP==prefer &&
				ReadWrite.ALLOW_NATIVE_BGZF==nativeBgzf, "Compression stream construction changed caller globals");
		}
	}

	/** Compare all decoded member bytes so compression cannot hide a routing/order regression. */
	private static void compareDecoded(Path left, Path right){
		final ByteFile a=ByteFile.makeByteFile1(left.toString(), false),
			b=ByteFile.makeByteFile1(right.toString(), false);
		try{
			for(byte[] row=a.nextLine(); row!=null; row=a.nextLine()){
				require(Arrays.equals(row, b.nextLine()), "Structural member rows differ from measured conversion");
			}
			require(b.nextLine()==null, "Measured conversion has extra member rows");
		}finally{
			final boolean failedA=a.close(), failedB=b.close();
			require(!failedA && !failedB, "Decoded bucket comparison failed");
		}
	}

	/** External gzip tools may stamp MTIME; all header flags, payload and trailer bytes must stay exact. */
	private static void compareGzipBytesExceptMtime(Path left, Path right) throws Exception{
		final byte[] a=Files.readAllBytes(left), b=Files.readAllBytes(right);
		if(Arrays.equals(a, b)){return;}
		require(a.length==b.length && a.length>=10, "Explicit default changed compressed length");
		require((a[0]&0xff)==0x1f && (a[1]&0xff)==0x8b && a[2]==8 &&
			(b[0]&0xff)==0x1f && (b[1]&0xff)==0x8b && b[2]==8,
			"Explicit default changed compressed codec");
		require(a[3]==0 && b[3]==0,
			"Only plain pigz output may differ in gzip MTIME; BGZF or optional gzip headers must stay exact");
		for(int i=0; i<a.length; i++){
			if(i>=4 && i<8){continue;}
			require(a[i]==b[i], "Explicit default changed compressed bytes outside gzip MTIME at offset "+i);
		}
	}

	/** Require the intended validation failure and absence of a successful completion marker. */
	private static void rejectStructural(String[] args, String message)throws Exception{
		boolean failed=false;
		try{runAdapter(args);}catch(IllegalArgumentException expected){
			require(expected.getMessage().contains(message), "Wrong structural failure: "+expected.getMessage());
			failed=true;
		}
		require(failed, "Invalid structural conversion passed");
		require(!Files.exists(Paths.get(args[3].substring(4)).resolve("PASS")),
			"Failed structural conversion wrote PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------    Literal adapter cases     ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks measured conversion against the literal three-family, five-member, two-bucket oracle. */
	private static void passCase(Path dir) throws Exception{
		Files.createDirectory(dir);
		writeInputs(dir, true, true);
		runAdapter(new String[]{
			"input="+dir.resolve("input.faa"),
			"members="+dir.resolve("members.tsv"),
			"codec="+dir.resolve("families.pddl"),
			"out="+dir.resolve("out"),
			"buckets=2",
			"expectedfamilies=3",
			"expectedmembers=5",
			"verification=sha80"
		});
		require(Files.isRegularFile(dir.resolve("out").resolve("PASS")), "Adapter did not write PASS");
		require(!Files.exists(dir.resolve("out").resolve("work")), "Default run retained work files");
		long[] bucket0=countBucket(dir.resolve("out").resolve("bucket_000.tsv.gz"), 0);
		long[] bucket1=countBucket(dir.resolve("out").resolve("bucket_001.tsv.gz"), 1);
		require(bucket0[0]==3 && bucket0[1]==300 && bucket1[0]==2 && bucket1[1]==200, "Wrong bucket counts");
		checkFamilyMap(dir.resolve("out").resolve("family_map.tsv"));
		checkModels(dir.resolve("out"));
		checkBucketSequence(dir.resolve("out").resolve("bucket_000.tsv.gz"), "g4|q6", repeat("XBZJ", 25));
	}

	/** Checks default resident conversion, its literal counts and explicit unmeasured provenance. */
	private static void defaultStructuralResidentCase(Path dir) throws Exception{
		Files.createDirectory(dir);
		writeInputs(dir, true, true);
		runAdapter(new String[]{
			"input="+dir.resolve("input.faa"),
			"members="+dir.resolve("members.tsv"),
			"codec="+dir.resolve("families.pddl"),
			"out="+dir.resolve("out"),
			"buckets=2"
		});
		final Path out=dir.resolve("out");
		require(Files.isRegularFile(out.resolve("PASS")), "Default structural resident run failed");
		require(!Files.exists(out.resolve("work")),
			"Default structural resident run must avoid stage rows");
		final long[] bucket0=countBucket(out.resolve("bucket_000.tsv.gz"), 0);
		final long[] bucket1=countBucket(out.resolve("bucket_001.tsv.gz"), 1);
		require(bucket0[0]==3 && bucket0[1]==300 && bucket1[0]==2 && bucket1[1]==200,
			"Default structural resident bucket counts changed");
		checkFamilyMap(out.resolve("family_map.tsv"));
		final String recipe=new String(Files.readAllBytes(out.resolve("recipe.tsv")), StandardCharsets.US_ASCII);
		require(recipe.contains("verification\tstructural\n") &&
			recipe.contains("storage\tbounded_resident\n") &&
			recipe.contains("maxresidentrows\t5\n") &&
			recipe.contains("maxresidentbytes\t"),
			"Default structural run must record bounded resident provenance");
		for(String row:Files.readAllLines(out.resolve("inputs.sha80.tsv"), StandardCharsets.US_ASCII)){
			require(row.startsWith("input\t") || row.endsWith("\tunmeasured"),
				"Default structural input identity must be unmeasured");
		}
	}

	/** Builds a mismatched ID/count/stable-ID case and requires rejection without PASS. */
	private static void failCase(Path dir, boolean matchingIds, boolean goodCounts,
		boolean missingStableId) throws Exception{
		Files.createDirectory(dir);
		writeInputs(dir, matchingIds, goodCounts, missingStableId);
		boolean failed=false;
		try{
			runAdapter(new String[]{
				"input="+dir.resolve("input.faa"),
				"members="+dir.resolve("members.tsv"),
				"codec="+dir.resolve("families.pddl"),
				"out="+dir.resolve("out"),
				"buckets=2"
			});
		}catch(IllegalArgumentException expected){failed=true;}
		require(failed, "Bad fixture should have failed: "+dir);
		require(!Files.exists(dir.resolve("out").resolve("PASS")), "Failed fixture wrote PASS");
	}

	/** Writes the shared inputs with all assigned stable IDs present; flags select ID/count corruption. */
	private static void writeInputs(Path dir, boolean matchingIds, boolean goodCounts) throws Exception{
		writeInputs(dir, matchingIds, goodCounts, false);
	}

	/** Writes the shared inputs with the literal X/B/Z/J ambiguity sequence in the last assigned row. */
	private static void writeInputs(Path dir, boolean matchingIds, boolean goodCounts,
		boolean missingStableId) throws Exception{
		writeInputs(dir, matchingIds, goodCounts, missingStableId, line(repeat("XBZJ", 25)));
	}

	/** Writes seven source records and three codec families beneath an existing directory.
	 * Five rows are assigned, one unassigned and one length-filtered. Flags deliberately
	 * corrupt the first ID, last family count or last assigned stable ID; q6 supplies
	 * the final FASTA sequence text including its terminator. Existing inputs are replaced. */
	private static void writeInputs(Path dir, boolean matchingIds, boolean goodCounts, boolean missingStableId,
		String q6) throws Exception{
		writeBytes(dir.resolve("input.faa"), (">g2|q0 desc\n"+line(100, 'A')+
			">g1|q1\n"+line(100, 'C')+
			">g1|q2\n"+line(100, 'D')+
			">g3|q3\n"+line(90, 'E')+
			">g2|q4\n"+line(100, 'F')+
			">g1|q5\n"+line(100, 'G')+
			">g4|q6\n"+q6).getBytes(StandardCharsets.US_ASCII));
		writeBytes(dir.resolve("members.tsv"), ("input_index\tinput_id\tfamily_id\tstatus\n"+
			"0\t"+(matchingIds ? "g2|q0" : "g2|WRONG")+"\t20\tassigned\n"+
			"1\tg1|q1\t7\tassigned\n"+
			"2\tg1|q2\t-1\tunassigned\n"+
			"3\tg3|q3\t-1\tfiltered_length\n"+
			"4\tg2|q4\t7\tassigned\n"+
			"5\tg1|q5\t20\tassigned\n"+
			"6\tg4|q6\t"+(missingStableId ? "999" : "100")+"\tassigned\n").getBytes(StandardCharsets.US_ASCII));
		final ProteinSequence[] sequences={
			new ProteinSequence("f20a", repeat(100, 'A')), new ProteinSequence("f20b", repeat(100, 'C')),
			new ProteinSequence("f7a", repeat(100, 'D')), new ProteinSequence("f7b", repeat(100, 'F')),
			new ProteinSequence("f100a", repeat(100, 'G')), new ProteinSequence("wrong_extra", repeat(100, 'H'))};
		final FamilySketchBuilder sketcher=new FamilySketchBuilder(16, 2, 1, false, 0, 0, null);
		final Family f20=new Family(20, "f20", 0, sequences[0].enc);
		f20.addMember(1);
		f20.refreshStats(sequences, sketcher);
		final Family f7=new Family(7, "f7", 2, sequences[2].enc);
		f7.addMember(3);
		f7.refreshStats(sequences, sketcher);
		final Family f100=new Family(100, "f100", 4, sequences[4].enc);
		if(!goodCounts){f100.addMember(5);}
		f100.refreshStats(sequences, sketcher);
		ClusterSketchCodec.write(dir.resolve("families.pddl").toString(), false, new Family[]{f20, f7, f100}, 2, 16, 1,
			false, 0, 0);
	}

	/** Requires rejection and no PASS when the last assigned sequence contains an unscorable residue. */
	private static void failResidueCase(Path dir, String q6) throws Exception{
		Files.createDirectory(dir);
		writeInputs(dir, true, true, false, q6);
		boolean failed=false;
		try{
			runAdapter(new String[]{
				"input="+dir.resolve("input.faa"),
				"members="+dir.resolve("members.tsv"),
				"codec="+dir.resolve("families.pddl"),
				"out="+dir.resolve("out"),
				"buckets=2"
			});
		}catch(IllegalArgumentException expected){failed=true;}
		require(failed, "Unscorable assigned residue should have failed: "+dir);
		require(!Files.exists(dir.resolve("out").resolve("PASS")), "Failed residue fixture wrote PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Output oracles        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads and closes the map, requiring the exact stable-ID order 20, 7, 100 and literal counts. */
	private static void checkFamilyMap(Path path){
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(new String(input.nextLine(), StandardCharsets.US_ASCII).equals(
				"stable_family_id\troster_rank\tfamily_id\tassigned\tresidues"),
				"Wrong family_map header");
			checkMapRow(input.nextLine(), lp, 20, 0, "f20", 2, 200);
			checkMapRow(input.nextLine(), lp, 7, 1, "f7", 2, 200);
			checkMapRow(input.nextLine(), lp, 100, 2, "f100", 1, 100);
			require(input.nextLine()==null, "Unexpected family_map row");
		}finally{require(!input.close(), "family_map read failed");}
	}

	/** Checks all five map fields against literal expectations, reusing the caller-owned parser. */
	private static void checkMapRow(byte[] line, LineParser1 lp, int stableId, int rank, String family, int assigned,
		int residues){
		lp.set(line);
		require(lp.terms()==5 && lp.parseInt(0)==stableId && lp.parseInt(1)==rank &&
			lp.termEquals(family, 2) && lp.parseInt(3)==assigned && lp.parseInt(4)==residues,
			"Wrong stable-ID to display-family mapping");
	}

	/** Checks the ordered consensus headers and exact model rows, including their shared FASTA binding. */
	private static void checkModels(Path out) throws Exception{
		final String pin=DigestSuffix.file(out.resolve("consensus.faa").toString());
		final java.util.List<String> fasta=Files.readAllLines(out.resolve("consensus.faa"), StandardCharsets.US_ASCII);
		require(fasta.size()==6 && fasta.get(0).equals(">f20") && fasta.get(2).equals(">f7") &&
			fasta.get(4).equals(">f100"),
			"Wrong shared consensus FASTA");
		final java.util.List<String> models=Files.readAllLines(out.resolve("models.tsv"), StandardCharsets.US_ASCII);
		require(models.size()==4 &&
			models.get(0).equals("rank\tfamily_id\tdirectory\tmembers\tresidues\tbundle_sha80\tconsensus_sha80\tprovenance_sha80\tretained_sha80"),
			"Wrong model roster header");
		require(models.get(1).equals("0\tf20\t"+out+"\t2\t200\t.\t"+pin+"\t.\t.") &&
			models.get(2).equals("1\tf7\t"+out+"\t2\t200\t.\t"+pin+"\t.\t.") &&
			models.get(3).equals("2\tf100\t"+out+"\t1\t100\t.\t"+pin+"\t.\t."),
			"Wrong shared consensus model roster");
	}

	/** Checks six-column rows, rank modulo two routing and increasing genome/ordinal order.
	 * Closes its reader and returns a new array containing member-row and amino-acid totals. */
	private static long[] countBucket(Path path, int bucket){
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		String lastGenome="";
		long rows=0, residues=0, lastOrdinal=-1;
		try{
			require(new String(input.nextLine(), StandardCharsets.US_ASCII).equals(
				"family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence"), "Wrong bucket header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==6 && lp.parseInt(0)%2==bucket, "Malformed bucket row");
				final String genome=lp.parseString(3);
				if(!genome.equals(lastGenome)){
					require(genome.compareTo(lastGenome)>0, "Bucket genome order failed");
					lastGenome=genome;
					lastOrdinal=-1;
				}
				require(lp.parseLong(4)>lastOrdinal, "Bucket ordinal order failed");
				lastOrdinal=lp.parseLong(4);
				rows++;
				residues+=lp.length(5);
			}
		}finally{require(!input.close(), "Bucket read failed");}
		return new long[]{rows, residues};
	}

	/** Requires exactly one row for id with the literal ambiguity sequence; closes its reader. */
	private static void checkBucketSequence(Path path, String id, String sequence){
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		boolean found=false;
		try{
			require(input.nextLine()!=null, "Missing bucket header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.termEquals(id, 2)){
					require(!found, "Duplicate row for ambiguous residue fixture");
					require(lp.termEquals(sequence, 5), "Assigned ambiguous residues were not preserved");
					found=true;
				}
			}
		}finally{require(!input.close(), "Bucket read failed");}
		require(found, "Assigned ambiguous residue fixture was not routed");
	}

	/*--------------------------------------------------------------*/
	/*----------------     Evidence collection      ----------------*/
	/*--------------------------------------------------------------*/

	/** Register declared fixture outputs, including files retained after expected failures. */
	private static void runAdapter(String[] args)throws Exception{
		Path out=null;
		int buckets=256;
		for(String arg:args){
			if(arg.startsWith("out=")){out=Paths.get(arg.substring(4));}
			if(arg.startsWith("buckets=")){buckets=Integer.parseInt(arg.substring(8));}
			if(arg.startsWith("codec=")){retainIfPresent(Paths.get(arg.substring(6)));}
		}
		try{DdlClustCalibrationAdapter.runForTest(args);}finally{
			if(out!=null){
				for(String name:new String[]{"families.tsv", "family_map.tsv", "consensus.faa", "models.tsv",
					"chunks.tsv", "buckets.tsv", "inputs.sha80.tsv", "summary.tsv", "recipe.tsv", "PASS"}){
					retainIfPresent(out.resolve(name));
				}
				for(int i=0; i<buckets; i++){
					retainIfPresent(out.resolve(String.format(java.util.Locale.ROOT, "bucket_%03d.tsv.gz", i)));
					retainIfPresent(out.resolve("work").resolve(String.format(java.util.Locale.ROOT, "stage_%03d.tsv",
						i)));
				}
			}
		}
	}

	/** Returns a normalized path relative to evidenceRoot, rejecting lexical paths outside that root. */
	private static String relative(Path path){
		final Path absolute=path.toAbsolutePath().normalize();
		require(absolute.startsWith(evidenceRoot), "Fixture output outside evidence root");
		return evidenceRoot.relativize(absolute).toString();
	}

	/** Registers an existing regular, non-symlink file in the sorted relative evidence roster. */
	private static void retainIfPresent(Path path){
		if(Files.isRegularFile(path) && !Files.isSymbolicLink(path)){evidence.add(relative(path));}
	}

	/** Writes fixture bytes with the supplied open options, registers the file and returns its path.
	 * With no options, the standard Files.write create/truncate behavior applies. */
	private static Path writeBytes(Path path, byte[] bytes,
		java.nio.file.OpenOption... options)throws java.io.IOException{
		Files.write(path, bytes, options);
		retainIfPresent(path);
		return path;
	}

	/** Creates a fixture symlink and records its relative link path and supplied target spelling. */
	private static Path symbolicLink(Path link, Path target)throws java.io.IOException{
		Files.createSymbolicLink(link, target);
		links.put(relative(link), target.toString());
		return link;
	}

	/** Writes PASS, the sorted symlink table and a new evidence.paths including itself.
	 * Reports success to stderr after those writes; failures leave partial evidence in place. */
	private static void finish()throws Exception{
		final String marker="DDLCLUST_CALIBRATION_ADAPTER_TEST_PASS";
		writeBytes(evidenceRoot.resolve("PASS"), (marker+"\n").getBytes(StandardCharsets.US_ASCII));
		final StringBuilder text=new StringBuilder();
		for(java.util.Map.Entry<String, String> link:links.entrySet()){
			text.append(link.getKey()).append('\t').append(link.getValue()).append('\n');
		}
		writeBytes(evidenceRoot.resolve("evidence.links.tsv"), text.toString().getBytes(StandardCharsets.UTF_8));
		evidence.add("evidence.paths");
		Files.write(evidenceRoot.resolve("evidence.paths"), evidence, StandardCharsets.UTF_8,
			java.nio.file.StandardOpenOption.CREATE_NEW);
		System.err.println(marker);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Evidence state        ----------------*/
	/*--------------------------------------------------------------*/

	/** Absolute normalized root for the active standalone invocation. */
	private static Path evidenceRoot;
	/** Sorted unique relative file paths; cleared when a mode starts. */
	private static final java.util.TreeSet<String> evidence=new java.util.TreeSet<String>();
	/** Sorted relative symlink paths and target spellings; kept separately from regular files. */
	private static final java.util.TreeMap<String, String> links=new java.util.TreeMap<String, String>();

	/*--------------------------------------------------------------*/
	/*----------------       Sequence helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a repeated-character sequence followed by one newline. */
	private static String line(int length, char c){return repeat(length, c)+"\n";}
	/** Appends one newline to the supplied literal sequence. */
	private static String line(String s){return s+"\n";}

	/** Builds a literal sequence containing length copies of c. */
	private static String repeat(int length, char c){
		final StringBuilder out=new StringBuilder(length+1);
		for(int i=0; i<length; i++){out.append(c);}
		return out.toString();
	}

	/** Concatenates copies repetitions of a literal motif. */
	private static String repeat(String s, int copies){
		final StringBuilder out=new StringBuilder(s.length()*copies);
		for(int i=0; i<copies; i++){out.append(s);}
		return out.toString();
	}

	/** Throws IllegalArgumentException on a failed fixture invariant, independent of -ea. */
	private static void require(boolean ok, String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}
}
