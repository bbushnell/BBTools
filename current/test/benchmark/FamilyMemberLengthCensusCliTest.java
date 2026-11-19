package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import benchmark.FamilyMemberLengthCensus;
import prot.DigestSuffix;
import structures.ByteBuilder;

/** Exercises the pinned census command with literal final-table oracles.
 * Run on Dori with out=NEW_DIRECTORY; retains the tiny inputs and outputs for review.
 * Expected statistics are written independently, never calculated by census helpers.
 * @author Brian, Keqing
 */
public final class FamilyMemberLengthCensusCliTest{

	/** Creates a fresh output directory and checks literal census tables, threshold states and input rejection. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1 && args[0].startsWith("out=") && args[0].length()>4,
			"Expected out=NEW_DIRECTORY for retained census CLI fixture evidence");
		final Path root=Paths.get(args[0].substring(4)).toAbsolutePath();
		check(!Files.exists(root), "Fixture root must be fresh; existing evidence must not be overwritten: "+root);
		Files.createDirectory(root);
		mixedRange(root.resolve("mixed"));
		legacyRange(root.resolve("legacy"));
		thresholdModes(root.resolve("thresholds"));
		modernThresholdModes(root.resolve("modern_thresholds"));
		System.out.println("FAMILY_MEMBER_LENGTH_CENSUS_CLI_TEST_PASS cases=6 exact_tables=16 stale_pin_rejected=true"
			+" structural_default=true unavailable_explicit=true malformed_state_rejected=true");
	}

	/** Two buckets and two batches split families across chunks; empty ranks remain visible.
	 * The 40/64/99 family has median64 but mean203/3; the next families straddle65/100.
	 */
	private static void mixedRange(final Path root) throws Exception{
		Files.createDirectory(root);
		write(root.resolve("families.tsv"), FAMILY_HEADER+
			"0\tf0\t3\n1\tf1\t3\n2\tf2\t3\n3\tf3\t2\n4\tf4\t2\n5\tempty\t0\n");
		// Deliberately unrelated consensus lengths catch accidental consensus-based statistics.
		write(root.resolve("consensus.faa"), ">f0\n"+sequence(120)+"\n>f1\n"+sequence(120)+
			"\n>f2\n"+sequence(120)+"\n>f3\n"+sequence(120)+"\n>f4\n"+sequence(120)+"\n>empty\nX\n");
		final String chunks=CHUNK_HEADER+
			chunk(root, "p00.tsv", "b00", 0, "g00", new int[]{0, 2, 4}, new int[]{40, 65, 100}, 205)+
			chunk(root, "p01.tsv", "b00", 1, "g01", new int[]{1, 3}, new int[]{64, 99}, 163)+
			chunk(root, "p10.tsv", "b01", 0, "g10", new int[]{0, 0, 2, 2, 4}, new int[]{64, 99, 65, 65, 100}, 393)+
			chunk(root, "p11.tsv", "b01", 1, "g11", new int[]{1, 1, 3}, new int[]{65, 65, 100}, 230);
		write(root.resolve("chunks.tsv"), chunks);
		final Path config=configuration(root, 6, 13, 2, 2, "minmemberaa=40\nmaxmemberaa=200\n");
		final Path out=root.resolve("result");
		run(config, out);
		same(out.resolve("family_mean_lengths.tsv"), MEAN_HEADER+
			"f0\t3\t203\t67.666666667\nf2\t3\t195\t65\nf4\t2\t200\t100\n"+
			"f1\t3\t194\t64.666666667\nf3\t2\t199\t99.500000000\nempty\t0\t0\tNA\n");
		same(out.resolve("family_lengths.tsv"), LENGTH_HEADER+
			"f0\t3\t40\t64\t99\t3\nf2\t3\t65\t65\t65\t3\nf4\t2\t100\t100\t100\t2\n"+
			"f1\t3\t64\t65\t65\t3\nf3\t2\t99\t99.5\t100\t2\nempty\t0\tNA\tNA\tNA\t0\n");
		same(out.resolve("summary.tsv"), "metric\tcount\nmembers_40_59\t1\nmembers_60_79\t7\nmembers_80_99\t2\n"+
			"members_100_124\t3\nmembers_125_149\t0\nmembers_150_199\t0\nmembers_200\t0\n"+
			"families_with_at_least_5_members_below_200\t0\nfamilies\t6\nmembers\t13\n");
		same(out.resolve("quantiles.tsv"), QUANTILES+"minmemberaa\t40\nmaxmemberaa\t200\n");
		same(out.resolve("PASS"), "FAMILY_MEMBER_LENGTH_CENSUS_PASS families=6 members=13\n");
		structural(root, config, out);
		final String pin=DigestSuffix.file(config.toString());
		final Path badOut=root.resolve("stale_pin_result");
		boolean rejected=false;
		try{
			FamilyMemberLengthCensus.main(new String[]{"in="+config, "insha80="+
				(pin.charAt(0)=='0' ? "1" : "0")+pin.substring(1), "verification=sha80", "out="+badOut});
		}catch(IllegalArgumentException expected){
			check("Preparation config fingerprint differs".equals(expected.getMessage()),
				"Stale pin must fail for its fingerprint, not an unrelated fixture error: "+expected.getMessage());
			rejected=true;
		}
		check(rejected && !Files.exists(badOut), "A stale config pin must fail before publishing any census directory");
	}

	/** Omitted member bounds preserve the historical100..200 reader and four-bin outputs. */
	private static void legacyRange(final Path root) throws Exception{
		Files.createDirectory(root);
		write(root.resolve("families.tsv"), FAMILY_HEADER+"0\tf0\t3\n1\tempty\t0\n");
		write(root.resolve("consensus.faa"), ">f0\n"+sequence(120)+"\n>empty\nX\n");
		write(root.resolve("chunks.tsv"), CHUNK_HEADER+
			chunk(root, "p00.tsv", "b00", 0, "legacy", new int[]{0, 0, 0}, new int[]{100, 150, 200}, 450));
		final Path config=configuration(root, 2, 3, 1, 1, "");
		final Path out=root.resolve("result");
		run(config, out);
		same(out.resolve("family_mean_lengths.tsv"), MEAN_HEADER+"f0\t3\t450\t150\nempty\t0\t0\tNA\n");
		same(out.resolve("family_lengths.tsv"), LENGTH_HEADER+"f0\t3\t100\t150\t200\t2\nempty\t0\tNA\tNA\tNA\t0\n");
		same(out.resolve("summary.tsv"), "metric\tcount\nmembers_100_124\t1\nmembers_125_149\t0\n"+
			"members_150_199\t1\nmembers_200\t1\nfamilies_with_at_least_5_members_below_200\t0\nfamilies\t2\nmembers\t3\n");
		same(out.resolve("quantiles.tsv"), QUANTILES);
		same(out.resolve("PASS"), "FAMILY_MEMBER_LENGTH_CENSUS_PASS families=2 members=3\n");
		structural(root, config, out);
	}

	/** Unknown digests are explicit input data, never forged hexadecimal placeholders. */
	private static void structural(final Path root, final Path original, final Path baseline) throws Exception{
		final Path chunks=root.resolve("structural_chunks.tsv"), config=root.resolve("structural.config");
		final StringBuilder rows=new StringBuilder(CHUNK_HEADER);
		final java.util.List<String> lines=Files.readAllLines(root.resolve("chunks.tsv"), StandardCharsets.US_ASCII);
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t", -1);
			check(fields.length==6, "Fixture chunk manifest must retain six fields");
			fields[3]="unmeasured";
			rows.append(String.join("\t", fields)).append('\n');
		}
		write(chunks, rows.toString());
		String recipe=new String(Files.readAllBytes(original), StandardCharsets.US_ASCII);
		recipe=recipe.replace("chunks=chunks.tsv\n", "chunks=structural_chunks.tsv\n")
			.replaceAll("(?m)^(familiessha80|chunkssha80|consensussha80)=[0-9a-f]{20}$", "$1=unmeasured");
		write(config, recipe);
		final Path out=root.resolve("structural_result");
		FamilyMemberLengthCensus.main(new String[]{"in="+config, "insha80=unmeasured", "verification=structural", "out="+out});
		final Path defaultOut=root.resolve("default_result"), inheritedOut=root.resolve("inherited_result");
		FamilyMemberLengthCensus.main(new String[]{"in="+config, "out="+defaultOut});
		FamilyMemberLengthCensus.main(new String[]{"in="+config, "insha80=00000000000000000000", "out="+inheritedOut});
		for(String file:new String[]{"family_mean_lengths.tsv", "family_lengths.tsv", "summary.tsv", "quantiles.tsv", "PASS"}){
			check(Arrays.equals(Files.readAllBytes(baseline.resolve(file)), Files.readAllBytes(out.resolve(file))),
				"Structural input handling must preserve the literal baseline table: "+file);
			check(Arrays.equals(Files.readAllBytes(baseline.resolve(file)), Files.readAllBytes(defaultOut.resolve(file))) &&
				Arrays.equals(Files.readAllBytes(baseline.resolve(file)), Files.readAllBytes(inheritedOut.resolve(file))),
				"Default structural mode must neither require nor rehash an inherited identity: "+file);
		}
		reject(config, "unmeasured", "sha80", root.resolve("unknown_hash"), "fingerprint differs");
		reject(config, "invalid", "structural", root.resolve("bad_label"), "Invalid structural preparation identity");
		reject(config, "unmeasured", "unknown", root.resolve("bad_policy"), "Unknown member census verification policy");
		final Path link=root.resolve("linked.config");
		Files.createSymbolicLink(link, config);
		reject(link, "unmeasured", "structural", root.resolve("linked_result"), "absolute nonsymlink");
	}
	/** Requires the named verification failure before any result directory is published. */
	private static void reject(final Path config, final String pin, final String policy, final Path out, final String reason) throws Exception{
		boolean failed=false;
		try{FamilyMemberLengthCensus.main(new String[]{"in="+config, "insha80="+pin, "verification="+policy, "out="+out});}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason), "Unexpected census rejection: "+expected);
			failed=true;
		}
		check(failed && !Files.exists(out), "Invalid census input must not publish a result: "+reason);
	}

	/** One legacy family has literal lower40/upper50 gates; verification must not change its summary. */
	private static void thresholdModes(final Path root) throws Exception{
		Files.createDirectory(root);
		final Path reference=root.resolve("references.tsv");
		final ByteBuilder bytes=new ByteBuilder("#family_reference\tfamily_reference_v7\t5\t1024\t0\t1\t1\t0\t"
			+"k5_identity30_blosum_uncapped_v4\ta48_f32\t3.8906\n"
			+"#gates\tcosine_lower\tidentity_lower\tblosum_lower\tcosine_accept\tidentity_upper\tblosum_final\tscalar_mask"
			+"\thydro_min\thydro_max\tcharge_min\tcharge_max\taromatic_min\taromatic_max\tsmall_min\tsmall_max\tpi_min\tpi_max\n"
			+"#consensus\tamino_acids\n#ddl\t1024\tA48\n#spectrum\t400\ta48_f32\n#family\tf0\t0\t1\t100\t96\n");
		for(float value:new float[]{0, 40, 0, 1, 50, 0}){bytes.appendFloatA48(value).tab();}
		bytes.append(0);
		for(int i=0; i<10; i++){bytes.tab().append(0);}
		bytes.nl().append(sequence(100)).nl();
		for(int i=0; i<1024; i++){
			if(i>0){bytes.tab();}
			bytes.appendA48(i==0 ? 1 : 0);
		}
		bytes.nl();
		for(int i=0; i<400; i++){
			if(i>0){bytes.tab();}
			bytes.appendFloatA48(i==0 ? 1f : 0f);
		}
		bytes.nl();
		write(reference, bytes.toString());
		final Path out=root.resolve("default"), pinned=root.resolve("pinned");
		FamilyMemberLengthCensus.main(new String[]{"reference="+reference, "marginscale=1", "out="+out});
		FamilyMemberLengthCensus.main(new String[]{"reference="+reference, "marginscale=1", "out="+pinned,
			"verification=sha80", "insha80="+DigestSuffix.file(reference.toString())});
		same(out.resolve("thresholds.tsv"), THRESHOLD_HEADER
			+"member_identity_lower\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t40\t40\t40\t40"
			+"\tglocal_full_query\t0\t0\tblosum62\n"
			+"family_identity_upper\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t50\t50\t50\t50"
			+"\tglocal_full_query\t0\t0\tblosum62\n"
			+"effective_identity_accept_floor30\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t50\t50\t50\t50"
			+"\tglocal_full_query\t0\t0\tblosum62\n");
		same(out.resolve("PASS"), "FAMILY_IDENTITY_THRESHOLDS_PASS families=1 measured=1 unavailable=0 empty=0"
			+" status=ALL_LIVE_FAMILIES_MEASURED quantile=nearest_rank floor=30\n");
		for(String file:new String[]{"thresholds.tsv", "PASS"}){
			check(Arrays.equals(Files.readAllBytes(out.resolve(file)), Files.readAllBytes(pinned.resolve(file))),
				"Explicit exact-match verification must preserve the threshold summary: "+file);
		}
		final Path rejected=root.resolve("missing_pin");
		boolean failed=false;
		try{FamilyMemberLengthCensus.main(new String[]{"reference="+reference, "marginscale=1",
			"verification=sha80", "out="+rejected});}
		catch(IllegalArgumentException expected){
			check("Reference fingerprint differs".equals(expected.getMessage()), "Explicit SHA80 requires a known pin");
			failed=true;
		}
		check(failed && !Files.exists(rejected), "Missing exact-match identity must fail before publishing thresholds");
	}

	/** Literal v16 inputs exercise state counts and even-population nearest-rank medians.
	 * The empty rank is separate from the unavailable live family; its zeros are not observations.
	 */
	private static void modernThresholdModes(final Path root) throws Exception{
		Files.createDirectory(root);
		final String mixed=modernReference(new String[]{"alignment_unavailable", "alignment_enabled", "alignment_enabled",
			"alignment_enabled", "alignment_enabled"}, new int[]{0, 20, 40, 60, 80}, new int[]{0, 40, 60, 80, 100}, true);
		final Path reference=root.resolve("mixed.tsv");
		write(reference, mixed);
		final Path out=root.resolve("mixed_result");
		thresholdRun(reference, out, "0.5", true);
		same(out.resolve("thresholds.tsv"), THRESHOLD_HEADER
			+"member_identity_lower\t6\t4\t1\t1\tUNEXPECTED_UNAVAILABLE\t0.500000\t20\t40\t80\t80"+MODERN_METRIC
			+"family_identity_upper\t6\t4\t1\t1\tUNEXPECTED_UNAVAILABLE\t0.500000\t40\t60\t100\t100"+MODERN_METRIC
			+"effective_identity_accept_floor30\t6\t4\t1\t1\tUNEXPECTED_UNAVAILABLE\t0.500000\t35\t50\t90\t90"+MODERN_METRIC);
		same(out.resolve("PASS"), "FAMILY_IDENTITY_THRESHOLDS_PASS families=5 measured=4 unavailable=1 empty=1"
			+" status=UNEXPECTED_UNAVAILABLE quantile=nearest_rank floor=30\n");
		final Path none=root.resolve("none.tsv"), noneOut=root.resolve("none_result");
		write(none, modernReference(new String[]{"alignment_unavailable"}, new int[]{0}, new int[]{0}, false));
		thresholdRun(none, noneOut, "1", true);
		same(noneOut.resolve("thresholds.tsv"), THRESHOLD_HEADER
			+"member_identity_lower\t1\t0\t1\t0\tUNEXPECTED_UNAVAILABLE\t1\tNA\tNA\tNA\tNA"+MODERN_METRIC
			+"family_identity_upper\t1\t0\t1\t0\tUNEXPECTED_UNAVAILABLE\t1\tNA\tNA\tNA\tNA"+MODERN_METRIC
			+"effective_identity_accept_floor30\t1\t0\t1\t0\tUNEXPECTED_UNAVAILABLE\t1\tNA\tNA\tNA\tNA"+MODERN_METRIC);
		final Path all=root.resolve("all.tsv"), allOut=root.resolve("all_result");
		write(all, modernReference(new String[]{"alignment_enabled"}, new int[]{40}, new int[]{50}, false));
		thresholdRun(all, allOut, "1", false);
		same(allOut.resolve("thresholds.tsv"), THRESHOLD_HEADER
			+"member_identity_lower\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t40\t40\t40\t40"+MODERN_METRIC
			+"family_identity_upper\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t50\t50\t50\t50"+MODERN_METRIC
			+"effective_identity_accept_floor30\t1\t1\t0\t0\tALL_LIVE_FAMILIES_MEASURED\t1\t50\t50\t50\t50"+MODERN_METRIC);
		thresholdReject(root, "unknown_state", mixed.replace("alignment_unavailable", "alignment_unknown"),
			"Unknown alignment calibration availability");
		thresholdReject(root, "inert_gate", modernReference(new String[]{"alignment_unavailable"},
			new int[]{10}, new int[]{0}, false), "inert zeros");
		thresholdReject(root, "missing_state", mixed.replace("\talignment_unavailable\t", "\t"),
			"Saved gate fields must match");
		thresholdReject(root, "count_mismatch", mixed.replace("\t6\t5\t1\t", "\t6\t4\t2\t"),
			"bounded live roster");
	}

	/** The fixture spells the wire schema independently; production codec writers are not its oracle. */
	private static String modernReference(final String[] states, final int[] lower, final int[] upper,
			final boolean empty){
		check(states.length>0 && states.length==lower.length && lower.length==upper.length,
			"Threshold fixture needs paired state and gate rows");
		final ByteBuilder bytes=new ByteBuilder("#family_reference\tfamily_reference_v16\t6\t512\t0\t");
		bytes.append(states.length+(empty ? 1 : 0)).tab().append(states.length).tab().append(empty ? 1 : 0)
			.append("\tk5_identity30_blosum_uncapped_v4\ta48_f32\t0\tpin_right_trim_left_aligned_query\t40\tblosum62\t20\n")
			.append("#gates\tcosine_lower\tidentity_lower\tblosum_lower\tcosine_accept\tidentity_upper\tblosum_final")
			.append("\tscalar_mask\talignment_state\thydro_min\thydro_max\tcharge_min\tcharge_max")
			.append("\taromatic_min\taromatic_max\tsmall_min\tsmall_max\tpi_min\tpi_max\n")
			.append("#consensus\tamino_acids\n#ddl\t512\tA48\n#spectrum\t20\ta48_f32\n");
		for(int rank=0; rank<states.length; rank++){
			bytes.append("#family\tf").append(rank).tab().append(rank).append("\t1\t100\t95\n");
			for(float value:new float[]{0, lower[rank], 0, 1, upper[rank], 0}){bytes.appendFloatA48(value).tab();}
			bytes.append(0).tab().append(states[rank]);
			for(int i=0; i<10; i++){bytes.tab().append(0);}
			bytes.nl().append(sequence(100)).nl();
			for(int i=0; i<512; i++){
				if(i>0){bytes.tab();}
				bytes.appendA48(i==0 ? 1 : 0);
			}
			bytes.nl();
			for(int i=0; i<20; i++){
				if(i>0){bytes.tab();}
				bytes.appendFloatA48(i==0 ? 1f : 0f);
			}
			bytes.nl();
		}
		if(empty){bytes.append("#family\tempty\t").append(states.length).append("\t0\t0\t0\nNA\nNA\nNA\nNA\n");}
		return bytes.toString();
	}

	/** Warnings are observable output: nonzero unavailable is unexpected; zero is explicitly reported. */
	private static void thresholdRun(final Path reference, final Path out, final String scale,
			final boolean warning) throws Exception{
		final java.io.PrintStream previous=System.err;
		final java.io.ByteArrayOutputStream captured=new java.io.ByteArrayOutputStream();
		try(java.io.PrintStream stream=new java.io.PrintStream(captured, true, "UTF-8")){
			System.setErr(stream);
			FamilyMemberLengthCensus.main(new String[]{"reference="+reference, "marginscale="+scale, "out="+out});
		}finally{System.setErr(previous);}
		final String diagnostic=new String(captured.toByteArray(), StandardCharsets.UTF_8);
		check(warning ? diagnostic.contains("WARNING:") && diagnostic.contains("this is unexpected") :
			diagnostic.contains("0 unavailable families") && !diagnostic.contains("WARNING:"),
			"Threshold availability must be reported plainly: "+diagnostic);
	}

	/** Writes a deliberately malformed threshold input and requires its specific rejection without output. */
	private static void thresholdReject(final Path root, final String name, final String input,
			final String reason) throws Exception{
		final Path reference=root.resolve(name+".tsv"), out=root.resolve(name+"_result");
		write(reference, input);
		boolean rejected=false;
		try{FamilyMemberLengthCensus.main(new String[]{"reference="+reference, "marginscale=1", "out="+out});}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason),
				"Malformed threshold record failed for the wrong reason: "+expected);
			rejected=true;
		}
		check(rejected && !Files.exists(out), "Malformed thresholds must fail before publishing output: "+name);
	}

	/** Explicit chunk totals are independent input assertions, not recomputed expected output. */
	private static String chunk(final Path root, final String name, final String batch, final int bucket,
			final String genome, final int[] ranks, final int[] lengths, final int residues) throws Exception{
		check(ranks.length==lengths.length && ranks.length>0, "Fixture chunk needs paired family and member-length entries");
		final StringBuilder rows=new StringBuilder(MEMBER_HEADER);
		for(int i=0; i<ranks.length; i++){
			rows.append(ranks[i]).append("\tf").append(ranks[i]).append('\t').append(genome).append('|').append(i)
				.append('\t').append(genome).append('\t').append(i).append('\t').append(sequence(lengths[i])).append('\n');
		}
		final Path part=root.resolve(name);
		write(part, rows.toString());
		return batch+"\t"+bucket+"\t"+part+"\t"+DigestSuffix.file(part.toString())+"\t"+ranks.length+"\t"+residues+"\n";
	}

	/** Relative manifest paths exercise the command's config-directory resolution. */
	private static Path configuration(final Path root, final int families, final int members,
			final int buckets, final int batches, final String bounds) throws Exception{
		final Path config=root.resolve("preparation.config");
		write(config, "schema=family_preparation_v1\nfamilies=families.tsv\nfamiliessha80="+
			DigestSuffix.file(root.resolve("families.tsv").toString())+"\nchunks=chunks.tsv\nchunkssha80="+
			DigestSuffix.file(root.resolve("chunks.tsv").toString())+"\nconsensus=consensus.faa\nconsensussha80="+
			DigestSuffix.file(root.resolve("consensus.faa").toString())+"\nexpectedfamilies="+families+
			"\nexpectedmembers="+members+"\nbuckets="+buckets+"\nbatches="+batches+
			"\nmaxreferenceaa=200\nmaxbucketmembers=32\nmaxbucketbytes=16384\n"+bounds);
		return config;
	}
	/** Executes the exact-verification census path against an existing config and a fresh output path. */
	private static void run(final Path config, final Path out) throws Exception{
		check(Files.isRegularFile(config) && !Files.exists(out), "CLI case needs pinned inputs and a fresh result path");
		FamilyMemberLengthCensus.main(new String[]{"in="+config, "insha80="+DigestSuffix.file(config.toString()),
			"verification=sha80", "out="+out});
	}
	/** Produces an alanine homopolymer within the fixture's supported 40–200-residue range. */
	private static String sequence(final int length){
		check(length>=40 && length<=200, "Literal fixture sequence length must fit the intended reader contract");
		final char[] residues=new char[length];
		Arrays.fill(residues, 'A');
		return new String(residues);
	}
	/** Writes a new ASCII fixture without replacing evidence from an earlier case. */
	private static void write(final Path path, final String text) throws Exception{
		check(!Files.exists(path), "Fixture inputs must not overwrite retained evidence: "+path);
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	/** Compares the entire output to a literal ASCII oracle, including row order and formatting. */
	private static void same(final Path path, final String expected) throws Exception{
		final byte[] actual=Files.readAllBytes(path);
		check(Arrays.equals(expected.getBytes(StandardCharsets.US_ASCII), actual),
			"Literal census output differs at "+path+"\nExpected:\n"+expected+"Actual:\n"+new String(actual, StandardCharsets.US_ASCII));
	}
	/** Fails immediately with the fixture-specific explanation when an expectation is violated. */
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}

	/** Literal roster schema, including empty-family rows. */
	private static final String FAMILY_HEADER="roster_rank\tfamily_id\tassigned\n";
	/** Literal threshold summary schema with availability counts and alignment metadata. */
	private static final String THRESHOLD_HEADER="metric\ttotal_families\tmeasured_families\tunavailable_families\tempty_families"
		+"\tavailability\tmargin_scale\tmin\tmedian\tp90\tmax\talignment_mode\tminimum_span\tminimum_paired\tmatrix_name\n";
	/** Alignment metadata appended to every modern threshold oracle row. */
	private static final String MODERN_METRIC="\tpin_right_trim_left_aligned_query\t40\t20\tblosum62\n";
	/** Chunk manifest schema carrying identities and independently supplied totals. */
	private static final String CHUNK_HEADER="batch_id\tbucket\tfile\tsha80\tmembers\tresidues\n";
	/** Literal retained-member schema consumed by the census fixture. */
	private static final String MEMBER_HEADER="family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence\n";
	/** Exact output schema for arithmetic member-length means. */
	private static final String MEAN_HEADER="family\tmembers\tlength_sum\tlength_mean\n";
	/** Exact output schema for family member-length quantiles and counts. */
	private static final String LENGTH_HEADER="family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n";
	/** Literal quantile-policy output; member peptides, not consensus lengths, supply observations. */
	private static final String QUANTILES="p10_p90\tnearest_rank\nmedian\tmean_of_two_middle_members_if_even\n"+
		"empty_family\tNA_quantiles_with_zero_members\nlength_source\tretained_member_peptides_not_consensus_or_original_gene\n";
}
