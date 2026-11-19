package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import ddl.DDLFamilySearch.StructuralInputs;
import prot.FamilyGcInputs;
import prot.FamilyGcPenalty;
import prot.FamilyReference;
import prot.FamilyCompactCodec;
import prot.ProteinSequence;

/**
 * Checks literal GC report/map joins, embedded compact metadata and rejected inputs.
 * Creates and retains evidence in one fresh directory; the literal tables define
 * three families and six members, including a singleton with no sample SD.
 * Run once per JVM: the check counter and evidence list are shared mutable state.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class FamilyGcInputsTest{

	/** Runs the file/compact joins and rejection cases, then writes PASS and evidence.paths.
	 * @param args One fresh, absolute, normalized evidence directory with an existing parent */
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply one fresh absolute evidence directory");
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && root.normalize().equals(root) && !Files.exists(root),
			"Fixture directory must be fresh and normalized");
		Files.createDirectory(root);
		final StructuralInputs tracked=new StructuralInputs();
		final FamilyGcPenalty p=load(root, "literal", REPORT, MAP, 3, 6, tracked, null);
		check(p.families()==3 && p.stableIdAt(1)==19, "Literal complete roster");
		p.requireRoster(new int[]{0, 1, 2}, new String[]{"f7", "f19", "f42"}, new long[]{3, 2, 1});
		check(p.score(10, 0, true, .25)==9, "Literal loaded GC penalty");
		check(Double.doubleToLongBits(p.meanAt(1))==Double.doubleToLongBits(34.0/55),
			"Exact double roundtrip from final report");
		check(p.meanAvailableAt(2) && p.score(8, 2, true, 0)==4,
			"Singleton measured mean remains usable despite unavailable sample SD");
		final FamilyReference[] compact={family("f7", 0, 3), family("f19", 1, 2), family("f42", 2, 1)};
		final FamilyGcPenalty embedded=FamilyGcInputs.fromCompact(compact, new int[]{7, 19, 42}, new long[]{3, 2, 1},
			new double[]{.5, 34.0/55, 1}, new double[]{.25, 0, 0}, 24, true, 8);
		embedded.requireRoster(new int[]{0, 1, 2}, new String[]{"f7", "f19", "f42"}, new long[]{3, 2, 1});
		check(embedded.families()==3 && embedded.stableIdAt(1)==19 &&
			Double.doubleToLongBits(embedded.meanAt(1))==Double.doubleToLongBits(34.0/55),
			"Embedded compact GC arrays bind the same live roster and exact means");
		final FamilyReference[] swapped={compact[1], compact[0], compact[2]};
		boolean shuffled=false;
		try{
			FamilyGcInputs.fromCompact(swapped, new int[]{7, 19, 42}, new long[]{3, 2, 1},
				new double[]{.5, 34.0/55, 1}, new double[]{.25, 0, 0}, 24, true, 8)
				.requireRoster(new int[]{0, 1, 2}, new String[]{"f7", "f19", "f42"}, new long[]{3, 2, 1});
		}catch(IllegalArgumentException expected){shuffled=true;}
		check(shuffled, "Compact GC rejects shuffled embedded live order");
		tracked.verify();
		Files.write(root.resolve("literal/map.tsv"), "\n".getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.APPEND);
		boolean changed=false;
		try{tracked.verify();}catch(IllegalStateException expected){changed=true;}
		check(changed, "Caller tracker retains the loaded files for later stability checks");

		load(root, "wrong_stable", REPORT, field(MAP, 1, 0, "8"), 3, 6, new StructuralInputs(),
			"differs from the accepted collection map");
		load(root, "wrong_rank", REPORT, field(MAP, 1, 1, "1"), 3, 6, new StructuralInputs(),
			"differs from the accepted collection map");
		load(root, "wrong_label", REPORT, field(MAP, 1, 2, "other"), 3, 6, new StructuralInputs(),
			"differs from the accepted collection map");
		load(root, "wrong_count", REPORT, field(MAP, 1, 3, "4"), 3, 6, new StructuralInputs(),
			"differs from the accepted collection map");
		load(root, "wrong_population", REPORT, MAP, 3, 7, new StructuralInputs(), "declared family/member cohort");
		load(root, "short_report", dropLast(REPORT), MAP, 3, 6, new StructuralInputs(),
			"families absent from the GC report");
		load(root, "short_map", REPORT, dropLast(MAP), 3, 6, new StructuralInputs(), "map ended before");
		load(root, "wrong_residues", REPORT, field(MAP, 1, 4, "119"), 3, 6, new StructuralInputs(), "SHORT population");
		load(root, "padded_mean", field(REPORT, 1, 5, "0.500000000000"), MAP, 3, 6, new StructuralInputs(),
			"Double.toString");
		load(root, "integral_sd", field(REPORT, 2, 6, "0"), MAP, 3, 6, new StructuralInputs(), "Double.toString");
		load(root, "nan_mean", field(REPORT, 1, 5, "NaN"), MAP, 3, 6, new StructuralInputs(), "Double.toString");
		load(root, "missing_mean", field(REPORT, 1, 4, "false"), MAP, 3, 6, new StructuralInputs(),
			"complete GC report");
		load(root, "singleton_sd", field(REPORT, 3, 6, "0.0"), MAP, 3, 6, new StructuralInputs(), "Singleton means");
		load(root, "bad_status", field(REPORT, 2, 8, "available"), MAP, 3, 6, new StructuralInputs(),
			"gate availability/status");
		load(root, "nonnested", field(REPORT, 1, 11, "3"), MAP, 3, 6, new StructuralInputs(), "bounded and nested");
		load(root, "constant_exclusion", field(REPORT, 2, 9, "1"), MAP, 3, 6, new StructuralInputs(),
			"cannot be excluded");
		load(root, "duplicate_stable", field(REPORT, 2, 0, "7"), field(MAP, 2, 0, "7"), 3, 6, new StructuralInputs(),
			"Duplicate or negative stable");
		final String marker="FAMILY_GC_INPUTS_TEST_PASS checks="+checks;
		write(root.resolve("PASS"), marker+"\n");
		evidence.add(root.resolve("evidence.paths").toString());
		Files.write(root.resolve("evidence.paths"), evidence, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		System.out.println(marker);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Case construction       ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes one case and checks its declared success or diagnostic substring.
	 * A null diagnostic requires success; an expected rejection returns null.
	 * The supplied tracker retains the successfully loaded inputs for later verification. */
	private static FamilyGcPenalty load(Path root, String name, String report, String map, int families, long members,
			StructuralInputs tracked, String diagnostic)throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		write(dir.resolve("report.tsv"), report);
		write(dir.resolve("map.tsv"), map);
		FamilyGcPenalty result=null;
		boolean rejected=false;
		try{
			result=FamilyGcInputs.load(dir.resolve("report.tsv"), dir.resolve("map.tsv"), families, members, 4,
				tracked);
		}catch(IllegalArgumentException e){
			check(diagnostic!=null && e.getMessage().contains(diagnostic), "Wrong rejection: "+e);
			rejected=true;
		}
		check(rejected==(diagnostic!=null), "Expected success/rejection for "+name);
		return result;
	}

	/** Replaces one zero-based cell, counting the header as row zero, and restores the final newline. */
	private static String field(String text, int row, int column, String value){
		final String[] rows=text.split("\n"), fields=rows[row].split("\t", -1);
		fields[column]=value;
		rows[row]=String.join("\t", fields);
		return String.join("\n", rows)+"\n";
	}

	/** Drops the final row of a newline-terminated literal table, retaining the preceding newline. */
	private static String dropLast(String text){return text.substring(0, text.lastIndexOf('\n', text.length()-2)+1);}
	/** Creates one ASCII evidence file without overwriting it and records its absolute path. */
	private static void write(Path path, String text)throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		evidence.add(path.toString());
	}

	/** Constructs a minimal compact reference with the requested roster identity and member count. */
	private static FamilyReference family(String id, int rank, long members){
		final byte[] raw="ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII);
		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		final float[] monomers=new float[FamilyCompactCodec.MONOMER_WIDTH];
		monomers[0]=1;
		return new FamilyReference(id, rank, members, raw.length*members, 36*members,
			new ProteinSequence(id, raw), ddl, monomers, 0, 0, 0, 0.72f, 30, 0, raw, 0, new double[10]);
	}

	/** Counts each attempted assertion and throws AssertionError on failure, independent of -ea. */
	private static void check(boolean valid, String message){
		checks++;
		if(!valid){throw new AssertionError(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Evidence and literals     ----------------*/
	/*--------------------------------------------------------------*/

	/** Attempted assertions in this JVM; main does not reset the count. */
	private static int checks;
	/** Created evidence paths in insertion order; main adds the list itself before writing it. */
	private static final ArrayList<String> evidence=new ArrayList<String>();
	/** Accepted literal report: ordinary variance, zero variance and singleton rows. */
	private static final String REPORT="stable_family_id\troster_rank\tfamily_id\tmembers\tmean_available\tmean_gc\tsample_sd\tgate_available\tgate_status"+
		"\texcluded_N0\texcluded_N0.5\texcluded_N1\texcluded_N1.5\texcluded_N2\texcluded_N2.5\texcluded_N3\texcluded_N4\texcluded_N5\texcluded_N6\n"+
		"7\t0\tf7\t3\ttrue\t0.5\t0.25\ttrue\tavailable\t2\t2\t0\t0\t0\t0\t0\t0\t0\t0\n"+
		"19\t1\tf19\t2\ttrue\t0.6181818181818182\t0.0\ttrue\tzero_sd\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"+
		"42\t2\tf42\t1\ttrue\t1.0\tNA\tfalse\tsingleton_neutral\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
	/** Matching collection roster with literal member and residue counts. */
	private static final String MAP="stable_family_id\troster_rank\tfamily_id\tassigned\tresidues\n"+
		"7\t0\tf7\t3\t120\n19\t1\tf19\t2\t198\n42\t2\tf42\t1\t40\n";
}
