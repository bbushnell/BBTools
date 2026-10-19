package test.prok;

import prok.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Literal label-table oracle: source order differs, a middle family is omitted,
 * and UTF-8 plus an unknown column must survive unchanged. No HBM claim is made.
 * @author Keqing
 */
public final class FamilyLabelsSubsetTest {
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Expected a fresh fixture directory");
		final Path root=Files.createDirectory(Paths.get(args[0]));
		FamilyLabelsSubset.validateOutput(root.resolve("new.tsv"));
		for(String name : new String[]{"stdout.tsv", "stderr.tsv", "/dev/null", "out.tsv.gz", "out.tsv.GZ", "out.tsv.bz2", "out.tsv.xz", "out|alias.tsv"}){
			final Path rejected=name.startsWith("out") ? root.resolve(name) : Paths.get(name);
			try{FamilyLabelsSubset.validateOutput(rejected);}
			catch(IllegalArgumentException expected){continue;}
			throw new AssertionError("Expected special/compressed output rejection: "+name);
		}
		final String header="family_id\tactive_index\trep_id\tlabel_version\tdecision_class\tproduct\textra\n";
		final String first="10\t0\ta\tv1\tspecific_full_length\tβ protein; A=B\tretained α\n";
		final String third="30\t2\tc\tv1\tunnamed\t\ttrailing\n";
		final Path input=root.resolve("parent.tsv"), output=root.resolve("subset.tsv");
		Files.write(input, (meta(OLD, OLD)+header+third
			+"20\t1\tb\tv1\tconflicting_evidence\t\tomitted\n"+first).getBytes(StandardCharsets.UTF_8));
		FamilyLabels.load(input.toString(), new int[]{10, 20, 30}, new String[]{"a", "b", "c"}, OLD, OLD);
		final byte[] before=Files.readAllBytes(input);
		final int rows=FamilyLabelsSubset.write(input, output, 2, IDENTITY, CONSENSUS, new int[]{0, -1, 1});
		final String expected=meta(IDENTITY, CONSENSUS)+header+first+"30\t1\tc\tv1\tunnamed\t\ttrailing\n";
		check(rows==2 && expected.equals(new String(Files.readAllBytes(output), StandardCharsets.UTF_8)),
			"Subset label bytes differ from the independent literal table");
		check(java.util.Arrays.equals(before, Files.readAllBytes(input)), "Source labels changed");
		final FamilyLabels labels=FamilyLabels.load(output.toString(), new int[]{10, 30}, new String[]{"a", "c"}, IDENTITY, CONSENSUS);
		check("β protein; A=B".equals(labels.get(0).product) && labels.get(1).product==null,
			"Transferred labels do not resolve to the expected local families");
		System.err.println("FAMILY_LABELS_SUBSET_TEST_PASS literal_bytes=true reordered=true omitted=true utf8=true unknown_column=true source_unchanged=true output_rejections=8");
	}
	private static String meta(final String identity, final String consensus){
		return "#schema_version\t1\n#kind\tbbgap_family_labels\n#label_version\tv1\n#family_identity_sha80\t"
			+identity+"\n#consensus_sha80\t"+consensus+"\n";
	}
	private static void check(final boolean ok, final String reason){if(!ok){throw new AssertionError(reason);}}
	private static final String OLD="00000000000000000000", IDENTITY="11111111111111111111", CONSENSUS="22222222222222222222";
}
