package prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/** Constructor-only reproduction of external model symlink and evidence hardlink
 * gaps. No annotation writer is invoked, even when checking vulnerable code.
 * @author Keqing
 */
public final class AnnotateGenesPackAliasTest {
	public static void main(final String[] args) throws Exception{
		check(args.length==2 && (args[1].equals("vulnerable") || args[1].equals("normalized") || args[1].equals("protected")),
			"Expected a fresh fixture directory and vulnerable|normalized|protected");
		AnnotateGenesPackInputTest.main(new String[]{args[0]});
		final Path base=Paths.get(args[0]), reference=base.resolve("reference");
		final Path pack=reference.resolve("packs/part3"), external=base.resolve("external_models");
		Files.move(pack.resolve("models"), external);
		Files.createSymbolicLink(pack.resolve("models"), external.toAbsolutePath());
		final Path outside=external.resolve("new.gff");
		final Path evidence=pack.resolve("evidence/decisions.tsv"), hardlink=base.resolve("evidence_alias.gff");
		Files.createLink(hardlink, evidence);
		final byte[] before=Files.readAllBytes(evidence);
		final boolean symlink=rejects(base, reference, outside, "immutable model package");
		final boolean alias=rejects(base, reference, hardlink, "reference pack payload");
		final boolean expected=!args[1].equals("vulnerable");
		check(symlink==expected && alias==expected, "Unexpected alias guards: model_symlink="+symlink+", evidence_hardlink="+alias);
		// A lexical normalize() before filesystem resolution changes the inode here.
		final Path deep=Files.createDirectories(pack.resolve("alt/deep"));
		final Path actual=Files.write(pack.resolve("alt/table.tsv"), before);
		Files.createSymbolicLink(pack.resolve("evidence/link"), deep.toAbsolutePath());
		final Path declared=pack.resolve("evidence/link/../table.tsv"), traversalAlias=base.resolve("traversal_alias.gff");
		Files.createLink(traversalAlias, actual);
		check(Files.isSameFile(declared, actual) && !Files.exists(declared.normalize()),
			"Fixture must distinguish filesystem symlink traversal from lexical normalization");
		Files.write(pack.resolve("manifest.tsv"), ("evidence\tevidence/link/../table.tsv\t1\t00000000000000000000\tfixture\n")
			.getBytes(java.nio.charset.StandardCharsets.US_ASCII), java.nio.file.StandardOpenOption.APPEND);
		final boolean traversal=rejects(base, reference, traversalAlias, "canonical relative");
		check(traversal==args[1].equals("protected"), "Unexpected symlink traversal guard: "+traversal);
		check(!Files.exists(outside) && java.util.Arrays.equals(before, Files.readAllBytes(evidence)),
			"Constructor-only regression must not modify any model/evidence bytes");
		check(java.util.Arrays.equals(before, Files.readAllBytes(actual)), "Traversal alias validation must preserve the actual payload");
		System.err.println("ANNOTATE_GENES_PACK_ALIAS_PASS expected="+args[1]+" model_symlink_rejected="+symlink
			+" evidence_hardlink_rejected="+alias+" symlink_traversal_rejected="+traversal+" inputs_unchanged=true");
	}
	private static boolean rejects(final Path base, final Path reference, final Path out, final String message) throws Exception{
		try{
			new AnnotateGenes(new String[]{base.resolve("genome.fna").toString(), base.resolve("genes.gff").toString(),
				"reference="+reference, "out="+out, "t=1"});
		}catch(IllegalArgumentException e){check(e.getMessage().contains(message), "Unexpected rejection: "+e); return true;}
		return false;
	}
	private static void check(final boolean ok, final String message){if(!ok){throw new AssertionError(message);}}
}
