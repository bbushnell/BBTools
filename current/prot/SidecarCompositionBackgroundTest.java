package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

/** Closed-form composition oracle and default-output preservation check.
 * @author Keqing
 */
public final class SidecarCompositionBackgroundTest {

	public static void main(final String[] args) throws Exception{
		check(args.length==2, "Expected fresh test root and pre-change selftest fixture");
		final Path root=Files.createDirectory(Paths.get(args[0]));
		final Path baseline=Paths.get(args[1]);
		FamilyShortlistSidecarBuilder.main(new String[]{"consensus="+baseline.resolve("consensus.fasta"),
			"familiesdir="+baseline.resolve("fam"), "kmersets="+baseline.resolve("sets.tsv"),
			"out="+root.resolve("default.tsv")});
		check(Arrays.equals(Files.readAllBytes(baseline.resolve("sidecar.tsv")),
			Files.readAllBytes(root.resolve("default.tsv"))), "Default sidecar changed from the pre-change output");

		final Path families=Files.createDirectory(root.resolve("families"));
		write(families.resolve("0.faa"), ">member\nAAAAAAAAAA\n");
		final Path consensus=root.resolve("consensus.faa"), kmers=root.resolve("kmers.tsv");
		write(consensus, ">rep active_index=0\nAAAAAAAAAA\n");
		write(kmers, "#schema_version\t1\n#kind\tcovering_set\n#bbtools_commit\ttest\n#alphabet\taa20\n#k\t5\nrep\tAAAAA\n");
		final ReducedAlphabet alphabet=ReducedAlphabet.named("aa20");
		final int aa=21*alphabet.code((byte)'A'), cc=21*alphabet.code((byte)'C');
		final StringBuilder vector=new StringBuilder("#background_400f\t");
		for(int i=0; i<400; i++){
			if(i>0){vector.append(',');}
			vector.append(i==aa ? "0.12345678912345678" : i==cc ? "0.8765432108765432" : "0");
		}
		vector.append('\n');
		final String content=HEADER+vector;
		final Path source=root.resolve("donor.tsv"), output=root.resolve("pinned.tsv");
		write(source, content);
		final String pin=DigestSuffix.file(source.toString());
		final SidecarCompositionBackground inherited=SidecarCompositionBackground.load(source.toString(), pin, "centered", 0, "families");
		check(inherited.frequencies[aa]==0.12345678912345678 && inherited.frequencies[cc]==0.8765432108765432,
			"Inherited probabilities differ from the consumer's decimal decoder");
		FamilyShortlistSidecarBuilder.main(new String[]{"consensus="+consensus, "familiesdir="+families,
			"kmersets="+kmers, "out="+output, "backgroundfrom="+source, "backgroundsha80="+pin});
		final String actual=new String(Files.readAllBytes(output), StandardCharsets.UTF_8);
		check(actual.contains(vector.toString()), "The donor background was rounded or rewritten");
		final FamilyShortlistSidecar loaded=FamilyShortlistSidecar.load(output.toString(), consensus.toString(), kmers.toString());
		check(loaded.background[aa]==0.12345678912345678 && loaded.background[cc]==0.8765432108765432,
			"Serialized background values changed on consumer load");
		// Pure A minus (p A + (1-p) C) normalizes to (+1,-1)/sqrt(2).
		for(int i=0; i<400; i++){
			final double expected=i==aa ? Math.sqrt(0.5) : i==cc ? -Math.sqrt(0.5) : 0;
			check(Math.abs(loaded.compositionMatrix[i]-expected)<1e-7, "Closed-form composition mismatch at "+i);
		}
		check(content.equals(new String(Files.readAllBytes(source), StandardCharsets.UTF_8)), "Donor source was changed");
		final Path precise=root.resolve("precise.tsv"), preciseOutput=root.resolve("precise_output.tsv");
		write(precise, content.replace("#composition_rawpseudo\t0", "#composition_rawpseudo\t0.123456789"));
		FamilyShortlistSidecarBuilder.main(new String[]{"consensus="+consensus, "familiesdir="+families,
			"kmersets="+kmers, "out="+preciseOutput, "backgroundfrom="+precise,
			"backgroundsha80="+DigestSuffix.file(precise.toString()), "rawpseudo=0.123456789"});
		check(FamilyShortlistSidecar.load(preciseOutput.toString(), consensus.toString(), kmers.toString()).rawPseudocount==0.123456789,
			"Inherited pseudocount was rounded during serialization");
		final Path infinite=root.resolve("infinite.tsv");
		write(infinite, content.replace("#composition_rawpseudo\t0", "#composition_rawpseudo\tInfinity"));
		try{SidecarCompositionBackground.load(infinite.toString(), DigestSuffix.file(infinite.toString()), "centered", Double.POSITIVE_INFINITY, "families");
			throw new AssertionError("Infinite pseudocount was accepted");
		}catch(IllegalArgumentException expected){}
		reject(source, "00000000000000000000", "centered");
		reject(source, pin, "raw");
		int rejected=3;
		for(String bad : new String[]{content.replace("#composition_dims\t400", "#composition_dims\t399"),
			content.replace("0.12345678912345678", "NaN"), content.replace("0.12345678912345678", "-0.1"),
			content.replace("0.8765432108765432", "0.1"), content+vector,
			HEADER, content.replace("#composition_rawpseudo\t0", "#composition_rawpseudo\t1"),
			content+"#background_400f\t0\textra\n", content+vector.toString().replace("#background_400f", "#background_other")}){
			final Path malformed=root.resolve("bad"+rejected+".tsv"); write(malformed, bad);
			reject(malformed, DigestSuffix.file(malformed.toString()), "centered"); rejected++;
		}
		System.err.println("SIDECAR_COMPOSITION_BACKGROUND_PASS default_bytes=true exact_header=true closed_form=true source_unchanged=true rejections="+rejected);
	}

	private static void reject(final Path source, final String pin, final String variant){
		try{SidecarCompositionBackground.load(source.toString(), pin, variant, 0, "families");}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Expected pinned-background rejection: "+source);
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}
	private static void check(final boolean ok, final String reason){if(!ok){throw new AssertionError(reason);}}
	private static final String HEADER="#schema_version\t1\n#kind\tfamily_shortlist_sidecar_v1\n"
		+"#composition_alphabet\taa20\n#composition_k\t2\n#composition_dims\t400\n"
		+"#composition_variant\tcentered\n#composition_rawpseudo\t0\n#composition_background\tfamilies\n";
}
