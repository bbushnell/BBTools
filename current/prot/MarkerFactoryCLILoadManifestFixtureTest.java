package prot;

import java.util.List;

/**
 * Disposable differential fixture for MarkerFactoryCLI.loadManifest's split-to-
 * LineParser1 conversion (RULE #12 baseline). Exercises comment/blank/whitespace-only
 * lines and leading/trailing whitespace within fields (incl. an empty lineage column)
 * to confirm the .trim() semantics survive the LineParser1 swap.
 */
public class MarkerFactoryCLILoadManifestFixtureTest {

	public static void main(String[] args){
		final List<GenomeProteins> genomes=MarkerFactoryCLI.loadManifest(args[0]);
		System.out.println("genomes="+genomes.size());
		for(GenomeProteins g : genomes){
			System.out.println("id=["+g.id+"] domain=["+g.domain+"] lineage=["+g.lineage+
				"] proteins="+g.size());
			for(ProteinSequence p : g.proteins){
				System.out.println("  protein id=["+p.id+"] len="+p.length());
			}
		}
	}
}
