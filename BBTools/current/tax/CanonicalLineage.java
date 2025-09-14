package tax;

import structures.ByteBuilder;

public class CanonicalLineage {

	/** 
	 * Produces e.g. k__something;p__something, etc
	 * kpcofgs levels must all be present.
	 * Missing values should use NA, e.g. k__NA;
	 * 
	 * @param TaxID Taxonomy ID to get lineage for
	 * @param tree Taxonomy tree to use
	 */
	public CanonicalLineage(int tid, TaxTree tree) {
		taxID=tid;
		if(tree==null || taxID<1) {
			// If invalid input, nodes will remain null (rendered as NA)
			levelsDefined=0;
			return;
		}

		// Start from the given taxID and traverse up the tree
		TaxNode current = tree.getNode(taxID);
		int levels=0;
		while(current != null) {
			int level = current.levelExtended;
			int index = -1;

			// Map taxonomy levels to array indices
			if(level == TaxTree.KINGDOM_E) {
				index = 0; // kingdom
			} else if(level == TaxTree.PHYLUM_E) {
				index = 1; // phylum
			} else if(level == TaxTree.CLASS_E) {
				index = 2; // class
			} else if(level == TaxTree.ORDER_E) {
				index = 3; // order
			} else if(level == TaxTree.FAMILY_E) {
				index = 4; // family
			} else if(level == TaxTree.GENUS_E) {
				index = 5; // genus
			} else if(level == TaxTree.SPECIES_E) {
				index = 6; // species
			}

			// If this is one of our target levels and not already set, save it
			if(index >= 0 && nodes[index] == null) {
				nodes[index] = current;
				levels++;
			}

			// Move up to parent
			if(current.pid == current.id) {
				break; // Avoid infinite loop at root
			}
			current = tree.getNode(current.pid);
		}
		levelsDefined=levels;
	}

	/**
	 * Converts the lineage to a ByteBuilder with DADA2 format.
	 * @return ByteBuilder containing the formatted lineage
	 */
	public ByteBuilder toBytes() {
		ByteBuilder bb = new ByteBuilder();

		// Standard prefixes for DADA2 taxonomy levels

		for(int level=0; level<nodes.length; level++) {
			TaxNode tn=nodes[level];
			if(level > 0) {bb.semi();}
			bb.append(prefixes[level]);
			bb.append(tn==null ? (level==6 ? "tid_"+taxID : "NA") : tn.name);
		}

		return bb;
	}

	public String toString() {
		return toBytes().toString();
	}
	
	public final int taxID;
	public final int levelsDefined;
	/** Holds nodes for the 7 standard taxonomic levels */
	final TaxNode[] nodes=new TaxNode[7];
	private static final String[] prefixes = {"k__", "p__", "c__", "o__", "f__", "g__", "s__"};
}