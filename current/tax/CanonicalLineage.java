package tax;

import structures.ByteBuilder;

/** Seven-rank lineage view containing borrowed nodes from one taxonomy tree.
 * Traversal keeps the nearest ancestor at each exact extended rank, from kingdom
 * through species. Domain and intermediate ranks are omitted. This is a snapshot
 * of node references, not a copy of the nodes or their mutable fields.
 * @author Brian Bushnell */
public class CanonicalLineage{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Collects canonical-rank ancestors, including the starting node.
	 * Null trees and nonpositive IDs produce an entirely missing lineage. A positive
	 * ID is resolved through {@link TaxTree#getNode(int)}, including its assertion
	 * policy for IDs outside the tree. Parent chains must be acyclic except at a
	 * self-parent root; a missing parent ends traversal.
	 * @param tid Starting TaxID, retained even when no node resolves
	 * @param tree Tree supplying nodes, or null for an empty lineage */
	public CanonicalLineage(int tid, TaxTree tree){
		taxID=tid;
		if(tree==null || taxID<1){
			levelsDefined=0;
			return;
		}
		TaxNode current=tree.getNode(taxID);
		int levels=0;
		while(current!=null){
			int level=current.levelExtended;
			int index=-1;
			if(level==TaxTree.KINGDOM_E){
				index=0;
			}else if(level==TaxTree.PHYLUM_E){
				index=1;
			}else if(level==TaxTree.CLASS_E){
				index=2;
			}else if(level==TaxTree.ORDER_E){
				index=3;
			}else if(level==TaxTree.FAMILY_E){
				index=4;
			}else if(level==TaxTree.GENUS_E){
				index=5;
			}else if(level==TaxTree.SPECIES_E){
				index=6;
			}
			if(index>=0 && nodes[index]==null){
				nodes[index]=current;
				levels++;
			}
			if(current.pid==current.id){break;}
			current=tree.getNode(current.pid);
		}
		levelsDefined=levels;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Renders all seven ranks with k__, p__, c__, o__, f__, g__, s__ prefixes.
	 * Missing ranks use NA, except missing species uses tid_ followed by the original
	 * TaxID. Names are appended verbatim without escaping embedded delimiters.
	 * @return Newly allocated semicolon-delimited lineage, without a final semicolon */
	public ByteBuilder toBytes(){
		ByteBuilder bb=new ByteBuilder();
		for(int level=0; level<nodes.length; level++){
			TaxNode tn=nodes[level];
			if(level>0){bb.semi();}
			bb.append(prefixes[level]);
			bb.append(tn==null ? (level==6 ? "tid_"+taxID : "NA") : tn.name);
		}
		return bb;
	}

	/** Returns the same seven-rank representation as {@link #toBytes()}.
	 * @return Newly allocated lineage string */
	@Override
	public String toString(){return toBytes().toString();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Original requested TaxID, including an invalid or unresolved value. */
	public final int taxID;
	/** Number of the seven rank slots filled during construction. */
	public final int levelsDefined;
	/** Borrowed nodes ordered kingdom, phylum, class, order, family, genus, species. */
	final TaxNode[] nodes=new TaxNode[7];
	/** Rendering prefixes, aligned with the seven node slots. */
	private static final String[] prefixes={"k__", "p__", "c__", "o__", "f__", "g__", "s__"};
}
