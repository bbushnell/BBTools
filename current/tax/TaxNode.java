package tax;

import java.io.Serializable;
import java.util.Comparator;

import shared.Tools;

/**
 * Represents a taxonomic identifier, such as a specific genus.
 * Includes the name, NCBI numeric id, parent id, and taxonomic level.
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxNode implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -4618526038942239246L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a TaxNode with ID and name, using default values for parent and level.
	 * @param id_ NCBI taxonomic identifier
	 * @param name_ Scientific or common name for this taxonomic unit
	 */
	public TaxNode(int id_, String name_){
		this(id_, -1, -1, -1, name_);
	}
	
	/**
	 * Creates a TaxNode with complete taxonomic information.
	 *
	 * @param id_ NCBI taxonomic identifier
	 * @param parent_ Parent node's taxonomic ID
	 * @param level_ Standard taxonomic level (kingdom, phylum, etc.)
	 * @param levelExtended_ Extended level including intermediate ranks
	 * @param name_ Scientific or common name for this taxonomic unit
	 */
	public TaxNode(int id_, int parent_, int level_, int levelExtended_, String name_){
		id=id_;
		pid=parent_;
		level=level_;
		levelExtended=levelExtended_;
		setOriginalLevel(levelExtended);
		name=name_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * @param split
	 * @param idx
	 * @return True if the node's name matches the 
	 */
	public boolean matchesName(String[] split, int idx, TaxTree tree) {
		if(idx<0){return true;}
		if(!split[idx].equalsIgnoreCase(name)){return false;}
		return tree.getNode(pid).matchesName(split, idx-1, tree);
	}
	
	@Override
	public String toString(){
		return "("+id+","+pid+","+countRaw+","+countSum+",'"+levelStringExtended(false)+"',"+(canonical() ? "T" : "F")+",'"+name+"')";
	}
	
	/**
	 * Compares this TaxNode with another for equality based on ID, parent,
	 * level, flag, and name fields.
	 * @param b TaxNode to compare against
	 * @return true if all key fields match between the two nodes
	 */
	public boolean equals(TaxNode b){
		if(id!=b.id || pid!=b.pid || levelExtended!=b.levelExtended || flag!=b.flag){return false;}
		if(name==b.name){return true;}
		if((name==null) != (b.name==null)){return false;}
		return name.equals(b.name);
	}
	
	/**
	 * Increments the raw count for this taxonomic node.
	 * Raw counts represent direct assignments to this specific taxon.
	 * @param amt Amount to add to the raw count
	 * @return Updated raw count value
	 */
	public long incrementRaw(long amt){
		if(amt==0){return countRaw;}
		if(verbose){System.err.println("incrementRaw("+amt+") node: "+this);}
		countRaw+=amt;
		assert(countRaw>=0) : "Overflow! "+countRaw+", "+amt;
		return countRaw;
	}
	
	/**
	 * Increments the sum count for this taxonomic node.
	 * Sum counts include assignments to this node and all descendants.
	 * @param amt Amount to add to the sum count (can be negative)
	 * @return Updated sum count value
	 */
	public long incrementSum(long amt){
		if(amt==0){return countSum;}
		if(verbose){System.err.println("incrementSum("+amt+") node: "+this);}
		countSum+=amt;
		assert(countSum>=0 || amt<0) : "Overflow! "+countSum+", "+amt;
		return countSum;
	}
	
	/** Checks if this node represents a simple taxonomic level.
	 * @return true if the extended level corresponds to a basic taxonomic rank */
	public boolean isSimple(){
		return TaxTree.isSimple(levelExtended);
	}
	
	/** Checks if this node represents a simple taxonomic level using alternative criteria.
	 * @return true if the extended level meets the simple2 classification */
	public boolean isSimple2(){
		return TaxTree.isSimple2(levelExtended);
	}
	
//	public String levelString(){return level<0 ? "unknown" : TaxTree.levelToString(level);}
	
	/**
	 * Returns the string representation of the taxonomic level.
	 * @param original If true, uses original level; if false, uses current extended level
	 * @return Human-readable taxonomic level name or "unknown"
	 */
	public String levelStringExtended(boolean original){
		int x=(original ? originalLevel() : levelExtended);
		return x<0 ? "unknown" : TaxTree.levelToStringExtended(x);
	}

	/** Returns abbreviated string representation of the taxonomic level.
	 * @return Single character or short abbreviation for the taxonomic level */
	public String levelToStringShort() {return level<0 ? "x" : TaxTree.levelToStringShort(level);}
	

	
	/** Checks if this node represents an unclassified taxonomic group.
	 * @return true if the name starts with "unclassified" */
	public boolean isUnclassified(){
		return name.startsWith("unclassified");
	}
	
	/** Checks if this node represents an environmental sample.
	 * @return true if the name starts with "environmental" */
	public boolean isEnvironmentalSample(){
		return name.startsWith("environmental");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Comparator for sorting TaxNode objects by abundance counts.
	 * Sorts primarily by sum count (descending), then by taxonomic level and ID. */
	public static class CountComparator implements Comparator<TaxNode>{
		
		@Override
		public int compare(TaxNode a, TaxNode b) {
			long x=b.countSum-a.countSum;
//			System.err.println("x="+x+" -> "+Tools.longToInt(x));
			if(x!=0){return Tools.longToInt(x);}
			return a.levelExtended==b.levelExtended ? a.id-b.id : a.levelExtended-b.levelExtended;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int hashCode(){return id;}
	
	/*--------------------------------------------------------------*/
	
	/** Checks if this node is marked as canonical in the taxonomic hierarchy.
	 * @return true if the canonical flag bit is set */
	public boolean canonical(){
		return (flag&CANON_MASK)==CANON_MASK;
	}
	
	/** Checks if the taxonomic level has been modified from its original value.
	 * @return true if current level differs from original level */
	public boolean levelChanged(){
		return originalLevel()!=levelExtended;
	}
	
	/** Retrieves the original taxonomic level before any modifications.
	 * @return Original taxonomic level, or -1 if not set */
	public int originalLevel(){
		int x=(int)(flag&ORIGINAL_LEVEL_MASK);
		return x==ORIGINAL_LEVEL_MASK ? -1 : x;
	}
	
	/** Checks if this node represents the cellular organisms taxonomic group.
	 * @return true if this node's ID matches the cellular organisms constant */
	public boolean cellularOrganisms(){
		return id==TaxTree.CELLULAR_ORGANISMS_ID;
	}
	
//	public int numChildren(){
//		return numChildren;
//	}
//
//	public int minParentLevelExtended(){
//		return minParentLevelExtended;
//	}
//
//	public int maxChildLevelExtended(){
//		return maxChildLevelExtended;
//	}
	
	/**
	 * Gets the minimum taxonomic level among this node and its ancestors.
	 * Uses parent level information if this node has no valid level.
	 * @return Minimum ancestor level including this node
	 */
	int minAncestorLevelIncludingSelf(){
		return levelExtended<1 ? minParentLevelExtended : levelExtended;
	}
	
	/**
	 * Gets the maximum taxonomic level among this node and its descendants.
	 * Uses child level information if this node has no valid level.
	 * @return Maximum descendant level including this node
	 */
	int maxDescendantLevelIncludingSelf(){
		return levelExtended<1 ? maxChildLevelExtended : levelExtended;
	}
	
	/**
	 * Creates a simplified version of the taxonomic name with only alphanumeric
	 * characters and underscores, suitable for file naming or identifiers.
	 * @return Simplified name with spaces replaced by underscores, or null if name is null
	 */
	public String simpleName(){
		if(name==null){return null;}
		StringBuilder sb=new StringBuilder();
		char last='?';
		for(int i=0; i<name.length(); i++){
			char c=name.charAt(i);
			if((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='1' && c<='0')){
				sb.append(c);
				last=c;
			}else{
				if(sb.length()>0 && last!=' '){sb.append(' ');}
				last=' ';
			}
		}
		String s=sb.toString().trim();
		return s.replace(' ', '_');
	}
	
	/** Checks if this node has a defined taxonomic rank.
	 * @return true if the extended level is not the "no rank" designation */
	public boolean isRanked() {return levelExtended!=TaxTree.NO_RANK_E;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Setters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sets the canonical status flag for this taxonomic node.
	 * @param b true to mark as canonical, false otherwise */
	public void setCanonical(boolean b){
		if(b){flag=flag|CANON_MASK;}
		else{flag=flag&~CANON_MASK;}
	}
	
	/** Sets the original taxonomic level in the flag bits.
	 * @param x Original taxonomic level to store */
	public void setOriginalLevel(int x){
		flag=(flag&~ORIGINAL_LEVEL_MASK)|(x&ORIGINAL_LEVEL_MASK);
	}
	
	/** Return true if changed */
	boolean discussWithParent(TaxNode parent){
		final int oldChildLevel=parent.maxChildLevelExtended;
		final int oldParentLevel=minParentLevelExtended;
		parent.maxChildLevelExtended=Tools.max(parent.maxChildLevelExtended, maxDescendantLevelIncludingSelf());
		minParentLevelExtended=Tools.min(parent.minAncestorLevelIncludingSelf(), minParentLevelExtended);
		return oldChildLevel!=parent.maxChildLevelExtended || oldParentLevel!=minParentLevelExtended;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** NCBI taxonomic identifier for this node */
	public final int id;
	/** Scientific or common name for this taxonomic unit */
	public final String name;
	/** Parent node's taxonomic identifier */
	public int pid;
	/** Standard taxonomic level (kingdom, phylum, class, etc.) */
	public int level;
	/** Extended taxonomic level including intermediate and non-standard ranks */
	public int levelExtended;
	
	/** Number of direct child nodes in the taxonomic tree */
	public int numChildren=0;
	/** Minimum extended level among all ancestor nodes */
	public int minParentLevelExtended=TaxTree.LIFE_E;
	/** Maximum extended level among all descendant nodes */
	public int maxChildLevelExtended=TaxTree.NO_RANK_E;
	
	private long flag=0;
	
	/** Direct count of sequences assigned specifically to this taxonomic node */
	public long countRaw=0;
	/** Total count including sequences assigned to this node and all descendants */
	public long countSum=0;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final long ORIGINAL_LEVEL_MASK=63; //bits 0-5
	private static final long CANON_MASK=64; //bit 6
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Controls verbose debug output during count operations */
	public static final boolean verbose=false;
	/** Comparator for sorting TaxNodes by count values in descending order */
	public static final CountComparator countComparator=new CountComparator();
	
	
}
