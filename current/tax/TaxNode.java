package tax;

import java.io.Serializable;
import java.util.Comparator;

import shared.Tools;

/**
 * Represents a taxonomic identifier with hierarchical relationships.
 * Includes NCBI numeric ID, parent relationships, taxonomic level information,
 * and counts for taxonomic abundance analysis. Parent links are TaxIDs resolved
 * through the owning tree, not necessarily array indices. Mutable fields and
 * counters are unsynchronized; callers must coordinate shared updates.
 *
 * @author Brian Bushnell
 * @date Mar 6, 2015
 */
public class TaxNode implements Serializable{

	/** Compatibility identifier for serialized taxonomy resources. */
	private static final long serialVersionUID=-4618526038942239246L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a named node whose parent and ranks are initially unknown.
	 * @param id_ Durable TaxID
	 * @param name_ Scientific name, retained by reference */
	public TaxNode(int id_, String name_){
		this(id_, -1, -1, -1, name_);
	}
	
	/** Creates a node with zero counters and an initially unset canonical flag.
	 * The original rank is saved in six bits; -1 denotes an unknown rank.
	 * @param id_ Durable TaxID
	 * @param parent_ Parent TaxID; a root names itself
	 * @param level_ Normal rank or -1 while constructing an incomplete node
	 * @param levelExtended_ Extended rank or -1 while constructing an incomplete node
	 * @param name_ Scientific name, retained by reference */
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

	/** Matches successive node names against tokens, walking backward through the lineage.
	 * Matching is case-insensitive. Every parent needed by the tokens must resolve.
	 * @param split Ancestor-to-descendant name tokens
	 * @param idx Token for this node; a negative index ends the match successfully
	 * @param tree Owning tree used to resolve parents
	 * @return Whether all tokens through idx match this node and its parents */
	public boolean matchesName(String[] split, int idx, TaxTree tree){
		if(idx<0){return true;}
		if(!split[idx].equalsIgnoreCase(name)){return false;}
		return tree.getNode(pid).matchesName(split, idx-1, tree);
	}
	
	/** Renders IDs, counters, current extended rank, canonical flag and name.
	 * @return Diagnostic tuple; the name is not escaped */
	@Override
	public String toString(){
		return "("+id+","+pid+","+countRaw+","+countSum+",'"+levelStringExtended(false)+"',"+(canonical() ? "T" : "F")+",'"+name+"')";
	}
	
	/** Compares IDs, parent, extended rank, flags and name, excluding counts.
	 * This overload does not override Object.equals; generic collections still use
	 * object identity. The argument must be non-null.
	 * @param b Node to compare
	 * @return Whether the selected structural fields match */
	public boolean equals(TaxNode b){
		if(id!=b.id || pid!=b.pid || levelExtended!=b.levelExtended || flag!=b.flag){return false;}
		if(name==b.name){return true;}
		if((name==null) != (b.name==null)){return false;}
		return name.equals(b.name);
	}
	
	/** Adds to the direct abundance count, asserting that the result is nonnegative.
	 * @param amt Signed adjustment; arithmetic uses ordinary long addition
	 * @return Updated direct count */
	public long incrementRaw(long amt){
		if(amt==0){return countRaw;}
		if(verbose){System.err.println("incrementRaw("+amt+") node: "+this);}
		countRaw+=amt;
		assert(countRaw>=0) : "Overflow! "+countRaw+", "+amt;
		return countRaw;
	}
	
	/** Adds to the accumulated count. Negative adjustments may produce negative totals.
	 * @param amt Signed adjustment, including subtraction during node gathering
	 * @return Updated accumulated count */
	public long incrementSum(long amt){
		if(amt==0){return countSum;}
		if(verbose){System.err.println("incrementSum("+amt+") node: "+this);}
		countSum+=amt;
		assert(countSum>=0 || amt<0) : "Overflow! "+countSum+", "+amt;
		return countSum;
	}
	
	/** Tests whether the current extended rank is a nonzero standard rank.
	 * @return Result of TaxTree.isSimple for this rank */
	public boolean isSimple(){
		return TaxTree.isSimple(levelExtended);
	}
	
	/** Tests standard ranks plus the strain/substrain ranks accepted by isSimple2.
	 * @return Whether the current rank belongs to that expanded set */
	public boolean isSimple2(){
		return TaxTree.isSimple2(levelExtended);
	}
	
//	public String levelString(){return level<0 ? "unknown" : TaxTree.levelToString(level);}
	
	/** Resolves the selected extended rank to its name.
	 * @param original Use the saved original rank instead of the current rank
	 * @return Rank name, or unknown for a negative rank */
	public String levelStringExtended(boolean original){
		int x=(original ? originalLevel() : levelExtended);
		return x<0 ? "unknown" : TaxTree.levelToStringExtended(x);
	}

	/** Returns the normal-rank abbreviation, or x for a negative rank.
	 * @return Short rank label */
	public String levelToStringShort(){return level<0 ? "x" : TaxTree.levelToStringShort(level);}
	

	
	/** Tests the case-sensitive unclassified prefix; the name must be non-null.
	 * @return Whether the name starts with unclassified */
	public boolean isUnclassified(){
		return name.startsWith("unclassified");
	}
	
	/** Tests the case-sensitive environmental prefix; the name must be non-null.
	 * @return Whether the name starts with environmental */
	public boolean isEnvironmentalSample(){
		return name.startsWith("environmental");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Orders nodes by descending cumulative count, then extended rank and TaxID.
	 * Count subtraction assumes its difference is representable as a long. */
	public static class CountComparator implements Comparator<TaxNode>{

		/** Creates a stateless count comparator. */
		public CountComparator(){}
		
		/**
		 * Compares two TaxNode objects for sorting by abundance.
		 * Primary sort: sum count (descending)
		 * Secondary sort: taxonomic level (ascending)
		 * Tertiary sort: taxonomic ID (ascending)
		 *
		 * @param a First TaxNode to compare
		 * @param b Second TaxNode to compare
		 * @return Negative if a sorts first, positive if b sorts first, zero on a tie
		 */
		@Override
		public int compare(TaxNode a, TaxNode b){
			long x=b.countSum-a.countSum;
//			System.err.println("x="+x+" -> "+Tools.longToInt(x));
			if(x!=0){return Tools.longToInt(x);}
			return a.levelExtended==b.levelExtended ? a.id-b.id : a.levelExtended-b.levelExtended;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns hash code based on the taxonomic ID.
	 * @return The NCBI taxonomic ID as hash code */
	@Override
	public final int hashCode(){return id;}
	
	/*--------------------------------------------------------------*/
	
	/** Reads the saved canonical flag, which is distinct from testing the current rank.
	 * @return Whether the canonical bit is set */
	public boolean canonical(){
		return (flag&CANON_MASK)==CANON_MASK;
	}
	
	/** Compares the saved original extended rank with the current rank.
	 * @return Whether the rank differs */
	public boolean levelChanged(){
		return originalLevel()!=levelExtended;
	}
	
	/** Decodes the original rank from its six-bit field.
	 * @return Saved extended rank, or -1 for the all-ones sentinel */
	public int originalLevel(){
		int x=(int)(flag&ORIGINAL_LEVEL_MASK);
		return x==ORIGINAL_LEVEL_MASK ? -1 : x;
	}

	/** Supplies the complete persisted flag word to TaxTreeText and overlays.
	 * @return Unmodified flag bits */
	final long rawFlag(){return flag;}
	
	/** Tests the durable ID for the cellular-organisms grouping node.
	 * @return Whether this is TaxTree.CELLULAR_ORGANISMS_ID */
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
	
	/** Returns this ranked node's rank, or the inferred ancestor bound if unranked.
	 * @return Extended-rank lower bound including this node */
	int minAncestorLevelIncludingSelf(){
		return levelExtended<1 ? minParentLevelExtended : levelExtended;
	}
	
	/** Returns this ranked node's rank, or the inferred descendant bound if unranked.
	 * @return Extended-rank upper bound including this node */
	int maxDescendantLevelIncludingSelf(){
		return levelExtended<1 ? maxChildLevelExtended : levelExtended;
	}
	
	/** Replaces runs of non-ASCII-alphanumeric characters with one underscore.
	 * Leading and trailing separator runs are removed; ASCII letter case is retained.
	 * @return Simplified name, null for a null name, or empty if no characters survive */
	public String simpleName(){
		if(name==null){return null;}
		StringBuilder sb=new StringBuilder();
		char last='?';
		for(int i=0; i<name.length(); i++){
			char c=name.charAt(i);
			if((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9')){
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
	
	/** Tests only inequality with no-rank; an incomplete negative rank also passes.
	 * @return Whether levelExtended differs from NO_RANK_E */
	public boolean isRanked(){return levelExtended!=TaxTree.NO_RANK_E;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Setters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sets or clears the canonical bit without altering the saved original rank.
	 * @param b New canonical status */
	public void setCanonical(boolean b){
		if(b){flag=flag|CANON_MASK;}
		else{flag=flag&~CANON_MASK;}
	}
	
	/** Stores the low six bits as original rank, preserving other flag bits.
	 * @param x Extended rank; -1 encodes the unknown sentinel */
	public void setOriginalLevel(int x){
		flag=(flag&~ORIGINAL_LEVEL_MASK)|(x&ORIGINAL_LEVEL_MASK);
	}

	/** Restores all persisted flags without interpreting or validating them.
	 * @param flag_ Complete flag word */
	final void setRawFlag(long flag_){flag=flag_;}
	
	/** Propagates descendant bounds upward and ancestor bounds downward for one edge.
	 * @param parent Resolved parent node; both nodes must be exclusively mutable
	 * @return Whether either propagated bound changed */
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
	
	/** Durable external TaxID, independent of a packed tree's array position. */
	public final int id;
	/** Scientific name, retained without normalization. */
	public final String name;
	/** Parent's durable TaxID; a root has pid equal to id. */
	public int pid;
	/** Current normal rank. */
	public int level;
	/** Current extended rank. */
	public int levelExtended;
	
	/** Direct child count, populated by tree construction. */
	public int numChildren=0;
	/** Lowest known ranked ancestor bound propagated during tree construction. */
	public int minParentLevelExtended=TaxTree.LIFE_E;
	/** Highest known ranked descendant bound propagated during tree construction. */
	public int maxChildLevelExtended=TaxTree.NO_RANK_E;
	
	/** Packed original-rank and canonical-status flags. */
	private long flag=0;
	
	/** Direct abundance before propagation to ancestors. */
	public long countRaw=0;
	/** Accumulated abundance, also mutated by destructive node gathering. */
	public long countSum=0;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Bits zero through five store the original rank; all ones mean unknown. */
	private static final long ORIGINAL_LEVEL_MASK=63;
	/** Bit six marks the node's saved canonical status. */
	private static final long CANON_MASK=64;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Compile-time diagnostic logging switch. */
	public static final boolean verbose=false;
	/** Shared stateless abundance-order comparator. */
	public static final CountComparator countComparator=new CountComparator();
	
	
}
