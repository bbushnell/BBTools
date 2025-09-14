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
	
	public TaxNode(int id_, String name_){
		this(id_, -1, -1, -1, name_);
	}
	
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
	
	public boolean equals(TaxNode b){
		if(id!=b.id || pid!=b.pid || levelExtended!=b.levelExtended || flag!=b.flag){return false;}
		if(name==b.name){return true;}
		if((name==null) != (b.name==null)){return false;}
		return name.equals(b.name);
	}
	
	public long incrementRaw(long amt){
		if(amt==0){return countRaw;}
		if(verbose){System.err.println("incrementRaw("+amt+") node: "+this);}
		countRaw+=amt;
		assert(countRaw>=0) : "Overflow! "+countRaw+", "+amt;
		return countRaw;
	}
	
	public long incrementSum(long amt){
		if(amt==0){return countSum;}
		if(verbose){System.err.println("incrementSum("+amt+") node: "+this);}
		countSum+=amt;
		assert(countSum>=0 || amt<0) : "Overflow! "+countSum+", "+amt;
		return countSum;
	}
	
	public boolean isSimple(){
		return TaxTree.isSimple(levelExtended);
	}
	
	public boolean isSimple2(){
		return TaxTree.isSimple2(levelExtended);
	}
	
//	public String levelString(){return level<0 ? "unknown" : TaxTree.levelToString(level);}
	
	public String levelStringExtended(boolean original){
		int x=(original ? originalLevel() : levelExtended);
		return x<0 ? "unknown" : TaxTree.levelToStringExtended(x);
	}

	public String levelToStringShort() {return level<0 ? "x" : TaxTree.levelToStringShort(level);}
	

	
	public boolean isUnclassified(){
		return name.startsWith("unclassified");
	}
	
	public boolean isEnvironmentalSample(){
		return name.startsWith("environmental");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public boolean canonical(){
		return (flag&CANON_MASK)==CANON_MASK;
	}
	
	public boolean levelChanged(){
		return originalLevel()!=levelExtended;
	}
	
	public int originalLevel(){
		int x=(int)(flag&ORIGINAL_LEVEL_MASK);
		return x==ORIGINAL_LEVEL_MASK ? -1 : x;
	}
	
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
	
	int minAncestorLevelIncludingSelf(){
		return levelExtended<1 ? minParentLevelExtended : levelExtended;
	}
	
	int maxDescendantLevelIncludingSelf(){
		return levelExtended<1 ? maxChildLevelExtended : levelExtended;
	}
	
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
	
	public boolean isRanked() {return levelExtended!=TaxTree.NO_RANK_E;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Setters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public void setCanonical(boolean b){
		if(b){flag=flag|CANON_MASK;}
		else{flag=flag&~CANON_MASK;}
	}
	
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
	
	public final int id;
	public final String name;
	public int pid;
	public int level;
	public int levelExtended;
	
	public int numChildren=0;
	public int minParentLevelExtended=TaxTree.LIFE_E;
	public int maxChildLevelExtended=TaxTree.NO_RANK_E;
	
	private long flag=0;
	
	public long countRaw=0;
	public long countSum=0;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final long ORIGINAL_LEVEL_MASK=63; //bits 0-5
	private static final long CANON_MASK=64; //bit 6
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final boolean verbose=false;
	public static final CountComparator countComparator=new CountComparator();
	
	
}
