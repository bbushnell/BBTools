package tax;

import java.util.Arrays;
import java.util.List;

import parse.LineParserS1;

/**
 * Mutable, rank-indexed view of nodes resolved through TaxTree's shared instance.
 * TaxID construction walks canonical ancestors; string parsing looks up names
 * independently and retains the last uniquely resolved field at each normal rank.
 * The shared tree must already be loaded; string lookup builds name indexes lazily.
 * Node references are borrowed, and this view is not safe for concurrent mutation.
 * @author Brian Bushnell
 */
public class Lineage{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a Lineage by parsing a taxonomic lineage string.
	 * @param line_ Semicolon-delimited lineage string in format like
	 * "k__Bacteria;p__Actinomycetota;c__Actinomycetes;..."
	 */
	public Lineage(String line_){
		parse(line_);
	}
	
	/** Constructs a Lineage from a taxonomy ID by traversing the taxonomic tree.
	 * @param tid_ NCBI taxonomy ID to build lineage from */
	public Lineage(int tid_){
		set(tid_);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Population          ----------------*/
	/*--------------------------------------------------------------*/

	/** Walks a resolved ID to the root and records canonical nodes by normal rank.
	 * Clears the node array and original text first. If several canonical ancestors
	 * have the same normal rank, the later, more ancestral node replaces the earlier
	 * node. An unresolved ID leaves the node array empty but is retained in tid.
	 * @param tid_ TaxID resolved by the shared tree's getNode method */
	public void set(int tid_){
		clear();
		tid=tid_;
		TaxTree tree=TaxTree.getTree();
		TaxNode tn=tree.getNode(tid);
		if(tn!=null){level=tn.level;}
		for(; tn!=null; tn=tree.getNode(tn.pid)){
			if(tn.canonical()){nodes[tn.level]=tn;}
			if(tn.id==tn.pid){break;}
		}
	}
	
	//sk__Bacteria;k__Bacillati;p__Actinomycetota;c__Actinomycetes;o__Micromonosporales;
	//f__Micromonosporaceae;g__Salinispora;s__Salinispora arenicola
	/**
	 * Parses a semicolon-delimited taxonomic lineage string.
	 * Expected format: "k__Bacteria;p__Actinomycetota;c__Actinomycetes;..."
	 * where each field has format "level__name". Both short rank prefixes and
	 * full rank names are accepted. Names with no match or more than
	 * one match are ignored. The actual resolved node determines the array rank;
	 * fields are not checked for mutual ancestry or sorted into rank order.
	 *
	 * @param line_ Semicolon-delimited lineage string to parse
	 * @return TaxID of the last uniquely resolved field, or -1 if none resolves
	 */
	public int parse(String line_){
		clear();
		line=line_;
		TaxTree tree=TaxTree.getTree();

		LineParserS1 lp=new LineParserS1(';');
		LineParserS1 lp_=new LineParserS1('_');
		lp.set(line);
		for(int i=0; i<lp.terms(); i++){
			String field=lp.parseString(i);
			lp_.set(field);
			String ss=lp_.parseString(0);
			//Validate the declared rank; the resolved node still determines its slot.
			int shortlevel=(ss.length()<=2 ? TaxTree.shortStringToLevel(ss) : TaxTree.stringToLevel(ss));
			String name=field.substring(ss.length()+2);
			if(name!=null && name.length()>0){
				List<TaxNode> list=tree.getNodesByName(name);
				if(list==null || list.size()>1){
					//ignore it
				}else{
					TaxNode tn=list.get(0);
					nodes[tn.level]=tn;
					tid=tn.id;
					level=tn.level;
				}
			}
		}
		return tid;
	}
	
	/** Clears original text and node references, and resets tid and level to -1. */
	public void clear(){
		line=null;
		Arrays.fill(nodes, null);
		tid=-1;
		level=-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Original parsed text, or null after clear or TaxID-based population. */
	public String line;
	/** Requested TaxID or last uniquely resolved string field, initially -1. */
	int tid=-1;
	/** Normal rank of the last resolved node, or -1 after clearing/unresolved input. */
	int level=-1;
	/** Borrowed nodes indexed by normal rank, from NO_RANK through LIFE. */
	public final TaxNode[] nodes=new TaxNode[TaxTree.LIFE+1];
	
}
