package tax;

import java.util.Arrays;
import java.util.List;

import shared.LineParserS1;

public class Lineage {
	
	public Lineage(String line_) {
		parse(line_);
	}
	
	public Lineage(int tid_) {
		set(tid_);
	}
	
	public void set(int tid_) {
		clear();
		tid=tid_;
		TaxTree tree=TaxTree.getTree();
		TaxNode tn=tree.getNode(tid);
		if(tn!=null) {level=tn.level;}
//		System.err.println("Looking for "+tid+"; found "+tn);
		for(; tn!=null; tn=tree.getNode(tn.pid)) {
//			System.err.println("Assigning "+tn+": canon="+tn.canonical()+", tid==pid="+(tn.id==tn.pid));
			if(tn.canonical()) {nodes[tn.level]=tn;}
			if(tn.id==tn.pid) {break;}
		}
	}
	
	//sk__Bacteria;k__Bacillati;p__Actinomycetota;c__Actinomycetes;o__Micromonosporales;
	//f__Micromonosporaceae;g__Salinispora;s__Salinispora arenicola
	public int parse(String line_) {
		clear();
		line=line_;
		TaxTree tree=TaxTree.getTree();

		LineParserS1 lp=new LineParserS1(';');
		LineParserS1 lp_=new LineParserS1('_');
		lp.set(line);
		for(int i=0; i<lp.terms(); i++) {
			String field=lp.parseString(i);
			lp_.set(field);
			String ss=lp.parseString(0);
			int shortlevel=TaxTree.stringToLevel(ss);
			String name=field.substring(ss.length()+2);
			if(name!=null && name.length()>0) {
				List<TaxNode> list=tree.getNodesByName(name);
				if(list==null || list.size()>1) {
					//ignore it
				}else {
					TaxNode tn=list.get(0);
					nodes[tn.level]=tn;
					tid=tn.id;
					level=tn.level;
				}
			}
		}
		return tid;
	}
	
	public void clear() {
		line=null;
		Arrays.fill(nodes, null);
		tid=-1;
	}
	
	public String line;
	int tid=-1;
	int level=-1;
	public final TaxNode[] nodes=new TaxNode[TaxTree.LIFE+1];
	
}
