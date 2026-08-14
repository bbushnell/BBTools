package prok;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import tax.TaxNode;
import tax.TaxTree;

/**
 * Maps taxIDs to phylum names using the BBTools taxonomy tree.
 * Reads taxIDs from stdin (one per line), writes taxID\tphylum to stdout.
 * @author Neptune, Brian Bushnell
 */
public class TaxIDToPhylum {

	public static void main(String[] args) throws Exception{
		TaxTree tree=TaxTree.loadTaxTree(System.err, false, false);
		BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
		String line;
		while((line=br.readLine())!=null){
			line=line.trim();
			if(line.isEmpty()){continue;}
			try{
				int tid=Integer.parseInt(line);
				TaxNode node=tree.getNode(tid, true);//skipAssertion: tids newer than the tree dump yield null, not a crash
				if(node==null){continue;}
				int phylumTid=tree.getIdAtLevelExtended(tid, TaxTree.PHYLUM_E);
				TaxNode phylumNode=(phylumTid>0 ? tree.getNode(phylumTid) : null);
				String phylum=(phylumNode!=null ? phylumNode.name : "Unknown");
				System.out.println(tid+"\t"+phylum);
			}catch(NumberFormatException e){}
		}
	}
}
