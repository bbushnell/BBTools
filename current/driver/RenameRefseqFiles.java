package driver;

import java.io.File;

import tax.TaxNode;
import tax.TaxTree;

public class RenameRefseqFiles {
	
	public static void main(String[] args){
		TaxTree tree=TaxTree.loadTaxTree(TaxTree.defaultTreeFile(), System.err, false, false);
		for(TaxNode tn : tree.nodes){
			if(tn!=null){
				String dir=tree.toDir(tn, args[0]);
				String path=dir+tn.id+".fa.gz";
				File f=new File(path);
				if(f.exists()){
					f.renameTo(new File(dir+"refseq_"+tn.id+".fa.gz"));
				}
			}
		}
	}
	
}
