package driver;

import java.io.File;

import tax.TaxNode;
import tax.TaxTree;

/**
 * Utility for standardizing RefSeq genome file naming conventions.
 * Renames RefSeq genome files by prefixing with "refseq_" and taxonomy ID.
 * @author Brian Bushnell
 */
public class RenameRefseqFiles {
	
	/**
	 * Renames RefSeq genome files by adding "refseq_" prefix to existing files.
	 * Loads taxonomy tree from default file, iterates through all taxonomy nodes,
	 * and renames files from "[id].fa.gz" to "refseq_[id].fa.gz" format.
	 * @param args Command-line arguments where args[0] is the base directory path
	 */
	public static void main(String[] args){
		TaxTree tree=TaxTree.loadTaxTree(TaxTree.defaultTreeFile(), System.err, false, false);
		for(TaxNode tn : tree.nodes){
			if(tn!=null){
				String dir=tree.toDir(tn, args[0]);
				String path=dir+tn.id+".fa.gz";
				File f=new File(path);
				if(f.exists()){
					//NOTE [driver/RenameRefseqFiles#001] LOW/dev (no .sh, no callers): File.renameTo() return ignored →
					//silent rename failure (target exists / cross-fs / perm), same as RenameByHeader#001. args[0] also
					//unguarded (AIOOBE 0 args). Dead one-off → LOW.
					f.renameTo(new File(dir+"refseq_"+tn.id+".fa.gz"));
				}
			}
		}
	}
	
}
