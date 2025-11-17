package driver;

import dna.Data;
import dna.Gene;

/**
 * Counts and categorizes genes across human chromosomes by type.
 * Processes gene data for chromosomes 1-24, classifying genes as pseudogenes,
 * untranslated, or coding, and outputs detailed gene type statistics.
 * @author Brian Bushnell
 */
public class CountRNAs {
	
	/**
	 * Main entry point that counts and categorizes genes by type across chromosomes.
	 * Iterates through chromosomes 1-24, examining each gene to classify it as
	 * pseudogene, untranslated (non-coding), or coding gene. Outputs summary
	 * statistics for each category.
	 *
	 * @param args Command line arguments: [0] genome build number, [1] gene map file path
	 */
	public static void main(String[] args){
		Data.GENOME_BUILD=Integer.parseInt(args[0]);
		Data.GENE_MAP=args[1];
		long coding=0;
		long noncoding=0;
		long pseudo=0;
		for(byte chrom=1; chrom<=24; chrom++){
			Gene[] genes=Data.getGenes(chrom);
			for(Gene g : genes){
				if(g.pseudo){
					pseudo++;
				}else if(g.untranslated){
					noncoding++;
				}else{
					coding++;
				}
			}
		}
		System.out.println("Gene map: "+Data.GENE_MAP);
		System.out.println("Pseudogenes: "+pseudo);
		System.out.println("Translated Genes: "+coding);
		System.out.println("Untranslated Genes: "+noncoding);
	}
	
}
