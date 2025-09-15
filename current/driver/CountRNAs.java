package driver;

import dna.Data;
import dna.Gene;

public class CountRNAs {
	
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
