package driver;

import dna.ChromosomeArray;
import dna.Data;

/**
 * Utility for comparing reference genome chromosomes stored as ChromosomeArray files.
 * Performs chromosome-by-chromosome comparison of two genome datasets.
 * Reports differences at the base level and index range mismatches.
 * @author Brian Bushnell
 */
public class CompareReferenceGenomes {
	
	/**
	 * Program entry point for comparing two reference genomes.
	 * Expects two arguments: file patterns with '#' placeholder for chromosome numbers.
	 * @param args Command-line arguments: [pattern1] [pattern2] where # is replaced by chromosome
	 */
	public static void main(String[] args){
		compareGenomes(args[0], args[1]);
	}
	
	/**
	 * Compares chromosomes 1-25 between two genome datasets using file patterns.
	 * Replaces '#' in patterns with chromosome numbers to generate filenames.
	 * Loads ChromosomeArray objects and performs base-by-base comparison.
	 *
	 * @param pattern1 File pattern for first genome (e.g., "genome1_chr#.dat")
	 * @param pattern2 File pattern for second genome (e.g., "genome2_chr#.dat")
	 */
	public static void compareGenomes(String pattern1, String pattern2){
		for(byte chrom=1; chrom<=25; chrom++){
			System.out.println("Comparing chromosome "+chrom);
			String fname1=pattern1.replace("#", ""+chrom);
			String fname2=pattern2.replace("#", ""+chrom);
			ChromosomeArray cha=ChromosomeArray.read(fname1);
			ChromosomeArray chb=ChromosomeArray.read(fname2);
			boolean result=compare(cha, chb);
			System.out.println("..."+(result ? "identical." : "different."));
		}
	}
	
	/**
	 * Performs detailed comparison between two ChromosomeArray objects.
	 * Checks index ranges for compatibility and compares bases within overlapping regions.
	 * Reports mismatches to stdout showing chromosome, position, and differing bases.
	 *
	 * @param cha First chromosome array to compare
	 * @param chb Second chromosome array to compare
	 * @return true if chromosomes are identical, false if any differences found
	 */
	public static boolean compare(ChromosomeArray cha, ChromosomeArray chb){
		boolean equal=true;
		if(cha.minIndex!=chb.minIndex || cha.maxIndex!=chb.maxIndex){
			System.out.println("Index mismatch in chrom "+cha.chromosome+":\n" +
					"("+cha.minIndex+" - "+cha.maxIndex+") vs ("+chb.minIndex+" - "+chb.maxIndex+")");
			equal=false;
		}
		int start=Data.max(cha.minIndex, chb.minIndex);
		int stop=Data.min(cha.maxIndex, chb.maxIndex);
		
		for(int i=start; i<=stop; i++){
			byte a=cha.get(i);
			byte b=chb.get(i);
			if(a!=b){
				System.out.println(((char)cha.chromosome)+"\t"+i+"\t"+((char)a)+" "+((char)b));
				equal=false;
			}
		}
		return equal;
		
	}
	
}
