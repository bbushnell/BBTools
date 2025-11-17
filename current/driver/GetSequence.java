package driver;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import shared.Shared;
import shared.Tools;
import structures.Range;

/**
 * Command-line utility for extracting and manipulating genomic sequences from chromosomes.
 * Retrieves genomic sequence segments by chromosome, location, and optional strand.
 * Supports flexible chromosome and location input parsing with genomic build selection.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class GetSequence {
	
	/**
	 * Main entry point for genomic sequence extraction.
	 * Parses command-line arguments for chromosome, build, strand, and coordinates.
	 * Extracts sequences from specified genomic regions and outputs original and
	 * reverse complement sequences with amino acid translations.
	 *
	 * Supports command-line format:
	 * - build=N (genome build selection)
	 * - chrN (chromosome specification)
	 * - + or - (strand specification, currently only + supported)
	 * - coordinate ranges in formats: start-end, [start,end], (start,end)
	 *
	 * @param args Command-line arguments including build, chromosome, strand, and coordinates
	 */
	public static void main(String[] args){
		
		int chrom=-1;
		byte strand=Shared.PLUS;
		
		/** Change base to zero or one for the coordinates mode */
		int base=0;
		
//		char c=args[1].charAt(0);
//		if(c=='+'){strand=Gene.PLUS;}
//		else if(c=='-'){strand=Gene.MINUS;}
//		else{assert(false) : "Invalid strand: "+args[1];}
		
		int firstLoc=-1;
		for(int i=0; i<args.length && firstLoc<0; i++){
			char x=args[i].charAt(0);
			if(args[i].startsWith("build")){
				final String arg=args[i];
				final String[] split=arg.split("=");
				Data.setGenome(Integer.parseInt(split[1]));
			}else if(Tools.isDigit(x) || x=='[' || x=='('){
				firstLoc=i;
			}else{
				if(args[i].startsWith("chr")){
					chrom=Gene.toChromosome(args[i]);
				}else if(x=='b'){
//					Data.GENOME_BUILD=Gene.toBuild(args[i]);
					Data.setGenome(Gene.toBuild(args[i]));
				}else if(x=='+' || x=='-'){
					strand=Gene.toStrand(args[i]);
				}else{
					assert(false) : "Bad parameter: "+args[i];
				}
			}
		}

//		System.out.println(Data.GENOME_BUILD);
//		System.out.println(chrom);
//		System.out.println(strand);
		
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray cha=Data.getChromosome(chrom);
		
		String[] array=new String[args.length];
		
//		System.out.println("firstLoc="+firstLoc+"/"+args.length);
		for(int i=firstLoc; i<args.length; i++){
//			System.out.println("Processing "+args[i]);
			args[i]=args[i].replace("[","").replace("]","").replace("(","").replace(")","").replace(",","").trim();
			
			Range r=Range.toRange(args[i]);
			array[i]=cha.getString((int)r.a-base, (int)r.b-base);
		}
		
		String combined="";
		
		System.out.println("Chrom Bounds: "+cha.minIndex+"-"+cha.maxIndex+" ("+cha.array.length+")");
		
		for(int i=firstLoc; i<array.length; i++){
			System.out.println("\nchr"+chrom+" ("+args[i]+") = \n"+array[i]);
//			System.out.println("AAs:\t"+AminoAcid.stringToAAs(array[i]));
			combined+=array[i];
			String s=AminoAcid.reverseComplementBases(array[i]);
			System.out.println("\n"+s+" (rcomp)");
//			System.out.println("AAs:\t"+AminoAcid.stringToAAs(s));
		}

		System.out.println("\nAAs:\n"+AminoAcid.stringToAAs(combined));
		System.out.println("\nAAs (reverse comp):\n"+AminoAcid.stringToAAs(AminoAcid.reverseComplementBases(combined)));
	}
	
	/**
	 * Retrieves a single base from a specified chromosome at a given position.
	 * @param chrom Chromosome number (1-based)
	 * @param a Genomic position on the chromosome
	 * @return Base value at the specified position
	 */
	public static byte get(int chrom, int a){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return cha.get(a);
	}
	
	/**
	 * Retrieves a genomic sequence from a chromosome between two positions.
	 * Uses positive strand by default.
	 *
	 * @param chrom Chromosome number (1-based)
	 * @param a Start position (inclusive)
	 * @param b End position (inclusive)
	 * @return Genomic sequence as a string between positions a and b
	 */
	public static String get(int chrom, int a, int b){
		return get(chrom, a, b, Shared.PLUS);
	}
	
	/**
	 * Retrieves a genomic sequence from a chromosome between two positions on specified strand.
	 * Currently only supports positive strand (Shared.PLUS).
	 *
	 * @param chrom Chromosome number (1-based)
	 * @param a Start position (inclusive)
	 * @param b End position (inclusive)
	 * @param strand Strand orientation (currently only Shared.PLUS supported)
	 * @return Genomic sequence as a string between positions a and b
	 */
	public static String get(int chrom, int a, int b, byte strand){
		assert(strand==Shared.PLUS) : "TODO";
		ChromosomeArray cha=Data.getChromosome(chrom);
		return cha.getString(a, b);
	}
	
	
	
	
	
}
