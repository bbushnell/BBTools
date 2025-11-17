package driver;

import dna.Data;
import dna.Gene;
import fileIO.ChainLine;
import shared.Shared;

/**
 * Genome coordinate translation utility for converting between different genome builds.
 * Supports translation between genome builds 36 and 37, handling chromosome and strand mapping
 * using ChainLine mapping for precise coordinate conversion.
 * @author Brian Bushnell
 */
public class Translator2 {
	
	
	/**
	 * Program entry point for coordinate translation.
	 * Translates genomic coordinates from command-line arguments between genome builds.
	 * Expected arguments: fromBuild toBuild chromosome location1 location2 ...
	 * @param args Command-line arguments: [fromBuild] [toBuild] [chromosome] [location...]
	 */
	public static void main(String[] args){
		
		int from=Gene.toBuild(args[0]);
		int to=Gene.toBuild(args[1]);

		if(from==18){from=36;}
		if(from==19){from=37;}
		if(to==18){to=36;}
		if(to==19){to=37;}
		assert(from!=to);
		assert(from==36 || from==37);
		assert(to==36 || to==37);
		
		int chrom=Gene.toChromosome(args[2]);
		
		ChainLine[][] lines=Data.getChainLines(from, to);

		for(int i=3; i<args.length; i++){
			int loc=Integer.parseInt(args[i]);
			int[] result=ChainLine.translate(loc, lines[chrom]);
			System.out.print("(build"+from+", chr"+chrom+", +, "+loc+")  ->  ");
			System.out.println(result==null ? "null" :
				"(build"+to+", chr"+result[0]+", "+Shared.strandCodes[result[1]]+", "+result[2]+")");
		}
		
//		Translator2 tr=new Translator2(from, to);
//
//		ChainLine[] array=lines[chrom];
//		int index=ChainLine.binarySearch(loc, array);
////		if(index<0){return null;}
//		ChainLine cl=array[index];
//
////		System.out.println(cl);
//
//		int[] dest=cl.translate(loc);
//
////		{qChrom, qStrand, qStart+loc-tStart};
//
//		System.out.println(chrom+", +, "+loc+"   ->   "+dest[0]+", "+Gene.strandCodes[dest[1]]+", "+dest[2]);
	}
	
	/** chrom, strand, loc */
	public static final int[] translate(int fromBuild, int toBuild, int chrom, int strand, int loc){
		ChainLine[][] lines=Data.getChainLines(fromBuild, toBuild);
		int[] result=ChainLine.translate(loc, lines[chrom]);
		if(result==null){return null;}
		int strand2=result[1];
		if(strand2==strand){
			result[1]=Shared.PLUS;
		}else{
			result[1]=Shared.MINUS;
		}
		return result;
	}
	
}
