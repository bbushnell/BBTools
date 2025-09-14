package driver;

import dna.ChromosomeArray;
import dna.Data;

public class CompareReferenceGenomes {
	
	public static void main(String[] args){
		compareGenomes(args[0], args[1]);
	}
	
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
