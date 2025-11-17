package driver;

import java.util.ArrayList;
import java.util.List;

import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;
import structures.CoverageArray2;
import var.VarLine;
import var.Variation;

/**
 * Generates no-call variant entries from coverage data.
 * Creates VarLine objects representing genomic regions with insufficient coverage
 * for confident variant calling, handling both haploid and diploid chromosomes.
 * @author Brian Bushnell
 */
public class GenerateNoCallsFromCoverage {
	
//	@Deprecated
//	public static ArrayList<VarLine> generateOld(byte chrom, CoverageArray ca, int build, char gender){
//
//		ArrayList<VarLine> lines=new ArrayList<VarLine>(256);
//
//		assert(Data.GENOME_BUILD==build);
//		ChromosomeArray chra=Data.getChromosome(chrom);
//
//		int start=-1;
//		int stop=-1;
//
//		for(int i=chra.minIndex; i<chra.maxIndex; i++){
//			boolean nc=(ca.get(i)<minCovered);
//
//			if(nc && start==-1){
//				start=i;
//			}
//
//			if(!nc && start>-1){
//				stop=i-1;
//
//				VarLine v1=new VarLine();
//				v1.ploidy=(chrom<=22 ? 2 : chrom>=24 ? 1 : (Byte)Variation.ploidyMap.get("?"));
//
//				v1.haplotype=1;
//				v1.chromosome=chrom;
//				v1.beginLoc=start;
//				v1.endLoc=stop;
//
//				v1.ref="=";
//				v1.call=null;
//
//				v1.totalScore=-1;
//				v1.xRef=-2;
//				v1.xRefArray=null;
//				v1.hapLink=-1;
//				v1.varType=Variation.NOCALL;
//
//				VarLine v2;
//				if((chrom==23 && gender=='M') || chrom==24 || chrom==25){
//					v2=null;
//				}else{
//					v2=new VarLine();
//					v2.ploidy=(chrom<=22 ? 2 : chrom>=24 ? 1 : (Byte)Variation.ploidyMap.get("?"));
//
//					v2.haplotype=2;
//					v2.chromosome=chrom;
//					v2.beginLoc=start;
//					v2.endLoc=stop;
//
//					v2.ref="=";
//					v2.call=null;
//
//					v2.totalScore=-1;
//					v2.xRef=-2;
//					v2.xRefArray=null;
//					v2.hapLink=-1;
//					v2.varType=Variation.NOCALL;
//				}
//
//
//				start=-1;
//				stop=-1;
//				lines.add(v1);
//				if(v2!=null){lines.add(v2);}
//			}
//
//
//		}
//
//		return lines;
//	}
	
	
	
	/**
	 * Generates no-call VarLine entries for regions with insufficient coverage.
	 * Analyzes coverage data to identify regions below the minimum thresholds and
	 * creates appropriate no-call entries for each haplotype based on chromosome type.
	 *
	 * @param chrom Chromosome number (1-25, where 23=X, 24=Y, 25=MT)
	 * @param ca Coverage array containing per-base coverage counts
	 * @param build Genome build number (must match Data.GENOME_BUILD)
	 * @param gender Sample gender ('M' or 'F'), affects haploid chromosome handling
	 * @return List of VarLine objects representing no-call regions
	 */
	public static ArrayList<VarLine> generate(byte chrom, CoverageArray2 ca, int build, char gender){
		
		assert(minCovered>=1);
		assert(minHalfCovered>=1);
		assert(minCovered>=minHalfCovered);
		
		
		ArrayList<VarLine> lines=new ArrayList<VarLine>(256);
		
		assert(Data.GENOME_BUILD==build);
		ChromosomeArray chra=Data.getChromosome(chrom);
		
		int start=-1;
		int stop=-1;
		
		
		byte level=-1;
		
		boolean haploid=(chrom==23 && gender=='M') || chrom==24 || chrom==25;
		
		for(int i=chra.minIndex; i<chra.maxIndex; i++){
			
			final byte newLevel;
			final int cov=ca.get(i);
			
			if(haploid){
				if(cov<minHalfCovered){
					newLevel=0;
				}else{
					newLevel=2;
				}
			}else{
				if(cov<minHalfCovered){
					newLevel=0;
				}else if(cov<minCovered){
					newLevel=1;
				}else{
					newLevel=2;
				}
			}
			
			
			if(level==-1){
				level=newLevel;
				start=i;
			}else if(level!=newLevel){ //The level changed; make VarLines
				
				stop=i-1;
				
				if(level==0){

					VarLine v1=new VarLine();
					v1.ploidy=(chrom<=22 ? 2 : chrom>=24 ? 1 : gender=='M' ? 1 : gender=='F' ? 2 : (Byte)Variation.ploidyMap.get("?"));
					v1.haplotype=1;
					v1.chromosome=chrom;
					v1.beginLoc=start;
					v1.endLoc=stop;

					v1.ref="=";
					v1.call=null;

					v1.totalScore=-1;
					v1.hapLink=-1;
					v1.varType=Variation.NOCALL;

					lines.add(v1);
				}
				
				if(level==0 || (level==1 && !haploid)){

					VarLine v2=new VarLine();
					v2.ploidy=(chrom<=22 ? 2 : chrom>=24 ? 1 : gender=='M' ? 1 : gender=='F' ? 2 : (Byte)Variation.ploidyMap.get("?"));

					v2.haplotype=2;
					v2.chromosome=chrom;
					v2.beginLoc=start;
					v2.endLoc=stop;

					v2.ref="=";
					v2.call=null;

					v2.totalScore=-1;
					v2.hapLink=-1;
					v2.varType=Variation.NOCALL;
					
					lines.add(v2);
				}


//				start=-1;
				stop=-1;
				level=newLevel;
				start=i;
			}
		}
		
		return lines;
	}
	
	
	/**
	 * Removes overlapping no-call entries by splitting haplotypes and processing separately.
	 * Splits the input by haplotype, removes duplicates within each haplotype,
	 * then recombines and sorts the results.
	 *
	 * @param input List of VarLine objects that may contain overlapping no-calls
	 * @param copies Number of chromosome copies (affects haplotype splitting logic)
	 * @return Deduplicated and sorted list of VarLine objects
	 */
	public static ArrayList<VarLine> removeDuplicateNocalls(List<VarLine> input, int copies){
		ArrayList<VarLine>[] haplo=splitHaplotypes(input, copies);

		ArrayList<VarLine> output=new ArrayList<VarLine>(256);
//		System.err.println("A: copies="+copies+"; input.size="+input.size()+"; haplo="+haplo[0].size()+", "+haplo[1].size());
		for(ArrayList<VarLine> alv : haplo){
			VarLine temp=alv.size()==0 ? null : alv.get(0);
			for(VarLine vl : alv){assert(vl.haplotype==temp.haplotype);}
			ArrayList<VarLine> alv2=removeDuplicateNocallsHaplotyped(alv);
//			assert(checkCopyCountHaplotyped(alv2)); //Very slow
			
//			output.addAll(removeDuplicateNocallsHaplotyped(alv2)); //This MUST be incorrect.
			
			output.addAll(alv2);
		}
		
		Shared.sort(output);
		
		return output;
	}
	
	/**
	 * Validates that no genomic position is covered by multiple VarLines within a haplotype.
	 * Checks for overlapping entries and nocall duplicates which would indicate
	 * processing errors in the variant calling pipeline.
	 *
	 * @param list List of VarLine objects from a single haplotype to validate
	 * @return true if no overlaps are found, false if validation fails
	 */
	public static boolean checkCopyCountHaplotyped(List<VarLine> list){
		
		int max=0;
		for(VarLine vl : list){
			if(vl.endLoc>max){max=vl.endLoc;}
		}

		byte[] sum=new byte[max+1];
//		byte[] vars=new byte[max+1];
		byte[] nocalls=new byte[max+1];
		
		for(VarLine vl : list){
			for(int i=vl.beginLoc; i<=vl.endLoc; i++){
				sum[i]++;
				if(vl.isNoCall()){nocalls[i]++;}
//				else{vars[i]++;}
			}
		}
		
		for(int i=0; i<sum.length; i++){
			if(nocalls[i]>1){
				assert(false) : "chr"+list.get(0).chromosome+", "+i;
				return false;
			}
			if(sum[i]>1){
				assert(false) : "chr"+list.get(0).chromosome+", "+i;
				return false;
			}
		}
		
		return true;
	}
	
	
	/** All elements of input should share haplotype */
	public static ArrayList<VarLine> removeDuplicateNocallsHaplotyped(ArrayList<VarLine> input){
		

//		System.err.println("B: input.size="+input.size());
		
		Shared.sort(input);
		
		ArrayList<VarLine> output=new ArrayList<VarLine>(256);

		boolean needToReprocess=false;

		VarLine prev=null;
		
		final boolean verbose=false;
		
		for(int i=0; i<input.size(); i++){
			VarLine current=input.get(i);
			
			assert(current.endLoc>=current.beginLoc) : current;
			
//			final VarLine current2=current;
			final VarLine prev2=prev;
			
//			if(current.chromosome==2 && (current.touches(8890433) || (prev!=null && prev.touches(8890433)))){
//				verbose=true;
//				System.err.println("current="+current);
//				System.err.println("touches? "+current.touches(8890433));
//				System.err.println("intersects? "+current.intersects(8890433));
//			}else if(prev==null && verbose){
//				System.err.println("current="+current);
//				System.err.println("touches? "+current.touches(8890433));
//				System.err.println("intersects? "+current.intersects(8890433));
//			}else{
//				verbose=false;
//			}
			
			boolean problem=prev!=null && prev.intersects(current);
			if(problem){
				if(prev.isPoint() && (current.endLoc==prev.beginLoc || current.beginLoc==prev.beginLoc)){
					problem=false;
				}
				if(current.isPoint() && (prev.endLoc==current.beginLoc || prev.beginLoc==current.beginLoc)){
					problem=false;
				}
			}
			
			if(problem){
				boolean ncc=current.isNoCall();
				boolean ncp=prev.isNoCall();
				boolean refc=current.isRef();
				boolean refp=prev.isRef();
				boolean varc=current.isTrueVariation();
				boolean varp=prev.isTrueVariation();
				if(!needToReprocess){
//					System.err.println("\nNeed to reprocess because:");
//					System.err.println("\n"+prev);
//					System.err.println("\n"+current);
				}
				needToReprocess=true;
				
				if((ncc && ncp) || (refc && refp) || (refc && ncp)/* || (refc && varp) || (ncp && varp)*/){ //Un-intersect them
					current=current.clone();
					{
						current.ref="=";
						if(refc){current.call="=";}
						else if(ncc){current.call=null;}
					}
					current.beginLoc=prev.endLoc+1;
					if(current.beginLoc>current.endLoc){current=null;}
					else{
						assert(!prev.intersects(current)
						|| (prev.isPoint() && (current.endLoc==prev.beginLoc || current.beginLoc==prev.beginLoc))
						|| (current.isPoint() && (prev.endLoc==current.beginLoc || prev.beginLoc==current.beginLoc))) :
							refp+", "+ncp+", "+refc+", "+ncc+"\n"+prev+"\n"+current;
					}
				}else if(ncc || refc){
					current=current.clone();
					{
						current.ref="=";
						if(refc){current.call="=";}
						else if(ncc){current.call=null;}
					}
					current.beginLoc=prev.endLoc+(prev.isPoint() ? 0 : 1);
					if(current.beginLoc>current.endLoc){current=null;}
					else{
						assert(!prev.intersects(current)
						|| (prev.isPoint() && (current.endLoc==prev.beginLoc || current.beginLoc==prev.beginLoc))
						|| (current.isPoint() && (prev.endLoc==current.beginLoc || prev.beginLoc==current.beginLoc))) :
							refp+", "+ncp+", "+refc+", "+ncc+"\n"+prev+"\n"+current;
					}
				}else if(ncp || refp){
					prev=prev.clone();
					{
						prev.ref="=";
						if(refp){prev.call="=";}
						else if(ncp){prev.call=null;}
					}
					prev.endLoc=current.beginLoc-1;
					if(prev.beginLoc>prev.endLoc){prev=null;}
					else{
						assert(!prev.intersects(current) ||
								(prev.isNoCall() && prev.lengthRef()==1 && current.isPoint()) //Corner case for intersection
								) : "\n"+prev+"\n\n"+current+"\n";
					}
					
					if(prev2.endLoc>current.endLoc || (prev2.endLoc==current.endLoc && current.isPoint())){
						VarLine temp=prev2.clone();
						{
							temp.ref="=";
							if(temp.isRef()){temp.call="=";}
							else if(temp.isNoCall()){temp.call=null;}
						}
						temp.beginLoc=current.endLoc+(current.isPoint() ? 0 : 1);
						if(temp.beginLoc<=temp.endLoc){

							assert(prev==null || !temp.intersects(prev));
							assert(!temp.intersects(current)
									|| (temp.isPoint() && (current.endLoc==temp.beginLoc || current.beginLoc==temp.beginLoc))
									|| (current.isPoint() && (temp.endLoc==current.beginLoc || temp.beginLoc==current.beginLoc))) :
										refp+", "+ncp+", "+refc+", "+ncc+"\n"+temp+"\n"+current;

							if(verbose){System.err.println("Current="+current+"\nprev="+prev+"\nAdding "+temp+"\n");}
							
							output.add(temp);
//							needToReprocess=true;
						}
					}
				}else{
					System.out.println("Warning: Deleted variation due to conflict! \n"+prev+"\n"+current+"\n");
//					assert(false) : "\n"+prev+"\n"+current+"\n";
					current=null;
				}
			}

			if(prev!=null){
				if(verbose){System.err.println("Current="+current+"\nAdding "+prev+"\n");}
				output.add(prev);
			}
			prev=current;
		}
		if(prev!=null){output.add(prev);}
		
		if(needToReprocess){return removeDuplicateNocallsHaplotyped(output);}
		
		Shared.sort(output);
		return output;
	}
	
	
	/**
	 * Separates VarLine objects into haplotype-specific arrays.
	 * Splits entries marked as haplotype 3 (both haplotypes) into separate
	 * entries for haplotypes 1 and 2 when multiple copies are present.
	 *
	 * @param input List of VarLine objects with mixed haplotype assignments
	 * @param copies Number of chromosome copies (affects splitting of haplotype 3)
	 * @return Array of two ArrayLists, [0] for haplotype 1, [1] for haplotype 2
	 */
	public static ArrayList<VarLine>[] splitHaplotypes(List<VarLine> input, int copies){
		ArrayList<VarLine>[] haplo=new ArrayList[2];
		for(int i=0; i<haplo.length; i++){
			haplo[i]=new ArrayList<VarLine>();
		}
		for(VarLine vl : input){
			if(vl.haplotype==1){
				haplo[0].add(vl);
			}else if(vl.haplotype==2){
				haplo[1].add(vl);
			}else{
				assert(vl.haplotype==3);
				if(copies>1){
					VarLine[] vl2=vl.splitLine();
					haplo[0].add(vl2[0]);
					haplo[1].add(vl2[1]);
				}else{
					haplo[0].add(vl);
				}
			}
		}
		
		return haplo;
	}
	

	/**
	 * Minimum coverage required for both haplotypes to be callable (diploid regions)
	 */
	public static int minCovered=2;
	/**
	 * Minimum coverage required for one haplotype to be callable (haploid regions)
	 */
	public static int minHalfCovered=1;
	
}
