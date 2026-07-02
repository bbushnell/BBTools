package driver;
import java.util.ArrayList;
import java.util.List;

import dna.Data;
import dna.Gene;
import structures.Range;


public class Search {
	
	public static List<Gene> findGenes(int p, Gene[] genes){
		ArrayList<Gene> list=new ArrayList<Gene>(16);
		
		for(int i=0; i<genes.length; i++){
			Gene g=genes[i];
			if(g.intersectsCode(p)){
				list.add(g);
			}
		}
		
		return list;
	}
	
	/**
	 * Finds genes that overlap with the specified point using binary search on ranges.
	 * More efficient than linear search for large datasets with pre-computed ranges.
	 * Can search for exact coding region overlaps or nearby genes within threshold.
	 *
	 * @param p The genomic position to search for
	 * @param ranges Array of genomic ranges containing gene references
	 * @param nearby If true, searches for nearby genes; if false, exact code overlaps only
	 * @return List of genes that overlap or are near the specified point, null if no range contains p
	 */
	//NOTE [driver/Search#002] LOW/dead: findGenesBinary/findGenesLinear/findPointLinear/overlaps/findGenes are DEAD — the only
	//external referents (dna.Data, var.Variation) call ONLY findPointBinary + containsPointBinary. So these carry the same
	//empty-array AIOOBE as #001 but unreachably: findGenesBinary→findPointBinary returns 0→ranges[0] on empty; findGenesLinear→
	//findPointLinear returns array.length-1 = -1 on empty→ranges[-1]. Dead-latent twins of #001; left as-is.
	public static List<Gene> findGenesBinary(int p, Range[] ranges, boolean nearby){
		ArrayList<Gene> list=null;
		int a=findPointBinary(p, ranges);

		Range r=ranges[a];
		
//		System.out.println("Searching for "+p+" in "+r+"; previous range was "+ranges[a-1]);
		if(!r.includes(p)){return list;}
		
		list=new ArrayList<Gene>(16);
		
		Gene[] genes2=(Gene[])r.obj1;
		assert(genes2.length>0);
		
//		System.out.println("Found "+genes2.length+" to consider.");
		
		
		//TODO: Specify whether tx or code (etc) coverage is needed.
		for(int i=0; i<genes2.length; i++){
			Gene g=genes2[i];
//			System.out.print("Does p overlap gene "+g.codeStart+" - "+g.codeEnd+"?");
			if(g.txStart>r.b+Data.NEAR){break;}
			if(nearby){
				if(g.intersectsNearby(p, p)){list.add(g);}
			}else{
				if(g.intersectsCode(p)){list.add(g);}
			}
		}
		
		return list;
	}
	
	/**
	 * Finds genes that overlap with the specified point using linear search on ranges.
	 * Searches for exact coding region overlaps only, not nearby genes.
	 * Less efficient than binary search but simpler implementation.
	 *
	 * @param p The genomic position to search for
	 * @param genes Array of genes (appears unused in current implementation)
	 * @param ranges Array of genomic ranges containing gene references
	 * @return List of genes that overlap the specified point, null if no range contains p
	 */
	public static List<Gene> findGenesLinear(int p, Gene[] genes, Range[] ranges){
		ArrayList<Gene> list=null;
		int a=findPointLinear(p, ranges);
		
		Range r=ranges[a];
		
//		System.out.println("Searching for "+p+" in "+r+"; previous range was "+ranges[a-1]);
		if(!r.includes(p)){return list;}
		
		list=new ArrayList<Gene>(16);
		
		Gene[] genes2=(Gene[])r.obj1;
		assert(genes2.length>0);
		
//		System.out.println("Found "+genes2.length+" to consider.");
		
		
		//TODO: Specify whether tx or code (etc) coverage is needed.
		for(int i=0; i<genes2.length; i++){
			Gene g=genes2[i];
//			System.out.print("Does p overlap gene "+g.codeStart+" - "+g.codeEnd+"?");
			if(g.txStart>r.b){break;}
			if(g.intersectsCode(p)){
				list.add(g);
//				System.out.println(" Yes.");
			}
		}
		
		return list;
	}

	public static int findPointLinear(int p, Range[] array){
		for(int i=0; i<array.length; i++){
			Range r=array[i];
			if(r.a>p){return i;} //Fail.
			if(r.includes(p)){return i;} //Success.
		}
		return array.length-1;
	}
	
	public static int findPointBinary(int p, Range[] array){
		assert(array!=null);
		//TODO: Possible bug [driver/Search#001] LOW→MEDIUM-if-reachable/ESCALATE: on an empty array this returns 0 — but 0
		//is NOT a valid index for a length-0 array, and every LIVE caller blindly indexes the result:
		//  · dna.Data.getNearestGeneSets(chrom,loc) L632: `r[index].obj1` — NO empty guard → r[0] AIOOBE.
		//  · var.Variation.filterCodingVariances L153: `ranges[rnum].intersects(...)` — NO empty guard → ranges[0] AIOOBE.
		//By contrast dna.Data.getNearestGeneSets(chrom,loc1,loc2) L656 DOES guard (`if(ranges==null||length==0)return null`),
		//and the sibling containsPointBinary below CORRECTLY guards empty→false. So the empty case is a shape the author handled
		//in 2 of 3 sibling sites but the findPointBinary contract + 2 callers miss it. Reachable iff geneSetRangeMatrix(chrom)/
		//geneCodeAndExon/geneNearby can be length 0 (a chrom with zero gene ranges) — the guarded twin implies the author thought
		//it could. Cross-package fix (driver+dna+var) & the return-on-empty is a design choice → ESCALATED, not fixed here.
		if(array.length==0){return 0;}
		int result=findPointBinary(p, 0, max(0, array.length-1), array);
		
		//TODO: Assertions
		
		return result;
	}
	
	public static boolean containsPointBinary(int p, Range[] array, int thresh){
		assert(array!=null);
		//n GOOD (studied praise / contrast to #001): THIS sibling guards empty correctly — returns false, never indexes.
		//n Same author, same file; the divergence with findPointBinary's empty→0 is what flags #001 as a real oversight.
		if(array.length==0){return false;}
		int rnum=findPointBinary(p, 0, max(0, array.length-1), array);
		
		int p1=p-thresh, p2=p+thresh;
		Range r=array[rnum];
		if(p2>=r.a && p1<=r.b){return true;}
		
		if(rnum==0 && p<r.a){return false;}
		
		assert(p>r.b) : "\n\n"+p+"\t"+rnum+"/"+array.length+"\t"+r+"\n\n"; //Otherwise, it violated the search contract.
		if(array.length<=rnum+1){return false;}
		
		Range r2=array[rnum+1];
		assert(r2.a>p) : "\n\n"+p+"\t"+rnum+"/"+array.length+"\t"+r+"\n\n"; //Otherwise, it violated the search contract.
		return (p2>=r.a && p1<=r.b);
	}
	
	public static int findPointBinary(int p, int a, int b, Range[] array){
		if(a>=b){
			
			//This line should ensure that p>array[a] when p is not within any range.
			//Except, of course, when p<(all ranges).
			//In other words, the return is strictly the LEFT (LOWER) index when p is between two ranges.
			if(a>0 && p<array[a].a){a--;}
			
			assert(a>=0);
			assert(a<array.length);
			assert(array[a].includes(p) || (a==0 && p<array[a].a) ||
					(p>array[a].b && (a==array.length-1 || p<array[a+1].a))) :
						"a="+a+", b="+b+", p="+p+", array[a]="+array[a];
			return a;
		}
		
		int mid=(a+b)/2;
		Range r=array[mid];
		
		if(r.a>p){
			return findPointBinary(p, a, mid-1, array);
		}else if(r.b<p){
			return findPointBinary(p, mid+1, b, array);
		}else{
			return mid;
		}
	}
	
	
	public static boolean overlaps(int a, Gene g){
		return a>=g.txStart && a<=g.txStop;
	}
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	
	
}
