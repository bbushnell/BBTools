package driver;
import java.util.ArrayList;
import java.util.List;

import dna.Data;
import dna.Gene;
import structures.Range;


/**
 * Utility class for searching genes within genomic ranges using various algorithms.
 * Provides linear and binary search implementations for efficient gene lookups by position.
 * Supports both exact overlaps and nearby gene searches with distance thresholds.
 * @author Brian Bushnell
 */
public class Search {
	
	/** Find genes in the array that overlap point "p" */
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
	
	/** Find genes in the array that overlap point "p" */
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
	
	/** Find genes in the array that overlap point "p" */
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

	/**
	 * Finds the index of the range containing the specified point using linear search.
	 * Returns the first range that includes the point or the first range starting after p.
	 *
	 * @param p The genomic position to locate
	 * @param array Array of genomic ranges to search
	 * @return Index of the range containing p, or last index if p is beyond all ranges
	 */
	public static int findPointLinear(int p, Range[] array){
		for(int i=0; i<array.length; i++){
			Range r=array[i];
			if(r.a>p){return i;} //Fail.
			if(r.includes(p)){return i;} //Success.
		}
		return array.length-1;
	}
	
	/**
	 * Finds the index of the range containing the specified point using binary search.
	 * More efficient than linear search for large range arrays.
	 *
	 * @param p The genomic position to locate
	 * @param array Array of genomic ranges to search
	 * @return Index of the range containing p or the appropriate boundary index
	 */
	public static int findPointBinary(int p, Range[] array){
		assert(array!=null);
		if(array.length==0){return 0;}
		int result=findPointBinary(p, 0, max(0, array.length-1), array);
		
		//TODO: Assertions
		
		return result;
	}
	
	/**
	 * Determines if any range contains the specified point within a threshold distance.
	 * Uses binary search to locate the appropriate range and checks threshold boundaries.
	 *
	 * @param p The genomic position to check
	 * @param array Array of genomic ranges to search
	 * @param thresh Distance threshold for proximity matching
	 * @return true if point is within threshold distance of any range, false otherwise
	 */
	public static boolean containsPointBinary(int p, Range[] array, int thresh){
		assert(array!=null);
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
	
	/**
	 * Recursive binary search implementation for finding point within range array.
	 * Returns the left (lower) index when point falls between ranges.
	 * Handles boundary conditions for points outside all ranges.
	 *
	 * @param p The genomic position to locate
	 * @param a Lower bound index for search
	 * @param b Upper bound index for search
	 * @param array Array of genomic ranges to search
	 * @return Index of range containing p or appropriate boundary index
	 */
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
	
	
	/**
	 * Checks if a genomic position overlaps with a gene's transcription boundaries.
	 * Tests if position falls within the gene's txStart and txStop coordinates.
	 *
	 * @param a The genomic position to test
	 * @param g The gene to check for overlap
	 * @return true if position overlaps gene transcription region, false otherwise
	 */
	public static boolean overlaps(int a, Gene g){
		return a>=g.txStart && a<=g.txStop;
	}
	/** Returns the smaller of two integer values */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integer values */
	private static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller of two long values */
	private static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger of two long values */
	private static final long max(long x, long y){return x>y ? x : y;}
	
	
}
