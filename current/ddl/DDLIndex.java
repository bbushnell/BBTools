package ddl;

import cardinality.DynamicDemiLog;
import structures.IntList;

/**
 * Inverted index from DynamicDemiLog bucket values to clade IDs.
 * For each (bucket, value) pair, stores the list of clades with that
 * value in that bucket. Enables O(buckets) lookup of all clades sharing
 * k-mers with a query DDL.
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLIndex {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DDLIndex(){this(BUCKETS, VALUES);}

	public DDLIndex(int buckets, int values){
		matrix=new IntList[buckets][values];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a clade's DDL to the index.
	 * @param cladeID Numeric identifier for this clade
	 * @param ddl The clade's DynamicDemiLog */
	public void add(int cladeID, DynamicDemiLog ddl){
		add(cladeID, ddl.maxArray(), ddl.filledBuckets());
	}

	/** Adds a clade's DDL values to the index.
	 * @param cladeID Numeric identifier for this clade
	 * @param maxArray The DDL's 16-bit bucket values
	 * @param filled Number of non-empty buckets in this DDL */
	public void add(int cladeID, char[] maxArray, int filled){
		if(cladeID>=filledBuckets.length){
			filledBuckets=java.util.Arrays.copyOf(filledBuckets, Math.max(cladeID+1, filledBuckets.length*2));
		}
		filledBuckets[cladeID]=filled;
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			IntList list=matrix[b][v];
			if(list==null){matrix[b][v]=list=new IntList(2);}
			list.add(cladeID);
		}
		numClades=Math.max(numClades, cladeID+1);
	}

	/** Counts matching buckets between a query DDL and all indexed clades.
	 * @param query The query DynamicDemiLog
	 * @return int[] of match counts indexed by cladeID */
	public int[] query(DynamicDemiLog query){
		return query(query.maxArray());
	}

	/** Counts matching buckets between query values and all indexed clades.
	 * @param maxArray The query DDL's 16-bit bucket values
	 * @return int[] of match counts indexed by cladeID */
	public int[] query(char[] maxArray){
		final int[] counts=new int[numClades];
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			final IntList list=matrix[b][v];
			if(list==null){continue;}
			for(int i=0; i<list.size; i++){
				counts[list.get(i)]++;
			}
		}
		return counts;
	}

	/** Finds the top N clades by match count.
	 * @param query The query DynamicDemiLog
	 * @param maxHits Maximum number of results to return
	 * @return int[][2] array of {cladeID, matchCount} pairs, sorted by matchCount descending */
	public int[][] topHits(DynamicDemiLog query, int maxHits){
		final int[] counts=query(query);
		final int[][] hits=new int[maxHits][2];
		int minVal=0, minIdx=0;
		int filled=0;
		for(int id=0; id<counts.length; id++){
			final int c=counts[id];
			if(c<=minVal && filled>=maxHits){continue;}
			if(filled<maxHits){
				hits[filled][0]=id;
				hits[filled][1]=c;
				filled++;
				if(filled==maxHits){
					minVal=Integer.MAX_VALUE;
					for(int i=0; i<filled; i++){
						if(hits[i][1]<minVal){minVal=hits[i][1]; minIdx=i;}
					}
				}
			}else{
				hits[minIdx][0]=id;
				hits[minIdx][1]=c;
				minVal=Integer.MAX_VALUE;
				for(int i=0; i<filled; i++){
					if(hits[i][1]<minVal){minVal=hits[i][1]; minIdx=i;}
				}
			}
		}
		java.util.Arrays.sort(hits, 0, filled, (a, b)->b[1]-a[1]);
		if(filled<maxHits){return java.util.Arrays.copyOf(hits, filled);}
		return hits;
	}

	/** Computes ANI between query and a specific indexed clade.
	 * @param queryMatches Match count from query() for this clade
	 * @param queryFilled Number of filled buckets in the query DDL
	 * @param cladeID The clade to compare against
	 * @param k K-mer length */
	public float ani(int queryMatches, int queryFilled, int cladeID, int k){
		if(queryMatches<=0){return 0;}
		final int refFilled=filledBuckets[cladeID];
		final int minDiv=Math.min(queryFilled, refFilled);
		if(minDiv<=0){return 0;}
		final double c=Math.min(1.0, (double)queryMatches/minDiv);
		return (float)Math.exp(Math.log(c)/k);
	}

	/** Returns the number of indexed clades. */
	public int numClades(){return numClades;}

	/** Returns the number of non-null IntLists (for memory diagnostics). */
	public long populatedCells(){
		long count=0;
		for(int b=0; b<matrix.length; b++){
			for(int v=0; v<matrix[b].length; v++){
				if(matrix[b][v]!=null){count++;}
			}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final IntList[][] matrix;
	private int[] filledBuckets=new int[1024];
	private int numClades=0;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int BUCKETS=2048;
	private static final int VALUES=65536;
}
