package ddl;

import java.util.ArrayList;
import java.util.Arrays;

import cardinality.DynamicDemiLog;

/**
 * Matrix-backed reference implementation of the DDL sketch index:
 * matrix[bucket][value][] holds the int[] of clade IDs with that value in that bucket.
 * Simple and exact, but memory-heavy: a values-wide reference array per bucket plus one
 * small int[] object per populated cell. See {@link DDLIndexBase}; the lean alternative
 * is {@link DDLIndexCSR}, which must match this class's query()/topHits() exactly.
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLIndex extends DDLIndexBase {

	public DDLIndex(){this(BUCKETS, VALUES);}

	public DDLIndex(int buckets){this(buckets, VALUES);}

	public DDLIndex(int buckets, int values){
		super(buckets, values);
		matrix=new int[buckets][values][];
	}

	/** Adds one clade's DDL incrementally. Not used by the bulk path (addAll); kept for
	 *  completeness and as the reference for what a single insertion means. */
	public void add(int cladeID, DynamicDemiLog ddl){
		add(cladeID, ddl.toAbsoluteArray(), ddl.filledBuckets());
	}

	/** Adds one clade's 16-bit bucket values incrementally. */
	public void add(int cladeID, char[] maxArray, int filled){
		if(cladeID>=filledBuckets.length){
			filledBuckets=Arrays.copyOf(filledBuckets, Math.max(cladeID+1, filledBuckets.length*2));
		}
		filledBuckets[cladeID]=filled;
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			int[] arr=matrix[b][v];
			if(arr==null){
				matrix[b][v]=new int[]{cladeID};
			}else{
				int[] grown=Arrays.copyOf(arr, arr.length+1);
				grown[arr.length]=cladeID;
				matrix[b][v]=grown;
			}
		}
		numClades=Math.max(numClades, cladeID+1);
	}

	@Override
	protected long accumulateCounts(char[] maxArray, int[] counts){
		long work=0;
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			final int[] arr=matrix[b][v];
			if(arr==null){continue;}
			work+=arr.length;
			for(int i=0; i<arr.length; i++){
				counts[arr[i]]++;
			}
		}
		return work;
	}

	@Override
	protected void buildBucketRange(ArrayList<DDLRecord> records, int bStart, int bEnd){
		final int numRecords=records.size();
		final int[] counts=new int[values];
		final int[] pos=new int[values];

		for(int b=bStart; b<bEnd; b++){
			Arrays.fill(counts, 0);

			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){counts[v]++;}
			}

			for(int v=0; v<values; v++){
				if(counts[v]>0){matrix[b][v]=new int[counts[v]];}
			}
			Arrays.fill(pos, 0);

			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){
					matrix[b][v][pos[v]++]=i;
				}
			}
		}
	}

	@Override
	public long populatedCells(){
		long count=0;
		for(int b=0; b<matrix.length; b++){
			for(int v=0; v<matrix[b].length; v++){
				if(matrix[b][v]!=null){count++;}
			}
		}
		return count;
	}

	private final int[][][] matrix;
}
