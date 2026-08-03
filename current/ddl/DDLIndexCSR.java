package ddl;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Compressed-sparse-row storage for the DDL sketch index. Per bucket it keeps two int[]:
 *   packed[b]  -- every clade ID in that bucket, grouped by value in value order, and
 *   offsets[b] -- int[values+1]; value v occupies packed[b][offsets[b][v] .. offsets[b][v+1]),
 *                 so count(b,v) = offsets[b][v+1] - offsets[b][v].
 * That is two arrays per bucket instead of a values-wide reference array plus one small int[]
 * per populated cell, so it eliminates the per-cell object header/padding overhead (the bulk of
 * the matrix index's footprint) and reads each bucket's IDs contiguously.
 *
 * Bulk-build only (no incremental add). Returns query()/topHits() results identical to
 * {@link DDLIndex}: within a value the IDs are stored in ascending record order, exactly as the
 * matrix builder fills them, and ranking only reads the per-clade counts. See {@link DDLIndexBase}.
 *
 * @author Noire
 * @date August 3, 2026
 */
public class DDLIndexCSR extends DDLIndexBase {

	public DDLIndexCSR(){this(BUCKETS, VALUES);}

	public DDLIndexCSR(int buckets){this(buckets, VALUES);}

	public DDLIndexCSR(int buckets, int values){
		super(buckets, values);
		packed=new int[buckets][];
		offsets=new int[buckets][];
	}

	@Override
	protected long accumulateCounts(char[] maxArray, int[] counts){
		long work=0;
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			final int[] off=offsets[b];
			if(off==null){continue;}
			final int lo=off[v], hi=off[v+1];
			if(hi<=lo){continue;}
			final int[] pk=packed[b];
			work+=(hi-lo);
			for(int j=lo; j<hi; j++){
				counts[pk[j]]++;
			}
		}
		return work;
	}

	@Override
	protected void buildBucketRange(ArrayList<DDLRecord> records, int bStart, int bEnd){
		final int numRecords=records.size();
		final int[] counts=new int[values];//scratch, reused each bucket

		for(int b=bStart; b<bEnd; b++){
			Arrays.fill(counts, 0);
			int total=0;

			//Pass 1: count occurrences of each value in this bucket.
			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){counts[v]++; total++;}
			}

			//Prefix-sum the counts into the offset array: off[v+1]=off[v]+counts[v].
			final int[] off=new int[values+1];
			for(int v=0; v<values; v++){off[v+1]=off[v]+counts[v];}

			//Pass 2: place each record's ID into its value's slice. Reuse counts as the
			//per-value running cursor (it is free once off is computed).
			final int[] pk=new int[total];
			Arrays.fill(counts, 0);
			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){pk[off[v]+counts[v]++]=i;}
			}

			offsets[b]=off;
			packed[b]=pk;
		}
	}

	@Override
	public long populatedCells(){
		long count=0;
		for(int b=0; b<offsets.length; b++){
			final int[] off=offsets[b];
			if(off==null){continue;}
			for(int v=0; v<values; v++){
				if(off[v+1]>off[v]){count++;}
			}
		}
		return count;
	}

	/** packed[b] = clade IDs for bucket b, grouped by value (ascending record order within a value). */
	private final int[][] packed;
	/** offsets[b] = int[values+1] prefix sums; value v is packed[b][offsets[b][v]..offsets[b][v+1]). */
	private final int[][] offsets;
}
