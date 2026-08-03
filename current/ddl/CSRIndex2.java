package ddl;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;

/**
 * Bit-packed CSR sketch index: a third {@link DDLIndexBase} storage that packs BOTH
 * the clade IDs and the value offsets into 21-bit fields (3 per 64-bit word) via
 * {@link CSRLine}. About a third less memory than {@link DDLIndexCSR}'s int[] arrays
 * (32k RefSeq DDL: ~24.6 GB -> ~16.4 GB), for one shift+mask per read.
 *
 * Ceiling: 2^21-1 = 2,097,151 sketches -- IDs and cumulative offsets are both 21-bit.
 * A larger set crashes loud in {@link #addAll} (use csr=t, the 32-bit CSR, instead).
 * Bulk-build only (no incremental add). Returns query()/topHits() results identical to
 * {@link DDLIndex}: within a value the IDs are stored in ascending record order, exactly
 * as the matrix builder fills them, and ranking reads only the per-clade counts.
 *
 * @author Brian Bushnell, Noire
 * @date August 3, 2026
 */
public class CSRIndex2 extends DDLIndexBase {

	public CSRIndex2(){this(BUCKETS, VALUES);}

	public CSRIndex2(int buckets){this(buckets, VALUES);}

	public CSRIndex2(int buckets, int values){
		super(buckets, values);
		lines=new CSRLine[buckets];
	}

	@Override
	public void addAll(ArrayList<DDLRecord> records, int threads){
		if(records!=null && records.size()>MAX_CLADES){
			KillSwitch.kill("CSRIndex2 supports at most "+MAX_CLADES+" sketches (21-bit IDs); got "
				+records.size()+".  Use csr=t for the 32-bit CSR index.");
		}
		super.addAll(records, threads);
	}

	@Override
	protected long accumulateCounts(char[] maxArray, int[] counts){
		long work=0;
		for(int b=0; b<maxArray.length; b++){
			final int v=maxArray[b];
			if(v==0){continue;}
			final CSRLine line=lines[b];
			if(line==null){continue;}
			work+=line.accumulate(v, counts);
		}
		return work;
	}

	@Override
	protected void buildBucketRange(ArrayList<DDLRecord> records, int bStart, int bEnd){
		final int numRecords=records.size();
		final int[] counts=new int[values];//per-value counter, reused each bucket
		final int[] off=new int[values+1];//cumulative offsets, reused each bucket

		for(int b=bStart; b<bEnd; b++){
			Arrays.fill(counts, 0);
			int total=0;

			//Pass 1: count each value's occurrences in this bucket (value 0 = empty, skipped).
			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){counts[v]++; total++;}
			}

			//Prefix-sum into off: off[0]=0, off[v+1]=off[v]+counts[v].
			off[0]=0;
			for(int v=0; v<values; v++){off[v+1]=off[v]+counts[v];}

			//Pack the offsets and allocate the (empty) packed ID array for this bucket.
			final CSRLine line=new CSRLine(off, total, values);

			//Pass 2: write each record's ID into its value's slice (reuse counts as cursor).
			Arrays.fill(counts, 0);
			for(int i=0; i<numRecords; i++){
				final int v=records.get(i).ddl.maxArray()[b];
				if(v!=0){line.setID(off[v]+counts[v]++, i);}
			}
			lines[b]=line;
		}
	}

	@Override
	public long populatedCells(){
		long count=0;
		for(CSRLine line : lines){
			if(line!=null){count+=line.populatedCells();}
		}
		return count;
	}

	/** Max sketches: IDs and offsets are 21-bit fields, so at most 2^21-1. */
	static final int MAX_CLADES=(1<<CSRLine.BITS)-1;//2,097,151

	/** lines[b] = the bit-packed CSR row for bucket b (null if never built). */
	private final CSRLine[] lines;
}
