package ddl;

/**
 * One CSR "row" (a single DDL bucket) with BOTH arrays bit-packed at 21 bits.
 * Where {@link DDLIndexCSR} keeps two int[] per bucket -- the clade IDs and the
 * value->slice offsets -- this keeps two long[] holding 21-bit fields, 3 to a
 * 64-bit word (fields at bits 0, 21, 42; the top bit is unused). That is 21.33
 * bits per entry instead of 32, about a third less memory, for one shift+mask on
 * each read. No field ever spans a word, so a read is one load, one shift, one
 * mask. Supports clade IDs and offsets up to 2^21-1 = 2,097,151 (see {@link CSRIndex2}).
 *
 * Layout for this bucket:
 *   packedOffsets -- (numValues+1) cumulative offsets; value v's IDs occupy slots
 *                    [offsetAt(v), offsetAt(v+1)) of packedIDs.
 *   packedIDs     -- the clade IDs, grouped by value in value order.
 * count(v) = offsetAt(v+1)-offsetAt(v). Value 0 (empty bucket) is never stored,
 * so offsetAt(0)==offsetAt(1)==0.
 *
 * @author Brian Bushnell, Noire
 * @date August 3, 2026
 */
public final class CSRLine {

	/** Packs plainOffsets (length numValues+1, each in [0, 2^21-1]) and allocates an
	 *  empty packedIDs sized for 'total' IDs; setID() fills it during the build. */
	CSRLine(int[] plainOffsets, int total, int numValues){
		this.numValues=numValues;
		final int nOff=numValues+1;
		packedOffsets=new long[words(nOff)];
		for(int i=0; i<nOff; i++){
			packedOffsets[i/3]|=((long)plainOffsets[i]&MASK)<<shift(i);
		}
		packedIDs=new long[words(total)];
	}

	/** Build-only: writes clade id into slot. Each slot is written exactly once and
	 *  packedIDs starts zeroed, so an OR is a set. */
	void setID(int slot, int id){
		packedIDs[slot/3]|=((long)id&MASK)<<shift(slot);
	}

	/** The cumulative offset stored at index i (i in [0, numValues]). */
	int offsetAt(int i){
		return (int)((packedOffsets[i/3]>>>shift(i))&MASK);
	}

	/** For value v, increments counts[id] for every clade id holding v in this bucket.
	 *  Returns the number of IDs visited (index work). */
	long accumulate(int v, int[] counts){
		final int lo=offsetAt(v), hi=offsetAt(v+1);
		for(int slot=lo; slot<hi; slot++){
			counts[(int)((packedIDs[slot/3]>>>shift(slot))&MASK)]++;
		}
		return hi-lo;
	}

	/** Non-empty (value) cells in this bucket, for memory diagnostics. */
	long populatedCells(){
		long c=0;
		int prev=offsetAt(0);
		for(int v=1; v<=numValues; v++){
			final int cur=offsetAt(v);
			if(cur>prev){c++;}
			prev=cur;
		}
		return c;
	}

	/** Longs needed to hold n 21-bit fields, 3 per word. */
	static int words(int n){return (n+2)/3;}
	/** Bit offset of field i within its word: 0, 21, or 42. */
	private static int shift(int i){return BITS*(i%3);}

	static final int BITS=21;
	static final long MASK=(1L<<BITS)-1;//2097151 = 2^21-1

	private final int numValues;
	/** (numValues+1) cumulative offsets, 21-bit fields, 3 per word. */
	private final long[] packedOffsets;
	/** clade IDs grouped by value, 21-bit fields, 3 per word. */
	private final long[] packedIDs;
}
