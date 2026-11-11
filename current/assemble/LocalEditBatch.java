package assemble;

import shared.Shared;
import stream.Read;
import structures.LongList;

/** Worker-reused original-coordinate edit list with one-pass materialization.
 * This class does not discover or validate biological corrections. Its caller
 * must validate complete dependency spans before recording independent edits.
 * Input arrays must remain immutable until commit or clear. Reference checks detect
 * replacement arrays, not in-place writes. One instance and its Read belong to one
 * worker throughout the transaction; installation is not a concurrent publication API.
 * @author Fischl */
public final class LocalEditBatch{

	/**
	 * Abandons any previous list, resets copiedBases and borrows an editable unpaired,
	 * alignment-free Read with matching or null qualities. Neither input array is copied.
	 * A failed reset leaves no active transaction and does not change either Read.
	 */
	public void reset(final Read input){
		clear(); copiedBases=0;
		HomopolymerIndelEdit.requireEditable(input);
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.bases.length)){
			throw new IllegalArgumentException("Edit batch requires bases and matching or null qualities.");
		}
		read=input; bases=input.bases; quality=input.quality;
	}
	/** Replaces one original base using the numeric allele code A/C/G/T=0/1/2/3. */
	public void addSubstitution(final int position, final int base){append(SUB, position, base);}
	/** Deletes one original base; recording does not yet mutate the Read. */
	public void addDeletion(final int position){append(DEL, position, 0);}
	/** Inserts a numeric A/C/G/T allele before an original base; position=length appends. */
	public void addInsertion(final int position, final int base){append(INS, position, base);}
	/** Appends only strictly increasing original coordinates, rejecting coincident edits
	 * before changing the list. Ordering alone does not prove dependency independence. */
	private void append(final int operation, final int position, final int base){
		if(read==null){throw new IllegalStateException("Reset the batch before recording an edit.");}
		if(position<0 || (operation==INS ? position>bases.length : position>=bases.length)){
			throw new IllegalArgumentException("Edit position outside original operation domain: "+position);
		}
		if(base<0 || base>3){throw new IllegalArgumentException("Replacement base code must be A/C/G/T=0/1/2/3: "+base);}
		if(edits.size>0 && position<=position(edits.array[edits.size-1])){
			throw new IllegalArgumentException("Original-coordinate edits must be strictly increasing; same-position edits require joint resolution.");
		}
		assert(operation>=SUB && operation<=INS) : "Only the three explicit edit methods may create packed operation codes.";
		// 31 position bits, two operation bits, two base bits. Cast BEFORE shifting.
		edits.add(((long)position<<4)|((long)operation<<2)|base);
		delta+=operation==INS ? 1 : operation==DEL ? -1 : 0;
	}
	/**
	 * Materializes all recorded operations, returning their count or -1 when over budget.
	 * The nonnegative budget counts records, even a substitution to the existing allele.
	 * An empty batch returns zero without replacing arrays. Otherwise installs new bases
	 * and optional qualities only after copying; old arrays remain untouched. Retained
	 * qualities are copied, substituted/inserted qualities are zero, and null stays null.
	 * Validation/allocation failures do not install replacements, including when another
	 * owner has already replaced the input arrays. Always closes an active batch; failed
	 * attempts can leave a partial copiedBases counter but cannot leave pending edits.
	 */
	public int commit(final int maxEdits){
		if(read==null){throw new IllegalStateException("No active edit batch.");}
		copiedBases=0;
		try{
			if(maxEdits<0){throw new IllegalArgumentException("Batch edit budget cannot be negative.");}
			HomopolymerIndelEdit.requireEditable(read);
			if(read.bases!=bases || read.quality!=quality){throw new IllegalStateException("Read arrays changed during the original-coordinate transaction.");}
			final int count=edits.size;
			if(count>maxEdits){return -1;}
			if(count==0){return 0;}
			final int length=resultLength(bases.length, delta);
			final byte[] output=new byte[length], outputQuality=quality==null ? null : new byte[length];
			int sourceEnd=bases.length, destinationEnd=length;
			for(int i=count-1; i>=0; i--){
				final long entry=edits.array[i]; final int position=position(entry), op=(int)((entry>>>2)&3);
				final int spanStart=position+(op==INS ? 0 : 1), spanLength=sourceEnd-spanStart;
				assert(op<=INS && spanLength>=0) : "Strictly ordered original coordinates partition source spans without overlap.";
				destinationEnd-=spanLength;
				System.arraycopy(bases, spanStart, output, destinationEnd, spanLength);
				if(quality!=null){System.arraycopy(quality, spanStart, outputQuality, destinationEnd, spanLength);}
				copiedBases+=spanLength;
				if(op!=DEL){
					output[--destinationEnd]=ALPHABET[(int)(entry&3)];
					// New arrays are zero-filled: substituted/inserted quality stays Q0.
				}
				sourceEnd=position;
			}
			destinationEnd-=sourceEnd;
			assert(destinationEnd==0) : "Net indel delta and backward span traversal must fill exactly the final array length.";
			System.arraycopy(bases, 0, output, 0, sourceEnd);
			if(quality!=null){System.arraycopy(quality, 0, outputQuality, 0, sourceEnd);}
			copiedBases+=sourceEnd;
			assert(copiedBases<=bases.length && (outputQuality==null || outputQuality.length==output.length)) :
				"Every retained original base is copied at most once; every output base has one quality when present.";
			read.bases=output; read.quality=outputQuality;
			return count;
		}finally{clear();}
	}
	/** Drops proposals/references while retaining the primitive buffer and copiedBases statistic. */
	public void clear(){edits.clear(); read=null; bases=null; quality=null; delta=0;}
	/** Returns pending operation records, not a count of independently validated corrections. */
	public int size(){return edits.size;}
	/** Checks nonnegative original length and a final length within Shared.MAX_ARRAY_LEN,
	 * including deletion to empty, without allocating an array or overflowing the sum. */
	public static int resultLength(final int originalLength, final long delta){
		if(originalLength<0 || delta<-(long)originalLength || delta>(long)Shared.MAX_ARRAY_LEN-originalLength){
			throw new IllegalArgumentException("Edit batch result exceeds supported array bounds: original="+originalLength+" delta="+delta);
		}
		return (int)(originalLength+delta);
	}
	/** Decodes the nonnegative original position from a packed edit. */
	private static int position(final long entry){
		assert((entry>>>35)==0) : "A packed edit uses exactly 31 nonnegative position bits and four operation/base bits.";
		return (int)(entry>>>4);
	}
	/** Original bases copied by the last commit attempt, excluding supplied substitution/
	 * insertion bases and quality copies; reset clears it, while clear preserves it. */
	public long copiedBases;
	/** Borrowed transaction target; null when no list can be committed. */
	private Read read;
	/** Original array identities, checked again before materialization. */
	private byte[] bases, quality;
	/** Net length change from recorded insertions minus deletions. */
	private long delta;
	/** Ordered packed position/operation/base records; capacity survives transaction closure. */
	private final LongList edits=new LongList(16);
	/** Two-bit operation codes packed above each two-bit replacement allele. */
	private static final int SUB=0, DEL=1, INS=2;
	/** Decodes the numeric allele accepted by substitution/insertion methods. */
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
