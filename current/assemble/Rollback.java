package assemble;

import stream.Read;
import structures.IntList;

/**
 * One-use snapshot of a read's string ID, flags, bases, qualities and optional counts.
 * Other read state, including numericID, mate and alignment fields, is not saved.
 * Correction callers own the read and counts exclusively while saving/restoring.
 * A restore may transfer snapshot arrays to the read; discard the snapshot afterward.
 * @author Brian Bushnell
 */
public class Rollback{

	/** Copies the saved read fields without a count profile.
	 * @param r Read with non-null bases and optional qualities */
	public Rollback(final Read r){
		this(r, null);
	}

	/**
	 * Copies bases, optional qualities and the count list's live prefix.
	 * The immutable string ID is retained; flags are copied by value.
	 * @param r Read with non-null bases and optional qualities
	 * @param counts Optional count profile; its capacity and unused tail are not saved
	 */
	public Rollback(final Read r, final IntList counts){
		assert(r!=null && r.bases!=null) : "Rollback needs read bases to snapshot the correction input.";
		id0=r.id;
		flags0=r.flags;
		bases0=r.bases.clone();
		quals0=(r.quality==null ? null : r.quality.clone());
		counts0=(counts==null ? null : counts.copy());
	}
	
	/** Restores saved read fields once, leaving any separately held counts unchanged.
	 * @param r Non-null read owned exclusively by the caller */
	public void rollback(final Read r){
		rollback(r, null);
	}
	
	/**
	 * Restores saved fields once, independently of current array lengths or count size.
	 * Reuses equal-length base/quality arrays and sufficient count capacity; otherwise
	 * transfers saved base/quality arrays or grows the destination count list.
	 * Saved null qualities clear any added qualities. Counts regain their saved size
	 * and values; unused capacity is not cleared. No other read fields are restored.
	 * @param r Non-null read owned exclusively by the caller
	 * @param counts Destination count list, or null to omit count restoration
	 * @throws IllegalArgumentException If counts are requested but none were saved
	 */
	public void rollback(final Read r, final IntList counts){
		assert(r!=null) : "Rollback needs a destination read for the saved correction input.";
		if(counts!=null && counts0==null){
			throw new IllegalArgumentException("Cannot restore counts: this Rollback saved no count profile.");
		}
		r.id=id0;
		r.flags=flags0;
		if(r.bases!=null && r.bases.length==bases0.length){
			System.arraycopy(bases0, 0, r.bases, 0, bases0.length);
		}else{r.bases=bases0;}
		if(quals0!=null && r.quality!=null && r.quality.length==quals0.length){
			System.arraycopy(quals0, 0, r.quality, 0, quals0.length);
		}else{r.quality=quals0;}
		if(counts!=null){
			if(counts.size!=counts0.size){counts.setSize(counts0.size);}
			System.arraycopy(counts0.array, 0, counts.array, 0, counts0.size);
		}
	}
	
	/** Saved read identifier. */
	final String id0;
	/** Saved read flags. */
	final int flags0;
	/** Saved bases, transferred to the read when its current storage cannot be reused. */
	final byte[] bases0;
	/** Saved optional qualities, likewise transferable on restore. */
	final byte[] quals0;
	/** Copied live count profile, or null. Exposed for correction comparisons;
	 * callers must not mutate it before restoration. */
	public final IntList counts0;
	
}
