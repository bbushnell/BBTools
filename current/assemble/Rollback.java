package assemble;

import stream.Read;
import structures.IntList;

/**
 * Saves and restores the state of a Read object for rollback operations, including ID, flags, bases, qualities, and optional counts.
 * Used to revert modifications during assembly or processing by cloning the original read state and restoring it later.
 * @author Brian Bushnell
 */
public class Rollback {

	/** Creates a rollback point for a read without count data by delegating to the full constructor.
	 * @param r The read to save state for */
	public Rollback(Read r){
		this(r, null);
	}

	/**
	 * Creates a rollback point for a read with optional count data.
	 * Saves the read's ID, flags, bases, quality scores, and count information.
	 * @param r The read to save state for
	 * @param counts Optional count data to preserve (may be null)
	 */
	public Rollback(Read r, IntList counts){
		id0=r.id;
		flags0=r.flags;
		bases0=r.bases.clone();
		quals0=(r.quality==null ? null : r.quality.clone());
		counts0=(counts==null ? null : counts.copy());
	}
	
	/** Restores a read to its saved state without count data by delegating to the full rollback method.
	 * @param r The read to restore */
	public void rollback(Read r){
		rollback(r, null);
	}
	
	/**
	 * Restores a read to its saved state including optional count data.
	 * Uses efficient array copying when sizes match, otherwise assigns new arrays.
	 * @param r The read to restore
	 * @param counts Optional count data to restore (may be null)
	 */
	public void rollback(Read r, IntList counts){
		r.id=id0;
		r.flags=flags0;
		if(r.length()==bases0.length){
			System.arraycopy(bases0, 0, r.bases, 0, bases0.length);
			if(quals0!=null){System.arraycopy(quals0, 0, r.quality, 0, quals0.length);}
			if(counts!=null){System.arraycopy(counts0.array, 0, counts.array, 0, counts0.size);}
		}else{
			r.bases=bases0;
			r.quality=quals0;
			if(counts!=null){
				counts.clear();
				counts.addAll(counts0);
			}
		}
	}
	
	/** Saved read identifier. */
	final String id0;
	/** Saved read flags. */
	final int flags0;
	/** Saved quality scores cloned from the original read (may be null). */
	final byte[] bases0, quals0;
	/** Saved count data (may be null). */
	public final IntList counts0;
	
}
