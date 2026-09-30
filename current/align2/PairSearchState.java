package align2;

import java.util.ArrayList;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Worker-owned shallow snapshot for a paired search attempt, before finalization.
 * Candidate lists and arrays are retained by reference, not copied: callers must
 * detach saved lists and keep saved object contents unchanged during another attempt.
 * BBMapThread restores entry state and detaches both site lists before a new
 * search. Its pseudo probe has a separate snapshot; the selected attempt alone
 * is finalized. Worker accounting and aligner scratch state are not captured.
 * Neither Read.insert nor neuralMapqs is saved. Search leaves insert unchanged;
 * BBMapS clears neural MAPQ before search when inference is enabled and fills it
 * after selection. This is not a general rollback facility for arbitrary Read
 * mutations. Restore is reusable until clear; neither method clones contents.
 * Not thread-safe, including the two-read restoration sequence.
 * @author Collei */
public final class PairSearchState{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Capture a reciprocal pair; reusing this holder replaces its previous snapshot. */
	public void capture(final Read read){
		assert(read!=null && read.mate!=null && read.mate.mate==read) : "Paired search requires reciprocal mate ownership";
		first.capture(read); second.capture(read.mate);
	}
	/** Restore both objects directly, avoiding setters that also mutate their mate. */
	public void restore(){
		assert(first.read!=null && second.read!=null) : "Capture must precede restore; a cleared snapshot cannot be restored";
		first.restore(); second.restore();
		assert(first.read.mate==second.read && second.read.mate==first.read) : "Restore must preserve captured pair ownership";
	}
	/** Release references after the attempt decision so a worker does not retain a prior pair. */
	public void clear(){first.clear(); second.clear();}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** One mate's scalar fields and retained references; private to this pair holder. */
	private static final class ReadState{
		/** Retains references; arrays and candidate objects remain owned by the caller. */
		void capture(final Read r){
			assert(r!=null) : "Cannot capture fields from a null Read";
			read=r; mate=r.mate; bases=r.bases; quality=r.quality; match=r.match; gaps=r.gaps;
			id=r.id; numericID=r.numericID; chrom=r.chrom; start=r.start; stop=r.stop;
			copies=r.copies; errors=r.errors; mapScore=r.mapScore; flags=r.flags; rand=r.rand;
			sites=r.sites; originalSite=r.originalSite; obj=r.obj; samline=r.samline;
		}
		/** Restores fields without Read setters that can also modify the other mate. */
		void restore(){
			assert(read!=null) : "A cleared holder has no restoration target";
			read.mate=mate; read.bases=bases; read.quality=quality; read.match=match; read.gaps=gaps;
			read.id=id; read.numericID=numericID; read.chrom=chrom; read.start=start; read.stop=stop;
			read.copies=copies; read.errors=errors; read.mapScore=mapScore; read.flags=flags; read.rand=rand;
			read.sites=sites; read.originalSite=originalSite; read.obj=obj; read.samline=samline;
			// Read.insert is private and is not changed by the search core. Do not call
			// setInsert here: it changes both mates and cannot preserve an invalid raw value.
		}
		/** Releases retained objects; stale scalars cannot be restored without a read. */
		void clear(){
			read=null; mate=null; bases=null; quality=null; match=null; gaps=null; id=null;
			sites=null; originalSite=null; obj=null; samline=null;
		}

		Read read, mate;
		byte[] bases, quality, match;
		int[] gaps;
		String id;
		long numericID;
		int chrom, start, stop, copies, errors, mapScore, flags;
		double rand;
		ArrayList<SiteScore> sites;
		SiteScore originalSite;
		Object obj;
		SamLine samline;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ReadState first=new ReadState(), second=new ReadState();
}
