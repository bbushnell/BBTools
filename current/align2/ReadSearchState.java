package align2;

import java.util.ArrayList;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Worker-owned shallow snapshot for one single-read search attempt.
 * Candidate lists, sites, arrays, and auxiliary objects are retained by reference;
 * capture does not clone their contents. Before a retry, BBMapThread restores
 * the entry state and detaches r.sites so the new search cannot modify the saved
 * low candidate list. Saved array contents must likewise remain unchanged.
 * Restores only the captured read, not fields in its mate or worker accounting.
 * Not a general Read rollback: insert and neuralMapqs are outside this snapshot.
 * Search does not update insert; BBMapS clears neural MAPQ before search when
 * inference is enabled and computes it after the selected attempt finishes.
 * Capture replaces an earlier snapshot; restore does not consume it. Clear
 * releases retained references at the attempt boundary. Not thread-safe.
 * @author Collei */
public final class ReadSearchState{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains one read's scalar values and object references without copying contents.
	 * @throws IllegalArgumentException If r is null; the prior snapshot is unchanged */
	public void capture(final Read r){
		if(r==null){throw new IllegalArgumentException("Cannot capture a null Read");}
		read=r; mate=r.mate; bases=r.bases; quality=r.quality; match=r.match; gaps=r.gaps;
		id=r.id; numericID=r.numericID; chrom=r.chrom; start=r.start; stop=r.stop;
		copies=r.copies; errors=r.errors; mapScore=r.mapScore; flags=r.flags; rand=r.rand;
		sites=r.sites; originalSite=r.originalSite; obj=r.obj; samline=r.samline;
	}

	/** Restores fields directly; Read.clearMapping/setInsert also mutate mate state.
	 * Does not copy saved objects or undo mutations to their contents.
	 * @throws IllegalStateException If capture has not succeeded since the last clear */
	public void restore(){
		if(read==null){throw new IllegalStateException("Capture must precede restore");}
		read.mate=mate; read.bases=bases; read.quality=quality; read.match=match; read.gaps=gaps;
		read.id=id; read.numericID=numericID; read.chrom=chrom; read.start=start; read.stop=stop;
		read.copies=copies; read.errors=errors; read.mapScore=mapScore; read.flags=flags; read.rand=rand;
		read.sites=sites; read.originalSite=originalSite; read.obj=obj; read.samline=samline;
	}

	/** Releases references and invalidates restore without modifying the captured read.
	 * Stale primitive values are harmless: the null read prevents their restoration. */
	public void clear(){
		read=null; mate=null; bases=null; quality=null; match=null; gaps=null; id=null;
		sites=null; originalSite=null; obj=null; samline=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private Read read, mate;
	private byte[] bases, quality, match;
	private int[] gaps;
	private String id;
	private long numericID;
	private int chrom, start, stop, copies, errors, mapScore, flags;
	private double rand;
	private ArrayList<SiteScore> sites;
	private SiteScore originalSite;
	private Object obj;
	private SamLine samline;
}
