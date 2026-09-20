package align2;

import java.util.ArrayList;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Reusable public-field snapshot for one mapping attempt.
 * Candidate lists and arrays are retained by reference; callers must detach
 * saved mutable structures before another attempt mutates them.
 * @author Collei
 */
public final class ReadSearchState {

	public void capture(Read r){
		if(r==null){throw new IllegalArgumentException("Cannot capture a null Read");}
		read=r;mate=r.mate;bases=r.bases;quality=r.quality;match=r.match;gaps=r.gaps;
		id=r.id;numericID=r.numericID;chrom=r.chrom;start=r.start;stop=r.stop;
		copies=r.copies;errors=r.errors;mapScore=r.mapScore;flags=r.flags;rand=r.rand;
		sites=r.sites;originalSite=r.originalSite;obj=r.obj;samline=r.samline;
	}

	/** Restore directly; Read setters can mutate mate state and are unsuitable here. */
	public void restore(){
		if(read==null){throw new IllegalStateException("Capture must precede restore");}
		read.mate=mate;read.bases=bases;read.quality=quality;read.match=match;read.gaps=gaps;
		read.id=id;read.numericID=numericID;read.chrom=chrom;read.start=start;read.stop=stop;
		read.copies=copies;read.errors=errors;read.mapScore=mapScore;read.flags=flags;read.rand=rand;
		read.sites=sites;read.originalSite=originalSite;read.obj=obj;read.samline=samline;
	}

	public void clear(){
		read=null;mate=null;bases=null;quality=null;match=null;gaps=null;id=null;
		sites=null;originalSite=null;obj=null;samline=null;
	}

	private Read read,mate;
	private byte[] bases,quality,match;
	private int[] gaps;
	private String id;
	private long numericID;
	private int chrom,start,stop,copies,errors,mapScore,flags;
	private double rand;
	private ArrayList<SiteScore> sites;
	private SiteScore originalSite;
	private Object obj;
	private SamLine samline;
}
