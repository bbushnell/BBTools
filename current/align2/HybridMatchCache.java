package align2;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;
import stream.SiteScore;

/** Two privately owned observation results for one hybrid attempt.
 * Experimental: only unchanged, gapless, initially unmatched sites qualify.
 * The complete genMatchString wrapper must still execute. @author Collei */
final class HybridMatchCache {

	void clear(){entries[0]=entries[1]=null;replaying=false;}
	void clear(int slot){
		assert(slot>=0 && slot<2) : "Each paired observation owns one of two end slots";
		entries[slot]=null;
	}
	void remember(int slot,long id,SiteScore original,SiteScore generated,
			byte[] plus,byte[] minus,int imperfect,int maximum){
		clear(slot);
		assert(!replaying) : "Observation cannot replace a cache during finalization";
		if(original==null || generated==null || original.match!=null || original.gaps!=null ||
				generated.gaps!=null || generated.match==null || generated.match.length==0 ||
				!sameState(original,generated)){return;}
		entries[slot]=new Entry(id,original,generated.match,plus,minus,imperfect,maximum);
	}
	boolean restore(long id,SiteScore site,byte[] plus,byte[] minus,
			int imperfect,int maximum,boolean secondary){
		if(!replaying || secondary){return false;}
		for(int i=0;i<entries.length;i++){
			final Entry e=entries[i];
			if(e==null || e.original!=site){continue;}
			entries[i]=null; // A retained observation is tried at most once.
			if(e.id!=id || site.match!=null || site.gaps!=null || e.imperfect!=imperfect ||
					e.maximum!=maximum || !sameState(e.before,site) ||
					!Arrays.equals(e.plus,plus) || !Arrays.equals(e.minus,minus)){return false;}
			site.match=e.match.clone();
			final AtomicLong audit=auditHits;
			if(audit!=null){audit.incrementAndGet();}
			return true;
		}
		return false;
	}
	private static boolean sameState(SiteScore a,SiteScore b){
		assert(a!=null && b!=null) : "Admission compares complete existing site states";
		return a.chrom==b.chrom && a.strand==b.strand && a.start==b.start && a.stop==b.stop &&
				a.quickScore==b.quickScore && a.score==b.score && a.slowScore==b.slowScore &&
				a.pairedScore==b.pairedScore && a.hits==b.hits && a.flags==b.flags &&
				a.rescued==b.rescued && a.perfect==b.perfect && a.semiperfect==b.semiperfect;
	}
	private static final class Entry {
		Entry(long id_,SiteScore site,byte[] result,byte[] plus_,byte[] minus_,int imperfect_,int maximum_){
			assert(site.match==null && site.gaps==null) : "Shallow scalar snapshot is safe only without arrays";
			id=id_;original=site;before=site.clone();match=result.clone();
			plus=plus_==null ? null : plus_.clone();minus=minus_==null ? null : minus_.clone();
			imperfect=imperfect_;maximum=maximum_;
		}
		final long id;
		final SiteScore original,before;
		final byte[] match,plus,minus;
		final int imperfect,maximum;
	}

	boolean replaying;
	/** Installed only by the functional test driver; ordinary runs leave it null. */
	static volatile AtomicLong auditHits;
	private final Entry[] entries=new Entry[2];
}
