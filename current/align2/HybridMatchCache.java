package align2;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;
import stream.SiteScore;

/** Two worker-owned native-match observations retained for one paired attempt.
 * Admission requires unchanged site state, no compressed-gap descriptor, and
 * no pre-existing match. The generated match may still contain short indels.
 * Only the inner genMatchStringForSite work is replaced; the outer wrapper must
 * execute to preserve score ordering and synchronization of the selected Read.
 * Not thread-safe: BBMapThread owns one cache, enables replay only during final
 * match generation, and clears it at attempt boundaries and in its finally block.
 * @author Collei */
final class HybridMatchCache {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	HybridMatchCache(){this(null);}
	/** The optional counter belongs to the same worker; null disables hit accounting. */
	HybridMatchCache(HybridMaxIndelStats stats_){stats=stats_;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Drops both observations and disables replay for the next attempt. */
	void clear(){entries[0]=entries[1]=null; replaying=false;}
	/** Drops one mate slot without changing whether the owner is replaying. */
	void clear(int slot){
		assert(slot>=0 && slot<2) : "Each paired observation owns one of two end slots";
		entries[slot]=null;
	}

	/** Replaces one slot only if native-match generation preserved the site's state.
	 * Both base orientations and the generated match are copied for later validation.
	 * @param slot Mate slot, 0 or 1
	 * @param id Read numeric ID, checked together with SiteScore object identity
	 * @param original Selected site before disposable match generation
	 * @param generated Disposable site after match generation
	 * @param plus Forward query bases, or null
	 * @param minus Reverse-complement query bases, or null
	 * @param imperfect Maximum imperfect score supplied to match generation
	 * @param maximum Maximum possible score supplied to match generation */
	void remember(int slot, long id, SiteScore original, SiteScore generated,
			byte[] plus, byte[] minus, int imperfect, int maximum){
		clear(slot);
		assert(!replaying) : "Observation cannot replace a cache during finalization";
		if(original==null || generated==null || original.match!=null || original.gaps!=null ||
				generated.gaps!=null || generated.match==null || generated.match.length==0 ||
				!sameState(original, generated)){return;}
		entries[slot]=new Entry(id, original, generated.match, plus, minus, imperfect, maximum);
	}

	/** Restores a private copy only during primary finalization with identical inputs.
	 * Once the site identity matches, the entry is consumed even if validation fails;
	 * a changed site cannot regain a stale observation later in the same attempt.
	 * @return True when site.match was populated; false leaves it untouched */
	boolean restore(long id, SiteScore site, byte[] plus, byte[] minus,
			int imperfect, int maximum, boolean secondary){
		if(!replaying || secondary){return false;}
		for(int i=0; i<entries.length; i++){
			final Entry e=entries[i];
			if(e==null || e.original!=site){continue;}
			entries[i]=null;// A retained observation is tried at most once.
			if(e.id!=id || site.match!=null || site.gaps!=null || e.imperfect!=imperfect ||
					e.maximum!=maximum || !sameState(e.before, site) ||
					!Arrays.equals(e.plus, plus) || !Arrays.equals(e.minus, minus)){return false;}
			site.match=e.match.clone();
			final AtomicLong audit=auditHits;
			if(audit!=null){audit.incrementAndGet();}
			if(stats!=null){stats.matchCacheHit();}
			return true;
		}
		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares the scalar candidate state used to admit and replay an observation.
	 * Array eligibility and query contents are checked separately by the caller. */
	private static boolean sameState(SiteScore a, SiteScore b){
		assert(a!=null && b!=null) : "Admission compares complete existing site states";
		return a.chrom==b.chrom && a.strand==b.strand && a.start==b.start && a.stop==b.stop &&
				a.quickScore==b.quickScore && a.score==b.score && a.slowScore==b.slowScore &&
				a.pairedScore==b.pairedScore && a.hits==b.hits && a.flags==b.flags &&
				a.rescued==b.rescued && a.perfect==b.perfect && a.semiperfect==b.semiperfect;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains site identity but owns all arrays used for validation and replay. */
	private static final class Entry {
		Entry(long id_, SiteScore site, byte[] result, byte[] plus_, byte[] minus_, int imperfect_, int maximum_){
			assert(site.match==null && site.gaps==null) : "Shallow scalar snapshot is safe only without arrays";
			id=id_; original=site; before=site.clone(); match=result.clone();
			plus=plus_==null ? null : plus_.clone(); minus=minus_==null ? null : minus_.clone();
			imperfect=imperfect_; maximum=maximum_;
		}
		final long id;
		final SiteScore original, before;
		final byte[] match, plus, minus;
		final int imperfect, maximum;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Enabled only around final match generation; always reset by the owner. */
	boolean replaying;
	/** Installed only by the functional test driver; ordinary runs leave it null. */
	static volatile AtomicLong auditHits;
	private final HybridMaxIndelStats stats;
	private final Entry[] entries=new Entry[2];
}
