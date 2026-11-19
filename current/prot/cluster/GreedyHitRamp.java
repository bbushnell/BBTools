package prot.cluster;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

/** Optional, monotonic greedy hit-floor ramp over actually published families.
 * One helper belongs to one greedy pass. The driver captures early() once per
 * query and reuses that boolean for its shortlist and competitive seed checks.
 * Model replacement never creates another qualifying stable family ID.
 * @author Brian Bushnell, Fischl
 */
public final class GreedyHitRamp{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Capacity bounds stable family IDs; the activation target must be explicitly positive.
	 * A target larger than capacity is valid and simply cannot be reached in this pass.
	 * The driver constructs this helper only when ramp mode is enabled.
	 */
	public GreedyHitRamp(final int capacity_, final int targetFamilies_, final int minimumMembers_){
		if(capacity_<0 || targetFamilies_<=0 || minimumMembers_<=0){
			throw new IllegalArgumentException("Greedy ramp requires nonnegative capacity and explicit positive family/member thresholds");
		}
		capacity=capacity_; targetFamilies=targetFamilies_; minimumMembers=minimumMembers_;
		marked=minimumMembers>1 ? new AtomicIntegerArray(capacity) : null;
		qualified=minimumMembers>1 ? new AtomicInteger() : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Admission / Capture    ----------------*/
	/*--------------------------------------------------------------*/

	/** Notify only after actual seed/join publication, using the published membership size.
	 * For minimumMembers=1 the published family count is sufficient, so this is a no-op.
	 * A winning mark immediately increments the count without intervening callbacks or
	 * allocation. A preempted marker can delay activation, but cannot cause early activation.
	 */
	public void admitted(final int familyId, final int publishedMembers){
		if(minimumMembers==1){return;}
		if(familyId<0 || familyId>=capacity || publishedMembers<0){
			throw new IllegalArgumentException("Invalid published family admission: id="+familyId+" size="+publishedMembers);
		}
		if(publishedMembers<minimumMembers || marked.get(familyId)!=0){return;}
		if(marked.compareAndSet(familyId,0,1)){
			final int count=qualified.incrementAndGet();
			assert(count>0 && count<=capacity) : "Each bounded stable family ID contributes at most one qualification";
		}
	}

	/** Capture once per query: true means the lower floor, false means the original floor.
	 * Once any query observes the threshold, later captures cannot return to the lower
	 * floor even if their supplied family-count snapshot was taken earlier.
	 */
	public boolean early(final int publishedFamilyCount){
		validatePublished(publishedFamilyCount);
		if(raised.get()){return false;}
		final int count=minimumMembers==1 ? publishedFamilyCount : qualified.get();
		assert(count>=0 && count<=capacity) : "Published families and once-per-ID qualifications cannot exceed stable-ID capacity";
		if(count>=targetFamilies){raised.compareAndSet(false,true);}
		return !raised.get();
	}

	/** A source snapshot may be old, but cannot be outside the allocated stable-ID population. */
	private void validatePublished(final int publishedFamilyCount){
		if(publishedFamilyCount<0 || publishedFamilyCount>capacity){
			throw new IllegalArgumentException("Published family count outside ramp capacity: "+publishedFamilyCount);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int capacity, targetFamilies, minimumMembers;
	private final AtomicIntegerArray marked;
	private final AtomicInteger qualified;
	private final AtomicBoolean raised=new AtomicBoolean();
}
