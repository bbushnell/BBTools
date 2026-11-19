package prot.cluster;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;

import ddl.MutableFamilyIndex;

/** Atomic family versions with a contiguous, monotonically published ID prefix.
 * A version is initialized completely before publication and never mutated afterward,
 * except through its shared append-only membership container. No slot is removed
 * during greedy admission. Pruning belongs after all query workers have joined.
 * @author Keqing
 */
public final class FamilySlots implements MutableFamilyIndex.SketchSource {

	/** At most one family can be seeded per eligible input. */
	public FamilySlots(final int capacity){
		if(capacity<0){throw new IllegalArgumentException("Negative eligible-input capacity: "+capacity);}
		slots=new AtomicReferenceArray<Family>(capacity);
	}

	/** Acquire-visible count: every slot below this boundary contains a complete version. */
	@Override public int size(){return published.get();}

	/** Returns the currently published version, including a just-seeded boundary slot.
	 * Only seed competitors inspect the boundary; rankers use IDs strictly below size().
	 */
	public Family get(final int id){
		if(id<0 || id>=slots.length()){throw new IllegalArgumentException("Family ID outside capacity: "+id);}
		return slots.get(id);
	}

	/** Borrows an immutable sketch for the snapshot builder or one ranking query. */
	@Override public char[] sketch(final int id){
		if(id<0 || id>=size()){throw new IllegalArgumentException("Sketch ID outside the published prefix: "+id);}
		final Family family=slots.get(id);
		assert(family!=null) : "Published prefix cannot contain holes: snapshot ranking reads every ID below size()";
		return family.sketchUnsafe();
	}

	/** Competes for the next empty slot, without reserving a hole in the visible prefix.
	 * A false result means the caller must compare its query with get(id), then join
	 * that family if admissible or proceed to the next slot. Both winners and losers
	 * help publish initialized slots, so a delayed winner cannot strand the prefix.
	 */
	public boolean seed(final int id, final Family candidate){
		validate(id,candidate);
		if(id>size()){
			// GreedyWorker can observe get(id-1) after its slot CAS but before
			// that winner advances published. Help the initialized prefix rather
			// than mistaking this legal observation for a hole or waiting on it.
			publishPrefix();
			if(id>size()){
				throw new IllegalStateException("Cannot skip unpublished family slots: id="+id+" published="+size());
			}
		}
		final boolean won=slots.compareAndSet(id,null,candidate);
		publishPrefix();
		return won;
	}

	/** Publishes a regenerated version only if the version it was built from is current.
	 * Membership ownership must be shared between expected and replacement; neither
	 * the losing candidate nor any published model may be mutated after this call.
	 */
	public boolean replace(final Family expected, final Family replacement){
		if(expected==null){throw new IllegalArgumentException("Missing previous family version");}
		validate(expected.familyId,replacement);
		if(expected.familyId!=replacement.familyId || expected.representativeIndex!=replacement.representativeIndex ||
			!expected.id.equals(replacement.id) || !expected.sharesMembers(replacement)){
			throw new IllegalArgumentException("Regeneration must preserve family identity and shared membership ownership");
		}
		return slots.compareAndSet(expected.familyId,expected,replacement);
	}

	/** Captures versions after the owner has joined all workers and the snapshot builder. */
	public Family[] completed(){
		publishPrefix();
		final int count=size();
		final Family[] result=new Family[count];
		for(int id=0; id<count; id++){
			result[id]=slots.get(id);
			assert(result[id]!=null && result[id].familyId==id) :
				"Frozen-pass handoff requires dense stable family IDs; missing/mismatched id="+id;
		}
		return result;
	}

	private void validate(final int id, final Family family){
		if(id<0 || id>=slots.length() || family==null || family.familyId!=id ||
			family.isEmpty() || family.sketchUnsafe()==null){
			throw new IllegalArgumentException("Publish a nonempty, fully initialized family at its own ID: "+id);
		}
	}

	/** Each successful count CAS extends only across an already initialized slot. */
	private void publishPrefix(){
		for(;;){
			final int count=published.get();
			if(count==slots.length() || slots.get(count)==null){return;}
			published.compareAndSet(count,count+1);
		}
	}

	private final AtomicReferenceArray<Family> slots;
	private final AtomicInteger published=new AtomicInteger();
}
