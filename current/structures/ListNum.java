package structures;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import shared.Random;

import shared.Shared;
import stream.HasID;
import stream.Read;

/**
 * Numbered list wrapper for multithreaded producer-consumer pipelines.
 * Borrows an ArrayList without copying it and carries a caller-supplied ID and
 * job type. This wrapper neither assigns sequential IDs nor synchronizes access
 * to the mutable list. Optional random-number assignment modifies each non-null
 * element whose exact class is Read, using one shared generator.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 2011
 */
public final class ListNum<K extends Serializable> implements Serializable, Iterable<K>, HasID{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private static final long serialVersionUID=-7509242172010729386L;

	/**
	 * Creates a numbered list with default job type.
	 * @param list_ Non-null ArrayList to wrap for this normal job.
	 * @param id_ Caller-supplied identifier
	 */
	public ListNum(ArrayList<K> list_, long id_){
		this(list_, id_, NORMAL);
	}

	/**
	 * Creates a numbered list with explicit poison/last flags.
	 * If both flags are true, the last-job type takes precedence.
	 *
	 * @param list_ ArrayList to wrap (null allowed for poison/last)
	 * @param id_ Caller-supplied identifier
	 * @param poison_ True if this is a poison pill
	 * @param last_ True if this is the last job
	 */
	public ListNum(ArrayList<K> list_, long id_, boolean poison_, boolean last_){
		this(list_, id_, last_ ? LAST : poison_ ? POISON : NORMAL);
	}

	/**
	 * Creates a numbered list with explicit job type.
	 * Optionally assigns random numbers to exact Read instances when enabled.
	 * Normal jobs require a non-null list, poison jobs require a null list, and
	 * last/prototype jobs may have either; these preconditions are assertion-checked.
	 *
	 * @param list_ ArrayList to wrap (null allowed for poison/last/proto)
	 * @param id_ Caller-supplied identifier
	 * @param type_ Job type flag (NORMAL/POISON/LAST/PROTO)
	 */
	public ListNum(ArrayList<K> list_, long id_, int type_){
		list=list_;
		id=id_;
		type=type_;
		if(GEN_RANDOM_NUMBERS && list!=null){
			for(K k : list){
				if(k!=null && k.getClass()==Read.class){
					((Read)k).rand=randy.nextDouble();
				}
			}
		}
		assert(list!=null || (poison() || last() || type==PROTO)); //Regular jobs may not be null (they can be empty though)
		assert(!poison() || list==null); //Poison should have a null list
		assert(!(poison() && last())); //There can only be one last but multiple poison
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the number of elements, or 0 if list is null.
	 * @return Element count */
	public final int size(){return list==null ? 0 : list.size();}

	/** Returns "ln.list=null" if list is null, otherwise list.toString(). */
	@Override
	public String toString(){return list==null ? "ln.list=null" : list.toString();}

	/** Returns true if the list is null or empty.
	 * @return true if empty */
	public final boolean isEmpty(){return list==null || list.isEmpty();}

	/**
	 * Returns the element at index i.
	 * @param i Index
	 * @return Element at index
	 */
	public final K get(int i){return list.get(i);}

	/**
	 * Replaces the element at index i.
	 * @param i Index to replace
	 * @param k New value
	 * @return Previous value
	 */
	public final K set(int i, K k){return list.set(i, k);}

	/**
	 * Removes and returns the element at index i.
	 * @param i Index to remove
	 * @return Removed value
	 */
	public final K remove(int i){return list.remove(i);}

	/** Appends an element to the list.
	 * @param k Element to add */
	public final void add(K k){list.add(k);}

	/** Removes all elements from the list. */
	public final void clear(){list.clear();}

	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns an iterator over the wrapped list (null if list is null).
	 * @return Iterator or null */
	@Override
	public Iterator<K> iterator(){return list==null ? null : list.iterator();}

	/** Returns the caller-supplied identifier for this list.
	 * @return ID value */
	@Override
	public final long id(){return id;}

	/** Returns true if this job is a poison pill.
	 * @return true if poison */
	@Override
	public final boolean poison(){return type==POISON;}

	/** Returns true if this job is marked as the last in sequence.
	 * @return true if last */
	@Override
	public final boolean last(){return type==LAST;}

	/** Classifies types at or above LAST as terminal; does not track processing completion.
	 * @return true if terminal */
	public final boolean finished(){return type>=LAST;}

	/**
	 * Creates a poison-pill ListNum with the given ID.
	 * @param id_ ID to assign
	 * @return New poison ListNum
	 */
	@Override
	public ListNum<K> makePoison(long id_){return new ListNum<K>(null, id_, POISON);}

	/**
	 * Creates a last-job ListNum with the given ID.
	 * @param id_ ID to assign
	 * @return New last-job ListNum
	 */
	@Override
	public ListNum<K> makeLast(long id_){return new ListNum<K>(null, id_, LAST);}

	/*--------------------------------------------------------------*/
	/*----------------            Random            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sets the seed used by the next call that enables random-number assignment.
	 * Does not reseed an already initialized generator.
	 * Negative seeds use system time plus a random offset.
	 * @param seed_ Seed value or negative for time-based seed
	 */
	public static synchronized void setDeterministicRandomSeed(long seed_){
		if(seed_>=0){seed=seed_;}
		else{seed=System.nanoTime()+(long)(Math.random()*10000000);}
	}

	/**
	 * Enables or disables random-number assignment for exact Read instances.
	 * Each enabling call replaces the shared generator using the current seed,
	 * then advances that seed. Configure this mode before constructing batches;
	 * reproducible assignment also depends on the order of generator calls.
	 * @param b true to enable deterministic mode
	 */
	public static synchronized void setDeterministicRandom(boolean b){
		GEN_RANDOM_NUMBERS=b;
		if(b){
			randy=Shared.threadLocalRandom(seed);
			seed++;
		}
	}

	/** Returns whether deterministic random number generation is enabled.
	 * @return true if deterministic mode is on */
	public static boolean deterministicRandom(){return GEN_RANDOM_NUMBERS;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed mutable list; null is permitted for poison, last, and prototype jobs. */
	public final ArrayList<K> list;
	/** Caller-supplied identifier, used for ordering by consumers that require it. */
	public final long id;
	/** Job type flag (NORMAL, POISON, LAST, PROTO). */
	public final int type;

	/** Optional first record number carried with this list. */
	public long firstRecordNum=-1;

	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/

	/** Controls whether Reads receive deterministic random numbers. */
	private static boolean GEN_RANDOM_NUMBERS=false;
	/** Shared random generator used by constructors while random assignment is enabled. */
	private static Random randy;
	/** Seed used to initialize deterministic random generator. */
	private static long seed=0;

	/** Prototype job type constant. */
	public static final int PROTO=-1;
	/** Normal job type constant. */
	public static final int NORMAL=0;
	/** Last-job type constant. */
	public static final int LAST=3;
	/** Poison-pill job type constant. */
	public static final int POISON=4;

}
