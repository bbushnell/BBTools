package structures;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import shared.Shared;
import stream.HasID;
import stream.Read;

/**
 * Numbered list wrapper for multithreaded producer-consumer patterns.
 * 
 * Wraps an ArrayList with a sequential ID to enable ordered processing in multithreaded
 * pipelines. Workers can process lists out of order while consumers reconstruct the
 * original sequence using the ID field.
 * 
 * Supports poison pill and last-job signaling for clean pipeline shutdown.
 * Optionally generates deterministic random numbers for Read objects to enable
 * reproducible subsampling in multithreaded contexts.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 2011
 * 
 * @param <K> Element type (must be Serializable)
 */
public final class ListNum<K extends Serializable> implements Serializable, Iterable<K>, HasID {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private static final long serialVersionUID=-7509242172010729386L;

	/** Standard constructor for data lists. Creates a numbered list with no special flags. */
	public ListNum(ArrayList<K> list_, long id_){
		this(list_, id_, NORMAL);
	}
	
	/**
	 * Full constructor with explicit poison and last flags.
	 * 
	 * @param list_ ArrayList to wrap (null allowed for poison pills)
	 * @param id_ Sequential identifier for ordering
	 * @param poison_ True if this is a poison pill signaling worker shutdown
	 * @param last_ True if this is the last job in the sequence
	 */
	public ListNum(ArrayList<K> list_, long id_, boolean poison_, boolean last_){
		this(list_, id_, last_ ? LAST : poison_ ? POISON : NORMAL);
	}
	
	/**
	 * Full constructor with explicit type flag.
	 * 
	 * @param list_ ArrayList to wrap (null allowed for poison pills)
	 * @param id_ Sequential identifier for ordering
	 * @param type_ Job type
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
	
	/** Returns number of elements in the list, or 0 if list is null. */
	public final int size(){return list==null ? 0 : list.size();}
	
	@Override
	public String toString(){return list==null ? "ln.list=null" : list.toString();}
	
	/** Returns true if list is null or empty. */
	public final boolean isEmpty(){return list==null || list.isEmpty();}

	/** Returns element at specified index. */
	public final K get(int i){return list.get(i);}
	
	/** Replaces element at specified index. */
	public final K set(int i, K k){return list.set(i, k);}
	
	/** Removes element at specified index. */
	public final K remove(int i){return list.remove(i);}
	
	/** Appends element to list. */
	public final void add(K k){list.add(k);}
	
	/** Removes all elements from list. */
	public final void clear(){list.clear();}	
	
	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Iterator<K> iterator(){return list==null ? null : list.iterator();}
	
	/** Returns sequential identifier for ordering. */
	@Override
	public final long id(){return id;}
	
	/** Returns true if this is a poison pill signaling worker shutdown. */
	@Override
	public final boolean poison(){return type==POISON;}
	
	/** Returns true if this is the last job in the sequence. */
	@Override
	public final boolean last(){return type==LAST;}
	
	/** Any terminal type */
	public final boolean finished(){return type>=LAST;}
	
	@Override
	public ListNum<K> makePoison(long id_) {return new ListNum<K>(null, id_, POISON);}
	
	@Override
	public ListNum<K> makeLast(long id_){return new ListNum<K>(null, id_, LAST);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Random            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Set the seed for deterministic random number generation.
	 * If seed is negative, uses current time plus random offset.
	 */
	public static synchronized void setDeterministicRandomSeed(long seed_){
		if(seed_>=0){seed=seed_;}
		else{seed=System.nanoTime()+(long)(Math.random()*10000000);}
	}
	
	/**
	 * Enable or disable deterministic random number generation for Read objects.
	 * When enabled, each Read in a ListNum gets a deterministic random value,
	 * allowing reproducible subsampling in multithreaded contexts.
	 */
	public static synchronized void setDeterministicRandom(boolean b){
		GEN_RANDOM_NUMBERS=b;
		if(b){
			randy=Shared.threadLocalRandom(seed);
			seed++;
		}
	}
	
	/** Returns true if deterministic random number generation is enabled. */
	public static boolean deterministicRandom(){return GEN_RANDOM_NUMBERS;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Wrapped list (may be null for poison pills). */
	public final ArrayList<K> list;
	/** Sequential identifier for ordering in multithreaded pipelines. */
	public final long id;
	/** Job type. */
	public final int type;
	
	public long firstRecordNum=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Enable deterministic random number generation for reproducible subsampling */
	private static boolean GEN_RANDOM_NUMBERS=false;
	/** Random number generator for deterministic mode */
	private static Random randy;
	/** Seed for deterministic random number generation */
	private static long seed=0;
	
	/** Prototype, for reflection. */
	public static final int PROTO=-1;
	/** Normal job. */
	public static final int NORMAL=0;
	/** True if this is the last job in the sequence. */
	public static final int LAST=3;
	/** True if this is a poison pill signaling worker shutdown. */
	public static final int POISON=4;
	
}