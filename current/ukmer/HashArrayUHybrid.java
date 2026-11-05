package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.IntList2;

/**
 * Stores long k-mer keys with positive scalar values or arrays of multiple values.
 * Primary cells hold one value directly or a negated index into an IntList2;
 * overflow entries use a two-dimensional HashForestU. Input arrays can be retained
 * and returned as live storage. Callers coordinate mutation and exclude resizing
 * during traversal/ownership use. Counting methods require primary entries to
 * remain direct scalars; victim HASH_COLLISIONs still delegate to the unresolved
 * UKM109 two-dimensional forest increment defect. Value-list promotion is for
 * set()/set(values) paths.
 *
 * @author Brian Bushnell
 * @date Oct 25, 2013
 */
public final class HashArrayUHybrid extends HashArrayU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the borrowed resize schedule; reserves list indices zero/one to keep tags below -1. */
	public HashArrayUHybrid(int[] schedule_, int k_, int kbig_){
		super(schedule_, k_, kbig_, true);
		values=allocInt1D(prime+extra);
		setList=new IntList2();
		setList.add(null);
		setList.add(null);
	}

//	public HashArrayUHybrid(int initialSize, int k_, int kbig_, boolean autoResize_){
//		super(initialSize, k_, kbig_, autoResize_, true);
//		values=allocInt1D(prime+extra);
//		setList=new IntList2();
//		setList.add(null);
//		setList.add(null);
//	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments the count for a k-mer, inserting if not present.
	 * Existing primary entries must still be direct scalars, not list tags.
	 * Victim hash collisions delegate to a two-dimensional HashForestU; repeated
	 * increments on an existing victim try to append its first value plus one,
	 * while retaining and returning that existing first value until UKM109 is fixed.
	 *
	 * @param kmer The k-mer to increment
	 * @return updated primary count, or unchanged first victim value under UKM109
	 */
	@Override
	public final int increment(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);

		if(cell==HASH_COLLISION){
			int x=victims.increment(kmer);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}else if(cellEmpty(cell)){
			setKmer(kmer.key(), cell);
			size++;
			values[cell]=1;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			values[cell]++;
			if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
			return values[cell];
		}
	}

	/**
	 * Increments a k-mer count and returns the number of new k-mers created.
	 * Returns 1 if the k-mer was newly inserted, 0 if it already existed.
	 *
	 * Existing primary entries must still be direct scalars, not list tags.
	 * Victim hash collisions delegate to a two-dimensional HashForestU; existing
	 * victim nodes run the unresolved UKM109 list-backed increment path before
	 * returning 0.
	 *
	 * @param kmer The k-mer to increment
	 * @return 1 if k-mer was newly created, 0 if it already existed
	 */
	@Override
	public final int incrementAndReturnNumCreated(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);

		if(cell==HASH_COLLISION){
			int x=victims.incrementAndReturnNumCreated(kmer);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}else if(cellEmpty(cell)){
			setKmer(kmer.key(), cell);
			size++;
			values[cell]=1;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			values[cell]++;
			if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
			return 0;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads the first value from a cell, following a list reference when needed.
	 * Encodings {@code >= -1} are returned directly; encodings {@code <= -2}
	 * are negated indices into setList, whose first element is returned.
	 *
	 * @param cell Primary array index to read from
	 * @return First list value or raw direct value; unused primary slots contain zero
	 */
	@Override
	protected final int readCellValue(int cell){
		final int x=values[cell];
		if(x>-2){return x;}
		return setList.get(0-x)[0];
	}

	/**
	 * Reads all values stored in a hash table cell.
	 * For direct cells, fills the provided singleton array. For list tags,
	 * returns the live setList array, including spare NOT_PRESENT slots.
	 * Callers must not mutate retained arrays independently.
	 *
	 * @param cell The cell index to read from
	 * @param singleton Single-element array for single-value returns
	 * @return Array containing all values for the cell
	 */
	@Override
	protected final int[] readCellValues(int cell, int[] singleton){
		final int x=values[cell];
		if(x>-2){
			singleton[0]=values[cell];
			return singleton;
		}
		return setList.get(0-x);
	}

	/** Returns true when no inline value or setList index is stored at this primary cell. */
	@Override
	protected final boolean cellEmpty(int cell){return values[cell]==0;}

	/**
	 * Inserts multiple values into a hash table cell for a k-mer.
	 * A one-element input uses scalar insertion. Existing lists merge until
	 * NOT_PRESENT; a new list retains vals directly or copies it with an
	 * existing scalar prepended. Retained input must remain stable.
	 *
	 * @param kmer The k-mer key being inserted
	 * @param vals Nonempty positive-value prefix with NOT_PRESENT padding, or a full positive array
	 * @param cell The target cell for insertion
	 */
	@Override
	protected final void insertValue(long[] kmer, int[] vals, int cell){
		if(verbose){System.err.println("insertValue("+kmer+", "+Arrays.toString(vals)+", "+cell+"); old="+values[cell]);}
		assert(matches(kmer, cell));
		if(vals.length==1){
			if(verbose){System.err.println("A: length=1");}
			insertValue(kmer, vals[0], cell);
			return;
		}
		final int old=values[cell];
		if(old==vals[0] && vals[1]==NOT_PRESENT){
			if(verbose){System.err.println("B: old==vals[0] && vals[1]==-1");}
			return; //Nothing to do
		}else if(old<-1){//An array already exists
			if(verbose){System.err.println("C: old<-1");}
			for(int i : vals){
				if(i==-1){break;}
				insertIntoList(i, -old);
			}
		}else{//Add the list
			final int[] temp;
			if(old>0){//Move the old value to a new array.  Note that this will probably never be used.
				if(verbose){System.err.println("D: old>0");}
				temp=allocInt1D(vals.length+1);
				temp[0]=old;
				for(int i=0; i<vals.length; i++){temp[i+1]=vals[i];}
			}else{
				if(verbose){System.err.println("E: old>0");}
				temp=vals;
			}
			values[cell]=-setList.size;
			setList.add(temp);
		}
	}

	/**
	 * Inserts a single value into a hash table cell for a k-mer.
	 * Empty cells store the positive value directly; direct cells promote to
	 * a four-slot setList array; tagged cells append to their existing list.
	 *
	 * @param kmer The k-mer key being inserted
	 * @param v The value to insert
	 * @param cell The target cell for insertion
	 */
	@Override
	protected final void insertValue(long[] kmer, int v, int cell){
		assert(matches(kmer, cell));
		assert(v>0);
		final int cc=values[cell];
		if(cc==v){
			return;
		}else if(cc<-1){
			insertIntoList(v, -cc);
		}else if(cc>0){
			values[cell]=-setList.size;
			setList.add(new int[] {cc, v, -1, -1});
		}else{
			values[cell]=v;
		}
	}

	/** Adds a value to one NOT_PRESENT-terminated setList slot, growing the local list if needed. */
	private final int insertIntoList(final int v, final int loc){

		if(loc>=setList.size){
			assert(loc==setList.size);
			setList.add(null);
		}

		int[] set=setList.get(loc);
		if(set==null){
			set=new int[] {-1, -1};
			setList.set(loc, set);
		}

		for(int i=0; i<set.length; i++){
			if(set[i]==v){return 0;}
			if(set[i]<0){set[i]=v;return 1;}
		}
		final int oldSize=set.length;
		final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
		assert(newSize>set.length) : "Overflow.";
		set=KillSwitch.copyOf(set, newSize);
		set[oldSize]=v;
		Arrays.fill(set, oldSize+1, newSize, -1);
		setList.set(loc, set);
		return 1;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns false as this implementation does not support rebalancing. */
	@Override
	public final boolean canRebalance(){return false;}

//	@Override
//	protected synchronized void resize(){
//
//		if(verbose){
//			System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
//		}
//
////		assert(TESTMODE);
////		if(TESTMODE){
////			for(int i=0; i<ll.size; i++){
////				assert(contains(ll.get(i), il.get(i)));
////				assert(!contains(ll.get(i), Integer.MAX_VALUE));
////			}
////		}
//
////		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
//		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
//		}
//
//		final long oldSize=size, oldVSize=victims.size;
//		final long totalSize=this.size+victims.size;
//
//		final long maxAllowedByLoadFactor=(long)(totalSize*minLoadMult);
//		final long minAllowedByLoadFactor=(long)(totalSize*maxLoadMult);
//
////		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);
//
//		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
//		if(maxAllowedByLoadFactor<prime){
//			sizeLimit=(long)(maxLoadFactor*prime);
//			return;
//		}
//
//		long x=10+(long)(prime*resizeMult);
//		x=Tools.max(x, minAllowedByLoadFactor);
//		x=Tools.min(x, maxAllowedByLoadFactor);
//
//		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
//
//		if(prime2<=prime){
//			sizeLimit=(long)(maxLoadFactor*prime);
//			assert(prime2==prime) : "Resizing to smaller array? "+totalSize+", "+prime+", "+x;
//			return;
//		}
//
//		prime=prime2;
////		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
//		long[][] oldk=arrays;
//		int[] oldc=values;
//		arrays=KillSwitch.allocLong2D(mult, prime2+extra);
//		for(int i=0; i<mult; i++){
//			Arrays.fill(arrays[i], NOT_PRESENT);
//		}
//		IntList2 oldList=setList;
//		setList=new IntList2();
//		setList.add(null);
//		setList.add(null);//TODO: May have to add 3 of them to avoid HASH_COLLISION at -2
//		values=allocInt1D(prime2+extra);
//		ArrayList<KmerNodeU> list=victims.toList();
//		victims.clear();
//		size=0;
//		sizeLimit=Long.MAX_VALUE;
//
//		final int[] singleton=new int[] {NOT_PRESENT};
//		final Kmer kmer=new Kmer(kbig);
//		{
//			for(int i=0; i<oldk.length; i++){
//				if(oldk[0][i]>NOT_PRESENT){
//					final int v=oldc[i];
//					fillKmer(i, kmer, oldk);
//					if(v>=0){
//						set(kmer, v);
//					}else{
//						set(kmer, oldList.get(-v));
//					}
//				}
//			}
//		}
//
//		for(KmerNodeU n : list){
//			if(n.pivot[0]>NOT_PRESENT){
//				kmer.setFrom(n.pivot());
//				if(n.numValues()>1){
//					set(kmer, n.values(singleton));
//				}else{
//					set(kmer, n.value());
//				}
//			}else{assert(false);}
//		}
//
//		assert(oldSize+oldVSize==size+victims.size) : oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size);
//
//		if(verbose){System.err.println("Resized to "+prime+". "+oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));}
//
//		sizeLimit=(long)(maxLoadFactor*prime);
//
////		assert(TESTMODE);
////		if(TESTMODE){
////			for(int i=0; i<ll.size; i++){
////				long[] kmer=ll.get(i);
////				int v=il.get(i);
////				assert(contains(kmer, v)) : i+", "+ll.size+", "+kmer+", "+v+", "+Arrays.toString(getValues(kmer, new int[1]));
////				assert(!contains(kmer, Integer.MAX_VALUE));
////			}
////		}
//	}


	/**
	 * Resizes the primary arrays after load growth exhausts the current limit.
	 * Uses the next scheduled size or retained dynamic sizing, then rehashes
	 * all primary cells and a victim-cache snapshot into fresh arrays.
	 * Existing setList arrays are reinserted by live tag. Primary owner markers
	 * are not remapped; clear ownership before growth. Callers must exclude all
	 * concurrent table access, including readers, while arrays are replaced.
	 */
	@Override
	protected synchronized void resize(){
		if(verbose){System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));}
		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
			KillSwitch.memKill(new OutOfMemoryError());
		}

		final int oldPrime=prime;
		final long oldSize=size, oldVSize=victims.size;
		final long totalSize=oldSize+oldVSize;

		if(schedule!=null){
			prime=nextScheduleSize();
			if(prime<=oldPrime){KillSwitch.memKill(new OutOfMemoryError());}
			sizeLimit=(long)((atMaxSize() ? maxLoadFactorFinal : maxLoadFactor)*prime);
		}else{//Old method
			final long maxAllowedByLoadFactor=(long)(totalSize*minLoadMult);
			final long minAllowedByLoadFactor=(long)(totalSize*maxLoadMult);

			//		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);

			assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
			if(maxAllowedByLoadFactor<prime){
				sizeLimit=(long)(maxLoadFactor*prime);
				return;
			}

			long x=10+(long)(prime*resizeMult);
			x=Tools.max(x, minAllowedByLoadFactor);
			x=Tools.min(x, maxAllowedByLoadFactor);

			int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));

			if(prime2<=prime){
				sizeLimit=(long)(maxLoadFactor*prime);
				assert(prime2==prime) : "Resizing to smaller array? "+totalSize+", "+prime+", "+x;
				return;
			}

			prime=prime2;
			sizeLimit=(long)(maxLoadFactor*prime);
		}
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		long[][] oldk=arrays;
		int[] oldc=values;
		arrays=KillSwitch.allocLong2D(mult, prime+extra);
		for(int i=0; i<mult; i++){
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		IntList2 oldList=setList;
		setList=new IntList2();
		setList.add(null);
		setList.add(null);//TODO: May have to add 3 of them to avoid HASH_COLLISION at -2
		values=allocInt1D(prime+extra);
		ArrayList<KmerNodeU> list=victims.toList();
		victims.clear();
		size=0;

		final int[] singleton=new int[] {NOT_PRESENT};
		final Kmer kmer=new Kmer(kbig);
		{
			for(int i=0; i<oldk[0].length; i++){
				//Values, not key bits, authoritatively encode occupancy; -1L is a valid packed word.
				if(oldc[i]!=0){
					final int v=oldc[i];
					fillKmer(i, kmer, oldk);
					if(v>=0){
						set(kmer, v);
					}else{
						set(kmer, oldList.get(-v));
					}
				}
			}
		}

		for(KmerNodeU n : list){
			//Victims never hold empties (see ukmer/HashArrayU#004)
			assert(n.pivot()!=null) : n;
			kmer.setFrom(n.pivot());
			if(n.numValues()>1){
				set(kmer, n.values(singleton));
			}else{
				set(kmer, n.value());
			}
		}

		assert(oldSize+oldVSize==size+victims.size) : KillSwitch.assertDie(oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));

		if(verbose){System.err.println("Resized to "+prime+". "+oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));}

		//TODO UKM107: scheduled final-size resizes should keep maxLoadFactorFinal after reinsertion.
		sizeLimit=(long)(maxLoadFactor*prime);

//		assert(TESTMODE);
//		if(TESTMODE){
//			for(int i=0; i<ll.size; i++){
//				long[] kmer=ll.get(i);
//				int v=il.get(i);
//				assert(contains(kmer, v)) : i+", "+ll.size+", "+kmer+", "+v+", "+Arrays.toString(getValues(kmer, new int[1]));
//				assert(!contains(kmer, Integer.MAX_VALUE));
//			}
//		}
	}


	/** @deprecated Rebalancing is not implemented for this hash table type. */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * @deprecated Regeneration is not supported for this hash table type.
	 * @param limit Unused threshold
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown
	 */
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Direct positive values, zero empties, and negative list tags aligned with key slots. */
	private int[] values;
	/** Live multi-value arrays addressed by negated tags; indices zero and one are reserved. */
	private IntList2 setList;

}
