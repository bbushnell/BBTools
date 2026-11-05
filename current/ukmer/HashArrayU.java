package ukmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Primes;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Open-addressed table for multiword {@link Kmer} keys with a bounded primary
 * probe window and a {@link HashForestU} victim cache for entries whose key
 * cannot be placed in the primary arrays. Subclasses provide either 1D counts
 * or 2D value-list storage; this base class owns packed-key placement, resize
 * scheduling, victim delegation, ownership cells and dump traversal. Callers
 * coordinate all shared access, finish structural changes before ownership
 * initialization and clear ownership before any later manual or automatic resize.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public abstract class HashArrayU extends AbstractKmerTableU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows a nonempty, stable prime schedule, allocates the first primary size
	 * and creates a victim cache with the same growth policy. {@code k_} is the
	 * per-word packing width; {@code kbig_} is the total decoded kmer length and
	 * determines the retained word count, including a partial tail under packed ukmers.
	 */
	HashArrayU(int[] schedule_, int k_, int kbig_, boolean twod_){
		schedule=schedule_;
		autoResize=schedule.length>1;
		prime=schedule[0];

		sizeLimit=(long)((schedule.length==1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
		k=k_;
		kbig=kbig_;
		//[ukmer/HashArrayU#003] (Amber, 2026-08-19) Was mult=kbig/k -- correct only under the
		//symmetric layout where kbig==k*mult by construction. Under Kmer.PACKED, k=FULL_WORD_K=32
		//and the last word is partial, so integer division rounds DOWN (62/32=1 vs true mult 2):
		//the table then allocated/stored/compared only word[0] of every kmer, silently truncating
		//keys. Benign-looking without coremask (unmasked hashing rarely co-locates same-word0
		//kmers), but MASK_CORE deliberately co-locates last-base extension variants -- which are
		//word0-identical -- so word0-only findKmer matching merged them (Tadpole2 packed=t
		//coremask=t 6330 vs 3421 errors). getMult is PACKED-aware; identical to kbig/k for every
		//symmetric caller. Caught by fillKey's temp.length==mult assert on the first resize.
		mult=Kmer.getMult(kbig);
		assert(k*mult>=kbig && mult==new Kmer(kbig).mult) :
			"Table geometry disagrees with Kmer: k="+k+", kbig="+kbig+", mult="+mult+", Kmer.mult="+new Kmer(kbig).mult;
		arrays=new long[mult][];
		for(int i=0; i<mult; i++){
			arrays[i]=allocLong1D(prime+extra);
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		victims=new HashForestU(Tools.max(10, prime/victimRatio), k, kbig, autoResize, twod_);
		TWOD=twod_;
	}

//	HashArrayU(int initialSize, int k_, int kbig_, boolean autoResize_, boolean twod){
//		if(initialSize>1){
//			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
//		}else{
//			initialSize=1;
//		}
//		schedule=null;
//		prime=initialSize;
//		sizeLimit=(long)(sizeLimit=(long)(maxLoadFactor*prime));
//		k=k_;
//		kbig=kbig_;
//		mult=kbig/k;
//		arrays=new long[mult][];
//		for(int i=0; i<mult; i++){
//			arrays[i]=allocLong1D(prime+extra);
//			Arrays.fill(arrays[i], NOT_PRESENT);
//		}
//		victims=new HashForestU(Tools.max(10, initialSize/victimRatio), k, autoResize_, twod);
//		autoResize=autoResize_;
//		TWOD=twod;
//	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, new int[1]);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Fetched "+Arrays.toString(old));}
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+v+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			final int[] singleton=new int[1];
//			int[] old=getValues(kmer, singleton);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Before: old="+Arrays.toString(old)+", v="+Arrays.toString(v));}
//			x=set0(kmer, v);
//			if(verbose){System.err.println("After:  old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", get()="+Arrays.toString(getValues(kmer, singleton)));}
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}

	/**
	 * Returns the first primary cell to probe for the kmer's current key/fingerprint.
	 */
	public final int kmerToCell(Kmer kmer){
		int cell=(int)(kmer.xor()%prime);
		return cell;
	}

	@Override
	public final int set(Kmer kmer, final int[] v){
		final int cell=findKmerOrEmpty(kmer);

		if(cell==HASH_COLLISION){
			if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
			final int x=victims.set(kmer, v);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
			return x;
		}
		final long[] key=kmer.key();

		assert(cell>=0);

		final boolean notpresent=cellEmpty(cell);
		if(notpresent){
			if(verbose){System.err.println("B2: Setting cell "+cell+" to kmer "+kmer);}
			setKmer(kmer.key(), cell);
		}

		if(verbose){System.err.println("A2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
		insertValue(key, v, cell);
		if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}

		if(notpresent){
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}

	/**
	 * Copies all retained packed key words into an already-selected empty primary cell.
	 * Occupancy remains subclass value-backed, so callers must install a value after
	 * copying the key words and before making the cell visible to lookup/walking.
	 */
	public final void setKmer(long[] key, int cell){
		if(verbose){System.err.println();}
		for(int i=0; i<mult; i++){
			arrays[i][cell]=key[i];
		}
	}

	@Override
	public final int set(final Kmer kmer, final int v){
		assert(kmer.mult==mult && kmer.len>=kmer.kbig);
		final int cell=findKmerOrEmpty(kmer);
//		assert(kmer.verify(false)); //123

		if(cell==HASH_COLLISION){
			if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
			final int x=victims.set(kmer, v);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
			return x;
		}
		assert(cell>=0);
		final long[] key=kmer.key();

		final boolean notpresent=cellEmpty(cell);
		if(notpresent){
			if(verbose){System.err.println("B2: Setting cell "+cell+" to kmer "+kmer);}
			setKmer(key, cell);
		}

		if(verbose){System.err.println("A2: Adding "+kmer+", "+v+", "+cell);}
		insertValue(key, v, cell);
		if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}

		if(notpresent){
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}


//	protected LongList ll=new LongList(); //123
//	protected IntList il=new IntList();

	@Override
	public final int setIfNotPresent(Kmer kmer, int value){
		final int cell=findKmerOrEmpty(kmer);

		if(cell==HASH_COLLISION){
			int x=victims.setIfNotPresent(kmer, value);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}
		assert(cell>=0);
		final long[] key=kmer.key();

		//[ukmer/HashArrayU#002] Presence must be read from the cell CONTENTS, not the index.
		//findKmerOrEmpty returns a valid index (>=0), or HASH_COLLISION (handled above) - it NEVER
		//returns the NOT_PRESENT (-1) sentinel. With the assert(cell>=0) above, 'cell==NOT_PRESENT'
		//was always false, so the insert branch was dead and this method never stored anything
		//(always returned 0). Occupancy is value-backed because every 64-bit key word is valid.
		if(cellEmpty(cell)){
			setKmer(key, cell);
			insertValue(key, value, cell);
			size++;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			return 0;
		}
	}

	@Override
	public final int getValue(Kmer kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}

	/**
	 * Looks up a kmer when the caller already knows the initial primary probe cell.
	 */
	public final int getValue(Kmer kmer, int startCell){
		int cell=findKmer(kmer, startCell);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}

	/**
	 * Alias for callers that treat the single stored value as a count.
	 */
	public final int getCount(Kmer kmer){
		return getValue(kmer);
	}

	/**
	 * Raw-key lookup is unsupported for this table family; callers use Kmer keys.
	 * @throws RuntimeException always from this base implementation
	 */
	@Override
	public int getValue(long[] key, long xor){
		throw new RuntimeException("Unimplemented");
	}

	/**
	 * Copies a live primary cell's key words into the caller-owned scratch array.
	 */
	protected final long[] fillKey(int cell, long[] temp){
		if(cellEmpty(cell)){return null;}
		return fillKey(cell, temp, arrays);
	}

	/**
	 * Copies a primary cell into the caller-owned Kmer and refreshes its alternate view.
	 */
	public final Kmer fillKmer(int cell, Kmer kmer){
		if(cellEmpty(cell)){return null;}
		return fillKmer(cell, kmer, arrays);
	}

	/**
	 * Copies a cell from an explicit matrix into a caller-owned Kmer; used while
	 * resizing reads old key arrays whose value arrays still mark occupied slots.
	 */
	public final Kmer fillKmer(int cell, Kmer kmer, long[][] matrix){
		long[] x=fillKey(cell, kmer.array1(), matrix);
//		assert(false) : x+"\ngetKmer("+cell+", kmer, matrix)"; //123
		if(x==null){return null;}
		kmer.fillArray2();
		if(verbose){System.err.println("Filled kmer "+kmer+": a1="+Arrays.toString(kmer.array1())+", a2="+Arrays.toString(kmer.array2())+", key="+Arrays.toString(kmer.key()));}
		return kmer;
	}

	/**
	 * Copies all key words from an explicit matrix; the caller is responsible for
	 * checking occupancy before reading because all long bit patterns are valid keys.
	 */
	protected final long[] fillKey(int cell, long[] temp, long[][] matrix){
		assert(temp.length==mult);
		//[ukmer/HashArrayU#004] Every 64-bit word is a valid packed key, including -1L
		//(32 T bases), so key data cannot encode occupancy.  Public/current-table callers
		//check cellEmpty(); resize callers use their saved value arrays before calling here.
		for(int i=0; i<temp.length; i++){
			temp[i]=matrix[i][cell];
		}
		if(verbose){System.err.println("cell="+cell+", matrix[0][cell]="+matrix[0][cell]+", temp="+Arrays.toString(temp)+"\nmatrix[0]="+Arrays.toString(matrix[0]));}
		return temp;
	}

	@Override
	public final int[] getValues(Kmer kmer, int[] singleton){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){
			singleton[0]=NOT_PRESENT;
			return singleton;
		}
		if(cell==HASH_COLLISION){return victims.getValues(kmer, singleton);}
		return readCellValues(cell, singleton);
	}

	@Override
	public final boolean contains(Kmer kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return false;}
		if(cell==HASH_COLLISION){return victims.contains(kmer);}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void initializeOwnership(){
		assert(owners==null);
		owners=allocAtomicInt(arrays[0].length);
		for(int i=0; i<arrays[0].length; i++){
			owners.set(i, NO_OWNER);
		}
		victims.initializeOwnership();
	}

	@Override
	public final void clearOwnership(){
		owners=null;
		victims.clearOwnership();
	}

	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.setOwner(kmer, newOwner);}
		return setOwner(kmer, newOwner, cell);
	}

	/**
	 * Raises the owner id for a known primary cell and leaves larger existing owners in
	 * place; callers must have called {@link #initializeOwnership()} first.
	 */
	public final int setOwner(final Kmer kmer, final int newOwner, final int cell){
//		kmer.verify(true);
		assert(matches(kmer.key(), cell)) : "cell="+cell+", key="+Arrays.toString(kmer.key())+", row="+Arrays.toString(cellToArray(cell))+"\n" +
				"kmer="+kmer+", array1="+Arrays.toString(kmer.array1())+", array2="+Arrays.toString(kmer.array2())+", row="+AbstractKmerTableU.toText(cellToArray(cell), kmer.kbig);
		final int original=owners.get(cell);
		int current=original;
		while(current<newOwner){
			boolean success=owners.compareAndSet(cell, current, newOwner);
			if(!success){current=owners.get(cell);}
			else{current=newOwner;}
		}
		assert(current>=original) : "original="+original+", current="+current+", newOwner="+newOwner+", re-read="+owners.get(cell);
		return current;
	}

	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.clearOwner(kmer, owner);}
		return clearOwner(kmer, owner, cell);
	}

	/**
	 * Releases a known primary cell only if its current owner still matches.
	 */
	public final boolean clearOwner(final Kmer kmer, final int owner, final int cell){
		assert(matches(kmer.key(), cell));
		boolean success=owners.compareAndSet(cell, owner, NO_OWNER);
		return success;
	}

	@Override
	public final int getOwner(final Kmer kmer){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.getOwner(kmer);}
		return getCellOwner(cell);
	}

	/**
	 * Returns the raw owner id for a primary cell.
	 */
	public final int getCellOwner(final int cell){
		return owners.get(cell);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores or merges a scalar value for a key already present at {@code cell}.
	 */
	protected abstract void insertValue(final long[] kmer, final int v, final int cell);

	/**
	 * Stores or merges a value list for a key already present at {@code cell}.
	 */
	protected abstract void insertValue(final long[] kmer, final int[] vals, final int cell);

	/** Reads the subclass's scalar value from an occupied primary cell. */
	protected abstract int readCellValue(int cell);

	/** Returns an occupied primary cell's values, using singleton for 1D storage. */
	protected abstract int[] readCellValues(int cell, int[] singleton);

	/** Returns true when the subclass's value storage marks this cell unused. */
	protected abstract boolean cellEmpty(int cell);

	@Override
	final Object get(long[] kmer){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Finds a key from its natural probe origin in the primary arrays.
	 * @return the matching primary cell, {@link #NOT_PRESENT}, or {@link #HASH_COLLISION}
	 */
	final int findKmer(Kmer kmer){
		return findKmer(kmer.key(), (int)(kmer.xor()%prime));
	}

	/**
	 * Finds a key from a caller-supplied primary probe origin.
	 */
	final int findKmer(Kmer kmer, int startCell){
		return findKmer(kmer.key(), startCell);
	}

//	final int findKmer(long[] key, long xor){
//		return findKmer(key, (int)(xor%prime));
//	}

//	final int findKmer(final long[] key, final int startCell){
//		int cell=startCell;
//		for(final int max=cell+extra; cell<max; cell++){
//			final long n=arrays[0][cell];
//			if(n==key[0]){
//				boolean success=true;
//				for(int i=1; i<mult && success; i++){
//					if(key[i]!=arrays[i][cell]){success=false;}
//				}
//				if(success){return cell;}
//			}else if(n==NOT_PRESENT){return NOT_PRESENT;}
//		}
//		return HASH_COLLISION;
//	}

	/**
	 * Probes up to {@link #extra} primary cells for a raw retained key.
	 * An empty value cell terminates the chain; an exhausted chain means the key, if
	 * present, must live in the victim cache.
	 */
	final int findKmer(final long[] key, final int startCell){
		int cell=startCell;

		final long[] array0=arrays[0];
		final long key0=key[0];
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array0[cell];
			if(cellEmpty(cell)){return NOT_PRESENT;}
			else if(n==key0){
				boolean success=true;
				for(int i=1; i<mult && success; i++){
					if(key[i]!=arrays[i][cell]){success=false;}
				}
				if(success){return cell;}
			}
		}
		return HASH_COLLISION;
	}

	/**
	 * Returns a matching primary cell or the first empty cell in the bounded probe
	 * window; {@link #HASH_COLLISION} sends the operation to the victim cache.
	 */
	final int findKmerOrEmpty(Kmer kmer){
		int cell=kmerToCell(kmer);
		if(verbose){System.err.println("Started at cell "+cell+" for "+kmer);}

		final long[] key=kmer.key();
		final long[] array0=arrays[0];
		final long key0=key[0];
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array0[cell];
			if(cellEmpty(cell)){
				if(verbose){System.err.println("Chose empty cell "+cell+" for "+kmer);}
				return cell;
			}else if(n==key0){
				boolean success=true;
				for(int i=1; i<mult && success; i++){
					if(key[i]!=arrays[i][cell]){success=false;}
				}
				if(success){
					if(verbose){System.err.println("Found cell "+cell+" containing "+kmer);}
					return cell;
				}
			}
		}
		return HASH_COLLISION;
	}

	/**
	 * Returns whether all retained key words match the supplied primary cell.
	 */
	final boolean matches(long[] key, int cell){
		assert(cell>=0);
		for(int i=0; i<mult; i++){
			if(key[i]!=arrays[i][cell]){return false;}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns true; concrete subclasses provide the actual synchronized resize. */
	@Override
	final boolean canResize(){return true;}

	/** Returns the primary-entry count; caller must add victims for a total size. */
	@Override
	final public long size(){return size;}

	/** Returns the allocated primary-array length, including the linear-probe tail. */
	@Override
	final public int arrayLength(){return arrays[0].length;}

	/**
	 * Resizes primary storage and reinserts existing primary and victim entries.
	 * Concrete subclasses replace the arrays instead of remapping ownership; callers
	 * must exclude all concurrent access and clear ownership before manual or
	 * automatic resizing.
	 */
	@Override
	protected abstract void resize();

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/** Unsupported base hook for converting a primary cell's key into a new array. */
	protected long[] cellToArray(int cell){throw new RuntimeException("Unimplemented");}

	/**
	 * Dumps occupied primary cells then victim entries as text.
	 * Scalar tables filter by mincount/maxcount; list-valued tables dump each
	 * retained value list without threshold filtering. {@code k} is the total kbig
	 * decoding length. A {@code true} return does not certify writer completion.
	 */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					tsw.print(toText(temp, readCellValues(i, singleton), k).append('\n'));
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount && readCellValue(i)<=maxcount){
					tsw.print(toText(temp, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}

	/**
	 * Dumps occupied primary cells then victim entries as text records encoded as bytes.
	 * Scalar tables filter by mincount/maxcount; list-valued tables dump each
	 * retained value list without threshold filtering. {@code k} is the total kbig
	 * decoding length. A {@code true} return may include a quota early-stop and
	 * does not certify writer completion.
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(temp, readCellValues(i, singleton), k);
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount && readCellValue(i)<=maxcount){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(temp, readCellValue(i), k);
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/**
	 * Dumps text-byte records while appending primary cells to the caller-thread's
	 * exclusive buffer. Scalar tables filter by mincount/maxcount; list-valued
	 * tables dump each retained value list without threshold filtering. {@code k}
	 * is the total kbig decoding length. A {@code true} return may include a quota
	 * early-stop; the caller submits the final tail, finishes the writer and checks
	 * writer state.
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, final int maxcount, AtomicLong remaining){
		final long[] key=new long[mult];
		final int alen=arrays[0].length;
		if(TWOD){
			final int[] singleton=new int[1];
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(temp, readCellValues(i, singleton), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}else{
			for(int i=0; i<alen; i++){
				long[] temp=fillKey(i, key);
				if(temp!=null && readCellValue(i)>=mincount && readCellValue(i)<=maxcount){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(temp, readCellValue(i), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/** Adds scalar primary counts and victims to a dense histogram capped at {@code max}. */
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			if(!cellEmpty(i)){
				int count=Tools.min(readCellValue(i), max);
				ca[count]++;
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}

	/** Adds scalar primary counts and victim counts to a sparse overflow-capable histogram. */
	@Override
	public void fillHistogram(SuperLongList sll){
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			if(!cellEmpty(i)){
				int count=readCellValue(i);
				sll.add(count);
			}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}

	/** Adds GC tallies for each retained primary word and delegates victim keys. */
	@Override
	public final void countGC(long[] gcCounts, int max){
		final int alen0=arrays.length;
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			if(!cellEmpty(i)){
				int count=Tools.min(readCellValue(i), max);
				for(int j=0; j<alen0; j++){
					gcCounts[count]+=gc(arrays[j][i]);
				}
			}
		}
		if(victims!=null){
			victims.countGC(gcCounts, max);
		}
	}

	/** Returns the overflow tree that owns keys outside the bounded primary probe window. */
	public HashForestU victims(){
		return victims;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional owner id per primary cell; null until initializeOwnership(), not remapped by resize. */
	AtomicIntegerArray owners;

	/** Retained packed Kmer words by word index, then primary cell index. */
	long[][] arrays;

	/** Number of hashable primary cells, excluding the extra linear-probe tail. */
	int prime;

	/** Number of entries stored directly in primary arrays; excludes victims. */
	long size=0;

	/** Combined primary-plus-victim size threshold for the next resize. */
	long sizeLimit;

	/** Overflow tree for keys that exhaust their bounded primary probe windows. */
	final HashForestU victims;

	/** True when constructed with a multi-prime schedule; remains true at the final prime. */
	final boolean autoResize;

	/** Number of bases held by a full packed word in this table layout. */
	final int k;

	/** Full kmer length in bases, possibly ending with a partial word. */
	final int kbig;

	/** Count of retained packed words in each key. */
	final int mult;//Length of Kmer arrays.

	/** True for list-valued tables; false for scalar-count tables. */
	public final boolean TWOD;

	/** Caller/shard lock exposed through getLock(); resize itself is synchronized. */
	private final Lock lock=new ReentrantLock();

	/** Returns the primary ownership cells; caller must handle a null uninitialized array. */
	public AtomicIntegerArray owners(){return owners;}

	/** Advances one schedule step when possible and otherwise keeps returning the final prime. */
	protected int nextScheduleSize(){
		if(schedulePos<schedule.length-1){schedulePos++;}
		return schedule[schedulePos];
	}

	/** Returns true after the resize schedule has reached its final prime. */
	protected boolean atMaxSize(){
		return schedulePos>=schedule.length-1;
	}

	/** Prime primary-table sizes, in resize order. */
	protected final int[] schedule;

	/** Index of the currently allocated prime in schedule. */
	private int schedulePos=0;

	/** Returns the caller/shard lock paired with this table. */
	@Override
	final Lock getLock(){return lock;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Initial divisor for victim-cache capacity; HashForestU can resize itself. */
	final static int victimRatio=16; //Initial divisor for victim cache size; it self-resizes.

	/**
	 * Primary overflow window after the prime hash range; probes past this point
	 * are delegated to victims.
	 */
	final static int extra=60; //Amazingly, increasing this gave increasing returns past 300.  Old default was 21.  Could allow higher maxLoadFactorFinal and smaller victim cache.

	/** Largest prime that still leaves room for the extra probe tail in an int index. */
	final static int maxPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra-20);

	/** Legacy geometric resize constant; scheduled resizes choose explicit primes. */
	final static float resizeMult=2f; //Resize by a minimum of this much; not needed for schedule

	/** Legacy minimum post-resize load; scheduled resizes choose explicit primes. */
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor; not needed for schedule

	/** Normal resize threshold for non-final scheduled tables. */
	final static float maxLoadFactor=0.88f; //Reaching this load triggers resizing

	/**
	 * Higher threshold used when a fixed or final scheduled table should stop growing.
	 */
	final static float maxLoadFactorFinal=0.95f; //Reaching this load triggers killing

	/** Reciprocal legacy minimum load factor. */
	final static float minLoadMult=1/minLoadFactor;

	/** Reciprocal normal resize load factor. */
	final static float maxLoadMult=1/maxLoadFactor;

}
