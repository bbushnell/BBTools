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
 * Stores kmers in a long[] and values in an int[][], with a victim cache.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public abstract class HashArrayU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs HashArrayU with configurable resizing schedule and k-mer parameters.
	 * Initializes arrays, victim cache, and size limits based on load factor settings.
	 *
	 * @param schedule_ Array of prime sizes for progressive resizing
	 * @param k_ Basic k-mer length
	 * @param kbig_ Extended k-mer length for multi-dimensional storage
	 * @param twod_ Whether to use two-dimensional value storage
	 */
	HashArrayU(int[] schedule_, int k_, int kbig_, boolean twod_){
		schedule=schedule_;
		autoResize=schedule.length>1;
		prime=schedule[0];
		
		sizeLimit=(long)((schedule.length==1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
		k=k_;
		kbig=kbig_;
		mult=kbig/k;
		arrays=new long[mult][];
		for(int i=0; i<mult; i++){
			arrays[i]=allocLong1D(prime+extra);
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		victims=new HashForestU(Tools.max(10, prime/victimRatio), k, autoResize, twod_);
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
	 * Computes hash table cell position for a k-mer.
	 * Uses XOR-based hash function modulo prime table size.
	 * @param kmer K-mer to hash
	 * @return Cell index in hash table
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
		
		final boolean notpresent=(arrays[0][cell]==NOT_PRESENT);
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
	 * Stores k-mer key in specified hash table cell.
	 * Copies all components of multi-dimensional k-mer key to arrays.
	 * @param key Multi-dimensional k-mer key array
	 * @param cell Target cell index in hash table
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
		
		final boolean notpresent=(arrays[0][cell]==NOT_PRESENT);
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
		
		if(cell==NOT_PRESENT){
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
	 * Retrieves value for k-mer starting search from specified cell.
	 * Allows optimized lookup when approximate position is known.
	 *
	 * @param kmer K-mer to look up
	 * @param startCell Starting cell for linear probe search
	 * @return Associated value, or NOT_PRESENT if k-mer not found
	 */
	public final int getValue(Kmer kmer, int startCell){
		int cell=findKmer(kmer, startCell);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}
	
	/**
	 * Gets count value for a k-mer (alias for getValue).
	 * @param kmer K-mer to count
	 * @return Count value associated with k-mer
	 */
	public final int getCount(Kmer kmer){
		return getValue(kmer);
	}

	/* (non-Javadoc)
	 * @see ukmer.AbstractKmerTableU#getValue(long[], long)
	 */
	@Override
	public int getValue(long[] key, long xor) {
		throw new RuntimeException("Unimplemented");
	}
	
	/**
	 * Fills temporary array with k-mer key from specified cell.
	 * Delegates to multi-array version using instance arrays.
	 *
	 * @param cell Cell index to read from
	 * @param temp Temporary array to fill with k-mer data
	 * @return Filled key array, or null if cell is empty
	 */
	protected final long[] fillKey(int cell, long[] temp) {
		return fillKey(cell, temp, arrays);
	}
	
	/**
	 * Fills Kmer object with data from specified hash table cell.
	 * @param cell Cell index to read from
	 * @param kmer Kmer object to fill
	 * @return Filled Kmer object, or null if cell is empty
	 */
	public final Kmer fillKmer(int cell, Kmer kmer) {
		return fillKmer(cell, kmer, arrays);
	}
	
	/**
	 * Fills Kmer object with data from cell in specified matrix.
	 * Allows reading from different array matrices during resizing operations.
	 *
	 * @param cell Cell index to read from
	 * @param kmer Kmer object to fill
	 * @param matrix Matrix to read k-mer data from
	 * @return Filled Kmer object, or null if cell is empty
	 */
	public final Kmer fillKmer(int cell, Kmer kmer, long[][] matrix) {
		long[] x=fillKey(cell, kmer.array1(), matrix);
//		assert(false) : x+"\ngetKmer("+cell+", kmer, matrix)"; //123
		if(x==null){return null;}
		kmer.fillArray2();
		if(verbose){System.err.println("Filled kmer "+kmer+": a1="+Arrays.toString(kmer.array1())+", a2="+Arrays.toString(kmer.array2())+", key="+Arrays.toString(kmer.key()));}
		return kmer;
	}
	
	/**
	 * Fills temporary array with k-mer key from cell in specified matrix.
	 * Used during resizing when reading from different array configurations.
	 *
	 * @param cell Cell index to read from
	 * @param temp Temporary array to fill
	 * @param matrix Matrix to read from
	 * @return Filled key array, or null if cell is empty
	 */
	protected final long[] fillKey(int cell, long[] temp, long[][] matrix) {
		assert(temp.length==mult);
		if(matrix[0][cell]<0){
//			assert(false) : matrix[0][cell]+"\ngetKmer("+cell+", kmer, matrix)\n"+Arrays.toString(matrix[0]); //123
			return null;
		}
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
	 * Sets thread ownership for k-mer at known cell position.
	 * Uses compare-and-swap to atomically update ownership when new owner has higher ID.
	 *
	 * @param kmer K-mer to claim ownership of
	 * @param newOwner Thread ID claiming ownership
	 * @param cell Known cell position containing the k-mer
	 * @return Current owner ID after operation
	 */
	public final int setOwner(final Kmer kmer, final int newOwner, final int cell){
//		kmer.verify(true);
		assert(matches(kmer.key(), cell)) : "cell="+cell+", key="+Arrays.toString(kmer.key())+", row="+Arrays.toString(cellToArray(cell))+"\n" +
				"kmer="+kmer+", array1="+Arrays.toString(kmer.array1())+", array2="+Arrays.toString(kmer.array2())+", row="+AbstractKmerTableU.toText(cellToArray(cell), kmer.k);
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
	 * Clears thread ownership for k-mer at known cell position.
	 *
	 * @param kmer K-mer to release ownership of
	 * @param owner Thread ID that should own the k-mer
	 * @param cell Known cell position containing the k-mer
	 * @return true if ownership was successfully cleared, false otherwise
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
	 * Gets thread owner of specified hash table cell.
	 * @param cell Cell index to check ownership of
	 * @return Thread ID of current owner
	 */
	public final int getCellOwner(final int cell){
		return owners.get(cell);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Inserts single value for k-mer at specified cell.
	 * Abstract method implemented by concrete subclasses for value storage.
	 *
	 * @param kmer K-mer key array
	 * @param v Value to insert
	 * @param cell Cell position for insertion
	 */
	protected abstract void insertValue(final long[] kmer, final int v, final int cell);
	
	/**
	 * Inserts value array for k-mer at specified cell.
	 * Abstract method implemented by concrete subclasses for multi-value storage.
	 *
	 * @param kmer K-mer key array
	 * @param vals Values to insert
	 * @param cell Cell position for insertion
	 */
	protected abstract void insertValue(final long[] kmer, final int[] vals, final int cell);

	/**
	 * Reads single value from specified cell.
	 * Abstract method implemented by subclasses based on storage strategy.
	 * @param cell Cell to read from
	 * @return Value stored in cell
	 */
	protected abstract int readCellValue(int cell);
	/**
	 * Reads all values from specified cell into array.
	 * Uses singleton for single values to minimize allocations.
	 *
	 * @param cell Cell to read from
	 * @param singleton Single-element array for return values
	 * @return Array of values from cell
	 */
	protected abstract int[] readCellValues(int cell, int[] singleton);

	@Override
	final Object get(long[] kmer){
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Finds cell containing specified k-mer.
	 * Uses k-mer's XOR hash to determine starting search position.
	 * @param kmer K-mer to find
	 * @return Cell index if found, NOT_PRESENT if not found, HASH_COLLISION if overflowed
	 */
	final int findKmer(Kmer kmer){
		return findKmer(kmer.key(), (int)(kmer.xor()%prime));
	}
		
	/**
	 * Finds k-mer starting search from specified cell.
	 * @param kmer K-mer to find
	 * @param startCell Starting position for linear probe search
	 * @return Cell index if found, NOT_PRESENT if not found, HASH_COLLISION if overflowed
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
	 * Finds k-mer using raw key starting from specified cell.
	 * Performs linear probing with early termination on empty cells.
	 * Optimized with local reference to first array for performance.
	 *
	 * @param key Multi-dimensional k-mer key to search for
	 * @param startCell Starting cell for linear probe search
	 * @return Cell index if found, NOT_PRESENT if not found, HASH_COLLISION if overflowed
	 */
	final int findKmer(final long[] key, final int startCell){
		int cell=startCell;
		
		final long[] array0=arrays[0];
		final long key0=key[0];
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array0[cell];
			if(n==key0){
				boolean success=true;
				for(int i=1; i<mult && success; i++){
					if(key[i]!=arrays[i][cell]){success=false;}
				}
				if(success){return cell;}
			}else if(n==NOT_PRESENT){return NOT_PRESENT;}
		}
		return HASH_COLLISION;
	}
	
	/**
	 * Finds k-mer or returns first empty cell for insertion.
	 * Used during insertion to locate existing k-mer or find insertion point.
	 * @param kmer K-mer to find or place
	 * @return Cell index containing k-mer or empty cell, HASH_COLLISION if table full
	 */
	final int findKmerOrEmpty(Kmer kmer){
		int cell=kmerToCell(kmer);
		if(verbose){System.err.println("Started at cell "+cell+" for "+kmer);}
		
		final long[] key=kmer.key();
		final long[] array0=arrays[0];
		final long key0=key[0];
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array0[cell];
			if(n==NOT_PRESENT){
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
	 * Checks if k-mer key matches content of specified cell.
	 * Compares all components of multi-dimensional key.
	 *
	 * @param key K-mer key to compare
	 * @param cell Cell to compare against
	 * @return true if key matches cell contents, false otherwise
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
	
	@Override
	final boolean canResize() {return true;}
	
	@Override
	final public long size() {return size;}
	
	@Override
	final public int arrayLength() {return arrays[0].length;}
	
	@Override
	protected abstract void resize();
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts cell contents to array representation.
	 * Currently unimplemented - throws RuntimeException.
	 *
	 * @param cell Cell to convert
	 * @return Array representation of cell contents
	 * @throws RuntimeException Always thrown as method is unimplemented
	 */
	protected long[] cellToArray(int cell){throw new RuntimeException("Unimplemented");}
	
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
				if(temp!=null && readCellValue(i)>=mincount){
					tsw.print(toText(temp, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
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
				if(temp!=null && readCellValue(i)>=mincount){
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
				if(temp!=null && readCellValue(i)>=mincount){
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
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int alen=arrays[0].length;
		final long[] array0=arrays[0];
		for(int i=0; i<alen; i++){
			long kmer=array0[i];
			if(kmer!=NOT_PRESENT){
				int count=Tools.min(readCellValue(i), max);
				ca[count]++;
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}
	
	@Override
	public void fillHistogram(SuperLongList sll){
		final int alen=arrays[0].length;
		final long[] array0=arrays[0];
		for(int i=0; i<alen; i++){
			long kmer=array0[i];
			if(kmer!=NOT_PRESENT){
				int count=readCellValue(i);
				sll.add(count);
			}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}
	
	@Override
	public final void countGC(long[] gcCounts, int max){
		final int alen0=arrays.length;
		final int alen=arrays[0].length;
		for(int i=0; i<alen; i++){
			long kmer0=arrays[0][i];
			if(kmer0!=NOT_PRESENT){
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
	
	/** Gets the victim cache used for collision handling.
	 * @return HashForestU instance serving as victim cache */
	public HashForestU victims(){
		return victims;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Thread ownership tracking array for concurrent access control */
	AtomicIntegerArray owners;
	/** Multi-dimensional storage arrays for k-mer keys */
	long[][] arrays;
	/** Current prime number used as hash table size */
	int prime;
	/** Current number of k-mers stored in the table */
	long size=0;
	/** Maximum k-mers allowed before resizing is triggered */
	long sizeLimit;
	/** Victim cache for handling hash table collisions */
	final HashForestU victims;
	/** Whether automatic resizing is enabled when size limit is reached */
	final boolean autoResize;
	/** Basic k-mer length */
	final int k;
	/** Extended k-mer length for multi-dimensional storage */
	final int kbig;
	/** Length of Kmer arrays (kbig/k multiplier) */
	final int mult;//Length of Kmer arrays.
	/** Whether two-dimensional value storage is enabled */
	public final boolean TWOD;
	/** Reentrant lock for synchronized access during resizing operations */
	private final Lock lock=new ReentrantLock();
	
	/** Gets the thread ownership tracking array.
	 * @return AtomicIntegerArray for ownership management */
	public AtomicIntegerArray owners() {return owners;}
	
	/** Advances to next size in resizing schedule.
	 * @return Next prime size for hash table resize */
	protected int nextScheduleSize(){
		if(schedulePos<schedule.length-1){schedulePos++;}
		return schedule[schedulePos];
	}
	
	/** Checks if hash table has reached maximum scheduled size.
	 * @return true if at final size in schedule, false if more sizes available */
	protected boolean atMaxSize(){
		return schedulePos>=schedule.length-1;
	}
	
	/** Array of prime sizes for progressive hash table resizing */
	protected final int[] schedule;
	/** Current position in the resizing schedule array */
	private int schedulePos=0;
	
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Initial divisor for victim cache size (16) - victim cache self-resizes */
	final static int victimRatio=16; //Initial divisor for victim cache size; it self-resizes.
	/**
	 * Extra slots for linear probing (60) - increasing this improved performance past 300
	 */
	final static int extra=60; //Amazingly, increasing this gave increasing returns past 300.  Old default was 21.  Could allow higher maxLoadFactorFinal and smaller victim cache.
	/** Maximum prime number that can be used as hash table size */
	final static int maxPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra-20);
	/** Minimum resize multiplier (2.0f) - not needed when using schedule */
	final static float resizeMult=2f; //Resize by a minimum of this much; not needed for schedule
	/** Minimum load factor after resize (0.58f) - not needed when using schedule */
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor; not needed for schedule
	/** Load factor that triggers resizing (0.88f) */
	final static float maxLoadFactor=0.88f; //Reaching this load triggers resizing
	/**
	 * Load factor that triggers termination when no more resizing possible (0.95f)
	 */
	final static float maxLoadFactorFinal=0.95f; //Reaching this load triggers killing
	/** Multiplier derived from minimum load factor (1/minLoadFactor) */
	final static float minLoadMult=1/minLoadFactor;
	/** Multiplier derived from maximum load factor (1/maxLoadFactor) */
	final static float maxLoadMult=1/maxLoadFactor;
	
}
