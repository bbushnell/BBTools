package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Tools;
import structures.SuperLongList;

/**
 * Single-value hash table for storing k-mer counts using parallel arrays.
 * Uses a long[] for k-mer storage and int[] for count values, with victim cache for collisions.
 * Provides memory-efficient k-mer counting with automatic resizing capabilities.
 * @author Brian Bushnell
 * @date Oct 25, 2013
 */
public final class HashArrayU1D extends HashArrayU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArrayU1D(int[] schedule_, int k_, int kbig_){
		super(schedule_, k_, kbig_, false);
		values=allocInt1D(prime+extra);
	}
	
//	public HashArrayU1D(int initialSize, int k_, int kbig_, boolean autoResize_){
//		super(initialSize, k_, kbig_, autoResize_, false);
//		values=allocInt1D(prime+extra);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);
		
		if(cell==HASH_COLLISION){
			int x=victims.increment(kmer);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}else if(arrays[0][cell]==NOT_PRESENT){
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
	 * Increments k-mer count and returns number of new entries created.
	 * Returns 1 for new k-mers, 0 for existing k-mers that were incremented.
	 * Triggers automatic resizing if size limits are exceeded.
	 * @param kmer The k-mer to increment
	 * @return 1 if a new k-mer was created, 0 if an existing k-mer was incremented
	 */
	@Override
	public final int incrementAndReturnNumCreated(final Kmer kmer){
//		assert(kmer.verify(false));
////		System.err.println("***");
////		System.err.println("Incrementing kmer "+kmer+"\n"+kmer.arraysToString());
////		System.err.println("Initial state:"+Arrays.toString(arrays[0])+"\n"+Arrays.toString(values)+"\nVictims.size: "+victims.size);
//		final int a=getValue(kmer);
//		final int x=incrementAndReturnNumCreated0(kmer);
//		final int b=getValue(kmer);
////		System.err.println("Kmer is now       "+kmer+"\n"+kmer.arraysToString());
//		assert(kmer.verify(false));
//		assert((a==-1 && b==1) || (a+1==b)) : a+", "+b+", "+kmer+"\n"+kmer.arraysToString()+"\n"+Arrays.toString(arrays[0])+"\n"+Arrays.toString(values);
//		return x;
//	}
//
//	public final int incrementAndReturnNumCreated0(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);
//		assert(victims.size<size+100);
//		System.err.println("size="+size+", victims="+victims.size+", sizeLimit="+sizeLimit+", autoResize="+autoResize);//123
		if(cell==HASH_COLLISION){
//			if(verbose || true){System.err.println("HASH_COLLISION - sending to victims.");}
			final int x=victims.incrementAndReturnNumCreated(kmer);
			if(autoResize && size+victims.size>sizeLimit){
				if(verbose){System.err.println("Exceeded size limit - resizing.");}
				resize();
			}
//			else{
				assert(!autoResize || size+victims.size<=sizeLimit+1) : sizeLimit+"<"+(size+victims.size)+", size="+size+", victims="+victims.size+", prime="+prime;
//			}
			return x;
		}else if(arrays[0][cell]==NOT_PRESENT){
			setKmer(kmer.key(), cell);
			size++;
			values[cell]=1;
			if(verbose){System.err.println("Added kmer "+kmer+", key "+Arrays.toString(kmer.key())+
					", a1 "+Arrays.toString(kmer.array1())+", a2 "+Arrays.toString(kmer.array2())+", xor "+kmer.xor()+", to cell "+cell+"\n" +
					"   array:"/*+Arrays.toString(arrays[0])*/);}
			if(autoResize && size+victims.size>sizeLimit){
				if(verbose){System.err.println("Exceeded size limit - resizing.");}
				resize();
			}
//			else{
				assert(!autoResize || size+victims.size<=sizeLimit+1) : sizeLimit+"<"+(size+victims.size)+", size="+size+", victims="+victims.size+", prime="+prime;
//			}
			return 1;
		}else{
			if(verbose){System.err.println("Already present - incrementing.");}
			assert(!autoResize || size+victims.size<=sizeLimit+1) : sizeLimit+"<"+(size+victims.size)+", size="+size+", victims="+victims.size+", prime="+prime;
			values[cell]++;
			if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
			return 0;
		}
	}
	
	/**
	 * Populates a histogram with all count values from this hash array.
	 * Includes counts from both main array and victim cache.
	 * @param sll The SuperLongList to fill with count values
	 */
	@Override
	public final void fillHistogram(SuperLongList sll){
		for(int i=0; i<values.length; i++){
			int count=values[i];
			if(count>0){sll.add(count);}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads the count value stored at a specific cell.
	 * @param cell The cell index to read
	 * @return The count value at the specified cell
	 */
	@Override
	public final int readCellValue(int cell) {
		return values[cell];
	}
	
	/**
	 * Reads cell value into a provided array for consistent interface.
	 * Used by superclass methods that expect array format.
	 * @param cell The cell index to read
	 * @param singleton Single-element array to store the value
	 * @return The singleton array with value at index 0
	 */
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		singleton[0]=values[cell];
		return singleton;
	}
	
	/**
	 * Inserts a single count value at the specified cell.
	 * Used during resizing and rehashing operations.
	 * @param kmer The k-mer key (for verification)
	 * @param v The count value to insert
	 * @param cell The target cell index
	 */
	@Override
	protected final void insertValue(long[] kmer, int v, int cell) {
		assert(matches(kmer, cell));
		values[cell]=v;
	}
	
	/**
	 * Inserts count value from array format at the specified cell.
	 * Array must contain exactly one element.
	 * @param kmer The k-mer key (for verification)
	 * @param vals Array containing the count value at index 0
	 * @param cell The target cell index
	 */
	@Override
	protected final void insertValue(long[] kmer, int[] vals, int cell) {
		assert(matches(kmer, cell));
		assert(vals.length==1);
		values[cell]=vals[0];
	}
	
	/**
	 * Converts a cell's k-mer data to array format.
	 * Creates array by copying from all k-mer storage arrays.
	 * @param cell The cell index to convert
	 * @return Array representation of the k-mer at the cell
	 */
	@Override
	protected long[] cellToArray(int cell){
		long[] r=new long[mult];
		for(int i=0; i<mult; i++){r[i]=arrays[i][cell];}
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Indicates whether this hash array supports rebalancing operations.
	 * @return false as this implementation does not support rebalancing */
	@Override
	public final boolean canRebalance() {return false;}
	
//	@Override
//	protected synchronized void resize(){
//		if(verbose){System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));}
//		final int oldPrime=prime;
//		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
//		}
//		
//		final long oldSize=size, oldVSize=victims.size;
//		final long totalSize=oldSize+oldVSize;
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
//		
//	}
	

	
	/**
	 * Resizes the hash table when load factor limits are exceeded.
	 * Allocates larger arrays, rehashes existing entries, and migrates victim cache.
	 * Uses scheduled sizing if available, otherwise calculates optimal size.
	 * Terminates program if maximum table size is exceeded.
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
//		KmerNodeU[] oldv=victims.array;
		arrays=allocLong2D(mult, prime+extra);
		for(int i=0; i<mult; i++){
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		values=allocInt1D(prime+extra);
		ArrayList<KmerNodeU> list=victims.toList();
		victims.clear();
		size=0;
//		assert(false);
		final Kmer kmer=new Kmer(kbig);
//		long kmersProcessed=0; //123
		{
			for(int i=0; i<oldk[0].length; i++){
//				assert(false) : oldk[0][i];
				if(oldk[0][i]>NOT_PRESENT){
//					kmersProcessed++;
//					assert(false) : oldk[0][i];
					Kmer temp=fillKmer(i, kmer, oldk);
					assert(temp==kmer);
					if(verbose){
						System.err.println("In cell "+i+", found kmer "+kmer+"; key="+Arrays.toString(kmer.key())+"; " +
								"a1="+Arrays.toString(kmer.array1())+"; a2="+Arrays.toString(kmer.array2()));
						System.err.println(Arrays.toString(oldk[0]));
						System.err.println(Arrays.toString(arrays[0]));
					}
					assert(temp!=null) : i+", "+kmer+", "+oldk[0][i];
					set(temp, oldc[i]);
					
//					assert(getValue(temp)==oldc[i]); //123
					
					if(verbose){
						System.err.println("prime="+prime+", xor="+kmer.xor()+", mod="+(kmer.xor()%prime));
						System.err.println("After set: kmer "+kmer+"; key="+Arrays.toString(kmer.key())+"; " +
								"a1="+Arrays.toString(kmer.array1())+"; a2="+Arrays.toString(kmer.array2()));
						System.err.println(Arrays.toString(arrays[0]));
					}
//					assert(kmer.verify(false)); //123
				}
			}
		}

		for(KmerNodeU n : list){
			if(n.pivot[0]>NOT_PRESENT){
				kmer.setFrom(n.pivot());
				set(kmer, n.value());
//				assert(getValue(kmer)==n.value()); //123 slow
			}
			else{assert(false) : "pivot="+n.pivot()+", n="+n;}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size+"; totalSize="+totalSize+", new total="+(size+victims.size)+
			"\noldPrime="+oldPrime+", prime="+prime+(prime<1000 ? (
			"\noldArray:"+Arrays.toString(oldk[0])+
			"\nnewArray:"+Arrays.toString(arrays[0])
			) : "");
	}
	
	/**
	 * Rebalancing operation is not supported in this implementation.
	 * @throws RuntimeException Always thrown as rebalancing is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Removes k-mers with counts at or below the specified limit.
	 * Processes both main arrays and victim cache, keeping only high-count k-mers.
	 * @param limit Count threshold; k-mers with counts <= limit are removed
	 * @return Number of k-mers that were removed
	 */
	@Override
	public long regenerate(final int limit){
		long sum=0;
		assert(owners==null) : "Clear ownership before regeneration.";
		final Kmer kmer=new Kmer(kbig);
		for(int pos=0; pos<values.length; pos++){
			Kmer key=fillKmer(pos, kmer);
			if(key!=null){
				final int value=values[pos];
				values[pos]=NOT_PRESENT;
				arrays[0][pos]=NOT_PRESENT;
				size--;
				if(value>limit){
					set(key, value);
				}else{
					sum++;
				}
			}
		}
		
		ArrayList<KmerNodeU> nodes=victims.toList();
		victims.clear();
		for(KmerNodeU node : nodes){
			int value=node.value();
			if(value<=limit){
				sum++;
			}else{
				kmer.setFrom(node.pivot());
				set(kmer, node.value());
			}
		}
		
		return sum;
	}
	
	public WalkerU1D walk(){
		return new WalkerU1D();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[] values;
	
	public int[] values(){return values;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Walker            ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO
	//Arrays are initialized to 0 instead of NOT_PRESENT
	
	/**
	 * Iterator for traversing all k-mer entries in the hash array.
	 * Walks through both main arrays and victim cache sequentially.
	 * Concurrent modification during iteration is not recommended.
	 */
	public class WalkerU1D extends WalkerU {
		
		WalkerU1D(){
			kmer=new Kmer(kbig);
			ha=HashArrayU1D.this;
			victims=ha.victims().toList();
		}
		
		public boolean next(){
			while(i<values.length && values[i]<=NOT_XPRESENT){i++;}
			if(i<values.length){
				fillKmer(i, kmer);
				value=values[i];
				assert(value!=NOT_XPRESENT);
				assert(kmer.len()>0) : kmer.len()+", "+value+", "+i+", "+values.length+"\n"
						+ "NOT_XPRESENT="+NOT_XPRESENT+", values[0]="+values[0]+", values[1]="+values[1];
				i++;
//				System.err.println("Y: "+kmer.len());
				return true;
			}
			if(i2<victims.size()){
				KmerNodeU kn=victims.get(i2);
				kn.fillKmer(kmer);
				value=kn.value();
				assert(value!=NOT_XPRESENT);
				i2++;
//				System.err.print("Z: "+kmer.len());
				return true;
			}
//			System.err.print("X: "+kmer.len());
			kmer.clearFast();
//			System.err.print("X2: "+kmer.len());
			value=NOT_XPRESENT;
			return false;
		}
		
		public Kmer kmer(){return kmer;}
		public int value(){return value;}
		
		private HashArrayU1D ha;
		private ArrayList<KmerNodeU> victims;
		
		private final Kmer kmer;
		private int value;

		/** Current index in main arrays during iteration; may point to an empty cell */
		private int i=0;
		private int i2=0;
	}
	
	//TODO: Remove after fixing array initialization
	private static final int NOT_XPRESENT=0;
}
