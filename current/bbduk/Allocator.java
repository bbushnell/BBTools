package bbduk;

import map.LongIntMapX;
import shared.Shared;
import shared.Tools;

public class Allocator{

	/**
	 * This allocates the data structures in multiple threads.  Unfortunately, it does not lead to any speedup, at least for ARRAY type.
	 * @param ways
	 * @param tableType
	 * @param schedule
	 * @param mask
	 * @return The preallocated table
	 */
	public static final LongIntMapX[] preallocate(int ways, int size){

		final LongIntMapX[] tables=new LongIntMapX[ways];
		
		{
			shared.Timer tm=new shared.Timer();
			final int t=Tools.max(1, Tools.min(Shared.threads(), 2, ways)); //More than 2 still improves allocation time, but only slightly; ~25% faster at t=4.
			final AllocThread[] allocators=new AllocThread[t];
			for(int i=0; i<t; i++){
				allocators[i]=new AllocThread(size, tables, i, t);
			}
			for(AllocThread at : allocators){at.start();}
			for(AllocThread at : allocators){
				while(at.getState()!=Thread.State.TERMINATED){
					try {
						at.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			tm.stop();
		}
		
		synchronized(tables){
			for(int i=0; i<tables.length; i++){
				final LongIntMapX akt=tables[i];
				if(akt==null){
					throw new RuntimeException("KmerTable allocation failed, probably due to lack of RAM: "+i+", "+tables.length);
				}
			}
		}
		
		return tables;
	}
	
	/**
	 * Creates a bitmask that zeros out the middle symbol(s) of a k-mer.
	 * Used for fuzzy k-mer matching by ignoring the center base(s).
	 *
	 * @param k K-mer length
	 * @param amino true for amino acid sequences, false for nucleotides
	 * @return Bitmask with middle bits zeroed
	 */
	public static final long makeMiddleMask(int k, boolean amino){
		final boolean odd=((k&1)==1);
		final int bitsPerSymbol=(amino ? 5 : 2);
		final int shift=bitsPerSymbol*((k-1)/2);
		final long middle=~((-1L)<<(odd ? bitsPerSymbol : 2*bitsPerSymbol));
		return ~(middle<<shift);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for parallel k-mer table allocation.
	 * Distributes table creation across multiple threads for faster initialization.
	 * @author Brian Bushnell
	 * @date Oct 23, 2013
	 */
	private static class AllocThread extends Thread{
		
		/**
		 * Creates allocation thread for specific table subset.
		 *
		 * @param type_ Table type to create
		 * @param schedule_ Size schedule for table growth
		 * @param mod_ Starting index modulo for this thread
		 * @param div_ Thread count divisor for work distribution
		 * @param mask_ Bit mask for table operations
		 * @param tables_ Shared array to store created tables
		 */
		AllocThread(int size_, LongIntMapX[] tables_, int mod_, int div_){
			size=Integer.highestOneBit(size_);
			mod=mod_;
			div=div_;
			tables=tables_;
		}
		
		@Override
		public void run(){
			//Initialize tables
			
//			Shared.printMemory();}
			for(int i=mod; i<tables.length; i+=div){
//				System.err.println("T"+i+" allocating "+i);
				final LongIntMapX akt=new LongIntMapX(size);
				synchronized(tables){
					tables[i]=akt;
				}
			}
		}
		
		/** Initial table size from schedule */
		private final int size;
		/** Starting index modulo for this thread */
		private final int mod;
		/** Thread count divisor for work distribution */
		private final int div;
		/** Shared array to store created tables */
		final LongIntMapX[] tables;
		
	}
	
}
