package ukmer;

import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import shared.Tools;
import structures.SuperLongList;

/**
 * Creates frequency histograms for k-mer tables with multi-threaded and single-threaded
 * processing options. Generates population frequency distribution for AbstractKmerTableU
 * collections by counting occurrences. Automatically selects between single-threaded and
 * multi-threaded execution based on thread availability.
 *
 * @author Brian Bushnell
 */
public final class HistogramMakerU {
	
	/**
	 * Creates a frequency histogram from k-mer tables using optimal threading strategy.
	 * Automatically chooses multi-threaded or single-threaded processing based on
	 * available threads.
	 * Uses multi-threading when more than 2 threads are available, otherwise uses
	 * single-threaded approach.
	 *
	 * @param tables Array of k-mer tables to process for histogram generation
	 * @param histMax Maximum histogram bucket size, determines array length
	 * @return Frequency histogram as long array where index represents count and
	 * value represents frequency
	 */
	public static long[] fillHistogram(final AbstractKmerTableU[] tables, final int histMax) {
		if(Shared.threads()>2){
			return fillHistogram_MT(tables, histMax);
		}else{
			return fillHistogram_ST(tables, histMax);
		}
	}
	
	/**
	 * Single-threaded histogram generation from k-mer tables.
	 * Sequentially processes each table and accumulates counts into histogram array.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Maximum histogram size
	 * @return Accumulated histogram array
	 */
	private static long[] fillHistogram_ST(final AbstractKmerTableU[] tables, final int histMax) {
		long[] ca=new long[histMax+1];
		for(AbstractKmerTableU set : tables){
			set.fillHistogram(ca, histMax);
		}
		return ca;
	}
	
	/**
	 * Multi-threaded histogram generation with dynamic thread allocation and load balancing.
	 * Calculates optimal thread count based on system threads, table count, and
	 * performance limits.
	 * Uses AtomicInteger for thread-safe work distribution and SuperLongList for
	 * per-thread accumulation.
	 *
	 * @param tables Array of k-mer tables to process in parallel
	 * @param histMax Maximum histogram bucket size
	 * @return Merged histogram array from all worker threads
	 */
	private static long[] fillHistogram_MT(final AbstractKmerTableU[] tables, final int histMax) {
		boolean errorState=false;
		int threads=Shared.threads();
		threads=Tools.min((threads>20 ? threads/2 : threads), (tables.length+1)/2, 32);
		if(threads<2){return fillHistogram_ST(tables, histMax);}
		
		final FillThread[] array=new FillThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<threads; i++){array[i]=new FillThread(tables, histMax, next);}
		for(int i=0; i<threads; i++){array[i].start();}
		
		//Wait for completion of all threads
		final long[] ca=new long[histMax+1];
		boolean success=true;
		for(FillThread pt : array){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			//Accumulate per-thread statistics
			
			pt.sll.addTo(ca);
			pt.sll=null;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		return ca;
	}
	
	/**
	 * Worker thread for parallel histogram generation from k-mer tables.
	 * Uses atomic work distribution to process tables and accumulates results in thread-local storage.
	 * Each thread processes tables until no more work is available.
	 */
	private static class FillThread extends Thread{
		
		/**
		 * Constructs a worker thread for histogram generation.
		 * Initializes thread-local SuperLongList with size optimization for efficient accumulation.
		 *
		 * @param tables_ Array of k-mer tables to process
		 * @param histMax_ Maximum histogram size for SuperLongList sizing
		 * @param next_ Atomic counter for thread-safe work distribution
		 */
		FillThread(final AbstractKmerTableU[] tables_, int histMax_, AtomicInteger next_){
			tables=tables_;
			next=next_;
			sll=new SuperLongList(Tools.mid(5000, histMax_, 100000));
		}
		
		@Override
		public void run(){
			for(int tnum=next.getAndIncrement(); tnum<tables.length; tnum=next.getAndIncrement()){
				tables[tnum].fillHistogram(sll);
			}
		}
		
		/** Array of k-mer tables to process for histogram generation */
		final AbstractKmerTableU[] tables;
		/** Atomic counter for thread-safe work distribution across k-mer tables */
		final AtomicInteger next;
		/** Thread-local storage for accumulating histogram counts during processing */
		SuperLongList sll;
		
	}
	
}
