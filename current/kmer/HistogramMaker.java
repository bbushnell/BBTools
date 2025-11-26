package kmer;

import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import shared.Tools;
import structures.SuperLongList;

/**
 * Creates k-mer frequency histograms from hash tables with multi-threaded and
 * single-threaded processing strategies. Generates frequency distribution of
 * k-mer occurrences across multiple hash tables for genomic analysis.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public final class HistogramMaker {

	
	/**
	 * Fills a histogram array with k-mer frequency counts from the provided tables.
	 * Automatically chooses between single-threaded and multi-threaded processing
	 * based on available CPU cores.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Maximum histogram bin value
	 * @return Array containing frequency counts for each histogram bin
	 */
	public static long[] fillHistogram(final AbstractKmerTable[] tables, final int histMax) {
		if(Shared.threads()>2){
			return fillHistogram_MT(tables, histMax);
		}else{
			return fillHistogram_ST(tables, histMax);
		}
	}
	
	/**
	 * Single-threaded histogram filling implementation.
	 * Processes tables sequentially to generate frequency distribution.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Maximum histogram bin value
	 * @return Array containing frequency counts for each histogram bin
	 */
	private static long[] fillHistogram_ST(final AbstractKmerTable[] tables, final int histMax) {
		long[] ca=new long[histMax+1];
		for(AbstractKmerTable set : tables){
			set.fillHistogram(ca, histMax);
		}
		return ca;
	}
	
	/**
	 * Multi-threaded histogram filling implementation.
	 * Distributes table processing across worker threads and accumulates results.
	 * Automatically adjusts thread count based on available cores and table count.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Maximum histogram bin value
	 * @return Array containing accumulated frequency counts from all threads
	 */
	private static long[] fillHistogram_MT(final AbstractKmerTable[] tables, final int histMax) {
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
	 * Worker thread for parallel histogram generation.
	 * Processes assigned k-mer tables and accumulates histogram data
	 * in a thread-local SuperLongList for later aggregation.
	 */
	private static class FillThread extends Thread{
		
		/**
		 * Constructs a FillThread with shared resources.
		 * @param tables_ Shared array of k-mer tables to process
		 * @param histMax_ Maximum histogram bin value
		 * @param next_ Atomic counter for work distribution among threads
		 */
		FillThread(final AbstractKmerTable[] tables_, int histMax_, AtomicInteger next_){
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
		
		/** Shared array of k-mer tables to process */
		final AbstractKmerTable[] tables;
		/** Atomic counter for thread-safe work distribution */
		final AtomicInteger next;
		/** Thread-local histogram accumulator for collecting frequency data */
		SuperLongList sll;
		
	}
	
}
