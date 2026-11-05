package kmer;

import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import shared.Tools;
import structures.SuperLongList;

/**
 * Aggregates the scalar-value histograms supplied by stable k-mer tables.
 * In count-valued tables, bin c counts distinct keys having count c; values at
 * or above histMax accumulate in the last bin. List-valued backends define their
 * own scalar view. The caller must exclude table mutations during traversal.
 * Worker failures are propagated after all workers terminate, without returning
 * a partial histogram.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public final class HistogramMaker{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Allocates a histogram and selects direct or worker processing from the configured
	 * thread count and table count. Input table references must be non-null and stable.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Nonnegative final bin index; histMax+1 must be an allocatable array length
	 * @return New histogram of length histMax+1, including overflow in its final bin
	 */
	public static long[] fillHistogram(final AbstractKmerTable[] tables, final int histMax){
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
	private static long[] fillHistogram_ST(final AbstractKmerTable[] tables, final int histMax){
		long[] ca=new long[histMax+1];
		for(AbstractKmerTable set : tables){
			set.fillHistogram(ca, histMax);
		}
		return ca;
	}

	/**
	 * Uses private accumulators and reduces them only after every worker succeeds.
	 * Caps workers at 32 and half the table count rounded up; configured counts above 20
	 * are halved first. Falls back to direct traversal if fewer than two workers result.
	 * Interrupted joins are logged and retried.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param histMax Maximum histogram bin value
	 * @return Array containing accumulated frequency counts from all threads
	 */
	private static long[] fillHistogram_MT(final AbstractKmerTable[] tables, final int histMax){
		int threads=Shared.threads();
		threads=Tools.min((threads>20 ? threads/2 : threads), (tables.length+1)/2, 32);
		if(threads<2){return fillHistogram_ST(tables, histMax);}

		final FillThread[] array=new FillThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<threads; i++){array[i]=new FillThread(tables, histMax, next);}
		for(int i=0; i<threads; i++){array[i].start();}

		//Join every worker before propagating a failure or reducing its private accumulator.
		final long[] ca=new long[histMax+1];
		Throwable failure=null;
		for(FillThread pt : array){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					pt.join();
				}catch(InterruptedException e){
					//Keep waiting so every owned accumulator is quiescent before reduction.
					e.printStackTrace();
				}
			}

			if(failure==null){failure=pt.failure;}
		}
		rethrowFailure(failure);
		for(FillThread pt : array){
			pt.sll.addTo(ca);
			pt.sll=null;
		}
		return ca;
	}

	/** Preserves the first worker exception/error after every owned worker has terminated. */
	private static void rethrowFailure(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException("K-mer histogram worker failed.", failure);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker thread for parallel histogram generation.
	 * Processes assigned k-mer tables and accumulates histogram data
	 * in a thread-local SuperLongList for later aggregation.
	 */
	private static class FillThread extends Thread{

		/**
		 * Borrows the tables/counter and allocates a private accumulator; does not start.
		 * @param tables_ Shared array of k-mer tables to process
		 * @param histMax_ Dense-accumulator sizing hint, clamped to 5000..100000
		 * @param next_ Atomic counter for work distribution among threads
		 */
		FillThread(final AbstractKmerTable[] tables_, int histMax_, AtomicInteger next_){
			tables=tables_;
			next=next_;
			sll=new SuperLongList(Tools.mid(5000, histMax_, 100000));
		}

		/** Claims tables and accumulates their values, retaining any failure for the caller. */
		@Override
		public void run(){
			try{
				for(int tnum=next.getAndIncrement(); tnum<tables.length; tnum=next.getAndIncrement()){
					tables[tnum].fillHistogram(sll);
				}
			}catch(Throwable t){failure=t;}
		}

		/** Shared array of k-mer tables to process */
		final AbstractKmerTable[] tables;
		/** Atomic counter for thread-safe work distribution */
		final AtomicInteger next;
		/** Thread-local histogram accumulator for collecting frequency data */
		SuperLongList sll;
		/** Failure published to the coordinator instead of returning a partial histogram. */
		volatile Throwable failure;

	}

}
