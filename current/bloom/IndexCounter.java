package bloom;

import java.lang.Thread.State;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;
import shared.Tools;

/**
 * Counts k-mers from indexed reference data using multi-threaded chromosome processing.
 * Optimized for processing large genomic reference databases by parallelizing across
 * chromosome segments. Supports canonical and non-canonical k-mer counting modes.
 *
 * @author Brian Bushnell
 * @date December 2, 2014
 */
public class IndexCounter extends KmerCountAbstract {
	
	/**
	 * Constructs an IndexCounter with specified k-mer parameters.
	 * Initializes bit manipulation constants for efficient k-mer encoding.
	 * @param k_ K-mer length (must be 1-32)
	 * @param rcomp_ True to use canonical k-mers (max of forward/reverse)
	 */
	public IndexCounter(final int k_, final boolean rcomp_){
		k=k_;
		rcomp=rcomp_;

		final int bitsPerChar=2;
		shift=bitsPerChar*k;
		shift2=shift-bitsPerChar;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		assert(k>=1 && k<33) : k;
	}
	
	/**
	 * Creates a new KCountArray and populates it by counting k-mers from the reference index.
	 * Convenience method that creates, populates, and finalizes the counting array.
	 *
	 * @param cells Number of cells in the counting array
	 * @param cbits Bits per cell for count storage
	 * @param hashes Number of hash functions to use
	 * @return Populated and finalized KCountArray
	 */
	public KCountArray makeKcaFromIndex(long cells, int cbits, int hashes){
		KCountArray kca=KCountArray.makeNew(cells, cbits, hashes, null, 0);
		try {
			countFromIndex(kca);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		kca.shutdown();
		return kca;
	}

	public KCountArray countFromIndex(KCountArray counts) throws Exception{
		
		final CountThread[] cta=new CountThread[Tools.min(Data.numChroms*THREADS_PER_CHROM, Shared.threads())];
		final AtomicInteger nextChrom=new AtomicInteger(0);
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(counts, nextChrom);
			cta[i].start();
		}
//		System.out.println("~1");
		for(int i=0; i<cta.length; i++){
//			System.out.println("~2");
			CountThread ct=cta[i];
			synchronized(ct){
//				System.out.println("~3");
				while(ct.getState()!=State.TERMINATED){
//					System.out.println("~4");
					try {
						ct.join(2000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
//					System.out.println("~5");
				}
			}
		}
		
		return counts;
	}
	
	/**
	 * Worker thread for parallel k-mer counting from reference chromosomes.
	 * Each thread processes multiple chromosome segments and maintains local counters
	 * before synchronizing with global statistics.
	 */
	private class CountThread extends Thread{
		
		/**
		 * Constructs a CountThread with shared counting array and chromosome coordinator.
		 * @param counts_ Shared KCountArray for k-mer counting
		 * @param nextChrom_ Atomic counter for thread-safe chromosome assignment
		 */
		CountThread(final KCountArray counts_, AtomicInteger nextChrom_){
			counts=counts_;
			nextChrom=nextChrom_;
		}
		
		@Override
		public void run(){
			count(counts);
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}
		
		/**
		 * Main counting logic that processes chromosomes assigned to this thread.
		 * Uses atomic counter to get next chromosome segment and processes it.
		 * Continues until all chromosome segments have been processed.
		 * @param counts KCountArray to increment with discovered k-mers
		 */
		private final void count(KCountArray counts){
			assert(k>=1 && counts!=null);
			final int maxCount=THREADS_PER_CHROM*Data.numChroms;
			for(int cnum=nextChrom.getAndIncrement(); cnum<maxCount; cnum=nextChrom.getAndIncrement()){
				ChromosomeArray ca=Data.getChromosome(cnum/THREADS_PER_CHROM+1);
				processChrom(ca, cnum%THREADS_PER_CHROM);
			}
		}
		
		/**
		 * Processes a specific segment of a chromosome for k-mer counting.
		 * Uses rolling hash to efficiently generate k-mers and their reverse complements.
		 * Handles canonical k-mer selection and increments the counting array.
		 *
		 * @param ca ChromosomeArray containing the sequence data
		 * @param segNum Segment number within the chromosome (0-3)
		 */
		private final void processChrom(ChromosomeArray ca, int segNum){
			assert(k<=maxShortKmerLength);
			assert(CANONICAL);

			final byte[] bases=ca.array;
			if(bases==null || bases.length<k){return;}
			final int segLength=bases.length/4;
			final int start=Tools.max(0, segNum*segLength-k);
			final int stop=Tools.min(bases.length, (segNum+1)*segLength);
			
			long kmer=0;
			long rkmer=0;
			int len=0;

			for(int i=start; i<stop; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						counts.increment(key);
						readsProcessedLocal++;
					}
				}
			}
		}
		/** Shared counting array for accumulating k-mer frequencies */
		private final KCountArray counts;
		/** Thread-safe counter for assigning chromosome segments to workers */
		private final AtomicInteger nextChrom;
		/** Thread-local counter for unique k-mers processed by this thread */
		private long keysCountedLocal=0;
		/** Thread-local counter for total k-mer instances processed by this thread */
		private long readsProcessedLocal=0;
	}
	
	/** K-mer length for counting operations */
	private final int k;
//	private final int cbits;
	/** Bit shift value for rolling hash operations (2 * k) */
	private final int shift;
	/** Bit shift value for reverse complement operations (shift - 2) */
	private final int shift2;
	/** Bit mask for k-mer extraction (handles k up to 32) */
	private final long mask;
	/** True to use canonical k-mers (maximum of forward and reverse complement) */
	private final boolean rcomp;
	
	/** Number of worker threads to assign per chromosome for parallel processing */
	private static final int THREADS_PER_CHROM=4;
	
}
