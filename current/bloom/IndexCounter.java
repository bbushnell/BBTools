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
 * Supports canonical and non-canonical counting and is optimized for large reference datasets.
 * @author Brian Bushnell
 * @date December 2, 2014
 */
public class IndexCounter extends KmerCountAbstract {
	
	/**
	 * Constructs an IndexCounter and initializes k-mer bit masks.
	 * @param k_ K-mer length (1-32)
	 * @param rcomp_ True to use canonical (max of forward/reverse) k-mers
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
	 * Creates and populates a new KCountArray by counting all reference chromosomes.
	 *
	 * @param cells Number of cells in the counting array
	 * @param cbits Bits per cell
	 * @param hashes Number of hash functions
	 * @return Populated and finalized KCountArray
	 */
	public KCountArray makeKcaFromIndex(long cells, int cbits, int hashes){
		//The 4-segment split double-counts each internal-boundary kmer and drops the last (length%4)
		//bases per chromosome (see IndexCounter#001). Correct ONLY for a presence bitmap (bits=1),
		//where the double-counts saturate; crash loud rather than emit inflated counts at bits>1.
		assert(cbits==1) : "IndexCounter is only correct for bits=1 (presence bitmap); see IndexCounter#001. cbits="+cbits;
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

	/**
	 * Populates an existing KCountArray by counting k-mers from all reference chromosomes in parallel.
	 * @param counts KCountArray to populate
	 * @return The populated KCountArray
	 * @throws Exception if worker threads fail
	 */
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
	
	private class CountThread extends Thread{
		
		/**
		 * Creates a counting thread with shared counting array and chromosome coordinator.
		 * @param counts_ Shared KCountArray
		 * @param nextChrom_ Atomic counter for assigning chromosome segments
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
		
		/** Processes assigned chromosome segments, counting k-mers into the shared array.
		 * @param counts Shared counting array */
		private final void count(KCountArray counts){
			assert(k>=1 && counts!=null);
			final int maxCount=THREADS_PER_CHROM*Data.numChroms;
			for(int cnum=nextChrom.getAndIncrement(); cnum<maxCount; cnum=nextChrom.getAndIncrement()){
				ChromosomeArray ca=Data.getChromosome(cnum/THREADS_PER_CHROM+1);
				processChrom(ca, cnum%THREADS_PER_CHROM);
			}
		}
		
		/**
		 * Counts k-mers for a chromosome segment using rolling hash and reverse complements.
		 * @param ca Chromosome data
		 * @param segNum Segment number (0-3) within the chromosome
		 */
		private final void processChrom(ChromosomeArray ca, int segNum){
			assert(k<=maxShortKmerLength);
			assert(CANONICAL);

			final byte[] bases=ca.array;
			if(bases==null || bases.length<k){return;}
			//TODO: Possible bug [bloom/IndexCounter#001] - the kmer ending at each internal segment
			//boundary (segNum*segLength-1) is counted by BOTH this segment (warmup) and segment segNum-1,
			//and segment 3 stops at 4*segLength so the last (length mod 4) bases' kmers are never counted.
			//Negligible in the only live caller (BBMap bloomfilter=t uses bits=1, so double-counts saturate).
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
		private final KCountArray counts;
		private final AtomicInteger nextChrom;
		private long keysCountedLocal=0;
		private long readsProcessedLocal=0;
	}
	
	private final int k;
//	private final int cbits;
	private final int shift;
	private final int shift2;
	private final long mask;
	private final boolean rcomp;
	
	private static final int THREADS_PER_CHROM=4;
	
}
