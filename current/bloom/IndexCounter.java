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
 * Supports canonical and non-canonical counting across partitioned chromosome arrays.
 * Loaded reference data and shared counting settings must remain stable while workers run.
 * The legacy readsProcessed tally counts accepted k-mer windows here, not input reads.
 * Failed-call tables and counters may be partial and must not be treated as completed results.
 * @author Brian Bushnell
 * @date December 2, 2014
 */
public class IndexCounter extends KmerCountAbstract{
	
	/** Builds a counter for short indexed-reference kmers. */
	public IndexCounter(final int k_, final boolean rcomp_){
		k=k_;
		rcomp=rcomp_;

		final int bitsPerChar=2;
		shift=bitsPerChar*k;
		shift2=shift-bitsPerChar;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		assert(k>=1 && k<33) : k;
	}
	
	/** Counts the loaded reference index into a finalized presence bitmap.
	 * Counting failures propagate instead of returning an incomplete table. */
	public KCountArray makeKcaFromIndex(long cells, int cbits, int hashes){
		//Retain the presence-only factory contract; caller-owned count tables use countFromIndex.
		assert(cbits==1) : "makeKcaFromIndex requires a presence bitmap (cbits=1); cbits="+cbits;
		KCountArray kca=KCountArray.makeNew(cells, cbits, hashes, null, 0);
		try{
			countFromIndex(kca);
		}catch(Exception e){
			throw new RuntimeException("Indexed reference counting failed; refusing to return incomplete counts.", e);
		}
		kca.shutdown();
		return kca;
	}

	/** Starts per-segment workers and, after successful startup, joins every worker.
	 * Worker failures propagate with their causes after joins; no table rollback is provided.
	 * Interrupted joins are retried without restoring the interrupt flag.
	 * @param counts Borrowed table supporting concurrent increments, or used with one worker
	 * @return The supplied table after all workers finish successfully
	 * @throws RuntimeException If a worker fails; discard the partial table and counters */
	public KCountArray countFromIndex(KCountArray counts) throws Exception{
		
		final CountThread[] cta=new CountThread[Tools.min(Data.numChroms*THREADS_PER_CHROM, Shared.threads())];
		final AtomicInteger nextChrom=new AtomicInteger(0);
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(counts, nextChrom);
			cta[i].start();
		}
		RuntimeException failure=null;
//		System.out.println("~1");
		for(int i=0; i<cta.length; i++){
//			System.out.println("~2");
			CountThread ct=cta[i];
			synchronized(ct){
//				System.out.println("~3");
				while(ct.getState()!=State.TERMINATED){
//					System.out.println("~4");
					try{
						ct.join(2000);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
//					System.out.println("~5");
				}
			}
			if(ct.failure!=null){
				if(failure==null){
					failure=new RuntimeException("Indexed reference worker failed; counts are incomplete.", ct.failure);
				}else{failure.addSuppressed(ct.failure);}
			}
		}
		
		if(failure!=null){throw failure;}
		return counts;
	}
	
	private class CountThread extends Thread{
		
		/** Binds this worker to shared output and a segment counter. */
		CountThread(final KCountArray counts_, AtomicInteger nextChrom_){
			counts=counts_;
			nextChrom=nextChrom_;
		}
		
		/** Counts assigned segments and folds local tallies into the shared totals. */
		@Override
		public void run(){
			try{
				count(counts);
				synchronized(getClass()){
					keysCounted+=keysCountedLocal;
					readsProcessed+=readsProcessedLocal;
					if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
					if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
				}
			}catch(Throwable t){failure=t;}
		}
		
		/** Processes all chromosome segments claimed by this worker. */
		private final void count(KCountArray counts){
			assert(k>=1 && counts!=null);
			final int maxCount=THREADS_PER_CHROM*Data.numChroms;
			for(int cnum=nextChrom.getAndIncrement(); cnum<maxCount; cnum=nextChrom.getAndIncrement()){
				ChromosomeArray ca=Data.getChromosome(cnum/THREADS_PER_CHROM+1);
				processChrom(ca, cnum%THREADS_PER_CHROM);
			}
		}
		
		/** Counts windows ending in this segment, warming up with at most k-1 preceding bases.
		 * The last segment includes the division remainder; each valid window is counted once. */
		private final void processChrom(ChromosomeArray ca, int segNum){
			assert(k<=maxShortKmerLength);
			assert(CANONICAL);

			final byte[] bases=ca.array;
			if(bases==null || bases.length<k){return;}
			//IndexCounter#001: k-1 warmup bases avoid duplicate boundary ends; retain the terminal remainder.
			final int segLength=bases.length/THREADS_PER_CHROM;
			final int start=Tools.max(0, segNum*segLength-(k-1));
			final int stop=(segNum==THREADS_PER_CHROM-1 ? bases.length : (segNum+1)*segLength);
			
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
		/** Failure published to the joining caller; null after successful completion. */
		private volatile Throwable failure;
	}
	
	private final int k;
//	private final int cbits;
	private final int shift;
	private final int shift2;
	private final long mask;
	private final boolean rcomp;
	
	private static final int THREADS_PER_CHROM=4;
	
}
