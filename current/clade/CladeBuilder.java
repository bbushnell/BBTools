package clade;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import bin.AdjustEntropy;
import fileIO.FileFormat;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ListNum;
import template.ThreadWaiter;
import tracker.EntropyTracker;

/** Constructs one query Clade from intact contigs with private per-worker accumulators.
 * No taxonomy lookup, gene calling, minimum-length filter or sequence mutation.
 * Callers must keep the Clade/sketch configuration fixed until the call returns.
 * @author Ganyu */
public final class CladeBuilder {
	/** Retains the existing whole-file Read loader; parallelizes per-contig signature construction. */
	public static Clade fromFile(String filename,int threads){
		require(filename!=null && threads>0,"File and positive worker count required for Clade construction");
		return fromReads(ReadInputStream.toReads(filename,FileFormat.FASTA,-1),filename,threads);
	}
	/** Read container overload; list order defines entropy-reduction order, not numericID. */
	public static Clade fromReads(ListNum<Read> batch,String name,int threads){
		require(batch!=null && batch.list!=null,"A real Read batch is required");
		return fromReads(batch.list,name,threads);
	}
	/** Returns a finished taxID0/level0 query with the supplied name, matching serial construction.
	 * Each listed Read contributes its bases once; mates are not implicitly traversed. */
	public static Clade fromReads(List<Read> reads,String name,int threads){
		require(reads!=null && threads>0,"Read list and positive worker count required");
		final Read[] input=reads.toArray(new Read[reads.size()]);
		for(Read read:input){require(read!=null && read.bases!=null,"Every listed contig needs bases");}
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150){AdjustEntropy.load(4,150);}
		final Clade result=new Clade(0,0,name);
		if(threads==1 || input.length<2){
			final EntropyTracker entropy=tracker();
			for(Read read:input){result.add(read.bases,entropy);}
			result.finish();return result;
		}
		final int count=Math.min(threads,input.length);
		final AtomicInteger next=new AtomicInteger();
		final float[] entropy=new float[input.length];
		final ArrayList<Worker> workers=new ArrayList<Worker>(count);
		for(int i=0;i<count;i++){workers.add(new Worker(input,entropy,next));}
		ThreadWaiter.startAndWait(workers);
		for(Worker worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Clade construction worker failed",worker.failure);}
			require(worker.success,"Clade worker did not finish; partial signatures cannot be used");
		}
		for(Worker worker:workers){
			final Clade partial=worker.partial;
			Tools.add(result.counts,partial.counts);
			result.incrementBases(partial.bases);result.contigs+=partial.contigs;result.addDDL(partial.ddl);
		}
		// Reproduce Clade.add's float arithmetic in original input order. Merging worker
		// weighted means would change rounding even though the integer spectra are exact.
		long bases=0;float mean=0;
		for(int i=0;i<input.length;i++){
			final int length=input[i].bases.length;
			mean=(mean*bases+entropy[i]*length)/(float)(bases+length);bases+=length;
		}
		result.entropy=mean;
		require(result.bases==bases && result.contigs==input.length,"Private accumulator merge must conserve all input contigs and bases");
		result.finish();return result;
	}
	private static EntropyTracker tracker(){
		return CladeObject.calcCladeEntropy ? new EntropyTracker(CladeObject.entropyK,CladeObject.entropyWindow,false) : null;
	}
	private static final class Worker extends Thread{
		Worker(Read[] input_,float[] entropy_,AtomicInteger next_){input=input_;entropy=entropy_;next=next_;}
		@Override public void run(){
			try{
				final EntropyTracker tracker=tracker();
				for(int i=next.getAndIncrement();i<input.length;i=next.getAndIncrement()){
					entropy[i]=partial.addSequence(input[i].bases,tracker);
				}
				success=true;
			}catch(Throwable e){failure=e;}
		}
		final Clade partial=new Clade(0,0,null);
		final Read[] input;
		final float[] entropy;
		final AtomicInteger next;
		boolean success;
		Throwable failure;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private CladeBuilder(){}
}
