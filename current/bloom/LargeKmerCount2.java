package bloom;

import java.util.ArrayList;
import shared.Random;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Legacy experiment combining forward rolling-hash counts with collision estimates.
 * Both mates contribute valid windows to one packed, saturating table. Repeated
 * windows and distinct keys colliding at one address are indistinguishable in the
 * table; the additional counters are heuristic estimates, not exact cardinalities.
 * Undefined bases break windows. Quality values and minQuality do not affect
 * counting. Per-table statistics reset on entry; readsProcessed accumulates roots.
 * Configuration/statistics are process-wide, so concurrent invocations are unsupported.
 *
 * @author Brian Bushnell
 * @date Jul 6, 2012
 */
public class LargeKmerCount2{
	
	/** Counts input and prints table occupancy plus experimental collision estimates.
	 * Empty tables report zero singleton/useful percentages.
	 * Sets the global FASTA minimum read length when the primary input has a FASTA suffix.
	 * @param args Input, optional mate input, index bits, cell bits and k */
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>4 || args[1].contains(".") ? args[1] : null);
		int indexbits=Integer.parseInt(args[args.length-3]);
		int cbits=Integer.parseInt(args[args.length-2]);
		int k=Integer.parseInt(args[args.length-1]);
		
		KCountArray2 count=null;
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		count=countFastq(fname1, fname2, indexbits, cbits, k);
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		long[] freq=count.transformToFrequency();

		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ", (100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			long x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ", (100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long estKmers=load+min(actualCollisions, (long)expectedCollisions);
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("probCollisions:\t         \t"+(long)probNewKeyCollisions);
		System.out.println("EstimateP:     \t         \t"+(sum2+(long)probNewKeyCollisions));
		System.out.println("expectedColl:  \t         \t"+(long)expectedCollisions);
		System.out.println("actualColl:    \t         \t"+(long)actualCollisions);
		System.out.println("estimateKmers: \t         \t"+estKmers);
		System.out.println();
		System.out.println("Singleton:     \t"+Tools.format("%.3f%%   ", (sum2==0 ? 0 : 100d*x/sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+Tools.format("%.3f%%   ", (sum2==0 ? 0 : 100d*x/sum2))+"\t"+x);
		
	}
	
	/**
	 * Builds a new forward-hash count table, consuming reads without changing them.
	 * maxReads is passed to CRIS: FASTQ roots/pairs or SAM nonheader records.
	 * readsProcessed accumulates roots, not individual mates.
	 * Resets per-table statistics. An undefined collision estimate stops counting;
	 * remaining input batches are drained and returned before closing and throwing.
	 * This avoids early cancellation of the shared reader; it can require reading
	 * the remaining input. Failed-call statistics are partial and must be discarded.
	 * Other counting/startup exceptions have no general finally cleanup of the owned reader.
	 * @param reads1 Primary sequence input, autodetected with FASTQ as fallback
	 * @param reads2 Optional separate mate input
	 * @param indexbits Log2 of allocated cells; at least two gives a positive primary modulus
	 * @param cbits Cell width supported by KCountArray2
	 * @param k Positive rolling-window length
	 * @return New saturated count table; occupancy is not an exact unique-key count
	 * @throws IllegalStateException If the estimator is undefined or input closure reports failure
	 */
	public static KCountArray2 countFastq(String reads1, String reads2, int indexbits, int cbits, int k){
		//TODO: Probable bug - k<=0 is not rejected before ring allocation/indexing.
		assert(indexbits>=1 && indexbits<40);
		final long cells=1L<<indexbits;
		final int kbits=ROTATE_DIST*k;
		final int xorShift=kbits%64;
		final long[] rotMasks=makeRotMasks(xorShift);
		final int[] buffer=new int[k];
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", indexbits="+indexbits+", cells="+cells+", cbits="+cbits);}
		if(verbose){System.err.println("xorShift="+xorShift+", rotMasks[3]="+Long.toHexString(rotMasks[3]));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
		expectedCollisions=0;
		actualCollisions=0;
		keysCounted=0;
		load=0;
		probNewKeyCollisions=0;
		invCells=1d/cells;
		invKmerSpace=Math.pow(0.5, 2*k);
		if(cells>=Math.pow(4, k)){invCells=0;}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cris.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		long kmer=0; //current kmer
		int len=0;  //distance since last contig start or ambiguous base
		IllegalStateException estimatorFailure=null;
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(Read r : reads){
					if(estimatorFailure!=null){break;}//Drain subsequent batches without counting.
					readsProcessed++;
					
					len=0;
					kmer=0;
					byte[] bases=r.bases;
					byte[] quals=r.quality;
					
					for(int i=0; i<bases.length; i++){
						byte b=bases[i];
						int x=AminoAcid.baseToNumber[b];
						int x2=buffer[len%buffer.length];
						buffer[len%buffer.length]=x;
						if(x<0){
							len=0;
							kmer=0;
						}else{
							kmer=(Long.rotateLeft(kmer, ROTATE_DIST)^x);
							len++;
							if(len>=k){
								keysCounted++;
								if(len>k){kmer=kmer^rotMasks[x2];}
								long hashcode=kmer&0x7fffffffffffffffL;
//								hashcode=randy.nextLong()&~((-1L)<<(2*k));
								long code1=hashcode%(cells-3);
//								long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
								int value=count.increment2(code1, 1);
								
								estimatorFailure=updateCollisionEstimate(value);
								if(estimatorFailure!=null){break;}
									
//									double estKmers=load+min(actualCollisions, expectedCollisions);
//									double probOldKmer=estKmers*invKmerSpace;
//									probNewKeyCollisions+=(prob*(1-prob2));
								
////								probCollisions+=(load*invCells);
//								if(value==0){load++;}
//								else{
////									long load2=keysCounted-load;
//									double prob=Math.sqrt(load*invCells);
//									double estKmers=load+probNewKeyCollisions;
//									double prob2=estKmers*invKmerSpace;
////									probCollisions+=(prob*(1-prob2));
////									probCollisions+=Math.sqrt(prob*(1-prob2));
//									probNewKeyCollisions+=Math.sqrt(prob*(1-prob2));
////									probCollisions+=min(prob, 1-prob2);
////									probCollisions+=(load*invCells);
//								}
							}
						}
					}
					
					if(r.mate!=null && estimatorFailure==null){
						len=0;
						kmer=0;
						bases=r.mate.bases;
						quals=r.mate.quality;
						for(int i=0; i<bases.length; i++){
							byte b=bases[i];
							int x=AminoAcid.baseToNumber[b];
							int x2=buffer[len%buffer.length];
							buffer[len%buffer.length]=x;
							if(x<0){
								len=0;
								kmer=0;
							}else{
								kmer=(Long.rotateLeft(kmer, ROTATE_DIST)^x);
								len++;
								if(len>=k){
									keysCounted++;
									if(len>k){kmer=kmer^rotMasks[x2];}
									long hashcode=kmer&0x7fffffffffffffffL;
//									hashcode=randy.nextLong()&~((-1L)<<(2*k));
									long code1=hashcode%(cells-3);
//									long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
									int value=count.increment2(code1, 1);
									
									estimatorFailure=updateCollisionEstimate(value);
									if(estimatorFailure!=null){break;}
										
//										double estKmers=load+min(actualCollisions, expectedCollisions);
//										double probOldKmer=estKmers*invKmerSpace;
//										probNewKeyCollisions+=(prob*(1-prob2));
									
////									probCollisions+=(load*invCells);
//									if(value==0){load++;}
//									else{
////										long load2=keysCounted-load;
//										double prob=Math.sqrt(load*invCells);
//										double estKmers=load+probNewKeyCollisions;
//										double prob2=estKmers*invKmerSpace;
////										probCollisions+=(prob*(1-prob2));
////										probCollisions+=Math.sqrt(prob*(1-prob2));
//										probNewKeyCollisions+=Math.sqrt(prob*(1-prob2));
////										probCollisions+=min(prob, 1-prob2);
////										probCollisions+=(load*invCells);
//									}
								}
							}
						}
					}
					
				}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			System.err.println("Finished reading");
			cris.returnList(ln);
			System.err.println("Returned list");
			final boolean errorState=ReadWrite.closeStream(cris);
			System.err.println("Closed stream");
			System.err.println("Processed "+readsProcessed+" reads.");
			if(errorState){
				final IllegalStateException inputFailure=new IllegalStateException("Input/close failure; large-kmer counts are incomplete: "+reads1);
				if(estimatorFailure==null){throw inputFailure;}
				estimatorFailure.addSuppressed(inputFailure);
			}
		}
		
		if(estimatorFailure!=null){throw estimatorFailure;}
		return count;
	}
	
	/** Updates the original collision heuristic for one increment, or returns its failure.
	 * Stages floating-point results before publication so undefined values never enter
	 * shared estimate fields. The caller drains input before throwing a returned failure.
	 * @param previousCount Nonnegative cell value returned by KCountArray2.increment2
	 * @return Null on success; diagnostic failure when the estimate is undefined */
	private static IllegalStateException updateCollisionEstimate(final int previousCount){
		assert(previousCount>=0) : "KCountArray2.increment2 returns a nonnegative prior cell count: "+previousCount;
		final double probCollision=load*invCells;
		final double expected=expectedCollisions+probCollision*(1-(load+min(expectedCollisions, actualCollisions))*invKmerSpace);
		if(!Double.isFinite(expected)){
			return new IllegalStateException("Nonfinite collision estimate: expected="+expected+", load="+load);
		}
		if(previousCount==0){
			expectedCollisions=expected;
			load++;
			return null;
		}
		actualCollisions++;
		final double denominator=min(expected, actualCollisions);
		if(!(denominator>0) || !Double.isFinite(denominator)){
			return new IllegalStateException("Undefined collision estimator: denominator must be positive and finite; denominator="+
					denominator+", expected="+expected+", actual="+actualCollisions+", load="+load+", invCells="+invCells);
		}
		final double probNewKey=(load*invCells)*expected/denominator;
		final double estKeys=load+probNewKeyCollisions;
		final double probOldKey=estKeys*invKmerSpace;
		final double updated=probNewKeyCollisions+probNewKey*(1-probOldKey);
		if(!Double.isFinite(updated)){
			return new IllegalStateException("Nonfinite collision estimate: newKeyCollisions="+updated+", expected="+expected);
		}
		expectedCollisions=expected;
		probNewKeyCollisions=updated;
		return null;
	}

	/** Returns rotated A/C/G/T codes for removing the departing base from the hash.
	 * Rotation follows Long.rotateLeft's modulo-64 distance semantics. */
	public static final long[] makeRotMasks(int rotDist){
		long[] masks=new long[4];
		for(long i=0; i<4; i++){
			masks[(int)i]=Long.rotateLeft(i, rotDist);
		}
		return masks;
	}
	
	/** Returns a new 2000-bin histogram of nonnegative values, clipping indices at 1999. */
	public static long[] transformToFrequency(int[] count){
		long[] freq=new long[2000];
		int max=freq.length-1;
		for(int i=0; i<count.length; i++){
			int x=count[i];
			x=min(x, max);
			freq[x]++;
		}
		return freq;
	}
	
	/** Sums a non-null integer array in a long accumulator. */
	public static long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	/** Sums a non-null long array with ordinary signed-long overflow semantics. */
	public static long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	/** Returns the smaller integer. */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger integer. */
	public static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller long. */
	public static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger long. */
	public static final long max(long x, long y){return x>y ? x : y;}
	/** Returns the smaller non-NaN value; otherwise follows the ordinary comparison expression. */
	public static final double min(double x, double y){return x<y ? x : y;}
	/** Returns the larger non-NaN value; otherwise follows the ordinary comparison expression. */
	public static final double max(double x, double y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Shared State          ----------------*/
	/*--------------------------------------------------------------*/

	/** Print table geometry during construction. */
	public static boolean verbose=true;
	/** Retained legacy setting; no active loop consults it. */
	public static byte minQuality=-5;
	/** Accumulated input roots; mates are processed within their root. */
	public static long readsProcessed=0;
	/** Reader-dependent limit: FASTQ roots/pairs or SAM nonheader records; negative is unlimited. */
	public static long maxReads=10000000L;
	/** Hash rotation in bits for each incorporated base. */
	public static final int ROTATE_DIST=2;
	
	/** Number of primary cells first occupied in the current invocation. */
	public static long load;
	/** Model-based collision accumulator for the current table. */
	public static double expectedCollisions;
	/** Increments hitting a previously occupied primary cell, including true repeats. */
	public static long actualCollisions;
	/** Heuristic new-key collision estimate, reset when counting starts. */
	public static double probNewKeyCollisions;
	/** Reciprocal allocated cells, set to zero when cells cover the nominal 4^k space. */
	public static double invCells;
	/** Reciprocal nominal nucleotide key space, 4^-k. */
	public static double invKmerSpace;
	/** Valid window occurrences from both mates in the current table. */
	public static long keysCounted;
	
	/** Seeded generator retained for the commented-out random-key experiment. */
	public static final Random randy=shared.Shared.random(1);
	
}
