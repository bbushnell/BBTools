package test.bloom;

import java.util.ArrayList;
import java.util.Arrays;

import bloom.KCountArray2;
import bloom.KmerCount3;
import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentGenericReadInputStream;
import stream.FastqReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Legacy experiment estimating long-window support from shorter forward k-mers.
 * For each quality-qualified k2-base window, bins the minimum constituent k-mer
 * count, clipped to 255. Counts window occurrences rather than distinct long keys;
 * the source table can also saturate. This is not an exact large-k-mer counter.
 * Input is FASTQ with non-null qualities and shared KmerCount3 limits/quality settings apply.
 *
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public class TestLargeKmer{
	
	/** Builds short counts, scans input again and prints all 256 support bins.
	 * @param args Input, optional mate input, k, cell bits and k2, in that order */
	public static void main(String[] args){
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>4 || args[1].contains(".") ? args[1] : null);
		int k=Integer.parseInt(args[args.length-3]);
		int cbits=Integer.parseInt(args[args.length-2]);
		int k2=Integer.parseInt(args[args.length-1]);
		
		KCountArray2 counts=KmerCount3.countFastq(fname1, fname2, k, cbits);
		long[] counts2=countK2(fname1, fname2, k, counts, k2);
		
		t.stop();
		System.out.println("Finished counting; time = "+t+"\n");
		
		for(int i=0; i<counts2.length; i++){
			System.out.println(i+":\t"+counts2[i]);
		}
	}
	
	/** Builds a new short-key table before scanning input for long-window support.
	 * @param fname1 Primary FASTQ input
	 * @param fname2 Optional second mate FASTQ
	 * @param k Short window length in 1..19
	 * @param cbits Cell width supported by KCountArray2
	 * @param k2 Long window length, at least k
	 * @return New 256-bin histogram of clipped minimum subwindow counts */
	public static long[] countK2(String fname1, String fname2, int k, int cbits, int k2){
		KCountArray2 counts=KmerCount3.countFastq(fname1, fname2, k, cbits);
		return countK2(fname1, fname2, k, counts, k2);
	}
	
	/** Scans both mates against a borrowed forward-key table without changing it.
	 * Maintains a sliding frequency histogram to find the least-supported subwindow.
	 * Undefined or low-quality bases interrupt qualifying windows.
	 * Normal completion closes the owned reader, but currently ignores its error result.
	 * Exceptions during counting bypass that close; incomplete-input detection is not guaranteed.
	 * @param fname1 Primary FASTQ input
	 * @param fname2 Optional second mate FASTQ
	 * @param k Short window length in 1..19, matching the table's key encoding
	 * @param counts1 Non-null table populated with forward k-base keys
	 * @param k2 Long window length, at least k; smaller values invalidate the ring
	 * @return New histogram counting each eligible k2-window occurrence once */
	public static long[] countK2(String fname1, String fname2, int k, KCountArray2 counts1, int k2){
		//TODO: Validate k2>=k before reader startup; smaller values create an invalid ring.
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		FastqReadInputStream fris1=new FastqReadInputStream(fname1, false);
		FastqReadInputStream fris2=(fname2==null ? null : new FastqReadInputStream(fname2, false));
		ConcurrentGenericReadInputStream cris=new ConcurrentGenericReadInputStream(fris1, fris2, KmerCount3.maxReads);
		
		cris.start();
		System.err.println("Started cris");
		boolean paired=cris.paired();
		
		long kmer=0; //current kmer
		int len=0;  //distance since last contig start or ambiguous base
		
		final long[] upperBound=new long[BOUND_LEN]; //Lowest upper bound provable of kmer count
		final int[] ring=new int[k2-k+1];
		final int[] subcount=new int[BOUND_LEN];
		final int maxValue=subcount.length-1;
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(Read r : reads){
					
					len=0;
					kmer=0;
					Arrays.fill(subcount, 0);
					
					byte[] bases=r.bases;
					byte[] quals=r.quality;
					for(int i=0; i<bases.length; i++){
						byte b=bases[i];
						int x=AminoAcid.baseToNumber[b];
						
						int ringpos=i%ring.length;
						int old=ring[ringpos];
						int value=0;
						
						if(x<0 || quals[i]<KmerCount3.minQuality){
							len=0;
							kmer=0;
						}else{
							kmer=((kmer<<2)|x)&mask;
							len++;
							
							if(len>=k){
								value=counts1.read(kmer);
							}
						}
						value=min(value, maxValue);
						
						ring[ringpos]=value;
						subcount[value]++;
						
						if(i>=ring.length){
							subcount[old]--;
						}
						
						if(len>=k2){
							int sub=0;
							while(sub<subcount.length && subcount[sub]==0){sub++;}
							assert(sub<subcount.length);
							upperBound[sub]++;
						}
						
					}
					
					if(r.mate!=null){
						bases=r.mate.bases;
						quals=r.mate.quality;

						len=0;
						kmer=0;
						Arrays.fill(subcount, 0);
						for(int i=0; i<bases.length; i++){
							byte b=bases[i];
							int x=AminoAcid.baseToNumber[b];
							
							int ringpos=i%ring.length;
							int old=ring[ringpos];
							int value=0;
							
							if(x<0 || quals[i]<KmerCount3.minQuality){
								len=0;
								kmer=0;
							}else{
								kmer=((kmer<<2)|x)&mask;
								len++;
								
								if(len>=k){
									value=counts1.read(kmer);
								}
							}
							value=min(value, maxValue);
							
							ring[ringpos]=value;
							subcount[value]++;
							
							if(i>=ring.length){
								subcount[old]--;
							}
							
							if(len>=k2){
								int sub=0;
								while(sub<subcount.length && subcount[sub]==0){sub++;}
								assert(sub<subcount.length);
								upperBound[sub]++;
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
			//TODO: Propagate the close error; otherwise malformed input can yield a partial histogram.
			ReadWrite.closeStreams(cris);
			System.err.println("Closed stream");
		}
		
		return upperBound;
	}
	
	/** Returns the smaller value. */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger value. */
	public static final int max(int x, int y){return x>y ? x : y;}
	
	/** Histogram width; values at least 255 share the last bin. */
	public static final int BOUND_LEN=256;
	
}
