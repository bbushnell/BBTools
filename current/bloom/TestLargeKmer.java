package bloom;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentGenericReadInputStream;
import stream.FastqReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class TestLargeKmer {
	
	public static void main(String args[]){
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
	
	public static long[] countK2(String fname1, String fname2, int k, int cbits, int k2){
		KCountArray2 counts=KmerCount3.countFastq(fname1, fname2, k, cbits);
		return countK2(fname1, fname2, k, counts, k2);
	}
	
	public static long[] countK2(String fname1, String fname2, int k, KCountArray2 counts1, int k2){
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
				//System.err.println("reads.size()="+reads.size());
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
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			System.err.println("Finished reading");
			cris.returnList(ln);
			System.err.println("Returned list");
			ReadWrite.closeStreams(cris);
			System.err.println("Closed stream");
		}
		
		return upperBound;
	}
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	
	public static final int BOUND_LEN=256;
	
}
