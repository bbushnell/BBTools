package bloom;

import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastqReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 6, 2012
 *
 */
public class LargeKmerCount {
	
public static void main(String[] args){
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>4 || args[1].contains(".") ? args[1] : null);
		int indexbits=Integer.parseInt(args[args.length-3]);
		int cbits=Integer.parseInt(args[args.length-2]);
		int k=Integer.parseInt(args[args.length-1]);
		
		KCountArray2 count=countFastq(fname1, fname2, indexbits, cbits, k);
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		long[] freq=count.transformToFrequency();

//		System.out.println(count+"\n");
//		System.out.println(Arrays.toString(freq)+"\n");
		
		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ",(100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			int x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ",(100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Unique:     \t         \t"+sum2);
		System.out.println("CollisionsA:\t         \t"+collisionsA);
		System.out.println("CollisionsB:\t         \t"+collisionsB);
		
		double modifier=(collisionsB)/(double)(32*collisionsA+8*collisionsB);
		
		System.out.println("Estimate:   \t         \t"+(sum2+collisionsA+collisionsB-(long)(collisionsA*modifier)));
		System.out.println();
		System.out.println("Singleton:  \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:     \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		
	}
	
	public static KCountArray2 countFastq(String reads1, String reads2, int indexbits, int cbits, int k){
		assert(indexbits>=1 && indexbits<40);
		collisionsA=0;
		collisionsB=0;
		final long cells=1L<<indexbits;
		final int kbits=ROTATE_DIST*k;
		final int xorShift=kbits%64;
		final long[] rotMasks=makeRotMasks(xorShift);
		final int[] buffer=new int[k];
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", indexbits="+indexbits+", cells="+cells+", cbits="+cbits);}
		if(verbose){System.err.println("xorShift="+xorShift+", rotMasks[3]="+Long.toHexString(rotMasks[3]));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
		
		FastqReadInputStream fris1=new FastqReadInputStream(reads1, false);
		FastqReadInputStream fris2=(reads2==null ? null : new FastqReadInputStream(reads2, false));
		ConcurrentGenericReadInputStream cris=new ConcurrentGenericReadInputStream(fris1, fris2, maxReads);
		
		cris.start();
		System.err.println("Started cris");
		boolean paired=cris.paired();
		System.err.println("Paired: "+paired);
		
		long kmer=0; //current kmer
		int len=0;  //distance since last contig start or ambiguous base
		
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
							kmer=(Long.rotateLeft(kmer,ROTATE_DIST)^x);
							len++;
							if(len>=k){
								if(len>k){kmer=kmer^rotMasks[x2];}
								long hashcode=kmer&0x7fffffffffffffffL;
								long code1=hashcode%(cells-3);
								long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
								int value=count.increment2(code1, 1);
								long temp=count.read(code2);
								if(temp>0){
									if(value==0){collisionsA++;}
									else{collisionsB++;}
								}
							}
						}
					}
					
					
					if(r.mate!=null){
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
								kmer=(Long.rotateLeft(kmer,ROTATE_DIST)^x);
								len++;
								if(len>=k){
									if(len>k){kmer=kmer^rotMasks[x2];}
									long hashcode=kmer&0x7fffffffffffffffL;
									long code1=hashcode%(cells-3);
									long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
									int value=count.increment2(code1, 1);
									long temp=count.read(code2);
									if(temp>0){
										if(value==0){collisionsA++;}
										else{collisionsB++;}
									}
								}
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
			ReadWrite.closeStream(cris);
			System.err.println("Closed stream");
			System.err.println("Processed "+readsProcessed+" reads.");
		}
		
		return count;
	}
	
	public static final long[] makeRotMasks(int rotDist){
		long[] masks=new long[4];
		for(long i=0; i<4; i++){
			masks[(int)i]=Long.rotateLeft(i, rotDist);
		}
		return masks;
	}
	
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
	
	public static long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	public static long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	
	public static boolean verbose=true;
	public static byte minQuality=-5;
	public static long readsProcessed=0;
	public static long maxReads=1000000L;
	public static final int ROTATE_DIST=2;

	public static long collisionsA=0;
	public static long collisionsB=0;
	
}
