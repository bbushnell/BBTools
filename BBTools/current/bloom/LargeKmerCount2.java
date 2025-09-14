package bloom;

import java.util.ArrayList;
import java.util.Random;

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
 * @author Brian Bushnell
 * @date Jul 6, 2012
 *
 */
public class LargeKmerCount2 {
	
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
		System.out.println("Singleton:     \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		
	}
	
	public static KCountArray2 countFastq(String reads1, String reads2, int indexbits, int cbits, int k){
		assert(indexbits>=1 && indexbits<40);
		final long cells=1L<<indexbits;
		final int kbits=ROTATE_DIST*k;
		final int xorShift=kbits%64;
		final long[] rotMasks=makeRotMasks(xorShift);
		final int[] buffer=new int[k];
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", indexbits="+indexbits+", cells="+cells+", cbits="+cbits);}
		if(verbose){System.err.println("xorShift="+xorShift+", rotMasks[3]="+Long.toHexString(rotMasks[3]));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
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
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
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
								keysCounted++;
								if(len>k){kmer=kmer^rotMasks[x2];}
								long hashcode=kmer&0x7fffffffffffffffL;
//								hashcode=randy.nextLong()&~((-1L)<<(2*k));
								long code1=hashcode%(cells-3);
//								long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
								int value=count.increment2(code1, 1);
								
								double probCollision=load*invCells;
//								expectedCollisions+=probCollision;
								expectedCollisions+=probCollision*(1-(load+min(expectedCollisions, actualCollisions))*invKmerSpace);
								if(value==0){load++;}
								else{
									actualCollisions++;
									double probNewKey=(load*invCells)*expectedCollisions/(min(expectedCollisions, actualCollisions));
									double estKeys=load+probNewKeyCollisions;
									double probOldKey=estKeys*invKmerSpace;
									probNewKeyCollisions+=probNewKey*(1-probOldKey);
									
//									double estKmers=load+min(actualCollisions, expectedCollisions);
//									double probOldKmer=estKmers*invKmerSpace;
//									probNewKeyCollisions+=(prob*(1-prob2));
								}
								
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
									keysCounted++;
									if(len>k){kmer=kmer^rotMasks[x2];}
									long hashcode=kmer&0x7fffffffffffffffL;
//									hashcode=randy.nextLong()&~((-1L)<<(2*k));
									long code1=hashcode%(cells-3);
//									long code2=((~hashcode)&0x7fffffffffffffffL)%(cells-5);
									int value=count.increment2(code1, 1);
									
									double probCollision=load*invCells;
//									expectedCollisions+=probCollision;
									expectedCollisions+=probCollision*(1-(load+min(expectedCollisions, actualCollisions))*invKmerSpace);
									if(value==0){load++;}
									else{
										actualCollisions++;
										double probNewKey=(load*invCells)*expectedCollisions/(min(expectedCollisions, actualCollisions));
										double estKeys=load+probNewKeyCollisions;
										double probOldKey=estKeys*invKmerSpace;
										probNewKeyCollisions+=probNewKey*(1-probOldKey);
										
//										double estKmers=load+min(actualCollisions, expectedCollisions);
//										double probOldKmer=estKmers*invKmerSpace;
//										probNewKeyCollisions+=(prob*(1-prob2));
									}
									
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
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	public static final double min(double x, double y){return x<y ? x : y;}
	public static final double max(double x, double y){return x>y ? x : y;}
	
	public static boolean verbose=true;
	public static byte minQuality=-5;
	public static long readsProcessed=0;
	public static long maxReads=10000000L;
	public static final int ROTATE_DIST=2;
	
	/** Non-empty cells in hash table */
	public static long load;
	/** Number of expected collisions */
	public static double expectedCollisions;
	/** Number of actual collisions (possibly by same value) */
	public static long actualCollisions;
	/** Number of probable collisions caused by new keys */
	public static double probNewKeyCollisions;
	/** Inverse of hash table size */
	public static double invCells;
	/** Inverse of number of potential kmers */
	public static double invKmerSpace;
	/** Inverse of number of potential kmers */
	public static long keysCounted;
	
	public static final Random randy=new Random(1);
	
}
