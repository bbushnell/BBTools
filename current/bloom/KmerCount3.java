package bloom;

import java.util.ArrayList;

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
 * @date Jul 5, 2012
 *
 */
public class KmerCount3 extends KmerCountAbstract {
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>3 || args[1].contains(".") ? args[1] : null);
		int k=Integer.parseInt(args[args.length-2]);
		int cbits=Integer.parseInt(args[args.length-1]);
		
		KCountArray2 count=null;
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		count=countFastq(fname1, fname2, k, cbits);
		
		
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
	}
	
	public static KCountArray2 countFastq(String reads1, String reads2, int k, int cbits){
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final long cells=mask+1;
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
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
						if(x<0 || quals[i]<minQuality){
							len=0;
							kmer=0;
						}else{
							kmer=((kmer<<2)|x)&mask;
							len++;
							if(len>=k){
//								System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
								count.increment(kmer, 1);
//								System.out.println(" -> "+count.read(kmer));
//								System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//								array[(int)kmer]++;
//								System.out.println(" -> "+array[(int)kmer]+"\n");
//								assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
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
							if(x<0 || quals[i]<minQuality){
								len=0;
								kmer=0;
							}else{
								kmer=((kmer<<2)|x)&mask;
								len++;
								if(len>=k){
									count.increment(kmer, 1);
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
	
}
