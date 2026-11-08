package bloom;

import java.util.ArrayList;
import java.util.BitSet;

import dna.AminoAcid;
import fileIO.FileFormat;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Short-kmer counter with reusable stream loops and trusted-mask recount paths.
 * Inherited static counters are cumulative until reset by the caller.
 * Ordinary stream loops count root fragments in readsProcessed; trusted-mask
 * loops do not advance it. keysCounted counts complete windows from both mates
 * and every requested strand pass.
 *
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public class KmerCount6 extends KmerCountAbstract{

	/** Counts kmers from command-line FASTQ inputs and prints a frequency summary. */
	public static void main(String[] args){

		Timer t=new Timer();

		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;

		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		KCountArray count=null;

		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}

		count=count(fname1, fname2, k, cbits, true, null);


		t.stop();
		System.out.println("Finished counting; time = "+t);

		printStatistics(count);

	}

	/** Prints a coarse count histogram and singleton/useful-key summary, reporting zero occupied ratios for an empty table. */
	public static void printStatistics(KCountArray count){
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
			long x=0;
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
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+Tools.format("%.3f    ",(sum2==0 ? 0 : keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+Tools.format("%.3f%%   ",(sum2==0 ? 0 : 100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+Tools.format("%.3f%%   ",(sum2==0 ? 0 : 100l*x/(double)sum2))+"\t"+x);
	}

	/**
	 * Opens, starts and closes FASTQ input while counting contiguous short kmers.
	 * A null table creates {@code KCountArray.makeNew(1L<<(2*k), cbits)} and therefore
	 * requires {@code k<20}; a supplied table receives increments directly and may have a
	 * hashed geometry for {@code 20<=k<32}.  rcomp destructively reverse-complements reads.
	 */
	public static KCountArray count(String reads1, String reads2, int k, int cbits, boolean rcomp, KCountArray count){
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));

		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits);
		}

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
		if(verbose){System.err.println("Paired: "+paired);}

		count(cris, k, rcomp, count);

		cris.close();
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}


		return count;
	}

	/**
	 * Consumes an already-started read stream and returns every borrowed list without
	 * closing the stream.  The destination table must be non-null; currently its
	 * contiguous path is always used and its legacy gap branch is assertion-disabled.
	 */
	public static void count(ConcurrentReadInputStream cris, int k, boolean rcomp, KCountArray count){
		assert(k<32 && k>=1 && (count!=null || k<20));

		assert(count!=null);

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);


		if(true /*count.gap==0*/){
			final int kbits=2*k;
			final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));


			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;

					addRead(r, count, k, mask, rcomp);
					if(r.mate!=null){
						addRead(r.mate, count, k, mask, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}else{
			final int k1=(k+1)/2;
			final int k2=k/2;
			final int kbits1=2*k1;
			final int kbits2=2*k2;
//			final int gap=count.gap;
			int gap=0; assert(false);
			final long mask1=~((-1L)<<(kbits1));
			final long mask2=~((-1L)<<(kbits2));
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;

					addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
					if(r.mate!=null){
						addReadSplit(r.mate, count, k1, k2, mask1, mask2, gap, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}

		if(verbose){System.err.println("Finished reading");}
		cris.returnList(ln);
		if(verbose){System.err.println("Returned list");}
	}



	/**
	 * Opens, starts and closes FASTQ input, optionally masks untrusted bases in place,
	 * and counts contiguous short kmers into a supplied or newly allocated table.
	 * The trusted table must support canonical queries at k. With conservative=true,
	 * retain bases covered by sampled windows reaching thresh; otherwise clear bases
	 * covered by sampled windows below thresh, allowing unexamined bases to remain.
	 * Cleared bases become N; existing qualities become zero. detectStepsize is the
	 * positive window-sampling step; the scanner modes have different origins.
	 * maxReads limits root fragments. This path does not increment readsProcessed.
	 */
	public static KCountArray count(final String reads1, final String reads2, final int k,
			final int cbits, final boolean rcomp, KCountArray count,
			final KCountArray trusted, final long maxReads, final int thresh,
			final int detectStepsize, final boolean conservative){

		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));

//		System.out.println("k="+k+", kbits="+kbits+", mask="+Long.toHexString(mask)+", thresh="+thresh);
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);

		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits);
		}

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
		if(verbose){System.err.println("Paired: "+paired);}

		count(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);

		cris.close();
		if(verbose){System.err.println("Closed stream");}

//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);

		return count;
	}



	/**
	 * Consumes but does not close an open stream while optionally masking untrusted
	 * bases. Conservative mode retains sampled high-count coverage; the other mode
	 * clears sampled low-count coverage. Qualities are zeroed only when present.
	 * Does not increment readsProcessed. Reads shorter than k are skipped unchanged,
	 * including metadata and orientation; their mate is still processed. Masks require
	 * a positive detectStepsize and canonical queries at k.
	 */
	public static void count(final ConcurrentReadInputStream cris, final int k,
			final boolean rcomp, final KCountArray count, final KCountArray trusted,
			final int thresh, final int detectStepsize, final boolean conservative){

		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//System.err.println("reads.size()="+reads.size());
			for(Read r : reads){

				Read r2=r.mate;
				if(r.length()>=k){
					if(trusted!=null){
						//BLM112: the length guard excludes null coverage and pointless strand mutation.
						BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) :
							ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
//						System.out.println("\n"+toString(bs, r.length()));
//						System.out.println(new String(r.bases));
						//FIXED 2026-07-05 [bloom/KmerCount6#001]: guard r.quality!=null (FASTA input has null
						//quals) - was an unguarded NPE. Matches the live twin ReadCounter.clearUntrustedBases.
						for(int i=bs.nextClearBit(0); i<r.length(); i=bs.nextClearBit(i+1)){
							r.bases[i]='N';
							if(r.quality!=null){r.quality[i]=0;}
						}
//						System.out.println(new String(r.bases));
//						System.out.println("used = "+Tools.format("%.3f%%",count.usedFraction()*100));
//						System.out.println("used = "+((KCountArray4)count).cellsUsed());
//						if(bs.length()<r.length()){r=null;}
					}
//					if(r!=null){addRead(r, count, k, mask, rcomp);}
					addRead(r, count, k, mask, rcomp);
				}
				if(r2!=null && r2.length()>=k){
					if(trusted!=null){
						BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
							ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
						for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
							r2.bases[i]='N';
							if(r2.quality!=null){r2.quality[i]=0;}//FIXED 2026-07-05 [bloom/KmerCount6#001]: FASTA null-quals guard.
						}
					}
					addRead(r2, count, k, mask, rcomp);
				}

			}
			//System.err.println("returning list");
			cris.returnList(ln);
			//System.err.println("fetching list");
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		if(verbose){System.err.println("Finished reading");}
		cris.returnList(ln);
		if(verbose){System.err.println("Returned list");}
	}

	/**
	 * Counts every complete contiguous valid kmer window.  Invalid bases reset the
	 * window; low-quality bases reset it only when a quality array exists, so absent
	 * qualities bypass {@link #minQuality}.  {@code mask} must cover exactly the low
	 * 2*k bits.  When requested, {@code rcomp} reverses {@code r} in place and runs
	 * a second strand pass rather than canonicalizing to one occurrence.
	 */
	public static void addRead(final Read r, final KCountArray count, final int k,
			final long mask, boolean rcomp){
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0 || (quals!=null && quals[i]<minQuality)){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				if(len>=k){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					count.increment(kmer);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addRead(r, count, k, mask, false);
		}
	}

	/**
	 * Counts split windows in a read using masks for the low 2*k1 and 2*k2 bits.
	 * Invalid bases or low qualities reset both arms; rcomp reverse-complements the read
	 * in place for the recursive second pass.  TODO: revisit odd-arm/gap geometry before
	 * accepting split-window behavior.
	 */
	public static void addReadSplit(final Read r, final KCountArray count, final int k1,
			final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] bases=r.bases;
		byte[] quals=r.quality;

		assert(kmer1>=kmer2);

//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;

		for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
			int x1=AminoAcid.baseToNumber[bases[i]];
			int x2=AminoAcid.baseToNumber[bases[j]];
			if(x1<0 || x2<0 || (quals!=null && (quals[i]<minQuality || quals[j]<minQuality))){
				len=0;
				kmer1=0;
				kmer2=0;
			}else{
				kmer1=((kmer1<<2)|x1)&mask1;
				kmer2=((kmer2<<2)|x2)&mask2;
				len++;
				if(len>=k1){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));

					long key=(kmer1<<shift)|kmer2;
//					System.err.println(Long.toHexString(key));
					count.increment(key);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addReadSplit(r, count, k1, k2, mask1, mask2, gap, false);
		}
	}

	/**
	 * Counts split windows over raw bases without qualities.  The masks must cover the
	 * low 2*k1 and 2*k2 bits; each counted key emits kmer1, kmer2 and merged key hex
	 * lines.  When requested, rcomp reverse-complements {@code bases} in place and
	 * recurses for the second pass without restoring it.  TODO: revisit odd-arm/gap
	 * geometry before accepting split-window behavior.
	 */
	public static void addReadSplit(final byte[] bases, final KCountArray count,
			final int k1, final int k2, final long mask1, final long mask2, final int gap,
			boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] quals=null;

		assert(kmer1>=kmer2);

//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;

		for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
			int x1=AminoAcid.baseToNumber[bases[i]];
			int x2=AminoAcid.baseToNumber[bases[j]];
			if(x1<0 || x2<0 || (quals!=null && (quals[i]<minQuality || quals[j]<minQuality))){
				len=0;
				kmer1=0;
				kmer2=0;
			}else{
				kmer1=((kmer1<<2)|x1)&mask1;
				kmer2=((kmer2<<2)|x2)&mask2;
				len++;
				if(len>=k1){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));

					long key=(kmer1<<shift)|kmer2;
					System.out.println(Long.toHexString(kmer1));
					System.out.println(Long.toHexString(kmer2));
					System.out.println(Long.toHexString(key));
					count.increment(key);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			Vector.reverseComplementInPlaceFast(bases);
			addReadSplit(bases, count, k1, k2, mask1, mask2, gap, false);
		}
	}

}
