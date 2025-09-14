package bloom;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import kmer.KmerTableSet;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SketchObject;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class ReadCounter_old extends KmerCountAbstract {
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
//			outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int matrixbits=-1;
		int hashes=1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("canonical")){
				CANONICAL=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		int kbits=Tools.min(2*k, 62);
		if(matrixbits<0){
			matrixbits=kbits;
		}
		matrixbits=Tools.min(kbits, matrixbits);
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		
		KCountArray counts=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
		ReadCounter_old rc=new ReadCounter_old(k, true, false, false, false);
		try {
			rc.count(fname1, fname2, counts);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		counts.shutdown();
		
//		verbose=true;
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		rc.printStatistics(counts);
		
	}
	
	/** Defaults for nucleotides. */
	public ReadCounter_old(final int k_){this(k_, true, false, false, false);}
	
	public ReadCounter_old(final int k_, final boolean rcomp_, 
			final boolean ecco_, final boolean merge_, final boolean amino_){
		k=k_;
		rcomp=rcomp_;
		ecco=ecco_;
		merge=merge_;
		amino=amino_;

		final int bitsPerChar=(amino ? AminoAcid.AMINO_SHIFT : 2);
		aminoShift=AminoAcid.AMINO_SHIFT;
		shift=bitsPerChar*k;
		shift2=shift-bitsPerChar;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		assert(!amino || k*bitsPerChar<64);
		assert(!amino || !rcomp);
		assert(k>0);
	}

	public void printStatistics(KCountArray counts){
		long[] freq=counts.transformToFrequency();

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
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+Tools.format("%.3f    ",(keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+Tools.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
	}
	
//	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, int cbits){
//		return makeKca(fname1, fname2, extraFiles, cbits, Tools.min(2*k, 35), 1, minQuality, 
//				maxReads, 1, 1, 1, 2);
//	}
	
//	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
//			int cbits, int matrixbits, int hashes, int minqual, long maxreads){
//		assert(matrixbits<63);
//		return makeKca(fname1, fname2, extraFiles, cbits, matrixbits, hashes, minqual, 
//				maxreads, 1, 1, 1, 2);
//	}
	
//	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
//			int cbits, int matrixbits, int hashes, int minqual,
//			long maxreads, int passes, int stepsize, int thresh1, int thresh2){
//		assert(matrixbits<63);
//		return makeKca(fname1, fname2, extraFiles,
//				cbits, 1L<<matrixbits, hashes, minqual,
//				maxreads, passes, stepsize, thresh1, thresh2, null, 0);
//	}
	
//	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
//			int cbits, long cells, int hashes, int minqual,
//			long maxreads, int passes, int stepsize, int thresh1, int thresh2){
//		return makeKca(fname1, fname2, extraFiles,
//				cbits, cells, hashes, minqual,
//				maxreads, passes, stepsize, thresh1, thresh2, null, 0);
//	}
	
	public KCountArray makeKca_als(ArrayList<String> fname1, ArrayList<String> fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int stepsize, int thresh1, int thresh2, 
			KCountArray prefilter, int prefilterLimit_){
		String a=null, b=null;
		ArrayList<String> list=new ArrayList<String>();
		if(fname1!=null){
			for(int i=0; i<fname1.size(); i++){
				if(i==0){a=fname1.get(i);}
				else{list.add(fname1.get(i));}
			}
		}
		if(fname2!=null){
			for(int i=0; i<fname2.size(); i++){
				if(i==0){b=fname2.get(i);}
				else{list.add(fname2.get(i));}
			}
		}
		if(extraFiles!=null){
			for(String s : extraFiles){
				list.add(s);
			}
		}
		return makeKca(a, b, list.isEmpty() ? null : list, cbits, cells, hashes, minqual, 
				maxreads, passes, stepsize, thresh1, thresh2,
				prefilter, prefilterLimit_);
	}
	
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes,
			KCountArray prefilter, int prefilterLimit){
		return makeKca(fname1, fname2, extraFiles, cbits, cells, hashes, 0, -1, 1, 1, 1, 1, prefilter, prefilterLimit);
	}
	
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int stepsize, int thresh1, int thresh2,
			KCountArray prefilter, int prefilterLimit_){
//		verbose=true;
		if(verbose){System.err.println("Making kca from ("+fname1+", "+fname2+")\nk="+k+", cells="+Tools.toKMG(cells)+", cbits="+cbits);}
		
		if(fname1==null && fname2==null && extraFiles==null){
			return KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);
		}
		
		boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		long oldmax=maxReads;
		byte oldq=minQuality;
		maxReads=maxreads;
		minQuality=(byte)minqual;
		//	System.out.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);
		
//		System.out.println("a");
		{//For processing input lists
			ArrayList<String> extra2=null;
			if(fname1!=null && fname1.contains(",")){
				String[] s=fname1.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname1=s[0];
			}
			if(fname2!=null && fname2.contains(",")){
				String[] s=fname2.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname2=s[0];
			}
			if(extra2!=null){
				if(extraFiles!=null){
					for(String s : extraFiles){
						extra2.add(s);
					}
				}
				extraFiles=extra2;
			}
		}
//		System.out.println("b");
		
		if(extraFiles!=null){
			for(String s : extraFiles){
				if(fileIO.FileFormat.hasFastaExtension(s)){
					assert(!FastaReadInputStream.SPLIT_READS);
				}
			}
		}
		
//		System.out.println("c");
		
		if(passes==1){
//			System.out.println("c1");
			if(fname1!=null){
				try {
					count(fname1, fname2, kca);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						count(s, null, kca);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();

		}else{
//			System.out.println("c2");
			assert(passes>1);
			KCountArray trusted=null;
			for(int i=1; i<passes; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
				//			if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}

				try {
					count(fname1, fname2, cbits, kca, trusted, maxreads, thresh1, step, conservative);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						try {
							count(s, null, cbits, kca, trusted, maxreads, thresh1, step, conservative);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
				kca.shutdown();
				
				System.out.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);

			}

			try {
				count(fname1, fname2, cbits, kca, trusted, maxreads, thresh2, stepsize, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						count(s, null, cbits, kca, trusted, maxreads, thresh2, stepsize, true);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();
		}
//		System.out.println("d");
		minQuality=oldq;
		maxReads=oldmax;
		FastaReadInputStream.SPLIT_READS=oldsplit;
		
		
		return kca;
	}
	
	public KCountArray count(String reads1, String reads2, KCountArray counts) throws Exception{
//		System.err.println("countFastq...  making a new cris");
		assert(counts!=null);
		
		{
			int pound=reads1.lastIndexOf('#');
			if(pound>=0 && reads2==null && !new File(reads1).exists()){
				String a=reads1.substring(0, pound);
				String b=reads1.substring(pound+1);
				reads1=a+1+b;
				reads2=a+2+b;
			}
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			if(ff2==null){ff1.preferShreds=true;}
//			if(ff2!=null){ //TODO - interleaved flag
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
//		countFastq(cris, count);
//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[Shared.threads()];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, counts);
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
//		System.out.println("~6");
		
		ReadWrite.closeStream(cris);
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return counts;
	}
	
	public KCountArray count(String reads1, String reads2, final int cbits, 
			KCountArray counts, final KCountArray trusted, final long maxReads, final int thresh, 
			final int detectStepsize, final boolean conservative)
			throws Exception{
		
		assert(counts!=null);
		
		{
			int pound=reads1.lastIndexOf('#');
			if(pound>=0 && reads2==null && !new File(reads1).exists()){
				String a=reads1.substring(0, pound);
				String b=reads1.substring(pound+1);
				reads1=a+1+b;
				reads2=a+2+b;
			}
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			if(ff2==null){ff1.preferShreds=true;}
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		

//		countFastq(cris, count, trusted, thresh, detectStepsize, conservative);

//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[Shared.threads()];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, counts, trusted, thresh, detectStepsize, conservative);
			cta[i].start();
		}
		
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		
//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		return counts;
	}
	
	private final int findOverlap(Read r1, Read r2, boolean ecc){
		if(vstrict){
			return BBMerge.findOverlapVStrict(r1, r2, ecc);
		}else{
			return BBMerge.findOverlapStrict(r1, r2, ecc);
		}
	}
	
	private class CountThread extends Thread{
		
		CountThread(final ConcurrentReadInputStream cris_, final KCountArray counts_){
			this(cris_, counts_, null, 2, 1, true);
		}
		
		CountThread(final ConcurrentReadInputStream cris_,
				final KCountArray counts_, final KCountArray trusted_, final int thresh_,
				final int detectStepsize_, final boolean conservative_){
			cris=cris_;
			counts=counts_;
			trusted=trusted_;
			thresh=thresh_;
			detectStepsize=detectStepsize_;
			conservative=conservative_;
		}
		
		@Override
		public void run(){
//			System.out.println("Running");
			if(trusted==null){
				count(cris);
			}else{
				count(cris, thresh, detectStepsize,  conservative);
			}
//			System.out.println("Finished: "+readsProcessedLocal);
			
			if(BUFFERED){dumpBuffer();}
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				increments+=incrementsLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}
		
		private void increment(final long key){
			if(BUFFERED){
				buffer.add(key);
				if(buffer.size>=BUFFERLEN){
					dumpBuffer();
				}
			}else{
				incrementByAmount(key, 1);
			}
		}
		
		
		private void dumpBuffer(){
			final int lim=buffer.size;
			if(lim<1){return;}
			final long[] array=buffer.array;
//			if(SORT_SERIAL){
//				buffer.sortSerial();
//			}else{
//				buffer.sort();
//			}
			buffer.sort();//Can be disabled via parallelsort flag but only affects arrays>10k
			long kmer=array[0]-1;//Ensures a nonmatch
			int streak=0;
			for(int i=0; i<lim; i++){
				long x=array[i];
				if(x==kmer){streak++;}
				else{
					if(streak>0){incrementByAmount(kmer, streak);}
					kmer=x;
					streak=1;
				}
			}
			assert(streak>0);
			incrementByAmount(kmer, streak);
			buffer.clear();
		}
		
		private void incrementByAmount(final long key, int amt){
			if(SKETCH_MODE){
				final long code=SketchObject.hash(key);
				if(code<minHashValue){return;}
				counts.increment(STORE_HASHED ? code : key, amt);
			}else{
				counts.increment(key, amt);
			}
			keysCountedLocal+=amt;
			incrementsLocal++;
		}
		
		private final void count(ConcurrentReadInputStream cris){
			assert(k>=1 && counts!=null);
//			System.out.println("Waiting for list");
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
//			System.out.println("Got list: "+(ln==null ? "null" : ln.id)+", "+(ln==null || ln.list==null ? "null" : ln.list.size()));

			long[] array=null;
			final Kmer kmer=new Kmer(k);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r1 : reads){

					Read r2=r1.mate;
					if((ecco || merge) && r2!=null){
						if(merge){
							final int insert=findOverlap(r1, r2, false);
							if(insert>0){
								r2.reverseComplement();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							findOverlap(r1, r2, true);
						}
					}
					readsProcessedLocal++;

					if(k<=maxShortKmerLength){
						array=addRead_Advanced(r1, array);
					}else{
						addReadBig(r1, kmer);
						addReadBig(r1.mate, kmer);
					}
					//						System.out.println(r);
					//						System.out.println("kmers hashed: "+keysCountedLocal);
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
		
		


		private final void count(final ConcurrentReadInputStream cris, final int thresh,
				final int detectStepsize, final boolean conservative){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			long[] array=null;
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r1 : reads){
					
					Read r2=r1.mate;
					if((ecco || merge) && r2!=null){
						if(merge){
							final int insert=findOverlap(r1, r2, false);
							if(insert>0){
								r2.reverseComplement();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							findOverlap(r1, r2, true);
						}
					}
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r1, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r1, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.length()));
//							System.out.println(new String(r.bases));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r1.length(); i=bs.nextClearBit(i+1)){
									r1.bases[i]='N';
									if(r1.quality!=null){r1.quality[i]=0;}
								}
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+Tools.format("%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.length()){r=null;}
						}
//						if(r!=null){addRead(r, count, rcomp);}
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
									r2.bases[i]='N';
									if(r2.quality!=null){r2.quality[i]=0;}
								}
							}
						}
					}
					array=addRead_Advanced(r1, array);

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
		 * Hash a read's kmers into the KCountArray.
		 * Advanced mode processes paired reads together and sorts kmers to eliminate spurious duplicates.
		 * @param r1
		 * @param countMatrix
		 * @param kbloom
		 * @param mask
		 * @param rcomp
		 */
		private final long[] addRead_Advanced(Read r1, long[] array){
			if(PREJOIN && r1.mate!=null && r1.insert()>0){
				r1.mate.reverseComplement();
				r1=r1.joinRead();
			}
			Read r2=r1.mate;
			final int len1=Tools.max(0, r1.length()-k+1);
			final int len2=(r2==null || r2.bases==null) ? 0 : Tools.max(0, r2.length()-k+1);
			final int len=len1+len2;
			if(len<1){return array;}
			
			if(!KEEP_DUPLICATE_KMERS){
				if(array==null || array.length!=len){array=new long[len];}//TODO: Should use LongList
				Arrays.fill(array, -1);
				fillKmerArray(r1, array, 0, len1);
				if(r2!=null){fillKmerArray(r2, array, len1, len);}
				
				Arrays.sort(array);
				long prev=-1;
				for(int i=0; i<array.length; i++){
					long kmer=array[i];
					if(kmer!=prev){
						increment(kmer);
						prev=kmer;
					}
				}
			}else{
				if(len1>0){addRead(r1);}
				if(len2>0){addRead(r2);}
			}
			return array;
		}
		
		private final void addReadBig(Read r, Kmer kmer){
			if(r==null || r.bases==null){return;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			int len=0;
			
			if(bases==null || bases.length<k){return;}
			kmer.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];

				//Update kmers
				kmer.addRight(b);
				
				if(minProb>0 && quals!=null){//Update probability
					prob=prob*KmerTableSet.PROB_CORRECT[quals[i]];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*KmerTableSet.PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					prob=1;
				}else{len++;}
				
				assert(len==kmer.len);
				
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
				if(len>=k && prob>=minProb){
//					System.err.println("Incrementing xor()="+kmer.xor());
					increment(kmer.xor());
//					counts.incrementAndReturnUnincremented(kmer.xor(), 1);
//					keysCountedLocal++;
				}
			}
		}
		
		private final void fillKmerArray(Read r, final long[] array, final int start, final int stop){
			if(amino){
				fillKmerArrayAmino(r, array, start, stop);
				return;
			}
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, array, start, stop);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k){return;}
			
			final int passes=(rcomp ? 2 : 1);
			for(int pass=0; pass<passes; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					assert(b>=0) : r.id+", "+bases.length+", "+i+", "+b+"\n"+(bases.length<20000 ? r.toFasta().toString() : "");
					int x=AminoAcid.baseToNumber[b];
					
//					int x=AminoAcid.baseToNumber[b<0 ? 'N' : b];

					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && prob>=minProb){
							array[idx]=Tools.max(array[idx], kmer);
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void fillKmerArrayAmino(Read r, final long[] array, final int start, final int stop){
			assert(false) : "TODO"; //TODO
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, array, start, stop);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k){return;}
			
			for(int pass=0; pass<2; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					assert(b>=0) : r.id+", "+bases.length+", "+i+", "+b+"\n"+(bases.length<20000 ? r.toFasta().toString() : "");
					int x=AminoAcid.baseToNumber[b];
					
//					int x=AminoAcid.baseToNumber[b<0 ? 'N' : b];

					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && prob>=minProb){
							array[idx]=Tools.max(array[idx], kmer);
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void addRead(Read r){
			if(amino){
				addReadAmino(r);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);

			final byte[] bases=r.bases;
			final byte[] quals=r.quality;

			if(bases==null || bases.length<k){return;}
			
			long kmer=0;
			long rkmer=0;
			int len=0;
			float prob=1;

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				final byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}

				if(x<0 || q<minQuality){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{
					len++;
					if(len>=k && prob>=minProb){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						increment(key);
					}
				}
			}
		}
		
		private final void addReadAmino(Read r){
			assert(!PREJOIN || r.mate==null);

			final byte[] bases=r.bases;

			if(bases==null || bases.length<k){return;}
			
			long kmer=0;
			int len=0;

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.acidToNumber[b];
				kmer=((kmer<<aminoShift)|x)&mask;

				if(x<0){
					len=0;
					kmer=0;
				}else{
					len++;
					if(len>=k){
						increment(kmer);
					}
				}
			}
		}
		
		private final void fillKmerArrayLong(Read r, final long[] array, final int start, final int stop){
			assert(k>maxShortKmerLength) : k;
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			Kmer kmer=new Kmer(k);
			
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			kmer.clear();
			
			for(int i=0, idx=start-k+1; i<bases.length; i++, idx++){
				byte b=bases[i];
				kmer.addRight(b);
				
				byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(kmer.len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}
				
				if(!AminoAcid.isFullyDefined(b) || q<minQuality){
					kmer.clear();
					prob=1;
				}
				if(kmer.len>=k && prob>=minProb){
					array[idx]=kmer.xor();
				}
			}
		}
		
		private final ConcurrentReadInputStream cris;
		
		private final KCountArray counts;
		private final KCountArray trusted;
		private final int thresh;
		private final int detectStepsize;
		private final boolean conservative;
		private long keysCountedLocal=0;
		private long incrementsLocal=0;
		private long readsProcessedLocal=0;
		private final long minHashValue=SketchObject.minHashValue;
		private final LongList buffer=(BUFFERED ? new LongList(BUFFERLEN) : null);
	}
	
	public static boolean vstrict=false;
	
	private final int k;
	private final int aminoShift;
	private final int shift;
	private final int shift2;
	private final long mask;
	private final boolean rcomp;
	private final boolean ecco;
	private final boolean merge;
	private final boolean amino;
	
}
