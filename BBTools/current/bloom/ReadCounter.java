package bloom;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.BitSet;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.Dedupe;
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
public class ReadCounter extends KmerCountAbstract {
	
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
		ReadCounter rc=new ReadCounter(k, true, false, false, false);
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
	public ReadCounter(final int k_){this(k_, true, false, false, false);}
	
	public ReadCounter(final int k_, final boolean rcomp_, 
			final boolean ecco_, final boolean merge_, final boolean amino_){
		k=k_;
		k2=k-1;
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
	
	public KCountArray makeKca_als(ArrayList<String> fname1, ArrayList<String> fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int thresh1, int thresh2, 
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
				maxreads, passes, thresh1, thresh2,
				prefilter, prefilterLimit_);
	}
	
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes,
			KCountArray prefilter, int prefilterLimit){
		return makeKca(fname1, fname2, extraFiles, cbits, cells, hashes, 0, -1, 1, 1, 1, prefilter, prefilterLimit);
	}
	
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int thresh1, int thresh2,
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
//				int step=(stepsize==1 ? 1 : stepsize+i%2);
//				//			if(!conservative){step=(step+3)/4;}
//				if(!conservative){step=Tools.min(3, (step+3)/4);}

				try {
					count(fname1, fname2, cbits, kca, trusted, maxreads, thresh1, conservative);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						try {
							count(s, null, cbits, kca, trusted, maxreads, thresh1, conservative);
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
				count(fname1, fname2, cbits, kca, trusted, maxreads, thresh2, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						count(s, null, cbits, kca, trusted, maxreads, thresh2, true);
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
		CountThread[] cta=new CountThread[Tools.min(Shared.threads(), MAX_COUNT_THREADS)];
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
			KCountArray counts, final KCountArray trusted, final long maxReads,
			final int thresh, final boolean conservative) throws Exception{
		
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

//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[Shared.threads()];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, counts, trusted, thresh, conservative);
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
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	public static final long toKmer(final byte[] bases, final int start, final int klen){
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	private class CountThread extends Thread{
		
		CountThread(final ConcurrentReadInputStream cris_, final KCountArray counts_){
			this(cris_, counts_, null, 2, true);
		}
		
		CountThread(final ConcurrentReadInputStream cris_,
				final KCountArray counts_, final KCountArray trusted_, final int thresh_,
				final boolean conservative_){
			cris=cris_;
			counts=counts_;
			trusted=trusted_;
			thresh=thresh_;
			conservative=conservative_;
		}
		
		@Override
		public void run(){
//			System.out.println("Running");
			count(cris);
//			System.out.println("Finished: "+readsProcessedLocal);
			
			if(BUFFERED){dumpBufferT();}
			
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
					dumpBufferT();
				}
			}else{
				incrementByAmountT(key, 1);
			}
		}
		
		private void incrementOneKmer(Read r, int knum) {
			if(r==null || r.length()<k) {return;}
			final long kmer=toKmer(r.bases, (int)((r.numericID+knum*k)%(r.length()-k2)), k);
			if(kmer>=0){increment(kmer);}
		}
		
		private void dumpBufferT(){
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
					if(streak>0){incrementByAmountT(kmer, streak);}
					kmer=x;
					streak=1;
				}
			}
			assert(streak>0);
			incrementByAmountT(kmer, streak);
			buffer.clear();
		}
		
		private void incrementByAmountT(final long key, int amt){
			if(SKETCH_MODE){
				final long code=SketchObject.hash(key);
				if(code<sketchMinHashValue){return;}
				counts.increment(STORE_HASHED ? code : key, amt);
			}else{
				counts.increment(key, amt);
			}
			keysCountedLocal+=amt;
			incrementsLocal++;
		}
		
		private final void count(final ConcurrentReadInputStream cris){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			final LongList buffer=new LongList(300);
			final Kmer kmer=(k>maxShortKmerLength ? new Kmer(k) : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r1 : reads){
					readsProcessedLocal+=r1.pairCount();
					
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
					
					if(trusted!=null){
						clearUntrustedBases(r1);
						clearUntrustedBases(r2);
					}

					if(amino) {
						addReadAmino(r1);
					}else if(k<=maxShortKmerLength){
						addRead_Advanced(r1, buffer);
					}else{
						addReadBig(r1, kmer);
						addReadBig(r1.mate, kmer);
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
		
		private int clearUntrustedBases(Read r) {
			if(r==null || trusted==null) {return 0;}
			BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) :
				ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
			if(bs==null){return 0;}
			int ns=0;
			for(int i=bs.nextClearBit(0); i<r.length(); i=bs.nextClearBit(i+1)){
				r.bases[i]='N';
				if(r.quality!=null){r.quality[i]=0;}
				ns++;
			}
			return ns;
		}
		
		/**
		 * Hash a read's kmers into the KCountArray.
		 * Advanced mode processes paired reads together,
		 * and sorts kmers to eliminate spurious duplicates.
		 * @param r1
		 * @param buffer
		 */
		private final void addRead_Advanced(Read r1, final LongList buffer){
			assert(k<=maxShortKmerLength);
			assert(!amino);
			
			buffer.clear();
			if(PREJOIN && r1.mate!=null && r1.insert()>0){
				r1.mate.reverseComplement();
				r1=r1.joinRead();
				assert(r1.mate==null);
			}
			Read r2=r1.mate;
			final int len1=Tools.max(0, r1.length()-k+1);
			final int len2=(r2==null || r2.bases==null) ? 0 : Tools.max(0, r2.length()-k+1);
			
			if(KMERS_PER_READ>0) {
				final int lim1=Tools.min(KMERS_PER_READ, r1.length()-k2);
				final int lim2=Tools.min(KMERS_PER_READ, r1.mateLength()-k2);
				for(int i=0; i<lim1; i++) {incrementOneKmer(r1, i);}
				for(int i=0; i<lim2; i++) {incrementOneKmer(r2, i);}
			}else if(!KEEP_DUPLICATE_KMERS){
				fillKmerArray(r1, buffer);
				if(r2!=null){fillKmerArray(r2, buffer);}
				if(buffer.size()<1) {return;}
				
				buffer.sort();
				long prev=-1;
				final long[] array=buffer.array;
				final int idmod=(int)(r1.numericID&IDMASK);
				for(int i=0, lim=buffer.size(); i<lim; i++){
					long kmer=array[i];
					if(kmer!=prev && (i&IDMASK)==idmod){
						increment(kmer);
						prev=kmer;
					}
				}
			}else{
				if(len1>0){addRead(r1);}
				if(len2>0){addRead(r2);}
			}
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
		
		private final void fillKmerArray(Read r, final LongList list){
			if(amino){
				fillKmerArrayAmino(r, list);
				return;
			}
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, list);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(list!=null);
			
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
						list.add(key);
					}
				}
			}
		}
		
		private final void fillKmerArrayAmino(Read r, final LongList array){
			assert(false) : "TODO"; //TODO
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
			final int idmod=(int)(r.numericID&IDMASK);
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
					if(len>=k && prob>=minProb && (i&IDMASK)==idmod){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						increment(key);
					}
				}
			}
		}
		
		private final void addReadAmino(Read r){
			if(r==null) {return;}
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
		
		private final void fillKmerArrayLong(Read r, final LongList list){
			assert(k>maxShortKmerLength) : k;
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(list!=null);
			Kmer kmer=new Kmer(k);
			
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			kmer.clear();
			
			for(int i=0; i<bases.length; i++){
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
					list.add(kmer.xor());
				}
			}
		}
		
		private long keysCountedLocal=0;
		private long incrementsLocal=0;
		private long readsProcessedLocal=0;
		
		private final ConcurrentReadInputStream cris;
		
		private final KCountArray counts;
		private final KCountArray trusted;
		private final int thresh;
		private final boolean conservative;
		private final long sketchMinHashValue=SketchObject.minHashValue;
		
		private final LongList buffer=(BUFFERED ? new LongList(BUFFERLEN) : null);
	}
	
	public int detectStepsize=1;//Jump this many bases when detecting errors
	
	private final int k;
	private final int k2;
	private final int aminoShift;
	private final int shift;
	private final int shift2;
	private final long mask;
	private final boolean rcomp;
	private final boolean ecco;
	private final boolean merge;
	private final boolean amino;
	
	public static boolean vstrict=false;
	
}
