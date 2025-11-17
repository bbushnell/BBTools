package clump;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import bloom.KCountArray;
import dna.AminoAcid;
import jgi.BBMerge;
import shared.Shared;
import shared.Tools;
import sketch.SketchTool;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 4, 2015
 *
 */
public class KmerComparator implements Comparator<Read> {
	
	/**
	 * Creates a KmerComparator with default seed, border, and hash parameters.
	 * @param k_ K-mer length for comparison
	 * @param addName_ Whether to append key information to read names
	 * @param rcomp_ Whether to reverse-complement reads when using minus strand k-mers
	 */
	public KmerComparator(int k_, boolean addName_, boolean rcomp_){
		this(k_, defaultSeed, defaultBorder, defaultHashes, addName_, rcomp_);
	}
	
	/**
	 * Creates a KmerComparator with fully customizable parameters.
	 *
	 * @param k_ K-mer length for comparison (must be >0 and <32)
	 * @param seed_ Random seed for hash function initialization
	 * @param border_ Number of bases to ignore at read ends (reduces end-effect bias)
	 * @param hashes_ Number of hash iterations (0-8, more iterations = better distribution)
	 * @param addName_ Whether to append key information to read names
	 * @param rcomp_ Whether to reverse-complement reads when using minus strand k-mers
	 */
	public KmerComparator(int k_, long seed_, int border_, int hashes_, boolean addName_, boolean rcomp_){
		k=k_;
		assert(k>0 && k<32);
		
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		seed=seed_;
		border=Tools.max(0, border_);
		hashes=Tools.mid(0, hashes_, 8);
		codes=SketchTool.makeCodes(8, 256, seed_, true);
		if(verbose){
			System.err.println("Made a comparator with k="+k+", seed="+seed+", border="+border+", hashes="+hashes);
		}
		addName=addName_;
		rcompReads=rcomp_;
	}
	
	/**
	 * Hashes reads using multiple threads for performance.
	 * Creates one HashThread per available CPU core to process reads in parallel.
	 *
	 * @param list List of reads to hash
	 * @param table K-mer count table for frequency filtering
	 * @param minCount Minimum k-mer count threshold for acceptance
	 */
	public void hashThreaded(ArrayList<Read> list, KCountArray table, int minCount){
		int threads=Shared.threads();
		ArrayList<HashThread> alt=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alt.add(new HashThread(i, threads, list, table, minCount));}
		for(HashThread ht : alt){ht.start();}
		
		/* Wait for threads to die */
		for(HashThread ht : alt){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Hashes all reads in a list sequentially.
	 *
	 * @param list List of reads to hash
	 * @param table K-mer count table for frequency filtering
	 * @param minCount Minimum k-mer count threshold for acceptance
	 * @param setObject Whether to store the ReadKey in the read's obj field
	 */
	public void hash(ArrayList<Read> list, KCountArray table, int minCount, boolean setObject) {
		for(Read r : list){hash(r, table, minCount, setObject);}
	}
	
	/**
	 * Hashes a single read and optionally its mate.
	 * If containment mode is enabled, both reads in a pair are processed.
	 *
	 * @param r1 The read to hash
	 * @param table K-mer count table for frequency filtering
	 * @param minCount Minimum k-mer count threshold for acceptance
	 * @param setObject Whether to store the ReadKey in the read's obj field
	 * @return The hash value of the selected k-mer
	 */
	public long hash(Read r1, KCountArray table, int minCount, boolean setObject){
		long x=hash_inner(r1, table, minCount, setObject);
		if(Clump.containment && r1.mate!=null){hash_inner(r1.mate, table, minCount, setObject);}
		return x;
	}
	
	/**
	 * Internal method that performs the actual k-mer hashing for a single read.
	 * Creates or reuses a ReadKey and calls fillMax to find the optimal k-mer.
	 *
	 * @param r1 The read to hash
	 * @param table K-mer count table for frequency filtering
	 * @param minCount Minimum k-mer count threshold for acceptance
	 * @param setObject Whether to store the ReadKey in the read's obj field
	 * @return The hash value of the selected k-mer
	 */
	private long hash_inner(Read r1, KCountArray table, int minCount, boolean setObject){
		ReadKey key;
		if(setObject){
			if(r1.obj==null){
				key=ReadKey.makeKey(r1, true);
			}else{
				key=(ReadKey)r1.obj;
				key.clear();
			}
		}else{
			key=getLocalKey();
		}
		assert(key.kmer==key.position && key.position==0);
		return fillMax(r1, key, table, minCount);
	}
	
//	public void fuse(Read r1){
//		Read r2=r1.mate;
//		if(r2==null){return;}
//		r1.mate=null;
//		final int len1=r1.length(), len2=r2.length();
//		int len=len1+len2+1;
//		byte[] bases=new byte[len];
//		for(int i=0; i<len1; i++){bases[i]=r1.bases[i];}
//		bases[len1]='N';
//		for(int i=0, j=len1+1; i<len2; i++){bases[j]=r2.bases[i];}
//	}
	
	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(Read a, Read b) {
		final ReadKey keyA, keyB;
		if(a.obj==null){
			keyA=ReadKey.makeKey(a, true);
			fillMax(a, keyA, null, 0);
		}else{keyA=(ReadKey)a.obj;}

		if(b.obj==null){
			keyB=ReadKey.makeKey(b, true);
			fillMax(b, keyB, null, 0);
		}else{keyB=(ReadKey)b.obj;}
		
		int x=keyA.compareTo(keyB);
		if(x==0 && compareSequence){
			x=KmerComparator2.compareSequence(a, b, 0);
		}
		if(x==0){
			float ea=a.expectedErrorsIncludingMate(true);
			float eb=b.expectedErrorsIncludingMate(true);
			if(ea!=eb){return ea>eb ? -1 : 1;}
		}
		return x==0 ? a.id.compareTo(b.id) : x;
	}
	
	/** Finds the global maximum */
	public long fillMax(Read r, ReadKey key, KCountArray table, int minCount){
		if(mergeFirst && r.pairnum()==0 && r.mate!=null){//This is probably unsafe in multithreaded mode unless the same thread handles both reads.
			int x=BBMerge.findOverlapStrict(r, r.mate, false);
			if(x>0){
				if(r.swapped()==r.mate.swapped()){r.mate.reverseComplementFast();}
				Read merged=r.joinRead(x);
				if(r.swapped()==r.mate.swapped()){r.mate.reverseComplementFast();}
				fillMax(merged, key, table, minCount);
				if(key.flipped){
					r.reverseComplementFast();
					r.setSwapped(true);
				}
				return key.kmer;
			}
		}
		if(r.length()<k){return fillShort(r, key);}
		assert(minCount>0 || table==null) : minCount;
		assert(table==null || minCount<=table.maxValue) : minCount;
		
		key.set(0, k-1, false); //TODO: Why is this k-1 rather than 0?
		final byte[] bases=r.bases, quals=r.quality;
		long kmer=0;
		long rkmer=0;
		int len=0;
		float prob=1;
		
		if(bases==null || bases.length<k){return -1;}
		
		long topCode=Long.MIN_VALUE;
		int topCount=-2;
		float topProb=0;
		final int localBorder=(bases.length>k+4*border ? border : 0);
		
		final int max=bases.length-localBorder;
		for(int i=localBorder; i<max; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(minProb>0 && quals!=null){//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>k){
					byte oldq=quals[i-k];
					prob=prob*PROB_CORRECT_INVERSE[oldq];
				}
			}
			
			if(x<0){
				len=0;
				prob=1;
			}else{len++;}
			
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				final long code=hash(kmax);
				boolean accept=false;
				int count=0;
				if(table!=null){
					assert(minCount>=1);
					if(code>topCode){
						count=table.read(kmax);
						accept=(count>=minCount || topCount<minCount);
					}else if(topCount<minCount){
						count=table.read(kmax);
						accept=count>=minCount;
					}
				}else{
					if(code>topCode){
						accept=(prob>=minProb || topProb<minProb);
					}else{
						accept=(prob>=minProb && topProb<minProb);
					}
				}
				
				if(accept){
					topCode=code;
					topCount=count;
					topProb=minProb;
					key.set(kmax, i, (kmax!=kmer));
				}
			}
		}
		if(topCode==Long.MIN_VALUE){
			return fillShort(r, key);
		}
//		if(bases.length<k){
//			final long kmax=Tools.max(kmer, rkmer);
//			key.set(kmax, bases.length-1, (kmax!=kmer));
//		}
		
//		assert(minCount<2) : minCount+", "+topCode+", "+topCount;
//		assert(minCount>0) : minCount+", "+topCode+", "+topCount;
		
//		if(topCode<0 && minCount>1){//Not needed
//			return fillMax(r, kmers, null, 0);
//		}
		
//		r.id+=" "+key.position+","+rcomp+","+(bases.length-key.position+k-2);
		if(key.kmerMinusStrand && rcompReads){
			key.flip(r, k);
		}
//		if(shortName){//This actually takes a lot of time!
//			r.id=r.numericID+" 1:"+(rcomp ? "t" : "f");
//			if(r.mate!=null){
//				r.mate.id=r.numericID+" 2:f";
//			}
		//		}else
		if(addName){r.id+=" "+key;}

		assert(key.kmer>=0 && key.position>=0) : key+"\n"+r;
		return key.kmer;
	}
	

	
	/* For teaching */
	/** Finds the global maximum, forward only */
	public long example1(Read r){
		//Handle degenerate reads shorter than kmer length
		if(r.length()<k){return fillShort(r, null);}
		
		final byte[] bases=r.bases;
		long kmer=0;//Stores the current kmer being built and analyzed
		int len=0;//Length bases since last N encountered
		
		long topCode=Long.MIN_VALUE;//Maximal hash code encountered
		long topKmer=Long.MIN_VALUE;//Kmer for maximal hash code
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];//The base, A, C, G, T, or N
			long x=AminoAcid.baseToNumber[b];//2-bit value, 0, 1, 2, 3, or -1
			kmer=((kmer<<2)|x)&mask;//Shift the new 2 bits into the kmer
			
			if(x<0){//If an N is encountered, reset the length
				len=0;
			}else{len++;}
			
			if(len>=k){//If the kmer is valid and full length, hash it
				final long code=hash(kmer);//Some arbitrary hash function
				
				if(code>topCode){//If this is a new maximum, retain it
					topKmer=kmer;
					topCode=code;
				}
			}
		}
		
		//If for whatever reason it failed (e.g., Ns), use the short method
		if(topCode==Long.MIN_VALUE){return fillShort(r, null);}
		
		//Return the pivot kmer.
		//It's also useful to store (and return) the position and orientation.
		return topKmer;
	}
	
	/* For teaching */
	/** Finds the global maximum, forward and reverse */
	public long example2(Read r){
		if(r.length()<k){return fillShort(r, null);}
		
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		if(bases==null || bases.length<k){return -1;}
		
		long topCode=Long.MIN_VALUE;
		long topKmer=Long.MIN_VALUE;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(x<0){
				len=0;
			}else{len++;}
			
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				final long code=hash(kmax);
				
				if(code>topCode){
					topKmer=kmax;
					topCode=code;
				}
			}
		}
		if(topCode==Long.MIN_VALUE){
			return fillShort(r, null);
		}
		
		return topKmer;
	}

	/* For teaching */
	/** Finds the global maximum, forward and reverse, with qualities and border */
	public long example3(Read r){
		if(r.length()<k){return fillShort(r, null);}
		
		final byte[] bases=r.bases, quals=r.quality;
		long kmer=0;
		long rkmer=0;
		int len=0;
		float prob=1;
		
		if(bases==null || bases.length<k){return -1;}
		
		long topCode=Long.MIN_VALUE;
		long topKmer=Long.MIN_VALUE;
		float topProb=0;
		final int localBorder=(bases.length>k+4*border ? border : 0);
		
		final int max=bases.length-localBorder;
		for(int i=localBorder; i<max; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(minProb>0 && quals!=null){//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>k){
					byte oldq=quals[i-k];
					prob=prob*PROB_CORRECT_INVERSE[oldq];
				}
			}
			
			if(x<0){
				len=0;
				prob=1;
			}else{len++;}
			
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				final long code=hash(kmax);
				boolean accept=false;
				
				if(code>topCode){
					accept=(prob>=minProb || topProb<minProb);
				}else{
					accept=(prob>=minProb && topProb<minProb);
				}
				
				if(accept){
					topKmer=kmax;
					topCode=code;
					topProb=minProb;
				}
			}
		}
		if(topCode==Long.MIN_VALUE){
			return fillShort(r, null);
		}
		
		return topKmer;
	}
	
	/** Generates a key when the read is shorter than k */
	public long fillShort(Read r, ReadKey key){
		final byte[] bases=r.bases;
		final int max=Tools.min(bases.length, k);
		key.set(0, max-1, false);
		long kmer=0;
		long rkmer=0;
		
		for(int i=0; i<max; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber0[b];
			long x2=AminoAcid.baseToComplementNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
		}

		final long kmax=Tools.max(kmer, rkmer);
		key.set(kmax, max-1, (kmax!=kmer));
		
		if(key.kmerMinusStrand && rcompReads){
			key.flip(r, k);
		}
		if(addName){r.id+=" "+key;}
		
		return key.kmer;
	}
		
	/** Sets the default number of hash iterations for new comparators.
	 * @param x Number of hash iterations (clamped to 0-8) */
	public static synchronized void setHashes(int x){
		defaultHashes=Tools.mid(0, x, 8);
	}
	
	//TODO:  This can be swapped with BBSketch hashing code.  Check speed.
	/**
	 * Hashes a k-mer using multiple rounds of XOR with lookup tables.
	 * Applies the configured number of hash iterations for better distribution.
	 * @param kmer The k-mer to hash
	 * @return Non-negative hash value (high bit cleared)
	 */
	public final long hash(long kmer){
		long code=kmer;
		for(int i=0; i<hashes; i++){//4 only half-hashes; 8 does full hashing
			int x=(int)(kmer&0xFF);
			kmer>>=8;
			code^=codes[i][x];
		}
		return code&Long.MAX_VALUE;
	}
	
	/** Worker thread for parallel k-mer hashing.
	 * Each thread processes a subset of reads using round-robin assignment. */
	private class HashThread extends Thread{
		
		/**
		 * Creates a hash worker thread with assigned parameters.
		 *
		 * @param id_ Thread identifier for round-robin assignment
		 * @param threads_ Total number of worker threads
		 * @param list_ List of reads to process
		 * @param table_ K-mer count table for frequency filtering
		 * @param minCount_ Minimum k-mer count threshold
		 */
		HashThread(int id_, int threads_, ArrayList<Read> list_, KCountArray table_, int minCount_){
			id=id_;
			threads=threads_;
			list=list_;
			table=table_;
			minCount=minCount_;
		}
		
		@Override
		public void run(){
			for(int i=id; i<list.size(); i+=threads){
				hash(list.get(i), table, minCount, true);
			}
		}
		
		/** Thread identifier for round-robin work assignment */
		final int id;
		/** Total number of worker threads */
		final int threads;
		/** List of reads to process */
		final ArrayList<Read> list;
		/** K-mer count table for frequency filtering */
		final KCountArray table;
		/** Minimum k-mer count threshold for acceptance */
		final int minCount;
	}
	
	/**
	 * Gets a thread-local ReadKey instance for temporary use.
	 * Creates a new key if none exists for the current thread.
	 * @return A cleared ReadKey ready for use
	 */
	static ReadKey getLocalKey(){
		ReadKey key=localReadKey.get();
		if(key==null){
			localReadKey.set(key=new ReadKey());
		}
		key.clear();
		return key;
	}
	
	/** K-mer length used for comparison */
	public final int k;

	/** Bit shift amount for k-mer operations (2*k) */
	final int shift;
	/** Bit shift amount for reverse complement operations (shift-2) */
	final int shift2;
	/** Bit mask for extracting k-mer bits */
	final long mask;
	
	/** Random seed for hash function initialization */
	final long seed;
	/** Number of bases to ignore at read ends */
	final int border;
	/** Number of hash iterations to perform */
	final int hashes;

	/** Whether to append key information to read names */
	final boolean addName;
	/** Whether to reverse-complement reads when using minus strand k-mers */
	final boolean rcompReads;
	
	/** Lookup tables for hash function operations */
	private final long[][] codes;
	
	/** Default random seed for new comparator instances */
	static long defaultSeed=1;
	/** Default number of hash iterations for new comparator instances */
	static int defaultHashes=4;
	/** Default border size for new comparator instances */
	static int defaultBorder=1;
	/**
	 * Minimum probability threshold for k-mer acceptance based on quality scores
	 */
	public static float minProb=0f;
	/** Whether to print debugging information during initialization */
	public static boolean verbose=true;

	/** Whether to attempt read merging before k-mer analysis */
	public static boolean mergeFirst=false;
	/**
	 * Whether to fall back to sequence comparison when k-mer signatures are identical
	 */
	public static boolean compareSequence=true;
	
	/**
	 * Thread-local storage for ReadKey instances to avoid object creation overhead
	 */
	public static ThreadLocal<ReadKey> localReadKey=new ThreadLocal<ReadKey>();
	
	/** Array mapping quality scores to probability of correct base call */
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 127);
	/** Array mapping quality scores to inverse probability of correct base call */
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 127);

}
