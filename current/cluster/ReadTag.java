package cluster;

import java.io.Serializable;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Mar 24, 2014
 *
 */
class ReadTag implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -6186366525723397478L;

	/**
	 * Constructs a ReadTag wrapper for the given read.
	 * Computes GC content and initializes strand information.
	 * @param r_ The read to wrap with clustering metadata
	 */
	public ReadTag(Read r_){
		r=r_;
		strand=r.strand();

		int gcCount_=0;
		for(byte b : r.bases){
			if(b=='G' || b=='C'){
				gcCount_++;
			}
		}
		gcCount=gcCount_;
		
		processHeader(r.id);
	}
	
	/**
	 * Processes read header to extract additional metadata.
	 * Currently unimplemented - sets default values for gc, depth, and cluster0.
	 * @param s The read header string to process
	 */
	private void processHeader(String s){
		assert(false) : "TODO";
		gc=-1;
		depth=-1;
		cluster0=-1;
	}

	/**
	 * Returns the first read in a pair.
	 * For single reads or forward reads (strand=0), returns the read itself.
	 * For reverse reads (strand=1), returns the mate.
	 * @return The first read in the pair
	 */
	Read r1(){
		return strand==0 ? r : r.mate;
	}
	
	/**
	 * Returns the second read in a pair.
	 * For reverse reads (strand=1), returns the read itself.
	 * For forward reads (strand=0), returns the mate.
	 * @return The second read in the pair, or null if single-ended
	 */
	Read r2(){
		return strand==1 ? r : r.mate;
	}
	
	/** Returns the ReadTag object for the first read in the pair.
	 * @return The ReadTag attached to the first read */
	ReadTag tag1(){
		return (ReadTag)r1().obj;
	}
	
	/** Returns the ReadTag object for the second read in the pair.
	 * @return The ReadTag attached to the second read, or null if single-ended */
	ReadTag tag2(){
		Read r2=r2();
		return r2==null ? null : (ReadTag)r2.obj;
	}
	
//	private int[] toKmers(final int k){
//		return ClusterTools.toKmers(r.bases, null, k);
//	}
	
	/**
	 * Returns cached sorted k-mer array for the specified k-mer length.
	 * Lazy initialization - computes and caches on first access.
	 * @param k1 The k-mer length to use
	 * @return Array of k-mers from the read sequence
	 */
	int[] kmerArray1(int k1){
		if(kmerArray1==null){kmerArray1=ClusterTools.toKmers(r.bases, null, k1);}
		return kmerArray1;
	}
	
	/**
	 * Returns cached k-mer count array for the specified k-mer length.
	 * Uses canonical k-mer ordering and lazy initialization.
	 * @param k2 The k-mer length to use for counting
	 * @return Array of k-mer counts in canonical order
	 */
	int[] kmerArray2(int k2){
		if(kmerArray2==null){kmerArray2=ClusterTools.toKmerCounts(r.bases, null, k2);}
		return kmerArray2;
	}
	
	/**
	 * Returns normalized k-mer frequency array for the specified k-mer length.
	 * Converts k-mer counts to frequencies with smoothing (95% of actual frequency
	 * plus 0.05% distributed equally among all k-mers).
	 *
	 * @param k2 The k-mer length to use
	 * @return Array of smoothed k-mer frequencies
	 */
	float[] kmerFreq2(int k2){
		if(kmerFreq2==null){
			int[] counts=kmerArray2(k2);
			if(counts!=null){
				long sum=shared.Vector.sum(counts);
				kmerFreq2=new float[counts.length];
				float extra=(0.05f/counts.length);
				float mult=0.95f/sum;
				for(int i=0; i<counts.length; i++){
					kmerFreq2[i]=counts[i]*mult+extra;
				}
			}
		}
		return kmerFreq2;
	}
	
	/** Sorted long kmers */
	private int[] kmerArray1;
	
	/** Canonically-ordered short kmer counts */
	private int[] kmerArray2;
	
	/** Cached normalized k-mer frequency array */
	private float[] kmerFreq2;
	
	/** The wrapped read object */
	final Read r;
	/** Strand orientation of the read (0=forward, 1=reverse) */
	final byte strand;
	/** Count of G and C bases in the read sequence */
	final int gcCount;
	
	/** Sequencing depth information for the read */
	int depth;
	/** Initial cluster assignment identifier */
	int cluster0=-1; //initial cluster
	/** Final cluster assignment identifier */
	int cluster1=-1; //final cluster

	/** GC content fraction of the read sequence */
	float gc;
	
}
