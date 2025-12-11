package cluster;

import java.io.Serializable;

import stream.Read;

/**
 * Wrapper class for sequence reads that provides clustering metadata and k-mer analysis.
 * Stores read information along with computed properties like GC content, cluster assignments,
 * and cached k-mer representations for efficient clustering operations.
 *
 * @author Brian Bushnell
 * @date Mar 24, 2014
 */
class ReadTag implements Serializable{
	
	/** Serialization version identifier for compatibility */
	private static final long serialVersionUID = -6186366525723397478L;

	/** Wraps a read, records strand/GC count, and processes header metadata for clustering.
	 * @param r_ Read to wrap */
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
	
	private void processHeader(String s){
		assert(false) : "TODO";
		gc=-1;
		depth=-1;
		cluster0=-1;
	}

	/** Returns the first read in the pair (forward read or this if unpaired).
	 * @return First read */
	Read r1(){
		return strand==0 ? r : r.mate;
	}
	
	/** Returns the second read in the pair (reverse read or mate), or null if single-ended.
	 * @return Second read or null */
	Read r2(){
		return strand==1 ? r : r.mate;
	}
	
	ReadTag tag1(){
		return (ReadTag)r1().obj;
	}
	
	ReadTag tag2(){
		Read r2=r2();
		return r2==null ? null : (ReadTag)r2.obj;
	}
	
//	private int[] toKmers(final int k){
//		return ClusterTools.toKmers(r.bases, null, k);
//	}
	
	/**
	 * Returns cached sorted k-mers of length k1, computing and caching on first access.
	 * @param k1 K-mer length
	 * @return Array of k-mers
	 */
	int[] kmerArray1(int k1){
		if(kmerArray1==null){kmerArray1=ClusterTools.toKmers(r.bases, null, k1);}
		return kmerArray1;
	}
	
	/**
	 * Returns cached canonical k-mer counts for length k2, computing on first access.
	 * @param k2 K-mer length
	 * @return Array of k-mer counts
	 */
	int[] kmerArray2(int k2){
		if(kmerArray2==null){kmerArray2=ClusterTools.toKmerCounts(r.bases, null, k2);}
		return kmerArray2;
	}
	
	/**
	 * Returns smoothed k-mer frequency array (95% observed, 5% evenly distributed) for length k2.
	 * @param k2 K-mer length
	 * @return Array of k-mer frequencies, or null if counts unavailable
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
	
	private int[] kmerArray1;
	
	private int[] kmerArray2;
	
	private float[] kmerFreq2;
	
	final Read r;
	final byte strand;
	final int gcCount;
	
	int depth;
	int cluster0=-1; //initial cluster
	int cluster1=-1; //final cluster

	float gc;
	
}
