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

	Read r1(){
		return strand==0 ? r : r.mate;
	}
	
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
	
	int[] kmerArray1(int k1){
		if(kmerArray1==null){kmerArray1=ClusterTools.toKmers(r.bases, null, k1);}
		return kmerArray1;
	}
	
	int[] kmerArray2(int k2){
		if(kmerArray2==null){kmerArray2=ClusterTools.toKmerCounts(r.bases, null, k2);}
		return kmerArray2;
	}
	
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
	
	private float[] kmerFreq2;
	
	final Read r;
	final byte strand;
	final int gcCount;
	
	int depth;
	int cluster0=-1; //initial cluster
	int cluster1=-1; //final cluster

	float gc;
	
}
