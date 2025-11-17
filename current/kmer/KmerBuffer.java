package kmer;

import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 30, 2015
 *
 */
public class KmerBuffer {
	
	/**
	 * Constructs a KmerBuffer with specified capacity and optional value tracking.
	 * @param buflen Initial capacity for the k-mer storage lists
	 * @param k_ Length of k-mers to be stored
	 * @param initValues If true, creates values list for storing associated integers
	 */
	public KmerBuffer(int buflen, int k_, boolean initValues){
		k=k_;
		kmers=new LongList(buflen);
		values=(initValues ? new IntList(buflen) : null);
	}
	
	/**
	 * Adds a k-mer to the buffer without an associated value.
	 * Only works when values tracking is disabled.
	 * @param kmer The k-mer encoded as a long
	 * @return Current size of the k-mer list after addition
	 */
	public int add(long kmer){
		assert(values==null);
		kmers.add(kmer);
		assert(values==null);
		return kmers.size;
	}
	
	/**
	 * Adds the same k-mer multiple times to the buffer.
	 * Only works when values tracking is disabled.
	 *
	 * @param kmer The k-mer encoded as a long
	 * @param times Number of times to add the k-mer
	 * @return Current size of the k-mer list after all additions
	 */
	public int addMulti(long kmer, int times){
		assert(values==null);
		for(int i=0; i<times; i++){kmers.add(kmer);}
		assert(values==null);
		return kmers.size;
	}
	
	/**
	 * Adds a k-mer with an associated integer value to the buffer.
	 * Requires values tracking to be enabled during construction.
	 *
	 * @param kmer The k-mer encoded as a long
	 * @param value Integer value to associate with the k-mer
	 * @return Current size of the k-mer list after addition
	 */
	public int add(long kmer, int value){
		kmers.add(kmer);
		values.add(value);
		assert(values.size==kmers.size);
		return kmers.size;
	}
	
	/** Removes all k-mers and values from the buffer, resetting size to zero.
	 * Clears both the k-mer list and values list if values tracking is enabled. */
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	//Returns raw size of kmers array, rather than actual number of kmers
	/**
	 * Returns the raw size of the k-mers array.
	 * Note: This returns the actual number of k-mers stored, not array capacity.
	 * @return Number of k-mers currently stored in the buffer
	 */
	final int size(){return kmers.size;}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<kmers.size; i++){
			if(i>0){bb.append(',');}
			bb.appendKmer(kmers.get(i), k);
		}
		return bb.toString();
	}
	
	/** Length of k-mers stored in this buffer */
	private final int k;
	/** Dynamic list storing k-mers encoded as long values */
	final LongList kmers;
	/** Optional dynamic list storing integer values associated with k-mers */
	final IntList values;
	
}
