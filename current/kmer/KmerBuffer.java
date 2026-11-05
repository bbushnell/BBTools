package kmer;

import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * Thread-confined buffer of packed keys, optionally paired with integer values.
 * Storage grows as needed; the caller decides when to flush. Clearing retains
 * capacity. Count-only and key/value insertion modes must not be mixed.
 *
 * @author Brian Bushnell
 * @date Jul 30, 2015
 */
public class KmerBuffer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Allocates key storage and optional value storage.
	 * @param buflen Initial list capacity, not a flush threshold
	 * @param k_ Key width used only by the diagnostic string representation
	 * @param initValues Whether every inserted key must have an associated value */
	public KmerBuffer(int buflen, int k_, boolean initValues){
		k=k_;
		kmers=new LongList(buflen);
		values=(initValues ? new IntList(buflen) : null);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Appends one key in count-only mode; returns the occupied key count. */
	public int add(long kmer){
		assert(values==null);
		kmers.add(kmer);
		assert(values==null);
		return kmers.size;
	}
	
	/** Appends {@code times} copies of a key in count-only mode; nonpositive counts add nothing.
	 * @return Occupied key count after insertion */
	public int addMulti(long kmer, int times){
		assert(values==null);
		for(int i=0; i<times; i++){kmers.add(kmer);}
		assert(values==null);
		return kmers.size;
	}
	
	/** Appends a key/value pair in value mode; returns the occupied key count. */
	public int add(long kmer, int value){
		kmers.add(kmer);
		values.add(value);
		assert(values.size==kmers.size);
		return kmers.size;
	}
	
	/** Clears logical contents while retaining both lists' backing storage. */
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	/** Returns the number of occupied keys, excluding unused backing capacity. */
	final int size(){return kmers.size;}
	
	/** Decodes occupied keys at the configured width, separated by commas; omits values. */
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<kmers.size; i++){
			if(i>0){bb.append(',');}
			bb.appendKmer(kmers.get(i), k);
		}
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Diagnostic decoding width; does not alter stored key bits. */
	private final int k;
	/** Reusable key storage, owned by this buffer. */
	final LongList kmers;
	/** Parallel value storage, or null in count-only mode. */
	final IntList values;
	
}
