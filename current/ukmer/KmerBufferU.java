package ukmer;

import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * Thread-confined buffer that copies consecutive key words into reusable storage.
 * A logical key occupies mult long entries and, in value mode, one integer value.
 * The caller controls flushing; this class neither writes tables nor enforces a
 * flush threshold. Packing settings must agree with the supplied keys.
 *
 * @author Brian Bushnell
 * @date Jul 9, 2015
 */
public class KmerBufferU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Allocates storage using the current Kmer packing mode.
	 * @param buflen Initial capacity in logical keys
	 * @param kbig Requested key length, normalized by Kmer's geometry helpers
	 * @param initValues Whether every inserted key must have a value */
	public KmerBufferU(int buflen, int kbig, boolean initValues){
		k=Kmer.getK(kbig);
		mult=Kmer.getMult(kbig);
		kmers=new LongList(buflen*mult);
		values=(initValues ? new IntList(buflen) : null);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Copies the selected key in count-only mode; returns occupied long words.
	 * Kmer.key() selects the canonical strand only when canonicalization is enabled. */
	public int add(Kmer kmer){
//		System.err.println("Adding "+kmer+"; this="+this+"; kmers.size="+kmers.size);
		add(kmer.key());
		return kmers.size;
//		System.err.println("Added "+kmer+"; this="+this+"; kmers.size="+kmers.size);
	}
	
	/** Copies the selected key and appends one value; requires value mode.
	 * Kmer.key() selects the canonical strand only when canonicalization is enabled. */
	public void add(Kmer kmer, int value){add(kmer.key(), value);}
	
	/** Copies exactly mult key words in count-only mode; does not retain the array. */
	public void add(long[] kmer){
		assert(values==null);
		assert(kmer.length==mult) : kmer.length+", "+mult+", "+k;
		kmers.append(kmer);
	}
	
	/** Copies exactly mult words and one value; requires value mode. */
	public void add(long[] kmer, int value){
		assert(kmer.length==mult);
		kmers.append(kmer);
		values.add(value);
		assert(values.size*mult==kmers.size);
	}
	
	/** Clears logical contents while retaining backing arrays. */
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	/** Returns occupied long words, not logical keys or allocated capacity. */
	final int size(){return kmers.size;}
	
	/** Decodes each occupied word at k bases, separated by commas.
	 * This diagnostic view pads a partial last word; it is not a full-key serializer. */
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
	
	/** Number of stored long words per logical key. */
	private final int mult;
	/** Width used to decode every word in the diagnostic view. */
	private final int k;
	/** Owned word storage with adjacent keys concatenated. */
	final LongList kmers;
	/** One value per logical key, or null in count-only mode. */
	final IntList values;
	
}
