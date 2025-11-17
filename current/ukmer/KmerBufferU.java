package ukmer;

import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 9, 2015
 *
 */
public class KmerBufferU {
	
	/**
	 * Constructs a k-mer buffer with specified capacity and configuration.
	 * Initializes internal storage based on k-mer length and value tracking requirements.
	 *
	 * @param buflen Buffer capacity in number of k-mers
	 * @param kbig K-mer bit length for determining storage multiplier
	 * @param initValues Whether to initialize value storage alongside k-mer storage
	 */
	public KmerBufferU(int buflen, int kbig, boolean initValues){
		k=Kmer.getK(kbig);
		mult=Kmer.getMult(kbig);
		kmers=new LongList(buflen*mult);
		values=(initValues ? new IntList(buflen) : null);
	}
	
	/**
	 * Adds a k-mer to the buffer using its key representation.
	 * @param kmer The k-mer object to add
	 * @return The new size of the k-mer storage after addition
	 */
	public int add(Kmer kmer){
//		System.err.println("Adding "+kmer+"; this="+this+"; kmers.size="+kmers.size);
		add(kmer.key());
		return kmers.size;
//		System.err.println("Added "+kmer+"; this="+this+"; kmers.size="+kmers.size);
	}
	
	/**
	 * Adds a k-mer with an associated value to the buffer.
	 * @param kmer The k-mer object to add
	 * @param value The integer value to associate with the k-mer
	 */
	public void add(Kmer kmer, int value){
		add(kmer.key(), value);
	}
	
	/**
	 * Adds a k-mer represented as a long array to the buffer.
	 * Used when values are not tracked (values must be null).
	 * @param kmer The k-mer as a long array representation
	 */
	public void add(long[] kmer){
		assert(values==null);
		assert(kmer.length==mult) : kmer.length+", "+mult+", "+k;
		kmers.append(kmer);
	}
	
	/**
	 * Adds a k-mer represented as a long array with an associated value.
	 * Maintains synchronization between k-mer and value storage arrays.
	 * @param kmer The k-mer as a long array representation
	 * @param value The integer value to associate with the k-mer
	 */
	public void add(long[] kmer, int value){
		assert(kmer.length==mult);
		kmers.append(kmer);
		values.add(value);
		assert(values.size*mult==kmers.size);
	}
	
	/** Clears all k-mers and values from the buffer.
	 * Resets both k-mer and value storage to empty state. */
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	//Returns raw size of kmers array, rather than actual number of kmers
	/**
	 * Returns the raw size of the k-mer storage array.
	 * Note: This is not the actual number of k-mers, but the array size
	 * which includes multiplier effects for long k-mer storage.
	 * @return Raw size of the underlying k-mer storage array
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
	
	/** Storage multiplier for k-mer representation based on k-mer length */
	private final int mult;
	/** K-mer length in bases */
	private final int k;
	/** Storage for k-mer data as long values */
	final LongList kmers;
	/**
	 * Optional storage for values associated with k-mers; null if not tracking values
	 */
	final IntList values;
	
}
