package bbduk;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLongArray;

import aligner.SideChannel3;
import fileIO.ByteStreamWriter;
import structures.IntList;
import var2.ScafMap;
import var2.VarMap;

/**
 * Abstract base class for BBDuk indices.
 * Defines the contract for kmer storage/retrieval and holds shared result state
 * (scaffold counts, reference names) populated during loading.
 * Implementations (e.g. Modulo, Mask) are responsible for managing their own
 * configuration parameters (k, dist, etc) and storage structures.
 * @author Brian Bushnell
 * @date November 19, 2025
 */
public abstract class BBDukIndex {
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** @return True if the index is fully loaded and ready for queries. */
	abstract boolean loaded();

	/**
	 * Loads assorted non-kmer data (reference sequences, sam headers, variants).
	 * @param in1_for_header Input file name for header parsing if needed.
	 */
	abstract void loadAssorted(String in1_for_header);

	/** Releases memory used by the index (kmers and scaffolds). */
	abstract void cleanup();

	/** Clear stored kmers to save memory. */
	abstract void unloadKmers();

	/** Clear stored sequence data to save memory. */
	abstract void unloadScaffolds();

	/** @return The number of tables/ways in the index. */
	abstract int ways();

	/**
	 * Adds a kmer to the table.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param k Kmer length
	 * @param extraBase Base added to end in case of deletions
	 * @param id Scaffold ID
	 * @param kmask Bitmask for kmer length
	 * @param hammingDistance Max subs allowed
	 * @param editDistance Max edits allowed
	 * @param tnum Thread/Table ID
	 * @return Number of kmers stored
	 */
	abstract long addToMap(long kmer, long rkmer, int k, long extraBase, int id, long kmask,
		int hammingDistance, int editDistance, int tnum);

	/** Adds short kmers on the right end of the read. */
	abstract long addToMapRightShift(long kmer, long rkmer, int id, int tnum);

	/** Adds short kmers on the left end of the read. */
	abstract long addToMapLeftShift(long kmer, long rkmer, long extraBase, int id, int tnum);

	/** Rebalances the specified table if necessary. */
	abstract void rebalance(int tnum);

	/** Dumps kmers to a ByteStreamWriter. */
	abstract void dump(ByteStreamWriter bsw, int i, int maxValue);

	/** Marks kmers as loaded. */
	abstract void setKmersLoaded();

	/** Fills scaffold names from reference file names. */
	abstract void toRefNames();

	/** @return The SideChannel3 instance for alignment, or null. */
	abstract SideChannel3 sidechannel();

	/**
	 * Queries the index for a kmer.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param l Kmer length mask
	 * @param i Query position
	 * @param len Kmer length
	 * @param qHammingDistance2 Hamming distance
	 * @return Scaffold ID or -1.
	 */
	abstract int getValue(long kmer, long rkmer, long l, int i, int len,
		int qHammingDistance2);

	/** Writes reference statistics to a file. */
	abstract void writeRefStats(String in1, String in2, long readsIn);
	
	/*--------------------------------------------------------------*/
	/*----------------        Shared Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	ArrayList<String> refNames;
	/** Alternate reference names (if primary ref fails). */
	ArrayList<String> altRefNames;
	/** Number of scaffolds per reference. */
	int[] refScafCounts;
	/** Array of reference files from which to load kmers */
	String[] ref;
	/** Alternate reference to be used if main reference has no kmers */
	String[] altref;
	/** Array of literal strings from which to load kmers */
	String[] literal;
	/** Optional reference for sam file */
	String samref;

	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	AtomicLongArray scaffoldReadCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	AtomicLongArray scaffoldBaseCounts;
	/** scaffoldLengths[id] stores the length of that scaffold */
	IntList scaffoldLengths=new IntList();
	
	/** Number of reference reads processed for kmer loading */
	long refReads;
	/** Number of reference bases processed for kmer loading */
	long refBases;
	/** Number of reference kmers encountered during loading */
	long refKmers;
	/** Number of unique kmers actually stored in hash tables */
	long storedKmers;
	
	/** Variant map (if loaded). */
	VarMap varMap;
	/** Scaffold map (if loaded). */
	ScafMap scafMap;
	/** True if variants should be fixed */
	boolean fixVariants;

	/** Has this class encountered errors while processing? */
	public boolean errorState;
	
}