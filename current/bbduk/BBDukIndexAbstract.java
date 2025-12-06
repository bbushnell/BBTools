package bbduk;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLongArray;

import aligner.SideChannel3;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.IntList;
import var2.ScafMap;
import var2.VarMap;

public abstract class BBDukIndexAbstract{

	abstract boolean loaded();

	abstract void loadAssorted(String in1_for_header);

	abstract void cleanup();

	/**
	 * Clear stored kmers.
	 */
	abstract void unloadKmers();

	/**
	 * Clear stored sequence data.
	 */
	abstract void unloadScaffolds();

	abstract int ways();

	abstract long addToMap(long kmer, long rkmer, int k, long extraBase, int id, long kmask,
		int hammingDistance, int editDistance, int tnum);

	abstract long addToMapRightShift(long kmer, long rkmer, int id, int tnum);

	abstract long addToMapLeftShift(long kmer, long rkmer, long extraBase, int id, int tnum);

	abstract void rebalance(int tnum);

	abstract void dump(ByteStreamWriter bsw, int i, int maxValue);
	

	abstract void setKmersLoaded();

	abstract void toRefNames();

	abstract void dump(TextStreamWriter tsw, int i, int maxValue);

	abstract SideChannel3 sidechannel();

	abstract int getValue(long kmer, long rkmer, long l, int i, int len,
		int qHammingDistance2);

	abstract void writeRefStats(String in1, String in2, long readsIn);
	
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	ArrayList<String> refNames;
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
	
	VarMap varMap;
	ScafMap scafMap;
	boolean fixVariants;

	/** Has this class encountered errors while processing? */
	public boolean errorState;
	
}