package cardinality;

import shared.Parser;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Manages multiple LogLog cardinality trackers for simultaneous k-mer size tracking
 * in genomic sequencing data. Creates and manages an array of CardinalityTracker
 * instances to estimate unique k-mer counts across different k-mer lengths
 * simultaneously. Validates and filters k-mer lengths dynamically, automatically
 * sorting and deduplicating k-mer lengths for efficiency.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MultiLogLog {
	
	/**
	 * Creates a MultiLogLog instance from parser configuration.
	 * Uses parser's loglog parameters including buckets, seed, minimum probability, and k-mer list.
	 * @param p Parser containing loglog configuration parameters
	 */
	public MultiLogLog(Parser p){
		this(p.loglogbuckets, p.loglogseed, p.loglogMinprob, p.loglogKlist);
	}
	
	/**
	 * Creates a MultiLogLog instance with explicit parameters.
	 * Filters k-mer list to keep only valid lengths, sorts and deduplicates them.
	 * Creates CardinalityTracker array with specified configuration for each k-mer length.
	 *
	 * @param buckets Number of buckets for each cardinality tracker
	 * @param seed Random seed for hash functions
	 * @param minProb Minimum probability threshold for k-mer quality filtering
	 * @param klist0 List of desired k-mer lengths
	 */
	public MultiLogLog(int buckets, long seed, float minProb, IntList klist0){
		assert(klist0.size>0) : "No valid kmer lengths specified.";
		IntList klist=new IntList(klist0.size);
		for(int i=0; i<klist0.size; i++){
			int x=klist0.get(i);
			int k=Kmer.getKbig(x);
			if(k>0){
				klist.add(k);
			}
		}
		klist.sort();
		klist.shrinkToUnique();
		assert(klist.size>0) : "No valid kmer lengths specified.";
		kArray=klist.toArray();
		counters=new LogLog[kArray.length];
		for(int i=0; i<kArray.length; i++){
			counters[i]=CardinalityTracker.makeTracker(buckets, kArray[i], seed, minProb);
		}
	}
	
	/**
	 * Hashes a read across all cardinality trackers.
	 * Processes the read through each tracker to update cardinality estimates for all k-mer lengths.
	 * @param r The read to hash and track
	 */
	public void hash(Read r){
		for(CardinalityTracker c : counters){
			c.hash(r);
		}
	}
	
	/** Array of valid k-mer lengths used by the cardinality trackers */
	public final int[] kArray;
	/** Array of cardinality trackers, one for each k-mer length */
	public final CardinalityTracker[] counters;
	
}
