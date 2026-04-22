package cardinality;

import parse.Parser;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Manages an array of CardinalityTrackers for simultaneous cardinality
 * estimation at multiple k-mer lengths. Validates, sorts, and deduplicates
 * the k-mer list on construction.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MultiLogLog {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Construct from parsed command-line arguments. */
	public MultiLogLog(Parser p){
		this(p.loglogbuckets, p.loglogseed, p.loglogMinprob, p.loglogKlist);
	}
	
	/**
	 * Full constructor. Filters klist0 to valid lengths, sorts, deduplicates,
	 * and creates one CardinalityTracker per k-mer length.
	 * @param buckets Number of buckets per tracker
	 * @param seed Random seed for hash functions
	 * @param minProb Minimum probability threshold
	 * @param klist0 List of desired k-mer lengths
	 */
	public MultiLogLog(int buckets, long seed, float minProb, IntList klist0){
		assert(klist0.size>0) : "No valid kmer lengths specified.";
		final IntList klist=new IntList(klist0.size);
		for(int i=0; i<klist0.size; i++){
			final int x=klist0.get(i);
			final int k=Kmer.getKbig(x);
			if(k>0){klist.add(k);}
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
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Hash a read through all trackers. */
	public void hash(Read r){
		for(CardinalityTracker c : counters){c.hash(r);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Valid k-mer lengths, sorted and deduplicated. */
	public final int[] kArray;
	/** One CardinalityTracker per k-mer length. */
	public final CardinalityTracker[] counters;
	
}
