package cardinality;

import shared.Parser;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

public class MultiLogLog {
	
	public MultiLogLog(Parser p){
		this(p.loglogbuckets, p.loglogseed, p.loglogMinprob, p.loglogKlist);
	}
	
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
	
	public void hash(Read r){
		for(CardinalityTracker c : counters){
			c.hash(r);
		}
	}
	
	public final int[] kArray;
	public final CardinalityTracker[] counters;
	
}
