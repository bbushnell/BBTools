package ddl;

import java.util.Arrays;

/** Optional per-query checkpoints for the unweighted inverted-index scan.
 * Constructed only by diagnostics; ordinary search does not allocate or update it.
 * One query worker owns this mutable scratch. Reset before each scan, capture the
 * initial state and every complete-key boundary, then requireComplete before reading
 * a finished result. Public arrays are borrowed scratch, overwritten by reuse.
 * @author Fischl
 */
public final class DDLIndexScanProgress{

	/** Allocate one reusable list of occupied bins and four fixed checkpoints. */
	public DDLIndexScanProgress(int width){
		if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS){throw new IllegalArgumentException("Invalid sketch width");}
		bins=new int[width];
	}

	/** Resets counters/results and records ascending occupied bins before scan permutation.
	 * Default checkpoints are the ceiling of each quarter of the occupied-key count.
	 */
	public void reset(char[] query){
		if(query==null || query.length!=bins.length){throw new IllegalArgumentException("Query/progress widths differ");}
		keys=0;
		nextCheckpoint=0;
		halfwayOne=halfwayTwo=0;
		for(int bin=0; bin<query.length; bin++){if(query[bin]!=0){bins[keys++]=bin;}}
		for(int i=0; i<4; i++){keysAt[i]=(keys*(i+1)+3)/4;}
		Arrays.fill(bestAt, 0);
		Arrays.fill(truthAt, 0);
		Arrays.fill(leaderAt, -1);
		Arrays.fill(postingsAt, 0);
	}

	/** Copies the first keys entries of order, requiring each occupied bin exactly once. */
	public void reset(char[] query, int[] order){
		reset(query);
		if(order==null || order.length<keys){throw new IllegalArgumentException("Missing occupied-bin permutation");}
		if(seen==null){seen=new boolean[bins.length];}
		Arrays.fill(seen, false);
		for(int i=0; i<keys; i++){
			final int bin=order[i];
			if(bin<0 || bin>=query.length || query[bin]==0 || seen[bin]){
				throw new IllegalArgumentException("Order must visit each occupied bin once");
			}
			seen[bin]=true;
			bins[i]=bin;
		}
	}

	/** Copies four cumulative key-count boundaries after reset and before capture.
	 * Boundaries may coincide when a quadrant contributes no selected keys.
	 */
	public void setBoundaries(int[] boundaries){
		if(boundaries==null || boundaries.length!=4 || boundaries[3]!=keys){
			throw new IllegalArgumentException("Four checkpoint boundaries must end at every occupied key");
		}
		int previous=0;
		for(int i=0; i<4; i++){
			if(boundaries[i]<previous || boundaries[i]>keys){throw new IllegalArgumentException("Checkpoint boundaries must be nested");}
			keysAt[i]=boundaries[i];
			previous=boundaries[i];
		}
	}

	/**
	 * Snapshots completed keys; coincident checkpoints share the same scan state.
	 * Call at consumed=0 and after every key, so no checkpoint can be skipped.
	 * @param work Current unweighted scan scores, touched-family count and posting count
	 * @param consumed Number of complete occupied query keys visited
	 * @param best Largest current hit count
	 * @param leader Family ID for best, or -1 when no family has a hit
	 * @param twoPlus Number of touched families with at least two hits
	 * @param truthFamily Valid score-array index, or a negative value when truth is absent
	 */
	public void capture(DDLFamilySearch.Workspace work, int consumed, int best, int leader, int twoPlus, int truthFamily){
		assert(work!=null && consumed>=0 && consumed<=keys && best>=0 && twoPlus>=0 && twoPlus<=work.size) :
			"Checkpoints require completed keys and nested positive-hit family sets";
		assert((best==0)==(leader<0)) : "An all-zero scan has no positive-score leader";
		while(nextCheckpoint<4 && consumed==keysAt[nextCheckpoint]){
			final int point=nextCheckpoint++;
			bestAt[point]=best;
			truthAt[point]=truthFamily<0 ? 0 : work.scores[truthFamily];
			leaderAt[point]=leader;
			postingsAt[point]=work.postings;
			if(point==1){halfwayOne=work.size; halfwayTwo=twoPlus;}
		}
	}

	/** Positive-score ties count as co-leadership; all-zero ties do not. */
	public boolean truthLeads(int checkpoint){
		if(checkpoint<0 || checkpoint>=4){throw new IllegalArgumentException("Invalid checkpoint");}
		return bestAt[checkpoint]>0 && truthAt[checkpoint]==bestAt[checkpoint];
	}

	/** Require every checkpoint before a caller emits a complete query row. */
	public void requireComplete(){
		if(nextCheckpoint!=4){throw new IllegalStateException("Incomplete diagnostic scan: "+nextCheckpoint+" checkpoints");}
	}

	/** Occupied-bin order; only the first keys entries are populated for this query. */
	public final int[] bins;
	/** Cumulative key counts, best/truth hit counts and leader IDs at four checkpoints. */
	public final int[] keysAt=new int[4], bestAt=new int[4], truthAt=new int[4], leaderAt=new int[4];
	/** Posting entries visited at each checkpoint, not distinct-family counts. */
	public final long[] postingsAt=new long[4];
	/** Occupied keys, then families with at least one/two hits at the second checkpoint. */
	public int keys, halfwayOne, halfwayTwo;
	/** Next checkpoint to capture, from zero through four. */
	private int nextCheckpoint;
	/** Reused duplicate-detection scratch for an explicit occupied-bin permutation. */
	private boolean[] seen;
}
