package assemble;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Observation-only trace of actual contig-walk junction checks. Supplied counts
 * are not modified and no return value is supplied to the walker. Repeated encounters
 * are retained explicitly; corpus construction must deduplicate canonical words.
 * Reference truth is deliberately unavailable in this process.
 * Callers select encounters with sample before record; record does not sample or deduplicate.
 * Nonbranch sampling is stateless and repeats the same decision for the same key or
 * key hash; it may run concurrently. Native stops bypass that sampling decision.
 * Record, completion and close serialize shared scratch/output access. Callers must keep
 * supplied neighborhoods and oriented Kmer state stable until record returns. Row order
 * follows observer-lock acquisition, not invocation order or canonical-word order across workers.
 * @author Fischl
 */
public final class LiveBranchObserver implements AutoCloseable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens a new trace with fixed K/depth/control settings; alias validation belongs to the driver.
	 * Starts the writer and supplies version/column headers before returning; output may still fail.
	 * The driver imposes its own supported K range; this constructor checks only K>=2
	 * and positive depth/sampling settings.
	 * @param path Trace destination; append and overwrite are disabled
	 * @param k_ Word length in bases
	 * @param minExtend_ Minimum extension count for the recorded policy screen
	 * @param controlEvery_ Hash-sampling denominator for native nonbranch encounters
	 */
	LiveBranchObserver(final String path, final int k_, final int minExtend_, final int controlEvery_){
		if(path==null || k_<2 || minExtend_<1 || controlEvery_<1){
			throw new IllegalArgumentException("Trace needs a path, K>=2, mce>=1 and positive control interval.");
		}
		k=k_; minExtend=minExtend_; controlEvery=controlEvery_;
		//TODO: Probable bug - startup/header failure after opening writer leaves
		//no constructed observer for the driver's try-with-resources to close.
		writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		writer.print("#livebranch-observation-v2\n");
		writer.print(HEADER+"\n");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * All native branch/hidden stops enter the trace; native nonbranch nodes are
	 * sampled deterministically by canonical table key. This is not an acceptance
	 * probability and must not be used as a model feature. Mix first so nearby
	 * packed DNA keys do not preferentially select one terminal nucleotide.
	 * @param canonicalKey Canonical short-key value or long-K table hash supplied by the walker
	 * @param stopped Whether either native branch flag or the hidden-branch flag is set
	 * @return True for every stop, or when the mixed nonnegative key is divisible by controlEvery
	 */
	boolean sample(final long canonicalKey, final boolean stopped){
		if(stopped || controlEvery==1){return true;}
		long hash=canonicalKey;
		hash=(hash^(hash>>>33))*0xff51afd7ed558ccdL;
		hash=(hash^(hash>>>33))*0xc4ceb9fe1a85ec53L;
		hash^=hash>>>33;
		return (hash&Long.MAX_VALUE)%controlEvery==0;
	}

	/**
	 * Decodes the oriented short key into reusable scratch and queues one canonicalized encounter.
	 * Borrows the four-element count arrays only for this call; their indices are A/C/G/T
	 * relative to forward. The caller must sample first and keep both arrays stable.
	 * @param forward Oriented packed DNA word, not its canonical table key; K must be at most 31
	 * @param count Current node's read-derived count
	 * @param left Left-extension counts in traversal orientation
	 * @param right Right-extension counts in traversal orientation
	 * @param lb Native left-branch flag before any model rescue
	 * @param rb Native right-branch flag before any model rescue
	 * @param arrival Evicted predecessor base (0=A through 3=T), or -1 at entry
	 * @param seed Whether this is a walk-entry or reentry observation
	 * @param hidden Native hidden-branch flag, including the walker's branch-policy conditions
	 */
	synchronized void record(final long forward, final int count, final int[] left, final int[] right,
			final boolean lb, final boolean rb, final int arrival, final boolean seed, final boolean hidden){
		assert(k<=31) : "Packed short-key decoding is only valid for the Tadpole1 representation.";
		word.clear().appendKmer(forward, k);
		write(count, left, right, lb, rb, arrival, seed, hidden);
	}

	/**
	 * Queues the same encounter schema from a long-K traversal word. Decodes array1 explicitly:
	 * AbstractKmerTableU.toText uses canonical key(), not traversal orientation. Uses each packed
	 * word's actual base count, including a shorter final word. Borrows forward and both
	 * four-element A/C/G/T count arrays without changing them; keep them stable for this call.
	 * @param forward Oriented Kmer whose kbig equals this observer's K
	 * @param count Current node's read-derived count
	 * @param left Left-extension counts in traversal orientation
	 * @param right Right-extension counts in traversal orientation
	 * @param lb Native left-branch flag before any model rescue
	 * @param rb Native right-branch flag before any model rescue
	 * @param arrival Evicted predecessor base (0=A through 3=T), or -1 at entry
	 * @param seed Whether this is a walk-entry or reentry observation
	 * @param hidden Native hidden-branch flag, including the walker's branch-policy conditions
	 */
	synchronized void record(final Kmer forward, final int count, final int[] left, final int[] right,
			final boolean lb, final boolean rb, final int arrival, final boolean seed, final boolean hidden){
		assert(forward.kbig==k) : "Observer and live count table must use the same K.";
		word.clear();
		final long[] words=forward.array1();
		for(int i=0; i<words.length; i++){word.appendKmer(words[i], forward.perWordK(i));}
		write(count, left, right, lb, rb, arrival, seed, hidden);
	}

	/**
	 * Canonicalizes the sequence together with base-indexed counts and arrival side.
	 * Seed/reentry and hidden flags are provenance, not direction-invariant model
	 * inputs. An entry record is not necessarily a length-K seed after left reversal.
	 * Reverse orientation swaps neighborhoods and complements base indices; branch flags
	 * are swapped. An interior arrival changes side 0 to 1 and complements its base;
	 * entry arrival stays -1. Hidden status is separate from the symmetric policy level.
	 * Requires the observer monitor and an already decoded word. The recorded policy uses
	 * minExtend and a fixed competitor fraction of 0.5; it does not change native flags.
	 * Stop rows carry sampling denominator 1; other rows carry controlEvery. A queued row
	 * increments rows, but asynchronous writer failure may still prevent durable output.
	 */
	private void write(final int count, final int[] left, final int[] right, final boolean lb,
			final boolean rb, final int arrival, final boolean seed, final boolean hidden){
		if(closed || completed){throw new IllegalStateException("Trace was used after its terminal marker.");}
		if(arrival< -1 || arrival>3){throw new IllegalArgumentException("Actual predecessor must be a base or -1 at entry.");}
		assert(word.length()==k) : "The trace must identify exactly the node whose neighbors were measured.";
		final int level=LiveBranchPolicy.level(count, left, right, lb, rb, minExtend, 0.5);
		boolean reverse=false;
		for(int i=0; i<k; i++){
			final int a=AminoAcid.baseToNumber[word.array[i]];
			final int b=3-AminoAcid.baseToNumber[word.array[k-1-i]];
			if(a<b){break;}
			if(a>b){reverse=true; break;}
		}
		row.clear();
		if(reverse){
			for(int i=k-1; i>=0; i--){row.append(AminoAcid.numberToBase[3-AminoAcid.baseToNumber[word.array[i]]]);}
		}else{row.append(word);}
		row.tab().append(k).tab().append(count);
		for(int i=0; i<4; i++){row.tab().append(reverse ? right[3-i] : left[i]);}
		for(int i=0; i<4; i++){row.tab().append(reverse ? left[3-i] : right[i]);}
		row.tab().append((reverse ? rb : lb) ? 1 : 0).tab().append((reverse ? lb : rb) ? 1 : 0);
		row.tab().append(arrival<0 ? -1 : (reverse ? 1 : 0));
		row.tab().append(arrival<0 ? -1 : (reverse ? 3-arrival : arrival));
		row.tab().append(seed ? 1 : 0).tab().append(hidden ? 1 : 0).tab().append(level);
		row.tab().append(lb || rb || hidden ? 1 : controlEvery).nl();
		writer.print(row);
		rows++;
	}

	/**
	 * Serializes supplied table moments and queued encounter count once. The caller must finish
	 * walk workers first; this method does not join them or recompute the moments.
	 * Enqueuing the trailer prevents further records, but successful close is still required.
	 * @param moments Borrowed final table moments, stable for this call; not the sampled-row population
	 */
	synchronized void complete(final SameKJoinFeatures.Moments moments){
		if(closed || completed || moments==null){throw new IllegalStateException("Trace completion must occur exactly once after workers finish.");}
		row.clear().append("#complete\t").append(rows).tab().append(moments.unique);
		row.tab().appendSlow(moments.observations).tab().appendSlow(moments.squares);
		row.tab().append(controlEvery).tab().append(k).tab().append(minExtend).tab().appendSlow(0.5).nl();
		writer.print(row);
		completed=true;
	}

	/**
	 * Drains/stops the writer without adding a completion trailer. Marks closed before waiting,
	 * so subsequent calls return even if the first close reported an output failure.
	 * A trailer already enqueued by complete does not by itself certify successful output.
	 */
	@Override
	public synchronized void close(){
		if(closed){return;}
		closed=true;
		if(writer.poisonAndWait()){throw new IllegalStateException("Live branch trace output failed.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed word length, rescue minimum depth and deterministic control-sampling denominator. */
	private final int k, minExtend, controlEvery;
	/** Owned output; successful construction must be paired with close. */
	private final ByteStreamWriter writer;
	/** Oriented-word and serialized-row scratch, protected by the observer monitor. */
	private final ByteBuilder word=new ByteBuilder(), row=new ByteBuilder();
	/** Queued encounters, including writer-buffered rows and repeat visits; I/O may still fail. */
	private long rows;
	/** Trailer queued and close attempted, respectively; neither substitutes for I/O success. */
	private boolean completed, closed;
	/** Exact 19-column trace schema, with provenance distinct from model features. */
	public static final String HEADER="word\tk\tcount\tleftA\tleftC\tleftG\tleftT\trightA\trightC\trightG\trightT\tleft_branch\tright_branch\tarrival_side\tarrival_base\tentry\thidden\tlevel\tcontrol_every";
}
