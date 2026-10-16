package assemble;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Observation-only trace of actual contig-walk junction checks. Counts remain
 * immutable and no return value is supplied to the walker. Repeated encounters
 * are retained explicitly; corpus construction must deduplicate canonical words.
 * Reference truth is deliberately unavailable in this process.
 * @author Fischl
 */
final class LiveBranchObserver implements AutoCloseable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a new diagnostic output only after the driver has checked all aliases. */
	LiveBranchObserver(final String path, final int k_, final int minExtend_, final int controlEvery_){
		if(path==null || k_<2 || minExtend_<1 || controlEvery_<1){
			throw new IllegalArgumentException("Trace needs a path, K>=2, mce>=1 and positive control interval.");
		}
		k=k_; minExtend=minExtend_; controlEvery=controlEvery_;
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
	 */
	boolean sample(final long canonicalKey, final boolean stopped){
		if(stopped || controlEvery==1){return true;}
		long hash=canonicalKey;
		hash=(hash^(hash>>>33))*0xff51afd7ed558ccdL;
		hash=(hash^(hash>>>33))*0xc4ceb9fe1a85ec53L;
		hash^=hash>>>33;
		return (hash&Long.MAX_VALUE)%controlEvery==0;
	}

	/** Decodes the oriented short key into reusable scratch, not its canonical table key. */
	synchronized void record(final long forward, final int count, final int[] left, final int[] right,
			final boolean lb, final boolean rb, final int arrival, final boolean seed, final boolean hidden){
		assert(k<=31) : "Packed short-key decoding is only valid for the Tadpole1 representation.";
		word.clear().appendKmer(forward, k);
		write(count, left, right, lb, rb, arrival, seed, hidden);
	}

	/** Decodes array1 explicitly: AbstractKmerTableU.toText uses canonical key(), not traversal orientation. */
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

	/** Writes exact table population after the successful walk; no second count table. */
	synchronized void complete(final SameKJoinFeatures.Moments moments){
		if(closed || completed || moments==null){throw new IllegalStateException("Trace completion must occur exactly once after workers finish.");}
		row.clear().append("#complete\t").append(rows).tab().append(moments.unique);
		row.tab().appendSlow(moments.observations).tab().appendSlow(moments.squares);
		row.tab().append(controlEvery).tab().append(k).tab().append(minExtend).tab().appendSlow(0.5).nl();
		writer.print(row);
		completed=true;
	}

	/** Shuts down the writer even when assembly fails; failed runs have no completion trailer. */
	@Override
	public synchronized void close(){
		if(closed){return;}
		closed=true;
		if(writer.poisonAndWait()){throw new IllegalStateException("Live branch trace output failed.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int k, minExtend, controlEvery;
	private final ByteStreamWriter writer;
	private final ByteBuilder word=new ByteBuilder(), row=new ByteBuilder();
	private long rows;
	private boolean completed, closed;
	static final String HEADER="word\tk\tcount\tleftA\tleftC\tleftG\tleftT\trightA\trightC\trightG\trightT\tleft_branch\tright_branch\tarrival_side\tarrival_base\tentry\thidden\tlevel\tcontrol_every";
}
