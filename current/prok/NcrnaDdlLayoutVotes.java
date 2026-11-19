package prok;

import map.IntHashMap;

/** Transient plurality voting for plain DDL profiles. Each member supplies at
 * most one greatest key per bin; equal compressed keys are the same candidate.
 * A consensus bonus applies once per distinct bin/key. Vote counts exist only
 * while building and never accompany the returned char[] profile.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlLayoutVotes{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlLayoutVotes(final int bins_){
		if(bins_<1 || bins_>32768){
			throw new IllegalArgumentException("Plain voting requires one to32768bins so bin15/key16 fits a nonnegative primitive map key");
		}
		bins=bins_;
		counts=new IntHashMap(bins);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Member Voting        ----------------*/
	/*--------------------------------------------------------------*/

	/** Begins another cluster while retaining the primitive map's capacity. */
	public void reset(){
		counts.clear();
		members=0;
		assert(counts.size()==0) : "A new cluster must never inherit the preceding cluster's member votes";
	}

	/** Consumes borrowed member scratch before the next sequence replaces it. */
	public void add(final char[] keys){
		if(keys==null || keys.length!=bins || members==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Each real voter needs one complete key array and a representable member count");
		}
		members++;
		for(int bin=0; bin<bins; bin++){
			final int key=keys[bin];
			if(key==0){continue;}
			final int support=counts.increment((bin<<16)|key);
			assert(support<=members)
				: "A real member supplies one vote per bin; raw-word multiplicity must not become additional voters";
		}
	}

	public Result finish(){return finish(null, 0, true);}

	/** Updates a current greedy profile after one real member, without rescanning
	 * the family's complete vote map. Its consensus prefix must be unchanged;
	 * consensus replacement requires finish() once before further incremental adds. */
	public boolean addAndSelect(final char[] keys, final KeyPrefix consensus, final int bonus, final char[] selected){
		if(consensus==null || consensus.buckets()!=bins || selected==null || selected.length!=bins || bonus<0){
			throw new IllegalArgumentException("Incremental profile updates require matching keys, current consensus prefix and nonnegative bonus");
		}
		add(keys);
		boolean changed=false;
		for(int bin=0; bin<bins; bin++){
			final int key=keys[bin];
			if(key==0){continue;}
			final int candidate=(bin<<16)|key, winner=(bin<<16)|selected[bin];
			final long support=counts.get(candidate)+(consensus.contains(candidate) ? (long)bonus : 0);
			final long old=Math.max(0, counts.get(winner))+(consensus.contains(winner) ? (long)bonus : 0);
			if(support>old || support==old && key>selected[bin]){
				selected[bin]=(char)key;
				changed=true;
			}
		}
		return changed;
	}

	/** Highest member support plus the declared bonus wins, then greatest key.
	 * Winner weights are temporary construction scratch, not profile metadata. */
	public Result finish(final NcrnaDdlLayoutKeys consensus, final int bonus){
		return finish(consensus, bonus, true);
	}

	/** Shared selector also retains the legacy smaller-key tie policy. The key
	 * prefix is borrowed read-only; weights never enter the real member counts. */
	public Result finish(final KeyPrefix consensus, final long bonus, final boolean greatestTie){
		if(members<1 || bonus<0 || bonus>Integer.MAX_VALUE+1L || bonus>0 && consensus==null
			|| consensus!=null && consensus.buckets()!=bins){
			throw new IllegalArgumentException("Voting needs real members, a nonnegative bonus and matching consensus bins");
		}
		final Result result=new Result(new char[bins], members);
		final long[] best=new long[bins];
		final int[] keys=counts.keys(), values=counts.values();
		final int invalid=counts.invalid();
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){select(result.keys, best, keys[i], values[i], consensus, bonus, greatestTie);}
		}
		if(bonus>0){
			for(int i=0; i<consensus.keyCount(); i++){
				final int packed=consensus.packedKey(i);
				select(result.keys, best, packed, Math.max(0, counts.get(packed)), consensus, bonus, greatestTie);
			}
		}
		return result;
	}

	private void select(final char[] winners, final long[] best, final int packed,
			final int support, final KeyPrefix consensus, final long bonus, final boolean greatestTie){
		final int bin=packed>>>16, key=packed&65535;
		assert(bin<bins && key>0 && support>=0 && support<=members)
			: "The transient vote map binds each candidate to a real bin, nonempty key and actual member support";
		final long votes=support+(bonus>0 && consensus.contains(packed) ? bonus : 0L);
		if(votes>best[bin] || votes==best[bin] && (greatestTie ? key>winners[bin] : key<winners[bin])){
			winners[bin]=(char)key;
			best[bin]=votes;
		}
	}

	public int buckets(){return bins;}
	public int members(){return members;}

	/** Read-only primitive consensus candidates used only during construction. */
	public interface KeyPrefix{
		int buckets();
		int keyCount();
		int packedKey(int index);
		boolean contains(int packed);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Result           ----------------*/
	/*--------------------------------------------------------------*/

	/** Construction result owns its keys; the scalar is a validation count only. */
	public static final class Result{
		private Result(final char[] keys_, final int members_){
			keys=keys_;
			members=members_;
		}
		public final char[] keys;
		public final int members;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int bins;
	private final IntHashMap counts;
	private int members;
}
