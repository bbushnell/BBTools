package prok;

import simd.Vector;

/** Plain member-voted DDL keys for one consensus model and sketch geometry.
 * One shared transient accumulator supplies plurality and consensus bonuses.
 * Legacy finish() uses smaller-key ties; the explicit weighted form uses greater
 * keys. Finished profiles retain only their keys and real member count, with
 * no support array or scoring metadata. Empty buckets never vote or match.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlProfile{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaDdlProfile(final NcrnaDdlLayoutVotes.Result result){
		assert(result!=null && result.members>0 && result.keys!=null)
			: "A finished profile must own keys selected from actual assigned members";
		members=result.members;
		keys=result.keys;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates four counters; use the scratch overload for repeated comparisons. */
	public int equalBuckets(final char[] query){return equalBuckets(query, new int[4]);}

	/** Borrows read-only query bins and caller-owned reusable comparison scratch. */
	public int equalBuckets(final char[] query, final int[] counts){
		if(query==null || query.length!=keys.length){
			throw new IllegalArgumentException("Query and member-voted profile must use the same bucket count");
		}
		final int matches=Vector.compareDDL(query, keys, counts)[1];
		assert(matches<=keys.length) : "A profile comparison counts each bucket at most once";
		return matches;
	}

	/** Borrowed read-only keys for the native one-key sketch comparator. */
	public char[] keys(){return keys;}
	public int buckets(){return keys.length;}

	/*--------------------------------------------------------------*/
	/*----------------        Member Voting         ----------------*/
	/*--------------------------------------------------------------*/

	/** Compatibility geometry/tie policy over the shared transient vote engine.
	 * Sources may be reused immediately; every finish owns independent plain keys. */
	public static final class Votes{

		public Votes(final int buckets){
			if(buckets<1 || buckets>32768 || Integer.bitCount(buckets)!=1){
				throw new IllegalArgumentException("Profile width must be a power of two in [1,32768]; packed votes require15bin bits and16key bits");
			}
			votes=new NcrnaDdlLayoutVotes(buckets);
		}

		public void add(final char[] sketch){votes.add(sketch);}

		/** Legacy plurality chooses the smaller nonempty key on equal support. */
		public NcrnaDdlProfile finish(){
			checkMembers();
			return new NcrnaDdlProfile(votes.finish(null, 0, false));
		}

		/** Every distinct raw consensus candidate gets its bonus once; greater
		 * keys break ties. CONSENSUS_WINS uses real members+1 without overflow. */
		public NcrnaDdlProfile finish(final NcrnaDdlQueryKeys consensus, final int weight){
			checkMembers();
			if(weight<CONSENSUS_WINS || weight!=0 && consensus==null
				|| consensus!=null && consensus.buckets()!=votes.buckets()){
				throw new IllegalArgumentException("Consensus voting needs matching raw-key width and nonnegative weight, or CONSENSUS_WINS");
			}
			final long bonus=weight==CONSENSUS_WINS ? (long)votes.members()+1 : weight;
			return new NcrnaDdlProfile(votes.finish(consensus, bonus, true));
		}

		/** Member plurality occupies the first three quarters; consensus maxima
		 * occupy only the last quarter and never add real construction members. */
		public NcrnaDdlProfile finishPartition(final char[] consensus){
			final int buckets=votes.buckets();
			if(consensus==null || consensus.length!=buckets || buckets<4){
				throw new IllegalArgumentException("Partition consensus must match the complete four-quarter profile width");
			}
			final NcrnaDdlProfile profile=finish(null, 0);
			final int memberBins=3*(buckets/4);
			for(int bucket=0; bucket<buckets; bucket++){
				if(bucket<memberBins){
					if(consensus[bucket]!=0){throw new IllegalArgumentException("Consensus keys cannot enter a member partition");}
				}else{
					if(profile.keys[bucket]!=0){throw new IllegalArgumentException("Member votes cannot enter the consensus quarter");}
					profile.keys[bucket]=consensus[bucket];
				}
			}
			return profile;
		}

		private void checkMembers(){
			if(votes.members()<1){throw new IllegalStateException("Cannot build a consensus-model profile without assigned members");}
		}

		private final NcrnaDdlLayoutVotes votes;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final int members;
	/** Bonus greater than all actual member votes, evaluated as a long. */
	public static final int CONSENSUS_WINS=-1;
	private final char[] keys;
}
