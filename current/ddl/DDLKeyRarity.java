package ddl;

import map.IntHashMap;

/** Immutable unweighted posting depths for one complete baseline family profile set.
 * Counts refer to stored (bucket,uint16 value) keys, not original protein words.
 * The caller must supply the complete baseline with the same population and sketch
 * geometry as the candidate profiles. No input rows are retained or changed.
 * @author Yaoyao
 */
public final class DDLKeyRarity{

	/*--------------------------------------------------------------*/
	/*----------------         Policy Names         ----------------*/
	/*--------------------------------------------------------------*/

	/** A construction policy, separate from query-time rarity scoring. */
	public enum Policy{
		SMALLEST("smallest"), RARITY("rarity");

		/** Binds the exact serialized/CLI spelling. */
		private Policy(String label_){
			label=label_;
		}
		/** Case-sensitive external spelling for this policy. */
		public final String label;

		/** Parses an exact label; null or an unknown value throws rather than selecting a default. */
		public static Policy parse(String value){
			for(Policy policy : values()){
				if(policy.label.equals(value)){return policy;}
			}
			throw new IllegalArgumentException("keytiebreak must be smallest or rarity: "+value);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds detached counts before any candidate row is selected. Zero is not a posting.
	 * Each row is one physical family; an all-zero row contributes no postings.
	 * This factory cannot validate external provenance; callers must bind that first.
	 * @param width Power-of-two sketch width in the supported family-profile range
	 * @param profiles Nonempty complete baseline; every row must have exactly width registers
	 * @return Immutable counts; input arrays may be reused after this method returns
	 */
	public static DDLKeyRarity fromProfiles(final int width, final char[][] profiles){
		if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1 ||
			profiles==null || profiles.length==0){
			throw new IllegalArgumentException("Rarity needs a nonempty complete baseline and a valid power-of-two width");
		}
		final IntHashMap counts=new IntHashMap(width*2);
		long postings=0;
		for(char[] row : profiles){
			if(row==null || row.length!=width){
				throw new IllegalArgumentException("Every baseline family must have the declared rarity width: "+width);
			}
			for(int bucket=0; bucket<width; bucket++){
				final char value=row[bucket];
				if(value!=0){
					final int depth=counts.increment((bucket<<16)|value);
					assert(depth>0 && depth<=profiles.length) : "One family contributes at most one posting per bucket/value";
					postings++;
				}
			}
		}
		return new DDLKeyRarity(width, profiles.length, postings, counts);
	}

	/** Takes exclusive ownership of the completed factory-local count map. */
	private DDLKeyRarity(int width_, int families_, long postings_, IntHashMap counts_){
		width=width_;
		families=families_;
		postings=postings_;
		counts=counts_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Frozen Counts         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the number of baseline rows carrying value in bucket.
	 * Absent nonzero keys have depth zero; zero is not a valid lookup key.
	 */
	public int depth(final int bucket, final char value){
		if(bucket<0 || bucket>=width || value==0){
			throw new IllegalArgumentException("Posting depth requires a valid bucket and a nonzero uint16 value");
		}
		final int count=counts.get((bucket<<16)|value);
		assert(count==-1 || count>0 && count<=families) : "Frozen posting counts must represent complete family rows";
		return count<0 ? 0 : count;
	}

	/** Chooses lower baseline depth, then the smaller unsigned key on an exact depth tie.
	 * Call only after equal member support is established; this method does not check support.
	 */
	public boolean prefers(final int bucket, final char candidate, final char incumbent){
		final int candidateDepth=depth(bucket, candidate), incumbentDepth=depth(bucket, incumbent);
		return candidateDepth<incumbentDepth || candidateDepth==incumbentDepth && candidate<incumbent;
	}

	/** Returns the number of register positions in each baseline row. */
	public int width(){return width;}
	/** Returns the number of baseline rows, including rows with no populated registers. */
	public int families(){return families;}
	/** Returns nonzero register occurrences across rows, not the number of distinct keys. */
	public long postings(){return postings;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Register width and total baseline row count. */
	private final int width, families;
	/** Total nonzero register occurrences counted during construction. */
	private final long postings;
	/** Owned frozen map from packed bucket/value keys to row counts; never exposed to callers. */
	private final IntHashMap counts;
}
