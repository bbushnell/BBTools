package ddl;

/** Reusable family sketches with exact member votes and configurable winner selection.
 * The default chooses nonempty modal keys, with smaller-key ties; explicit policies
 * can change ties or add a selection-only consensus bonus without changing member counts.
 * A builder owns mutable profile, member and selection scratch and is not thread-safe.
 * Use fork() for independent workers sharing the immutable configuration and key policy.
 * Sequence inputs are consumed during each call and are not retained by the builder.
 * <p>
 * Start each family with reset(), add its members, then write into caller-owned output.
 * Query sketch methods overwrite their output without adding family votes or word counts.
 * Valid-word totals count valid word starts before blacklist exclusion, not distinct keys.
 * All input residues are consumed in plain mode, including proteins longer than 200 aa.
 * Existing boolean constructors retain fixed stop-relative tail geometry.
 * Relative coordinates require an explicit position mode on both model and query builders.
 * @author Brian Bushnell, Fischl
 */
public final class FamilySketchBuilder{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Configuration is explicit; random input-order seeds are unrelated to sketchSeed. */
	public FamilySketchBuilder(final int width_, final int k_, final long sketchSeed_,
			final boolean positional_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_){
		this(width_, ProteinKmerGeometry.exact(k_), sketchSeed_, positional_, partitions_, overlap_, blacklist_);
	}

	/** Exact-k counterpart for readers that retain an explicit stored position mode. */
	public FamilySketchBuilder(final int width_, final int k_, final long sketchSeed_,
			final ProteinPositionMode positionMode_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_){
		this(width_,ProteinKmerGeometry.exact(k_),sketchSeed_,positionMode_,partitions_,overlap_,blacklist_);
	}

	/** Configuration is explicit; random input-order seeds are unrelated to sketchSeed. */
	public FamilySketchBuilder(final int width_, final ProteinKmerGeometry geometry_, final long sketchSeed_,
			final boolean positional_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_){
		this(width_,geometry_,sketchSeed_,positional_ ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,
			partitions_,overlap_,blacklist_);
	}

	/** Explicit coordinate identity is immutable and follows each worker fork. */
	public FamilySketchBuilder(final int width_, final ProteinKmerGeometry geometry_, final long sketchSeed_,
			final ProteinPositionMode positionMode_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_){
		this(width_,geometry_,sketchSeed_,positionMode_,partitions_,overlap_,blacklist_,null);
	}

	/** A nonnull immutable context selects packed finalization; null retains the default policy.
	 * The final constructor may also force packed selection for a positive consensus bonus. */
	public FamilySketchBuilder(final int width_, final ProteinKmerGeometry geometry_, final long sketchSeed_,
			final ProteinPositionMode positionMode_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_,
			final DDLKeyTieBreak.PolicyContext keyPolicy_){
		this(width_,geometry_,sketchSeed_,positionMode_,partitions_,overlap_,blacklist_,keyPolicy_,keyPolicy_!=null);
	}

	/**
	 * Resolves or validates the official blacklist policy, checks the kernel binding and
	 * allocates private member/profile scratch. An omitted key policy uses SMALLEST.
	 * Nondefault tie policies require packed selection; a positive consensus bonus forces it.
	 * Geometry, blacklist and explicit policy are shared with forks; mutable scratch is not.
	 */
	public FamilySketchBuilder(final int width_, final ProteinKmerGeometry geometry_, final long sketchSeed_,
			final ProteinPositionMode positionMode_, final int partitions_, final int overlap_, final DDLPentamerBlacklist blacklist_,
			final DDLKeyTieBreak.PolicyContext keyPolicy_,final boolean packedSelection_){
		require(width_>0 && width_<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width_)==1,
			"Sketch width must be a bounded power of two");
		require(geometry_!=null, "Family sketches require an explicit protein geometry");
		final DDLPentamerBlacklist selected=ProteinBlacklistPolicy.resolve(geometry_,blacklist_);
		require(positionMode_!=null, "Family/query sketches require an explicit position mode");
		if(positionMode_.positional()){
			DDLFamilyProfiles.validatePositionalWidth(width_);
			DDLFamilyProfiles.validatePositionalLayout(partitions_,overlap_,geometry_,positionMode_);
		}
		width=width_; geometry=geometry_; k=geometry.k(); sketchSeed=sketchSeed_; positionMode=positionMode_;
		partitions=partitions_; overlap=overlap_; blacklist=selected;
		member=new char[width]; profile=new DDLFamilyProfiles.Profile(width);
		final String pin=blacklist==null ? null : blacklist.sha80();
		keyPolicy=keyPolicy_==null ? new DDLKeyTieBreak.PolicyContext(DDLKeyTieBreak.Policy.SMALLEST,null,
			width,geometry,sketchSeed,positionMode,partitions,overlap,pin,null) : keyPolicy_;
		keyPolicy.requireGeometry(width,geometry,sketchSeed,positionMode,partitions,overlap,pin);
		// Positive consensus bonuses require packed candidate selection so
		// consensus-only keys can enter the modal vote. The tested W1..W10000
		// lookup-only arms all lost raw objective; see
		// results/family_sketch_phase_a_20261116.
		packedSelection=packedSelection_ || keyPolicy.consensusMemberBonus>0;
		require(packedSelection || keyPolicy.policy==DDLKeyTieBreak.Policy.SMALLEST,
			"Nondefault policies require packed selection");
		keyScratch=packedSelection ? new DDLKeyTieBreak.Workspace(width) : null;
		consensusKeys=keyPolicy.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY || keyPolicy.consensusMemberBonus>0 ?
			new DDLKeyTieBreak.KeySet(width) : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Shares configuration, resolved blacklist and key policy, with fresh worker-owned scratch.
	 * Existing family votes and word counters are not copied. Construction revalidates the binding.
	 */
	public FamilySketchBuilder fork(){
		return new FamilySketchBuilder(width, geometry, sketchSeed, positionMode, partitions, overlap, blacklist,keyPolicy,packedSelection);
	}

	/** Clears family votes, member count and valid-word total while retaining scratch capacity.
	 * Member/selection scratch is overwritten by the next sketch/finalization, not zeroed here. */
	public void reset(){
		profile.clear(); validWords=0;
		assert(profile.members()==0) : "New family votes must start from an empty modal profile";
	}

	/** Adds one standard-code member (0..19 or encoded X), with one vote per nonzero register.
	 * An all-zero sketch still increments the member count; valid words include blacklisted words.
	 * Encoded X breaks words; other nonstandard codes are rejected by the shared kernel. */
	public void addEncoded(final byte[] sequence){
		final long words=sketchEncoded(sequence,member);
		profile.add(member); validWords=Math.addExact(validWords,words);
		assert(profile.members()>0) : "A completed member sketch must contribute one family vote";
	}

	/** Adds one ASCII member; ambiguity/stops break words in the shared kernel.
	 * Member counting and pre-blacklist valid-word totals match addEncoded. */
	public void addRaw(final byte[] sequence){
		final long words=sketchRaw(sequence,member);
		profile.add(member); validWords=Math.addExact(validWords,words);
		assert(profile.members()>0) : "A completed raw member sketch must contribute one family vote";
	}

	/** Overwrites exact-width output with nonempty modal keys, choosing smaller keys on ties.
	 * Empty votes do not compete with nonzero keys. Requires at least one added member,
	 * SMALLEST policy and zero consensus bonus; output is not retained and votes are unchanged. */
	public void writeProfile(final char[] out){
		require(keyPolicy.policy==DDLKeyTieBreak.Policy.SMALLEST && keyPolicy.consensusMemberBonus==0,
			"Nondefault family selection requires explicit identity and consensus");
		require(out!=null && out.length==width && profile.members()>0, "A family profile needs members and an exact-width output");
		for(int bucket=0; bucket<width; bucket++){out[bucket]=profile.value(bucket,false);}
	}

	/**
	 * Writes the selected profile for this rebuild without changing accumulated member votes.
	 * Nonpacked selection delegates to the default overload and ignores label/consensus.
	 * Packed diversity policies require a nonempty familyLabel; consensus-diversity or
	 * a positive bonus also requires encodedConsensus, enumerated afresh on every call.
	 * Each argument is ignored when its corresponding condition does not apply.
	 * Output must have exact width and is not retained.
	 */
	public void writeProfile(final char[] out,final String familyLabel,final byte[] encodedConsensus){
		require(out!=null && out.length==width && profile.members()>0,"Require populated profile and exact-width output");
		if(!packedSelection){writeProfile(out);return;}
		final boolean hash=keyPolicy.policy==DDLKeyTieBreak.Policy.DIVERSITY || keyPolicy.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY;
		final long seed=hash ? DDLKeyTieBreak.familySeed(familyLabel) : 0;
		if(consensusKeys!=null){
			require(encodedConsensus!=null,"Consensus policy requires the current encoded consensus");
			// Adds all consensus keys to the selectable set before modal reduction.
			// Phase A measured this as lookup-only, not a family-rebuild result.
			// See results/family_sketch_phase_a_20261116.
			DDLFamilyProfiles.enumerateKeysEncoded(encodedConsensus,geometry,sketchSeed,positionMode,partitions,overlap,blacklist,consensusKeys);
		}
		profile.writePackedValues(keyPolicy,seed,consensusKeys,keyScratch,null,out);
	}

	/** Returns the immutable policy context shared with forks; exposes no builder scratch. */
	public DDLKeyTieBreak.PolicyContext policyContext(){return keyPolicy;}
	/** Returns the effective finalization strategy, including a bonus-forced packed path. */
	public boolean packedSelection(){return packedSelection;}

	/** Overwrites caller-owned output from standard encoded residues without updating family state.
	 * Returns valid word starts before blacklist exclusion, not distinct keys or overlap-tag copies.
	 * Output must have exact width; neither input nor output is retained. */
	public long sketchEncoded(final byte[] sequence, final char[] out){
		require(out!=null && out.length==width, "Encoded sketch output must have the configured width");
		return positionMode.positional() ? DDLFamilyProfiles.sketchPositionalEncoded(sequence,geometry,sketchSeed,partitions,overlap,out,blacklist,positionMode) :
			DDLFamilyProfiles.sketchEncoded(sequence,geometry,sketchSeed,out,blacklist);
	}

	/** ASCII counterpart of sketchEncoded; ambiguity/stops break words.
	 * Overwrites output and returns pre-blacklist valid word starts without adding a member. */
	public long sketchRaw(final byte[] sequence, final char[] out){
		require(out!=null && out.length==width, "Raw sketch output must have the configured width");
		return positionMode.positional() ? DDLFamilyProfiles.sketchPositional(sequence,geometry,sketchSeed,partitions,overlap,positionMode,out,blacklist) :
			DDLFamilyProfiles.sketch(sequence,geometry,sketchSeed,out,blacklist);
	}

	/** Number of members added since reset, including members with no retained register. */
	public int members(){return profile.members();}
	/** Sum of pre-blacklist valid word starts from added members since reset. */
	public long validWords(){return validWords;}
	/** Configured number of registers, including empty ones. */
	public int width(){return width;}
	/** Configured coordinate mode shared by family and query sketches. */
	public ProteinPositionMode positionMode(){return positionMode;}
	/** Fixed seeds and new query sketches must use the same hash geometry and blacklist. */
	public boolean matchesGeometry(int width_,int k_,long seed_,boolean positional_,
			int partitions_,int overlap_,String blacklistSha80){
		return matchesGeometry(width_,k_,seed_,positional_ ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,
			partitions_,overlap_,blacklistSha80);
	}

	/** A relative builder must not accept a fixed-coordinate seed merely because both use tags. */
	public boolean matchesGeometry(int width_,int k_,long seed_,ProteinPositionMode positionMode_,
			int partitions_,int overlap_,String blacklistSha80){
		return matchesGeometry(width_,ProteinKmerGeometry.exact(k_),seed_,positionMode_,partitions_,overlap_,blacklistSha80);
	}

	/** Compares kernel settings and the declared blacklist identity, including plain-mode layout fields.
	 * Does not compare selection policy, consensus bonus or raw-tail counts, and performs no file check. */
	public boolean matchesGeometry(int width_,ProteinKmerGeometry geometry_,long seed_,ProteinPositionMode positionMode_,
			int partitions_,int overlap_,String blacklistSha80){
		return width==width_ && geometry.equals(geometry_) && sketchSeed==seed_ && positionMode==positionMode_ &&
			partitions==partitions_ && overlap==overlap_ &&
			(blacklist==null ? blacklistSha80==null : blacklist.sha80().equals(blacklistSha80));
	}

	private static void require(final boolean ok, final String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Register count, cached geometry k and positional layout parameters. */
	private final int width, k, partitions, overlap;
	/** Immutable word projection and hash seed shared with query sketches. */
	private final ProteinKmerGeometry geometry;
	private final long sketchSeed;
	/** Coordinate identity and resolved shared blacklist. */
	private final ProteinPositionMode positionMode;
	private final DDLPentamerBlacklist blacklist;
	/** Reused one-member register buffer and this family's vote accumulator. */
	private final char[] member;
	private final DDLFamilyProfiles.Profile profile;
	/** Immutable selection configuration and effective strategy. */
	private final DDLKeyTieBreak.PolicyContext keyPolicy;
	private final boolean packedSelection;
	/** Private packed-sort scratch, or null for direct modal finalization. */
	private final DDLKeyTieBreak.Workspace keyScratch;
	/** Private all-consensus-key set when required by the policy/bonus; otherwise null. */
	private final DDLKeyTieBreak.KeySet consensusKeys;
	/** Running pre-blacklist word-start total for the current family's added members. */
	private long validWords;
}
