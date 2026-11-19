package ddl;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import map.IntHashMap;

/** Fixed 40-bit selection scores and 24-bit construction tie signals.
 * Selection scores normally equal member support; an explicit consensus weight
 * can increase them without changing the underlying support or occupancy.
 * Raw rarity counts held input-tail registers, never family posting lists.
 * Immutable policy/count objects may be shared. Diagnostics, key sets, counters
 * and workspaces are mutable and require one owner or external serialization.
 * @author Yaoyao, Fischl
 */
public final class DDLKeyTieBreak{

	private DDLKeyTieBreak(){}

	/** Low-order bits reserved for the secondary tie signal. */
	public static final int SIGNAL_BITS=24;
	/** Largest representable secondary signal. */
	public static final int MAX_SIGNAL=0xffffff;
	/** Largest primary selection score, including any explicit member bonus. */
	public static final long MAX_VOTES=(1L<<40)-1;
	/** Low 23 hash bits leave one bit for consensus-preferred ties. */
	public static final int HASH_MASK=0x7fffff;
	/** Secondary-signal bit; distinct from the optional primary member bonus. */
	public static final int CONSENSUS_BONUS=1<<23;
	/** Stable recipe identity for family-specific diversity ordering. */
	public static final String HASH_RECIPE="fnv1a64-utf8-family-xor-uint16-splitmix64-finalizer-low23-v1";

	/** Secondary ordering among candidates with equal primary selection scores. */
	public enum Policy{
		SMALLEST("smallest"), RAW_TAIL_RARITY("rawtail"), DIVERSITY("diversity"), CONSENSUS_DIVERSITY("consensus");
		private Policy(String label_){label=label_;}
		/** Exact case-sensitive command-line and stored policy name. */
		public final String label;
		/** Parses one exact policy label; null and unknown labels are rejected. */
		public static Policy parse(String label){
			for(Policy value:values()){if(value.label.equals(label)){return value;}}
			throw new IllegalArgumentException("keytiebreak must be smallest, rawtail, diversity or consensus: "+label);
		}
	}

	/** Wire identity is centralized for builders, codecs and preparation receipts. */
	public static String encodingVersion(Policy policy){
		if(policy==null){throw new IllegalArgumentException("Missing key policy");}
		switch(policy){
			case SMALLEST: return "member40-key24-smallest-v1";
			case RAW_TAIL_RARITY: return "member40-distinct-count-rank24-v1";
			case DIVERSITY: return "member40-fnv1a-splitmix23-v1";
			case CONSENSUS_DIVERSITY: return "member40-consensus1-fnv1a-splitmix23-v1";
			default: throw new AssertionError(policy);
		}
	}

	/** Immutable full-view policy/kernel binding; never owns worker scratch. */
	public static final class PolicyContext{
		/** Binds a policy with no primary consensus member bonus. */
		public PolicyContext(Policy policy_,RawTailCounts rawCounts_,int width_,ProteinKmerGeometry geometry_,
				long sketchSeed_,ProteinPositionMode positionMode_,int partitions_,int overlap_,
				String blacklistSha80_,String rawInputBinding_){
			this(policy_,rawCounts_,width_,geometry_,sketchSeed_,positionMode_,partitions_,overlap_,
				blacklistSha80_,rawInputBinding_,0);
		}

		/** A family-only bonus is applied once per distinct consensus bin/key, not per occurrence.
		 * Phase A tested W1/W10/W1000/W10000 as fixed-family lookup arms; all
		 * reduced ALL/40-99 raw objective and stayed off by default. See
		 * results/family_sketch_phase_a_20261116.
		 * The old constructor remains the byte-compatible W0 recipe.
		 */
		public PolicyContext(Policy policy_,RawTailCounts rawCounts_,int width_,ProteinKmerGeometry geometry_,
				long sketchSeed_,ProteinPositionMode positionMode_,int partitions_,int overlap_,
				String blacklistSha80_,String rawInputBinding_,long consensusMemberBonus_){
			if(consensusMemberBonus_<0 || consensusMemberBonus_>MAX_VOTES){
				throw new IllegalArgumentException("Consensus member bonus must fit nonnegative40-bit selection scores");
			}
			if(policy_==null || geometry_==null || positionMode_==null || width_<1 ||
				width_>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width_)!=1 || partitions_<0 || overlap_<0){
				throw new IllegalArgumentException("Policy context requires a complete valid kernel binding");
			}
			if(positionMode_.positional()){
				DDLFamilyProfiles.validatePositionalWidth(width_);
				DDLFamilyProfiles.validatePositionalLayout(partitions_,overlap_,geometry_,positionMode_);
			}
			if(policy_==Policy.RAW_TAIL_RARITY){
				if(rawCounts_==null || rawCounts_.width()!=width_ || rawInputBinding_==null || rawInputBinding_.trim().isEmpty()){
					throw new IllegalArgumentException("Rawtail needs frozen full-view counts and a canonical input binding");
				}
			}else if(rawCounts_!=null || rawInputBinding_!=null){
				throw new IllegalArgumentException("Only rawtail consumes a raw count/input binding");
			}
			policy=policy_; rawCounts=rawCounts_; width=width_; geometry=geometry_; sketchSeed=sketchSeed_;
			positionMode=positionMode_; partitions=partitions_; overlap=overlap_; blacklistSha80=blacklistSha80_;
			rawInputBinding=rawInputBinding_; encodingVersion=DDLKeyTieBreak.encodingVersion(policy);
			consensusMemberBonus=consensusMemberBonus_;
		}
		/** Requires the same width, encoding, hash seed, position layout and blacklist binding.
		 * Checks declared identity; it does not recount raw inputs or inspect a blacklist file.
		 */
		public void requireGeometry(int width_,ProteinKmerGeometry geometry_,long sketchSeed_,
				ProteinPositionMode positionMode_,int partitions_,int overlap_,String blacklistSha80_){
			if(width!=width_ || !geometry.equals(geometry_) || sketchSeed!=sketchSeed_ || positionMode!=positionMode_ ||
				partitions!=partitions_ || overlap!=overlap_ || !java.util.Objects.equals(blacklistSha80,blacklistSha80_)){
				throw new IllegalArgumentException("Key policy context does not match the family/query kernel");
			}
		}
		/** Secondary selection policy shared by worker builders. */
		public final Policy policy;
		/** Immutable full-input counts for raw-tail rarity; null for other policies. */
		public final RawTailCounts rawCounts;
		/** Policy encoding identity and optional canonical raw-input identity. */
		public final String encodingVersion,rawInputBinding;
		/** Selection-only weight; actual member support/occupancy are never increased. */
		public final long consensusMemberBonus;
		final int width;
		private final int partitions,overlap;
		private final long sketchSeed;
		private final ProteinKmerGeometry geometry;
		private final ProteinPositionMode positionMode;
		private final String blacklistSha80;
	}

	/** Packs a primary score in the high 40 bits and a tie signal in the low 24.
	 * The votes argument may include an explicit selection-only consensus bonus.
	 * Callers must not compare the result as an ordinary signed long.
	 */
	public static long pack(long votes,int signal){
		if(votes<0 || votes>MAX_VOTES || signal<0 || signal>MAX_SIGNAL){
			throw new IllegalArgumentException("Pivot requires unsigned40-bit votes and24-bit signal");
		}
		return (votes<<SIGNAL_BITS)|(long)signal;
	}

	/** Monotone signed representation, useful for independent checks and diagnostics. */
	public static long signedOrder(long pivot){return pivot^Long.MIN_VALUE;}
	/** Extracts the unsigned primary score; it is not necessarily actual support. */
	public static long votes(long pivot){return pivot>>>SIGNAL_BITS;}
	/** Extracts the low 24-bit secondary ordering signal. */
	public static int signal(long pivot){return (int)(pivot&MAX_SIGNAL);}

	/** Family labels are exact retained IDs, including shard and split suffixes.
	 * Deliberate modulo2^64 arithmetic belongs only to this named hash recipe.
	 */
	public static long familySeed(String familyId){
		if(familyId==null || familyId.isEmpty()){throw new IllegalArgumentException("Require a retained family ID");}
		long hash=0xcbf29ce484222325L;
		for(byte value:familyId.getBytes(StandardCharsets.UTF_8)){
			hash^=value&255L;
			hash*=0x100000001b3L;
		}
		return hash;
	}

	/** Diversity and consensus-diversity use the same 23-bit hash; collisions remain possible. */
	public static int diversity(long familySeed,int key){
		checkKey(key);
		long hash=familySeed^(long)key;
		hash=(hash^(hash>>>30))*0xbf58476d1ce4e5b9L;
		hash=(hash^(hash>>>27))*0x94d049bb133111ebL;
		hash^=hash>>>31;
		return (int)(hash&HASH_MASK);
	}

	/** Returns a larger-is-preferred secondary signal for one nonzero key.
	 * Consensus membership precedes per-bin max reduction. This raw-tail branch
	 * uses a direct count limited to 24 bits; live PolicyContext selection instead
	 * calls RawTailCounts.raritySignal for an exact distinct-count rank.
	 */
	public static int tieSignal(Policy policy,int key,int rawTailCount,long familySeed,boolean inConsensus){
		if(policy==null){throw new IllegalArgumentException("An explicit pivot policy is required");}
		checkKey(key);
		switch(policy){
			case SMALLEST: return 65535-key;
			case RAW_TAIL_RARITY:
				if(rawTailCount<0 || rawTailCount>MAX_SIGNAL){
					throw new IllegalArgumentException("Raw-tail key count exceeds24 bits: "+rawTailCount);
				}
				return MAX_SIGNAL-rawTailCount;
			case DIVERSITY: return diversity(familySeed,key);
			case CONSENSUS_DIVERSITY: return diversity(familySeed,key)|(inConsensus ? CONSENSUS_BONUS : 0);
			default: throw new AssertionError(policy);
		}
	}

	/** One stable LSD radix sort, descending unsigned64-bit pivot, with paired keys.
	 * Input keys must be unique and ascending, supplied by canonical bitset traversal.
	 * Equal complete pivots keep that order; no secondary comparator or tie scan.
	 * All scratch is caller-owned and reusable across family buckets.
	 */
	public static void sortDescending(long[] pivots,char[] keys,int size,long[] pivotScratch,
			char[] keyScratch,int[] bins){
		if(pivots==null || keys==null || pivotScratch==null || keyScratch==null || bins==null ||
			pivots==pivotScratch || keys==keyScratch || size<0 || size>pivots.length || size>keys.length ||
			size>pivotScratch.length || size>keyScratch.length || bins.length!=256){
			throw new IllegalArgumentException("Stable pivot sort requires distinct adequate scratch and256 bins");
		}
		for(int i=0; i<size; i++){
			if(keys[i]==0 || i>0 && keys[i-1]>=keys[i]){
				throw new IllegalArgumentException("Candidate keys must be nonzero, unique and ascending before sorting");
			}
		}
		if(size<2){return;}
		long[] source=pivots, target=pivotScratch;
		char[] sourceKeys=keys, targetKeys=keyScratch;
		for(int shift=0; shift<64; shift+=8){
			Arrays.fill(bins,0);
			for(int i=0; i<size; i++){bins[255-(int)((source[i]>>>shift)&255)]++;}
			int offset=0;
			for(int digit=0; digit<256; digit++){
				final int count=bins[digit]; bins[digit]=offset; offset+=count;
			}
			for(int i=0; i<size; i++){
				final int digit=255-(int)((source[i]>>>shift)&255), destination=bins[digit]++;
				target[destination]=source[i]; targetKeys[destination]=sourceKeys[i];
			}
			final long[] oldSource=source; source=target; target=oldSource;
			final char[] oldKeys=sourceKeys; sourceKeys=targetKeys; targetKeys=oldKeys;
		}
		assert(source==pivots && sourceKeys==keys) : "Eight stable passes return output to caller arrays";
	}

	/** Single-owner counters; never changes sorted arrays or the selected key.
	 * Collisions are equal complete pivots among maximal-primary-score candidates.
	 * With consensus-diversity this requires equal consensus bits and 23-bit hashes.
	 * A primary member bonus can make that score differ from actual member support.
	 */
	public static final class Diagnostics{
		/** Column order returned by summary(). */
		public static final String HEADER="occupied_bins\ttied_bins\tmaximal_candidates\tcollision_bins\tcollision_pairs\twinning_collision_bins\tdistinct_tied_winner_keys";
		/** Adds one bin from arrays already sorted by descending unsigned pivot. */
		public void observe(long[] pivots,char[] keys,int size){
			if(pivots==null || keys==null || size<0 || size>pivots.length || size>keys.length){
				throw new IllegalArgumentException("Diagnostics require valid sorted candidate arrays");
			}
			if(size==0){return;}
			occupiedBins++;
			final long maximum=votes(pivots[0]);
			int tied=1;
			while(tied<size && votes(pivots[tied])==maximum){tied++;}
			maximalCandidates+=tied;
			if(tied<2){return;}
			tiedBins++;
			checkKey(keys[0]);
			tiedWinners[keys[0]]++;
			long pairs=0;
			int first=0;
			while(first<tied){
				int end=first+1;
				while(end<tied && pivots[end]==pivots[first]){end++;}
				final long count=end-first;
				pairs+=count*(count-1)/2;
				assert(end==tied || Long.compareUnsigned(pivots[end-1],pivots[end])>0) : "Diagnostics require descending pivots";
				first=end;
			}
			if(pairs>0){collisionBins++;collisionPairs+=pairs;}
			if(pivots[0]==pivots[1]){winningCollisionBins++;}
		}
		/** Adds another owner's quiescent counters without clearing either object. */
		public void add(Diagnostics other){
			if(other==null || other==this){throw new IllegalArgumentException("Merge distinct worker diagnostics");}
			occupiedBins=Math.addExact(occupiedBins,other.occupiedBins);
			tiedBins=Math.addExact(tiedBins,other.tiedBins);
			maximalCandidates=Math.addExact(maximalCandidates,other.maximalCandidates);
			collisionBins=Math.addExact(collisionBins,other.collisionBins);
			collisionPairs=Math.addExact(collisionPairs,other.collisionPairs);
			winningCollisionBins=Math.addExact(winningCollisionBins,other.winningCollisionBins);
			for(int key=1;key<tiedWinners.length;key++){tiedWinners[key]=Math.addExact(tiedWinners[key],other.tiedWinners[key]);}
		}
		/** Bins won by this key where at least two candidates shared the top primary score. */
		public long tiedWins(int key){checkKey(key);return tiedWinners[key];}
		/** Tab-separated cumulative counts in HEADER order; does not reset this object. */
		public String summary(){
			int distinct=0;
			for(long count:tiedWinners){if(count>0){distinct++;}}
			return occupiedBins+"\t"+tiedBins+"\t"+maximalCandidates+"\t"+collisionBins+"\t"+collisionPairs+"\t"+winningCollisionBins+"\t"+distinct;
		}
		private long occupiedBins,tiedBins,maximalCandidates,collisionBins,collisionPairs,winningCollisionBins;
		/** Winning key histogram for bins tied on primary score; zero is unused. */
		private final long[] tiedWinners=new long[65536];
	}

	/** Single-owner reusable membership before per-bin maximum reduction.
	 * Stores each (bucket, nonzero key) pair once; different buckets remain distinct.
	 */
	public static final class KeySet{
		/** Creates an empty set for a bounded power-of-two sketch width. */
		public KeySet(int width_){
			if(width_<1 || width_>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width_)!=1){
				throw new IllegalArgumentException("All-key membership requires a valid profile width");
			}
			width=width_;
			keys=new IntHashMap(width*2);
		}
		/** Removes all pairs while retaining map capacity. */
		public void clear(){keys.clear();}
		/** Adds a validated pair idempotently; no occurrence count is retained. */
		public void add(int bucket,int key){
			checkKey(key);
			checkBucket(bucket);
			keys.put((bucket<<16)|key,1);
		}
		/** Tests one pair after validating both coordinates. */
		public boolean contains(int bucket,int key){
			checkKey(key);
			checkBucket(bucket);
			return keys.contains((bucket<<16)|key);
		}
		/** Configured number of buckets, including empty buckets. */
		public int width(){return width;}
		/** Number of distinct bucket/key pairs, not residues or occupied buckets. */
		public int size(){return keys.size();}
		/** Overwrites caller-owned output with per-bin maxima, using zero for empty bins.
		 * This permits direct parity checks against the ordinary sketch.
		 */
		public void maxima(char[] output){
			if(output==null || output.length!=width){throw new IllegalArgumentException("Membership/maxima width differs");}
			Arrays.fill(output,(char)0);
			for(int packed:keys.keys()){
				if(packed==keys.invalid()){continue;}
				final int bucket=packed>>>16,key=packed&65535;
				if(key>output[bucket]){output[bucket]=(char)key;}
			}
		}
		/** Requires an index inside this set's declared bucket domain. */
		private void checkBucket(int bucket){
			if(bucket<0 || bucket>=width){throw new IllegalArgumentException("All-key bucket outside width");}
		}
		private final int width;
		private final IntHashMap keys;
	}

	/** Thread-owned candidate enumeration and sort scratch; holds keys, not member sequences.
	 * Collection clears the previous family; capacity growth discards old sort contents.
	 */
	public static final class Workspace{
		/** Creates empty per-bin sets and lazily grown candidate/sort arrays. */
		public Workspace(int width){
			if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1){
				throw new IllegalArgumentException("Candidate workspace requires valid width");
			}
			candidates=new BitSet[width];
			for(int i=0;i<width;i++){candidates[i]=new BitSet();}
		}
		/** Replaces candidate membership from a trusted packed bucket/key count map. */
		void collect(IntHashMap counts){
			for(BitSet bucket:candidates){bucket.clear();}
			for(int packed:counts.keys()){
				if(packed==counts.invalid()){continue;}
				candidates[packed>>>16].set(packed&65535);
			}
		}
		/** Positive bonuses admit every consensus key, even if no member retained it as its bin winner.
		 * This made W1000/W10000 identical in the Phase A fixed-family lookup
		 * study; both lost -344/-218 ALL/40-99 raw objective and were not adopted.
		 */
		void collect(IntHashMap counts,KeySet consensus){
			if(consensus==null || consensus.width()!=candidates.length){
				throw new IllegalArgumentException("Consensus candidate union requires matching bin domains");
			}
			collect(counts);
			for(int packed:consensus.keys.keys()){
				if(packed!=consensus.keys.invalid()){candidates[packed>>>16].set(packed&65535);}
			}
		}
		/** Ensures each candidate array can hold size entries; old contents are scratch. */
		void capacity(int size){
			if(keys.length>=size){return;}
			keys=new char[size];
			keyScratch=new char[size];
			pivots=new long[size];
			pivotScratch=new long[size];
		}
		/** Canonical ascending enumeration source for each bucket. */
		final BitSet[] candidates;
		/** Paired key buffers reused by the stable radix sort. */
		char[] keys=new char[0],keyScratch=new char[0];
		/** Paired unsigned pivot buffers, with the same active length as keys. */
		long[] pivots=new long[0],pivotScratch=new long[0];
		/** Reused radix histogram and destination offsets. */
		final int[] bins=new int[256];
	}

	/** Rejects zero (the empty-register sentinel) and values outside uint16. */
	private static void checkKey(int key){
		if(key<1 || key>65535){throw new IllegalArgumentException("Require a nonzero uint16 stored key: "+key);}
	}

	/** A frozen complete original-input count table, shared by every comparison arm.
	 * Construction takes ownership of the counter's storage without copying it;
	 * callers never receive mutable arrays. Each bucket ranks all 65,535 nonzero
	 * keys by distinct observed count, including absent keys with count zero.
	 */
	public static final class RawTailCounts{
		/** Takes the completed count array and derives exact per-bucket count ranks. */
		private RawTailCounts(int width_,long tails_,long postings_,int[] counts_){
			width=width_; tails=tails_; postings=postings_; counts=counts_;
			ranks=new char[counts.length];
			final int[] sorted=new int[65535];
			int maximum=0;
			for(int bucket=0;bucket<width;bucket++){
				final int base=bucket<<16;
				System.arraycopy(counts,base+1,sorted,0,sorted.length);
				Arrays.sort(sorted);
				int distinct=0;
				for(int value:sorted){if(distinct==0 || sorted[distinct-1]!=value){sorted[distinct++]=value;}}
				maximum=Math.max(maximum,sorted[distinct-1]);
				for(int key=1;key<=65535;key++){
					final int rank=Arrays.binarySearch(sorted,0,distinct,counts[base|key]);
					assert(rank>=0 && rank<65535) : "Every exact count has one shared distinct-value rank";
					ranks[base|key]=(char)rank;
				}
			}
			maximumCount=maximum;
		}
		/** Number of input-tail sketches retaining this nonzero key in this bucket. */
		public int count(int bucket,int key){
			checkKey(key);
			if(bucket<0 || bucket>=width){throw new IllegalArgumentException("Raw-tail bucket outside declared width");}
			return counts[(bucket<<16)|key];
		}
		/** Larger-is-rarer signal from the bucket's distinct-count rank, not a raw count. */
		public int raritySignal(int bucket,int key){
			count(bucket,key);// Validate the same nonzero-key domain.
			return MAX_SIGNAL-ranks[(bucket<<16)|key];
		}
		/** Largest single bucket/key count across the complete input view. */
		public int maximumCount(){return maximumCount;}
		/** Number of registers in each counted tail sketch. */
		public int width(){return width;}
		/** Number of sketches added, including sketches with no occupied register. */
		public long tails(){return tails;}
		/** Total nonzero registers observed across all added sketches. */
		public long postings(){return postings;}
		private final int width;
		private final long tails,postings;
		/** Owned flattened bucket/key counts; key zero is unused in every bucket. */
		private final int[] counts;
		/** Zero-based distinct-count rank at each nonzero key. */
		private final char[] ranks;
		private final int maximumCount;
	}

	/** Single-owner streaming accumulator. Add one already-reduced DDL per raw tail.
	 * At512 bins counts use128MiB; freezing retains another64MiB for exact ranks.
	 * No serialized intermediate.
	 */
	public static final class RawTailCounter{
		/** Requires an exact positive input count; freeze rejects an incomplete view. */
		public RawTailCounter(int width_,long expectedTails_){
			if(width_<1 || width_>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width_)!=1 ||
				expectedTails_<1){
				throw new IllegalArgumentException("Require valid width and a positive complete input population");
			}
			width=width_; expectedTails=expectedTails_;
			counts=new int[Math.multiplyExact(width,65536)];
		}
		/** Consumes one sketch synchronously without retaining its array; zeros do not vote. */
		public void add(char[] sketch){
			if(counts==null || sketch==null || sketch.length!=width || tails>=expectedTails){
				throw new IllegalArgumentException("Raw-tail input must match width and exact expected population before freeze");
			}
			for(int bucket=0; bucket<width; bucket++){
				final int key=sketch[bucket];
				if(key==0){continue;}
				final int index=(bucket<<16)|key;
				if(counts[index]==Integer.MAX_VALUE){throw new IllegalArgumentException("Raw-tail count overflow");}
				counts[index]++; postings=Math.addExact(postings,1L);
			}
			tails++;
		}
		/** Transfers the completed counts to one immutable result and disables this counter.
		 * A second freeze or any subsequent add is rejected after successful transfer.
		 */
		public RawTailCounts freeze(){
			if(counts==null || tails!=expectedTails){
				throw new IllegalArgumentException("Freeze requires the complete expected raw input population, exactly once");
			}
			final RawTailCounts result=new RawTailCounts(width,tails,postings,counts);
			counts=null;
			return result;
		}
		private final int width;
		private final long expectedTails;
		private long tails,postings;
		/** Mutable owned counts until successful freeze; null after ownership transfer. */
		private int[] counts;
	}
}
