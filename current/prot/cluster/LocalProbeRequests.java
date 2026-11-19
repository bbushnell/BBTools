package prot.cluster;

import java.util.Arrays;
import java.util.Random;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;
import structures.LongList;

/** Reproducible uniform or opt-in size-rank-stratified local-probe requests.
 * Size stratification does not stratify on cohesion or outcomes. Sampling uses
 * current frozen memberships, cheap structural eligibility and independent
 * operation-specific random streams. It never constructs a proposal, evaluates
 * its objective or replaces a request that later proves infeasible.
 * @author Brian Bushnell, Fischl
 */
public final class LocalProbeRequests{

	private LocalProbeRequests(){}

	/*--------------------------------------------------------------*/
	/*----------------           Sampling           ----------------*/
	/*--------------------------------------------------------------*/

	/** Samples stable parent IDs and reciprocal nearest pairs without replacement.
	 * The table must have been built from this SAME frozen baseline and hash
	 * configuration. ID/count/width checks cannot certify unchanged model bytes
	 * or blacklist/seed identity; the driver retains those input bindings.
	 * The fusion table must apply the configured current-size ratio before its
	 * shortlist. Merge cutoffs are checked in BOTH cached directions, after
	 * nearest selection, without falling back to a second-best neighbor.
	 * A count may be zero, but their sum must be positive. Asking for more than
	 * either eligible population fails before any output is opened; no clamping.
	 */
	public static Selection sample(final LocalProbeBaseline baseline, final FamilyNeighborTable table,
			final LayerCakeClustering.Options options, final TailAlignmentGate mergeGate,
			final int fissionCount, final int fusionCount, final long seed){
		return sample(baseline, table, options, mergeGate, fissionCount, fusionCount, seed, false);
	}

	/** Optional quintiles rank parents by current size and pairs by combined size.
	 * Canonical stable IDs break ties, even when a size tie crosses a band boundary.
	 * Population and requested-count remainders go to the earliest bands. Each band
	 * uses its own deterministic stream; an unavailable quota fails without backfill.
	 */
	public static Selection sample(final LocalProbeBaseline baseline, final FamilyNeighborTable table,
			final LayerCakeClustering.Options options, final TailAlignmentGate mergeGate,
			final int fissionCount, final int fusionCount, final long seed, final boolean stratified){
		require(fissionCount>=0 && fusionCount>=0 && (long)fissionCount+fusionCount>0 &&
			(long)fissionCount+fusionCount<=Integer.MAX_VALUE-8,
			"Probe counts must be nonnegative with a positive bounded combined count");
		validate(baseline, table, options, mergeGate);
		final LongList fissions=new LongList(), fusions=new LongList();
		final long floor=Math.max((long)options.fissionMinSize, 2L*options.minChildSize);
		int reciprocal=0;
		// Table rows are already stable-ID ordered. One partner per row makes
		// canonical left<right pair order independent of caller roster ordering.
		for(int row=0; row<table.size(); row++){
			final int id=table.familyId(row);
			if(table.currentMembers(row)>=floor){fissions.add(pair(id, -1));}
			final int edge=table.nearestEdge(row);
			if(edge<0){continue;}
			final int other=table.targetRow(edge), partner=table.familyId(other);
			if(id>=partner){continue;}
			final int reverse=table.nearestEdge(other);
			if(reverse<0 || table.targetRow(reverse)!=row){continue;}
			assert(table.currentMembers(row)>0 && table.currentMembers(other)>0 &&
				table.reason(edge)==TailAlignmentDecision.Reason.ACCEPT &&
				table.reason(reverse)==TailAlignmentDecision.Reason.ACCEPT) :
				"A filtered reciprocal nearest pair requires two positive-count, geometry-valid endpoints";
			reciprocal++;
			if(passes(table, edge, mergeGate) && passes(table, reverse, mergeGate)){fusions.add(pair(id, partner));}
		}
		if(fissionCount>fissions.size || fusionCount>fusions.size){
			throw new IllegalArgumentException("Requested fission="+fissionCount+", fusion="+fusionCount+
				" but eligible fission="+fissions.size+", fusion="+fusions.size+
				" (reciprocal fusion before cheap gates="+reciprocal+"); counts are not clamped or backfilled");
		}
		final int[] left=new int[fissionCount+fusionCount], right=new int[left.length];
		final Bands splitBands, mergeBands;
		if(stratified){
			splitBands=drawStratified(baseline, fissions, fissionCount, seed^FISSION_SALT, left, right, 0);
			mergeBands=drawStratified(baseline, fusions, fusionCount, seed^FUSION_SALT, left, right, fissionCount);
		}else{
			splitBands=mergeBands=null;
			draw(fissions, fissionCount, new Random(seed^FISSION_SALT), left, right, 0);
			draw(fusions, fusionCount, new Random(seed^FUSION_SALT), left, right, fissionCount);
		}
		return new Selection(left, right, seed, fissionCount, fusionCount, fissions.size,
			fusions.size, reciprocal, floor, options, mergeGate, splitBands, mergeBands);
	}

	/** Sorts primitive size/rank keys, then samples independently inside five fixed rank bands. */
	private static Bands drawStratified(final LocalProbeBaseline baseline, final LongList pool,
			final int count, final long operationSeed, final int[] left, final int[] right, final int offset){
		final long[] order=new long[pool.size];
		for(int i=0; i<pool.size; i++){
			final long size=candidateSize(baseline, pool.array[i]);
			assert(size>0 && size<=0xfffffffeL) : "At most two positive signed-int memberships fit unsigned 32-bit size";
			// The pool is already in canonical ID order. Flipping the high sign bit
			// lets signed Arrays.sort order the unsigned size word before that rank.
			order[i]=((size<<32)|(i&0xffffffffL))^Long.MIN_VALUE;
		}
		Arrays.sort(order);
		final Bands bands=new Bands();
		for(int band=0, start=0; band<STRATA; band++){
			final int population=pool.size/STRATA+(band<pool.size%STRATA ? 1 : 0);
			final int quota=count/STRATA+(band<count%STRATA ? 1 : 0);
			require(quota<=population, "Size stratum "+band+" has "+population+" candidates for quota "+quota+"; no backfill");
			bands.population[band]=population; bands.quota[band]=quota;
			bands.minSize[band]=population==0 ? -1 : candidateSize(baseline, pool.array[(int)order[start]]);
			bands.maxSize[band]=population==0 ? -1 : candidateSize(baseline, pool.array[(int)order[start+population-1]]);
			// Overflow is intentional seed mixing, never a count or array index.
			bands.seed[band]=operationSeed^(STRATUM_SALT*(band+1));
			start+=population;
		}
		int out=offset;
		for(int band=0, start=0; band<STRATA; band++){
			final int population=bands.population[band];
			final Random random=new Random(bands.seed[band]);
			for(int i=0; i<bands.quota[band]; i++){
				final int at=start+i+random.nextInt(population-i);
				final long selected=order[at]; order[at]=order[start+i]; order[start+i]=selected;
				final long chosen=pool.array[(int)selected];
				left[out]=(int)(chosen>>>32); right[out]=(int)chosen; out++;
			}
			start+=population;
		}
		assert(out==offset+count) : "Each rank band contributes exactly its declared quota";
		return bands;
	}

	/** Sizes come from current assigned memberships, never original codec counts. */
	private static long candidateSize(final LocalProbeBaseline baseline, final long candidate){
		final int left=(int)(candidate>>>32), right=(int)candidate;
		final long size=baseline.selected[baseline.selectedIndex(left)].size();
		return size+(right<0 ? 0L : baseline.selected[baseline.selectedIndex(right)].size());
	}

	/** Only ten bands of scalar metadata survive sampling; sorted candidate scratch does not. */
	private static final class Bands{
		final int[] population=new int[STRATA], quota=new int[STRATA];
		final long[] minSize=new long[STRATA], maxSize=new long[STRATA], seed=new long[STRATA];
	}

	/** Checks the cached selector's population and geometry before any random draws. */
	private static void validate(final LocalProbeBaseline baseline, final FamilyNeighborTable table,
			final LayerCakeClustering.Options options, final TailAlignmentGate gate){
		require(baseline!=null && table!=null && options!=null && gate!=null,
			"Probe sampling requires the frozen baseline, fusion table, structural options and merge gate");
		require(table.size()==baseline.selected.length && table.width()==baseline.sourceHeader.width &&
			table.maxSizeRatio()>0 && table.maxSizeRatio()==options.fusionSizeRatio &&
			table.shortlistLimit()==options.mergeShortlist,
			"Fusion sampling requires the matching roster/width and size-prefiltered shortlist, never an unfiltered cache");
		final TailAlignmentGate measured=table.geometry();
		require(measured.mode==gate.mode && measured.denominator==gate.denominator &&
			measured.maxLength==gate.maxLength && measured.minSpan==gate.minSpan && measured.minPaired==gate.minPaired,
			"Cached native alignment geometry must match the probe merge gate");
		for(int row=0; row<table.size(); row++){
			final int id=table.familyId(row), at=baseline.selectedIndex(id);
			require(at>=0 && at<baseline.selected.length && baseline.selected[at]!=null &&
				baseline.selected[at].familyId==id && baseline.selected[at].size()==table.currentMembers(row),
				"Neighbor table must contain exactly the selected stable IDs and their current membership counts");
		}
	}

	/** Native acceptance cutoffs apply to cached raw metrics, not a clamped quality score. */
	private static boolean passes(final FamilyNeighborTable table, final int edge, final TailAlignmentGate gate){
		return table.reason(edge)==TailAlignmentDecision.Reason.ACCEPT &&
			table.identity(edge)>=gate.minIdentity && table.scoreRatio(edge)>=gate.minScoreRatio;
	}

	/** Partial Fisher-Yates gives a uniform ordered sample and stable prefixes as count grows.
	 * Pools are private scratch. Separate streams ensure one operation's requested
	 * count cannot change the other's selected endpoints. Random.nextInt uses its
	 * bounded sampler rather than a biased remainder of an arbitrary random word.
	 */
	private static void draw(final LongList pool, final int count, final Random random,
			final int[] left, final int[] right, final int offset){
		assert(count>=0 && count<=pool.size && offset>=0 && offset+count<=left.length && left.length==right.length) :
			"Validated no-replacement draw must fit both its eligible pool and retained request arrays";
		for(int i=0; i<count; i++){
			final int at=i+random.nextInt(pool.size-i);
			final long chosen=pool.array[at]; pool.array[at]=pool.array[i]; pool.array[i]=chosen;
			left[offset+i]=(int)(chosen>>>32); right[offset+i]=(int)chosen;
		}
	}
	private static long pair(final int left, final int right){return ((long)left<<32)|(right&0xffffffffL);}

	/*--------------------------------------------------------------*/
	/*----------------         Retained plan        ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable requests and their sampling provenance; no response-dependent update API.
	 * The caller records this metadata with its existing run configuration/input
	 * pins. Candidate IDs are zero-based row positions: fissions in draw order,
	 * followed by fusions in draw order. A family may occur in both operation
	 * arms because each probe independently restores the same baseline.
	 */
	public static final class Selection{
		/** Takes exclusive ownership of sampled primitive IDs and records the uncensored pool census. */
		private Selection(final int[] left_, final int[] right_, final long seed_,
				final int fissions, final int fusions, final int eligibleFissions,
				final int eligibleFusions, final int reciprocals, final long floor,
				final LayerCakeClustering.Options options, final TailAlignmentGate gate,
				final Bands splitBands_, final Bands mergeBands_){
			assert(left_.length==right_.length && left_.length==(long)fissions+fusions &&
				fissions<=eligibleFissions && fusions<=eligibleFusions && eligibleFusions<=reciprocals) :
				"Retained requests must match exact requested counts and fit the precomputed eligible populations";
			left=left_; right=right_; seed=seed_;
			fissionCount=fissions; fusionCount=fusions;
			fissionPoolSize=eligibleFissions; fusionPoolSize=eligibleFusions; reciprocalFusionPoolSize=reciprocals;
			fissionFloor=floor; fissionMinSize=options.fissionMinSize; minChildSize=options.minChildSize;
			fusionSizeRatio=options.fusionSizeRatio; mergeShortlist=options.mergeShortlist;
			mergeMinIdentity=gate.minIdentity; mergeMinScoreRatio=gate.minScoreRatio;
			splitBands=splitBands_; mergeBands=mergeBands_;
		}
		public int size(){return left.length;}
		public boolean isFission(final int index){checkIndex(index); return index<fissionCount;}
		public int leftId(final int index){checkIndex(index); return left[index];}
		public int rightId(final int index){checkIndex(index); return right[index];}
		public String samplingMethod(){return stratified() ? STRATIFIED_METHOD : SAMPLING_METHOD;}
		public long fissionSeed(){return seed^FISSION_SALT;}
		public long fusionSeed(){return seed^FUSION_SALT;}
		public boolean stratified(){return splitBands!=null;}
		public int strata(){return stratified() ? STRATA : 0;}
		public int stratumPopulation(final boolean fission, final int band){return bands(fission, band).population[band];}
		public int stratumQuota(final boolean fission, final int band){return bands(fission, band).quota[band];}
		/** Inclusive size extrema; -1 for an empty band. */
		public long stratumMinSize(final boolean fission, final int band){return bands(fission, band).minSize[band];}
		public long stratumMaxSize(final boolean fission, final int band){return bands(fission, band).maxSize[band];}
		public long stratumSeed(final boolean fission, final int band){return bands(fission, band).seed[band];}

		/** Zero-based band of a retained request, or -1 for the legacy uniform method. */
		public int stratum(final int index){
			checkIndex(index);
			if(!stratified()){return -1;}
			final boolean fission=index<fissionCount;
			int position=fission ? index : index-fissionCount;
			for(int band=0; band<STRATA; band++){
				final int quota=stratumQuota(fission, band);
				if(position<quota){return band;}
				position-=quota;
			}
			throw new IllegalStateException("Retained request must lie in one sampled rank band");
		}

		private Bands bands(final boolean fission, final int band){
			require(stratified() && band>=0 && band<STRATA, "Stratum metadata requires a stratified selection and band 0..4");
			return fission ? splitBands : mergeBands;
		}

		/** Emits exactly the established five-column format. sample_seed records the
		 * master seed, not a per-proposal RNG seed. The driver preflights output paths
		 * against inputs; overwrite is explicit. Sampling already completed before
		 * this method opens anything, so over-request failures cannot truncate output.
		 */
		public void writeRequests(final String path, final boolean overwrite){
			require(path!=null && !path.isEmpty(), "A request table needs an output path");
			final ByteStreamWriter writer=new ByteStreamWriter(path, overwrite, false, true);
			writer.start();
			try{
				writer.println(StructuralProbeBatch.REQUEST_HEADER);
				final ByteBuilder row=new ByteBuilder();
				for(int i=0; i<left.length; i++){
					row.clear().append(i).tab().append(i<fissionCount ? "fission" : "fusion").tab()
						.append(left[i]).tab().append(right[i]).tab().append(seed);
					writer.println(row);
				}
			}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Probe request write failed: "+path);}}
		}
		private void checkIndex(final int index){
			if(index<0 || index>=left.length){throw new IllegalArgumentException("Unknown sampled request row: "+index);}
		}

		private final int[] left, right;
		private final Bands splitBands, mergeBands;
		public final long seed, fissionFloor;
		public final int fissionCount, fusionCount, fissionPoolSize, fusionPoolSize, reciprocalFusionPoolSize;
		public final int fissionMinSize, minChildSize, mergeShortlist;
		public final double fusionSizeRatio, mergeMinIdentity, mergeMinScoreRatio;
	}

	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	/** Fixed salts spell the operation names; retain them with this algorithm version. */
	public static final long FISSION_SALT=0x46495353494F4EL, FUSION_SALT=0x465553494F4EL;
	public static final String SAMPLING_METHOD="uniform_within_operation_partial_fisher_yates_v1";
	public static final int STRATA=5;
	public static final long STRATUM_SALT=0x9e3779b97f4a7c15L;
	public static final String STRATIFIED_METHOD="size_rank_quintiles_partial_fisher_yates_v1";
}
