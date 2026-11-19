package prot;

import java.util.Arrays;
import java.util.HashSet;

/** Immutable resident admission vector and its initial fitting-positive support.
 * Standalone values bind reference identities, including their frozen gates.
 * Classifier-created values also bind the originating resident search policy;
 * attachment checks the released geometry. This is not a file format or a
 * certificate of measured false-positive rates.
 * @author Brian, Keqing
 */
public final class FamilyHitFloorVector{

	/** All-one measurement probe. It has no measured support and cannot be selected
	 * for held-out evaluation or deployment until withInitialSupport is called. */
	public static FamilyHitFloorVector probe(final FamilyReference[] families){
		require(families!=null && families.length>0 && families.length<=Integer.MAX_VALUE/BINS,
			"A floor probe needs a nonempty live roster with seven-bin offsets fitting int");
		final FamilyReference[] roster=families.clone();
		final HashSet<String> ids=new HashSet<String>();
		final int[] ranks=new int[roster.length];
		for(int i=0; i<roster.length; i++){
			final FamilyReference family=roster[i];
			if(family==null || family.members<=0 || !ids.add(family.id)){
				throw new IllegalArgumentException("A floor vector requires unique live family IDs; missing, empty or repeated row="+i);
			}
			ranks[i]=family.rank;
		}
		Arrays.sort(ranks);
		for(int i=0; i<ranks.length; i++){
			require(ranks[i]>=0 && (i==0 || ranks[i]>ranks[i-1]),
				"Physical family ranks must be unique; live ordinals need not equal physical ranks");
		}
		final byte[] floors=new byte[roster.length]; Arrays.fill(floors, (byte)1);
		return new FamilyHitFloorVector(roster, floors, null, null);
	}

	/** Only the classifier creates an executable probe; standalone roster values
	 * intentionally lack this exact resident policy binding. */
	static FamilyHitFloorVector forClassifier(final FamilyTailCascade.SearchConfig config){
		require(config!=null,"An executable floor probe needs its originating classifier policy");
		final FamilyHitFloorVector probe=probe(config.families);
		return new FamilyHitFloorVector(probe.roster,probe.floors,null,config);
	}
	private FamilyHitFloorVector(final FamilyReference[] roster_, final byte[] floors_, final boolean[] support_,
			final FamilyTailCascade.SearchConfig context_){
		roster=roster_; floors=floors_; support=support_; context=context_;
		int min=6, max=1;
		for(byte floor:floors){min=Math.min(min, floor); max=Math.max(max, floor);}
		minimum=min; maximum=max;
	}

	/** Freeze eligibility once from the initial floor-one terminal winners, not
	 * any later replay. Unsupported families immediately return to baseline six.
	 * Each accepted query contributes to one winner; rejected queries contribute none. */
	public FamilyHitFloorVector withInitialSupport(final long positiveQueries, final long[] winnerBins){
		require(support==null, "Initial fitting support is frozen once; later replay winners must not redefine eligibility");
		require(positiveQueries>0 && winnerBins!=null && winnerBins.length==roster.length*BINS,
			"Initial support requires a positive query denominator and seven winner bins per live family");
		final boolean[] measured=new boolean[roster.length];
		final byte[] selected=floors.clone();
		long total=0;
		for(int f=0; f<roster.length; f++){
			for(int hit=0; hit<BINS; hit++){
				final long count=winnerBins[f*BINS+hit];
				if(count<0 || hit==0 && count!=0 || count>positiveQueries-total){
					throw new IllegalArgumentException("Initial floor-one winners must be nonnegative, have positive hits and count each query at most once; family="+f+", bin="+hit);
				}
				total+=count; measured[f]|=count>0;
			}
			if(!measured[f]){selected[f]=6;}
		}
		return new FamilyHitFloorVector(roster, selected, measured, context);
	}

	/** Copy a proposed monotone raise while preserving the initial support mask.
	 * The result still needs exact whole-classifier replay before scientific selection. */
	public FamilyHitFloorVector withFloors(final byte[] proposed){
		require(support!=null, "Proposed selectable floors require measured initial fitting support first");
		require(proposed!=null && proposed.length==floors.length,
			"Proposed floors must cover the complete live reference roster");
		final byte[] copy=proposed.clone();
		for(int i=0; i<copy.length; i++){
			if(copy[i]<floors[i] || copy[i]>6 || !support[i] && copy[i]!=6){
				throw new IllegalArgumentException("Proposals only raise floors within one through six, preserving six for unsupported families; ordinal="+i);
			}
		}
		return new FamilyHitFloorVector(roster, copy, support, context);
	}

	/** Fitting and replay must share the actual index, gates, ranking and margin
	 * objects. Matching IDs alone cannot authorize transfer to a changed policy. */
	boolean boundToPolicy(final FamilyTailCascade.SearchConfig c){
		return context!=null && c!=null && boundTo(c.families) && c.index==context.index &&
			c.geometry==context.geometry && c.width==context.width && c.seed==context.seed &&
			c.top==context.top && c.ddlReject==context.ddlReject && c.ddlAccept==context.ddlAccept &&
			c.positionMode==context.positionMode && c.partitions==context.partitions && c.overlap==context.overlap &&
			c.blacklist==context.blacklist && c.alignmentMetric.same(context.alignmentMetric) &&
			c.scalarBounds==context.scalarBounds && c.scalarGates==context.scalarGates &&
			c.identityFloor==context.identityFloor && c.blosumFloor==context.blosumFloor &&
			c.margins==context.margins && c.heuristicRanking==context.heuristicRanking &&
			c.features==context.features && c.cosineMode==context.cosineMode &&
			c.memberScalarMeans==context.memberScalarMeans && c.lengthRanking==context.lengthRanking && c.lengthMode==context.lengthMode &&
			c.familyAdmission==context.familyAdmission && c.lengthScaledShortAdmission==context.lengthScaledShortAdmission &&
			c.truth==context.truth && c.union==context.union && c.candidateOnly==context.candidateOnly && c.floorAudit==context.floorAudit;
	}

	/** A new array with the same immutable references is allowed; altered gates,
	 * membership, order or any replacement reference require a fresh binding. */
	public boolean boundTo(final FamilyReference[] families){
		if(families==null || families.length!=roster.length){return false;}
		for(int i=0; i<roster.length; i++){if(families[i]!=roster[i]){return false;}}
		return true;
	}

	public int size(){return floors.length;}
	public boolean hasInitialSupport(){return support!=null;}
	public int floor(final int ordinal){
		require(ordinal>=0 && ordinal<floors.length, "A raw-hit floor requires a valid live-family ordinal");
		return floors[ordinal];
	}
	/** Detached selector input; callers cannot mutate this vector. */
	public byte[] copyFloors(){return floors.clone();}
	public boolean[] copyInitialSupport(){
		require(support!=null, "Unmeasured initial support must not become an all-false measured mask");
		return support.clone();
	}
	public int[] copyPhysicalRanks(){
		final int[] ranks=new int[roster.length];
		for(int i=0; i<ranks.length; i++){ranks[i]=roster[i].rank;}
		return ranks;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public final int minimum, maximum;
	private final FamilyReference[] roster;
	private final FamilyTailCascade.SearchConfig context;
	private final byte[] floors;
	/** Null means not measured; a nonnull all-false mask is a real measured result. */
	private final boolean[] support;
	private static final int BINS=FamilyHitFloorSelector.BINS;
}
