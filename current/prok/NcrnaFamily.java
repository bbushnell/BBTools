package prok;

import consensus.BaseGraph;
import map.LongHashSet;
import ml.CellNet;

/**
 * One conserved-ncRNA family's loaded resources + per-family scavenger config,
 * as a shared (read-only after loading), family-agnostic bundle. Populated by
 * CallGenes.loadNcrnaResources(); consumed by GeneCaller, which builds one
 * NcrnaScavenger per bundle per instance (mirroring TrnaCaller's per-instance
 * construction -- NcrnaScavenger's TrnaKmerIndex is per-thread-mutable, so the
 * bundle itself is shared but each GeneCaller instance gets its own scavenger).
 *
 * <p>Two bundles may share the same kmerSet/kLong (e.g. srp_small and
 * srp_large both seed from one combined srp_17mers.fa, per Noire, 2026-08-23:
 * "mixing 120bp and 300bp sequences in one consensus doesn't hold together"
 * -- so the CONSENSUS libraries are separate, but the conserved-kmer SEED set
 * is shared across both size classes).
 * @author Neptune
 */
public class NcrnaFamily {

	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				7, 60, true, 11f, 0.48f, 0.072f, 12);
	}

	/** Forward-ported from Noire's ncRNA-family-loading tree (2026-08-28, C3 merge): scoreA/
	 * scoreB/idPass/idBorderline replace the flat bestId*100 orfScore with a per-family-tuned
	 * scoreA+scoreB*orfLen*bestId*bestId formula in NcrnaScavenger.alignWindow -- Noire's #1
	 * recall lever (+11.5pp rnasep, per her findings). Defaults here (0f,20f,0.75f,0.65f)
	 * reproduce the OLD flat bestId*100-equivalent shape only in the trivial sense that no
	 * caller in this tree constructs a family via this delegation path with real ncRNA data
	 * yet -- callers that DO load real families must pass real per-family values (see
	 * CallGenes.loadNcrnaResources). Does NOT touch TrnaCaller's tRNA orfScore (bestId*100,
	 * unrelated, untouched) -- this formula is NcrnaScavenger-only. */
	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				0f, 20f, 0.75f, 0.65f);
	}

	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_,
				null, null, null, null, -1, -1, -1, -1, 0f);
	}

	/** Full constructor, adding C3's boundary-precision-NN resources (Noire's spec,
	 * plans/c3_ncrnaboundaryscorer_spec.md; G11, 2026-08-28). OFF-by-default is
	 * STRUCTURAL, not a separate boolean on this class: boundary5NetTemplate==null means
	 * off -- every shorter constructor above delegates here with all-null/sentinel
	 * boundary args, so a family built any other way than this exact overload has boundary
	 * refinement permanently off. NcrnaScavenger checks boundary5NetTemplate!=null (not a
	 * flag on this class) before ever cloning a net or calling refineBoundaries -- see its
	 * constructor and trimToAlignmentExtent. meanLen is REQUIRED whenever a net is given
	 * (asserted in NcrnaScavenger's constructor) -- families differ hugely in length
	 * (rnasep ~380bp vs srp_small ~95bp), so there is no safe shared default. */
	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate_, CellNet boundary3NetTemplate_,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_, 0.75f, 0.9f,
				boundary5NetTemplate_, boundary3NetTemplate_, boundaryStartTable_, boundaryStopTable_,
				boundaryStartInside_, boundaryStartOutside_, boundaryStopInside_, boundaryStopOutside_, boundaryMeanLen_);
	}

	/** Full family configuration, including the exact boundary candidates used by both
	 * vector generation and inference. Arrays are validated and cloned so family bundles
	 * remain read-only after loading. */
	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			float hbmPass_, float collapseFrac_,
			CellNet boundary5NetTemplate_, CellNet boundary3NetTemplate_,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_, int[] boundaryStartOffsets_, int[] boundaryStopOffsets_){
		name=name_;
		library=library_;
		models=models_;
		modelNames=modelNames_;
		kmerSet=kmerSet_;
		kLong=kLong_;
		minLen=minLen_;
		windowPad=windowPad_;
		indexK=indexK_;
		indexTopN=indexTopN_;
		adaptive=adaptive_;
		adaptFloor=adaptFloor_;
		adaptTopFrac=adaptTopFrac_;
		adaptQFrac=adaptQFrac_;
		fixedMinHits=fixedMinHits_;
		scoreA=scoreA_;
		scoreB=scoreB_;
		idPass=idPass_;
		idBorderline=idBorderline_;
		hbmPass=hbmPass_;
		collapseFrac=collapseFrac_;
		boundary5NetTemplate=boundary5NetTemplate_;
		boundary3NetTemplate=boundary3NetTemplate_;
		boundaryStartTable=boundaryStartTable_;
		boundaryStopTable=boundaryStopTable_;
		boundaryStartInside=boundaryStartInside_;
		boundaryStartOutside=boundaryStartOutside_;
		boundaryStopInside=boundaryStopInside_;
		boundaryStopOutside=boundaryStopOutside_;
		boundaryMeanLen=boundaryMeanLen_;
		boundaryStartOffsets=validatedOffsets(boundaryStartOffsets_, "start");
		boundaryStopOffsets=validatedOffsets(boundaryStopOffsets_, "stop");
	}

	/** Full family configuration including the two remaining alignment/candidate tunables. */
	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			float hbmPass_, float collapseFrac_,
			CellNet boundary5NetTemplate_, CellNet boundary3NetTemplate_,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_, hbmPass_, collapseFrac_,
				boundary5NetTemplate_, boundary3NetTemplate_, boundaryStartTable_, boundaryStopTable_,
				boundaryStartInside_, boundaryStartOutside_, boundaryStopInside_, boundaryStopOutside_, boundaryMeanLen_,
				LEGACY_START_OFFSETS, LEGACY_STOP_OFFSETS);
	}

	private static int[] validatedOffsets(int[] offsets, String label){
		if(offsets==null || offsets.length<1){throw new IllegalArgumentException("Empty ncRNA boundary "+label+" offsets");}
		final int[] copy=offsets.clone();
		boolean hasZero=false;
		for(int i=0; i<copy.length; i++){
			if(copy[i]==0){hasZero=true;}
			if(i>0 && copy[i]<=copy[i-1]){throw new IllegalArgumentException(
				"ncRNA boundary "+label+" offsets must be strictly increasing: "+java.util.Arrays.toString(copy));}
		}
		if(!hasZero){throw new IllegalArgumentException("ncRNA boundary "+label+" offsets must include 0: "
			+java.util.Arrays.toString(copy));}
		return copy;
	}

	static final int[] LEGACY_START_OFFSETS={-3,-2,-1,0,1,2};
	static final int[] LEGACY_STOP_OFFSETS={-4,-3,-2,-1,0,1};

	public final String name;
	public final byte[][] library;
	public final BaseGraph[] models;
	public final String[] modelNames;
	public final LongHashSet kmerSet;
	public final int kLong;
	public final int minLen;
	public final int windowPad;

	public final int indexK;
	public final int indexTopN;
	public final boolean adaptive;
	public final float adaptFloor;
	public final float adaptTopFrac;
	public final float adaptQFrac;
	public final int fixedMinHits;

	/** Forward-ported from Noire's tree (2026-08-28) -- see the 15-arg constructor's javadoc. */
	public final float scoreA;
	public final float scoreB;
	public final float idPass;
	public final float idBorderline;
	public final float hbmPass;
	public final float collapseFrac;

	/** C3 boundary-precision-NN resources (G11, 2026-08-28) -- null/-1/0f (this class' default
	 * for every constructor except the full one) means OFF. Read-only shared templates; each
	 * NcrnaScavenger instance clones its own net copies (thread safety), mirroring
	 * TrnaCaller's BOUNDARY_5_NET_TEMPLATE/BOUNDARY_3_NET_TEMPLATE pattern. */
	public final CellNet boundary5NetTemplate;
	public final CellNet boundary3NetTemplate;
	public final TrnaBoundaryFeatures.NinemerTable boundaryStartTable;
	public final TrnaBoundaryFeatures.NinemerTable boundaryStopTable;
	public final int boundaryStartInside;
	public final int boundaryStartOutside;
	public final int boundaryStopInside;
	public final int boundaryStopOutside;
	public final float boundaryMeanLen;
	final int[] boundaryStartOffsets;
	final int[] boundaryStopOffsets;
}
