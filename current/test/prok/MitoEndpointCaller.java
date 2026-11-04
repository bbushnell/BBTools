package test.prok;

import prok.NcrnaBoundaryRefiner;
import prok.MitoSsuModelPolicy;
import prok.NcrnaBoundaryScorer;
import prok.NcrnaFamily;
import prok.Orf;
import prok.RrnaEndpointCallerFeatures;
import prok.SeedOffsetTable;

import java.util.Arrays;
import ml.CellNet;

/** Explicit mitochondrial sliding refinement alongside unchanged prok endpoint nets.
 * Resources are shared read-only; every caller worker owns its network state.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointCaller implements NcrnaBoundaryRefiner {
	public static final class Resources implements NcrnaBoundaryRefiner.Factory {
		public Resources(String[] names_,byte[][] refs_,CellNet net_,MitoEndpointTables tables_,int radius_){
			require(names_!=null && refs_!=null && names_.length==158 && refs_.length==names_.length
				&& MitoSsuModelPolicy.coverage("16S",names_),"Mito endpoint resources require the unchanged original five plus153named coverage models");
			require(net_!=null && net_.numInputs()==10 && net_.numOutputs()==1 && tables_!=null && tables_.models==refs_.length,
				"One shared10-input/1-output net must bind the complete model-indexed count resource");
			require(radius_==128 || radius_==192 || radius_==256,"Mito inference uses the selected192or the explicit128/256comparison controls");
			require(tables_.meanLength==MitoEndpointTables.medianLength(refs_,5),"Training and caller must share one mitochondrial length denominator");
			names=names_.clone();refs=refs_.clone();net=net_;tables=tables_;radius=radius_;
			offsets=new int[2*radius+1];for(int i=0;i<offsets.length;i++){offsets[i]=i-radius;}
		}
		public static Resources load(NcrnaFamily family,String net,String tables,String pin,int radius){
			require(family!=null && "16S".equals(family.name) && SeedOffsetTable.validSha80(pin),"Explicit experimental resources must bind the selected16S library and recorded sha80");
			return new Resources(family.modelNames,family.library,NcrnaBoundaryScorer.load(net),MitoEndpointTables.load(tables,family.library.length,pin),radius);
		}
		/** Validate the historical experiment's raw-anchor policy once per worker.
		 * Original prokaryotic endpoint resources may coexist; applies() leaves their
		 * models on that path. Mutable network state is copied for this worker only. */
		@Override
		public MitoEndpointCaller create(final NcrnaFamily family){
			require(family!=null && family.reuseConsensusAlignment && !family.trimAlignmentExtent
				&& !family.voteEnds && !family.voteWindows && family.models==null
				&& !family.pacBioConsensusAlignment && !family.modelEndClipping && !family.modelClipRescue
				&& family.boundaryNetsByModel==null && family.boundary5NetTemplate==null
				&& family.boundary3NetTemplate==null,
				"Mito sliding refinement requires raw full-model Quantum anchors and no competing legacy boundary policy");
			final MitoEndpointCaller caller=new MitoEndpointCaller(this);
			caller.validateLibrary(family.name, family.modelNames, family.library);
			return caller;
		}
		public final String[] names;public final byte[][] refs;public final CellNet net;
		public final MitoEndpointTables tables;public final int radius;final int[] offsets;
	}
	public MitoEndpointCaller(Resources resources_){
		require(resources_!=null,"Per-worker mitochondrial inference requires explicit resources");resources=resources_;net=resources.net.copy(false);
	}
	@Override
	public boolean applies(int model){
		require(model>=0 && model<resources.refs.length,"Winning model must belong to the bound16S library");return model>=5;
	}
	public void validateLibrary(String family,String[] names,byte[][] refs){
		require("16S".equals(family) && Arrays.equals(names,resources.names) && refs.length==resources.refs.length,
			"Mito endpoint model order must equal the actual caller order");
		for(int i=0;i<refs.length;i++){require(Arrays.equals(refs[i],resources.refs[i]),"Mito endpoint consensus bases differ from caller model "+i);}
	}
	/** Bases are the complete sense-oriented source contig, not a cropped window.
	 * Original five models, partial calls and non-Quantum anchors remain unchanged. */
	@Override
	public boolean refine(Orf orf,byte[] bases,int model,float gc,float identity,boolean quantumIdentity,int fixedMask,int minLen,int maxLen){
		if(!applies(model)){return false;}
		require(orf!=null && bases!=null && "16S".equals(orf.ncrnaFamily) && resources.names[model].equals(orf.trnaModel)
			&& (fixedMask&~3)==0,"Mito refinement requires the caller's strand-oriented model and physical-end provenance");
		if(fixedMask!=0 || orf.ncrnaPartial!=0 || !quantumIdentity){return false;}
		require(orf.start>=0 && orf.stop<bases.length && orf.stop-orf.start>=15 && Float.isFinite(identity) && identity>=0 && identity<=1,
			"Complete mitochondrial refinement requires a valid full-model Quantum anchor");
		final MitoEndpointTables t=resources.tables;
		final NcrnaBoundaryScorer.RefinementResult r=NcrnaBoundaryScorer.refineBoundariesWithCounts(net,net,bases,orf.start,orf.stop,
			resources.refs[model],null,t.table(model,0),t.table(model,1),t.insideCount(0),9-t.insideCount(0),t.insideCount(1),9-t.insideCount(1),
			gc,t.meanLength,resources.offsets,resources.offsets,0,0,true,true,NcrnaBoundaryScorer.FEATURES_V3,identity,null,null,t.table(-1,0),t.table(-1,1));
		require(Float.isFinite(r.startChosenScore) && Float.isFinite(r.stopChosenScore),"Nonfinite learned scores cannot become plausible unchanged calls");
		// Both window edges above are genuine contig edges because no subwindow is copied.
		// The established rRNA fallback retains the raw call if a joint span is invalid.
		return RrnaEndpointCallerFeatures.apply(orf,orf.start+r.startOffset,orf.stop+r.stopOffset,bases.length,minLen,maxLen);
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final Resources resources;private final CellNet net;
}
