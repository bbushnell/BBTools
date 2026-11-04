package prok;

import java.util.Arrays;

import consensus.BaseGraph;
import idaligner.QuantumAligner;
import prok.NcrnaEndpointTable.LeaveOneOut;
import prok.TrnaBoundaryFeatures.NinemerTable;

/** Worker-local model selection and candidate-specific V1 boundary features.
 * Uses the accepted family's mapping shortlist policy and Quantum for model
 * selection. Feature alignment calls the exact shared inference function; it
 * must not be substituted independently of the shipped caller.
 * Each prepare borrows one context until the next prepare. The ten-value array
 * is reusable scratch and must be consumed before the next fill call.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryExamples{

	/** Borrows the accepted family and completed tables; owns shortlist and leave-one-out scratch. */
	public NcrnaBoundaryExamples(final NcrnaFamily family_, final NcrnaEndpointTable start_,
			final NcrnaEndpointTable stop_){
		if(family_==null || family_.library==null || family_.library.length==0
				|| family_.models==null || family_.models.length!=family_.library.length
				|| start_==null || stop_==null || start_.type!=TrnaBoundaryFeatures.BoundaryType.START
				|| stop_.type!=TrnaBoundaryFeatures.BoundaryType.STOP){
			throw new IllegalArgumentException("Boundary examples require an accepted consensus/HBM pair and both endpoint tables");
		}
		family=family_;
		start=start_;
		stop=stop_;
		startLoo=start.leaveOneOut();
		stopLoo=stop.leaveOneOut();
		index=new TrnaKmerIndex(family.library, family.indexK, family.adaptive,
			family.adaptFloor, family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits);
		plateau=new float[family.library.length];
		upperMedianLength=PackagedBoundaryConfig.upperMedianLength(family.library);
		for(int model=0; model<family.library.length; model++){
			final BaseGraph graph=family.models[model];
			if(graph==null || !Arrays.equals(graph.original, family.library[model])){
				throw new IllegalArgumentException("Accepted consensus/HBM sequences disagree at model "+model);
			}
			plateau[model]=TrnaBoundaryFeatures.modelMaxCoverage(graph);
			if(plateau[model]<=0){throw new IllegalArgumentException("Boundary models need positive training coverage");}
		}
	}

	/** Selects the highest-identity shortlisted model; empty shortlist is an explicit skip.
	 * Caller supplies only eligible loci, and enables LOO only for records counted
	 * in these tables (training partition for tuning, all eligible loci for refit).
	 */
	public boolean prepare(final NcrnaLearningContext context_, final boolean leaveOneOut){
		if(context_==null || !context_.hasFeatures()){throw new IllegalArgumentException("Vector preparation requires valid span geometry and whole-contig GC");}
		context=context_;
		model=-1;
		final byte[] gene=Arrays.copyOfRange(context.read.bases, context.start, context.stop+1);
		final int[] shortlist=index.shortlist(gene, family.indexTopN);
		float best=-1;
		for(int candidate:shortlist){
			final float identity=QuantumAligner.alignStatic(gene, family.library[candidate], null);
			if(!(identity>=0 && identity<=1)){throw new IllegalStateException("Invalid Quantum identity during boundary model selection");}
			if(identity>best){
				best=identity;
				model=candidate;
			}
		}
		if(model<0){return false;}
		if(leaveOneOut){
			startLoo.set(context.read.bases, context.start);
			stopLoo.set(context.read.bases, context.stop);
		}
		startCounts=leaveOneOut ? startLoo : start;
		stopCounts=leaveOneOut ? stopLoo : stop;
		return true;
	}

	/** Fills one candidate; returns a named limitation instead of emitting an invalid span. */
	public int fill(final boolean isStop, final int offset){
		if(model<0){throw new IllegalStateException("Select a model before generating boundary candidates");}
		final long left=(long)context.start+(isStop ? 0 : offset);
		final long right=(long)context.stop+(isStop ? offset : 0);
		if(left<0){return BEFORE_CONTEXT;}
		if(right>=context.read.length()){return AFTER_CONTEXT;}
		if(right-left+1<16){return TOO_SHORT;}
		final int s=(int)left, e=(int)right;
		final NcrnaEndpointTable table=isStop ? stop : start;
		final NinemerTable counts=isStop ? stopCounts : startCounts;
		final byte[] candidate=Arrays.copyOfRange(context.read.bases, s, e+1);
		TrnaBoundaryFeatures.aniAndFuzzinessFeature(candidate,
			family.library[model], family.models[model], !isStop, plateau[model], aniFuzz);
		TrnaBoundaryFeatures.enrichmentProfile(context.read.bases,
			isStop ? e : s, table.type, table.inside, table.outside, counts, profile);
		values[0]=aniFuzz[0];
		values[1]=profile[0];
		values[2]=profile[1];
		values[3]=profile[2];
		values[4]=isStop ? 1 : 0;
		values[5]=aniFuzz[1];
		values[6]=aniFuzz[2];
		values[7]=aniFuzz[3];
		values[8]=(e-s+1)/upperMedianLength;
		values[9]=context.gc;
		label=offset==0 ? 1 : 0;
		assert(finiteFeatures(offset));
		return VALID;
	}

	/** Runs the complete feature check only when assertions are enabled. */
	private boolean finiteFeatures(final int offset){
		for(int i=0; i<values.length; i++){
			assert(Float.isFinite(values[i])) : "Training requires finite features: index="+i+", value="+values[i]+", model="+model+", offset="+offset;
		}
		return true;
	}

	/** Selected accepted-library model index, or -1 before a successful preparation. */
	public int model(){return model;}
	/** Upper median accepted-library length, preserving the inference normalization contract. */
	public final float upperMedianLength;
	/** Ten candidate features, overwritten on every successful fill; consume before reuse. */
	public final float[] values=new float[10];
	/** One for the original endpoint, zero for a shifted endpoint; valid after successful fill. */
	public int label;
	/** Worker-local feature outputs, reused across both ends and successive loci. */
	private final float[] aniFuzz=new float[4], profile=new float[3];
	private final NcrnaFamily family;
	private final TrnaKmerIndex index;
	private final NcrnaEndpointTable start, stop;
	private final LeaveOneOut startLoo, stopLoo;
	private final float[] plateau;
	private NinemerTable startCounts, stopCounts;
	private NcrnaLearningContext context;
	private int model=-1;
	/** fill completed and values/label describe the requested candidate. */
	public static final int VALID=0;
	/** Candidate start would precede the available context. */
	public static final int BEFORE_CONTEXT=1;
	/** Candidate stop would follow the available context. */
	public static final int AFTER_CONTEXT=2;
	/** Candidate span would be shorter than the inference minimum of sixteen bases. */
	public static final int TOO_SHORT=3;
}
