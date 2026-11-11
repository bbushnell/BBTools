package assemble;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Optional serial veto for frozen reciprocal fusion proposals. The model sees
 * the same retained context and original-tip entropy as the v2 training corpus.
 * It borrows one phase's counts and never ranks partners or overrides a guard.
 * Instances reuse mutable network, feature and sequence scratch; all calls are
 * serial, with stable contigs and counts until scoring completes. Model dimensions
 * alone do not verify feature ordering, training provenance or calibration.
 * @author Fischl
 */
public final class FusionNeuralGate{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads a fresh inference-only network from text; the caller supplies its chosen cutoff. */
	public FusionNeuralGate(final String path, final float cutoff_){
		this(CellNetParser.loadInferenceFromLines(ByteFile.toLines(path)), cutoff_);
	}

	/**
	 * Checks the v2 input/output dimensions and finite cutoff in [0,1]. Retains the
	 * supplied network without copying and enables its per-network simdFF flag;
	 * callers must not concurrently use or mutate that network's inference scratch.
	 * This constructor does not change the legacy global network representation mode.
	 */
	public FusionNeuralGate(final CellNet net_, final float cutoff_){
		if(net_==null || net_.numInputs()!=INPUTS || net_.numOutputs()!=1){
			throw new IllegalArgumentException("Fusion neural gate requires 109 inputs and one output.");
		}
		if(!Float.isFinite(cutoff_) || cutoff_<0 || cutoff_>1){
			throw new IllegalArgumentException("fusencutoff must be an explicit finite value in [0,1].");
		}
		net=net_;
		net.simdFF=true;
		cutoff=cutoff_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows a stable count provider at the supplied K and resets per-phase counters.
	 * The caller must end the previous phase first and keep the table live until
	 * endPhase; this method neither loads a table nor enforces ownership elsewhere.
	 */
	public void beginPhase(final int k_, final TadpoleGraph.Counts counts){
		assert(measurer==null) : "Detach the previous fusion phase before borrowing another table.";
		k=k_;
		measurer=new FusionJoinDiagnostic(k, counts);
		evaluations=rejected=featureNanos=inferenceNanos=0;
	}

	/**
	 * Reports completed evaluations and drops the borrowed provider; repeated calls
	 * without a new phase do nothing. The caller clears its table afterward. Network
	 * and scratch buffers remain reusable, and counters remain readable until beginPhase.
	 */
	public void endPhase(){
		if(measurer==null){return;}
		System.err.println("Fusion neural: k="+k+", cutoff="+cutoff+", evaluations="+evaluations+
				", rejected="+rejected+", featureSeconds="+(featureNanos*1e-9)+
				", inferenceSeconds="+(inferenceNanos*1e-9)+".");
		measurer=null;
	}

	/**
	 * Scores one proposal and accepts at score>=cutoff, including equality. The caller
	 * supplies a frozen proposal after its other guards; this method does not deduplicate
	 * calls or apply graph guards. Only completed feature/inference calls enter counters
	 * and timings; an exception propagates without recording an evaluated candidate.
	 */
	boolean accepts(final Contig source, final boolean sr, final int st,
			final Contig dest, final boolean dr, final int dt, final int overlap){
		final long start=System.nanoTime();
		features(source, sr, st, dest, dr, dt, overlap);
		final long measured=System.nanoTime();
		final float score=score(vector);
		featureNanos+=measured-start;
		inferenceNanos+=System.nanoTime()-measured;
		evaluations++;
		final boolean accepted=score>=cutoff;
		if(!accepted){rejected++;}
		return accepted;
	}

	/**
	 * Runs native inference with configured input normalization, independent of phase
	 * attachment and candidate counters. Inputs must contain exactly109 finite values;
	 * CellNet copies them into its mutable scratch. Rejects nonfinite/out-of-range output.
	 */
	public float score(final float[] inputs){
		assert(inputs.length==INPUTS) : "Fusion v2 inference must consume exactly its 109 trained inputs.";
		for(float input : inputs){
			if(!Float.isFinite(input)){throw new IllegalArgumentException("Nonfinite fusion network input.");}
		}
		final float score=net.applyInput(inputs).feedForward();
		if(!Float.isFinite(score) || score<0 || score>1){
			throw new IllegalStateException("Fusion network output must be finite and in [0,1]: "+score);
		}
		return score;
	}

	/**
	 * Returns a borrowed vector overwritten by subsequent feature calls. Matches
	 * FusionJoinCollector's 200-base retained flanks and FusionEntropyCorpus's
	 * original, pre-trim terminal100-base tips.
	 * Nine-decimal quantization reproduces the persisted training input floats.
	 * Reverse flags select oriented contig views; trims and overlap are base counts in
	 * those views. Requires a positive exact overlap and at least one retained flanking
	 * base on each side. Input contigs are read, not changed; scoring/counters are separate.
	 */
	public float[] features(final Contig source, final boolean sr, final int st,
			final Contig dest, final boolean dr, final int dt, final int overlap){
		assert(measurer!=null && source!=null && dest!=null) : "Live fusion features need contigs and phase counts.";
		final int sourceEnd=source.length()-st, destStart=dt+overlap;
		if(st<0 || dt<0 || overlap<1 || sourceEnd<=overlap || destStart>=dest.length()){
			throw new IllegalArgumentException("Invalid retained fusion geometry.");
		}
		final int left=Math.min(200, sourceEnd-overlap), right=Math.min(200, dest.length()-destStart);
		final int length=left+overlap+right, start=sourceEnd-overlap-left;
		if(product.length<length){product=Arrays.copyOf(product, length);}
		for(int i=0; i<left+overlap; i++){product[i]=baseAt(source, sr, start+i);}
		for(int i=0; i<overlap; i++){
			if(product[left+i]!=baseAt(dest, dr, dt+i)){
				throw new IllegalArgumentException("Fusion network received a nonexact overlap.");
			}
		}
		for(int i=0; i<right; i++){product[left+overlap+i]=baseAt(dest, dr, destStart+i);}
		final float[] core=measurer.measure(product, length, left, overlap, st, dt, null);
		assert(core.length==INPUTS-2) : "The v2 entropy inputs extend, rather than replace, all107 v1 values.";
		for(int i=0; i<core.length; i++){vector[i]=quantize(core[i]);}
		sourceTip=tip(source, sr, true, sourceTip);
		destTip=tip(dest, dr, false, destTip);
		final float a=entropy.mean(sourceTip, true), b=entropy.mean(destTip, false);
		vector[INPUTS-2]=quantize(Math.min(a, b));
		vector[INPUTS-1]=quantize(Math.max(a, b));
		return vector;
	}

	/** Reuses exact-sized tip buffers, including rare contigs shorter than100 bases. */
	private static byte[] tip(final Contig contig, final boolean reverse, final boolean right, byte[] buffer){
		assert(contig.length()>0) : "Only nonempty contigs have a fusion tip.";
		final int length=Math.min(FusionTipEntropy.TIP_BASES, contig.length());
		if(buffer.length!=length){buffer=new byte[length];}
		final int start=right ? contig.length()-length : 0;
		for(int i=0; i<length; i++){buffer[i]=baseAt(contig, reverse, start+i);}
		return buffer;
	}

	/** Matches native vector serialization/parsing without a temporary String per input. */
	private float quantize(final float value){
		decimal.clear().append(value, 9);
		return Parse.parseFloat(decimal.array, 0, decimal.length());
	}

	/** Reads an oriented base without reversing or copying the complete contig. */
	private static byte baseAt(final Contig contig, final boolean reverse, final int pos){
		return reverse ? AminoAcid.baseToComplementExtended[contig.bases[contig.length()-1-pos]] : contig.bases[pos];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained mutable inference state, owned here or supplied for exclusive serial use. */
	private final CellNet net;
	/** Inclusive acceptance threshold. */
	private final float cutoff;
	/** Reused original-tip measurement scratch. */
	private final FusionTipEntropy entropy=new FusionTipEntropy();
	/** Reused nine-decimal serialization scratch. */
	private final ByteBuilder decimal=new ByteBuilder(32);
	/** Borrowed109-value feature result; never returned as a fresh snapshot. */
	private final float[] vector=new float[INPUTS];
	/** Product capacity may grow; tip buffers are resized to each exact terminal span. */
	private byte[] product=new byte[0], sourceTip=new byte[0], destTip=new byte[0];
	/** Current borrowed count provider and feature scratch, or null outside a phase. */
	private FusionJoinDiagnostic measurer;
	/** Count-table K for the current/most recently ended phase. */
	private int k;
	/** Completed accepts calls in the current or most recently ended phase. */
	public long evaluations;
	/** Phase veto count and accumulated feature/inference time in nanoseconds. */
	long rejected, featureNanos, inferenceNanos;
	/**107 core features followed by symmetric original-tip entropy minimum/maximum. */
	static final int INPUTS=109;
}
