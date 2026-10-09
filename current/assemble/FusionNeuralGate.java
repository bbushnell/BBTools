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
 * @author Fischl
 */
final class FusionNeuralGate {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads isolated inference state; the caller must supply its chosen cutoff. */
	FusionNeuralGate(final String path, final float cutoff_){
		this(CellNetParser.loadInferenceFromLines(ByteFile.toLines(path)), cutoff_);
	}

	/** Validates the fixed v2 schema without changing legacy global network mode. */
	FusionNeuralGate(final CellNet net_, final float cutoff_){
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

	/** Attaches the only current table and resets per-phase cost/candidate counters. */
	void beginPhase(final int k_, final TadpoleGraph.Counts counts){
		assert(measurer==null) : "Detach the previous fusion phase before borrowing another table.";
		k=k_;
		measurer=new FusionJoinDiagnostic(k, counts);
		evaluations=rejected=featureNanos=inferenceNanos=0;
	}

	/** Drops borrowed counts before the owning table is cleared; never clears a table itself. */
	void endPhase(){
		if(measurer==null){return;}
		System.err.println("Fusion neural: k="+k+", cutoff="+cutoff+", evaluations="+evaluations+
				", rejected="+rejected+", featureSeconds="+(featureNanos*1e-9)+
				", inferenceSeconds="+(inferenceNanos*1e-9)+".");
		measurer=null;
	}

	/** Scores an unchanged proposal once, before any contigs in its frozen pass are merged. */
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

	/** Native inference, including any serialized input normalization. */
	float score(final float[] inputs){
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
	 * Returns a borrowed vector. Matches FusionJoinCollector's 200-base retained
	 * flanks and FusionEntropyCorpus's original, pre-trim terminal100-base tips.
	 * Nine-decimal quantization reproduces the persisted training input floats.
	 */
	float[] features(final Contig source, final boolean sr, final int st,
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

	private final CellNet net;
	private final float cutoff;
	private final FusionTipEntropy entropy=new FusionTipEntropy();
	private final ByteBuilder decimal=new ByteBuilder(32);
	private final float[] vector=new float[INPUTS];
	private byte[] product=new byte[0], sourceTip=new byte[0], destTip=new byte[0];
	private FusionJoinDiagnostic measurer;
	private int k;
	long evaluations, rejected, featureNanos, inferenceNanos;
	static final int INPUTS=109;
}
