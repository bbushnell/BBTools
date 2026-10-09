package aligner;

import java.util.Arrays;

/** Exact CYK decision on a caller's already refined, strand-oriented
 * candidate extent. The extent and caller's path score are not changed here.
 * Rfam GA is a provisional engineering threshold, not a calibrated global-CYK
 * operating point; retain candidate score distributions for later selection.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelGate {
	public CovarianceModelGate(CovarianceModel model_, float threshold_, long maxScoreCells){
		this(model_, threshold_, maxScoreCells, false);
	}
	/** Local scoring must be explicitly selected; the supplied threshold remains
	 * an experimental choice, never recalibrated by switching scoring modes. */
	public CovarianceModelGate(CovarianceModel model_, float threshold_, long maxScoreCells, boolean local_){
		require(model_!=null && Float.isFinite(threshold_), "A CM decision requires an explicit model and finite bit-score threshold");
		model=model_;threshold=threshold_;local=local_;scorer=new CovarianceModelScoreOnly(model, maxScoreCells, local);
	}
	/** Read the actual pinned model header; never substitute a guessed family cutoff. */
	public static float gatheringThreshold(CovarianceModel model){
		require(model!=null, "A gathering threshold belongs to a specific parsed model");
		final String text=model.header("GA");require(text!=null, "Model lacks GA; declare an explicit provisional threshold: "+model.accession);
		final float value=Float.parseFloat(text.trim());require(Float.isFinite(value), "Model GA must be a finite bit score: "+model.accession);return value;
	}
	/** Inclusive coordinates in a strand-oriented array. Unsupported alphabets
 * and resource errors propagate; neither may masquerade as a scored rejection. */
	public Decision evaluate(byte[] bases,int start,int stop){
		require(bases!=null && start>=0 && start<=stop && stop<bases.length,
			"CM candidate bounds must select a nonempty existing interval; start="+start+", stop="+stop);
		final long begin=System.nanoTime();
		final byte[] sequence=Arrays.copyOfRange(bases,start,stop+1);
		final CovarianceModelScoreOnly.Result result=scorer.align(sequence);
		final double seconds=(System.nanoTime()-begin)*1e-9;
		require(Float.isFinite(result.score) || result.score==Float.NEGATIVE_INFINITY,
			"Invalid CM arithmetic cannot decide a biological candidate");
		return new Decision(result.score,Float.isFinite(result.score) && result.score>=threshold,sequence.length,seconds,result.peakScoreCells);
	}
	public static final class Decision {
		Decision(float bits_,boolean accepted_,int length_,double seconds_,long cells_){
			assert(length_>0 && seconds_>=0 && cells_>0):"Decisions describe one completed nonempty exact-scoring attempt";
			bits=bits_;accepted=accepted_;length=length_;seconds=seconds_;peakScoreCells=cells_;
		}
		public final float bits;
		public final boolean accepted;
		public final int length;
		public final double seconds;
		public final long peakScoreCells;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final CovarianceModel model;
	public final float threshold;
	public final boolean local;
	private final CovarianceModelScoreOnly scorer;
}
