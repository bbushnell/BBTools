package aligner;

import java.util.Arrays;
import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;

/** Placement policy with explicit experimental clipped-anchor and below-cutoff
 * rescue options. Existing constructors retain the D32/D34 policies.
 * No saved exact score or oracle decides fallback.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelPlacementGate {
	public CovarianceModelPlacementGate(CovarianceModel model_, float threshold_, int radius_, long scoreLimit, long traceLimit){
		this(model_, threshold_, radius_, scoreLimit, traceLimit, false);
	}
	public CovarianceModelPlacementGate(CovarianceModel model_, float threshold_, int radius_, long scoreLimit, long traceLimit, boolean adaptive_){
		this(model_, threshold_, radius_, scoreLimit, traceLimit, adaptive_, false);
	}
	/** Local scoring is explicit and applies to both banded and exact paths.
	 * Anchor availability and the D32 widening/fallback policy are unchanged. */
	public CovarianceModelPlacementGate(CovarianceModel model_, float threshold_, int radius_, long scoreLimit, long traceLimit, boolean adaptive_, boolean local_){
		this(model_, threshold_, radius_, scoreLimit, traceLimit, adaptive_, local_, false);
	}
	/** D34 extension is explicit until its biological comparison is accepted.
	 * It changes only the placement corridors, never the candidate sequence. */
	public CovarianceModelPlacementGate(CovarianceModel model_, float threshold_, int radius_, long scoreLimit, long traceLimit, boolean adaptive_, boolean local_, boolean extendClipped_){
		this(model_, threshold_, radius_, scoreLimit, traceLimit, adaptive_, local_, extendClipped_, false);
	}
	/** D34b tries extended bands even at physical input edges, then rescues a
	 * below-cutoff score on the identical candidate. This is not a default. */
	public CovarianceModelPlacementGate(CovarianceModel model_, float threshold_, int radius_, long scoreLimit, long traceLimit, boolean adaptive_, boolean local_, boolean extendClipped_, boolean rescueBelowThreshold_){
		require(model_!=null && Float.isFinite(threshold_) && radius_>=0,
			"Placement decisions require a bound model, finite cutoff and declared radius");
		require(!rescueBelowThreshold_ || extendClipped_, "D34b must try the extended clipped-query corridors before score-based rescue");
		model=model_;threshold=threshold_;adaptive=adaptive_;local=local_;radius=adaptive ? 32 : radius_;
		extendClipped=extendClipped_;rescueBelowThreshold=rescueBelowThreshold_;
		banded=new CovarianceModelBandedCyk(model, new CovarianceModelLengthBands(model, 2), scoreLimit, traceLimit, local);
		exact=new CovarianceModelScoreOnly(model, scoreLimit, local);
	}
	/** Coordinates refer to one strand-oriented chromosome; the retained trace
	 * uses the consensus as query and the window beginning at windowStart as ref. */
	public Decision evaluate(byte[] bases, int start, int stop, CovarianceModelConsensusMap map,
			AlignmentStats stats, int windowStart){
		require(bases!=null && start>=0 && start<=stop && stop<bases.length && windowStart>=0,
			"CM placement must consume the existing inclusive refined candidate");
		require(map==null || map.belongsTo(model), "Consensus anchors must belong to this exact loaded CM instance");
		require(stats==null || stats.rLen>=0 && (long)windowStart+stats.rLen<=bases.length,
			"Retained reference-window coordinates must fit the same oriented chromosome");
		final long begin=System.nanoTime();final byte[] sequence=Arrays.copyOfRange(bases, start, stop+1);
		final boolean physicalEdge=extendClipped && !rescueBelowThreshold && clippedAtPhysicalEdge(stats, windowStart, bases.length);
		final CovarianceModelPlacementBands placement=map==null || physicalEdge ? null
			: map.compose(stats, windowStart, start, stop, radius, extendClipped);
		final long composed=System.nanoTime();
		CovarianceModelBandedCyk.Result raw=null;CovarianceModelLengthBands.Layout[] layouts=null;
		int finalRadius=0, attempts=0;long peakScore=0, peakTrace=0;
		if(placement!=null){
			float previous=Float.NaN;
			for(int width=radius; ; width*=2){
				final CovarianceModelPlacementBands current=width==radius ? placement
					: new CovarianceModelPlacementBands(model, sequence.length, placement.low, placement.high, width);
				layouts=current.layouts();raw=banded.align(sequence, layouts, Double.NaN);attempts++;finalRadius=width;
				final float score=raw.alignment.score;
				require(Float.isFinite(score) || score==Float.NEGATIVE_INFINITY, "Invalid arithmetic cannot decide widening or fallback");
				require(Float.isNaN(previous) || score>=previous,
					"Nested placement corridors retain every prior parse, so widening cannot lower the optimum");
				peakScore=Math.max(peakScore, raw.peakScoreCapacity);peakTrace=Math.max(peakTrace, raw.traceCells);
				if(!adaptive || !widen(width, score, previous)){break;}previous=score;
			}
		}
		final CovarianceModelPlacementBands.EdgeContacts contacts=raw==null ? CovarianceModelPlacementBands.NO_CONTACTS
			: CovarianceModelPlacementBands.edgeContacts(raw.alignment, layouts, local);
		final long bandEnd=System.nanoTime();
		final float rawBits=raw==null ? Float.NaN : raw.alignment.score;
		require(raw==null || Float.isFinite(rawBits) || rawBits==Float.NEGATIVE_INFINITY,
			"Invalid arithmetic is not a legitimate missing-parse fallback");
		final String reason=map==null ? "MISSING_MAP" : placement==null
			? (stats==null || stats.matchString==null ? "MISSING_TRACE" : physicalEdge ? "PHYSICAL_SEQUENCE_EDGE" : "CLIPPED_QUERY")
			: rawBits==Float.NEGATIVE_INFINITY ? "IMPOSSIBLE_BAND" : rescueBelowThreshold && rawBits<threshold ? "BELOW_THRESHOLD" : "NONE";
		final CovarianceModelScoreOnly.Result repair=raw==null || rawBits==Float.NEGATIVE_INFINITY || rescueBelowThreshold && rawBits<threshold ? exact.align(sequence) : null;
		final long end=System.nanoTime();final float bits=repair==null ? rawBits : repair.score;
		require(Float.isFinite(bits) || bits==Float.NEGATIVE_INFINITY, "Invalid effective CM score cannot decide a candidate");
		return new Decision(bits, rawBits, Float.isFinite(bits) && bits>=threshold, sequence.length, reason,
			contacts, finalRadius, attempts, (composed-begin)*1e-9, raw==null ? 0 : (bandEnd-composed)*1e-9, repair==null ? 0 : (end-bandEnd)*1e-9, (end-begin)*1e-9,
			peakScore, peakTrace, repair==null ? 0 : repair.peakScoreCells);
	}
	static boolean clippedAtPhysicalEdge(AlignmentStats stats, int windowStart, int sequenceLength){
		require(windowStart>=0 && sequenceLength>0, "Physical-edge attribution uses the supplied oriented sequence, not a guessed original contig");
		if(!(stats instanceof EndClippedAligner.Result)){return false;}
		final EndClippedAligner.Result r=(EndClippedAligner.Result)stats;
		return r.qStart>0 && (long)windowStart+r.rStart==0
			|| r.qStop<r.qLen-1 && (long)windowStart+r.rStop==sequenceLength-1;
	}
	/** D32: compare successive band scores, never a saved exact reference.
	 * Equal impossible scores cannot terminate widening before radius128. */
	static boolean widen(int width, float score, float previous){
		require(width==32 || width==64 || width==128, "Adaptive placement uses the reviewed32/64/128 schedule");
		return width<128 && (!Float.isFinite(score) || !Float.isFinite(previous) || score!=previous);
	}
	public static final class Decision {
		Decision(float bits_, float raw_, boolean accepted_, int length_, String reason_, CovarianceModelPlacementBands.EdgeContacts contacts_, int radius_, int attempts_, double anchor_, double band_, double exact_, double total_,
				long bandCells_, long traceBytes_, long exactCells_){
			assert(length_>0 && total_>=0 && bandCells_>=0 && traceBytes_>=0 && exactCells_>=0):
				"Placement decision reports one completed candidate and stage-specific storage capacities";
			bits=bits_;rawBits=raw_;accepted=accepted_;length=length_;fallbackReason=reason_;contacts=contacts_;
			finalRadius=radius_;bandAttempts=attempts_;
			anchorSeconds=anchor_;bandSeconds=band_;exactSeconds=exact_;seconds=total_;
			bandScoreCells=bandCells_;traceBytes=traceBytes_;exactScoreCells=exactCells_;
		}
		public boolean usedFallback(){return !fallbackReason.equals("NONE");}
		public final float bits, rawBits;
		public final boolean accepted;
		public final int length;
		public final String fallbackReason;
		/** Radius0/attempts0 means no band was run; contacts belong to the final trace. */
		public final int finalRadius, bandAttempts;
		public final CovarianceModelPlacementBands.EdgeContacts contacts;
		public final double anchorSeconds, bandSeconds, exactSeconds, seconds;
		/** Stage capacities, not simultaneous JVM memory or minimum heap. */
		public final long bandScoreCells, traceBytes, exactScoreCells;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final CovarianceModel model;
	public final float threshold;
	public final int radius;
	public final boolean adaptive;
	public final boolean local;
	public final boolean extendClipped;
	public final boolean rescueBelowThreshold;
	private final CovarianceModelBandedCyk banded;
	private final CovarianceModelScoreOnly exact;
}
