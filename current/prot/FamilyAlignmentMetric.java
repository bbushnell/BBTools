package prot;

import prot.cluster.TailAlignmentDecision;
import prot.cluster.TailAlignmentGate;

/** Explicit geometry and normalization contract for a calibrated family package.
 * A descriptor alone does not select new gates or enable a production mode.
 * Workers and observations are single-threaded scratch; the descriptor is immutable.
 * @author Brian, Keqing
 */
public final class FamilyAlignmentMetric{

	/** Historical packages use full-query glocal scores without a separate span gate. */
	public static final FamilyAlignmentMetric LEGACY=new FamilyAlignmentMetric(Mode.GLOCAL_FULL_QUERY,0);

	public FamilyAlignmentMetric(Mode mode_,int minimumSpan_){
		this(mode_,minimumSpan_,ProteinSubstitutionMatrix.BLOSUM62);
	}
	/** Custom matrices are supported by the pinned family-evaluation geometries.
	 * Historical glocal packages remain explicitly BLOSUM62, never silently rebound.
	 */
	public FamilyAlignmentMetric(Mode mode_,int minimumSpan_,ProteinSubstitutionMatrix matrix_){
		this(mode_,minimumSpan_,matrix_,0);
	}
	/** Paired columns consume a residue on both sides; gaps on either side never count. */
	public FamilyAlignmentMetric(Mode mode_,int minimumSpan_,ProteinSubstitutionMatrix matrix_,int minimumPaired_){
		require(minimumPaired_>=0 && minimumPaired_<=1000 && (mode_!=Mode.GLOCAL_FULL_QUERY || minimumPaired_==0),
			"Paired minimum must be in0..1000 and requires a pinned metric");
		minimumPaired=minimumPaired_;
		require(mode_!=null,"A calibrated alignment must name its geometry and denominator");
		require(matrix_!=null,"A calibrated alignment must bind its substitution matrix");
		require(mode_!=Mode.GLOCAL_FULL_QUERY || matrix_.same(ProteinSubstitutionMatrix.BLOSUM62),
			"Custom matrices require a pinned family metric; legacy glocal packages are BLOSUM62-only");
		require(mode_==Mode.GLOCAL_FULL_QUERY ? minimumSpan_==0 : minimumSpan_>0,
			"Legacy alignment has span0; a pinned metric requires an explicit positive span on both sequences");
		mode=mode_;minimumSpan=minimumSpan_;matrix=matrix_;
	}

	public Worker newWorker(){return new Worker(this,null);}
	public Worker newWorker(FamilyMatchInformation information){return new Worker(this,information);}
	public boolean same(FamilyAlignmentMetric other){return other!=null && mode==other.mode && minimumSpan==other.minimumSpan && minimumPaired==other.minimumPaired && matrix.same(other.matrix);}
	public String name(){return mode.label;}

	/** One mutable scorer per consuming thread; never shares observations or active consensuses. */
	public static final class Worker{
		private Worker(FamilyAlignmentMetric metric_,FamilyMatchInformation information_){
			require(information_==null || metric_.mode!=Mode.GLOCAL_FULL_QUERY,
				"Information measurement requires a pinned recorded traceback, not legacy aggregate metrics");
			metric=metric_;information=information_;consensus=new HbmPositionModel.ConsensusScorer(metric.matrix);
		}

		/** Measures without applying family acceptance cutoffs. Only MEASURED has usable normalized
		 * metrics. Other states are explicit failures, never missing values disguised as NaN.
		 * fullQuerySelf is the historical BLOSUM62 cache. Custom workers recompute
		 * their matrix's full-query self score, so a shared legacy query cache cannot
		 * contaminate custom normalization or reported provenance.
		 */
		public void measure(ProteinSequence query,ProteinSequence reference,int fullQuerySelf,Observation out){
			measure(query,reference,fullQuerySelf,out,null);
		}
		/** Optional development counts share the exact existing path and scratch lifetime. */
		public void measure(ProteinSequence query,ProteinSequence reference,int fullQuerySelf,Observation out,
				FamilyAlignmentCounts counts){
			if(counts!=null){counts.clear();}
			require(counts==null || metric.mode!=Mode.GLOCAL_FULL_QUERY,
				"Residue counts require a pinned recorded traceback; legacy aggregate metrics are insufficient");
			require(query!=null && reference!=null && query.length()>0 && reference.length()>0 && out!=null,
				"Family alignment requires nonempty encoded sequences and caller-owned output");
			require(fullQuerySelf>0,"Full-query self score must be positive, as required by the historical metric");
			final int selectedSelf="blosum62".equals(metric.matrix.name) ? fullQuerySelf : metric.matrix.selfScore(query.enc);
			require(selectedSelf>0,"Selected matrix requires positive full-query self normalization");
			out.clear();out.fullQuerySelf=selectedSelf;
			if(metric.mode==Mode.GLOCAL_FULL_QUERY){
				final AAAlignment a=GlocalAminoLinear.align(query.enc,reference.enc,false);
				out.rawScore=a.rawScore;out.denominatorSelf=fullQuerySelf;
				out.queryStart=a.qStart;out.queryEnd=a.qStop;out.referenceStart=a.tStart;out.referenceEnd=a.tStop;
				out.columns=a.length;out.identities=a.identities;out.pairedColumns=a.identities+a.mismatches;
				out.identityPct=a.pident();out.blosum=a.rawScore/(double)fullQuerySelf;
				out.state=State.MEASURED;
			}else if(metric.mode==Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY ||
					metric.mode==Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY){
				if(query.length()<metric.minimumSpan || reference.length()<metric.minimumSpan){
					out.state=State.INPUT_TOO_SHORT;return;
				}
				final TailAlignmentGate activeGate=gate(Math.max(query.length(),reference.length()));
				final TailAlignmentDecision a=activeGate.decide(query.enc,reference.enc,consensus.resetConsensus(reference.enc),counts);
				copy(a,out);
			}else if(metric.mode==Mode.PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY){
				if(query.length()<metric.minimumSpan || reference.length()<metric.minimumSpan){
					out.state=State.INPUT_TOO_SHORT;return;
				}
				final TailAlignmentGate activeGate=gate(Math.max(query.length(),reference.length()));
				final byte[] queryRev=reverse(query.enc,true), referenceRev=reverse(reference.enc,false);
				final TailAlignmentDecision a=activeGate.decide(queryRev,referenceRev,consensus.resetConsensus(referenceRev),counts);
				copy(a,out);
				if(a.reason!=TailAlignmentDecision.Reason.NO_ALIGNMENT){
					out.queryStart=query.length()-1-a.queryEnd;out.queryEnd=query.length()-1-a.queryStart;
					out.referenceStart=reference.length()-1-a.referenceEnd;out.referenceEnd=reference.length()-1-a.referenceStart;
				}
			}else{
				throw new AssertionError("Unhandled family metric: "+metric.mode);
			}
			require(Double.isFinite(out.identityPct) && out.identityPct>=0 && out.identityPct<=100 && Double.isFinite(out.blosum),
				"Measured family metrics must be finite and identity must remain in0..100");
		}

		/** Diagnostic-only baseline score on the same pair/matrix; never an acceptance value.
		 * Invoke only after a real overhang measurement, with inputs meeting its span floor. */
		public int maximalTrimScore(ProteinSequence query,ProteinSequence reference){
			require(metric.mode==Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY &&
				query.length()>=metric.minimumSpan && reference.length()>=metric.minimumSpan,
				"Paired trimming diagnostics require span-eligible overhang-mode inputs");
			final HbmPositionModel.TailResult a=consensus.resetConsensus(reference.enc)
				.alignTail(query.enc,Math.max(query.length(),reference.length()),false,true);
			return a==null ? 0 : a.score/HbmPositionModel.SCALE;
		}

		private TailAlignmentGate gate(final int inputLength){
			final int needed=Math.max(inputLength,metric.minimumPaired);
			if(gate==null || gate.maxLength<needed){
				// The capacity only grows, and always covers both full inputs: no hidden tail cropping.
				final boolean overhangOnly=metric.mode==Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY;
				gate=new TailAlignmentGate(overhangOnly ? TailAlignmentGate.Mode.PIN_RIGHT_OVERHANG_ONLY :
					TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT,TailAlignmentGate.Denominator.ALIGNED_QUERY,
					needed,metric.minimumSpan,0,0,information,metric.matrix,metric.minimumPaired);
			}
			return gate;
		}

		private byte[] reverse(final byte[] input,final boolean query){
			byte[] out=query ? reversedQuery : reversedReference;
			if(out==null || out.length!=input.length){
				out=new byte[input.length];
				if(query){reversedQuery=out;}else{reversedReference=out;}
			}
			for(int i=0,j=input.length-1;i<input.length;i++,j--){out[i]=input[j];}
			return out;
		}

		private static void copy(final TailAlignmentDecision a,final Observation out){
			if(a.reason==TailAlignmentDecision.Reason.NO_ALIGNMENT){out.state=State.NO_ALIGNMENT;return;}
			require(a.score%HbmPositionModel.SCALE==0,
				"Direct substitution emissions and gap costs must have integral matrix-unit scores");
			out.rawScore=a.score/HbmPositionModel.SCALE;
			out.queryStart=a.queryStart;out.queryEnd=a.queryEnd;
			out.referenceStart=a.referenceStart;out.referenceEnd=a.referenceEnd;
			if(a.reason==TailAlignmentDecision.Reason.QUERY_SPAN_SHORT || a.reason==TailAlignmentDecision.Reason.REFERENCE_SPAN_SHORT ||
					a.reason==TailAlignmentDecision.Reason.PAIRED_SPAN_SHORT){
				out.columns=a.columns;out.identities=a.identities;out.pairedColumns=a.paired;
				out.state=State.SPAN_TOO_SHORT;return;
			}
			require(a.accepted() || a.reason==TailAlignmentDecision.Reason.NONPOSITIVE_SCORE,
				"Measurement-only probe gates must retain every signed span-qualified alignment");
			require(a.denominatorSelfScore>0 && a.denominatorSelfScore%HbmPositionModel.SCALE==0,
				"Aligned-query normalization requires a positive integral selected-matrix self score");
			out.denominatorSelf=a.denominatorSelfScore/HbmPositionModel.SCALE;
			out.columns=a.columns;out.identities=a.identities;out.pairedColumns=a.paired;
			// Use the same percent-identity arithmetic as AAAlignment.pident, including gap columns.
			out.identityPct=100.0*a.identities/a.columns;
			out.blosum=out.rawScore/(double)out.denominatorSelf;out.state=State.MEASURED;
			out.information=a.information;out.informationMeasured=a.informationMeasured;
		}

		private final FamilyAlignmentMetric metric;
		private final FamilyMatchInformation information;
		private final HbmPositionModel.ConsensusScorer consensus;
		private byte[] reversedQuery,reversedReference;
		private TailAlignmentGate gate;
	}

	/** Reusable result. Numeric zeros in non-MEASURED states are inert, not calibrated scores.
	 * A span-rejected alignment retains its signed raw score and coordinates for diagnosis.
	 */
	public static final class Observation{
		public boolean measured(){return state==State.MEASURED;}
		public double identityPct(){
			if(!measured()){throw new IllegalArgumentException("No calibrated identity is available for alignment state "+state);}
			return identityPct;
		}
		public double blosum(){
			if(!measured()){throw new IllegalArgumentException("No calibrated BLOSUM ratio is available for alignment state "+state);}
			return blosum;
		}
		private void clear(){
			state=State.NO_ALIGNMENT;rawScore=fullQuerySelf=denominatorSelf=columns=identities=pairedColumns=0;
			queryStart=queryEnd=referenceStart=referenceEnd=-1;identityPct=blosum=0;
			information=0;informationMeasured=false;
		}
		public State state=State.NO_ALIGNMENT;
		public double information;
		public boolean informationMeasured;
		public int rawScore,fullQuerySelf,denominatorSelf,columns,identities;
		/** Residue-residue traceback columns; full-query span alone includes gaps. */
		public int pairedColumns;
		public int queryStart=-1,queryEnd=-1,referenceStart=-1,referenceEnd=-1;
		private double identityPct,blosum;
	}

	public enum State{MEASURED,INPUT_TOO_SHORT,NO_ALIGNMENT,SPAN_TOO_SHORT}
	public enum Mode{
		GLOCAL_FULL_QUERY("glocal_full_query"),
		PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY("pin_right_trim_left_aligned_query"),
		PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY("pin_left_trim_right_aligned_query"),
		PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY("pin_right_overhang_only_aligned_query");
		Mode(String label_){label=label_;}
		public final String label;
	}
	public final Mode mode;
	public final int minimumSpan,minimumPaired;
	public final ProteinSubstitutionMatrix matrix;
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
}
