package prot;

/** Stop-facing evidence measurements, separate from full-protein family assignment.
 * Thresholds are explicit experimental settings, not calibrated probabilities.
 * The caller owns mutable Metrics; profiles and Policy may be shared by workers.
 * @author Keqing
 */
public final class HbmTailGate{

	/** Replays an original-coordinate tail trace and returns its acceptance class.
	 * Only paired columns support coverage. Consumed terminal gap columns do not.
	 * The core denominator is its intersection with the capped reference window.
	 * This layer checks sequence/geometry; HbmCompetitiveGate adds the profile
	 * policy's terminal-score condition using its bound positional model.
	 */
	public static int evaluate(final byte[] query, final byte[] reference,
			final HbmPositionModel.TailResult alignment, final int coreStart, final int coreEnd,
			final Policy policy, final Metrics out){
		require(query!=null && query.length>0 && reference!=null && reference.length>0
			&& policy!=null && out!=null, "Tail evidence requires query, reference, policy and worker metrics");
		require(coreStart>=0 && coreEnd>=coreStart && coreEnd<reference.length, "Tail family core lies outside its reference");
		out.clear();
		if(alignment==null){return out.status=NO_ALIGNMENT;}
		require(alignment.path!=null && alignment.path.length>0 && alignment.score>0,
			"Tail evidence requires a recorded positive-scoring alignment");
		require(alignment.queryStart>=Math.max(0, query.length-policy.maxLength)
			&& alignment.queryStart<query.length && alignment.queryEnd==query.length-1
			&& alignment.start>=Math.max(0, reference.length-policy.maxLength)
			&& alignment.start<reference.length && alignment.end==reference.length-1,
			"Tail trace must lie inside capped windows and consume both stop-facing ends");
		final int cappedCoreStart=Math.max(coreStart, reference.length-policy.maxLength);
		out.coreLength=Math.max(0, coreEnd-cappedCoreStart+1);
		out.queryStart=alignment.queryStart; out.referenceStart=alignment.start;
		int q=alignment.queryStart, r=alignment.start, identities=0;
		int lastQ=-1, lastR=-1;
		for(byte op : alignment.path){
			if(op=='m'){
				require(q<query.length && r<reference.length, "Paired tail column exceeds sequence bounds");
				out.paired++;
				if(query[q]==reference[r] && query[q]>=0 && query[q]<20){identities++;}
				if(r>=cappedCoreStart && r<=coreEnd){out.pairedCore++;}
				//TODO: Paired terminal columns include mismatches; they do not establish
				//conservation. The 20261018 robustness assay accepts 199/200 consensus
				//queries with every final20 residue replaced. Calibrate terminal signal
				//against real supported and unsupported ends before production use.
				if(q>=query.length-policy.terminalWindow && r>=reference.length-policy.terminalWindow){out.terminalPaired++;}
				if(q>=query.length-policy.terminalWindow){out.queryTerminalPaired++;}
				lastQ=q++; lastR=r++;
			}else if(op=='I'){q++;}
			else if(op=='D'){r++;}
			else{throw new IllegalArgumentException("Unknown tail path operation: "+(char)op);}
			require(q<=query.length && r<=reference.length, "Tail path overruns sequence bounds");
		}
		require(q==query.length && r==reference.length && out.paired>0,
			"Tail evidence must pair residues and consume both declared terminal coordinates");
		out.querySuffix=query.length-1-lastQ; out.referenceSuffix=reference.length-1-lastR;
		out.identity=(float)(100.0*identities/(policy.profileEvidence ? out.paired : alignment.path.length));
		out.pairedFraction=(float)(out.paired/(double)Math.max(query.length-alignment.queryStart, reference.length-alignment.start));
		out.coreCoverage=out.coreLength==0 ? 0 : (float)(out.pairedCore/(double)out.coreLength);
		assert(out.pairedCore<=out.coreLength && out.pairedFraction<=1 && out.coreCoverage<=1)
			: "Paired tail coverage cannot exceed the available query/reference/core columns";
		if(out.querySuffix>policy.maxTerminalGap || out.referenceSuffix>policy.maxTerminalGap
			|| out.terminalPaired<policy.minTerminalPaired
			|| (policy.profileEvidence && out.queryTerminalPaired<policy.minTerminalPaired)){return out.status=UNSUPPORTED_END;}
		if(out.paired<policy.minPaired || out.pairedCore<policy.minPaired){return out.status=SHORT_SUPPORT;}
		if(out.identity<policy.minIdentity){return out.status=LOW_IDENTITY;}
		if(out.pairedFraction<policy.minCoverage || out.coreCoverage<policy.minCoverage){return out.status=LOW_COVERAGE;}
		return out.status=SUPPORTED;
	}

	/** Immutable explicit policy; defaults are a conservative starting experiment. */
	public static final class Policy{
		public Policy(int maxLength_, int minPaired_, int maxTerminalGap_, int terminalWindow_,
				int minTerminalPaired_, float minIdentity_, float minCoverage_){
			this(maxLength_, minPaired_, maxTerminalGap_, terminalWindow_, minTerminalPaired_, minIdentity_, minCoverage_, false);
		}
		private Policy(int maxLength_, int minPaired_, int maxTerminalGap_, int terminalWindow_,
				int minTerminalPaired_, float minIdentity_, float minCoverage_, boolean profileEvidence_){
			require(maxLength_>0 && minPaired_>0 && minPaired_<=maxLength_
				&& maxTerminalGap_>=0 && maxTerminalGap_<=maxLength_
				&& terminalWindow_>0 && terminalWindow_<=maxLength_
				&& minTerminalPaired_>0 && minTerminalPaired_<=terminalWindow_
				&& Float.isFinite(minIdentity_) && minIdentity_>=0 && minIdentity_<=100
				&& Float.isFinite(minCoverage_) && minCoverage_>0 && minCoverage_<=1,
				"Invalid tail length, terminal-support, identity or coverage thresholds");
			maxLength=maxLength_; minPaired=minPaired_; maxTerminalGap=maxTerminalGap_;
			terminalWindow=terminalWindow_; minTerminalPaired=minTerminalPaired_;
			minIdentity=minIdentity_; minCoverage=minCoverage_;
			profileEvidence=profileEvidence_;
		}
		public final int maxLength, minPaired, maxTerminalGap, terminalWindow, minTerminalPaired;
		public final float minIdentity, minCoverage;
		/** Paired identity plus positive terminal profile score; joint terminal pairing is still required. */
		public final boolean profileEvidence;
	}

	/** Counts refer to this one alignment, never accumulated across family candidates. */
	public static final class Metrics{
		public void copyFrom(Metrics other){
			assert(other!=null) : "Winner measurements must come from evaluated tail evidence";
			status=other.status; paired=other.paired; pairedCore=other.pairedCore; terminalPaired=other.terminalPaired;
			queryTerminalPaired=other.queryTerminalPaired; terminalScore64=other.terminalScore64;
			coreLength=other.coreLength; queryStart=other.queryStart; referenceStart=other.referenceStart;
			querySuffix=other.querySuffix; referenceSuffix=other.referenceSuffix;
			identity=other.identity; pairedFraction=other.pairedFraction; coreCoverage=other.coreCoverage;
		}
		void clear(){
			status=NO_ALIGNMENT; paired=pairedCore=terminalPaired=coreLength=0;
			queryTerminalPaired=0; terminalScore64=Integer.MIN_VALUE;
			queryStart=referenceStart=querySuffix=referenceSuffix=-1;
			identity=pairedFraction=coreCoverage=Float.NaN;
		}
		public int status, paired, pairedCore, terminalPaired, coreLength;
		public int queryTerminalPaired, terminalScore64;
		public int queryStart, referenceStart, querySuffix, referenceSuffix;
		public float identity, pairedFraction, coreCoverage;
	}
	public static String statusName(final int status){
		if(status<0 || status>=NAMES.length){throw new IllegalArgumentException("Unknown tail evidence class: "+status);}
		return NAMES[status];
	}
	private static void require(boolean value, String message){if(!value){throw new IllegalArgumentException(message);}}
	private HbmTailGate(){}
	public static final int NO_ALIGNMENT=0, SUPPORTED=1, UNSUPPORTED_END=2, SHORT_SUPPORT=3, LOW_IDENTITY=4, LOW_COVERAGE=5;
	public static final Policy EXPERIMENTAL=new Policy(500, 30, 10, 20, 10, 40.612846f, 0.8f);
	public static final Policy PROFILE_JOINT=new Policy(500, 30, 10, 20, 10, 40.612846f, 0.8f, true);
	public static Policy forName(String name){
		if("geometry".equalsIgnoreCase(name)){return EXPERIMENTAL;}
		if("profile".equalsIgnoreCase(name)){return PROFILE_JOINT;}
		throw new IllegalArgumentException("hbmtailrule must be geometry or profile: "+name);
	}
	/** Change only the window cap; retain every named gate threshold and profile rule. */
	public static Policy forName(String name, int maxLength){
		final Policy base=forName(name);
		return maxLength==base.maxLength ? base : new Policy(maxLength, base.minPaired, base.maxTerminalGap,
			base.terminalWindow, base.minTerminalPaired, base.minIdentity, base.minCoverage, base.profileEvidence);
	}
	private static final String[] NAMES={"no_alignment", "stop_supported", "unsupported_end", "short_support", "low_identity", "low_coverage"};
	public static final int STATUS_COUNT=NAMES.length;
}
