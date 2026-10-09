package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import align2.MultiStateAligner9PacBio;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import stream.Read;
import stream.SamLine;

/** Worker-local adapter for the unchanged PacBio MSA. Experimental 18S only.
 * Native endpoints and padding are evidence, never instructions to expand or clamp.
 * @author Raiden
 */
final class Euk18sPacBioAligner {
	Euk18sPacBioAligner(){this(align2.PacBioScoreParameters.DEFAULT);}
	Euk18sPacBioAligner(align2.PacBioScoreParameters costs_){this(costs_, false);}
	Euk18sPacBioAligner(align2.PacBioScoreParameters costs_, boolean rolling_){
		check(costs_!=null,"Each worker must bind immutable PacBio costs before its first allocation");costs=costs_;
		rolling=rolling_;
	}

	Result align(final byte[] query, final byte[] ref, final boolean compareQuantum){
		check(query!=null && ref!=null && query.length>0 && ref.length>0,
			"PacBio fillUnlimited requires nonempty query and reference");
		check(SamLine.INTRON_LIMIT==Integer.MAX_VALUE,
			"D37 identity discounts all deletion runs by sqrt; a changed global intron limit would silently omit some");
		final long begin=System.nanoTime(), cpuBegin=cpuTime();
		long allocationNanos=0;
		final int[] score;
		final byte[] trace;
		if(rolling){
			if(rollingMsa==null){rollingMsa=new align2.PacBioRollingAligner(costs);}
			rollingMsa.fillUnlimited(query, ref, 0, ref.length-1);
			score=rollingMsa.score();trace=rollingMsa.traceback();
			allocationNanos=rollingMsa.allocationNanos();
			maxQuery=Math.max(query.length, maxQuery);maxRef=Math.max(ref.length, maxRef);
		}else{
			if(msa==null || query.length>maxQuery || ref.length>maxRef){
				final long allocationBegin=System.nanoTime();
				maxQuery=Math.max(query.length, maxQuery);maxRef=Math.max(ref.length, maxRef);
				msa=null;// Permit collection of the old matrix before a larger allocation.
				msa=new MultiStateAligner9PacBio(maxQuery,maxRef,costs);
				allocationNanos=System.nanoTime()-allocationBegin;
			}
			final int[] end=msa.fillUnlimited(query, ref, 0, ref.length-1, null);
			check(end!=null && end.length>=3, "Unlimited PacBio fill must supply the endpoint consumed by score and traceback");
			score=msa.score(query, ref, 0, ref.length-1, end[0], end[1], end[2], false);
			trace=msa.traceback(query, ref, 0, ref.length-1, end[0], end[1], end[2], false);
		}
		final Result r=new Result(trace, score, query.length, ref.length);
		r.maxQuery=maxQuery;r.maxRef=maxRef;r.worker=Thread.currentThread().getId();
		r.workingMatrixBytes=rolling ? rollingMsa.workingBytes() : matrixBytes(query.length, ref.length);
		r.retainedMatrixBytes=rolling ? rollingMsa.retainedBytes() : matrixBytes(maxQuery, maxRef);
		r.allocationNanos=allocationNanos;r.primaryCpuNanos=elapsedCpu(cpuBegin);
		r.primaryNanos=System.nanoTime()-begin;
		if(compareQuantum){
			final long qBegin=System.nanoTime(), qCpuBegin=cpuTime();
			final AlignmentStats q=new AlignmentStats(true);
			r.quantumIdentity=QuantumAligner.alignAndTraceStatic(query, ref, q);
			r.quantumStart=q.rStart;r.quantumStop=q.rStop;
			r.quantumCpuNanos=elapsedCpu(qCpuBegin);r.quantumNanos=System.nanoTime()-qBegin;
			r.comparedQuantum=true;
		}
		return r;
	}

	/** Primitive payload of packed[3][maxRows+1][maxColumns+2] in the native constructor.
	 * Excludes array headers, limits, trace scratch, GC-retained old arrays and all other heap. */
	static long matrixBytes(final int rows, final int columns){
		assert(rows>0 && columns>0) : "Matrix payload describes a real nonempty alignment allocation";
		return 12L*(rows+1L)*(columns+2L);
	}
	private static long cpuTime(){
		return CPU.isCurrentThreadCpuTimeSupported() && CPU.isThreadCpuTimeEnabled() ? CPU.getCurrentThreadCpuTime() : -1;
	}
	private static long elapsedCpu(final long before){
		final long after=cpuTime();return before<0 || after<0 ? -1 : after-before;
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new IllegalStateException(why);}}

	/** Expanded native trace. X consumes virtual left reference; Y consumes query
	 * only (MultiStateAligner9PacBio.traceback2/score2). Both count as identity errors. */
	static final class Result extends AlignmentStats {
		Result(final byte[] trace, final int[] nativeScore, final int queryLength, final int referenceLength){
			super(true);
			check(nativeScore!=null && (nativeScore.length==6 || nativeScore.length==8),
				"PacBio score2 returns six values, or eight with padding; preserve that contract");
			check(trace!=null && trace.length>0, "Native traceback is required for D37 gap-aware identity");
			matchString=trace;score=nativeScore[0];rStart=nativeScore[1];rStop=nativeScore[2];
			qLen=queryLength;rLen=referenceLength;
			padLeft=nativeScore.length==8 ? nativeScore[6] : 0;
			padRight=nativeScore.length==8 ? nativeScore[7] : 0;
			int q=0, r=rStart;
			for(final byte op:trace){
				switch(op){
					case 'm': matches++;q++;r++;break;
					case 'S': subs++;q++;r++;break;
					case 'N': ns++;q++;r++;break;
					case 'D': dels++;r++;break;
					case 'I': ins++;q++;break;
					case 'X': leftOverhang++;q++;r++;break;
					case 'Y': rightOverhang++;q++;break;
					default: throw new IllegalStateException("Unsupported expanded PacBio operation: "+(char)op);
				}
			}
			check(q==qLen && r-1==rStop, "Native trace must agree with full query and score2 reference endpoints: q="
				+q+" expected="+qLen+" stop="+(r-1)+" expected="+rStop);
			check(trace[trace.length-1]!='D', "Read.identitySkewed does not support a terminal deletion; native winner requires investigation");
			identity=Read.identitySkewed(trace, true, true, false, false);
			flatIdentity=Read.identityFlat(trace, true);
			validBounds=rStart>=0 && rStop>=rStart && rStop<rLen;
		}
		String outcome(final float cutoff, final int maximumLength){
			assert(Float.isFinite(cutoff) && cutoff>=0 && cutoff<=1) : "Diagnostic outcome must use the caller's actual identity gate";
			return !validBounds ? "REJECT_NATIVE_BOUNDS" : rStop-rStart+1>maximumLength ? "REJECT_MAXLEN"
				: identity<cutoff ? "REJECT_IDENTITY" : "PRIMARY_ELIGIBLE";
		}
		final int padLeft, padRight;
		final boolean validBounds;
		final float flatIdentity;
		int leftOverhang, rightOverhang, maxQuery, maxRef;
		long worker, workingMatrixBytes, retainedMatrixBytes, allocationNanos, primaryNanos, primaryCpuNanos;
		boolean comparedQuantum;
		float quantumIdentity=Float.NaN;
		int quantumStart=-1, quantumStop=-1;
		long quantumNanos=0, quantumCpuNanos=-1;
	}
	private MultiStateAligner9PacBio msa;
	private align2.PacBioRollingAligner rollingMsa;
	final boolean rolling;
	final align2.PacBioScoreParameters costs;
	private int maxQuery, maxRef;
	private static final ThreadMXBean CPU=ManagementFactory.getThreadMXBean();
}
