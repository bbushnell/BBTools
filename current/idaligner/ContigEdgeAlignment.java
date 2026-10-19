package idaligner;

import java.util.Arrays;

/** Bounded alternatives for physical-edge calls. Normal alignment operates on
 * a full model and a reference padded only beyond actual contig boundaries.
 * Verification statistics always describe observed genomic columns, never the
 * artificial padding. Interior refinement preserves the original physical end;
 * the padded-alignment alternative is evaluated separately before selection.
 * @author Raiden */
public final class ContigEdgeAlignment {
	public static EndClippedAligner.Result padded(final byte[] model, final byte[] ref,
			final boolean leftEdge, final boolean rightEdge, final int quantumCap, final boolean scrabble){
		if(model==null || ref==null || model.length==0 || ref.length==0 || quantumCap<1){
			throw new IllegalArgumentException("Padded alignment requires nonempty model/reference and a positive Quantum cap");
		}
		final int left=leftEdge?model.length:0, right=rightEdge?model.length:0;
		final long total=(long)left+ref.length+right;
		if(total>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Padded genomic window exceeds Java array capacity");}
		final byte[] padded=new byte[(int)total];Arrays.fill(padded,(byte)'N');
		System.arraycopy(ref,0,padded,left,ref.length);
		final AlignmentStats stats=new AlignmentStats(true);// No stale Quantum start hint.
		if(scrabble){ScrabbleAligner.alignAndTraceStatic(model,padded,stats);}
		else{QuantumAligner.alignAndTraceStatic(model,padded,stats,quantumCap);}
		return observed(stats,model.length,left,ref.length,leftEdge,rightEdge);
	}

	/** Remove only columns outside the real reference interval. Real genomic Ns
	 * stay in the trace and retain the caller's explicit half-identity credit. */
	public static EndClippedAligner.Result observed(final AlignmentStats raw, final int modelLength,
			final int paddingLeft, final int referenceLength, final boolean leftEdge, final boolean rightEdge){
		assert(raw!=null && modelLength>0 && paddingLeft>=0 && referenceLength>0)
			:"Virtual-to-genomic conversion requires the exact padded reference geometry";
		if(raw.matchString==null || raw.rStart<0 || raw.rStop<raw.rStart || raw.rStop>=raw.rLen){return null;}
		final long observedEnd=(long)paddingLeft+referenceLength;
		if(observedEnd>raw.rLen){throw new IllegalArgumentException("Observed interval exceeds padded reference");}
		int q=0,r=raw.rStart,from=-1,to=-1,qFrom=-1,rFrom=-1;
		for(int k=0; k<raw.matchString.length; k++){
			final byte op=raw.matchString[k];
			assert(op=='m' || op=='S' || op=='I' || op=='D' || op=='N')
				:"Padding removal must use the same five trace operations as AlignmentStats";
			final boolean observed=op=='I'
				? r>=paddingLeft && r<=observedEnd && !(leftEdge && r==paddingLeft) && !(rightEdge && r==observedEnd)
				: r>=paddingLeft && r<observedEnd;
			if(observed){
				if(from<0){from=k;qFrom=q;rFrom=r;}
				to=k+1;
			}
			if(op!='D'){q++;}if(op!='I'){r++;}
		}
		assert(q==modelLength && r==raw.rStop+1)
			:"Normal full-model traceback must conserve query and padded reference consumption before cropping";
		if(from<0){return null;}// Padding-only alignments contain no genomic evidence.
		final EndClippedAligner.Result out=new EndClippedAligner.Result();
		out.qLen=modelLength;out.rLen=referenceLength;
		out.qStart=qFrom;out.rStart=rFrom-paddingLeft;
		out.setFromMatchString(Arrays.copyOfRange(raw.matchString,from,to));
		final int qConsumed=out.matches+out.subs+out.ins+out.ns;
		final int rConsumed=out.matches+out.subs+out.dels+out.ns;
		if(qConsumed==0 || rConsumed==0){return null;}
		out.qStop=out.qStart+qConsumed-1;out.rStop=out.rStart+rConsumed-1;
		out.identity=(out.matches+.5f*out.ns)/out.matchString.length;
		out.optimizationScore=out.score;// Observed elementary score, not padded search score.
		out.physicalClips=(out.qStart>0 && leftEdge && out.rStart==0?1:0)
			|(out.qStop<modelLength-1 && rightEdge && out.rStop==referenceLength-1?2:0);
		assert(out.rStart>=0 && out.rStop<referenceLength && out.qStart>=0 && out.qStop<modelLength)
			:"Stripping virtual columns must return a span wholly within the original model and real genome";
		assert((out.qStart==0 || (out.physicalClips&1)!=0) && (out.qStop==modelLength-1 || (out.physicalClips&2)!=0))
			:"A model end removed by virtual padding must still be anchored at its actual physical contig side";
		return out;
	}

	/** Normal Scrabble endpoint re-anchoring on the observed model slice. The
	 * accepted identity/trace remain detection evidence; returned bounds are only
	 * refinement proposals. Physical sides stay fixed throughout interpretation.
	 * Reference is the original candidate window, so context is bounded by it. */
	public static Bounds refineInterior(final byte[] model, final byte[] ref,
			final EndClippedAligner.Result accepted, final int extension){
		return refineInterior(model,ref,accepted,extension,0);
	}
	/** Whole-contig form used by callers: accepted coordinates are relative to
	 * the candidate window at referenceOffset; proposals use full-reference coordinates. */
	public static Bounds refineInterior(final byte[] model, final byte[] ref,
			final EndClippedAligner.Result accepted, final int extension, final int referenceOffset){
		if(accepted==null){return null;}
		if(extension<0 || model.length!=accepted.qLen || referenceOffset<0 || (long)referenceOffset+accepted.rLen>ref.length){
			throw new IllegalArgumentException("Interior refinement must preserve the accepted model/reference frame");
		}
		assert(accepted.qStart>=0 && accepted.qStop>=accepted.qStart && accepted.qStop<model.length
			&& accepted.rStart>=0 && accepted.rStop>=accepted.rStart && accepted.rStop<accepted.rLen)
			:"Refinement must slice only coordinates retained by the accepted alignment in its original frames";
		final int mask=accepted.physicalClips;
		final int originalStart=referenceOffset+accepted.rStart,originalStop=referenceOffset+accepted.rStop;
		assert(((mask&1)==0 || originalStart==0) && ((mask&2)==0 || originalStop==ref.length-1))
			:"A candidate-window boundary cannot substitute for the original physical contig boundary during refinement";
		if(mask==0 || mask==3){return new Bounds(originalStart,originalStop,false);}
		final int from=Math.max(0,originalStart-extension),to=Math.min(ref.length,originalStop+extension+1);
		final byte[] query=Arrays.copyOfRange(model,accepted.qStart,accepted.qStop+1);
		final byte[] reference=Arrays.copyOfRange(ref,from,to);
		shared.Tools.toUpperCase(reference);// Match the caller's private uppercase trim window.
		final AlignmentStats stats=new AlignmentStats(true);
		// Keep model/query orientation even when indels make the query longer.
		ScrabbleAligner.alignAndTraceStatic(query,reference,stats);
		if(stats.matchString==null || stats.rStart<0 || stats.rStop<stats.rStart || stats.rStop>=reference.length){
			return new Bounds(originalStart,originalStop,false);
		}
		final int start=(mask&1)!=0?originalStart:from+stats.rStart;
		final int stop=(mask&2)!=0?originalStop:from+stats.rStop;
		if(start>stop){return new Bounds(originalStart,originalStop,false);}
		assert(((mask&1)==0 || start==originalStart) && ((mask&2)==0 || stop==originalStop))
			:"Only the interior-facing endpoint may move; the original physical endpoint is exempt";
		return new Bounds(start,stop,start!=originalStart || stop!=originalStop);
	}
	public static final class Bounds{
		Bounds(int start_,int stop_,boolean changed_){start=start_;stop=stop_;changed=changed_;}
		public final int start,stop;public final boolean changed;
	}
	private ContigEdgeAlignment(){}
}
