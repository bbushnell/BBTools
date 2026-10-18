package idaligner;

import java.util.Arrays;
import structures.ByteBuilder;

/** Exact local alignment for the experimental consensus-end-clipping caller.
 * Both sequences have free ends; internal gaps and substitutions cost one,
 * matches score one, and N scores zero, as in Quantum's elementary scores.
 * Unlike Quantum this explores the full matrix without sparse pruning or its
 * position-dependent penalty. Thus a comparison changes search as well as ends.
 * Scratch belongs to one worker and is reused. Trace storage is one byte/cell.
 * @author Brian Bushnell, Raiden
 */
public final class EndClippedAligner {

	/** End-to-end model alignment, except model bases physically off a contig end.
	 * Reference flanks are free. A model prefix may start at reference position zero
	 * only when leftEdge is true; a model suffix may end at the last reference base
	 * only when rightEdge is true. There is no internal local restart.
	 * The caller must pass physical contig boundaries, not candidate-window edges. */
	public Result alignContigEdges(final byte[] query, final byte[] ref, final boolean leftEdge, final boolean rightEdge){
		return alignContigEdges(query, ref, leftEdge, rightEdge, 1, 1);
	}

	/** Experimental weighted edge objective; existing callers retain the 1/-1 overload.
	 * Substitutions and either gap direction pay penalty; N remains neutral.
	 * AlignmentStats.score retains its elementary definition. optimizationScore
	 * records this weighted objective, so identity and HBM consumers are unchanged. */
	public Result alignContigEdges(final byte[] query, final byte[] ref, final boolean leftEdge, final boolean rightEdge,
			final int matchReward, final int penalty){
		if(query==null || ref==null){throw new IllegalArgumentException("Edge alignment needs both sequences");}
		if(query.length==0 || ref.length==0){return null;}
		if(matchReward<1 || penalty<1 || ((long)query.length+ref.length)*Math.max(matchReward, penalty)>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Positive edge scores must fit the signed DP range");
		}
		final int rows=query.length, cols=ref.length+1;
		final long cells=(rows+1L)*cols;
		if(cols<1 || cells>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Edge traceback exceeds a Java array");}
		if(prev.length<cols){prev=new int[cols];curr=new int[cols];}
		if(trace.length<cells){trace=new byte[(int)cells];}
		Arrays.fill(prev, 0, cols, 0);// Free reference prefix with the full model still unconsumed.
		int best=Integer.MIN_VALUE, bestI=0, bestJ=0;
		for(int i=1; i<=rows; i++){
			curr[0]=leftEdge ? 0 : -i*penalty;
			trace[i*cols]=leftEdge ? (byte)0 : (byte)2;
			final int offset=i*cols;
			final byte q=query[i-1];
			for(int j=1; j<cols; j++){
				final byte r=ref[j-1];
				final int add=q=='N' || r=='N' ? 0 : q==r ? matchReward : -penalty;
				int score=prev[j-1]+add;
				byte direction=1;
				if(prev[j]-penalty>score){score=prev[j]-penalty;direction=2;}
				if(curr[j-1]-penalty>score){score=curr[j-1]-penalty;direction=3;}
				curr[j]=score;trace[offset+j]=direction;
				if((i==rows || (rightEdge && j==cols-1)) && score>best){best=score;bestI=i;bestJ=j;}
			}
			final int[] swap=prev;prev=curr;curr=swap;
		}
		if(best<=0){return null;}
		final Result result=new Result();
		result.qLen=rows;result.rLen=ref.length;
		result.qStop=bestI-1;result.rStop=bestJ-1;
		int i=bestI, j=bestJ;
		ops.clear();
		while(i>0){
			if(j==0 && leftEdge){break;}
			final byte direction=trace[i*cols+j];
			if(direction==1){
				final byte q=query[--i], r=ref[--j];
				ops.append((byte)(q=='N' || r=='N' ? 'N' : q==r ? 'm' : 'S'));
			}else if(direction==2){ops.append((byte)'I');i--;}
			else{assert(direction==3 && j>0):"Edge traceback must reach a permitted model/reference boundary";ops.append((byte)'D');j--;}
		}
		result.qStart=i;result.rStart=j;
		final byte[] match=new byte[ops.length()];
		for(int p=0; p<match.length; p++){match[p]=ops.array[match.length-1-p];}
		result.setFromMatchString(match);
		result.identity=(result.matches+0.5f*result.ns)/match.length;
		result.optimizationScore=best;
		result.matchReward=matchReward;result.nonmatchPenalty=penalty;
		result.physicalClips=(result.qStart>0 ? 1 : 0)|(result.qStop<rows-1 ? 2 : 0);
		assert(result.matches*matchReward-(result.subs+result.ins+result.dels)*penalty==best
			&& result.matches+result.subs+result.ins+result.ns==bestI-i
			&& result.matches+result.subs+result.dels+result.ns==bestJ-j)
			:"Edge traceback must conserve score and both consumed spans before caller coordinates are trusted";
		assert((result.qStart==0 || (leftEdge && result.rStart==0))
			&& (result.qStop==rows-1 || (rightEdge && result.rStop==ref.length-1)))
			:"Clipping a model end requires that alignment end to reach the physical contig boundary";
		return result;
	}

	/** Diagnostic boundary-local support filter. Walk inward from a genuinely
	 * clipped physical end; retain only the connected segment whose sliding
	 * windows meet the supplied identity floor. This never joins across a failed
	 * window. N has the same half credit as reported identity, not positive
	 * evidence hidden by excluding it from the denominator.
	 * Null means no supported anchor or equally scoring disconnected anchors;
	 * the latter are rejected as ambiguous, never resolved by strand orientation.
	 * No production caller uses this until its positive controls are reviewed. */
	public static Result trimEdgeSupport(final Result source, final boolean leftEdge, final boolean rightEdge,
			final int window, final float minimumIdentity){
		if(window<1 || !Float.isFinite(minimumIdentity) || minimumIdentity<=.5f || minimumIdentity>1){
			throw new IllegalArgumentException("Local edge support needs a positive window and identity floor in (.5,1]");
		}
		if(source==null){return null;}
		assert(source.matchString!=null && source.matchString.length>0)
			:"Boundary support consumes the existing verified operation trace, never reconstructed identity";
		final boolean left=(source.physicalClips&1)!=0 && leftEdge && source.rStart==0;
		final boolean right=(source.physicalClips&2)!=0 && rightEdge && source.rStop==source.rLen-1;
		if(!left && !right){return source;}// Full models and interior local rescues are unchanged.
		final int n=source.matchString.length;
		final int leftLength=left ? supportedLength(source.matchString, false, window, minimumIdentity) : 0;
		final int rightLength=right ? supportedLength(source.matchString, true, window, minimumIdentity) : 0;
		if(leftLength==n || rightLength==n){return source;}
		final Result a=leftLength==0 ? null : crop(source, 0, leftLength);
		final Result b=rightLength==0 ? null : crop(source, n-rightLength, n);
		if(a==null){return b;}if(b==null){return a;}
		// Two supported ends separated by a failed window cannot become one call.
		// Reject equal evidence rather than introducing an oriented-left preference.
		if(a.optimizationScore==b.optimizationScore){return null;}
		return a.optimizationScore>b.optimizationScore ? a : b;
	}
	private static int supportedLength(final byte[] match, final boolean reverse, final int window, final float floor){
		assert(match.length>0 && window>0):"An edge-connected support window must contain at least one trace column";
		final int n=match.length, width=Math.min(n,window);
		int twiceMatches=0;
		for(int k=0; k<width; k++){twiceMatches+=identityCredit(match[reverse ? n-1-k : k]);}
		if(twiceMatches/(2f*width)<floor){return 0;}
		int retained=width;
		for(int k=width; k<n; k++){
			twiceMatches+=identityCredit(match[reverse ? n-1-k : k]);
			twiceMatches-=identityCredit(match[reverse ? n-1-(k-width) : k-width]);
			if(twiceMatches/(2f*width)<floor){break;}
			retained=k+1;
		}
		return retained;
	}
	private static int identityCredit(final byte op){
		assert(op=='m' || op=='S' || op=='I' || op=='D' || op=='N')
			:"Local identity must use the same five operations as AlignmentStats.setFromMatchString";
		return op=='m' ? 2 : op=='N' ? 1 : 0;
	}
	private static Result crop(final Result source, final int from, final int to){
		assert(from>=0 && from<to && to<=source.matchString.length && (from==0 || to==source.matchString.length))
			:"Local support may remove only the interior-facing trace end while preserving a physical anchor";
		int prefixQ=0,prefixR=0;
		for(int k=0; k<from; k++){final byte op=source.matchString[k];prefixQ+=op=='D'?0:1;prefixR+=op=='I'?0:1;}
		final Result out=new Result();
		out.qLen=source.qLen;out.rLen=source.rLen;
		out.qStart=source.qStart+prefixQ;out.rStart=source.rStart+prefixR;
		out.setFromMatchString(Arrays.copyOfRange(source.matchString,from,to));
		final int q=out.matches+out.subs+out.ins+out.ns, r=out.matches+out.subs+out.dels+out.ns;
		if(q==0 || r==0){return null;}
		out.qStop=out.qStart+q-1;out.rStop=out.rStart+r-1;
		out.identity=(out.matches+.5f*out.ns)/(to-from);
		out.matchReward=source.matchReward;out.nonmatchPenalty=source.nonmatchPenalty;
		out.optimizationScore=out.matches*out.matchReward-(out.subs+out.ins+out.dels)*out.nonmatchPenalty;
		out.physicalClips=source.physicalClips&((from==0 ? 1 : 0)|(to==source.matchString.length ? 2 : 0));
		out.locallyTrimmed=true;
		assert(out.physicalClips!=0 && out.qStart>=source.qStart && out.qStop<=source.qStop
			&& out.rStart>=source.rStart && out.rStop<=source.rStop)
			:"Cropping must conserve both coordinate axes and the original physical truncation provenance";
		return out;
	}

	/** Highest positive local score; ties keep the first row-major endpoint and
	 * prefer diagonal, then insertion, then deletion. Null means no positive path.
	 * No minimum span is imposed here: callers must filter the returned alignment.
	 * Inputs are uppercase, strand-oriented; neither is swapped or modified. */
	public Result align(final byte[] query, final byte[] ref){
		if(query==null || ref==null){throw new IllegalArgumentException("Local alignment needs both sequences");}
		if(query.length==0 || ref.length==0){return null;}
		final int rows=query.length, cols=ref.length+1;
		final long cells=(rows+1L)*cols;
		if(cols<1 || cells>Integer.MAX_VALUE-8){
			throw new IllegalArgumentException("Local traceback exceeds a Java array: query="+rows+", reference="+ref.length);
		}
		if(prev.length<cols){prev=new int[cols];curr=new int[cols];}
		if(trace.length<cells){trace=new byte[(int)cells];}
		Arrays.fill(prev, 0, cols, 0);
		int best=0, bestI=0, bestJ=0;
		for(int i=1; i<=rows; i++){
			curr[0]=0;
			final int offset=i*cols;
			final byte q=query[i-1];
			for(int j=1; j<cols; j++){
				final byte r=ref[j-1];
				final int add=q=='N' || r=='N' ? 0 : q==r ? 1 : -1;
				int score=prev[j-1]+add;
				byte direction=1;
				if(prev[j]-1>score){score=prev[j]-1;direction=2;}
				if(curr[j-1]-1>score){score=curr[j-1]-1;direction=3;}
				if(score<=0){score=0;direction=0;}
				curr[j]=score;trace[offset+j]=direction;
				if(score>best){best=score;bestI=i;bestJ=j;}
			}
			final int[] swap=prev;prev=curr;curr=swap;
		}
		if(best==0){return null;}
		final Result result=new Result();
		result.qLen=rows;result.rLen=ref.length;
		result.qStop=bestI-1;result.rStop=bestJ-1;
		int i=bestI, j=bestJ;
		ops.clear();
		while(i>0 && j>0){
			final byte direction=trace[i*cols+j];
			if(direction==0){break;}
			if(direction==1){
				final byte q=query[--i], r=ref[--j];
				ops.append((byte)(q=='N' || r=='N' ? 'N' : q==r ? 'm' : 'S'));
			}else if(direction==2){ops.append((byte)'I');i--;}
			else{assert(direction==3) : "Trace direction must select one DP predecessor";ops.append((byte)'D');j--;}
		}
		result.qStart=i;result.rStart=j;
		final byte[] match=new byte[ops.length()];
		for(int p=0; p<match.length; p++){match[p]=ops.array[match.length-1-p];}
		result.setFromMatchString(match);
		result.optimizationScore=best;
		assert(result.score==best && result.matches+result.subs+result.ins+result.ns==bestI-i
			&& result.matches+result.subs+result.dels+result.ns==bestJ-j)
			: "Trace must reproduce local DP score and both consumed spans before caller coordinates are trusted";
		// Quantum's return identity gives a neutral N half credit; keep that
		// convention explicitly while excluding clipped bases from the denominator.
		result.identity=(result.matches+0.5f*result.ns)/match.length;
		return result;
	}

	/** Query/model bounds are inclusive and refer to the original full model.
	 * matchString covers only these bounds, never the clipped terminal sequence. */
	public static final class Result extends AlignmentStats {
		Result(){super(true);}
		public int qStart, qStop;
		public int optimizationScore;
		/** Physical model truncations before any optional interior-facing crop. */
		public int physicalClips;
		public boolean locallyTrimmed;
		public int matchReward=1, nonmatchPenalty=1;
	}

	private int[] prev=new int[0], curr=new int[0];
	private byte[] trace=new byte[0];
	private final ByteBuilder ops=new ByteBuilder();
}
