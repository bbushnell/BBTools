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
		if(query==null || ref==null){throw new IllegalArgumentException("Edge alignment needs both sequences");}
		if(query.length==0 || ref.length==0){return null;}
		final int rows=query.length, cols=ref.length+1;
		final long cells=(rows+1L)*cols;
		if(cols<1 || cells>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Edge traceback exceeds a Java array");}
		if(prev.length<cols){prev=new int[cols];curr=new int[cols];}
		if(trace.length<cells){trace=new byte[(int)cells];}
		Arrays.fill(prev, 0, cols, 0);// Free reference prefix with the full model still unconsumed.
		int best=Integer.MIN_VALUE, bestI=0, bestJ=0;
		for(int i=1; i<=rows; i++){
			curr[0]=leftEdge ? 0 : -i;
			trace[i*cols]=leftEdge ? (byte)0 : (byte)2;
			final int offset=i*cols;
			final byte q=query[i-1];
			for(int j=1; j<cols; j++){
				final byte r=ref[j-1];
				final int add=q=='N' || r=='N' ? 0 : q==r ? 1 : -1;
				int score=prev[j-1]+add;
				byte direction=1;
				if(prev[j]-1>score){score=prev[j]-1;direction=2;}
				if(curr[j-1]-1>score){score=curr[j-1]-1;direction=3;}
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
		assert(result.score==best && result.matches+result.subs+result.ins+result.ns==bestI-i
			&& result.matches+result.subs+result.dels+result.ns==bestJ-j)
			:"Edge traceback must conserve score and both consumed spans before caller coordinates are trusted";
		assert((result.qStart==0 || (leftEdge && result.rStart==0))
			&& (result.qStop==rows-1 || (rightEdge && result.rStop==ref.length-1)))
			:"Clipping a model end requires that alignment end to reach the physical contig boundary";
		return result;
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
	}

	private int[] prev=new int[0], curr=new int[0];
	private byte[] trace=new byte[0];
	private final ByteBuilder ops=new ByteBuilder();
}
