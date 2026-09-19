package align2;

import java.util.Arrays;

import idaligner.IDAlignerStatics;

/**
 * Reusable int implementation of QuantumAligner's whole-top sparse ranking.
 * Explores competing high-scoring ridges and can bridge unexplored deletions.
 * This is a ranking prealigner, not BBMap's authoritative variable-gap MSA.
 *
 * @author Brian Bushnell, Collei
 * @date September 18, 2026
 */
public class QuantumRanker {

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public QuantumRanker(){this(512);}

	public QuantumRanker(int initialColumns){
		assert(initialColumns>0) : "Scratch capacity must be positive; requested "+initialColumns;
		ensureCapacity(initialColumns);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Aligns a query glocally to an inclusive reference window. Result is reused. */
	public Result align(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart){
		return align(query, ref, refStart, refStop, expectedStart, false, false);
	}

	/** Scalar/SIMD-selectable entry point for qualification and mapper dispatch. */
	public Result align(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart, final boolean useSIMD){
		return align(query, ref, refStart, refStop, expectedStart, useSIMD, false,
				Integer.MAX_VALUE);
	}

	/** Score-only alignment constrained to paths with at most maxEdits edits. */
	public Result align(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart, final boolean useSIMD,
			final int maxEdits){
		return align(query, ref, refStart, refStop, expectedStart, useSIMD, false,
				maxEdits);
	}

	/** Optional compact sparse traceback; score-only calls record no ancestry. */
	public Result align(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart, final boolean useSIMD,
			final boolean traceback){
		return align(query, ref, refStart, refStop, expectedStart, useSIMD, traceback,
				Integer.MAX_VALUE);
	}

	/** Optional traceback plus a hard per-path edit budget. */
	public Result align(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart, final boolean useSIMD,
			final boolean traceback, final int maxEdits){
		return alignInternal(query, ref, refStart, refStop, expectedStart, useSIMD, traceback,
				maxEdits);
	}

	private Result alignInternal(final byte[] query, final byte[] ref, final int refStart,
			final int refStop, final int expectedStart, final boolean useSIMD,
			final boolean traceback, final int maxEdits){
		result.clear();
		assert(query!=null && ref!=null) : "QuantumRanker requires non-null query and reference";
		assert(maxEdits>=0) : "QuantumRanker maxEdits must be nonnegative: "+maxEdits;
		assert(refStart>=0 && refStop<ref.length && refStart<=refStop) :
			"Window must be inside reference: "+refStart+"-"+refStop+" of "+ref.length;
		final int qLen=query.length, rLen=refStop-refStart+1;
		if(qLen<1 || rLen<1){result.failureCode=FAIL_EMPTY;return result;}
		if(qLen>rLen){result.failureCode=FAIL_QUERY_LONGER;return result;}
		if(rLen>MAX_WINDOW){result.failureCode=FAIL_WINDOW_TOO_LONG;return result;}
		final int bandWidth=decideBandwidth(query, ref, refStart, rLen, expectedStart,
				useSIMD);
		final int effectiveMaxEdits=(maxEdits==Integer.MAX_VALUE ? Integer.MAX_VALUE :
				Math.min(maxEdits, bandWidth));
		result.bandwidth=bandWidth;
		result.editBudget=effectiveMaxEdits;
		ensureCapacity(rLen+4);
		if(useSIMD){
			for(int j=0; j<rLen; j++){refCodes[j]=ref[refStart+j];}
		}
		final int topWidth=Math.min(qLen, bandWidth*2);
		final int insPad=-16-rLen/128-(int)Math.sqrt(rLen)-(5*Math.max(0, rLen-qLen))/4;
		final int sideWidth0=1, sideWidthMax=Math.min(qLen, rLen);
		final int scoreWidth0=bandWidth+1;

		for(int j=0; j<=rLen; j++){
			prevScore[j]=0;
			prevMeta[j]=pack(j, 0);
			prevEdits[j]=0;
			active[j]=j;
		}
		if(traceback){
			Arrays.fill(prevNode, 0, rLen+1, -1);
			Arrays.fill(currNode, 0, rLen+1, -1);
			traceSize=0;
		}
		int activeSize=rLen+1;
		int nextSize=0;
		for(int j=0; (j<=sideWidth0 || j<=topWidth*2) && j<qLen; j++){next[nextSize++]=j;}

		int bridgeTime=BRIDGE_PERIOD;
		int maxPos=0, maxMeta=0, maxNode=-1, maxRowEdits=EDIT_BAD;
		int maxScore=BAD, prevRowScore=BAD;
		long cells=0;
		boolean ambiguous=false, prunedForEdits=false;

		for(int i=1; i<=qLen; i++){
			currScore[0]=(i<=effectiveMaxEdits ? i*INS : BAD);
			currMeta[0]=0;
			currEdits[0]=(i<=effectiveMaxEdits ? i : EDIT_BAD);
			if(traceback){currNode[0]=(i<=effectiveMaxEdits ?
					addTraceNode(prevNode[0], 1, (byte)'I') : -1);}
			prunedForEdits|=(i>effectiveMaxEdits);
			while(activeSize>0 && active[activeSize-1]>rLen){activeSize--;}
			if(activeSize<2){
				result.cells=cells;result.editBudgetExceeded=prunedForEdits;
				result.failureCode=(prunedForEdits ? FAIL_EDIT_BUDGET : FAIL_NO_PATH);
				return result;
			}

			final byte q=query[i-1];
			ambiguous|=!canonical(q);
			final boolean nextMatch=(q==ref[refStart+Math.min(rLen-1, maxPos)] && q!='N');
			final boolean bridge=(bridgeTime<1 && !nextMatch);
			final int extra=nextMatch ? 2 : bridge ? 35 : 3;
			bridgeTime=(bridge ? BRIDGE_PERIOD : bridgeTime-1);
			final int oldLast=active[activeSize-1];
			final int extendLimit=Math.min(oldLast+extra, rLen);
			for(int e=oldLast+1; e<extendLimit; e++){active[activeSize++]=e;}
			cells+=activeSize-1;
			if(active[activeSize-1]<rLen){
				prevScore[rLen]=BAD; prevMeta[rLen]=0;
				if(traceback){prevNode[rLen]=-1;}
			}

			prevRowScore=maxScore;
			maxScore=BAD;
			maxPos=0;
			maxMeta=0;
			maxNode=-1;
			maxRowEdits=EDIT_BAD;

			final int sideWidth=mid(sideWidth0, topWidth*2-i, sideWidthMax);
			assert(nextSize>=sideWidth) : "Sparse frontier lost required side band: "+nextSize+" < "+sideWidth;
			nextSize=sideWidth;
			assert(next[nextSize-1]+1==sideWidth || rLen<sideWidth) :
				"Quantum side band must remain consecutive through "+sideWidth;
			final int scoreWidth=scoreWidth0+Math.max(0, topWidth-i);
			final int forcedInsOffset=insPad-i;
			if(useSIMD){
				simd.QuantumSIMD.diagonalUp((int)q, refCodes, active, 1, activeSize,
						prevScore, prevMeta, diagonalUpScore, diagonalUpMeta);
			}

			for(int idx=1; idx<activeSize; idx++){
				final int j=active[idx];
				final byte r=ref[refStart+j-1];
				final boolean hasN=(q=='N' || r=='N');
				ambiguous|=(!canonical(q) || !canonical(r));
				final boolean match=(q==r && q!='N');
				final int add=match ? MATCH : (hasN ? N_SCORE : SUB);

				final int diagonal, up;
				if(useSIMD){
					diagonal=diagonalUpScore[idx];
					up=Integer.MIN_VALUE;
				}else{
					diagonal=prevScore[j-1]+add;
					up=prevScore[j]+INS;
				}
				final int directLeft=currScore[j-1]+DEL;
				final int jumpLeft=maxScore+DEL*(j-maxPos);
				final boolean useDirectLeft=(directLeft>=jumpLeft);
				final int left=(useDirectLeft ? directLeft : jumpLeft);
				int value=diagonal;
				int meta=(useSIMD ? diagonalUpMeta[idx] : prevMeta[j-1]);
				final boolean useUp;
				if(useSIMD){
					final int diagonalScore=prevScore[j-1]+add;
					final int upScore=prevScore[j]+INS;
					useUp=(upScore>diagonalScore);
				}else{
					useUp=(up>value);
					if(useUp){value=up; meta=prevMeta[j];}
				}
				int editCount=(useUp ? prevEdits[j]+1 :
						prevEdits[j-1]+(match ? 0 : 1));
				int parentNode=-1, run=1;
				byte op=0;
				if(traceback){
					parentNode=(useUp ? prevNode[j] : prevNode[j-1]);
					op=(byte)(useUp ? 'I' : match ? 'm' : hasN ? 'N' : 'S');
				}
				if(left>value){
					value=left;
					meta=(useDirectLeft ? incrementDeletion(currMeta[j-1], 1) :
							incrementDeletion(maxMeta, j-maxPos));
					editCount=(useDirectLeft ? currEdits[j-1]+1 :
							maxRowEdits+(j-maxPos));
					if(traceback){
						parentNode=(useDirectLeft ? currNode[j-1] : maxNode);
						run=(useDirectLeft ? 1 : j-maxPos);
						op='D';
					}
				}
				value-=Math.max(0, (j+forcedInsOffset)/2);
				final boolean withinBudget=(editCount<=effectiveMaxEdits && value>BAD/2);
				prunedForEdits|=!withinBudget && editCount>effectiveMaxEdits && editCount<EDIT_BAD;
				final int node=(traceback && withinBudget ?
						addTraceNode(parentNode, run, op) : -1);

				final int scoreDifference=prevRowScore-value;
				final int last=next[nextSize-1];
				final boolean addPosition=withinBudget &&
						(j<sideWidth || scoreDifference<scoreWidth);
				final boolean live=(withinBudget && match && last<j+1);
				final boolean retain=(addPosition || live);
				currScore[j]=(retain ? value : BAD);
				currMeta[j]=(retain ? meta : 0);
				currEdits[j]=(retain ? editCount : EDIT_BAD);
				if(traceback){currNode[j]=(retain ? node : -1);}
				prevScore[j-1]=BAD;
				prevMeta[j-1]=0;
				prevEdits[j-1]=EDIT_BAD;
				if(traceback){prevNode[j-1]=-1;}

				if(addPosition){
					final int jp1=j+1, jp2=j+2, jp3=j+3;
					if(last==jp2 && jp2<rLen){next[nextSize++]=jp3;}
					else{
						int tail=last;
						if(tail<j){next[nextSize++]=j; tail=j;}
						if(tail<jp1){next[nextSize++]=jp1; tail=jp1;}
						if(tail<jp2){next[nextSize++]=jp2; tail=jp2;}
						if(tail<jp3){next[nextSize++]=jp3;}
					}
				}else if(live){next[nextSize++]=j+1;}

				if(withinBudget && value>maxScore){
					maxScore=value; maxPos=j; maxMeta=meta; maxRowEdits=editCount;
					if(traceback){maxNode=node;}
				}
			}
			if(maxScore<=BAD/2){
				result.cells=cells;result.editBudgetExceeded=prunedForEdits;
				result.failureCode=(prunedForEdits ? FAIL_EDIT_BUDGET : FAIL_NO_PATH);
				return result;
			}

			final int oldActiveSize=activeSize;
			int[] temp=prevScore; prevScore=currScore; currScore=temp;
			temp=prevMeta; prevMeta=currMeta; currMeta=temp;
			temp=prevEdits; prevEdits=currEdits; currEdits=temp;
			if(traceback){temp=prevNode; prevNode=currNode; currNode=temp;}
			temp=active; active=next; next=temp;
			activeSize=nextSize;
			nextSize=oldActiveSize;
		}

		final int origin=origin(maxMeta), deletions=deletions(maxMeta);
		final int refAlnLength=maxPos-origin;
		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final int matches=Math.max(0, (maxScore+qLen+deletions)/2);
		final int substitutions=Math.max(0, qLen-matches-insertions);
		final int denominator=matches+substitutions+insertions+deletions;

		result.supported=true;
		result.vectorized=useSIMD;
		result.uncertain=ambiguous;
		result.ambiguous=ambiguous;
		result.score=maxScore;
		result.rStart=refStart+origin;
		result.rStop=refStart+maxPos-1;
		result.matches=matches;
		result.substitutions=substitutions;
		result.insertions=insertions;
		result.deletions=deletions;
		result.identity=(denominator<1 ? 0 : matches/(float)denominator);
		result.cells=cells;
		result.editBudgetExceeded=false;
		if(traceback){
			result.match=traceback(maxNode);
			installTraceCounts(result, qLen, maxPos-origin);
		}
		assert(result.rStart>=refStart && result.rStop<=refStop) :
			"Result must remain inside supplied window: "+result;
		assert(result.substitutions+result.insertions+result.deletions<=effectiveMaxEdits) :
			"Supported Quantum result exceeds edit budget: "+result+" > "+effectiveMaxEdits;
		return result;
	}

	private void ensureCapacity(final int columns){
		if(prevScore!=null && prevScore.length>=columns){return;}
		int capacity=1;
		while(capacity<columns){capacity<<=1;}
		prevScore=new int[capacity];
		currScore=new int[capacity];
		prevMeta=new int[capacity];
		currMeta=new int[capacity];
		prevEdits=new int[capacity];
		currEdits=new int[capacity];
		active=new int[capacity];
		next=new int[capacity];
		refCodes=new int[capacity];
		diagonalUpScore=new int[capacity];
		diagonalUpMeta=new int[capacity];
		prevNode=new int[capacity];
		currNode=new int[capacity];
	}

	private int addTraceNode(final int parent, final int run, final byte op){
		assert(parent<traceSize) : "Trace parent must precede child: "+parent+" >= "+traceSize;
		assert(run>0 && run<=TRACE_RUN_MASK) : "Trace run exceeds 24 bits: "+run;
		if(traceSize>=traceNodes.length){
			final int next=Math.max(1024, traceNodes.length*2);
			traceNodes=Arrays.copyOf(traceNodes, next);
		}
		traceNodes[traceSize]=(((long)parent)<<32)|(((long)run)<<8)|(op&0xFFL);
		return traceSize++;
	}

	private byte[] traceback(final int finalNode){
		assert(finalNode>=0 && finalNode<traceSize) :
				"A supported traced alignment requires a terminal node: "+finalNode+" / "+traceSize;
		int length=0, node=finalNode, remaining=traceSize+1;
		while(node>=0){
			final long packed=traceNodes[node];
			length+=(int)((packed>>>8)&TRACE_RUN_MASK);
			node=(int)(packed>>32);
			assert(remaining-->0) : "Cycle while sizing QuantumRanker traceback";
		}
		final byte[] match=new byte[length];
		int position=length;
		node=finalNode;
		remaining=traceSize+1;
		while(node>=0){
			final long packed=traceNodes[node];
			final int run=(int)((packed>>>8)&TRACE_RUN_MASK);
			final byte op=(byte)packed;
			for(int i=0; i<run; i++){match[--position]=op;}
			node=(int)(packed>>32);
			assert(remaining-->0) : "Cycle while materializing QuantumRanker traceback";
		}
		assert(position==0) : "Traceback must fill the allocated match array: "+position;
		return match;
	}

	private static void installTraceCounts(final Result result, final int qLen,
			final int refLen){
		int m=0, s=0, ins=0, del=0, n=0;
		for(byte op : result.match){
			if(op=='m'){m++;}else if(op=='S'){s++;}else if(op=='I'){ins++;}
			else if(op=='D'){del++;}else if(op=='N'){n++;}
			else{assert(false) : "Unknown QuantumRanker trace operation: "+(char)op;}
		}
		assert(m+s+ins+n==qLen) : "Trace query consumption differs: "+(m+s+ins+n)+" != "+qLen;
		assert(m+s+del+n==refLen) : "Trace reference consumption differs: "+(m+s+del+n)+" != "+refLen;
		result.matches=m;
		result.substitutions=s;
		result.insertions=ins;
		result.deletions=del;
		result.score=m-s-ins-del;
		final int denominator=m+s+ins+del+n;
		result.identity=(denominator<1 ? 0 : m/(float)denominator);
	}

	private int decideBandwidth(final byte[] query, final byte[] ref,
			final int refStart, final int rLen, final int expectedStart,
			final boolean useSIMD){
		final int qLen=query.length, maxLen=Math.max(qLen, rLen);
		final int log2=(int)(Math.log(maxLen+256)*INV_LN_2);
		int maxBandwidth=Math.min(qLen/4+2, Math.min(maxLen/32, log2+2));
		maxBandwidth=Math.max(2, maxBandwidth)+3;
		if(useSIMD){
			prealignPos[0]=expectedStart;
			prealignPos[1]=0;
			int substitutions=IDAlignerStatics.countSubs(query, ref, prealignPos,
					maxBandwidth);
			if(substitutions<maxBandwidth){
				return Math.min(substitutions+1, maxBandwidth);
			}
		}
		final int expectedOffset=mid(0, expectedStart-refStart, rLen-1);
		int substitutions=0, q=0;
		for(; q<qLen && substitutions<maxBandwidth; q++){
			final int r=q+expectedOffset;
			substitutions+=(r<0 || r>=rLen || query[q]!=ref[refStart+r] ? 1 : 0);
		}
		return Math.min(substitutions+1, maxBandwidth);
	}

	private static int pack(final int origin, final int deletions){
		assert(origin>=0 && origin<=META_MASK) : "Origin exceeds 16-bit window coordinate: "+origin;
		assert(deletions>=0 && deletions<=META_MASK) : "Deletion count exceeds 16 bits: "+deletions;
		return (origin<<META_BITS)|deletions;
	}

	private static int incrementDeletion(final int meta, final int amount){
		final int count=deletions(meta);
		return (count>=META_MASK-amount ? (meta|META_MASK) : meta+amount);
	}

	private static int origin(final int meta){return meta>>>META_BITS;}
	private static int deletions(final int meta){return meta&META_MASK;}
	private static int mid(final int low, final int value, final int high){return Math.max(low, Math.min(value, high));}
	private static boolean canonical(final byte base){
		return base=='A' || base=='C' || base=='G' || base=='T';
	}

	/*--------------------------------------------------------------*/
	/*----------------            Result            ----------------*/
	/*--------------------------------------------------------------*/

	public static class Result {
		private void clear(){
			supported=uncertain=ambiguous=vectorized=editBudgetExceeded=prealignmentOnly=false;
			score=rStart=rStop=matches=substitutions=insertions=deletions=failureCode=0;
			bandwidth=editBudget=0;
			identity=0;
			cells=0;
			match=null;
		}
		@Override
		public String toString(){
			return "Result(score="+score+", start="+rStart+", stop="+rStop+
					", identity="+identity+", M="+matches+", S="+substitutions+
					", I="+insertions+", D="+deletions+", cells="+cells+")";
		}
		public boolean supported;
		public boolean uncertain;
		public boolean ambiguous;
		public boolean vectorized;
		public boolean editBudgetExceeded;
		public boolean prealignmentOnly;
		public int score;
		public int rStart;
		public int rStop;
		public int matches;
		public int substitutions;
		public int insertions;
		public int deletions;
		public int bandwidth;
		public int editBudget;
		public int failureCode;
		public float identity;
		public long cells;
		public byte[] match;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int[] prevScore;
	private int[] currScore;
	private int[] prevMeta;
	private int[] currMeta;
	private int[] prevEdits;
	private int[] currEdits;
	private int[] active;
	private int[] next;
	private int[] refCodes;
	private int[] diagonalUpScore;
	private int[] diagonalUpMeta;
	private int[] prevNode;
	private int[] currNode;
	private final int[] prealignPos=new int[2];
	private long[] traceNodes=new long[1024];
	private int traceSize;
	private final Result result=new Result();

	private static final int META_BITS=16;
	private static final int META_MASK=(1<<META_BITS)-1;
	private static final int MAX_WINDOW=META_MASK;
	private static final int MATCH=1;
	private static final int SUB=-1;
	private static final int INS=-1;
	private static final int DEL=-1;
	private static final int N_SCORE=0;
	private static final int BAD=Integer.MIN_VALUE/4;
	private static final int EDIT_BAD=1<<28;
	private static final int BRIDGE_PERIOD=16;
	private static final int TRACE_RUN_MASK=0xFFFFFF;
	private static final double INV_LN_2=1.4426950408889634;

	public static final int FAIL_NONE=0;
	public static final int FAIL_EMPTY=1;
	public static final int FAIL_QUERY_LONGER=2;
	public static final int FAIL_WINDOW_TOO_LONG=3;
	public static final int FAIL_EDIT_BUDGET=4;
	public static final int FAIL_NO_PATH=5;
}
