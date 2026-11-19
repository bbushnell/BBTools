package prot;

import java.util.Arrays;

/** Immutable position scores; integer DP uses selected score units times64.
 * Derived legacy profiles retain BLOSUM62 units and the historical gap cost.
 */
public final class HbmPositionModel {
	public static final int SCALE=64, GAP=GlocalAminoLinear.GAP*SCALE;
	private final int[][] scores;
	private final int length;
	private final byte[] consensus;
	private final int gap;
	private HbmPositionModel(int[][] scores_){
		this(scores_,null);
	}
	private HbmPositionModel(int[][] scores_,byte[] consensus_){
		this(scores_,consensus_,GAP);
	}
	private HbmPositionModel(int[][] scores_,byte[] consensus_,int gap_){
		if(gap_<1 || gap_>6400){throw new IllegalArgumentException("Gap exceeds bounded fixed-point DP range");}
		gap=gap_;
		if(scores_.length!=22 || scores_[0].length<1){throw new IllegalArgumentException("Position matrix must have 22 rows and nonzero length");}
		length=scores_[0].length; scores=new int[22][];
		if(consensus_!=null && consensus_.length!=length){throw new IllegalArgumentException("Profile consensus length differs from score columns");}
		consensus=consensus_==null ? null : consensus_.clone();
		for(int a=0; a<22; a++){
			if(scores_[a].length!=length){throw new IllegalArgumentException("Ragged position matrix");}
			if(a!=20){for(int value : scores_[a]){if(value< -6400 || value>6400){throw new IllegalArgumentException("Position score exceeds DP bound");}}}
			scores[a]=scores_[a].clone();
		}
	}
	/** Derives detached score tables; neither graphs nor their arrays escape this call. */
	static HbmPositionModel[] derive(AAGraph[] graphs,String kind,double beta,boolean clip){
		return derive(graphs,kind,beta,clip,-4);
	}
	/** The floor is in half-bit score units before fixed-point rounding; the ceiling remains 11. */
	static HbmPositionModel[] derive(AAGraph[] graphs,String kind,double beta,boolean clip,double clipMin){
		validateFloor(clipMin);
		if(!kind.equals("expected") && !kind.equals("logodds") && !kind.equals("consensus")){throw new IllegalArgumentException("Unknown position score kind: "+kind);}
		if(!(beta>0 && beta<1)){throw new IllegalArgumentException("Smoothing mixture must be in (0,1)");}
		final double[] background=backgroundProbabilities(graphs);
		final HbmPositionModel[] models=new HbmPositionModel[graphs.length];
		for(int i=0; i<graphs.length; i++){models[i]=deriveOne(graphs[i],kind,beta,clip,background,clipMin);}
		return models;
	}
	/** Fresh global-background vector, using exactly the original assay's REF-count rule. */
	static double[] backgroundProbabilities(AAGraph[] graphs){
		final double[] background=new double[20]; double total=20;
		Arrays.fill(background,1);//Finite background even if a training alphabet omits a residue.
		for(AAGraph g : graphs){for(AAGraphNode node : g.ref){for(int a=0; a<20; a++){background[a]+=node.count[a]; total+=node.count[a];}}}
		for(int a=0; a<20; a++){background[a]/=total;}
		return background;
	}
	static HbmPositionModel deriveOne(AAGraph graph,String kind,double beta,boolean clip,double[] background){
		return deriveOne(graph,kind,beta,clip,background,-4);
	}
	public static HbmPositionModel deriveOne(AAGraph graph,String kind,double beta,boolean clip,double[] background,double clipMin){
		final int[][] counts=new int[graph.ref.length][];
		for(int j=0; j<counts.length; j++){counts[j]=graph.ref[j].count;}
		return deriveColumns(graph.pivot,counts,kind,beta,clip,background,clipMin);
	}
	/** Carries selected counts once; it never adds another scaffold observation. */
	static HbmPositionModel deriveTraversal(AAGraph.Traversal traversal,double[] background){
		final int[][] counts=new int[traversal.length()][22];
		for(int j=0; j<counts.length; j++){for(int a=0; a<22; a++){counts[j][a]=traversal.count(j,a);}}
		return deriveColumns(traversal.consensus(),counts,"logodds",0.01,true,background,-4);
	}
	private static HbmPositionModel deriveColumns(byte[] pivot,int[][] counts,String kind,double beta,boolean clip,double[] background,double clipMin){
		validateFloor(clipMin);
		if(!(beta>0 && beta<1) || (!kind.equals("expected") && !kind.equals("logodds") && !kind.equals("consensus"))){throw new IllegalArgumentException("Invalid profile kind or smoothing mixture");}
		if(background.length!=20){throw new IllegalArgumentException("Profile background needs 20 real-residue probabilities");}
		double sum=0;
		for(double p : background){if(!(p>0 && p<=1) || !Double.isFinite(p)){throw new IllegalArgumentException("Invalid background probability: "+p);}sum+=p;}
		if(Math.abs(sum-1)>1e-9){throw new IllegalArgumentException("Profile background probabilities must sum to one: "+sum);}
		assert(pivot.length==counts.length) : "Derived emission columns must share the pivot coordinate frame";
		final int[][] table=new int[22][counts.length];
		Arrays.fill(table[20],Integer.MIN_VALUE/4);//Encoded stop is not a query residue.
		for(int j=0; j<counts.length; j++){
			long count=0; for(int a=0; a<20; a++){count+=counts[j][a];}
			for(int a=0; a<22; a++){
				if(a==20){continue;}
				double value=0;
				if(kind.equals("consensus")){value=Blosum62.score((byte)a,pivot[j]);}
				else if(kind.equals("expected")){
					for(int b=0; b<20; b++){
						final double p=count==0 ? background[b] : counts[j][b]/(double)count;
						value+=p*Blosum62.score((byte)a,(byte)b);
					}
				}else if(a!=Blosum62.X_CODE){
					final double p=count==0 ? background[a] : counts[j][a]/(double)count;
					value=2*Math.log(((1-beta)*p+beta*background[a])/background[a])/Math.log(2);
					if(clip){value=Math.max(clipMin,Math.min(11,value));}
				}
				if(!Double.isFinite(value) || Math.abs(value)>100){throw new IllegalArgumentException("Invalid position score: "+value);}
				table[a][j]=(int)Math.round(SCALE*value);
			}
		}
		return new HbmPositionModel(table,pivot);
	}
	/** Fails before graph accumulation if scores belong to a different reference frame. */
	void requireConsensus(byte[] expected){
		if(consensus==null || !Arrays.equals(consensus,expected)){throw new IllegalArgumentException("Position profile is not bound to the supplied consensus");}
	}
	/** Adds neutral background columns for the consensus-growth pass, matching an X-padded pivot. */
	HbmPositionModel padded(int pad){
		if(pad<0 || (long)length+2L*pad>Integer.MAX_VALUE-8 || consensus==null){throw new IllegalArgumentException("Invalid padding or unbound profile");}
		if(pad==0){return this;}
		final int n=length+2*pad;
		final byte[] pivot=new byte[n]; Arrays.fill(pivot,Blosum62.X_CODE);
		System.arraycopy(consensus,0,pivot,pad,length);
		final int[][] table=new int[22][n]; Arrays.fill(table[20],Integer.MIN_VALUE/4);
		for(int a=0; a<22; a++){System.arraycopy(scores[a],0,table[a],pad,length);}
		return new HbmPositionModel(table,pivot,gap);
	}
	private static void validateFloor(double clipMin){
		if(!Double.isFinite(clipMin) || clipMin< -100 || clipMin>0){
			throw new IllegalArgumentException("clipmin must be finite and in [-100,0] to preserve the bounded DP score range: "+clipMin);
		}
	}
	/** Fixture factory copies every row, preserving the model's immutability. */
	public static HbmPositionModel fromScores(int[][] table){return new HbmPositionModel(table);}
	/** Direct BLOSUM62 consensus scorer, matching deriveOne(..., "consensus", ...). */
	public static HbmPositionModel fromConsensus(byte[] consensus){
		return fromConsensus(consensus,ProteinSubstitutionMatrix.BLOSUM62);
	}
	/** Direct consensus scores with an explicit immutable substitution/gap binding. */
	public static HbmPositionModel fromConsensus(byte[] consensus,ProteinSubstitutionMatrix matrix){
		if(matrix==null){throw new IllegalArgumentException("Consensus model requires a substitution matrix");}
		final int[][] table=consensusTable(consensus,matrix);
		return new HbmPositionModel(table,consensus,matrix.linearGap*SCALE);
	}
	public int scoreAt(byte aa,int position){return scores[aa][position];}
	public int length(){return length;}
	private static int[][] consensusTable(byte[] consensus,ProteinSubstitutionMatrix matrix){
		if(consensus==null || consensus.length<1){throw new IllegalArgumentException("Consensus profile needs at least one residue");}
		final int[][] table=new int[22][consensus.length];
		fillConsensus(table,consensus,matrix);
		return table;
	}
	private static void fillConsensus(int[][] table,byte[] consensus,ProteinSubstitutionMatrix matrix){
		assert(table.length==22) : "Consensus scorers use the full 22-symbol encoded amino-acid alphabet";
		Arrays.fill(table[20], 0, consensus.length, Integer.MIN_VALUE/4);
		for(int j=0; j<consensus.length; j++){
			final byte pivot=consensus[j];
			if(pivot<0 || pivot>21 || pivot==20){throw new IllegalArgumentException("Invalid consensus residue: "+pivot);}
			for(int a=0; a<22; a++){
				if(a!=20){table[a][j]=SCALE*matrix.score((byte)a,pivot);}
			}
		}
	}
	private static final class Scratch{
		int[] matrix=new int[0];
		void ensure(int n){if(matrix.length<n){matrix=new int[n];}}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);
	/** Same glocal recurrence, end scan and tie order as GlocalAminoLinear; only emissions differ. */
	public Result align(byte[] query,boolean record){
		return align(query,record,scores,length,gap);
	}
	private static Result align(byte[] query,boolean record,int[][] scores,int length,int gap){
		final int m=query.length,n=length,stride=n+1;
		if(m<1 || (long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Empty/oversized profile alignment: "+m+" x "+n);
		}
		final Scratch scratch=SCRATCH.get(); scratch.ensure((m+1)*stride);
		final int[] h=scratch.matrix;
		Arrays.fill(h,0,stride,0);
		for(int i=1; i<=m; i++){
			final byte a=query[i-1];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded query residue: "+a);}
			final int[] row=scores[a];
			final int cur=i*stride,prev=cur-stride;
			int diag=h[prev],left=h[cur]=-i*gap;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=diag+row[j-1];
				final int u=up-gap; if(u>value){value=u;}
				final int l=left-gap; if(l>value){value=l;}
				h[cur+j]=value; left=value; diag=up;
			}
		}
		int bestJ=1,best=h[m*stride+1];
		for(int j=2; j<=n; j++){if(h[m*stride+j]>best){best=h[m*stride+j];bestJ=j;}}
		int i=m,j=bestJ,start=bestJ-1;
		final byte[] ops=record ? new byte[m+n] : null; int at=record ? ops.length : 0;
		while(i>0){
			if(j>0 && h[i*stride+j]==h[(i-1)*stride+j-1]+scores[query[i-1]][j-1]){
				start=j-1; if(record){ops[--at]='m';} i--; j--;
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-gap){
				if(record){ops[--at]='I';} i--;
			}else{
				if(j<1 || h[i*stride+j]!=h[i*stride+j-1]-gap){throw new AssertionError("Profile traceback has no valid predecessor");}
				start=j-1; if(record){ops[--at]='D';} j--;
			}
		}
		return new Result(best,start,bestJ-1,record ? Arrays.copyOfRange(ops,at,ops.length) : null);
	}
	/**
	 * Aligns bounded C-terminal windows, with free left clipping on both sequences.
	 * Both right ends are consumed; terminal gaps retain their ordinary penalties.
	 * Returns null when no positive-scoring suffix reaches both ends. This is an
	 * alignment primitive, not a stop-validation gate: callers must check aligned
	 * span and terminal support. Coordinates refer to the original query/profile.
	 */
	public TailResult alignTail(byte[] query, int maxLength, boolean record){
		return alignTail(query, maxLength, record, true);
	}
	/** Right-pinned bounded tail alignment, with optional query left clipping. */
	public TailResult alignTail(byte[] query, int maxLength, boolean record, boolean trimLeft){
		return alignTail(query,maxLength,record,trimLeft,scores,length,gap);
	}
	/** Right-pinned overlap: only top/left boundary starts are free; interior scores
	 * never restart at zero. Returns the signed endpoint score and traces to a boundary,
	 * even through nonpositive cells. A gap-only optimum has zero span on one side;
	 * consumers must apply their span gates rather than invent a paired overlap.
	 */
	public TailResult alignTailOverhang(byte[] query,int maxLength,boolean record){
		return alignTail(query,maxLength,record,true,scores,length,gap,true);
	}

	/** Left-pinned bounded head alignment, with optional query right clipping. */
	public TailResult alignHead(byte[] query, int maxLength, boolean record, boolean trimRight){
		return alignHead(query,maxLength,record,trimRight,scores,length,gap);
	}

	private static TailResult alignTail(byte[] query, int maxLength, boolean record, boolean trimLeft,
			int[][] scores,int length,int gap){
		return alignTail(query,maxLength,record,trimLeft,scores,length,gap,false);
	}
	private static TailResult alignTail(byte[] query,int maxLength,boolean record,boolean trimLeft,
			int[][] scores,int length,int gap,boolean overhangOnly){
		assert(!overhangOnly || trimLeft) : "Boundary-only overlap requires free starts on both top and left edges";
		if(query==null || query.length<1 || maxLength<1){
			throw new IllegalArgumentException("Tail alignment requires a nonempty query and positive window limit");
		}
		final int m=Math.min(query.length, maxLength), n=Math.min(length, maxLength);
		final int qOffset=query.length-m, rOffset=length-n, stride=n+1;
		if((long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Oversized tail alignment: "+m+" x "+n);
		}
		final Scratch scratch=SCRATCH.get(); scratch.ensure((m+1)*stride);
		final int[] h=scratch.matrix;
		Arrays.fill(h, 0, stride, 0);
		for(int i=1; i<=m; i++){
			final byte a=query[qOffset+i-1];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded tail residue: "+a);}
			final int[] row=scores[a];
			final int cur=i*stride, prev=cur-stride;
			int diag=h[prev], left=h[cur]=trimLeft ? 0 : -i*gap;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=trimLeft && !overhangOnly ? Math.max(0, diag+row[rOffset+j-1]) : diag+row[rOffset+j-1];
				value=Math.max(value, up-gap);
				value=Math.max(value, left-gap);
				h[cur+j]=value; left=value; diag=up;
			}
		}
		final int score=h[m*stride+n];
		if(!overhangOnly && score<=0){return null;}
		int i=m, j=n;
		final byte[] ops=record ? new byte[m+n] : null;
		int at=record ? ops.length : 0;
		while(i>0 && (trimLeft ? j>0 && (overhangOnly || h[i*stride+j]>0) : true)){
			if(j>0 && h[i*stride+j]==h[(i-1)*stride+j-1]+scores[query[qOffset+i-1]][rOffset+j-1]){
				if(record){ops[--at]='m';} i--; j--;
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-gap){
				if(record){ops[--at]='I';} i--;
			}else{
				if(j<1 || h[i*stride+j]!=h[i*stride+j-1]-gap){throw new AssertionError("Tail traceback has no valid predecessor");}
				if(record){ops[--at]='D';} j--;
			}
		}
		assert(overhangOnly || i<m && j<n) : "Positive tail evidence must consume residues from both sequences";
		assert(!overhangOnly || i==0 || j==0) : "Overlap traceback must reach a free boundary, never an interior restart";
		return new TailResult(score, qOffset+i, query.length-1, rOffset+j, length-1,
			record ? Arrays.copyOfRange(ops, at, ops.length) : null);
	}
	private static TailResult alignHead(byte[] query, int maxLength, boolean record, boolean trimRight,
			int[][] scores,int length,int gap){
		if(query==null || query.length<1 || maxLength<1){
			throw new IllegalArgumentException("Head alignment requires a nonempty query and positive window limit");
		}
		final int m=Math.min(query.length, maxLength), n=Math.min(length, maxLength), stride=n+1;
		if((long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Oversized head alignment: "+m+" x "+n);
		}
		final Scratch scratch=SCRATCH.get(); scratch.ensure((m+1)*stride);
		final int[] h=scratch.matrix;
		Arrays.fill(h, 0, stride, 0);
		for(int i=1; i<=m; i++){
			final byte a=query[m-i];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded head residue: "+a);}
			final int[] row=scores[a];
			final int cur=i*stride, prev=cur-stride;
			int diag=h[prev], left=h[cur]=trimRight ? 0 : -i*gap;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=trimRight ? Math.max(0, diag+row[n-j]) : diag+row[n-j];
				value=Math.max(value, up-gap);
				value=Math.max(value, left-gap);
				h[cur+j]=value; left=value; diag=up;
			}
		}
		final int score=h[m*stride+n];
		if(score<=0){return null;}
		int i=m, j=n;
		final byte[] ops=record ? new byte[m+n] : null;
		int at=record ? ops.length : 0;
		while(i>0 && (trimRight ? j>0 && h[i*stride+j]>0 : true)){
			if(j>0 && h[i*stride+j]==h[(i-1)*stride+j-1]+scores[query[m-i]][n-j]){
				if(record){ops[--at]='m';} i--; j--;
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-gap){
				if(record){ops[--at]='I';} i--;
			}else{
				if(j<1 || h[i*stride+j]!=h[i*stride+j-1]-gap){throw new AssertionError("Head traceback has no valid predecessor");}
				if(record){ops[--at]='D';} j--;
			}
		}
		assert(i<m && j<n) : "Positive head evidence must consume residues from both sequences";
		final byte[] path=record ? Arrays.copyOfRange(ops, at, ops.length) : null;
		if(path!=null){reverseInPlace(path);}
		return new TailResult(score, 0, m-1-i, 0, n-1-j, path);
	}
	private static void reverseInPlace(byte[] array){
		for(int left=0,right=array.length-1; left<right; left++,right--){
			final byte b=array[left]; array[left]=array[right]; array[right]=b;
		}
	}
	/**
	 * Exact bottom-right score of {@link #alignTail}, with zero replacing its null result.
	 * Uses the same suffix offsets, zero boundaries/floor, integer emissions and gap penalty.
	 * Only the preceding row and current left cell affect a fill; no geometry is returned.
	 * For raw-score ranking ONLY: callers needing spans/density/trace must keep alignTail.
	 * Reuses one thread-local row, with no result or per-candidate allocation after growth.
	 */
	public int scoreTail(byte[] query,int maxLength){
		if(query==null || query.length<1 || maxLength<1){
			throw new IllegalArgumentException("Tail alignment requires a nonempty query and positive window limit");
		}
		final int m=Math.min(query.length, maxLength), n=Math.min(length, maxLength);
		final int qOffset=query.length-m, rOffset=length-n, stride=n+1;
		// Keep alignTail's domain even though rolling storage alone could support larger matrices.
		if((long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Oversized tail alignment: "+m+" x "+n);
		}
		final ScoreScratch scratch=SCORE_SCRATCH.get();
		if(scratch.row.length<stride){scratch.row=new int[stride];}
		final int[] h=scratch.row;
		Arrays.fill(h, 0, stride, 0);
		for(int i=1; i<=m; i++){
			final byte a=query[qOffset+i-1];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded tail residue: "+a);}
			final int[] row=scores[a];
			int diag=0, left=0;
			for(int j=1; j<=n; j++){
				final int up=h[j];
				int value=Math.max(0, diag+row[rOffset+j-1]);
				value=Math.max(value, up-gap);
				value=Math.max(value, left-gap);
				h[j]=value; left=value; diag=up;
			}
		}
		assert(h[n]>=0 && (long)h[n]<=6400L*Math.min(m,n)) : "Zero floor and bounded emissions cap a suffix score; gap steps cannot increase it";
		return h[n];
	}
	private static final class ScoreScratch{int[] row=new int[0];}
	private static final ThreadLocal<ScoreScratch> SCORE_SCRATCH=ThreadLocal.withInitial(ScoreScratch::new);

	/** Worker-local reusable consensus scorer with immutable substitution/gap units. */
	public static final class ConsensusScorer {
		public ConsensusScorer(){
			this(ProteinSubstitutionMatrix.BLOSUM62);
		}
		/** Workers may use different immutable matrices concurrently. */
		public ConsensusScorer(ProteinSubstitutionMatrix matrix_){
			if(matrix_==null){throw new IllegalArgumentException("Consensus scorer requires a substitution matrix");}
			matrix=matrix_;gap=matrix.linearGap*SCALE;
			scores=new int[22][];
			for(int i=0; i<scores.length; i++){scores[i]=new int[0];}
		}

		/** Replaces the active consensus in-place without reallocating existing large rows. */
		public ConsensusScorer resetConsensus(byte[] consensus){
			length=0;
			if(consensus==null || consensus.length<1){
				throw new IllegalArgumentException("Consensus profile needs at least one residue");
			}
			for(int i=0; i<scores.length; i++){
				if(scores[i].length<consensus.length){scores[i]=new int[consensus.length];}
			}
			fillConsensus(scores,consensus,matrix);
			if(active.length<consensus.length){active=new byte[consensus.length];}
			System.arraycopy(consensus, 0, active, 0, consensus.length);
			length=consensus.length;
			return this;
		}

		public Result align(byte[] query,boolean record){
			requireReady();
			return HbmPositionModel.align(query,record,scores,length,gap);
		}

		public TailResult alignTail(byte[] query,int maxLength,boolean record,boolean trimLeft){
			requireReady();
			return HbmPositionModel.alignTail(query,maxLength,record,trimLeft,scores,length,gap);
		}

		/** Signed boundary-only overlap with the same immutable matrix and reusable rows. */
		public TailResult alignTailOverhang(byte[] query,int maxLength,boolean record){
			requireReady();
			return HbmPositionModel.alignTail(query,maxLength,record,true,scores,length,gap,true);
		}

		public TailResult alignHead(byte[] query,int maxLength,boolean record,boolean trimRight){
			requireReady();
			return HbmPositionModel.alignHead(query,maxLength,record,trimRight,scores,length,gap);
		}

		public int length(){return length;}

		public boolean matchesConsensus(byte[] consensus){
			if(consensus==null || consensus.length!=length){return false;}
			for(int i=0; i<length; i++){if(consensus[i]!=active[i]){return false;}}
			return true;
		}

		private void requireReady(){
			if(length<1){throw new IllegalStateException("Reusable consensus scorer is not initialized");}
		}

		public final ProteinSubstitutionMatrix matrix;
		private final int gap;
		private final int[][] scores;
		private byte[] active=new byte[0];
		private int length;
	}

	/** Diagnostic-only alternative to ask whether an unpaired terminal optimum had a paired-end substitute.
	 * The DP window, free left clipping and scoring are unchanged, but the final traceback step is forced
	 * to pair the last query residue to the last model column. Production tail assignment does not use this.
	 */
	public TailResult alignTailTerminalPair(byte[] query, int maxLength, boolean record){
		if(query==null || query.length<1 || maxLength<1){
			throw new IllegalArgumentException("Tail alignment requires a nonempty query and positive window limit");
		}
		final int m=Math.min(query.length, maxLength), n=Math.min(length, maxLength);
		final int qOffset=query.length-m, rOffset=length-n, stride=n+1;
		if((long)(m+1)*stride>Integer.MAX_VALUE-8 || ((long)m+n)*6400>Integer.MAX_VALUE/2){
			throw new IllegalArgumentException("Oversized tail alignment: "+m+" x "+n);
		}
		final Scratch scratch=SCRATCH.get(); scratch.ensure((m+1)*stride);
		final int[] h=scratch.matrix;
		Arrays.fill(h, 0, stride, 0);
		for(int i=1; i<=m; i++){
			final byte a=query[qOffset+i-1];
			if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded tail residue: "+a);}
			final int[] row=scores[a];
			final int cur=i*stride, prev=cur-stride;
			int diag=h[prev], left=h[cur]=0;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=Math.max(0, diag+row[rOffset+j-1]);
				value=Math.max(value, up-gap);
				value=Math.max(value, left-gap);
				h[cur+j]=value; left=value; diag=up;
			}
		}
		final int score=h[(m-1)*stride+n-1]+scores[query[query.length-1]][length-1];
		if(score<=0){return null;}
		int i=m-1, j=n-1;
		final byte[] ops=record ? new byte[m+n] : null;
		int at=record ? ops.length : 0;
		if(record){ops[--at]='m';}
		while(i>0 && j>0 && h[i*stride+j]>0){
			if(h[i*stride+j]==h[(i-1)*stride+j-1]+scores[query[qOffset+i-1]][rOffset+j-1]){
				if(record){ops[--at]='m';} i--; j--;
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-gap){
				if(record){ops[--at]='I';} i--;
			}else{
				if(h[i*stride+j]!=h[i*stride+j-1]-gap){throw new AssertionError("Tail traceback has no valid predecessor");}
				if(record){ops[--at]='D';} j--;
			}
		}
		return new TailResult(score, qOffset+i, query.length-1, rOffset+j, length-1,
			record ? Arrays.copyOfRange(ops, at, ops.length) : null);
	}
	/** Scores every operation, including terminal insertions, on an existing full-query trace. */
	public int scorePath(byte[] query,int start,byte[] path){
		if(path==null){throw new IllegalArgumentException("Profile rescoring requires a recorded path");}
		for(byte a : query){if(a<0 || a>21 || a==20){throw new IllegalArgumentException("Invalid encoded query residue: "+a);}}
		int q=0,r=start; long sum=0;
		for(byte op : path){
			if(op=='m'){
				if(q>=query.length || r<0 || r>=length){throw new IllegalArgumentException("Profile path overran a sequence");}
				sum+=scores[query[q++]][r++];
			}else if(op=='I'){sum-=gap;q++;}
			else if(op=='D'){sum-=gap;r++;}
			else{throw new IllegalArgumentException("Unknown traceback operation: "+op);}
		}
		if(q!=query.length || r>length){throw new IllegalArgumentException("Trace does not consume the full query within the model");}
		if(sum<Integer.MIN_VALUE || sum>Integer.MAX_VALUE){throw new ArithmeticException("Profile path score overflow");}
		return (int)sum;
	}
	/** Scores are fixed-point profile units, with no BLOSUM-specific E-value interpretation. */
	public static final class TailResult{
		public TailResult(int score_, int queryStart_, int queryEnd_, int start_, int end_, byte[] path_){
			score=score_; queryStart=queryStart_; queryEnd=queryEnd_; start=start_; end=end_; path=path_;
		}
		public final int score, queryStart, queryEnd, start, end;
		public final byte[] path;
	}
	/** Scores are fixed-point profile units, with no BLOSUM-specific E-value interpretation. */
	public static final class Result{
		public Result(int score_,int start_,int end_,byte[] path_){score=score_;start=start_;end=end_;path=path_;}
		public final int score,start,end;
		public final byte[] path;
	}
}
