package prot;

import java.util.Arrays;

/** Immutable experimental position scores; the integer DP uses BLOSUM units times 64. */
public final class HbmPositionModel {
	public static final int SCALE=64, GAP=GlocalAminoLinear.GAP*SCALE;
	private final int[][] scores;
	private final int length;
	private final byte[] consensus;
	private HbmPositionModel(int[][] scores_){
		this(scores_,null);
	}
	private HbmPositionModel(int[][] scores_,byte[] consensus_){
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
	static HbmPositionModel deriveOne(AAGraph graph,String kind,double beta,boolean clip,double[] background,double clipMin){
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
		return new HbmPositionModel(table,pivot);
	}
	private static void validateFloor(double clipMin){
		if(!Double.isFinite(clipMin) || clipMin< -100 || clipMin>0){
			throw new IllegalArgumentException("clipmin must be finite and in [-100,0] to preserve the bounded DP score range: "+clipMin);
		}
	}
	/** Fixture factory copies every row, preserving the model's immutability. */
	static HbmPositionModel fromScores(int[][] table){return new HbmPositionModel(table);}
	int scoreAt(byte aa,int position){return scores[aa][position];}
	public int length(){return length;}
	private static final class Scratch{
		int[] matrix=new int[0];
		void ensure(int n){if(matrix.length<n){matrix=new int[n];}}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);
	/** Same glocal recurrence, end scan and tie order as GlocalAminoLinear; only emissions differ. */
	public Result align(byte[] query,boolean record){
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
			int diag=h[prev],left=h[cur]=-i*GAP;
			for(int j=1; j<=n; j++){
				final int up=h[prev+j];
				int value=diag+row[j-1];
				final int u=up-GAP; if(u>value){value=u;}
				final int l=left-GAP; if(l>value){value=l;}
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
			}else if(h[i*stride+j]==h[(i-1)*stride+j]-GAP){
				if(record){ops[--at]='I';} i--;
			}else{
				if(j<1 || h[i*stride+j]!=h[i*stride+j-1]-GAP){throw new AssertionError("Profile traceback has no valid predecessor");}
				start=j-1; if(record){ops[--at]='D';} j--;
			}
		}
		return new Result(best,start,bestJ-1,record ? Arrays.copyOfRange(ops,at,ops.length) : null);
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
			}else if(op=='I'){sum-=GAP;q++;}
			else if(op=='D'){sum-=GAP;r++;}
			else{throw new IllegalArgumentException("Unknown traceback operation: "+op);}
		}
		if(q!=query.length || r>length){throw new IllegalArgumentException("Trace does not consume the full query within the model");}
		if(sum<Integer.MIN_VALUE || sum>Integer.MAX_VALUE){throw new ArithmeticException("Profile path score overflow");}
		return (int)sum;
	}
	/** Scores are fixed-point profile units, with no BLOSUM-specific E-value interpretation. */
	public static final class Result{
		Result(int score_,int start_,int end_,byte[] path_){score=score_;start=start_;end=end_;path=path_;}
		public final int score,start,end;
		public final byte[] path;
	}
}
