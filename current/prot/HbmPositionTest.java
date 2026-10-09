package prot;

import java.util.Arrays;
import java.util.Random;

/** Exact BLOSUM parity, exhaustive tiny-path optimality and emission arithmetic oracles. */
public final class HbmPositionTest {
	public static void main(String[] args){
		try{test();}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	private static void test(){
		final Random random=new Random(20261008);
		final double[] bg=new double[20]; Arrays.fill(bg,0.05);
		for(int trial=0; trial<400; trial++){
			final byte[] query=sequence(random,1+random.nextInt(12)),ref=sequence(random,1+random.nextInt(12));
			final HbmPositionModel model=HbmPositionModel.deriveOne(new AAGraph(ref,0),"consensus",0.1,false,bg);
			final AAAlignment old=GlocalAminoLinear.align(query,ref,true);
			final HbmPositionModel.Result newer=model.align(query,true);
			check(newer.score==old.rawScore*64 && newer.start==old.tStart && newer.end==old.tStop && Arrays.equals(newer.path,old.match),"Consensus table must reproduce exact BLOSUM score, endpoints and ties");
			check(model.scorePath(query,newer.start,newer.path)==newer.score,"Optimal trace must reproduce its DP score");
		}
		for(int trial=0; trial<300; trial++){
			final int n=1+random.nextInt(3); final byte[] query=sequence(random,1+random.nextInt(3));
			final int[][] table=new int[22][n];
			for(int a=0; a<22; a++){for(int j=0; j<n; j++){table[a][j]=(random.nextInt(16)-4)*64;}}
			final HbmPositionModel model=HbmPositionModel.fromScores(table);
			int best=Integer.MIN_VALUE;
			for(int start=0; start<=n; start++){best=Math.max(best,enumerate(query,table,0,start,0));}
			final HbmPositionModel.Result result=model.align(query,true);
			check(result.score==best,"Exhaustive path oracle disagrees with profile DP");
			check(model.scorePath(query,result.start,result.path)==best,"Profile path failed to reproduce exhaustive optimum");
			final int before=model.scoreAt((byte)0,0);table[0][0]++;
			check(model.scoreAt((byte)0,0)==before,"Caller matrix mutation must not change immutable model");
		}
		final byte a=Blosum62.encode(new byte[]{'A'},"A")[0],r=Blosum62.encode(new byte[]{'R'},"R")[0];
		final AAGraph graph=new AAGraph(new byte[]{a},0);graph.ref[0].add(r,1);
		final HbmPositionModel expected=HbmPositionModel.deriveOne(graph,"expected",0.1,false,bg);
		check(expected.scoreAt(a,0)==Math.round(32.0*(Blosum62.score(a,a)+Blosum62.score(a,r))),"Expected BLOSUM must average the two equally observed residues");
		final HbmPositionModel log=HbmPositionModel.deriveOne(graph,"logodds",0.1,false,bg);
		check(log.scoreAt(a,0)==Math.round(128*Math.log(9.1)/Math.log(2)),"Half-bit log odds and smoothing arithmetic differ");
		check(log.scoreAt(Blosum62.X_CODE,0)==0,"Unknown residue must add zero log-odds evidence");
		final byte absent=(byte)((a!=2 && r!=2) ? 2 : 3);
		check(HbmPositionModel.deriveOne(graph,"logodds",0.1,true,bg).scoreAt(absent,0)== -4*64,"Clipped unseen residue must reach the declared floor");
		final HbmPositionModel floor4=HbmPositionModel.deriveOne(graph,"logodds",0.01,true,bg,-4);
		final HbmPositionModel floor5=HbmPositionModel.deriveOne(graph,"logodds",0.01,true,bg,-5);
		check(floor5.scoreAt(absent,0)== -5*64,"Unseen residue at beta .01 must reach the requested -5 floor");
		check(floor5.scoreAt(a,0)==floor4.scoreAt(a,0),"Changing the floor must preserve positive emission scores");
		check(floor5.scoreAt(Blosum62.X_CODE,0)==0,"Changing the floor must preserve neutral unknown-residue scores");
		for(int residue=0; residue<20; residue++){
			check(floor4.scoreAt((byte)residue,0)==HbmPositionModel.deriveOne(graph,"logodds",0.01,true,bg).scoreAt((byte)residue,0),"Explicit -4 must preserve the original default");
		}
		for(double bad : new double[]{Double.NaN,Double.NEGATIVE_INFINITY,-101,1}){
			boolean failed=false;
			try{HbmPositionModel.deriveOne(graph,"logodds",0.01,true,bg,bad);}catch(IllegalArgumentException good){failed=true;}
			check(failed,"Invalid clipping floor must fail before building scores");
		}
		check(log.scorePath(new byte[]{a,a},0,new byte[]{'m','I'})==log.scoreAt(a,0)-HbmPositionModel.GAP,"Terminal insertions must be scored");
		Arrays.fill(graph.ref[0].count,0);
		check(HbmPositionModel.deriveOne(graph,"logodds",0.1,false,bg).scoreAt(a,0)==0,"Empty column must supply background, not infinite evidence");
		boolean rejected=false;try{log.align(new byte[]{20},true);}catch(IllegalArgumentException good){rejected=true;}
		check(rejected,"Encoded stop must fail loudly");
		System.err.println("POSITION_MODEL_TEST_PASS parity=400 exhaustive=300 emissions=PASS");
	}
	/** Enumerates paths directly, with no DP table or recurrence memoization. */
	private static int enumerate(byte[] query,int[][] table,int q,int r,int score){
		if(q==query.length){return r>0 ? score : Integer.MIN_VALUE;}
		int best=enumerate(query,table,q+1,r,score-HbmPositionModel.GAP);
		if(r<table[0].length){
			best=Math.max(best,enumerate(query,table,q+1,r+1,score+table[query[q]][r]));
			best=Math.max(best,enumerate(query,table,q,r+1,score-HbmPositionModel.GAP));
		}
		return best;
	}
	private static byte[] sequence(Random random,int length){
		final byte[] out=new byte[length];for(int i=0; i<length; i++){out[i]=(byte)random.nextInt(20);}return out;
	}
	private static void check(boolean condition,String reason){if(!condition){throw new AssertionError(reason);}}
}
