package aligner;

import java.util.Arrays;
import static aligner.CovarianceModel.*;

/** Immutable scoring configuration over a parsed global CM.
 * Local setup follows Infernal1.1.5 cm_modelconfig.c:397-543 and CMLogoddsify.
 * Topology/emissions are unchanged; no parser-owned probability array is mutated.
 * This configures local model boundaries, not target search or truncation.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelConfiguration {
	public CovarianceModelConfiguration(CovarianceModel model_, boolean local_){
		require(model_!=null, "Scoring configuration requires the actual parsed model");
		model=model_;local=local_;final int states=model.states();
		beginScore=new float[states];endScore=new float[states];
		Arrays.fill(beginScore, Float.NEGATIVE_INFINITY);Arrays.fill(endScore, Float.NEGATIVE_INFINITY);
		if(!local){transitionScore=model.transitionScore;elSelf=Float.NEGATIVE_INFINITY;pBegin=pEnd=0;return;}
		require(model.nodes()>1 && model.nodeType[0]==ROOT, "Local entry requires the CM root and its first model node");
		pBegin=probability(model, "PBEGIN", .05f);pEnd=probability(model, "PEND", .05f);
		final String el=model.header("ELSELF");require(el!=null, "Local-end length cost must come from the required ELSELF header");
		elSelf=Float.parseFloat(el);require(Float.isFinite(elSelf) && elSelf<=0, "ELSELF is a finite log probability, never a positive reward");
		transitionScore=new float[states][];
		for(int v=0; v<states; v++){
			require(model.transitionScore[v]!=null && model.transition[v]!=null && model.transition[v].length==model.transitionScore[v].length,
				"Local rescaling needs the parsed probabilities corresponding to each transition score, state="+v);
			transitionScore[v]=model.transitionScore[v].clone();
		}
		int starts=0, exits=0;
		for(int n=2; n<model.nodes(); n++){if(canBegin(model.nodeType[n])){starts++;}}
		for(int n=1; n<model.nodes(); n++){if(canEnd(model, n)){exits++;}}
		beginScore[principal(model, 1)]=bits(1-pBegin);
		for(int n=2; n<model.nodes(); n++){if(canBegin(model.nodeType[n])){beginScore[principal(model, n)]=bits(pBegin/starts);}}
		Arrays.fill(transitionScore[0], Float.NEGATIVE_INFINITY);
		for(int n=1; n<model.nodes(); n++){
			if(!canEnd(model, n)){continue;}
			final int v=principal(model, n);final float end=pEnd/exits;endScore[v]=bits(end);
			float sum=0;for(float p:model.transition[v]){require(Float.isFinite(p) && p>=0, "Local transition rescaling cannot repair invalid probabilities");sum+=p;}
			final float denominator=sum+end;require(Float.isFinite(denominator) && denominator>0, "Eligible local exits need a positive outgoing probability mass");
			final float scale=(float)(1.0/denominator);
			// Infernal rescales ordinary transitions but does NOT rescale end[v].
			// Normalizing both here would implement a different scoring model.
			for(int c=0; c<transitionScore[v].length; c++){transitionScore[v][c]=bits(model.transition[v][c]*scale);}
		}
	}
	private static boolean canBegin(byte node){return node==MATP || node==MATL || node==MATR || node==BIF;}
	private static boolean canEnd(CovarianceModel model, int n){
		assert(n>0 && n<model.nodes()):"Local exits are considered only after the root node";
		final byte node=model.nodeType[n];
		return (node==MATP || node==MATL || node==MATR || node==BEGL || node==BEGR)
			&& n+1<model.nodes() && model.nodeType[n+1]!=END;
	}
	private static int principal(CovarianceModel model, int n){
		final int v=model.nodeState[n];require(v>=0 && v<model.states() && model.node[v]==n,
			"Local begin/end probabilities belong only to the declared principal state of node "+n);return v;
	}
	private static float probability(CovarianceModel model, String key, float fallback){
		final String text=model.header(key);final float p=text==null ? fallback : Float.parseFloat(text);
		require(Float.isFinite(p) && p>=0 && p<=1, key+" must be a finite probability mass in [0,1]");return p;
	}
	private static float bits(float p){
		assert(Float.isFinite(p) && p>=0 && p<=1):"CMLogoddsify converts probability mass, not an unnormalized score";
		return p==0 ? Float.NEGATIVE_INFINITY : (float)(Math.log(p)/Math.log(2));
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final CovarianceModel model;
	public final boolean local;
	public final float elSelf, pBegin, pEnd;
	final float[][] transitionScore;
	final float[] beginScore, endScore;
}
