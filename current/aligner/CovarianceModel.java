package aligner;

import map.ObjectMap;

/** Parsed INFERNAL1/a model. Arrays are package-private and read-only after parsing.
 * Scores are reconstructed from normalized probabilities, not copied from rounded text.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModel {

	CovarianceModel(String name_, String accession_, int clen_, int window_, int states, int nodes,
			ObjectMap<String,String> headers_, float[] null_){
		assert(states>0 && nodes>0 && clen_>0):"CM arrays require positive declared dimensions";
		name=name_;accession=accession_;clen=clen_;window=window_;headers=headers_;nullModel=null_;
		type=new byte[states];node=new int[states];parentLast=new int[states];parentCount=new int[states];
		childFirst=new int[states];childCount=new int[states];qdb=new int[states][4];
		transition=new float[states][];emission=new float[states][];
		transitionScore=new float[states][];emissionScore=new float[states][];
		nodeType=new byte[nodes];nodeState=new int[nodes];mapLeft=new int[nodes];mapRight=new int[nodes];
		consLeft=new byte[nodes];consRight=new byte[nodes];rfLeft=new byte[nodes];rfRight=new byte[nodes];
		consensusLeft=new int[nodes];consensusRight=new int[nodes];
	}

	public int states(){return type.length;}
	public int nodes(){return nodeType.length;}
	public String header(String key){return headers.get(key);}
	/** Infernal reserves the state after all stored CM states for local-end EL. */
	byte traceType(int state){
		if(state<0 || state>states()){throw new IllegalArgumentException("Trace state leaves the CM and its virtual EL state: "+state);}
		return state==states() ? EL : type[state];
	}
	public String traceStateName(int state){return STATE_NAMES[traceType(state)];}
	public int stateCount(String state){
		assert(state!=null):"State counting requires an explicit type name";
		int n=0;for(byte t:type){if(STATE_NAMES[t].equals(state)){n++;}}return n;
	}

	public final String name, accession;
	public final int clen, window;
	private final ObjectMap<String,String> headers;
	final float[] nullModel;
	final byte[] type, nodeType, consLeft, consRight, rfLeft, rfRight;
	final int[] node, parentLast, parentCount, childFirst, childCount, nodeState;
	final int[] mapLeft, mapRight, consensusLeft, consensusRight;
	final int[][] qdb;
	final float[][] transition, emission, transitionScore, emissionScore;
	static final byte D=0, MP=1, ML=2, MR=3, IL=4, IR=5, S=6, E=7, B=8, EL=9;
	static final byte BIF=0, MATP=1, MATL=2, MATR=3, BEGL=4, BEGR=5, ROOT=6, END=7;
	static final String[] STATE_NAMES={"D", "MP", "ML", "MR", "IL", "IR", "S", "E", "B", "EL"};
	static final String[] NODE_NAMES={"BIF", "MATP", "MATL", "MATR", "BEGL", "BEGR", "ROOT", "END"};
}
