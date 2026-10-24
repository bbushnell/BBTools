package test.aligner;

import aligner.*;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import map.ObjectMap;
import structures.IntList;
import static aligner.CovarianceModel.*;

/** Independent derivation enumeration versus interval DP on tiny grammars.
 * Enumeration builds emitted strings rather than filling interval cells.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelCykTest {
	public static void main(String[] args){
		assert(args.length==0):"Exhaustive fixtures use explicit small grammars";
		final CovarianceModel linear=linear();compareAll(linear, 3);
		final CovarianceModel branch=branch();compareAll(branch, 3);
		final CovarianceModelCyk.Result b=new CovarianceModelCyk(branch, 10000).align(bytes("AT"));
		check(b.score==5 && b.modelPosition[0]==1 && b.modelPosition[1]==2, "Bifurcation must preserve left/right consensus order");
		check(!Float.isFinite(new CovarianceModelCyk(branch, 10000).align(bytes("A")).score), "Impossible lengths cannot acquire a trace");
		final CovarianceModel pair=pair();final CovarianceModelCyk scorer=new CovarianceModelCyk(pair, 10000);
		final CovarianceModelCyk.Result p=scorer.align(bytes("AT"));
		check(p.score==8 && p.matched==2 && p.inserted==0 && p.deleted==0, "Paired emission must use the joint score");
		check(p.modelPosition[0]==1 && p.modelPosition[1]==2, "MP must assign two ordered model coordinates");
		check(scorer.align(bytes("aU")).score==8, "Lowercase and T/U share canonical RNA coding");
		boolean unsupported=false;try{scorer.align(bytes("A?"));}catch(CovarianceModelCyk.UnsupportedSequenceException expected){unsupported=true;}
		check(unsupported, "Illegal input must be labeled unsupported, never silently substituted");
		boolean limited=false;try{new CovarianceModelCyk(pair, 1).align(bytes("AT"));}catch(IllegalArgumentException expected){limited=true;}
		check(limited, "Matrix sizing must fail before an excessive allocation");
		final CovarianceModel tie=tie();final CovarianceModelCyk.Result t=new CovarianceModelCyk(tie, 10000).align(bytes("A"));
		check(t.choice.get(0)==0 && t.states.get(1)==1, "Equal alternatives retain first child order");
		searchFixtures();
		System.out.println("CM_CYK_TEST_PASS exhaustiveStrings=170");
		check(searchChecks==458, "Both grammars must retain the declared exhaustive search population");
		System.out.println("CM_SEARCH_TEST_PASS exhaustiveSearches="+searchChecks);
	}
	private static void compareAll(CovarianceModel m, int length){
		assert(length<=3):"Derivation enumeration is deliberately exponential and bounded";
		final ArrayList<Derivation> all=enumerate(m, 0, length);final CovarianceModelCyk scorer=new CovarianceModelCyk(m, 100000);
		final CovarianceModelScoreOnly compact=new CovarianceModelScoreOnly(m, 100000);
		for(int n=0; n<=length; n++){for(int word=0; word<(1<<(2*n)); word++){
			final byte[] sequence=new byte[n];int code=word;for(int i=0; i<n; i++){sequence[i]=ALPHABET[code&3];code>>>=2;}
			final String text=new String(sequence, StandardCharsets.US_ASCII);float expected=Float.NEGATIVE_INFINITY;
			for(Derivation d:all){if(d.sequence.equals(text) && d.score>expected){expected=d.score;}}
			final CovarianceModelCyk.Result actual=scorer.align(sequence);
			check(actual.score==expected, "Enumeration disagrees for "+m.name+":"+text+" expected="+expected+" actual="+actual.score);
			check(compact.align(sequence).score==expected, "Reused triangular decks disagree with independent derivations for "+m.name+":"+text);
			if(Float.isFinite(actual.score)){check(actual.matched+actual.inserted==n, "Trace emissions must conserve even an empty input");}
			for(int limit=1; limit<=Math.max(1, n); limit++){compareSearch(scorer, sequence, text, all, limit);}
		}}
	}
	private static void compareSearch(CovarianceModelCyk scorer, byte[] sequence, String text, ArrayList<Derivation> all, int limit){
		assert(limit>0):"Interval-length bounds belong to the search comparison, not CM emissions";
		float best=Float.NEGATIVE_INFINITY;final IntList starts=new IntList(), stops=new IntList();
		for(int stop=1; stop<=sequence.length; stop++){for(int width=1; width<=Math.min(stop, limit); width++){
			final int start=stop-width;final String sub=text.substring(start, stop);float score=Float.NEGATIVE_INFINITY;
			for(Derivation d:all){if(d.sequence.equals(sub) && d.score>score){score=d.score;}}
			if(!Float.isFinite(score) || score<best){continue;}if(score>best){best=score;starts.clear();stops.clear();}
			starts.add(start+1);stops.add(stop);
		}}
		final CovarianceModelCyk.Result r=scorer.search(sequence, limit);searchChecks++;
		check(r.score==best, "Search differs from independent derivation enumeration: "+text+" limit="+limit);
		check(r.rootStarts.size==starts.size && r.rootStops.size==stops.size, "Every equal best root interval must be visible");
		for(int i=0; i<starts.size; i++){check(r.rootStarts.get(i)==starts.get(i) && r.rootStops.get(i)==stops.get(i), "Tie coordinates/order differ from exact enumeration");}
		if(Float.isFinite(best)){check(r.targetFrom==starts.get(0) && r.targetTo==stops.get(0), "Selected interval must be the first reported tie");
			check(r.matched+r.inserted==r.targetTo-r.targetFrom+1, "Selected trace must consume only its chosen interval");
			for(int i=1; i<=sequence.length; i++){check((r.stateAtBase[i-1]>=0)==(i>=r.targetFrom && i<=r.targetTo), "Search cannot assign a model state to an unaligned flank");}}
		else{check(r.targetFrom==0 && r.targetTo==0 && r.states.size==0, "Impossible search cannot invent endpoints or a trace");}
	}
	private static void searchFixtures(){
		final CovarianceModel m=pair();m.emissionScore[1][1]=9;final CovarianceModelCyk scorer=new CovarianceModelCyk(m, 100000);
		final CovarianceModelCyk.Result plus=scorer.search(bytes("TTACGG"), 6);
		check(plus.score==9 && plus.targetFrom==3 && plus.targetTo==4, "Embedded plus-strand pair must exclude both flanks");
		final byte[] original=bytes("TTGTGG");final CovarianceModelCyk.Result minus=scorer.search(CovarianceModelSearch.reverseComplement(original), 6);
		check(minus.score==9 && CovarianceModelSearch.forwardPosition(minus.targetTo, 6, true)==3 && CovarianceModelSearch.forwardPosition(minus.targetFrom, 6, true)==4, "Minus-strand bounds must map to original forward coordinates");
		for(int i=1; i<=6; i++){check(CovarianceModelSearch.forwardPosition(CovarianceModelSearch.forwardPosition(i, 6, true), 6, true)==i, "Strand conversion must be an involution at both endpoints and every interior base");}
		check(Arrays.equals(CovarianceModelSearch.reverseComplement(CovarianceModelSearch.reverseComplement(original)), original), "Reverse complement must preserve the original canonical sequence on round trip");
		final CovarianceModelCyk.Result ties=scorer.search(bytes("ACAC"), 4);
		check(ties.score==9 && ties.rootStarts.size==2 && ties.rootStarts.get(0)==1 && ties.rootStarts.get(1)==3, "Repeated exact copies must retain both tied best intervals");
		check(!Float.isFinite(scorer.search(bytes("AC"), 1).score), "Explicit maximum hit length must be enforced");
		boolean rejected=false;try{scorer.search(bytes("?AC"), 3);}catch(CovarianceModelCyk.UnsupportedSequenceException expected){rejected=true;}
		check(rejected, "Illegal flanks remain explicitly unsupported, even when a canonical subinterval exists");
	}
	private static ArrayList<Derivation> enumerate(CovarianceModel m, int v, int budget){
		assert(budget>=0):"Generative enumeration tracks remaining emitted length";
		final ArrayList<Derivation> out=new ArrayList<Derivation>();final byte type=m.type[v];
		if(type==E){out.add(new Derivation("", 0));return out;}
		if(type==B){final ArrayList<Derivation> left=enumerate(m, m.childFirst[v], budget), right=enumerate(m, m.childCount[v], budget);
			for(Derivation a:left){for(Derivation b:right){if(a.sequence.length()+b.sequence.length()<=budget){out.add(new Derivation(a.sequence+b.sequence, a.score+b.score));}}}return out;}
		final int emitted=type==MP ? 2 : (type==ML || type==MR || type==IL || type==IR ? 1 : 0);if(emitted>budget){return out;}
		for(int c=0; c<m.childCount[v]; c++){final float trans=m.transitionScore[v][c];if(!Float.isFinite(trans)){continue;}
			for(Derivation child:enumerate(m, m.childFirst[v]+c, budget-emitted)){
				if(emitted==0){out.add(new Derivation(child.sequence, child.score+trans));continue;}
				for(int a=0; a<4; a++){if(type==MP){for(int b=0; b<4; b++){
					final float score=child.score+trans+m.emissionScore[v][a*4+b];if(Float.isFinite(score)){out.add(new Derivation((char)ALPHABET[a]+child.sequence+(char)ALPHABET[b], score));}}}
					else{final String sequence=type==ML || type==IL ? (char)ALPHABET[a]+child.sequence : child.sequence+(char)ALPHABET[a];
						final float score=child.score+trans+m.emissionScore[v][a];if(Float.isFinite(score)){out.add(new Derivation(sequence, score));}}}
			}
		}return out;
	}
	static CovarianceModel linear(){
		final CovarianceModel m=model("all_states", 10, 3, 2);
		state(m, 0, S, 0, 1, new float[]{-4,-4,-1,-3,-3,-5});
		state(m, 1, IL, 0, 1, new float[]{-2,-4,-1,-3,-3,-5});
		state(m, 2, IR, 0, 2, new float[]{-2,-1,-3,-3,-5});
		state(m, 3, MP, 1, 7, new float[]{-4,-4,0});state(m, 4, ML, 1, 7, new float[]{-4,-4,0});
		state(m, 5, MR, 1, 7, new float[]{-4,-4,0});state(m, 6, D, 1, 7, new float[]{-4,-4,0});
		state(m, 7, IL, 1, 7, new float[]{-2,-4,0});state(m, 8, IR, 1, 8, new float[]{-2,0});state(m, 9, E, 2, -1, new float[0]);
		m.consensusLeft[0]=0;m.consensusRight[0]=3;m.consensusLeft[1]=1;m.consensusRight[1]=2;
		m.emissionScore[3][3]=8;m.emissionScore[4][0]=2;m.emissionScore[5][3]=3;return m;
	}
	static CovarianceModel branch(){
		final CovarianceModel m=model("bifurcation", 6, 6, 2);state(m, 0, S, 0, 1, new float[]{0});state(m, 1, B, 1, 2, new float[0]);m.childCount[1]=4;
		state(m, 2, ML, 2, 3, new float[]{0});state(m, 3, E, 3, -1, new float[0]);state(m, 4, MR, 4, 5, new float[]{0});state(m, 5, E, 5, -1, new float[0]);
		m.consensusLeft[2]=1;m.consensusRight[4]=2;m.emissionScore[2][0]=2;m.emissionScore[4][3]=3;return m;
	}
	static CovarianceModel pair(){
		final CovarianceModel m=model("pair", 3, 3, 2);state(m, 0, S, 0, 1, new float[]{0});state(m, 1, MP, 1, 2, new float[]{0});state(m, 2, E, 2, -1, new float[0]);
		m.consensusLeft[1]=1;m.consensusRight[1]=2;m.emissionScore[1][3]=8;return m;
	}
	static CovarianceModel tie(){
		final CovarianceModel m=model("tie", 5, 5, 1);state(m, 0, S, 0, 1, new float[]{0,0});state(m, 1, D, 1, 3, new float[]{0});state(m, 2, D, 2, 3, new float[]{0});
		state(m, 3, ML, 3, 4, new float[]{0});state(m, 4, E, 4, -1, new float[0]);m.consensusLeft[3]=1;return m;
	}
	private static CovarianceModel model(String name, int states, int nodes, int clen){
		assert(states>=nodes):"Test grammar dimensions must fit its explicit node mapping";
		final ObjectMap<String,String> headers=new ObjectMap<String,String>(String.class, String.class);
		headers.put("QDBBETA1", "1e-7");headers.put("QDBBETA2", "1e-15");
		return new CovarianceModel(name, "TEST", clen, 20, states, nodes, headers, new float[]{.25f,.25f,.25f,.25f});
	}
	private static void state(CovarianceModel m, int v, byte type, int node, int child, float[] scores){
		assert(v<m.states() && node<m.nodes()):"Fixture states must stay inside declared dimensions";
		m.type[v]=type;m.node[v]=node;m.childFirst[v]=child;m.childCount[v]=scores.length;m.transitionScore[v]=scores;
		m.qdb[v]=type==E ? new int[]{0,0,0,0} : new int[]{0,0,20,20};
		final int n=type==MP ? 16 : (type==ML || type==MR || type==IL || type==IR ? 4 : 0);m.emissionScore[v]=new float[n];Arrays.fill(m.emissionScore[v], -1);
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final class Derivation {Derivation(String s,float x){sequence=s;score=x;}final String sequence;final float score;}
	private static final byte[] ALPHABET={'A','C','G','T'};
	private static int searchChecks=0;
}
