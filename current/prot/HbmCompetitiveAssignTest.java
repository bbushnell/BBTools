package prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

/** Exact input-repair and independent sorted-shortlist tests. @author Keqing */
public final class HbmCompetitiveAssignTest {
	public static void main(String[] args){
		final byte[] normal=bytes("ACDE*"), lowercase=bytes("acde");
		check(HbmCompetitiveAssign.trimEdges(normal, "normal")==normal, "Conventional terminal marker was rewritten");
		check(HbmCompetitiveAssign.trimEdges(lowercase, "lower")==lowercase, "Case was rewritten");
		check(Arrays.equals(HbmCompetitiveAssign.trimEdges(bytes("**ACDE**"), "edges"), bytes("ACDE")), "Malformed edge repair differs");
		check(HbmCompetitiveAssign.queryId("id|0 description").equals("id|0"), "Query ID tokenization differs");
		reject(()->HbmCompetitiveAssign.queryId(" blank"), "Empty query ID");
		reject(()->HbmCompetitiveAssign.trimEdges(bytes("***"), "empty"), "Empty query");
		reject(()->HbmCompetitiveAssign.trimEdges(bytes("AC*DE"), "internal"), "Internal query stop");
		boolean rejected=false;
		try{new ProteinSequence("unsupported", HbmCompetitiveAssign.trimEdges(bytes("ACUDE"), "unsupported"));}
		catch(RuntimeException expected){rejected=true;}
		check(rejected, "Unsupported residues were silently accepted");
		final Random random=new Random(20261010);
		for(int trial=0; trial<1000; trial++){
			final int n=50+random.nextInt(100); final float[] scores=new float[n]; final Integer[] order=new Integer[n];
			for(int i=0; i<n; i++){scores[i]=random.nextInt(13)-6; order[i]=i;}
			Arrays.sort(order, (a,b)->{final int c=Float.compare(scores[b], scores[a]); return c==0 ? Integer.compare(a,b) : c;});
			final int[] top=new int[50]; HbmCompetitiveAssign.top(scores, top);
			for(int i=0; i<top.length; i++){check(top[i]==order[i], "Top50 differs from complete independent sort");}
			final int limit=50+random.nextInt(n-49);
			final Integer[] prefix=new Integer[limit];
			for(int i=0; i<limit; i++){prefix[i]=i;}
			Arrays.sort(prefix, (a,b)->{final int c=Float.compare(scores[b], scores[a]); return c==0 ? Integer.compare(a,b) : c;});
			HbmCompetitiveAssign.top(scores, top, limit);
			for(int i=0; i<top.length; i++){check(top[i]==prefix[i], "Restricted top50 differs from independent prefix sort");}
		}
		final int[] restricted=new int[2];
		HbmCompetitiveAssign.top(new float[]{1, 3, 2, 100, 99}, restricted, 3);
		check(Arrays.equals(restricted, new int[]{1, 2}), "High-scoring excluded families displaced eligible candidates");
		reject(()->HbmCompetitiveAssign.top(new float[]{1, 2}, new int[2], 1), "Invalid shortlist dimensions");
		reject(()->HbmCompetitiveAssign.top(new float[]{1, 2}, new int[1], 3), "Invalid shortlist dimensions");
		reject(()->HbmCompetitiveAssign.top(new float[]{0, Float.NaN}, new int[1]), "Nonfinite");
		reject(()->HbmCompetitiveAssign.top(new float[]{Float.NEGATIVE_INFINITY, 1}, new int[1]), "Nonfinite");
		System.err.println("HBM_COMPETITIVE_INPUT_TEST_PASS shortlist_trials=1000 edge_bytes=true terminal_marker=true internal_rejected=true unsupported_rejected=true");
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean ok, String reason){if(!ok){throw new AssertionError(reason);}}
	private static void reject(Runnable action, String text){
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(text), "Wrong rejection: "+e); return;}
		throw new AssertionError("Invalid input accepted: "+text);
	}
}
