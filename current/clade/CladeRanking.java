package clade;

import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;
import dna.Data;
import tax.TaxTree;

/**
 * Ranking neural network: replaces compositeScore as the display ordering
 * for QuickClade results. Loaded from resources/ranking.bbnet; if absent,
 * compositeScore ordering is retained (graceful degradation).
 *
 * The net takes a 35-dim input vector per hit (same as RankingVectorizer)
 * and outputs a single score. Higher = better match. Ranking builds its vector
 * on the composite (heuristic) order it was trained on, then re-sorts the
 * display. Confidence is cached AFTERWARD, on this ranking-sorted order -- the
 * order the ranking-regenerated confidence bundle was trained on, and the order
 * whose ranking score its slow (49-dim) variant reads as an input dimension.
 *
 * Not used in fast mode (no DDL sketch data).
 *
 * @author Neptune
 * @date July 2026
 */
public class CladeRanking{

	/** Runtime on/off for the ranking net (set via useranking=t/f). ready() also requires a loaded net. */
	public static boolean enabled=true;
	public static boolean ready(){return net!=null && enabled;}

	/**
	 * Score each hit in the display list with the ranking net, then re-sort
	 * by ranking score (descending). Called AFTER the composite sort (so the
	 * ranking vector sees the heuristic order it was trained on) and BEFORE
	 * confidence caching (so confidence sees this ranking-sorted order, which
	 * its bundle was trained on, and its slow variant can read rankingScore).
	 *
	 * @param display The composite-sorted hit list for one query
	 * @param buckets The DDL bucket count (4096 or 32768)
	 */
	public static void scoreAndResort(ArrayList<Comparison> display, int buckets){
		if(!ready() || display.size()<2){return;}
		final CellNet cn=getThreadNet();
		final int n=display.size();
		final int top=Math.min(TOP_COUNT, n);
		final float m=display.get(0).composite;
		final float t=display.get(Math.min(TOP_COUNT-1, n-1)).composite;
		final float spread=m-t;
		final float log2b=(float)(Math.log(Math.max(buckets, 1))/LN2/16.0);

		// Precompute ancestor KEYS for sharing counts, from the lineage STRING (no TaxTree). Each key is the FULL
		// lineage prefix up to a SHARE_LEVEL, so homonyms (same rank name under different ancestors) never collide
		// -- matching the global uniqueness of the node IDs this replaces.
		final String[][] anc=new String[n][SHARE_LEVELS.length];
		for(int i=0; i<n; i++){
			final CharSequence lin=display.get(i).ref.lineage();
			for(int L=0; L<SHARE_LEVELS.length; L++){
				anc[i][L]=ancestorKey(lin, SHARE_LEVELS[L]);
			}
		}

		// Score each hit
		for(int i=0; i<n; i++){
			float[] v=buildVector(i, display, m, t, spread, anc, top, log2b);
			cn.applyInput(v);
			display.get(i).rankingScore=(float)cn.feedForward();
		}

		// Re-sort by ranking score
		display.sort((a, b)->Float.compare(b.rankingScore, a.rankingScore));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Vector building       ----------------*/
	/*--------------------------------------------------------------*/

	/** Build the 35-dim ranking vector for hit i. Same construction as RankingVectorizer.buildVector. */
	private static float[] buildVector(int i, ArrayList<Comparison> display,
			float m, float t, float spread, String[][] anc, int top, float log2b){
		final Comparison h=display.get(i);
		final float[] v=new float[VECTOR_DIMS];
		int p=0;

		// Main (26): same as ConfidenceVectorizer/RankingVectorizer dims 1-26
		final boolean spectraValid=spectraPresent(h);
		final boolean sketchValid=(h.kid>=0);
		final float len=Math.max(h.query.bases, 1);
		v[p++]=(float)(Math.log(len)/LN2);
		v[p++]=(float)(0.001*Math.sqrt(len));
		v[p++]=spectraValid ? r5(h.gcdif) : 0;
		v[p++]=spectraValid ? r5(h.strdif) : 0;
		v[p++]=spectraValid ? r5(h.hhdif) : 0;
		v[p++]=spectraValid ? r5(h.cagadif) : 0;
		v[p++]=spectraValid ? r5(h.k3dif) : 0;
		v[p++]=spectraValid ? r5(h.k4dif) : 0;
		v[p++]=spectraValid ? r5(h.k5dif) : 0;
		v[p++]=spectraValid ? r5(h.k3dif/(h.k5dif+0.01f)) : 0;
		v[p++]=sketchValid ? r4(h.kid) : 0;
		v[p++]=sketchValid ? r4(h.wkid) : 0;
		v[p++]=spectraValid ? 1 : 0;
		v[p++]=sketchValid ? 1 : 0;
		final float ssuANI=1-h.ssudif;
		final boolean ssuValid=(ssuANI>0);
		v[p++]=ssuValid ? r4(ssuANI) : 0;
		v[p++]=ssuValid ? 1 : 0;
		p=domainOneHot(h.ref.domain, v, p);
		final boolean lcaValid=(h.sketchLCA>=0);
		v[p++]=lcaValid ? encodeLevel9(h.sketchLCA) : 0;
		v[p++]=lcaValid ? 1 : 0;
		v[p++]=log2b;

		// Ranking block (9)
		final float denom=spread+EPS;
		float norm=(denom<=0 ? 0 : (h.composite-t)/denom);
		norm=Math.max(-CLAMP, Math.min(CLAMP, norm));
		v[p++]=norm;
		v[p++]=Math.min(1f, 0.1f*spread);
		v[p++]=0.1f*Math.min(display.size(), TOP_COUNT);
		for(int L=0; L<SHARE_LEVELS.length; L++){
			v[p++]=0.1f*sharingCount(i, L, anc, top);
		}
		v[p++]=(float)(0.1*Math.log(Math.max(h.ref.bases, 1))/LN2);
		final long card=(h.ref.ddl!=null ? h.ref.ddl.cardinality() : 0);
		v[p++]=(card>0 ? (float)(0.1*Math.log(card)/LN2) : 0);

		return v;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean spectraPresent(Comparison c){
		return c.k3dif<1 || c.k4dif<1 || c.k5dif<1;
	}

	private static float r5(float x){return Math.round(x*100000)/100000f;}
	private static float r4(float x){return Math.round(x*10000)/10000f;}

	private static float encodeLevel9(int level){
		float x=(11-level)/9f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	private static int sharingCount(int i, int L, String[][] anc, int top){
		final String self=anc[i][L];
		if(self==null){return 0;}
		int count=0;
		for(int j=0; j<top; j++){
			if(self.equals(anc[j][L])){count++;}
		}
		return count;
	}

	/** The FULL lineage prefix up to (and including) the token at `level` (s__/g__/f__/p__), or null when that
	 *  rank is absent from the lineage. Uses the whole path -- not the bare name -- so a homonym (same rank name
	 *  under a different ancestor) yields a different key, matching the global uniqueness of the tree node IDs
	 *  this replaces. Prefixes are anchored to a ';' token boundary so e.g. k__ cannot match inside sk__. */
	private static String ancestorKey(CharSequence lineage, int level){
		if(lineage==null){return null;}
		final String s=";"+lineage;
		final String prefix=";"+levelShort(level)+"__";
		final int pos=s.indexOf(prefix);
		if(pos<0){return null;}
		final int end=s.indexOf(';', pos+prefix.length());
		return end<0 ? s : s.substring(0, end);
	}

	/** Short rank-prefix letter for the SHARE_LEVELS, matching the lineage-string tokens. */
	private static String levelShort(int level){
		if(level==TaxTree.SPECIES){return "s";}
		if(level==TaxTree.GENUS){return "g";}
		if(level==TaxTree.FAMILY){return "f";}
		if(level==TaxTree.PHYLUM){return "p";}
		return null;
	}

	/** 7-way domain one-hot from the record's STORED domain category (0-6); all-zero when unknown. Embedded at
	 *  DB build (cladeloader adddomain=t) so no TaxTree is loaded at query time. Reproduces the tree-based 1-hot
	 *  the ranking net trained on (adddomain uses the identical superkingdom/kingdom logic). */
	private static int domainOneHot(int domain, float[] v, int p){
		for(int k=0; k<7; k++){v[p+k]=(k==domain ? 1 : 0);}
		return p+7;
	}

	private static CellNet getThreadNet(){
		CellNet cn=threadNet.get();
		if(cn==null){cn=net.copy(true); threadNet.set(cn);}
		return cn;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Static loading         ----------------*/
	/*--------------------------------------------------------------*/

	private static CellNet loadNet(){
		String path=Data.findPath("?ranking.bbnet", true);
		if(path==null){return null;}
		try{return CellNetParser.load(path);}
		catch(Exception e){System.err.println("WARNING: failed to load ranking net: "+e.getMessage()); return null;}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private static final CellNet net=loadNet();
	private static final ThreadLocal<CellNet> threadNet=new ThreadLocal<>();

	static final double LN2=Math.log(2);
	static final int VECTOR_DIMS=35;
	static final int TOP_COUNT=10;
	static final float CLAMP=2f;
	static final float EPS=0.02f;
	static final int[] SHARE_LEVELS={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY, TaxTree.PHYLUM};
}
