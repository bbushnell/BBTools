package clade;

import java.util.ArrayList;

import ml.CellNet;
import dna.Data;
import shared.Tools;
import structures.FloatList;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Estimates the probability that a QuickClade hit is taxonomically correct at a
 * given level.
 *
 * Three models are supported, selected automatically from what
 * confidence.bbnets.gz actually contains:
 *
 * 1. V2 (preferred) -- one network per level for 9 levels (species..domain,
 *    including superkingdom), taking the 48-dimension feature vector defined by
 *    {@link ConfidenceVectorizer}. Because two of those dimensions describe OTHER
 *    hits of the same query (the strongest competitor, and the most taxonomically
 *    remote hit among the top ten), a V2 score cannot be computed from a lone
 *    Comparison -- it needs the query's hit group. See profile().
 * 2. V1 (legacy) -- the 8-level, 10-feature hybrid: all-length networks for
 *    order..domain and per-length-bin networks for species..family.
 * 3. No model file -- the 5-parameter PROB_K4/PROB_K5 sigmoid tables.
 *
 * Calibration converts a network's raw score into a probability. Two forms ship
 * in the bundle: the 4-parameter K*sigmoid(a*logit(x)+b)^c, and an isotonic
 * lookup table. The table is PREFERRED when present -- measured on held-out data
 * it beat the parametric form at every level (ECE 7x-27x lower), because the
 * parametric failure is structured rather than noisy: entire confidence bands
 * come out biased in one direction, which a scalar ECE hides when most of the
 * mass sits in the two extreme bins. The parametric constants are retained as a
 * fallback and for bundles that carry no table.
 *
 * @author Brian Bushnell, Ady, Noire
 */
public class CladeConfidence {

	/*--------------------------------------------------------------*/
	/*----------------          V2 entry point      ----------------*/
	/*--------------------------------------------------------------*/

	/** True when a V2 bundle is loaded and complete (per-level nets OR a single multi-output net). */
	public static boolean v2Ready(){return v2 || multiV2 || multiV2Slow;}

	/**
	 * Probability that group.get(index) is correct, for all 9 levels in LEVELS
	 * order (species first). Returns null when no V2 model is loaded.
	 *
	 * The 48-dim vector is built ONCE and reused across the 9 networks; the
	 * taxonomy lookups it needs (domain one-hot, pairwise LCAs across the top
	 * hits) dominate its cost, so building it per level would be wasteful.
	 */
	public static float[] profile(ArrayList<Comparison> group, int index){
		if((!v2 && !multiV2 && !multiV2Slow) || group==null || index<0 || index>=group.size()){return null;}
		final Comparison self=group.get(index);
		if(self.query==null || self.ref==null){return null;}
		final FloatList fl=getThreadInput();
		buildVector(fl, group, index);
		//Slow (49-dim) net: used ONLY when ranking actually ran for this hit -- slow mode sets rankingScore,
		//fast mode leaves it -Inf. Append the raw ranking score as dim 49, matching ConfidenceVectorizer,
		//whose training dim came through the 6-decimal RankingScore column, so round to 6 here too.
		if(multiV2Slow && Float.isFinite(self.rankingScore)){
			fl.add(r6(self.rankingScore));
			final CellNet net=getThreadMultiNetSlow();
			net.applyInput(fl);
			net.feedForward();
			final float[] raw=net.getOutput();
			final float[] out=new float[levelsMultiSlow.length];
			for(int i=0; i<levelsMultiSlow.length; i++){out[i]=applyMultiCalSlow(raw[i], i);}
			return out;
		}
		if(multiV2){ //ONE net, N output heads, one forward pass, N per-output calibrations
			final CellNet net=getThreadMultiNet();
			net.applyInput(fl);
			net.feedForward();
			final float[] raw=net.getOutput();
			final float[] out=new float[levelsMulti.length];
			for(int i=0; i<levelsMulti.length; i++){out[i]=applyMultiCal(raw[i], i);}
			return out;
		}
		//No slow-net hit and no 48-dim multi net -> use the per-level V2 nets if present, else return null so
		//the caller (cacheConfidence) takes its legacy/sigmoid fallback. Without this, a slow-only bundle would
		//return a non-null EMPTY array for a fast-mode hit, suppressing that fallback and crashing downstream.
		if(!v2){return null;}
		final float[] out=new float[levelsV2.length];
		for(int i=0; i<levelsV2.length; i++){
			final CellNet net=getThreadNetV2(i);
			if(net==null){out[i]=-1; continue;}
			net.applyInput(fl);
			out[i]=applyCal(net.feedForward(), i);
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------      V2 feature assembly     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds the 48-dim vector for group.get(i). This MUST stay dimension-for-
	 * dimension identical to ConfidenceVectorizer.buildVector/smallVector,
	 * including the invalid-fill=0 convention -- that class defines the training
	 * distribution and is the authority here.
	 *
	 * Values are rounded to the precision appendResultMachine printed, because
	 * the training vectors were parsed back out of that text: the difs went
	 * through 5 decimals and ssuID/kid/wkid through 4. Feeding full precision at
	 * inference would present the networks a subtly different distribution than
	 * the one they were fit on.
	 */
	private static void buildVector(FloatList fl, ArrayList<Comparison> group, int i){
		final Comparison h=group.get(i);
		fl.size=0;

		final float gcdif=r5(h.gcdif), strdif=r5(h.strdif), hhdif=r5(h.hhdif), cagadif=r5(h.cagadif);
		final float k3dif=r5(h.k3dif), k4dif=r5(h.k4dif), k5dif=r5(h.k5dif);
		//The validity tests MUST run on the ROUNDED values, not the raw ones. Training read
		//these back from printed text, so a dif of 0.999996 reached the vectorizer as 1.00000
		//and counted as "spectra absent". Testing the raw value here would call it present and
		//fill dims 3-10 that the training distribution has as zeros -- an eight-dimension
		//divergence on exactly the rows sitting against the 1.0 boundary.
		final boolean spectraValid=(k3dif<1 || k4dif<1 || k5dif<1);
		final boolean sketchValid=(h.kid>=0);

		//--- main (26) ---
		final float len=Math.max(h.query.bases, 1);
		fl.add((float)(Math.log(len)*INV_LN2));            //1  log2(len)
		fl.add((float)(0.001*Math.sqrt(len)));             //2  0.001*sqrt(len)
		fl.add(spectraValid ? gcdif : 0);                  //3
		fl.add(spectraValid ? strdif : 0);                 //4
		fl.add(spectraValid ? hhdif : 0);                  //5
		fl.add(spectraValid ? cagadif : 0);                //6
		fl.add(spectraValid ? k3dif : 0);                  //7
		fl.add(spectraValid ? k4dif : 0);                  //8
		fl.add(spectraValid ? k5dif : 0);                  //9
		fl.add(spectraValid ? k3dif/(k5dif+0.01f) : 0);    //10
		fl.add(sketchValid ? r4(h.kid) : 0);               //11 kid
		fl.add(sketchValid ? r4(h.wkid) : 0);              //12 wkid
		fl.add(spectraValid ? 1 : 0);                      //13 spectra-present
		fl.add(sketchValid ? 1 : 0);                       //14 sketch-present
		final float ssuID=r4(1f-h.ssudif);
		final boolean ssuValid=(ssuID>0);
		fl.add(ssuValid ? ssuID : 0);                      //15 SSU-ANI
		fl.add(ssuValid ? 1 : 0);                          //16 SSU-valid
		domainOneHot(h.ref.domain, fl);                    //17-23 domain 1-hot (7); reads the stored domain (no tree)
		final boolean lcaValid=(h.sketchLCA>0);
		fl.add(lcaValid ? encodeLevel(h.sketchLCA) : 0);   //24 cross-method LCA
		fl.add(lcaValid ? 1 : 0);                          //25 cross-method LCA valid
		fl.add(log2Buckets());                             //26 log2(#buckets)/16

		//--- small #1: strongest competitor; small #2: most remote among the top ten ---
		smallVector(altTopOther(i, group), h, group, fl);
		smallVector(altMostRemote(i, group), h, group, fl);
		assert(fl.size==VECTOR_DIMS) : fl.size+" != "+VECTOR_DIMS;
	}

	/** 11-dim alt-hit block, all zero when the alt does not exist. */
	private static void smallVector(int alt, Comparison self, ArrayList<Comparison> group, FloatList fl){
		if(alt<0){
			for(int k=0; k<SMALL_DIMS; k++){fl.add(0);}
			return;
		}
		final Comparison a=group.get(alt);
		final boolean aSketch=(a.kid>=0);
		fl.add(1);                                         //1 valid
		fl.add(r5(a.k4dif));                               //2 k4dif
		fl.add(aSketch ? r4(a.kid) : 0);                   //3 kid
		fl.add(aSketch ? r4(a.wkid) : 0);                  //4 wkid
		fl.add(aSketch ? 1 : 0);                           //5 sketch-present(alt)
		fl.add(encodeLevel(lcaLevel(a, self)));            //6 LCA(alt,this) -- always valid
		final boolean altLcaValid=(a.sketchLCA>0);
		fl.add(altLcaValid ? encodeLevel(a.sketchLCA) : 0);//7 alt's cross-method LCA
		fl.add(altLcaValid ? 1 : 0);                       //8
		final float aSsu=r4(1f-a.ssudif);
		final boolean altSsuValid=(aSsu>0);
		fl.add(altSsuValid ? aSsu : 0);                    //9 SSU-ANI(alt)
		fl.add(altSsuValid ? 1 : 0);                       //10
		fl.add(r5(a.gcdif));                               //11 gcdif(query, alt-ref)
	}

	/** Top hit other than i: hit 0, or hit 1 when i is itself hit 0. -1 if the query has one hit. */
	private static int altTopOther(int i, ArrayList<Comparison> group){
		if(group.size()<2){return -1;}
		return (i==0) ? 1 : 0;
	}

	/** Among the top TOP_EXAMINE hits, the one whose LCA with hit i is shallowest (most remote);
	 *  ties resolve to the highest-ranked, which is the strongest. */
	private static int altMostRemote(int i, ArrayList<Comparison> group){
		final Comparison self=group.get(i);
		int best=-1, bestLevel=-1;
		final int examine=Math.min(TOP_EXAMINE, group.size());
		for(int j=0; j<examine; j++){
			if(j==i){continue;}
			final int lvl=lcaLevel(group.get(j), self);
			if(lvl>bestLevel){bestLevel=lvl; best=j;} //higher level = shallower = more remote
		}
		return best;
	}

	/** Formal-rank level of the LCA between two hits' references, from their lineage strings (no TaxTree). The
	 *  lineages now carry d__ (domain), so cross-kingdom-same-domain pairs resolve correctly. Fallback: viruses
	 *  carry no d__ token (not a domain in NCBI), so when the lineages share no rank, use the stored domain field
	 *  -- same domain -> domain-level LCA, else unrelated (LIFE). Reproduces the tree-based LCA the net trained on. */
	private static int lcaLevel(Comparison a, Comparison b){
		if(a.ref==null || b.ref==null){return TaxTree.LIFE;}
		final int lca=CladeIndex.lineageLCA(a.ref.lineage(), b.ref.lineage());
		if(lca>0){return lca;}
		return (a.ref.domain>=0 && a.ref.domain==b.ref.domain) ? TaxTree.DOMAIN : TaxTree.LIFE;
	}

	/** Promotes a node up the parent chain to the next formal rank (>=SPECIES). */
	private static int formalLevel(int tid, TaxTree tree){
		if(tid<1 || tree==null){return TaxTree.LIFE;}
		TaxNode n=tree.getNode(tid);
		while(n!=null && n.level<TaxTree.SPECIES){n=(n.pid>0 && n.pid!=n.id) ? tree.getNode(n.pid) : null;}
		return n==null ? TaxTree.LIFE : n.level;
	}

	/** Continuous 0-1 LCA encoding: life(11)->0 ... species(2)->1.0, constant step 1/9. */
	private static float encodeLevel(int level){
		final float x=(11-level)/9f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	/** 7-way domain one-hot from the record's STORED domain category (0-6): bacteria, archaea, virus, animal,
	 *  plant, fungi, other-euk; all-zero when unknown (domain<0 or >6). The domain is embedded at DB build time
	 *  (cladeloader adddomain=t, which computes it from the tree ONCE), so no TaxTree is loaded at query time.
	 *  All-zero-when-unknown reproduces the training distribution exactly (the tree-based builder gave all-zero
	 *  for an unclassified superkingdom). */
	private static void domainOneHot(int domain, FloatList fl){
		for(int k=0; k<7; k++){fl.add(k==domain ? 1 : 0);}
	}

	/** log2(bucket count)/16, matching ConfidenceVectorizer's buckets= parameter. 0 in no-sketch (fast) mode. */
	private static float log2Buckets(){
		final int b=CladeIndex.USE_SKETCH_INDEX ? Clade.DDL_BUCKETS : 0;
		return (float)(Math.log(Math.max(b, 1))/Math.log(2)/16.0);
	}

	/** Round to the 5 decimals appendResultMachine printed for the dif columns. */
	private static float r5(float f){return Math.round(f*100000f)/100000f;}
	/** Round to the 4 decimals appendResultMachine printed for ssuID/kid/wkid. */
	private static float r4(float f){return Math.round(f*10000f)/10000f;}
	/** Round to the 6 decimals appendResultMachine printed for the RankingScore column. */
	private static float r6(float f){return Math.round(f*1000000f)/1000000f;}

	/*--------------------------------------------------------------*/
	/*----------------           Calibration        ----------------*/
	/*--------------------------------------------------------------*/

	/** Isotonic table when the bundle carries one, else the 4-parameter form. */
	private static float applyCal(float raw, int idx){
		if(loaded!=null && loaded.allLenLutX!=null && loaded.allLenLutX[idx]!=null){
			return lookup(raw, loaded.allLenLutX[idx], loaded.allLenLutY[idx]);
		}
		return calibrate(raw, loaded.allLenCal[idx]);
	}

	/** Monotone table lookup, linearly interpolated between knots and clamped at the ends. */
	static float lookup(float x, float[] lx, float[] ly){
		final int last=lx.length-1;
		if(x<=lx[0]){return ly[0];}
		if(x>=lx[last]){return ly[last];}
		int lo=0, hi=last;
		while(lo+1<hi){
			final int mid=(lo+hi)>>>1;
			if(lx[mid]<=x){lo=mid;}else{hi=mid;}
		}
		final float dx=lx[hi]-lx[lo];
		final float f=(dx<=0 ? 0 : (x-lx[lo])/dx);
		return ly[lo]+f*(ly[hi]-ly[lo]);
	}

	/** p = K * sigmoid(a * logit(x) + b) ^ c */
	private static float calibrate(float x, float[] p) {
		x=Tools.mid(0.0001f, x, 0.9999f);
		float K=p[0], a=p[1], b=p[2], c=p[3];
		double lx=Math.log(x / (1.0 - x));
		double s=1.0 / (1.0 + Math.exp(-(a * lx + b)));
		return (float)Math.min(1.0, K * Math.pow(s, c));
	}

	/*--------------------------------------------------------------*/
	/*----------------      V1 / fallback entries   ----------------*/
	/*--------------------------------------------------------------*/

	public static float probCorrect(int length, float gcdif, float strdif,
			float hhdif, float cagadif,
			float k3dif, float k4dif, float k5dif, int taxLevel) {
		int idx = levelToIndexV1(taxLevel);
		if(idx<0){return -1;}
		//A V2 or multi-output bundle carries no V1-shaped allLenNets (the multi net lives in
		//loaded.multiNet), and a lone hit cannot supply the alt-hit dimensions, so this legacy path
		//degrades to the sigmoid tables rather than dereference a null allLenNets. The real V2/multi
		//confidence is computed group-aware via profile() and cached upstream.
		if(!v2 && !multiV2 && !multiV2Slow && USE_NN && loaded!=null && loaded.allLenNets!=null && k5dif<1.0f){
			CellNet net;
			float[] cal;
			if(loaded.allLenNets[idx]!=null){
				net=loaded.allLenNets[idx];
				cal=loaded.allLenCal[idx];
			}else{
				int bin=binIndex(length);
				net=loaded.binNets[idx][bin];
				cal=loaded.binCal[idx][bin];
			}
			if(net!=null){
				return predictNN(net, cal, idx, length, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif);
			}
		}
		return predictSigmoid(idx, length, k3dif, k4dif, k5dif);
	}

	/** Backwards-compatible call without gcdif/strdif/hhdif/cagadif; always uses sigmoid. */
	public static float probCorrect(int length, float k3dif, float k4dif, float k5dif, int taxLevel) {
		int idx = levelToIndexV1(taxLevel);
		if(idx<0){return -1;}
		return predictSigmoid(idx, length, k3dif, k4dif, k5dif);
	}

	private static float predictNN(CellNet master, float[] cal, int idx, int length,
			float gcdif, float strdif, float hhdif, float cagadif,
			float k3dif, float k4dif, float k5dif) {
		CellNet net = getThreadNet(master, idx, length);
		FloatList fl = getThreadInput();
		fl.size=0;
		fl.add((float)(Math.log(Math.max(length, 1)) * INV_LN2));
		fl.add((float)(0.001 * Math.sqrt(Math.max(length, 1))));
		fl.add(gcdif);
		fl.add(strdif);
		fl.add(hhdif);
		fl.add(cagadif);
		fl.add(k3dif);
		fl.add(k4dif);
		fl.add(k5dif);
		fl.add((float)(k3dif / (k5dif + 0.01)));
		net.applyInput(fl);
		float raw = net.feedForward();
		return calibrate(raw, cal);
	}

	private static float predictSigmoid(int idx, int length, float k3dif, float k4dif, float k5dif) {
		double[] params = (k5dif < 1.0f) ? PROB_K5[idx] : PROB_K4[idx];
		float kdif = (k5dif < 1.0f) ? k5dif : k4dif;
		double floor = params[0];
		double log2len = Math.log(Math.max(length, 1)) * INV_LN2 - LOG2_REF;
		double z = params[1] + params[2] * log2len + params[3] * k3dif + params[4] * kdif;
		double sigmoid = 1.0 / (1.0 + Math.exp(z));
		return (float)(floor + (1.0 - floor) * sigmoid);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Thread-local nets     ----------------*/
	/*--------------------------------------------------------------*/

	private static CellNet getThreadNetV2(int idx){
		CellNet[] local=threadNetsV2.get();
		if(local==null){
			local=new CellNet[levelsV2.length];
			threadNetsV2.set(local);
		}
		if(local[idx]==null){
			final CellNet master=loaded.allLenNets[idx];
			if(master==null){return null;}
			local[idx]=master.copy(false);
		}
		return local[idx];
	}

	/** Thread-local copy of the single multi-output net. */
	private static CellNet getThreadMultiNet(){
		CellNet n=threadMultiNet.get();
		if(n==null){n=loaded.multiNet.copy(false); threadMultiNet.set(n);}
		return n;
	}

	/** Per-output calibration for the multi-output net: isotonic table when present, else 4-param. */
	private static float applyMultiCal(float raw, int idx){
		if(loaded.multiLutX!=null && loaded.multiLutX[idx]!=null){
			return lookup(raw, loaded.multiLutX[idx], loaded.multiLutY[idx]);
		}
		return calibrate(raw, loaded.multiCal[idx]);
	}

	/** Thread-local copy of the optional slow (49-input) multi-output net. */
	private static CellNet getThreadMultiNetSlow(){
		CellNet n=threadMultiNetSlow.get();
		if(n==null){n=loadedSlow.multiNet.copy(false); threadMultiNetSlow.set(n);}
		return n;
	}

	/** Per-output calibration for the slow net's heads: isotonic table when present, else 4-param. */
	private static float applyMultiCalSlow(float raw, int idx){
		if(loadedSlow.multiLutX!=null && loadedSlow.multiLutX[idx]!=null){
			return lookup(raw, loadedSlow.multiLutX[idx], loadedSlow.multiLutY[idx]);
		}
		return calibrate(raw, loadedSlow.multiCal[idx]);
	}

	private static CellNet getThreadNet(CellNet master, int idx, int length) {
		CellNet[][] local = threadNets.get();
		if(local==null){
			local = new CellNet[NUM_LEVELS][NUM_BINS+1];
			threadNets.set(local);
		}
		int slot;
		if(loaded.allLenNets[idx]!=null){
			slot=NUM_BINS;
		}else{
			slot=binIndex(length);
		}
		if(local[idx][slot]==null){
			local[idx][slot] = master.copy(false);
		}
		return local[idx][slot];
	}

	private static FloatList getThreadInput() {
		FloatList fl = threadInput.get();
		if(fl==null){
			fl = new FloatList(VECTOR_DIMS);
			threadInput.set(fl);
		}
		return fl;
	}

	private static final ThreadLocal<CellNet[]> threadNetsV2 = new ThreadLocal<>();
	private static final ThreadLocal<CellNet> threadMultiNet = new ThreadLocal<>();
	private static final ThreadLocal<CellNet> threadMultiNetSlow = new ThreadLocal<>();
	private static final ThreadLocal<CellNet[][]> threadNets = new ThreadLocal<>();
	private static final ThreadLocal<FloatList> threadInput = new ThreadLocal<>();

	/*--------------------------------------------------------------*/
	/*----------------            Lookups           ----------------*/
	/*--------------------------------------------------------------*/

	static int binIndex(int length) {
		int raw = 63 - Long.numberOfLeadingZeros(Math.max(1, (long)length / 2500));
		return Tools.mid(0, raw, NUM_BINS-1);
	}

	/** V1 level->index; the legacy 8-level model has no superkingdom slot. */
	private static int levelToIndexV1(int taxLevel) {
		switch (taxLevel) {
			case TaxTree.SPECIES: return 0;
			case TaxTree.GENUS:   return 1;
			case TaxTree.FAMILY:  return 2;
			case TaxTree.ORDER:   return 3;
			case TaxTree.CLASS:   return 4;
			case TaxTree.PHYLUM:  return 5;
			case TaxTree.KINGDOM: return 6;
			case TaxTree.DOMAIN:  return 7;
			default: return -1;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	private static SerialNNLoader.LoadedNets loadNets() {
		String path = Data.findPath("?confidence.bbnets.gz", true);
		if(path==null){return null;}
		return SerialNNLoader.load(path);
	}

	/** Optional 49-input slow bundle (ranking dim). Null when absent -> the fast bundle covers every mode.
	 *  warn=false: this file is OPTIONAL, so its absence must be silent (no scary "Cannot find" + stack
	 *  trace). A deployment with only confidence.bbnets.gz is a supported configuration, not an error. */
	private static SerialNNLoader.LoadedNets loadSlowNets() {
		String path = Data.findPath("?confidence49.bbnets.gz", false);
		if(path==null){return null;}
		return SerialNNLoader.load(path);
	}

	/**
	 * A bundle is V2 only if it supplies a complete set: 9 levels, every network
	 * present, every network taking the 48-dim vector. Anything short of that is
	 * treated as V1 (or as absent), so a partial or older file degrades to a
	 * working model instead of producing confident nonsense.
	 */
	private static boolean detectV2(SerialNNLoader.LoadedNets ln){
		if(ln==null || ln.levels<1 || ln.allLenNets==null){return false;}
		if(buildLevels(ln).length!=ln.levels){return false;} //every net must declare a known level
		for(int i=0; i<ln.levels; i++){
			final CellNet net=ln.allLenNets[i];
			if(net==null || net.numInputs()!=VECTOR_DIMS || net.numOutputs()!=1){return false;}
		}
		return true;
	}

	/**
	 * A single-net multi-output bundle: `#multioutput 1`, one 48-input net with `levels` output heads,
	 * every head declaring a known level. Backward-compatible: the per-level V2 path is untouched, and a
	 * bundle that is not multi-output makes this return false.
	 */
	private static boolean detectMultiV2(SerialNNLoader.LoadedNets ln){
		if(ln==null || !ln.multioutput || ln.multiNet==null || ln.levels<1){return false;}
		if(ln.multiNet.numInputs()!=VECTOR_DIMS || ln.multiNet.numOutputs()!=ln.levels){return false;}
		return buildLevelsFromLabels(ln.multiLabels).length==ln.levels;
	}

	/**
	 * The OPTIONAL slow bundle: a single multi-output net taking the 49-dim vector (the 48 features plus
	 * Neptune's ranking score as dim 49). Loaded from confidence49.bbnets.gz and used per-hit only in slow
	 * mode, where the ranking net has run. Requires exactly VECTOR_DIMS+1 inputs so it can never be confused
	 * with the 48-dim fast bundle.
	 */
	private static boolean detectMultiV2Slow(SerialNNLoader.LoadedNets ln){
		if(ln==null || !ln.multioutput || ln.multiNet==null || ln.levels<1){return false;}
		if(ln.multiNet.numInputs()!=VECTOR_DIMS+1 || ln.multiNet.numOutputs()!=ln.levels){return false;}
		return buildLevelsFromLabels(ln.multiLabels).length==ln.levels;
	}

	/**
	 * A slow net's outputs get LABELED by activeLevels(), which returns the PRIMARY (fast) model's level
	 * array whenever a fast model is loaded. So the slow bundle is usable in the hybrid only if its levels
	 * are element-wise identical to the primary's; otherwise a slow-scored hit would be paired with the
	 * wrong (or wrong-length) level list -- silent-wrong output, or an AIOOBE in Comparison.confidence. A
	 * mismatched pair is a packaging error, so disable the slow net and use the 48-dim net everywhere.
	 * (Slow-only, no fast model: the primary IS the slow list, so this is trivially satisfied.)
	 */
	private static boolean slowLevelsMatchPrimary(SerialNNLoader.LoadedNets slow){
		if(slow==null || slow.multiLabels==null){return false;}
		final int[] s=buildLevelsFromLabels(slow.multiLabels);
		if(s.length==0){return false;}
		final int[] primary=multiV2 ? levelsMulti : (v2 ? levelsV2 : s);
		return java.util.Arrays.equals(s, primary);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final boolean USE_NN = true;
	private static final int NUM_LEVELS = 8;
	private static final int NUM_BINS = 11;
	private static final double INV_LN2 = 1.0 / Math.log(2);
	private static final double LOG2_REF = Math.log(2500) * INV_LN2;

	/** Levels the loaded bundle actually supplies, in file order. Empty when no V2 model is active. */
	public static int[] activeLevels(){return multiV2 ? levelsMulti : (multiV2Slow ? levelsMultiSlow : levelsV2);}
	static final int MAIN_DIMS=26, SMALL_DIMS=11, VECTOR_DIMS=MAIN_DIMS+2*SMALL_DIMS; //48
	static final int TOP_EXAMINE=10;
	/**
	 * Maps a #label name to its TaxTree level. Deliberately NOT a fixed nine-element list:
	 * the shipped bundle omits superkingdom (no node in the taxonomy carries that rank, so
	 * its network merely duplicates kingdom), and a level absent from the file must be
	 * absent from inference too rather than silently shifting every index after it.
	 */
	private static int levelFromLabel(String label){
		if(label==null){return -1;}
		final String s=label.trim().toLowerCase();
		if(s.equals("species")){return TaxTree.SPECIES;}
		if(s.equals("genus")){return TaxTree.GENUS;}
		if(s.equals("family")){return TaxTree.FAMILY;}
		if(s.equals("order")){return TaxTree.ORDER;}
		if(s.equals("class")){return TaxTree.CLASS;}
		if(s.equals("phylum")){return TaxTree.PHYLUM;}
		if(s.equals("kingdom")){return TaxTree.KINGDOM;}
		if(s.equals("superkingdom")){return TaxTree.SUPERKINGDOM;}
		if(s.equals("domain")){return TaxTree.DOMAIN;}
		return -1;
	}

	/** Level list built from the bundle; null/empty means no usable V2 model. */
	private static int[] buildLevels(SerialNNLoader.LoadedNets ln){
		if(ln==null){return new int[0];}
		return buildLevelsFromLabels(ln.allLenLabels);
	}

	/** Maps an array of #label names to TaxTree level codes; empty if any is unlabeled/unknown. */
	private static int[] buildLevelsFromLabels(String[] labels){
		if(labels==null || labels.length<1){return new int[0];}
		final int[] out=new int[labels.length];
		for(int i=0; i<labels.length; i++){
			out[i]=levelFromLabel(labels[i]);
			if(out[i]<0){return new int[0];} //unlabeled or unknown -> refuse the whole bundle
		}
		return out;
	}

	//CLEVER [verified in-file]: graceful NN->sigmoid degradation. loaded=loadNets() returns null if confidence.bbnets.gz is absent (loadNets null on missing path); probCorrect then falls through to predictSigmoid (5-param PROB_K4/K5 tables) -- calibrated confidence even with no model file. Per-thread net copies (getThreadNet/getThreadNetV2) avoid CellNet contention; calibrate() clamps x to [1e-4,0.9999] so logit never hits +/-Inf.
	private static final SerialNNLoader.LoadedNets loaded = loadNets();
	private static final boolean v2 = detectV2(loaded);
	private static final boolean multiV2 = detectMultiV2(loaded);
	private static final int[] levelsV2 = v2 ? buildLevels(loaded) : new int[0];
	private static final int[] levelsMulti = multiV2 ? buildLevelsFromLabels(loaded.multiLabels) : new int[0];
	private static final SerialNNLoader.LoadedNets loadedSlow = loadSlowNets();
	private static final boolean multiV2Slow = detectMultiV2Slow(loadedSlow) && slowLevelsMatchPrimary(loadedSlow);
	private static final int[] levelsMultiSlow = multiV2Slow ? buildLevelsFromLabels(loadedSlow.multiLabels) : new int[0];

	// Sigmoid fallback parameters
	static final double[][] PROB_K4 = {
		{0.0113446, -3.5630508, 0.6403488, 51.9226674, 13.6112786},
		{0.0795275, -5.0429731, 0.4555063, 59.6982223, 10.8967730},
		{0.2127537, -4.4862228, 0.2292117, 48.4028932,  7.7315173},
		{0.3256019, -5.2138624, 0.1925119, 49.9547608,  7.6156285},
		{0.5084922, -5.4920396, 0.0667606, 63.4318036,  0.5094376},
		{0.6071749, -6.5750378, 0.0950648, 57.9167172,  6.3128603},
		{0.6801771, -7.1822199, 0.1554127, 43.4878516, 15.1889807},
		{0.8182326, -8.3517371, 0.1448564, 32.9901037, 21.5209495},
	};

	static final double[][] PROB_K5 = {
		{0.0058409, -2.5413033, 0.4705238, 70.8753575, -0.7640568},
		{0.0556486, -3.6911274, 0.3030928, 71.3175058, -2.1744591},
		{0.1795990, -3.3621772, 0.1272263, 57.1529654, -2.5656541},
		{0.2788433, -4.2430087, 0.1167844, 57.7251810, -2.0374220},
		{0.4376369, -4.2088439,-0.0017810, 63.6425910, -5.5906412},
		{0.5352429, -5.0486971, 0.0025206, 64.2620784, -4.7840435},
		{0.6134413, -5.4231686, 0.0501707, 55.7766015, -2.2099322},
		{0.7916569, -7.1943760, 0.0995638, 49.9422788,  2.2346173},
	};
}
