package clade;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import structures.FloatList;

/**
 * Head-to-head: does ONE continuous depth-regression net (RegressionTrainer output) calibrate as well
 * as the shipped 8-net boolean confidence bundle, on the SAME data and the SAME metric? This is the
 * "1 continuous net vs 9 boolean nets" question, put onto a common axis so the two are comparable.
 *
 * Both models are scored by per-level CONDITIONAL calibration -- the same idea as {@link ConfidenceEval}:
 * mean |mean_confidence - accuracy| over occupied (length x KID-similarity) slices, so a model only wins
 * by being honest WITHIN each length/similarity regime, not by exploiting the base rate. Because val_CS
 * carries no clade tag, we recover the 7-way domain one-hot from the vector itself (dims 17-23) and
 * macro-average over domains (rare domains count equally), and also report the pooled number.
 *
 * Boolean side: the bundle's per-level net feedForwards the raw 48-dim vector, then applyCal (isotonic
 * LUT if present, else the 4-param K*sigmoid) -- replicated verbatim from CladeConfidence V2 inference.
 * Continuous side: the single net feedForwards the same vector to a predicted DEPTH in [0,1]; that depth
 * is mapped to P(correct at level L) by a per-level ISOTONIC (PAV) calibrator FIT ON A DISJOINT SET
 * (isofit=, default a stride-sample of reg_train) and evaluated on val -- so the continuous side gets no
 * in-sample advantage, matching how the boolean bundle's own calibration was fit on training data.
 *
 * "Correct at level L" is derived from the continuous depth label: encodeLevel(level)=(11-level)/9, so
 * true_depth >= threshold(L) means the hit's LCA is at least as deep as L. Levels species..domain (8),
 * superkingdom skipped to match the shipped bundle.
 *
 * NOTE: the pooled/macro-domain numbers here are NOT the 0.0129 macro-over-CLADE live-audit figure; this
 * is an apples-to-apples model-vs-model comparison on val_CS. Lower = better for both.
 *
 * Usage: confidenceheadtohead.sh bundle=<confidence.bbnets.gz> cont=<net.bbnet> val=<val_CS.tsv>
 *        isofit=<reg_train.tsv> [maxrows=1000000] [isomax=400000] [isostride=30] [minn=40]
 *
 * @author Noire
 */
public class ConfidenceHeadToHead {

	// The 8 boolean levels, in a fixed order, with TaxTree level codes and depth thresholds.
	static final String[] LNAME={"species","genus","family","order","class","phylum","kingdom","domain"};
	static final int[] LCODE={2,3,4,5,6,7,8,10};                 //TaxTree codes (superkingdom=9 skipped)
	static final int NL=LNAME.length;
	static final String[] DNAME={"bacteria","archaea","virus","animal","plant","fungi","otherEuk","unclass"};
	static final int[] SIZES={1000,2500,5000,10000,25000,50000,100000,250000,500000,1000000};
	static final float EPS=1e-4f;

	public static void main(String[] args){
		String bundlePath=null, contPath=null, valPath=null, isoPath=null;
		int maxrows=1000000, isomax=400000, isostride=30, minn=40;
		for(String a : args){
			final String[] s=a.split("=", 2);
			final String k=s[0].toLowerCase(), v=(s.length>1 ? s[1] : null);
			if(k.equals("bundle")){bundlePath=v;}
			else if(k.equals("cont") || k.equals("net")){contPath=v;}
			else if(k.equals("val") || k.equals("in")){valPath=v;}
			else if(k.equals("isofit") || k.equals("iso")){isoPath=v;}
			else if(k.equals("maxrows")){maxrows=Integer.parseInt(v);}
			else if(k.equals("isomax")){isomax=Integer.parseInt(v);}
			else if(k.equals("isostride")){isostride=Integer.parseInt(v);}
			else if(k.equals("minn")){minn=Integer.parseInt(v);}
			else{throw new RuntimeException("Unknown parameter: "+a);}
		}
		if(bundlePath==null||contPath==null||valPath==null||isoPath==null){
			throw new RuntimeException("Required: bundle= cont= val= isofit=");
		}

		//--- load the boolean bundle and align its per-file levels to LNAME order ---
		final SerialNNLoader.LoadedNets ln=SerialNNLoader.load(bundlePath);
		if(ln==null){throw new RuntimeException("bundle load failed: "+bundlePath);}
		final CellNet[] boolNet=new CellNet[NL];
		final float[][] boolCal=new float[NL][], boolLutX=new float[NL][], boolLutY=new float[NL][];
		for(int i=0; i<ln.levels; i++){
			final int my=nameIndex(ln.allLenLabels==null?null:ln.allLenLabels[i]);
			if(my<0){continue;}                                  //e.g. superkingdom -> not in the 8
			boolNet[my]=ln.allLenNets[i];
			boolCal[my]=ln.allLenCal!=null?ln.allLenCal[i]:null;
			boolLutX[my]=ln.allLenLutX!=null?ln.allLenLutX[i]:null;
			boolLutY[my]=ln.allLenLutY!=null?ln.allLenLutY[i]:null;
		}
		for(int L=0; L<NL; L++){
			if(boolNet[L]==null){throw new RuntimeException("bundle missing level "+LNAME[L]);}
		}
		final CellNet cont=CellNetParser.load(contPath);
		if(cont==null){throw new RuntimeException("continuous net load failed: "+contPath);}
		final boolean multi=(cont.numOutputs()==NL);
		if(cont.numOutputs()!=1 && !multi){
			throw new RuntimeException("cont net has "+cont.numOutputs()+" outputs; expected 1 (scalar depth) or "+NL+" (per-level heads)");
		}
		System.err.println("cont mode: "+(multi ? ("MULTI ("+NL+" per-level output heads)") : "SCALAR (one depth + per-level isotonic)"));
		final int indims=cont.numInputs();   //48, or 49 with the ranking dim; bundle stays 48 (reads first 48 cols)

		//--- FIT per-level isotonic for the continuous net on the DISJOINT isofit sample ---
		final FloatList fl=new FloatList(48);
		final ArrayList<float[]> fitScore=new ArrayList<float[]>();   //per-level score (float[NL]) per kept row
		final ArrayList<byte[]> fitLab=new ArrayList<byte[]>();       //per-level label (byte[NL]) per kept row
		{
			final ByteFile bf=ByteFile.makeByteFile(isoPath, true);
			byte[] line; long seen=0;
			while((line=bf.nextLine())!=null && fitScore.size()<isomax){
				if(line.length==0 || line[0]=='#'){continue;}
				if((seen++ % isostride)!=0){continue;}
				final float[] row=parse(line, indims+1);
				if(row==null){continue;}
				final float[] sc=scores(cont, fl, row, multi);       //scalar->replicated; multi->8 heads
				final byte[] lab=new byte[NL];
				for(int L=0; L<NL; L++){lab[L]=(byte)(row[indims]>=thresh(L)-EPS?1:0);}
				fitScore.add(sc); fitLab.add(lab);
			}
			bf.close();
		}
		final int nf=fitScore.size();
		System.err.println("isofit rows="+nf);
		final float[][] isoX=new float[NL][], isoY=new float[NL][];
		{
			final float[] sc=new float[nf];
			final byte[] lb=new byte[nf];
			for(int L=0; L<NL; L++){
				for(int i=0; i<nf; i++){sc[i]=fitScore.get(i)[L]; lb[i]=fitLab.get(i)[L];}
				final float[][] knots=pav(sc, lb);
				isoX[L]=knots[0]; isoY[L]=knots[1];
			}
		}

		//--- EVAL: stream val, accumulate conditional-calibration cells for both models ---
		//key = domain*10^7 + level*10^6 + sizeIdx*10 + kidBucket   (domain 8, level 8, size 10, kid 7)
		final LinkedHashMap<Long,Cell> boolCells=new LinkedHashMap<Long,Cell>();
		final LinkedHashMap<Long,Cell> contCells=new LinkedHashMap<Long,Cell>();
		double boolBrier=0, contBrier=0; long brierN=0;
		{
			final ByteFile bf=ByteFile.makeByteFile(valPath, true);
			byte[] line; int n=0;
			while((line=bf.nextLine())!=null && n<maxrows){
				if(line.length==0 || line[0]=='#'){continue;}
				final float[] row=parse(line, indims+1);
				if(row==null){continue;}
				n++;
				final int dom=domainOf(row), sz=sizeIdx(row[0]), kb=kidBucket(row[10]);
				final float[] cs=scores(cont, fl, row, multi);   //per-level score (scalar->replicated; multi->8 heads)
				for(int L=0; L<NL; L++){
					final int label=(row[indims]>=thresh(L)-EPS)?1:0;
					final float bp=applyCal(feed(boolNet[L], fl, row), boolCal[L], boolLutX[L], boolLutY[L]);
					final float cp=lookup(cs[L], isoX[L], isoY[L]);
					final long key=((long)dom)*10000000L+((long)L)*1000000L+((long)sz)*10L+kb;
					add(boolCells, key, bp, label);
					add(contCells, key, cp, label);
					boolBrier+=(bp-label)*(bp-label); contBrier+=(cp-label)*(cp-label); brierN++;
				}
			}
			bf.close();
			System.err.println("val eval rows="+n);
		}

		//--- report ---
		final PrintStream o=System.out;
		o.println("HEAD-TO-HEAD: continuous depth net vs boolean 8-net bundle (val_CS, conditional calibration).");
		o.println("condCalErr = mean |mean_conf - accuracy| over (length x KID) slices >= "+minn+" samples. Lower=better.");
		o.println("NOT the 0.0129 macro-CLADE live figure; this is model-vs-model on val_CS.\n");

		o.println(String.format("%-10s %12s %12s", "level", "BOOL cCE", "CONT cCE"));
		for(int L=0; L<NL; L++){
			o.println(String.format("%-10s %12.5f %12.5f", LNAME[L],
					levelCCE(boolCells, L, minn), levelCCE(contCells, L, minn)));
		}
		o.println();
		final double bMacro=macroDomain(boolCells, minn), cMacro=macroDomain(contCells, minn);
		final double bPool=pooled(boolCells, minn), cPool=pooled(contCells, minn);
		o.println(String.format("%-18s BOOL=%.5f   CONT=%.5f   (diff %+.5f)", "MACRO-DOMAIN cCE",
				bMacro, cMacro, cMacro-bMacro));
		o.println(String.format("%-18s BOOL=%.5f   CONT=%.5f   (diff %+.5f)", "POOLED cCE",
				bPool, cPool, cPool-bPool));
		o.println(String.format("%-18s BOOL=%.5f   CONT=%.5f", "BRIER",
				boolBrier/Math.max(brierN,1), contBrier/Math.max(brierN,1)));
		o.println();
		final String v=(cMacro<bMacro-0.0005 ? "CONTINUOUS calibrates better"
				: (cMacro>bMacro+0.0005 ? "BOOLEAN calibrates better" : "~EQUAL calibration"));
		o.println(">>> "+v+" (by macro-domain cCE; +diff = continuous worse).");
	}

	/*----------------------------------------------------------------*/

	/** One (model, slice) accumulator. */
	static class Cell{ double sumConf=0, sumLab=0; long n=0;
		void add(float c, int l){sumConf+=c; sumLab+=l; n++;}
		double gap(){return Math.abs(sumConf/n - sumLab/n);}
	}
	static void add(LinkedHashMap<Long,Cell> m, long key, float conf, int label){
		Cell c=m.get(key); if(c==null){c=new Cell(); m.put(key, c);} c.add(conf, label);
	}

	/** Mean gap over this level's occupied slices (pooled over domain). */
	static double levelCCE(LinkedHashMap<Long,Cell> m, int level, int minn){
		double s=0; int k=0;
		//re-slice ignoring domain: key by (level,size,kid)
		final LinkedHashMap<Long,Cell> agg=new LinkedHashMap<Long,Cell>();
		for(java.util.Map.Entry<Long,Cell> e : m.entrySet()){
			final long key=e.getKey();
			if((key/1000000L)%10!=level){continue;}
			final long sub=key%1000000L;               //size*10+kid, drop domain+level
			Cell c=agg.get(sub); if(c==null){c=new Cell(); agg.put(sub, c);}
			c.sumConf+=e.getValue().sumConf; c.sumLab+=e.getValue().sumLab; c.n+=e.getValue().n;
		}
		for(Cell c : agg.values()){if(c.n>=minn){s+=c.gap(); k++;}}
		return k>0?s/k:Double.NaN;
	}

	/** Pooled cCE: all (level,size,kid) slices (domain collapsed), one mean. */
	static double pooled(LinkedHashMap<Long,Cell> m, int minn){
		final LinkedHashMap<Long,Cell> agg=new LinkedHashMap<Long,Cell>();
		for(java.util.Map.Entry<Long,Cell> e : m.entrySet()){
			final long sub=e.getKey()%10000000L;       //drop domain, keep level+size+kid
			Cell c=agg.get(sub); if(c==null){c=new Cell(); agg.put(sub, c);}
			c.sumConf+=e.getValue().sumConf; c.sumLab+=e.getValue().sumLab; c.n+=e.getValue().n;
		}
		double s=0; int k=0;
		for(Cell c : agg.values()){if(c.n>=minn){s+=c.gap(); k++;}}
		return k>0?s/k:Double.NaN;
	}

	/** Macro over domains: per-domain mean gap (over its level x slice cells), averaged across domains. */
	static double macroDomain(LinkedHashMap<Long,Cell> m, int minn){
		final double[] sum=new double[DNAME.length]; final int[] cnt=new int[DNAME.length];
		for(java.util.Map.Entry<Long,Cell> e : m.entrySet()){
			if(e.getValue().n<minn){continue;}
			final int dom=(int)(e.getKey()/10000000L);
			if(dom<0||dom>=DNAME.length){continue;}
			sum[dom]+=e.getValue().gap(); cnt[dom]++;
		}
		double s=0; int k=0;
		for(int d=0; d<DNAME.length; d++){if(cnt[d]>0){s+=sum[d]/cnt[d]; k++;}}
		return k>0?s/k:Double.NaN;
	}

	/*----------------------------------------------------------------*/

	static int nameIndex(String label){
		if(label==null){return -1;}
		final String s=label.trim().toLowerCase();
		for(int i=0; i<LNAME.length; i++){if(LNAME[i].equals(s)){return i;}}
		return -1;
	}
	/** encodeLevel threshold for level index L: (11-code)/9. */
	static float thresh(int L){return (11-LCODE[L])/9f;}

	/** Domain from the 7-way one-hot at cols 16..22 (dims 17-23); 7=unclass when all zero. */
	static int domainOf(float[] r){
		for(int i=0; i<7; i++){if(r[16+i]>0.5f){return i;}}
		return 7;
	}
	/** Snap log2len (col 0) to the nearest of the 10 shred sizes; returns its index 0..9. */
	static int sizeIdx(float log2len){
		int best=0; double bd=1e18;
		for(int i=0; i<SIZES.length; i++){
			final double d=Math.abs(Math.log(SIZES[i])/Math.log(2)-log2len);
			if(d<bd){bd=d; best=i;}
		}
		return best;
	}
	/** ConfidenceEval's KID buckets: <0 noSketch(0), ==0 k=0(1), then ascending ranges. col 10 = kid. */
	static int kidBucket(float k){
		if(k<0){return 0;} if(k==0){return 1;}
		if(k<=0.6f){return 2;} if(k<=0.8f){return 3;} if(k<=0.9f){return 4;} if(k<=0.97f){return 5;}
		return 6;
	}

	/** Parse a data line into float[49] (48 features + depth label); null if malformed. */
	static float[] parse(byte[] line, int ncol){
		final String[] f=new String(line).split("\t");
		if(f.length!=ncol){return null;}
		final float[] r=new float[ncol];
		try{for(int i=0; i<ncol; i++){r[i]=Float.parseFloat(f[i]);}}catch(Exception e){return null;}
		return r;
	}
	/** Feed a net's own numInputs raw features (from the front of the row); returns its scalar output. */
	static float feed(CellNet net, FloatList fl, float[] row){
		final int ni=net.numInputs();
		fl.size=0; for(int i=0; i<ni; i++){fl.add(row[i]);}
		net.applyInput(fl); return net.feedForward();
	}
	/** Per-level continuous score. multi=8 output heads -> one score per level; scalar=one depth
	 *  replicated to every level. Reads the net's own numInputs features (48 bundle-compat, or 49 w/ ranking). */
	static float[] scores(CellNet net, FloatList fl, float[] row, boolean multi){
		final int ni=net.numInputs();
		fl.size=0; for(int i=0; i<ni; i++){fl.add(row[i]);}
		net.applyInput(fl); final float f=net.feedForward();
		if(multi){return net.getOutput();}                 //getOutput() allocates fresh (length NL)
		final float[] s=new float[NL]; java.util.Arrays.fill(s, f); return s;
	}

	/*------- calibration replicated verbatim from CladeConfidence -------*/
	static float applyCal(float raw, float[] cal, float[] lx, float[] ly){
		if(lx!=null){return lookup(raw, lx, ly);}
		return calibrate(raw, cal);
	}
	static float lookup(float x, float[] lx, float[] ly){
		if(lx==null || lx.length==0){return x;}
		final int last=lx.length-1;
		if(x<=lx[0]){return ly[0];}
		if(x>=lx[last]){return ly[last];}
		int lo=0, hi=last;
		while(lo+1<hi){final int mid=(lo+hi)>>>1; if(lx[mid]<=x){lo=mid;}else{hi=mid;}}
		final float dx=lx[hi]-lx[lo];
		final float f=(dx<=0?0:(x-lx[lo])/dx);
		return ly[lo]+f*(ly[hi]-ly[lo]);
	}
	static float calibrate(float x, float[] p){
		if(p==null){return x;}
		x=Math.max(0.0001f, Math.min(0.9999f, x));
		final float K=p[0], a=p[1], b=p[2], c=p[3];
		final double lx=Math.log(x/(1.0-x));
		final double s=1.0/(1.0+Math.exp(-(a*lx+b)));
		return (float)Math.min(1.0, K*Math.pow(s, c));
	}

	/*------- pool-adjacent-violators isotonic fit: score->P, monotone non-decreasing -------*/
	/** Returns {knotX[], knotY[]} ascending in x; y is the pooled label rate. */
	static float[][] pav(float[] score, byte[] label){
		final int n=score.length;
		final Integer[] ord=new Integer[n];
		for(int i=0; i<n; i++){ord[i]=i;}
		java.util.Arrays.sort(ord, (a,b)->Float.compare(score[a], score[b]));
		final double[] val=new double[n], w=new double[n], ss=new double[n]; int nb=0;
		for(int k=0; k<n; k++){
			final int i=ord[k];
			val[nb]=label[i]; w[nb]=1; ss[nb]=score[i]; nb++;
			while(nb>=2 && val[nb-1]<val[nb-2]-1e-12){
				final double nw=w[nb-1]+w[nb-2];
				val[nb-2]=(val[nb-1]*w[nb-1]+val[nb-2]*w[nb-2])/nw;
				ss[nb-2]=ss[nb-1]+ss[nb-2]; w[nb-2]=nw; nb--;
			}
		}
		final float[] x=new float[nb], y=new float[nb];
		for(int b=0; b<nb; b++){x[b]=(float)(ss[b]/w[b]); y[b]=(float)val[b];}
		return new float[][]{x, y};
	}
}
