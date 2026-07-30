package clade;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;

/**
 * Evaluates the calibration of QuickClade per-hit confidence CONDITIONALLY on the two drivers that
 * should govern it -- query length and similarity-to-reference (KID). Consumes the (qlen, kid, conf,
 * label) pair files emitted per (bundle, clade, level) by {@link ConfidenceAudit}, and reports:
 *
 *   - Conditional calibration error: mean |mean_conf - accuracy| over occupied length x similarity
 *     slices (>= minn samples), macro-averaged across clades so every clade counts equally.
 *   - Brier score (a proper scoring rule) per clade and overall -- rewards calibration AND sharpness,
 *     so a useless base-rate predictor cannot win the way it can on aggregate ECE alone.
 *   - Monotonicity at the fine levels: mean confidence must RISE with length and with similarity.
 *
 * A well-calibrated net gives longer, closer shreds higher fine-level confidence BECAUSE they are
 * genuinely more often correct, so conditional calibration operationalizes "the results feel right"
 * with no tuned weights. Compare two bundles by whichever holds calibration across ALL slices and
 * shows the correct monotonic trends.
 *
 * Usage: confidenceeval.sh dir=<pairdir> bundles=NEW,OLD clades=archaea,bacteria,... [levels=...] [minn=40]
 * Pair files are named &lt;bundle&gt;_&lt;clade&gt;_pairs_&lt;level&gt;.tsv (cols: qlen kid conf label).
 *
 * @author Noire
 */
public class ConfidenceEval {

	public static void main(String[] args){
		ConfidenceEval x=new ConfidenceEval(args);
		x.process();
	}

	public ConfidenceEval(String[] args){
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("dir") || a.equals("in") || a.equals("path")){dir=b;}
			else if(a.equals("bundles")){bundles=b.split(",");}
			else if(a.equals("clades")){clades=b.split(",");}
			else if(a.equals("levels")){levels=b.split(",");}
			else if(a.equals("minn")){minn=Integer.parseInt(b);}
			else if(a.equals("tol")){tol=Float.parseFloat(b);}
			else if(a.equals("vtol")){vtol=Float.parseFloat(b);}
			else{throw new RuntimeException("Unknown parameter: "+arg);}
		}
		if(dir==null){throw new RuntimeException("dir= (directory of pair files) is required");}
		if(!dir.endsWith("/")){dir=dir+"/";}
		if(clades==null){throw new RuntimeException("clades= is required");}
	}

	/*--------------------------------------------------------------*/

	/** Ordered similarity buckets. kid<0 = no sketch comparison; kid==0 = sketch present but ZERO
	 *  identity, a large and meaningful class (genuinely novel / low-similarity queries) that earns
	 *  its own bucket rather than being merged upward or dropped through a bin-edge gap; positive kid
	 *  falls in ascending ranges. Monotonicity skips only "noSketch"; "k=0" is the bottom rung. */
	static final String[] KID_NAME={"noSketch", "k=0", "0-.6", ".6-.8", ".8-.9", ".9-.97", ">.97"};

	static int kidBucket(float k){
		if(k<0){return 0;}
		if(k==0){return 1;}
		if(k<=0.6f){return 2;}
		if(k<=0.8f){return 3;}
		if(k<=0.9f){return 4;}
		if(k<=0.97f){return 5;}
		return 6;
	}

	/** One (length, kidBucket) slice: running sums of confidence and label. */
	static class Cell{
		double sumConf=0, sumLabel=0; long n=0;
		void add(float conf, int label){sumConf+=conf; sumLabel+=label; n++;}
		double meanConf(){return sumConf/Math.max(n,1);}
		double accuracy(){return sumLabel/Math.max(n,1);}
		double absGap(){return Math.abs(meanConf()-accuracy());}
	}

	/** Reads a pair file; returns rows as float[]{qlen,kid,conf,label}. Missing file -> empty. */
	ArrayList<float[]> readPairs(String bundle, String clade, String level){
		final ArrayList<float[]> rows=new ArrayList<float[]>();
		final String path=dir+bundle+"_"+clade+"_pairs_"+level+".tsv";
		final FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, true, false);
		if(ff==null || !ff.exists()){return rows;}
		final ByteFile bf=ByteFile.makeByteFile(ff);
		byte[] line;
		while((line=bf.nextLine())!=null){
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] f=new String(line).split("\t");
			if(f.length<4){continue;}
			try{
				rows.add(new float[]{Float.parseFloat(f[0]), Float.parseFloat(f[1]),
						Float.parseFloat(f[2]), Float.parseFloat(f[3])});
			}catch(NumberFormatException e){/* skip malformed */}
		}
		bf.close();
		return rows;
	}

	/*--------------------------------------------------------------*/

	void process(){
		final PrintStream out=System.out;
		out.println("CONDITIONAL CALIBRATION EVAL -- confidence honesty within length x similarity slices.");
		out.println("Scalar = mean |mean_conf - accuracy| over slices (>= "+minn+" samples), macro over clades.");
		out.println("Brier = mean (conf-label)^2 (proper score). Both lower = better.\n");

		final double[] macroCCE=new double[bundles.length];
		final double[] macroBrier=new double[bundles.length];

		for(int bi=0; bi<bundles.length; bi++){
			final String bundle=bundles[bi];
			out.println("### "+bundle);
			double clSum=0, brSum=0; int clN=0;
			for(String clade : clades){
				final ArrayList<Double> cellErrs=new ArrayList<Double>();
				double brierSum=0; long brierN=0;
				for(String level : levels){
					final LinkedHashMap<Long,Cell> cells=new LinkedHashMap<Long,Cell>();
					for(float[] r : readPairs(bundle, clade, level)){
						final float qlen=r[0], kid=r[1], conf=r[2], label=r[3];
						if(conf<0){continue;}
						final long key=(((long)qlen)<<8)|kidBucket(kid);
						Cell c=cells.get(key);
						if(c==null){c=new Cell(); cells.put(key, c);}
						c.add(conf, (int)label);
						brierSum+=(conf-label)*(conf-label); brierN++;
					}
					for(Cell c : cells.values()){
						if(c.n>=minn){cellErrs.add(c.absGap());}
					}
				}
				final double cce=mean(cellErrs);
				final double brier=brierSum/Math.max(brierN,1);
				out.println(String.format("    %-22s  condCalErr=%.4f  brier=%.4f  (%d slices)",
						clade, cce, brier, cellErrs.size()));
				if(!Double.isNaN(cce)){clSum+=cce; brSum+=brier; clN++;}
			}
			macroCCE[bi]=clN>0 ? clSum/clN : Double.NaN;
			macroBrier[bi]=clN>0 ? brSum/clN : Double.NaN;
			out.println(String.format("  >> %s MACRO  condCalErr=%.4f  brier=%.4f\n",
					bundle, macroCCE[bi], macroBrier[bi]));
		}

		if(bundles.length==2){
			final double d=macroCCE[0]-macroCCE[1];
			final String v=(d<-vtol ? bundles[0]+" better" : (d>vtol ? bundles[1]+" better" : "~equal"));
			out.println(String.format(">>> %s condCalErr=%.4f vs %s=%.4f  (diff %+.4f -> %s)\n",
					bundles[0], macroCCE[0], bundles[1], macroCCE[1], d, v));
		}

		out.println("MONOTONICITY (fine levels): mean confidence must RISE with length and similarity.");
		for(String bundle : bundles){
			out.println("### "+bundle);
			for(String level : new String[]{"species", "genus"}){
				for(String clade : clades){monotonic(out, bundle, clade, level);}
			}
		}
	}

	/** Prints mean confidence by ascending length and by ascending KID bucket, flagging drops. */
	void monotonic(PrintStream out, String bundle, String clade, String level){
		final LinkedHashMap<Integer,Cell> byLen=new LinkedHashMap<Integer,Cell>();
		final Cell[] byKid=new Cell[KID_NAME.length];
		for(int i=0; i<byKid.length; i++){byKid[i]=new Cell();}
		for(float[] r : readPairs(bundle, clade, level)){
			final float conf=r[2];
			if(conf<0){continue;}
			final int len=(int)r[0];
			Cell c=byLen.get(len);
			if(c==null){c=new Cell(); byLen.put(len, c);}
			c.add(conf, (int)r[3]);
			byKid[kidBucket(r[1])].add(conf, (int)r[3]);
		}
		final Integer[] lens=byLen.keySet().toArray(new Integer[0]);
		Arrays.sort(lens);
		final StringBuilder ls=new StringBuilder();
		double prev=-1; int lenViol=0;
		for(Integer L : lens){
			final Cell c=byLen.get(L);
			if(c.n<minn){continue;}
			final double m=c.meanConf();
			if(prev>=0 && m<prev-tol){lenViol++;}
			prev=m;
			ls.append(L).append(':').append(String.format("%.2f", m)).append(' ');
		}
		final StringBuilder ks=new StringBuilder();
		prev=-1; int kidViol=0;
		for(int i=1; i<byKid.length; i++){ //skip noSketch bucket 0
			final Cell c=byKid[i];
			if(c.n<minn){continue;}
			final double m=c.meanConf();
			if(prev>=0 && m<prev-tol){kidViol++;}
			prev=m;
			ks.append(KID_NAME[i]).append(':').append(String.format("%.2f", m)).append(' ');
		}
		String flag="";
		if(lenViol>0){flag+=" LEN-VIOL("+lenViol+")";}
		if(kidViol>0){flag+=" KID-VIOL("+kidViol+")";}
		out.println(String.format("  %-22s %-7s  len[%s] kid[%s]%s", clade, level, ls.toString().trim(), ks.toString().trim(), flag));
	}

	static double mean(ArrayList<Double> a){
		if(a.isEmpty()){return Double.NaN;}
		double s=0; for(Double d : a){s+=d;} return s/a.size();
	}

	/*--------------------------------------------------------------*/

	String dir=null;
	String[] bundles={"NEW", "OLD"};
	String[] clades=null;
	String[] levels={"species","genus","family","order","class","phylum","kingdom","domain"};
	int minn=40;
	float tol=0.01f;     //monotonicity tolerance (a drop bigger than this flags a violation)
	float vtol=0.0005f;  //NEW-vs-OLD verdict tolerance (distinct from monotonicity)
}
