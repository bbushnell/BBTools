package ml;

import java.util.ArrayList;

import fileIO.ByteFile;

/**
 * Evaluate ranking nets against compositeScore on benchmark data.
 * Reports top-1 selection accuracy, position-weighted ranking score (WDCG),
 * and within-query concordance, for each net and the composite baseline.
 *
 * Usage: java ml.RankingEval vectors=<ranking_vectors.tsv> nets=<dir_of_bbnets>
 *
 * The vector file is RankingVectorizer output (#dims 33 1).
 * Rows are grouped by query (maxemit=10 consecutive rows per query).
 * Dim 29 (0-indexed 28) = 0.1 * min(nHits, 10); a change in this value OR
 * reaching 10 consecutive rows marks a query boundary.
 * Dim 27 (0-indexed 26) = normalized composite; the net re-ranks by its own score.
 * The label (column 34, 0-indexed 33) = LCA depth truth (1.0=exact, 0=life).
 *
 * @author Neptune
 * @date July 2026
 */
public class RankingEval{

	public static void main(String[] args) throws Exception{
		String vecFile=null, netDir=null;
		for(String arg : args){
			final String[] kv=arg.split("=", 2);
			final String k=kv[0].toLowerCase();
			final String v=(kv.length>1) ? kv[1] : "";
			if(k.equals("vectors") || k.equals("in")){vecFile=v;}
			else if(k.equals("nets") || k.equals("netdir")){netDir=v;}
		}
		if(vecFile==null || netDir==null){
			throw new IllegalArgumentException("required: vectors=<file> nets=<dir>");
		}

		// Load vectors and labels
		final ArrayList<float[]> inputs=new ArrayList<>();
		final ArrayList<Float> labels=new ArrayList<>();
		final ArrayList<Float> composites=new ArrayList<>();
		final ArrayList<Integer> qstarts=new ArrayList<>();

		int NI=33; // default, overridden by #dims header
		final ByteFile bf=ByteFile.makeByteFile(vecFile, true);
		byte[] line;
		int rowInQuery=0;
		float prevDim29=-1;
		while((line=bf.nextLine())!=null){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String h=new String(line);
				if(h.startsWith("#dims")){
					final String[] hf=h.split("\t");
					if(hf.length>=2){NI=Integer.parseInt(hf[1].trim());}
				}
				continue;
			}
			final String s=new String(line);
			final String[] f=s.split("\t");
			if(f.length<NI+1){continue;}
			final float[] x=new float[NI];
			for(int i=0; i<NI; i++){x[i]=Float.parseFloat(f[i]);}
			final float label=Float.parseFloat(f[NI]);
			final float dim29=x[28]; // 0.1 * nHits (same position in both 33 and 35 dim vectors)
			final float comp=x[26]; // normalized composite

			// Detect query boundary: dim29 changes or we hit 10 rows
			if(inputs.isEmpty() || rowInQuery>=10 || (dim29!=prevDim29 && !inputs.isEmpty())){
				qstarts.add(inputs.size());
				rowInQuery=0;
			}
			inputs.add(x);
			labels.add(label);
			composites.add(comp);
			prevDim29=dim29;
			rowInQuery++;
		}
		bf.close();

		final int n=inputs.size();
		final int nq=qstarts.size();
		System.err.println("Loaded "+n+" hits, "+nq+" queries from "+vecFile);

		// Load all .bbnet files from netDir
		final java.io.File dir=new java.io.File(netDir);
		final java.io.File[] netFiles=dir.listFiles((d,name)->name.endsWith(".bbnet"));
		if(netFiles==null || netFiles.length==0){
			throw new IllegalArgumentException("no .bbnet files in "+netDir);
		}
		java.util.Arrays.sort(netFiles);

		// Print net sizes
		System.err.println("\n=== NET SIZES ===");
		for(java.io.File nf : netFiles){
			System.err.println(String.format("%-24s %8d bytes  (%,.1f KB)",
				nf.getName().replace(".bbnet",""), nf.length(), nf.length()/1024.0));
		}

		// Run inference for all nets
		final ArrayList<String> netNames=new ArrayList<>();
		final ArrayList<double[]> netPreds=new ArrayList<>();
		for(java.io.File nf : netFiles){
			final CellNet net=CellNetParser.load(nf.getAbsolutePath());
			final double[] preds=new double[n];
			for(int i=0; i<n; i++){
				net.applyInput(inputs.get(i));
				preds[i]=net.feedForward();
			}
			netNames.add(nf.getName().replace(".bbnet",""));
			netPreds.add(preds);
		}

		// === TOP-1 SELECTION ===
		System.err.println("\n=== TOP-1 SELECTION (which hit does each ranker pick?) ===");
		System.err.println(String.format("%-24s %10s %10s %12s", "ranker", "exact%", ">=genus%", "meanDepth"));
		evalSelection("composite(incumbent)", null, composites, labels, qstarts, n);
		for(int ni=0; ni<netNames.size(); ni++){
			evalSelection(netNames.get(ni), netPreds.get(ni), null, labels, qstarts, n);
		}
		evalSelection("ORACLE", null, null, labels, qstarts, n);

		// === POSITION-WEIGHTED RANKING SCORE (WDCG) ===
		// Weights: position 1 gets 10, position i gets 10/i.
		// Re-ranks hits by each ranker's score, then sums truth[rank_i] * 10/i.
		// Normalized by the oracle's score so 1.0 = perfect ranking.
		System.err.println("\n=== POSITION-WEIGHTED SCORE (10/i weights, normalized by oracle) ===");
		System.err.println(String.format("%-24s %10s %12s %12s", "ranker", "nWDCG", "rawWDCG", "oracleWDCG"));
		evalWDCG("composite(incumbent)", null, composites, labels, qstarts, n);
		for(int ni=0; ni<netNames.size(); ni++){
			evalWDCG(netNames.get(ni), netPreds.get(ni), null, labels, qstarts, n);
		}

		// === WITHIN-QUERY CONCORDANCE ===
		System.err.println("\n=== WITHIN-QUERY CONCORDANCE vs LCA depth (0.5=random) ===");
		System.err.println(String.format("%-24s %10s %10s %8s %10s", "ranker", "withTies", "noTies", "tie%", "pairs"));
		evalConcordance("composite", null, composites, labels, qstarts, n);
		for(int ni=0; ni<netNames.size(); ni++){
			evalConcordance(netNames.get(ni), netPreds.get(ni), null, labels, qstarts, n);
		}

		// === COMPOSITE GRADE: 10*top1_depth + concordance ===
		System.err.println("\n=== COMPOSITE GRADE (10*meanTop1Depth + concordance, higher=better) ===");
		System.err.println(String.format("%-24s %8s %10s", "ranker", "KB", "grade"));
		{
			double[] top1d=getTop1Depths(null, composites, labels, qstarts, n);
			double conc=getConcordance(null, composites, labels, qstarts, n);
			System.err.println(String.format("%-24s %8s %10.4f", "composite", "-", 10*top1d[0]+conc));
		}
		for(int ni=0; ni<netNames.size(); ni++){
			double[] top1d=getTop1Depths(netPreds.get(ni), null, labels, qstarts, n);
			double conc=getConcordance(netPreds.get(ni), null, labels, qstarts, n);
			long bytes=netFiles[ni].length();
			System.err.println(String.format("%-24s %7.1f %10.4f", netNames.get(ni), bytes/1024.0, 10*top1d[0]+conc));
		}
	}

	static void evalSelection(String name, double[] preds, ArrayList<Float> composites,
			ArrayList<Float> labels, ArrayList<Integer> qstarts, int n){
		int exact=0, genus=0, queries=0;
		double sumDepth=0;
		final boolean oracle=name.equals("ORACLE");
		for(int qi=0; qi<qstarts.size(); qi++){
			final int s=qstarts.get(qi);
			final int e=(qi+1<qstarts.size() ? qstarts.get(qi+1) : n);
			if(s>=n || e<=s){continue;}
			int pick=s;
			if(oracle){
				for(int i=s; i<e; i++){if(labels.get(i)>labels.get(pick)){pick=i;}}
			}else if(preds!=null){
				for(int i=s; i<e; i++){if(preds[i]>preds[pick]){pick=i;}}
			}else if(composites!=null){
				for(int i=s; i<e; i++){if(composites.get(i)>composites.get(pick)){pick=i;}}
			}
			final float d=labels.get(pick);
			sumDepth+=d;
			if(d>=0.999f){exact++;}
			if(d>=0.799f){genus++;}
			queries++;
		}
		System.err.println(String.format("%-24s %10.2f %10.2f %12.4f",
			name, 100.0*exact/Math.max(queries,1), 100.0*genus/Math.max(queries,1),
			sumDepth/Math.max(queries,1)));
	}

	static void evalWDCG(String name, double[] preds, ArrayList<Float> composites,
			ArrayList<Float> labels, ArrayList<Integer> qstarts, int n){
		double sumScore=0, sumOracle=0;
		int queries=0;
		for(int qi=0; qi<qstarts.size(); qi++){
			final int s=qstarts.get(qi);
			final int e=(qi+1<qstarts.size() ? qstarts.get(qi+1) : n);
			if(s>=n || e<=s){continue;}
			final int len=e-s;
			// Get indices sorted by ranker score (descending)
			final Integer[] idx=new Integer[len];
			for(int i=0; i<len; i++){idx[i]=i;}
			java.util.Arrays.sort(idx, (a,b)->{
				final double va=(preds!=null ? preds[s+a] : composites.get(s+a));
				final double vb=(preds!=null ? preds[s+b] : composites.get(s+b));
				return Double.compare(vb, va);
			});
			// Get indices sorted by truth (descending) for oracle
			final Integer[] oidx=new Integer[len];
			for(int i=0; i<len; i++){oidx[i]=i;}
			java.util.Arrays.sort(oidx, (a,b)->Float.compare(labels.get(s+b), labels.get(s+a)));

			double score=0, oracle=0;
			for(int i=0; i<len; i++){
				final double w=10.0/(i+1);
				score+=labels.get(s+idx[i])*w;
				oracle+=labels.get(s+oidx[i])*w;
			}
			sumScore+=score;
			sumOracle+=oracle;
			queries++;
		}
		final double norm=(sumOracle>0) ? sumScore/sumOracle : 0;
		System.err.println(String.format("%-24s %10.4f %12.4f %12.4f",
			name, norm, sumScore/Math.max(queries,1), sumOracle/Math.max(queries,1)));
	}

	static double[] getTop1Depths(double[] preds, ArrayList<Float> composites,
			ArrayList<Float> labels, ArrayList<Integer> qstarts, int n){
		int queries=0;
		double sumDepth=0;
		for(int qi=0; qi<qstarts.size(); qi++){
			final int s=qstarts.get(qi);
			final int e=(qi+1<qstarts.size() ? qstarts.get(qi+1) : n);
			if(s>=n || e<=s){continue;}
			int pick=s;
			if(preds!=null){
				for(int i=s; i<e; i++){if(preds[i]>preds[pick]){pick=i;}}
			}else if(composites!=null){
				for(int i=s; i<e; i++){if(composites.get(i)>composites.get(pick)){pick=i;}}
			}
			sumDepth+=labels.get(pick);
			queries++;
		}
		return new double[]{sumDepth/Math.max(queries,1)};
	}

	static double getConcordance(double[] preds, ArrayList<Float> composites,
			ArrayList<Float> labels, ArrayList<Integer> qstarts, int n){
		long conc=0, disc=0, tiePred=0;
		for(int qi=0; qi<qstarts.size(); qi++){
			final int s=qstarts.get(qi);
			final int e=Math.min(qi+1<qstarts.size() ? qstarts.get(qi+1) : n, n);
			if(s>=n){break;}
			for(int a=s; a<e; a++){
				for(int b=a+1; b<e; b++){
					final float la=labels.get(a), lb=labels.get(b);
					if(la==lb){continue;}
					final double pa=(preds!=null ? preds[a] : composites.get(a));
					final double pb=(preds!=null ? preds[b] : composites.get(b));
					if(pa==pb){tiePred++; continue;}
					if((pa>pb)==(la>lb)){conc++;}else{disc++;}
				}
			}
		}
		final double total=conc+disc+tiePred;
		return (conc+0.5*tiePred)/Math.max(total,1);
	}

	static void evalConcordance(String name, double[] preds, ArrayList<Float> composites,
			ArrayList<Float> labels, ArrayList<Integer> qstarts, int n){
		long conc=0, disc=0, tiePred=0;
		for(int qi=0; qi<qstarts.size(); qi++){
			final int s=qstarts.get(qi);
			final int e=Math.min(qi+1<qstarts.size() ? qstarts.get(qi+1) : n, n);
			if(s>=n){break;}
			for(int a=s; a<e; a++){
				for(int b=a+1; b<e; b++){
					final float la=labels.get(a), lb=labels.get(b);
					if(la==lb){continue;}
					final double pa=(preds!=null ? preds[a] : composites.get(a));
					final double pb=(preds!=null ? preds[b] : composites.get(b));
					if(pa==pb){tiePred++; continue;}
					if((pa>pb)==(la>lb)){conc++;}else{disc++;}
				}
			}
		}
		final double total=conc+disc+tiePred;
		final double withTies=(conc+0.5*tiePred)/Math.max(total,1);
		final double noTies=conc/Math.max((double)(conc+disc),1);
		System.err.println(String.format("%-24s %10.4f %10.4f %7.1f%% %10d",
			name, withTies, noTies, 100.0*tiePred/Math.max(total,1), (long)total));
	}
}
