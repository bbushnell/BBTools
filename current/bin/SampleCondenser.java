package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import shared.Tools;
import simd.Vector;
import structures.FloatList;
import structures.IntList;

/**
 * Estimates the number of LOGICAL samples in multi-sample depth data and
 * condenses correlated samples into logical groups by summation.
 * Extracted from CovMaker so that both CovMaker (covmaker.sh preprocessing)
 * and QuickBin (in-pipeline autocondense) share ONE implementation.
 * Correlated samples inflate coverage evidence: near-duplicate libraries make
 * independence-assuming statistics count one signal N times (measured:
 * 8-sample trap = 2 logical x 4 near-dupes collapsed warpbin 976 to 108 Total
 * and bent QuickBin 938 to 873; condensing to 2 healed both, qbwb_bench
 * 2026-09-09).  Summing grouped columns also averages noise, which measured
 * BETTER than an equal number of independent raw samples.
 *
 * @author Brian Bushnell
 * @author Amber
 * @date September 9, 2026
 */
public class SampleCondenser {

	public SampleCondenser(DataLoader loader_, PrintStream outstream_){
		loader=loader_;
		outstream=outstream_;
	}

	/**
	 * Estimates the number of LOGICAL samples for condense=auto: Pearson correlation
	 * of log-depths over the proxy set, single-linkage merge distances (d=1-r), and
	 * the largest relative jump (elbow) between within-group and between-group merges.
	 * Guard band: pairs with r>alwaysMergeR always condense; pairs with r<neverMergeR
	 * never do. Validated 32/32 on ground-truth grids (group structures x jitter<=0.3
	 * x zeros<=90%) plus real correlated arms (qbwb_bench 2026-09-09); spectral
	 * alternatives (PR/erank) rejected — they underestimate uneven group structures.
	 * @return Target sample count in [1, numDepths]; numDepths means no condensation.
	 */
	public int chooseTarget(ArrayList<Contig> contigs){
		final int n=loader.numDepths;
		if(n<2){return n;}
		final int proxyCount=Math.min(contigs.size(), maxContigsToCompare);
		ArrayList<Contig> proxy=(contigs.size()>proxyCount ?
			(ArrayList<Contig>)contigs.clone() : contigs);
		if(contigs.size()>proxyCount){Collections.sort(proxy);}

		//Log-depth matrix over the proxy set; presence tracks depth>0 because rows
		//where BOTH samples are absent are EXCLUDED from that pair's correlation:
		//shared absence is not evidence of similarity (NEON-like data is mostly
		//zeros in every library, and co-absence would otherwise read as high r and
		//wrongly condense independent sparse samples). Present-in-one-absent-in-
		//other rows ARE included - informative discordance.
		final double[][] x=new double[n][proxyCount];
		final boolean[][] present=new boolean[n][proxyCount];
		for(int c=0; c<proxyCount; c++){
			Contig ctg=proxy.get(c);
			for(int s=0; s<n; s++){
				final float d=ctg.depth(s);
				x[s][c]=Math.log(d+0.5);
				present[s][c]=(d>0);
			}
		}
		//Pairwise Pearson r on co-informative rows; single-linkage edges (distance, i, j).
		final double[][] edges=new double[n*(n-1)/2][3];
		int e=0;
		for(int i=0; i<n; i++){
			for(int j=i+1; j<n; j++){
				double sw=0, sx=0, sy=0, sxx=0, syy=0, sxy=0;
				for(int c=0; c<proxyCount; c++){
					if(!present[i][c] && !present[j][c]){continue;}
					final double xi=x[i][c], yj=x[j][c];
					sw++; sx+=xi; sy+=yj; sxx+=xi*xi; syy+=yj*yj; sxy+=xi*yj;
				}
				double r=0;
				if(sw>0){
					double cov=sxy/sw-(sx/sw)*(sy/sw);
					double vx=sxx/sw-(sx/sw)*(sx/sw), vy=syy/sw-(sy/sw)*(sy/sw);
					double denom=Math.sqrt(vx*vy);
					r=(denom>0 ? cov/denom : 0);
				}
				edges[e++]=new double[] {1-r, i, j};
			}
		}
		Arrays.sort(edges, (p, q) -> Double.compare(p[0], q[0]));

		//Kruskal-style single linkage: record the n-1 component-joining merge distances.
		final int[] parent=new int[n];
		for(int i=0; i<n; i++){parent[i]=i;}
		final double[] merges=new double[n-1];
		int m=0;
		for(double[] ed : edges){
			int a=findRoot(parent, (int)ed[1]), b=findRoot(parent, (int)ed[2]);
			if(a!=b){parent[a]=b; merges[m++]=ed[0];}
		}
		assert(m==n-1) : m+", "+n;

		//Elbow: largest relative jump; guards clamp the cut into the allowed band.
		final double EPS=1e-4;
		int cut=-1;
		double best=gapRatio;
		for(int i=0; i<m-1; i++){
			double ratio=(merges[i+1]+EPS)/(merges[i]+EPS);
			if(ratio>best){best=ratio; cut=i;}
		}
		int k=(cut<0 ? 0 : cut+1);//Number of merges accepted
		int mandatory=0, allowed=0;
		while(mandatory<m && merges[mandatory]<=1-alwaysMergeR){mandatory++;}
		while(allowed<m && merges[allowed]<1-neverMergeR){allowed++;}
		k=Tools.mid(mandatory, k, allowed);
		final int target=n-k;

		StringBuilder sb=new StringBuilder("Merge distances (1-r):");
		for(int i=0; i<m; i++){sb.append(String.format(" %.4f", merges[i]));}
		outstream.println(sb.toString());
		outstream.println("Auto condense target: "+target+" of "+n+" samples"
			+(cut<0 ? " (no elbow found)" : String.format(" (elbow ratio %.1f)", best)));
		return target;
	}

	/**
	 * Merges samples down to target using Log-Normalized Cosine Similarity on a
	 * proxy set, then replays the merges (column SUMS) on the full contig data.
	 * Optimized to pre-calculate logs and update normalization dynamically.
	 * Requires BinObject.invSampleDepthSum (fill via loader.calcDepthSum first);
	 * re-runs calcDepthSum on the condensed columns before returning.
	 */
	public void condense(ArrayList<Contig> contigs, int target){
		final int current=loader.numDepths;
		if(current<=target){return;}

		outstream.println("Condensing "+current+" samples to "+target+" using Log-Weighted Metrics.");
		shared.Timer t=new shared.Timer();

		//1. Build Proxy Matrix (Linear)
		final int proxyCount=Math.min(contigs.size(), maxContigsToCompare);
		ArrayList<Contig> proxyContigs=(contigs.size()>proxyCount ?
			(ArrayList<Contig>)contigs.clone() : contigs);
		if(contigs.size()>proxyCount){Collections.sort(proxyContigs);}

		float[][] proxyMatrix=new float[current][proxyCount];
		float[] magnitudes=new float[current];
		float[] entropy=new float[current];

		//Track dynamic normalization factors
		double[] currentInvNorms=Arrays.copyOf(BinObject.invSampleDepthSum, current);

		for(int sIdx=0; sIdx<current; sIdx++){
			for(int cIdx=0; cIdx<proxyCount; cIdx++){
				float d=proxyContigs.get(cIdx).depth(sIdx);
				proxyMatrix[sIdx][cIdx]=d;
				magnitudes[sIdx]+=d;
			}
			entropy[sIdx]=calcSampleEntropy(proxyMatrix[sIdx]);
		}

		//2. Build Log-Proxy Matrix (Pre-calculation optimization)
		float[][] logProxyMatrix=new float[current][];
		if(logNormCosine){
			for(int i=0; i<current; i++){
				logProxyMatrix[i]=transformToLog(proxyMatrix[i], currentInvNorms[i]);
			}
		}
		t.stopAndStart("Built proxy matrix:");

		//3. State Tracking
		IntList[] groups=new IntList[current];
		boolean[] dead=new boolean[current];
		for(int i=0; i<current; i++){
			groups[i]=new IntList();
			groups[i].add(i);
		}

		//4. Compute Cost Matrix
		float[][] costMatrix=new float[current][current];
		for(int i=0; i<current; i++){
			for(int j=i+1; j<current; j++){
				//Pass the pre-calculated log matrix
				costMatrix[i][j]=calculateCost((logNormCosine ? logProxyMatrix : proxyMatrix), magnitudes, entropy, i, j);
			}
		}
		t.stopAndStart("Built cost matrix:");

		//5. Greedy Merge Loop
		int active=current;
		while(active>target){
			int bestI=-1, bestJ=-1;
			float minCost=Float.MAX_VALUE;

			for(int i=0; i<current; i++){
				if(dead[i]){continue;}
				for(int j=i+1; j<current; j++){
					if(dead[j]){continue;}
					if(costMatrix[i][j]<minCost){
						minCost=costMatrix[i][j];
						bestI=i;
						bestJ=j;
					}
				}
			}

			if(bestI==-1){break;}

			//--- PERFORM MERGE ---

			//A. Update Linear Proxy Matrix
			for(int k=0; k<proxyCount; k++){
				proxyMatrix[bestI][k]+=proxyMatrix[bestJ][k];
			}

			//B. Update Magnitudes & Entropy
			magnitudes[bestI]+=magnitudes[bestJ];
			entropy[bestI]=calcSampleEntropy(proxyMatrix[bestI]);

			//C. Update Normalization Factor (Harmonic sum logic for inverse)
			//NewInv = (InvI * InvJ) / (InvI + InvJ) which equals 1/(SumI+SumJ)
			double invI=currentInvNorms[bestI];
			double invJ=currentInvNorms[bestJ];
			currentInvNorms[bestI]=(invI*invJ)/(invI+invJ);

			//D. Update Log-Proxy Matrix for the NEW merged row only
			if(logNormCosine){
				logProxyMatrix[bestI]=transformToLog(proxyMatrix[bestI], currentInvNorms[bestI]);
			}

			//E. Update Groups
			groups[bestI].addAll(groups[bestJ]);
			groups[bestJ]=null;
			dead[bestJ]=true;
			active--;

			//F. Update costs
			for(int k=0; k<current; k++){
				if(k!=bestI && !dead[k]){
					int row=Math.min(bestI, k);
					int col=Math.max(bestI, k);
					costMatrix[row][col]=calculateCost((logNormCosine ? logProxyMatrix : proxyMatrix), magnitudes, entropy, bestI, k);
				}
			}
		}
		t.stopAndStart("Calculated merge order:");

		applyMergesToContigs(contigs, groups, dead, target);
		t.stopAndStart("Applied merges:");
	}

	private static int findRoot(int[] p, int x){
		while(p[x]!=x){p[x]=p[p[x]]; x=p[x];}
		return x;
	}

	private float calculateCost(float[][] matrix, float[] mags, float[] ents, int i, int j){
		//Use the pre-calculated vectors directly!
		float cos=(useCosine ? (1.0f-Vector.cosineSimilarity(matrix[i], matrix[j])) : 1.0f);
		if(negCosine){cos=1.0001f-cos;}

		float magWeight=(useMagnitude ? (mags[i]+mags[j]) : 1.0f);
		if(magnitudePower!=1.0){magWeight=(float)Math.pow(magWeight, magnitudePower);}
		if(rootMagnitude){magWeight=(float)Math.sqrt(magWeight);}

		float entWeight=(useEntropy ? (ents[i]+ents[j]) : 1.0f);
		if(entropyPower!=1.0){entWeight=(float)Math.pow(entWeight, entropyPower);}

		return cos*magWeight*entWeight;
	}

	private float[] transformToLog(float[] rawDepths, double normalizationFactor){
		float[] transformed=new float[rawDepths.length];
		for(int k=0; k<rawDepths.length; k++){
			float d=(float)(rawDepths[k]*normalizationFactor+0.25f);
			transformed[k]=(float)Math.log(d);
		}
		return transformed;
	}

	private void applyMergesToContigs(ArrayList<Contig> contigs, IntList[] groups, boolean[] dead, int target) {
		//Replay Merges on Full Data
		IntList[] finalGroups=new IntList[target];
		int next=0;
		for(int i=0; i<dead.length; i++){
			if(!dead[i]){
				finalGroups[next++]=groups[i];
			}
		}

		outstream.println("Applying merges to "+contigs.size()+" contigs...");

		for(Contig c : contigs){
			FloatList oldDepths=c.depthList();
			float[] newDepths=new float[target];
			for(int g=0; g<target; g++){
				IntList cols=finalGroups[g];
				float sum=0;
				for(int k=0; k<cols.size(); k++){
					sum+=oldDepths.get(cols.get(k));
				}
				newDepths[g]=sum;
			}

			c.clearDepth();
			for(float f : newDepths){c.appendDepth(f);}
		}

		loader.numDepths=target;

		System.err.println("After merging:");
		double totalEntropy=loader.calcDepthSum(contigs);
		System.err.println("Depth Entropy:      \t"+String.format("%.4f", totalEntropy));
		System.err.println("Samples Equivalent: \t"+BinObject.samplesEquivalent);
		System.err.println("Samples:            \t"+loader.numDepths);
	}

	// Helper for the entropy calculation: sum of log(d+1)
	private float calcSampleEntropy(float[] proxyDepths) {
		double ent = 0;
		for (float d : proxyDepths) {
			if (d > 0) { ent += Math.log(d + 1); }
		}
		return (float)ent/proxyDepths.length;
	}

	/** Stop merging when the next merge distance exceeds gapRatio x the previous. */
	public float gapRatio=3.0f;
	/** Pairs correlated above this always condense, below neverMergeR never do. */
	public float alwaysMergeR=0.95f, neverMergeR=0.70f;
	public int maxContigsToCompare=100000;
	public boolean useCosine=true;
	public boolean negCosine=false;
	public boolean useMagnitude=true;
	public boolean useEntropy=false;
	public boolean rootMagnitude=false;
	public boolean logNormCosine=false;
	public float magnitudePower=1f;
	public float entropyPower=0.25f;

	private final DataLoader loader;
	private final PrintStream outstream;
}
