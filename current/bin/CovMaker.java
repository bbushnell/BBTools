package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import structures.FloatList;
import structures.IntList;

/**
 * Utility for converting, condensing, and optimizing coverage data.
 * Performs SIMD-accelerated Similarity Merging and Entropy-based sorting.
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 21, 2026
 */
public class CovMaker {

	public static void main(String[] args){
		Timer t=new Timer();
		CovMaker x=new CovMaker(args);
		t.outstream=x.outstream;
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CovMaker(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLoader3.MAX_SAMPLES=10000; //Override memory limits
		//		SamLoader3.MAX_CONCURRENT_FILES=16;

		loader=new DataLoader(outstream);

		{//Parse the arguments
			final Parser parser=new Parser();
			for(int i=0; i<args.length; i++){
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}

				if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
				}else if(a.equals("condense") || a.equals("samples")){
					if(b!=null && Tools.startsWithLetter(b)) {
						boolean x=Parse.parseBoolean(b);
						assert(!x) : arg;
						if(!x) {condense=-1;}
					}
					condense=Integer.parseInt(b);
				}else if(a.equals("minseed")){
					Binner.minSizeToCompare=Binner.minSizeToMerge=Parse.parseIntKMG(b);
				}else if(a.equals("permute") || a.equals("sort") || a.equals("reorder")){
					reorder=Parse.parseBoolean(b);
				}else if(a.equals("cosine") || a.equals("cos")){
					useCosine=Parse.parseBoolean(b);
				}else if(a.equals("negcosine") || a.equals("negcos")){
					negCosine=Parse.parseBoolean(b);
				}else if(a.equals("lognormcosine") || a.equals("lognormcos") || a.equals("lognorm")){
					logNormCosine=Parse.parseBoolean(b);
				}else if(a.equals("magnitude") || a.equals("mag")){
					useMagnitude=Parse.parseBoolean(b);
				}else if(a.equals("entropy") || a.equals("ent")){
					useEntropy=Parse.parseBoolean(b);
				}else if(a.equals("rootmagnitude") || a.equals("rootmag")){
					rootMagnitude=Parse.parseBoolean(b);
				}else if(a.equals("magpower") || a.equals("magnitudepower")){
					magnitudePower=Float.parseFloat(b);
				}else if(a.equals("entpower") || a.equals("entropypower")){
					entropyPower=Float.parseFloat(b);
				}else if(a.equals("compare")){
					maxContigsToCompare=Parse.parseIntKMG(b);
				}else if(a.equals("out") || a.equals("outcov") || a.equals("covout")){
					out=b;
				}else if(a.equalsIgnoreCase("trackcardinality") || a.equalsIgnoreCase("loglog")) {
					SamLoader3.CARDINALITY=Parse.parseBoolean(b);
				}else if(a.equalsIgnoreCase("cardinality")) {
					cardinality=Parse.parseKMG(b);
				}else if(a.equalsIgnoreCase("maxsamples")) {
					SamLoader3.MAX_SAMPLES=Integer.parseInt(b);
				}else if(a.equalsIgnoreCase("maxconcurrentfiles") || a.equalsIgnoreCase("concurrentfiles")
					|| a.equalsIgnoreCase("readthreads")) {
					SamLoader3.MAX_CONCURRENT_FILES=Integer.parseInt(b);
				}else if(a.equals("in")){
					parser.in1=b;
				}else if(a.equals("ref")){
					outstream.println("Adding ref file "+b);
					ref=b;
				}else if(a.equals("refout") || a.equals("outref") || a.equals("outr")){
					outRef=b;
				}else if(loader.parse(arg, a, b)){
					//do nothing
				}else if(parser.parse(arg, a, b)){
					//do nothing
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}

			overwrite=parser.overwrite;
			append=parser.append;

			if(parser.in1!=null) {
				if(DataLoader.looksLikeCovFile(parser.in1)){
					loader.covIn=parser.in1;
				}else{
					for(String s : parser.in1.split(",")) {loader.readFiles.add(s);}
				}
			}
		}

		loader.checkInput();
	}

	void process(Timer t){
		outstream.println("Loading data...");
		Timer t2=new Timer(outstream);
		ArrayList<Contig> contigs=null;
		
		if(ref!=null && loader.covIn!=null) {
			contigs=loader.loadContigSubsetFromCov(ref, loader.covIn, loader.minContigToLoad);
		}else {
			if(ref!=null) {
				contigs=DataLoader.loadContigsFromFasta(ref, loader.minContigToLoad, false);
				t.stopAndStart("Loaded contigs:");
			}else if(loader.covIn==null) {
				contigs=loader.loadContigsFromSam(loader.readFiles.get(0));
				t.stopAndStart("Loaded header:");
			}
			contigs=loader.loadDepth(contigs, false, false);
		}

		if(condense>0 && loader.numDepths>condense){
			condenseSamples(contigs, condense);
			t.stopAndStart("Condensed samples:");
		}

		if(reorder){
			reorderSamples(contigs);
			t.stopAndStart("Reordered:");
		}

		if(out!=null){
			DataLoader.writeCov(out, contigs, loader.numDepths, outstream);
		}
		t2.stop("Total time:");
	}

	/**
	 * Merges samples using Log-Normalized Cosine Similarity on a proxy set.
	 * Optimized to pre-calculate logs and update normalization dynamically.
	 */
	private void condenseSamples(ArrayList<Contig> contigs, int target){
		final int current=loader.numDepths;
		if(current<=target){return;}

		outstream.println("Condensing "+current+" samples to "+target+" using Log-Weighted Metrics.");
		Timer t=new Timer();

		//1. Build Proxy Matrix (Linear)
		final int proxyCount=Math.min(contigs.size(), maxContigsToCompare);
		ArrayList<Contig> proxyContigs=(contigs.size()>proxyCount ? 
			(ArrayList<Contig>)contigs.clone() : contigs);
		if(contigs.size()>proxyCount){Collections.sort(proxyContigs);}

		float[][] proxyMatrix=new float[current][proxyCount];
		float[] magnitudes=new float[current];
		float[] entropy=new float[current];

		//New: Track dynamic normalization factors
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
//			System.err.println(Arrays.toString(currentInvNorms));

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
		//5. Replay Merges on Full Data
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

	private float entropy(ArrayList<Contig> contigs, int sample, int limit) {
		double entropy=0;
		int count=0;
		for(int i=0; i<limit; i++){
			Contig c=contigs.get(i);
			if(c.size()>=Binner.minSizeToCompare) {
				float d=c.depth(sample);
				if(d>0){
					entropy+=Math.log(d+1);
					count++;
				}
			}
		}
		return (float)entropy;
	}

	/**
	 * Sorts samples by Entropy to put high-info samples at indices 0, 1, 2.
	 * Fixes indexing performance in Binner.
	 */
	private void reorderSamples(ArrayList<Contig> contigs){
		outstream.println("Permuting samples by Entropy...");
		int samples=loader.numDepths;
		double[] entropy=new double[samples];

		//Calculate Entropy using Proxy (Top 100k) is sufficient and faster
		int limit=Math.min(contigs.size(), maxContigsToEntropy);
		for(int i=0; i<limit; i++){
			Contig c=contigs.get(i);
			if(c.size()>=Binner.minSizeToCompare) {
				for(int s=0; s<samples; s++){
					float d=c.depth(s);
					if(d>0){
						entropy[s]+=Math.log(d+1);
					}
				}
			}
		}

		Integer[] idxs=new Integer[samples];
		for(int i=0; i<samples; i++){idxs[i]=i;}
		Arrays.sort(idxs, (a, b) -> Double.compare(entropy[b], entropy[a])); //Desc

		outstream.println("Sample Priority: "+Arrays.toString(idxs));

		//Apply permutation
		for(Contig c : contigs){
			FloatList old=new FloatList(c.depthList());
			c.clearDepth();
			for(int i=0; i<samples; i++){
				c.appendDepth(old.get(idxs[i]));
			}
		}
	}

	private DataLoader loader;
	private String out=null;
	private String ref=null;
	private String outRef=null;
	private int condense=-1;
	private int maxContigsToCompare=100000;
	private int maxContigsToEntropy=100000;
	private boolean reorder=true;
	private boolean useCosine=true;
	private boolean negCosine=false;
	private boolean useMagnitude=true;
	private boolean useEntropy=false;
	private boolean rootMagnitude=false;
	private boolean logNormCosine=false;
	private float magnitudePower=1f;
	private float entropyPower=0.25f;
	long cardinality=0;

	private boolean overwrite=true;
	private boolean append=false;

	public static boolean verbose=false;
	private PrintStream outstream=System.err;
}