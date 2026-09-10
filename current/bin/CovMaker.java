package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.FloatList;

/**
 * Utility for converting, condensing, and optimizing coverage data.
 * Performs SIMD-accelerated Similarity Merging and Entropy-based sorting.
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 21, 2026
 */
public class CovMaker {
	//REVIEW (Eru 2026-06-19): standalone coverage utility (covmaker.sh) — condense/reorder/convert cov data. NOT
	//NN-input (its cosine merge clusters SAMPLES, not contigs). Fix #001: condense=f NumberFormatException (parse).
	//Dead: entropy(ArrayList,int,int) (361). +verified-clever: harmonic inverse-norm O(1) merge update (262).

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
		condenser=new SampleCondenser(loader, outstream);//Shared estimator+merge engine (also used by QuickBin autocondense)

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
					if(b!=null && b.equalsIgnoreCase("auto")){
						condense=AUTO;
					}else if(b!=null && Tools.startsWithLetter(b)) {
						boolean x=Parse.parseBoolean(b);
						assert(!x) : arg;
						if(!x) {condense=-1;}
					}else{
						//FIX [bin/CovMaker#001]: was an UNCONDITIONAL condense=Integer.parseInt(b) after the boolean
						//block, so the intended `condense=f`/`samples=f` (disable -> condense=-1) threw NumberFormatException
						//on the "f". Now guarded by else: letter -> boolean path, number -> parseInt. (Verified predict-then-run.)
						condense=Integer.parseInt(b);
					}
				}else if(a.equals("gapratio")){
					condenser.gapRatio=Float.parseFloat(b);
				}else if(a.equals("alwaysmerge") || a.equals("alwaysmerger")){
					condenser.alwaysMergeR=Float.parseFloat(b);
				}else if(a.equals("nevermerge") || a.equals("nevermerger")){
					condenser.neverMergeR=Float.parseFloat(b);
				}else if(a.equals("minseed")){
					Binner.minSizeToCompare=Binner.minSizeToMerge=Parse.parseIntKMG(b);
				}else if(a.equals("permute") || a.equals("sort") || a.equals("reorder")){
					reorder=Parse.parseBoolean(b);
				}else if(a.equals("cosine") || a.equals("cos")){
					condenser.useCosine=Parse.parseBoolean(b);
				}else if(a.equals("negcosine") || a.equals("negcos")){
					condenser.negCosine=Parse.parseBoolean(b);
				}else if(a.equals("lognormcosine") || a.equals("lognormcos") || a.equals("lognorm")){
					condenser.logNormCosine=Parse.parseBoolean(b);
				}else if(a.equals("magnitude") || a.equals("mag")){
					condenser.useMagnitude=Parse.parseBoolean(b);
				}else if(a.equals("entropy") || a.equals("ent")){
					condenser.useEntropy=Parse.parseBoolean(b);
				}else if(a.equals("rootmagnitude") || a.equals("rootmag")){
					condenser.rootMagnitude=Parse.parseBoolean(b);
				}else if(a.equals("magpower") || a.equals("magnitudepower")){
					condenser.magnitudePower=Float.parseFloat(b);
				}else if(a.equals("entpower") || a.equals("entropypower")){
					condenser.entropyPower=Float.parseFloat(b);
				}else if(a.equals("compare")){
					condenser.maxContigsToCompare=Parse.parseIntKMG(b);
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

		if(condense==AUTO){
			condense=condenser.chooseTarget(contigs);
			t.stopAndStart("Chose condense target:");
		}
		if(condense>0 && loader.numDepths>condense){
			condenser.condense(contigs, condense);
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

	//condenseSamples/chooseCondenseTarget and their helpers moved VERBATIM to bin.SampleCondenser
	//(2026-09-09) so QuickBin's in-pipeline autocondense shares the exact implementation; this class
	//now only parses flags into the condenser and delegates in process().

	//DEAD: no callers (reorderSamples has its own inline entropy loop at 388-398). Deletion candidate.
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
	/** condense=auto sentinel; resolved to a real target by SampleCondenser.chooseTarget. */
	private static final int AUTO=-2;
	/** Estimator + merge engine shared with QuickBin's autocondense; parse writes its knobs directly. */
	private final SampleCondenser condenser;
	private int maxContigsToEntropy=100000;
	private boolean reorder=true;
	long cardinality=0;

	private boolean overwrite=true;
	private boolean append=false;

	public static boolean verbose=false;
	private PrintStream outstream=System.err;
}