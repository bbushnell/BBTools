package bin;

import java.util.ArrayList;
import shared.Random;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.FloatList;
import tracker.EntropyTracker;

/**
 * Calculates and compensates for entropy variations across different GC content
 * levels in genomic sequences. Generates random sequences with specified GC content
 * and measures their k-mer entropy to create compensation references.
 * @author Brian Bushnell
 */
public class AdjustEntropy {
	
	/**
	 * Program entry point for entropy analysis and compensation.
	 * Generates random sequences across GC content range, measures entropy,
	 * and optionally processes input file with compensation metrics.
	 * @param args Command-line arguments: k, window, step, length, trials, [filename]
	 */
	public static void main(String[] args) {

		int k=Integer.parseInt(args[0]);
		int window=Integer.parseInt(args[1]);
		float step=Float.parseFloat(args[2]);
		int length=(args.length<4 ? 100000 : Integer.parseInt(args[3]));
		int trials=(args.length<5 ? 9 : Integer.parseInt(args[4]));
		assert(step>0 && step<1);
		assert(length>window);
		assert(window>k);
		assert(k>0 && k<=10);
		
		AdjustEntropy ea=new AdjustEntropy();
		EntropyTracker et=new EntropyTracker(k, window, false);
		FloatList fl=new FloatList();
		System.out.println("#K\t"+k);
		System.out.println("#Window\t"+window);
		System.out.println("#Length\t"+length);
		System.out.println("#Step\t"+step);
		System.out.println("#Trials\t"+trials);
		System.out.println("#GC\tAvgEnt\tMaxEnt");
		for(double gc=0; gc<1 || Tools.absdif(gc, 1)<0.5*step; gc+=step) {
			float entropy=ea.randomSequenceEntropy(length, (float)gc, et, trials);
			fl.add(entropy);
			System.out.println(String.format("%.4f\t%.6f\t%.6f", gc, entropy, ea.max));
		}
		
		String fname=args.length<6 ? null : args[5];
		entropyArray=fl.toArray();
		if(fname!=null) {
			processFile(fname, et);
		}
	}
	
	/**
	 * Processes FASTA file to calculate GC-compensated entropy and strandedness
	 * metrics for each sequence read.
	 * @param fname Input FASTA filename to process
	 * @param et EntropyTracker configured with appropriate k-mer parameters
	 */
	static void processFile(String fname, EntropyTracker et) {
		ArrayList<Read> reads=ReadInputStream.toReads(fname, FileFormat.FASTA, 2000);
		System.out.println("#GC\tgcCompEntropy\tstrandedness");
		int[] counts=new int[1<<(2*et.k())];
		for(Read r : reads) {
			float gc=r.gc();
			float entropy=et.averageEntropy(r.bases, false);
			float comp=compensate(gc, entropy);
			float strandedness=EntropyTracker.strandedness(r.bases, counts, et.k());
			System.out.println(String.format("%.4f\t%.6f\t%6f", gc, comp, strandedness));
		}
	}
	
	/**
	 * Calculates maximum expected entropy for a given GC content using linear
	 * interpolation between adjacent values in the reference entropy array.
	 *
	 * @param gc GC content fraction (0.0 to 1.0)
	 * @param array Reference entropy values indexed by GC content steps
	 * @return Interpolated maximum entropy value for the specified GC content
	 */
	//⚠️ NN-INPUT (compensate→SpectraCounter:223/Clade:255→CellNet) - DOCUMENT+FLAG ONLY, do NOT patch.
	//TODO: Possible bug [bin/AdjustEntropy#001] - NN-FROZEN: at an exact GC bin boundary (gc*steps integer → lowBin==highBin) both interpolation fractions are 0, so max=(0*lowE+0*highE)*steps=0 instead of array[lowBin]. ALWAYS for gc=0 and gc=1; also interior boundaries (e.g. gc=0.5 when steps is even) → compensate then returns min(1,1+entropy)=1. Wrong NN input for boundary-GC contigs. Fix (retrain only): if(lowBin==highBin) return array[lowBin]. See report.
	static float maxEntropy(float gc, float[] array) {
		int steps=array.length-1;
		float stepSize=1f/steps;
		int lowBin=(int)Math.floor(gc*steps);
		int highBin=(int)Math.ceil(gc*steps);
		float lowE=array[lowBin], highE=array[highBin];
		float lowGC=stepSize*lowBin, highGC=stepSize*highBin;
		float lowFraction=highGC-gc, highFraction=gc-lowGC;
		float max=(lowFraction*lowE+highFraction*highE)*steps;
//		assert(false) : "GC "+gc+" -> "+max;
		return max;
	}
	
	/** Loads default entropy compensation data using k=4 and window=150. */
	public static synchronized void load() {
		load(4, 150);
	}
	
	/**
	 * Loads entropy compensation data for specified k-mer length and window size.
	 * Constructs filename pattern and locates entropy reference file.
	 * @param k K-mer length for entropy calculation
	 * @param window Sliding window size for entropy measurement
	 */
	public static synchronized void load(int k, int window) {
		String fname="?entropy_k"+k+"_w"+window+".tsv";
		fname=Data.findPath(fname);
		setEntropyFile(fname);
	}
	
	/**
	 * Sets the entropy reference file if different from currently loaded file.
	 * Prevents redundant file loading by checking filename against cached value.
	 * @param fname Entropy reference file path to load
	 */
	private static synchronized void setEntropyFile(String fname) {
		assert(fname!=null);
		if(fnameLoaded==null || !fnameLoaded.equals(fname)) {
			entropyArray=loadEntropyFile(fname);
		}
	}
	
	/**
	 * Parses entropy reference file containing GC-entropy mapping data.
	 * Extracts k-mer and window parameters from header lines and entropy
	 * values from data lines.
	 *
	 * @param fname Path to entropy reference file
	 * @return Array of entropy values indexed by GC content steps
	 */
	private static float[] loadEntropyFile(String fname) {
		ArrayList<byte[]> lines=ByteFile.toLines(fname);
		FloatList floats=new FloatList(lines.size());
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line : lines) {
			lp.set(line);
			if(lp.startsWith('#')){
				if(lp.startsWith("#K\t")) {
					kLoaded=lp.parseInt(1);
				}else if(lp.startsWith("#Window\t")) {
					wLoaded=lp.parseInt(1);
				}
			}else {
				//TODO: Possible bug [bin/AdjustEntropy#002] - QUESTION (NN-frozen, don't patch): loads column 1 = AvgEnt (main() writes a tab-separated #GC/AvgEnt/MaxEnt header = cols 0,1,2), yet maxEntropy()/its local call the loaded array "max"; the real MaxEnt is column 2. Either intentional (compensate vs avg-random, misleading name) or a column off-by-one. Float to Brian - changing the column changes NN input. See report.
				floats.add(lp.parseFloat(1));
			}
		}
		return floats.toArray();
	}
	
	/** Returns entropy as a fraction of random entropy for this GC level */
	//⚠️ NN-INPUT (feeds c.entropy → Oracle vector via SpectraCounter; gcCompEntropy in clade/Clade) - DOCUMENT ONLY, do NOT patch the formula (the net trained on it).
	//claim: returns 1-(max-entropy) clamped to <=1 - how close the contig entropy is to the GC-matched random reference; entropy>=max clamps to 1. Requires entropyArray pre-loaded (callers call load() first, else NPE in maxEntropy).
	public static float compensate(float gc, float entropy) {
		float max=maxEntropy(gc, entropyArray);
		return Tools.min(1, 1-(max-entropy));//entropy/max;
	}
	
	/**
	 * Generates multiple random sequences and calculates statistical entropy
	 * measures. Uses top quartile samples to compute robust entropy estimate
	 * that avoids low outliers.
	 *
	 * @param len Length of random sequences to generate
	 * @param gc Target GC content fraction
	 * @param et EntropyTracker for calculating k-mer entropy
	 * @param trials Number of random sequences to generate for statistics
	 * @return Average entropy from top-performing samples
	 */
	//claim: returns the mean of the top max(sqrt(trials), trials/4) of the trials sorted samples (robust upper estimate, drops low outliers); asserts min<=mid<=max. Offline (table generation), not NN runtime.
	float randomSequenceEntropy(int len, float gc, EntropyTracker et, int trials) {
		FloatList fl=new FloatList(trials);
		for(int i=0; i<trials; i++) {fl.add(randomSequenceEntropy(len, gc, et));}
		fl.sort();
		min=fl.get(0);
		max=fl.lastElement();
		assert(min<=max) : fl;
		
		double sum=0;
		int count=0;
		int samples=Tools.max((int)Math.sqrt(trials), trials/4);
		for(int i=0; i<samples; i++) {
			sum+=fl.get(fl.size()-i-1);
			count++;
		}
		assert(count==samples);
		mid=(float)(sum/samples);
			
//		else if((trials&1)==1) {//odd
//			double sum=0;
//			int count=0;
//			int center=trials/2;
//			for(int i=0; i<=trials/4; i++) {
//				sum+=fl.get(center+i);
//				sum+=fl.get(center-i);
//				count+=2;
//			}
//			mid=(float)(sum/count);
//		}else {
//			double sum=0;
//			int count=0;
//			int center=trials/2-1;
//			for(int i=0; i<=trials/4; i++) {
//				sum+=fl.get(center+i+1);
//				sum+=fl.get(center-i);
//				count+=2;
//			}
//			mid=(float)(sum/count);
//		}
		assert(min<=max) : fl;
		assert(min<=mid) : fl;
		assert(mid<=max) : fl;
		return mid;
	}
	
	/**
	 * Generates single random sequence with specified GC content and measures
	 * its k-mer entropy using the provided EntropyTracker.
	 *
	 * @param len Length of random sequence to generate
	 * @param gc Target GC content fraction
	 * @param et EntropyTracker for calculating k-mer entropy
	 * @return Entropy value of the generated random sequence
	 */
	static float randomSequenceEntropy(int len, float gc, EntropyTracker et) {
		byte[] bases=randomSequence(len, gc);
		return et.averageEntropy(bases, false);
	}
	
	/**
	 * Generates random DNA sequence with specified length and GC content.
	 * Uses probabilistic base selection where high GC bases (C,G) are chosen
	 * with probability 'gc' and low GC bases (A,T) with probability '1-gc'.
	 *
	 * @param len Length of sequence to generate
	 * @param gc Target GC content fraction (0.0 to 1.0)
	 * @return Byte array containing random DNA sequence
	 */
	//claim: resulting GC == gc. high=2 (C/G at atcg[2,3]) chosen with prob gc (nextFloat()>=gc has prob 1-gc → high=0=A/T); low bit picks within the pair. Correct at gc=0 (all AT) / gc=1 (all CG). nextInt()&1 unbiased even for negative ints.
	static byte[] randomSequence(int len, float gc) {
		byte[] bases=new byte[len];
		Random randy=Shared.threadLocalRandom();
		byte[] atcg={'A','T','C','G'};
		for(int i=0; i<len; i++) {
			int high=(randy.nextFloat()>=gc ? 0 : 2);
			int low=(randy.nextInt()&1);
			byte b=atcg[high+low];
			bases[i]=b;
		}
		return bases;
	}
	
	float min, mid, max;
	/** Filename of currently loaded entropy reference file, null if none loaded */
	private static String fnameLoaded=null;
	/** K-mer length parameter from loaded entropy reference file */
	public static int kLoaded=0;
	/** Window size parameter from loaded entropy reference file */
	public static int wLoaded=0;
	/**
	 * Array of entropy values indexed by GC content steps for compensation calculations
	 */
	private static float[] entropyArray=null;
	
}
