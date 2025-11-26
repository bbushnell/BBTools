package bin;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;

import clade.CladeObject;
import dna.AminoAcid;
import ml.CellNet;
import shared.Tools;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntHashSet;
import tax.TaxTree;

/**
 * Superclass for binner classes
 * @author Brian Bushnell
 * @date Feb 4, 2025
 *
 */
public class BinObject {
	
	/**
	 * Sets the quantization factor for k-mer sampling.
	 * Higher values reduce memory usage by only considering every nth canonical k-mer.
	 * @param x Quantization factor; must be positive
	 */
	public static void setQuant(int x) {
		quant=x;
		assert(quant>0);
		initialize();
	}
	
//	public static void setK(int x) {
//		assert(k>0 && k<16);
//		k=x;
//	}
	
	/** Initializes static data structures for k-mer analysis.
	 * Creates remap matrices, canonical k-mer mappings, and GC content mappings. */
	private static void initialize() {
		remapMatrix=makeRemapMatrix(2, 5, true);
		//K=1 is ACGTN for use in GC calcs.
		//K=2 is noncanonical for use in strandedness calcs.
		canonicalKmers=makeCanonicalKmers();
		invCanonicalKmers=makeInvCanonicalKmers();
		gcmapMatrix=makeGCMapMatrix();
	}
	
	/**
	 * Creates a matrix mapping k-mers to canonical indices for different k values.
	 * Each k-value gets its own remapping array to convert raw k-mers to canonical indices.
	 *
	 * @param mink Minimum k-mer length to process
	 * @param maxk Maximum k-mer length to process
	 * @param specialCase2 If true, uses non-canonical mapping for k=2 (for strandedness calculations)
	 * @return Matrix where [k][kmer] gives the canonical index for that k-mer
	 */
	private static synchronized int[][] makeRemapMatrix(int mink, int maxk, boolean specialCase2){
		int[][] matrix=new int[maxk+1][];
		for(int i=mink; i<=maxk; i++) {
			matrix[i]=makeRemap(i);
		}
		if(specialCase2 && 2<=maxk && 2>=mink) {
			matrix[2]=new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		}
		return matrix;
	}
	
	/**
	 * Creates a matrix mapping canonical k-mer indices to their GC content.
	 * Used for GC-based comparisons between sequences.
	 * @return Matrix where [k][canonical_index] gives GC content for that canonical k-mer
	 */
	private static synchronized int[][] makeGCMapMatrix(){
		int[][] matrix=new int[remapMatrix.length][];
		for(int i=0; i<matrix.length; i++) {
			int[] remap=remapMatrix[i];
			if(remap!=null) {
				matrix[i]=gcmap(i, remap);
			}
		}
		return matrix;
	}
	
	/**
	 * Calculates the number of canonical k-mers for each k value.
	 * Used to determine array sizes for frequency counting.
	 * @return Array where index k gives the number of canonical k-mers for that k
	 */
	private static synchronized int[] makeCanonicalKmers() {
		int[] array=new int[remapMatrix.length];
		for(int i=0; i<array.length; i++) {
			int[] remap=remapMatrix[i];
			int max=(remap==null ? 1 : Tools.max(remap)+1);
			array[i]=max;
		}
		return array;
	}
	
	/**
	 * Calculates the inverse of canonical k-mer counts for normalization.
	 * Pre-computed to avoid repeated division operations.
	 * @return Array where index k gives 1/canonicalKmers[k]
	 */
	private static synchronized float[] makeInvCanonicalKmers() {
		float[] array=new float[canonicalKmers.length];
		for(int i=0; i<array.length; i++) {
			array[i]=1f/canonicalKmers[i];
		}
		return array;
	}
	
	/**
	 * Creates a remapping array for a specific k-mer length.
	 * Maps each possible k-mer to its canonical representative index,
	 * considering quantization and reverse complement equivalence.
	 *
	 * @param k K-mer length
	 * @return Array where index is k-mer value and value is canonical index
	 */
	public static int[] makeRemap(int k){
		final int bits=2*k;
		final int max=(int)((1L<<bits)-1);
		int count=0;
		IntHashMap canonMap=new IntHashMap();
		IntHashMap kmerMap=new IntHashMap();
		for(int kmer=0; kmer<=max; kmer++){
//			int ungapped=ungap(kmer, k, gap);
			int canon=Tools.min(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k));
			if(canon%quant==0 && !canonMap.containsKey(canon)) {
				canonMap.put(canon, count);
				count++;
			}
			int idx=canonMap.get(canon);
			kmerMap.put(kmer, idx);
		}
		int[] remap=new int[max+1];
		Arrays.fill(remap, -1);
		for(int kmer=0; kmer<=max; kmer++){
			remap[kmer]=kmerMap.get(kmer);
//			System.err.println(AminoAcid.kmerToString(kmer, k2)+" -> "+AminoAcid.kmerToString(remap[kmer], k));
		}
		return remap;
	}
	
	/**
	 * Removes gap bases from the middle of a k-mer to create an ungapped k-mer.
	 * Used for gapped k-mer analysis where middle bases are ignored.
	 *
	 * @param kmer The gapped k-mer
	 * @param k Total k-mer length including gap
	 * @param gap Number of middle bases to remove
	 * @return Ungapped k-mer with gap bases removed
	 */
	public static int ungap(int kmer, int k, int gap) {
		if(gap<1) {return kmer;}
		int half=k/2;
		int halfbits=half*2;
		int gapbits=2*gap;
		int mask=~((-1)<<halfbits);
		int ungapped=(kmer&mask)|((kmer>>gapbits)&~mask);
		return ungapped;
	}
	
	/**
	 * Creates a mapping from canonical k-mer indices to their GC content.
	 * Counts G and C bases in each canonical k-mer.
	 *
	 * @param k K-mer length
	 * @param remap Array mapping k-mers to canonical indices
	 * @return Array where index is canonical k-mer index and value is GC count
	 */
	public static int[] gcmap(int k, int[] remap){
		int[] gcContent=new int[] {0, 1, 1, 0};
		final int bits=2*k;
		final int max=(int)((1L<<bits)-1);
		int[] gcmap=new int[canonicalKmers[k]];
		for(int kmer=0; kmer<=max; kmer++){
			int gc=0;
			for(int i=0, kmer2=kmer; i<k; i++) {
				gc+=gcContent[kmer2&3];
				kmer2>>=2;
			}
			int idx=remap[kmer];
			gcmap[idx]=gc;
		}
		return gcmap;
	}
	
	/**
	 * Counts k-mer frequencies in a sequence and updates the counts array.
	 * Uses canonical k-mer mapping to avoid counting reverse complements separately.
	 *
	 * @param bases Sequence bases to analyze
	 * @param counts Array to store k-mer frequencies (indexed by canonical k-mer)
	 * @param k K-mer length
	 * @return Number of valid k-mers processed
	 */
	public static int countKmers(final byte[] bases, final int[] counts, int k){
		if(quant>1) {return countKmers_quantized(bases, counts, k);}
		if(bases==null || bases.length<k){return 0;}
		
		final int shift=2*k;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
//		int rkmer=0;
		int len=0;
		int valid=0;
		int[] remap=remapMatrix[k];
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
//			int x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
//			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x>=0){
				len++;
				if(len>=k) {
					valid++;
					counts[remap[kmer]]++;
				}
			}else{len=kmer=0;}
		}
		return valid;
	}
	
	/**
	 * Counts k-mers of multiple lengths simultaneously in a single pass.
	 * More efficient than multiple separate counting operations.
	 *
	 * @param bases Sequence bases to analyze
	 * @param counts 2D array where counts[k] stores frequencies for k-mers of length k
	 * @param kmax Maximum k-mer length to count
	 */
	public static void countKmersMulti(final byte[] bases, final long[][] counts, int kmax){
		if(bases==null || bases.length<1){return;}
		
		int kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumberACGTother[b];
			counts[1][x]++;//monomers user 0,1,2,3 plus 4 for undefined
			kmer=((kmer<<2)|x);
			if(x>=0){
				len++;
				for(int k=2; k<=kmax && k<=len; k++) {
					int masked=kmer&masks[k];
					int canon=remapMatrix[k][masked];
					counts[k][canon]++;
				}
			}else{len=kmer=0;}
		}
	}
	
	/**
	 * Counts k-mer frequencies with quantization to reduce memory usage.
	 * Only counts k-mers that are divisible by the quantization factor.
	 *
	 * @param bases Sequence bases to analyze
	 * @param counts Array to store k-mer frequencies (indexed by canonical k-mer)
	 * @param k K-mer length
	 * @return Number of valid k-mers processed
	 */
	public static int countKmers_quantized(final byte[] bases, final int[] counts, int k){
		if(bases==null || bases.length<k){return 0;}
//		counts=(counts!=null ? counts : new int[canonicalKmers]);
		
		final int shift=2*k;
		final int mask=~((-1)<<shift);
		int[] remap=remapMatrix[k];
		
		int kmer=0;
		int len=0;
		int valid=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k) {
					valid++;
					int pos=remap[kmer];
					if(pos>=0) {counts[pos]++;}
				}
			}else{len=kmer=0;}
		}
		return valid;
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}

		return (float)sum;
	}
	
	/**
	 * Calculates normalized absolute difference between two k-mer count arrays.
	 * Automatically normalizes by the sum of each array.
	 *
	 * @param a First k-mer count array
	 * @param b Second k-mer count array
	 * @return Normalized sum of absolute differences
	 */
	static final float absDif(int[] a, int[] b){
		return absDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * @param a Contig kmer counts
	 * @param b Cluster kmer counts
	 * @return Score
	 */
	static final float absDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float rmsDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
//			double d=Tools.absdif((double)a[i], (double)b[i]);
			double d=(a[i])-(b[i]);
			sum+=d*d;
		}

		return (float)Math.sqrt(sum/a.length);
	}
	
	/**
	 * Calculates normalized root mean square difference between two k-mer count arrays.
	 * Automatically normalizes by the sum of each array.
	 *
	 * @param a First k-mer count array
	 * @param b Second k-mer count array
	 * @return Normalized root mean square of differences
	 */
	static final float rmsDif(int[] a, int[] b){
		return rmsDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * @param a Contig kmer counts
	 * @param b Cluster kmer counts
	 * @return Score
	 */
	static final float rmsDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		long sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			float d=(ai-bi);
			sum+=d*d;
		}
		return (float)Math.sqrt(sum/a.length);
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float ksFunction(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			double ai=a[i]+0.0005;
			double bi=b[i]+0.0005;
			double d=(double)ai*Math.log(ai/bi);
			sum+=d;
		}
		
		return (float)sum;
	}
	
	/**
	 * Validates a collection of bins for consistency and correctness.
	 * Checks that clusters contain valid contigs and contigs have proper cluster relationships.
	 *
	 * @param list Collection of bins to validate
	 * @param allowLeafContigs Whether to allow contigs without cluster assignments
	 * @return true if all bins are valid
	 */
	static boolean isValid(Collection<? extends Bin> list, boolean allowLeafContigs) {
		for(Bin b : list) {
			if(b.isCluster()) {
				Cluster c=(Cluster)b;
				assert(c.isValid());
				for(Contig x : c.contigs) {assert(x.isValid());}
			}else {
				Contig c=(Contig)b;
				assert(c.isValid());
				assert(allowLeafContigs || c.cluster()==null);
//				assert(c.cluster()==null || c.cluster().isValid()); //This is too slow
			}
		}
		return true;
	}
	
//	public static float calculateShannonEntropy(float[] depths) {
//		IntHashMap map=new IntHashMap(depths.length*2);
//		for (int i=0; i<depths.length; i++) {
//			float depth=depths[i]+0.125f;
//			int log=(int)Math.round(4*Math.log(depth));
//			map.increment(log);
//		}
//		final float invTotal=1f/depths.length;
//		final int[] values=map.values();
//		float entropy=0;
//		for(int count : values) {
//			if(count>0) {
//				float probability=(float)count*invTotal;
//				entropy-=probability*Math.log(probability)*Tools.invlog2;
//			}
//		}
//		return entropy;
//	}
	
	//Not useful
	/**
	 * Calculates Shannon entropy of depth distribution for diversity measurement.
	 * Uses logarithmic binning of depth values to reduce noise.
	 * Higher entropy indicates more diverse depth distribution.
	 *
	 * @param depths List of depth values
	 * @param limit Maximum number of depths to consider
	 * @return Shannon entropy of the depth distribution
	 */
	static float calculateShannonEntropy(FloatList depths, int limit) {
		final int numDepths=Math.min(limit, depths.size);
		IntHashMap map=new IntHashMap(numDepths*2);
		for(int i=0; i<numDepths; i++) {
			float depth=(float)(depths.get(i)*invSampleDepthSum[i]+0.25f);
			int log=(int)Math.round(8*Math.log(depth));
			map.increment(log);
		}
		final float invTotal=1f/numDepths;
		final int[] values=map.values();
		float entropy=0;
		for(int count : values) {
			if(count>0) {
				float probability=(float)count*invTotal;
				entropy-=probability*Math.log(probability)*Tools.invlog2;
			}
		}
//		System.err.println("depths="+depths);
//		System.err.println("values="+Arrays.toString(values));
//		System.err.println("entropy="+entropy);
//		assert(Math.random()>0.1);
		return entropy;
	}
	
	/**
	 * Counts the number of distinct depth bins in the distribution.
	 * Uses logarithmic binning to group similar depth values.
	 * @param depths List of depth values
	 * @return Number of distinct depth bins
	 */
	static int calculateDistinctValues(FloatList depths) {
		final int numDepths=depths.size;
		IntHashSet set=new IntHashSet(numDepths);//Could alternately use an IntList and sort it
		for(int i=0; i<numDepths; i++) {
			float depth=(float)(depths.get(i)*invSampleDepthSum[i]+0.25f);
			int log=(int)Math.round(8*Math.log(depth));
			set.add(log);
		}
		return set.size();
	}
	
	/**
	 * Loads the taxonomic tree for taxonomic classification.
	 * Uses default tree file if treePath is set to "auto".
	 * @return Loaded taxonomic tree or null if loading fails
	 */
	static TaxTree loadTree() {
		if("auto".equals(treePath)){treePath=TaxTree.defaultTreeFile();}
		if(treePath!=null) {
			tree=CladeObject.tree=TaxTree.loadTaxTree(treePath, System.err, false, false);
		}
		return tree;
	}
	
	/** Looks for tid_1234 or tid|1234, with any delimiters */
	public static int parseTaxID(String line) {
		if(!parseTaxid) {return -1;}
		String term="tid_";
		int pos=line.indexOf(term);
		if(pos<0) {pos=line.indexOf("tid|");}
		if(pos<0) {return -1;}
		long id=0;
		for(int i=pos+4; i<line.length(); i++) {
			char c=line.charAt(i);
			if(c<'0' || c>'9') {break;}
			id=id*10+(c-'0');
		}
		assert(id>0 && id<Integer.MAX_VALUE) : id+"\n"+line+"\n";
		return (int)id;
	}
	
	/**
	 * Resolves a taxonomic ID using the loaded taxonomic tree.
	 * Converts sequence header to taxonomic ID and resolves to canonical form.
	 * @param s Sequence header string
	 * @return Resolved taxonomic ID or original parsed ID if tree unavailable
	 */
	public static int resolveTaxID(String s) {
		int tid=parseTaxID(s);
		if(tid<1 || tree==null) {return tid;}
		return tree.resolveID(tid);
	}
	
//	public static int k() {return k;}
////	public static int gap() {return gap;}
//
//	/** Kmer length for frequencies */
//	private static int k=4;
//	/** Kmer gap length */
////	private static int gap=0;
	
	/** Quantization factor determining how many k-mers to use for comparisons */
	private static int quant=1;//Determines how many tetramers to use for comparisons
	/** Maps a kmer to index in frequency array */
	public static int[][] remapMatrix=makeRemapMatrix(2, 5, true);
	/** Number of canonical kmers; frequency array length */
	public static int[] canonicalKmers=makeCanonicalKmers();
	/** Inverse of canonical k-mer counts for normalization efficiency */
	public static float[] invCanonicalKmers=makeInvCanonicalKmers();
	/** Maps a kmer to index in gc content array */
	public static int[][] gcmapMatrix=makeGCMapMatrix();
	
	/** Bit masks for extracting k-mers of different lengths */
	private static final int[] masks={0, 3, 15, 63, 255, 1023, 4095};
	
	/** Print status messages to this output stream */
	static PrintStream outstream=System.err;
	/** Taxonomic tree for sequence classification */
	static TaxTree tree=null;
	/** Path to taxonomic tree file */
	static String treePath="auto";
	
	/** Minimum size threshold for forming clusters */
	static int minClusterSize=50000;
	/** Minimum number of contigs required per cluster */
	static int minContigsPerCluster=1;
	/** Boost factor for depth-based scoring */
	static float depthBoost=0.25f;
	/** Method number for calculating depth ratios */
	static int depthRatioMethod=4;
	
	/** Whether to include Euclidean distance in comparison calculations */
	static boolean addEuclidian=false;
	/** Whether to include Hellinger distance in comparison calculations */
	static boolean addHellinger=true;
	/** Whether to include 3-mer Hellinger distance in comparisons */
	static boolean addHellinger3=true;
	/** Whether to include 5-mer Hellinger distance in comparisons */
	static boolean addHellinger5=true;
	/** Whether to include absolute difference in comparison calculations */
	static boolean addAbsDif=true;
	/** Whether to include Jensen-Shannon divergence in comparisons */
	static boolean addJsDiv=true;
	/** Whether to include entropy measures in comparison calculations */
	static boolean addEntropy=true;
	/** Whether to include strandedness analysis in comparisons */
	static boolean addStrandedness=true;
	/** Whether to include GC composition in comparison calculations */
	static boolean addGCComp=true;
	/** Multiplier for small number handling in vector calculations */
	static float vectorSmallNumberMult=5f;
	/** Whether to apply square root transformation to small numbers */
	static boolean vectorSmallNumberRoot=false;
	/** Flag indicating if currently creating bin mappings */
	static boolean makingBinMap=false;
	
	/** Whether to count 3-mer frequencies */
	public static boolean countTrimers=true;
	/** Whether to count 5-mer frequencies */
	public static boolean countPentamers=true;
	/** Minimum sequence size for counting pentamers */
	public static int minPentamerSizeCount=2000;
	/** Minimum sequence size for pentamer-based comparisons */
	public static int minPentamerSizeCompare=40000;
	
	/** Whether to print verbose output messages */
	static boolean loud=false;
	/** Whether to print detailed processing information */
	static boolean verbose;
	/** Whether to print stepwise correlation coefficient information */
	static boolean printStepwiseCC=false;
	
	public static float sketchDensity=1/100f;
	/** Whether to create sketches for individual contigs */
	static boolean sketchContigs=false;
	/** Whether to create sketches for clusters */
	static boolean sketchClusters=false;
	/** Whether to output sketch information */
	static boolean sketchOutput=false;
	public static boolean sketchInBulk=true;
	/** Size parameter for sketch generation */
	static int sketchSize=20000;
	
	/** Whether validation mode is enabled */
	static boolean validation=false;
	/** Whether grading mode is enabled */
	static boolean grading=false;
	/** Whether to parse taxonomic IDs from sequence headers */
	static boolean parseTaxid=true;
	/** Whether to use proxy values for zero depth calculations */
	static boolean depthZeroProxy=true;
	/** Global time counter for processing operations */
	static int globalTime=0;

	/** Sum of depths for each sample in multi-sample analysis */
	static double[] sampleDepthSum;
	/** Inverse of sample depth sums for normalization */
	static double[] invSampleDepthSum;
	/** Entropy measure across samples */
	static double sampleEntropy=1;
	/** Number of equivalent samples for statistical calculations */
	static int samplesEquivalent=1;
	/** Neural network for small sequence classification */
	static CellNet net0small=null;
	/** Neural network for medium sequence classification */
	static CellNet net0mid=null;
	/** Neural network for large sequence classification */
	static CellNet net0large=null;

	/** K-mer length for entropy calculations */
	static int entropyK=4;
	/** Window size for entropy calculations */
	static int entropyWindow=150;
	/** Whether to calculate clade-based entropy */
	static boolean calcCladeEntropy=false;//Currently this just affects queries, not ref.
	/** Minimum lineage level for entropy calculations */
	static int MIN_LINEAGE_LEVEL_E=0;
	
}
