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
 * Superclass for binner classes that provides k-mer frequency analysis and comparison functions.
 * Contains shared functionality for sequence binning operations including k-mer counting,
 * canonical k-mer mapping, GC content analysis, and distance calculations between bins.
 *
 * @author Brian Bushnell
 * @date Feb 4, 2025
 */
public class BinObject {
	
	public static void setQuant(int x) {
		quant=x;
		assert(quant>0);
		initialize();
	}
	
//	public static void setK(int x) {
//		assert(k>0 && k<16);
//		k=x;
//	}
	
	private static void initialize() {
		remapMatrix=makeRemapMatrix(2, 5, true);
		//K=1 is ACGTN for use in GC calcs.
		//K=2 is noncanonical for use in strandedness calcs.
		canonicalKmers=makeCanonicalKmers();
		invCanonicalKmers=makeInvCanonicalKmers();
		gcmapMatrix=makeGCMapMatrix();
	}
	
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
	
	private static synchronized int[] makeCanonicalKmers() {
		int[] array=new int[remapMatrix.length];
		for(int i=0; i<array.length; i++) {
			int[] remap=remapMatrix[i];
			int max=(remap==null ? 1 : Tools.max(remap)+1);
			array[i]=max;
		}
		return array;
	}
	
	private static synchronized float[] makeInvCanonicalKmers() {
		float[] array=new float[canonicalKmers.length];
		for(int i=0; i<array.length; i++) {
			array[i]=1f/canonicalKmers[i];
		}
		return array;
	}
	
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
	
	public static int ungap(int kmer, int k, int gap) {
		if(gap<1) {return kmer;}
		int half=k/2;
		int halfbits=half*2;
		int gapbits=2*gap;
		int mask=~((-1)<<halfbits);
		int ungapped=(kmer&mask)|((kmer>>gapbits)&~mask);
		return ungapped;
	}
	
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
	 * Calculates absolute difference between two k-mer frequency arrays.
	 * Sums the absolute differences of corresponding elements.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return Sum of absolute differences between corresponding elements
	 */
	static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}

		return (float)sum;
	}
	
	static final float absDif(int[] a, int[] b){
		return absDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * Calculates normalized absolute difference between two k-mer count arrays.
	 * Uses pre-computed inverse normalization factors for efficiency.
	 *
	 * @param a First k-mer count array
	 * @param b Second k-mer count array
	 * @param inva Inverse of sum of array a (1/sum(a))
	 * @param invb Inverse of sum of array b (1/sum(b))
	 * @return Normalized sum of absolute differences
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
	 * Calculates root mean square difference between two k-mer frequency arrays.
	 * More sensitive to large differences than absolute difference.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return Root mean square of differences between corresponding elements
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
	
	static final float rmsDif(int[] a, int[] b){
		return rmsDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * Calculates normalized root mean square difference between two k-mer count arrays.
	 * Uses pre-computed inverse normalization factors for efficiency.
	 *
	 * @param a First k-mer count array
	 * @param b Second k-mer count array
	 * @param inva Inverse of sum of array a (1/sum(a))
	 * @param invb Inverse of sum of array b (1/sum(b))
	 * @return Normalized root mean square of differences
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
	 * Calculates Kullback-Leibler divergence between two k-mer frequency distributions.
	 * Measures how one probability distribution diverges from another.
	 * Adds small pseudocount to avoid log(0) issues.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return KL divergence from b to a
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
	
	static TaxTree loadTree() {
		if("auto".equals(treePath)){treePath=TaxTree.defaultTreeFile();}
		if(treePath!=null) {
			tree=CladeObject.tree=TaxTree.loadTaxTree(treePath, System.err, false, false);
		}
		return tree;
	}
	
	/**
	 * Parses taxonomic ID from a sequence header line.
	 * Looks for patterns like "tid_1234" or "tid|1234" with any delimiters.
	 * @param line Sequence header line to parse
	 * @return Taxonomic ID if found, -1 otherwise
	 */
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
	
	private static int quant=1;//Determines how many tetramers to use for comparisons
	/** Maps k-mers to canonical indices for different k values */
	public static int[][] remapMatrix=makeRemapMatrix(2, 5, true);
	/**
	 * Number of canonical k-mers for each k value; determines frequency array length
	 */
	public static int[] canonicalKmers=makeCanonicalKmers();
	public static float[] invCanonicalKmers=makeInvCanonicalKmers();
	/** Maps canonical k-mer indices to their GC content for different k values */
	public static int[][] gcmapMatrix=makeGCMapMatrix();
	
	private static final int[] masks={0, 3, 15, 63, 255, 1023, 4095};
	
	/** Output stream for status messages */
	static PrintStream outstream=System.err;
	static TaxTree tree=null;
	static String treePath="auto";
	
	static int minClusterSize=50000;
	static int minContigsPerCluster=1;
	static float depthBoost=0.25f;
	static int depthRatioMethod=4;
	
	static boolean addEuclidian=false;
	static boolean addHellinger=true;
	static boolean addHellinger3=true;
	static boolean addHellinger5=true;
	static boolean addAbsDif=true;
	static boolean addJsDiv=true;
	static boolean addEntropy=true;
	static boolean addStrandedness=true;
	static boolean addGCComp=true;
	static float vectorSmallNumberMult=5f;
	static boolean vectorSmallNumberRoot=false;
	static boolean makingBinMap=false;
	
	public static boolean countTrimers=true;
	public static boolean countPentamers=true;
	public static int minPentamerSizeCount=2000;
	public static int minPentamerSizeCompare=40000;
	
	static boolean loud=false;
	static boolean verbose;
	static boolean printStepwiseCC=false;
	
	/** Density parameter for sketch-based comparisons */
	public static float sketchDensity=1/100f;
	static boolean sketchContigs=false;
	static boolean sketchClusters=false;
	static boolean sketchOutput=false;
	/** Whether to process sketches in bulk operations */
	public static boolean sketchInBulk=true;
	static int sketchSize=20000;
	
	static boolean validation=false;
	static boolean grading=false;
	static boolean parseTaxid=true;
	static boolean depthZeroProxy=true;
	static int globalTime=0;

	static double[] sampleDepthSum;
	static double[] invSampleDepthSum;
	static double sampleEntropy=1;
	static int samplesEquivalent=1;
	static CellNet net0small=null;
	static CellNet net0mid=null;
	static CellNet net0large=null;

	static int entropyK=4;
	static int entropyWindow=150;
	static boolean calcCladeEntropy=false;//Currently this just affects queries, not ref.
	static int MIN_LINEAGE_LEVEL_E=0;
	
}
