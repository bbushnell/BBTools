package bin;

import java.util.ArrayList;
import java.util.HashMap;
import shared.Random;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntHashSet;
import structures.LongHashSet;

/**
 * @author Brian Bushnell
 * @date Feb 23, 2025
 *
 */
public class AllToAllVectorMaker extends BinObject {

	/** Program entry point for generating all-to-all comparison vectors.
	 * @param args Command-line arguments for configuration */
	public static void main(String[] args){
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AllToAllVectorMaker x=new AllToAllVectorMaker(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs an AllToAllVectorMaker with command-line arguments.
	 * Parses configuration parameters and initializes data loader.
	 * @param args Command-line arguments containing configuration options
	 */
	public AllToAllVectorMaker(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		loader=new DataLoader(outstream);
		loader.netFileLarge=loader.netFileMid=loader.netFileSmall=null;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("rate") || a.equals("positivity")){
				positiveRate=Float.parseFloat(b);
			}else if(a.equals("edgefraction")){
				edgeFraction=Float.parseFloat(b);
			}else if(a.equals("gcdif") || a.equals("maxgcdif")){
				maxGCDif=Float.parseFloat(b);
			}else if(a.equals("maxkmerdif")){
				maxKmerDif=Float.parseFloat(b);
			}else if(a.equals("maxdepthratio")){
				maxDepthRatio=Float.parseFloat(b);
			}else if(a.equals("lines")){
				lines=Parse.parseKMG(b);
			}else if(a.equals("rolls")){
				baseRolls=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxClusterContigs") || a.equalsIgnoreCase("mcc")){
				maxClusterContigs=Integer.parseInt(b);
			}else if(a.equals("kmerdif") || a.equals("outkmerdif")){
				outKmerDif=b;
			}else if(a.equals("kmerfraction") || a.equals("outkmerfraction")){
				outKmerFraction=b;
			}else if(a.equals("minlen")){
				minlen=Parse.parseIntKMG(b);
			}else if(a.equals("maxlen")){
				maxlen=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("printSizeInVector")){
				Oracle.printSizeInVector=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printNetOutputInVector")){
				Oracle.printNetOutputInVector=Parse.parseBoolean(b);
			}else if(loader.parse(arg, a, b)){
				//do nothing
			}else if(SimilarityMeasures.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			out1=parser.out1;
		}
		maxProduct=maxKmerDif*maxDepthRatio*0.75f;
		KmerProb.load();
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
	}
	
	/**
	 * Main processing method that generates comparison vectors.
	 * Loads contigs, groups by taxonomy, and outputs training data.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		validation=true;
		ArrayList<Contig> contigs=allContigs=loader.loadData();
		if(verbose){outstream.println("Finished reading data.");}
		loader.trimEdges(contigs, Binner.maxEdges, Binner.minEdgeWeight, true);
		
		HashMap<Integer, ArrayList<Contig>> map=new HashMap<Integer, ArrayList<Contig>>();
		for(Contig c : contigs) {
			if(c.labelTaxid>0) {
				ArrayList<Contig> list=map.get(c.labelTaxid);
				if(list==null) {map.put(c.labelTaxid, list=new ArrayList<Contig>());}
				list.add(c);
			}
		}
		
		outputResults(contigs, map);
		if(outKmerDif!=null) {
			outputKmerDifs(outKmerDif, 0);
			outputKmerDifs(outKmerDif, 1);
		}
		if(outKmerFraction!=null) {
			outputKmerDifFraction(outKmerFraction, 1.0/1024, 0.25);
		}
		
		t.stop();
		
		outstream.println("dif3 avg: \t"+dif3good/count3good);
		outstream.println("dif34 avg:\t"+dif34good/count3good);
		outstream.println("dif45 avg:\t"+dif45good/count5good);
		outstream.println("dif5 avg: \t"+dif5good/count5good);
		outstream.println("dif3/dif4:\t"+dif3good/dif34good);
		outstream.println("dif4/dif5:\t"+dif45good/dif5good);

		outstream.println("Positive: \t"+positiveLines);
		outstream.println("Negative: \t"+negativeLines);
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+loader.contigsLoaded+
				" \t"+Tools.format("%.2fk bases/sec", (loader.basesLoaded/(double)(t.elapsed))*1000000));
		assert(!errorState) : "An error was encountered.";
	}
	
	/**
	 * Generates and outputs comparison vectors for machine learning training.
	 * Creates positive and negative pairs according to specified rate.
	 * @param contigs List of all contigs to compare
	 * @param map Contigs grouped by taxonomic ID
	 */
	private void outputResults(ArrayList<Contig> contigs, HashMap<Integer, ArrayList<Contig>> map){
		LongHashSet used=new LongHashSet();
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
		
		randy=Shared.threadLocalRandom(seed);
		Oracle oracle=new Oracle(999999, 0);
		if(bsw!=null) {//Print header
			vecBuffer.clear();
			oracle.toVector(contigs.get(0), contigs.get(1), vecBuffer, true);
			int weights=Oracle.printWeightInVector>0 ? 1 : 0;
			bsw.print("#dims\t").print(vecBuffer.size-(1+weights)).tab().print(1).tab().println(weights);
		}
		ArrayList<ArrayList<Contig>> clusters=new ArrayList<ArrayList<Contig>>(map.values());
		
		while(linesOut<lines) {
			final boolean positive=(randy.nextFloat()<=positiveRate);
			ByteBuilder bb=null;
			while(bb==null) {bb=makeLine(clusters, oracle, positive);}
			if(bb!=null) {
				if(bsw!=null) {
					bsw.print(bb.nl());
					bb.clear();
				}
				linesOut++;
			}
		}

		if(bsw!=null) {
			errorState=bsw.poisonAndWait() | errorState;
		}
	}
	
	/**
	 * Outputs k-mer difference distributions to file.
	 * Generates percentile distributions for different contig lengths.
	 * @param fname Output filename (% replaced with sign value)
	 * @param sign Classification sign (0=different taxa, 1=same taxa)
	 */
	private void outputKmerDifs(String fname, int sign) {
		fname=fname.replaceFirst("%", sign+"");
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, true, false, true);
		FloatList[] matrix=kmerDifMatrix[sign];
		assert(matrix!=null);
		for(int lenIdx=0; lenIdx<matrix.length; lenIdx++) {
			FloatList list=matrix[lenIdx];
//			System.err.println(lenIdx+": "+(list==null ? "null" : list.size()+""));
			if(list!=null && list.size()>=100) {
				int length=KmerProb.dequantizeLength(lenIdx);
				list.sort();
				bsw.print(length).tab().print(list.size());
				for(int i=0; i<=100; i++) {
					bsw.tab().print(list.percentile(i*0.01), 8);
				}
				bsw.println();
			}
		}
		
		int x=Tools.binarySearch(new float[1], 1);
		//TODO: use this to make an array of percentiles indexed by kmer dif.
		
		bsw.poisonAndWait();
	}
	
	/**
	 * Outputs probability fractions indexed by k-mer difference values.
	 * Creates lookup table for classification probability by k-mer similarity.
	 *
	 * @param fname Output filename
	 * @param incr Increment size for k-mer difference bins
	 * @param max Maximum k-mer difference value
	 */
	private void outputKmerDifFraction(String fname, double incr, double max) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, true, false, true);
		bsw.println("#ceil(size)\tcount\tprobs");
		FloatList[] matrix1=kmerDifMatrix[1];
		FloatList[] matrix0=kmerDifMatrix[0];
		
		bsw.print(" \t ");
		float key=0;
		for(int i=0; key<=max; i++) {
			key=(float)(i*incr);
			bsw.tab().print(key, 10);
		}
		
		for(int lenIdx=0; lenIdx<matrix1.length; lenIdx++) {
			FloatList list1=matrix1[lenIdx];
			FloatList list0=matrix0[lenIdx];
			int size=Tools.min(list1==null ? 0 : list1.size(), list0==null ? 0 : list0.size());
			if(size>100) {
				int length=KmerProb.dequantizeLength(lenIdx);
				FloatList fractions=fractionIndexedByKmerDif(list1, list0, incr, max);
				bsw.print(length).tab().print(size);
				for(int i=0; i<fractions.size(); i++) {
					bsw.tab().print(fractions.get(i), 8);
				}
				bsw.println();
//				bsw.poisonAndWait();
//				assert(false);
			}
		}
		
		int x=Tools.binarySearch(new float[1], 1);
		//TODO: use this to make an array of percentiles indexed by kmer dif.
		
		bsw.poisonAndWait();
	}
	
	/**
	 * Calculates probability fractions from positive and negative k-mer difference lists.
	 * Creates smoothed probability curves for classification decisions.
	 *
	 * @param plus K-mer differences for same-taxa pairs
	 * @param minus K-mer differences for different-taxa pairs
	 * @param incr Bin size for difference values
	 * @param max Maximum difference to analyze
	 * @return List of probability fractions indexed by k-mer difference
	 */
	private FloatList fractionIndexedByKmerDif(FloatList plus, FloatList minus, double incr, double max) {
		plus.shrink().sort();
		minus.shrink().sort();
		float invPlus=1f/Math.max(1, plus.size());
		float invMinus=1f/Math.max(1, minus.size());
		FloatList fractions=new FloatList(1+(int)Math.ceil(max/incr));
		float key=0;
		for(int i=0; key<=max; i++) {
			key=(float)(i*incr);
			int idxPlus=Tools.max(1, Tools.binarySearch(plus.array, key));
			int idxMinus=Tools.binarySearch(minus.array, key);
			float fractionPlus=idxPlus*invPlus;
			float fractionMinus=idxMinus*invMinus;
			float fraction=fractionPlus/(fractionPlus+fractionMinus);
			fractions.add(fraction);
		}
		for(int i=fractions.size()-2; i>=0; i--) {//Fix low sample size weirdness
			fractions.set(i, Tools.max(fractions.get(i), fractions.get(i+1)));
		}
		return fractions;
	}
	
//	private ByteBuilder makeLine(ArrayList<Contig> contigs, HashMap<Integer, ArrayList<Contig>> map, 
//			LongHashSet used, Oracle oracle) {
//		Contig a=null;
//		while(a==null || a.labelTaxid<1 || a.size()<minlen) {
//			int idx=randomIndex(randy, contigs.size(), baseRolls+1);
//			a=contigs.get(idx);
//		}
//		assert(a.labelTaxid>0);
//		boolean positive=(randy.nextFloat()<=positiveRate);
//		ArrayList<Contig> list=(positive ? map.get(a.labelTaxid) : contigs);
//		Contig b=findOther(a, list, used, null, randy, positive);
//		if(b==null) {return null;}
//		assert(b.labelTaxid>0) : a.name()+", "+b.name()+", "+positive;
//		vecBuffer.clear();
//		oracle.toVector(a, b, vecBuffer, true);
//		if(outKmerDif!=null || outKmerFraction!=null) {
//			int same=(a.labelTaxid==b.labelTaxid) ? 1 : 0;
//			float dif=SimilarityMeasures.calculateDifferenceAverage(a.counts, b.counts);
//			int size=(int)Tools.min(a.size(), b.size());
//			FloatList difs=getDifList(size, same);
//			difs.add(dif);
//			assert(getDifList(size, same).size>0);
//		}
////		assert(false) : Arrays.toString(kmerDifMatrix)+", "+
////			Arrays.toString(kmerDifMatrix[0])+", "+Arrays.toString(kmerDifMatrix[1]);
//		
//		return toLine(vecBuffer);
//	}
	
	/**
	 * Creates a single training line by comparing contigs or clusters.
	 * Randomly selects comparison type and generates feature vector.
	 *
	 * @param clusters Contig clusters grouped by taxonomy
	 * @param oracle Feature extraction engine
	 * @param positive Whether to create positive (same taxa) or negative pair
	 * @return Formatted training line or null if no valid pair found
	 */
	private ByteBuilder makeLine(ArrayList<ArrayList<Contig>> clusters, Oracle oracle, final boolean positive) {
		ArrayList<Contig> alist=clusters.get(randy.nextInt(clusters.size()));
		ArrayList<Contig> blist=alist;
		while(!positive && alist==blist) {blist=clusters.get(randy.nextInt(clusters.size()));}
		FloatList vector=null;
		for(int i=0; i<9 && vector==null; i++) {
			vector=makeVector(alist, blist, minlen, maxlen, oracle);
		}
//		System.err.println(vector==null ? "fail" : "success");
		return vector==null ? null : toLine(vector);
	}
	
	/**
	 * Creates feature vector by comparing elements from two contig lists.
	 * Randomly selects individual contigs or clusters for comparison.
	 *
	 * @param alist First contig list
	 * @param blist Second contig list
	 * @param minSize Minimum element size in bases
	 * @param maxSize Maximum element size in bases
	 * @param oracle Feature extraction engine
	 * @return Feature vector or null if no valid comparison possible
	 */
	private FloatList makeVector(ArrayList<Contig> alist, ArrayList<Contig> blist, 
			int minSize, int maxSize, Oracle oracle) {
		int numClusters=randy.nextInt3();
//		System.err.println(numClusters+", "+(alist==blist));
//		System.err.println("numClusters="+numClusters);
		if(numClusters==0) {
			IntHashSet used=new IntHashSet(7);
			Contig a=selectContig(alist, minSize, maxSize, used);
			Contig b=selectContig(blist, minSize, Integer.MAX_VALUE, used);
			if(!passesFilter(a, b)) {return null;}
			assert(a!=b);
			vecBuffer.clear();
//			System.err.println("size="+a.size()+", "+a.numContigs()+", "+b.size()+", "+b.numContigs());
			trackRatio(a, b);
			return oracle.toVector(a, b, vecBuffer, true);
		}else if(numClusters==1) {
			IntHashSet used=new IntHashSet(7);
			Contig a=selectContig(alist, minSize, maxSize, used);
			if(a==null) {return null;}
			Cluster b=selectCluster(blist, 2+randy.nextInt(maxClusterContigs-1), minSize, Integer.MAX_VALUE, used, 3);
			if(!passesFilter(a, b)) {return null;}
			vecBuffer.clear();
//			System.err.println("size="+a.size()+", "+a.numContigs()+", "+b.size()+", "+b.numContigs());
			trackRatio(a, b);
			FloatList fl=oracle.toVector(a, b, vecBuffer, true);
			decluster(b);
			return fl;
		}else {
			Cluster a=selectCluster(alist, 2+randy.nextInt(maxClusterContigs-1), minSize, maxSize, null, 3);
			if(a==null) {return null;}
			Cluster b=selectCluster(blist, 2+randy.nextInt(maxClusterContigs-1), minSize, Integer.MAX_VALUE, a.contigSet, 3);
			if(!passesFilter(a, b)) {
				decluster(a);
				return null;
			}
			vecBuffer.clear();
//			System.err.println("size="+a.size()+", "+a.numContigs()+", "+b.size()+", "+b.numContigs());
			trackRatio(a, b);
			FloatList fl=oracle.toVector(a, b, vecBuffer, true);
			decluster(a);
			decluster(b);
			return fl;
		}
	}
	
	/**
	 * Resets cluster state by clearing contig assignments.
	 * Allows contigs to be reassigned to different clusters.
	 * @param clust Cluster to reset
	 */
	private void decluster(Cluster clust) {
		for(Contig c : clust) {c.cluster=null; c.dest=-1;}
		clust.clear();
	}
	
	/**
	 * Tracks k-mer difference statistics for same-taxon pairs.
	 * Accumulates cosine differences for 3-mer, 4-mer, and 5-mer counts.
	 * @param a First bin for comparison
	 * @param b Second bin for comparison
	 */
	private void trackRatio(Bin a, Bin b) {
		if(a.labelTaxid!=b.labelTaxid) {return;}
		double dif3=SimilarityMeasures.cosineDifference(a.trimers, b.trimers);
		double dif4=SimilarityMeasures.cosineDifference(a.tetramers, b.tetramers);
		double dif5=(a.numPentamers<BinObject.minPentamerSizeCompare ||
				b.numPentamers<BinObject.minPentamerSizeCompare ? -1 :
					SimilarityMeasures.cosineDifference(a.pentamers, b.pentamers));
		
		dif3good+=dif3;
		dif34good+=dif4;
		count3good++;
		
		if(dif5<0) {return;}
		dif45good+=dif4;
		dif5good+=dif5;
		count5good++;
	}
	
	/**
	 * Selects either a single contig or cluster from list.
	 * Wrapper method that delegates to appropriate selection function.
	 *
	 * @param list Contig list to select from
	 * @param maxElements Maximum elements in selection (1=contig, >1=cluster)
	 * @param minSize Minimum size in bases
	 * @param maxSize Maximum size in bases
	 * @param used Set of already-used contig IDs
	 * @return Selected bin (contig or cluster)
	 */
	private Bin selectBin(ArrayList<Contig> list, int maxElements, int minSize, int maxSize, IntHashSet used) {
		if(maxElements==1) {return selectContig(list, minSize, maxSize, used);}
		return selectCluster(list, maxElements, minSize, maxSize, used, 1);
	}
	
	/**
	 * Creates cluster by randomly selecting contigs from list.
	 * Attempts to meet size constraints through iterative selection.
	 *
	 * @param list Contig list to select from
	 * @param maxElements Maximum contigs in cluster
	 * @param minSize Minimum total cluster size in bases
	 * @param maxSize Maximum total cluster size in bases
	 * @param used Set of contig IDs to avoid
	 * @param tries Number of retry attempts with larger clusters
	 * @return Created cluster or null if constraints cannot be met
	 */
	private Cluster selectCluster(ArrayList<Contig> list, int maxElements, int minSize, int maxSize, IntHashSet used, int tries) {
		IntHashSet set=new IntHashSet(7);
		long size=0;
		for(int i=0; i<100; i++) {
			Contig c=list.get(randy.nextInt(list.size()));
			long size2=size+c.size();
			if(size2>=minSize && size2<=maxSize && !set.contains(c.id()) 
					&& (used==null || !used.contains(c.id()))) {
				set.add(c.id());
				size=size2;
			}
			if(set.size()>=maxElements) {break;}
			if(size>minSize) {
				if(i>20 && set.size()>=2) {break;}
				if(randy.nextFloat()<0.05f) {break;}
			}
		}
		if(size<minSize || size>maxSize) {//fail
			if(size<minSize && tries>1 && set.size()>=maxElements) {
				return selectCluster(list, maxElements+maxClusterContigs, minSize, maxSize, used, tries-1);
			}
			return null;
		}
		Cluster clust=new Cluster(0);
		for(int i : set.toArray()) {
			clust.add(allContigs.get(i));
		}
		return clust;
	}
	
	/**
	 * Randomly selects single contig meeting size constraints.
	 * Tries multiple attempts to find valid contig.
	 *
	 * @param list Contig list to select from
	 * @param minSize Minimum contig size in bases
	 * @param maxSize Maximum contig size in bases
	 * @param used Set of contig IDs to avoid
	 * @return Selected contig or null if no valid contig found
	 */
	private Contig selectContig(ArrayList<Contig> list, int minSize, int maxSize, IntHashSet used) {
		for(int i=0; i<40; i++) {
			Contig c=list.get(randy.nextInt(list.size()));
			if(c.size()>=minSize && c.size()<=maxSize && (used==null || !used.contains(c.id()))) {
				used.add(c.id());
				return c;
			}
		}
//		System.err.println("Can't find contig in range ("+minSize+", "+maxSize+") in list: ");
//		for(int i=0; i<list.size() && i<1000; i++) {
//			System.err.print(list.get(i).size()+", ");
//		}
		return null;
	}
	
	/**
	 * Formats feature vector as tab-separated output line.
	 * Tracks positive and negative line counts.
	 * @param vector Feature values to format
	 * @return Formatted line for output
	 */
	private ByteBuilder toLine(FloatList vector) {
		lineBuffer.clear();
		for(int i=0; i<vector.size(); i++) {
			if(i>0) {lineBuffer.tab();}
			lineBuffer.append(vector.get(i), 7, true);
		}
		if(vector.lastElement()==1) {positiveLines++;}
		else {negativeLines++;}
		return lineBuffer;
	}
	
//	private Contig findOther(final Contig a, ArrayList<Contig> contigs, 
//			LongHashSet used, Oracle oracle, Random randy, boolean positive) {
//		for(int i=0; i<100; i++) {
//			int idx=randomIndex(randy, contigs.size(), baseRolls);
//			Contig b=contigs.get(idx);
//			if(a.pairMap!=null && randy.nextFloat()<edgeFraction) {
//				ArrayList<KeyValue> edges=KeyValue.toList(a.pairMap);
//				KeyValue kv=edges.get(randy.nextInt(Tools.min(edges.size(), 4)));
//				if(kv.key<allContigs.size()) {b=allContigs.get(kv.key);}
//				positive=(a.labelTaxid==b.labelTaxid);//Keep it either way
////				System.err.print('.');
//			}
//			boolean same=(a.labelTaxid==b.labelTaxid);
//			if(a!=b && b.labelTaxid>0 && (same || Math.abs(a.gc()-b.gc())<=maxGCDif) &&
//					(a.size()<=maxlen && b.size()<=maxlen) && b.size()>=minlen) {
//				final long key=toKey(a.id(), b.id());
//				if((a.labelTaxid==b.labelTaxid)==positive && !used.contains(key)) {
//					if(same || oracle==null || oracle.similarity(a, b, 1)>=0) {
//						used.add(key);
//						return b;
//					}
//				}
//			}
//		}
//		return null;
//	}
	
	/**
	 * Determines if bin pair passes quality filters for comparison.
	 * Checks GC content difference, depth ratio, and k-mer difference.
	 *
	 * @param a First bin to compare
	 * @param b Second bin to compare
	 * @return true if pair passes all filters, false otherwise
	 */
	private boolean passesFilter(Bin a, Bin b) {
		if(a==null || b==null || a==b) {return false;}
		final float gcDif=Tools.absdif(a.gc(), b.gc());
		final boolean same=a.labelTaxid==b.labelTaxid;
		if(gcDif>maxGCDif) {
//			System.err.println("Failed filter: "+same+", gcDif="+gcDif);
			return false;
		}
		final float depthRatio=a.depthRatio(b);
		if(depthRatio>maxDepthRatio) {
//			System.err.println("Failed filter: "+same+", depthRatio="+depthRatio);
			return false;
		}
		final float kmerDif=SimilarityMeasures.calculateDifferenceAverage(a.tetramers, b.tetramers);
		if(kmerDif>maxKmerDif) {
//			System.err.println("Failed filter: "+same+", kmerDif="+kmerDif);
			return false;
		}
		final float product=kmerDif*depthRatio;
		if(product>maxProduct) {
//			System.err.println("Failed filter: "+same+", product="+product);
			return false;
		}
		return true;
	}
	
	/**
	 * Generates biased random index favoring smaller values.
	 * Takes minimum of multiple random rolls to bias selection.
	 *
	 * @param randy Random number generator
	 * @param max Maximum index value (exclusive)
	 * @param rolls Number of rolls to take minimum from
	 * @return Biased random index
	 */
	private int randomIndex(Random randy, int max, int rolls) {
		int idx=randy.nextInt(max);
		for(int i=randy.nextInt(rolls+1); i>0; i--) {
			idx=Math.min(idx, randy.nextInt(max));
		}
		return idx;
	}
	
	/**
	 * Creates unique key from two contig IDs.
	 * Ensures consistent key regardless of parameter order.
	 *
	 * @param a First contig ID
	 * @param b Second contig ID
	 * @return Unique 64-bit key for the pair
	 */
	private static long toKey(int a, int b) {
		return (((long)Math.min(a, b))<<32)|(long)Math.max(a, b);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Gets k-mer difference list for given size and classification.
	 * Creates new list if none exists for this size category.
	 *
	 * @param size Contig size to quantize into bins
	 * @param sameGenome Classification (0=different taxa, 1=same taxa)
	 * @return FloatList for accumulating k-mer differences
	 */
	FloatList getDifList(int size, int sameGenome) {
		int idx=KmerProb.quantizeLength(size);
		FloatList[] matrix=kmerDifMatrix[sameGenome];
		if(matrix[idx]==null) {matrix[idx]=new FloatList();}
		return matrix[idx];
	}
	
	/** Matrix storing k-mer difference distributions by size and classification */
	FloatList[][] kmerDifMatrix=new FloatList[2][38];
	
	/*--------------------------------------------------------------*/
	
	/** Primary output filename for training vectors */
	private String out1=null;

	/** Output filename for k-mer difference distributions */
	private String outKmerDif=null;
	/** Output filename for k-mer fraction probability tables */
	private String outKmerFraction=null;
	
	/** File format handler for primary output */
	private final FileFormat ffout1;
	
	/** Data loading component for reading contig data */
	DataLoader loader=null;
	/** Random seed for reproducible pair selection */
	long seed=-1;
	/** Target number of training lines to generate */
	long lines=1000000;
	/** Number of training lines generated so far */
	long linesOut=0;
	long posCount=0;
	long negCount=0;
	/** Fraction of training pairs that should be positive (same taxon) */
	float positiveRate=0.5f;
	/** Fraction of pairs selected from assembly graph edges */
	float edgeFraction=0.1f;
	/** Number of random rolls for biased index selection */
	int baseRolls=2;
	/** Count of positive training examples generated */
	long positiveLines=0;
	/** Count of negative training examples generated */
	long negativeLines=0;
	/** Maximum number of contigs allowed in a cluster */
	int maxClusterContigs=9;
	/** Random number generator for pair selection */
	Random randy;

	/** Accumulated 3-mer cosine differences for same-taxon pairs */
	double dif3good=0;
	/** Accumulated 4-mer cosine differences for same-taxon pairs */
	double dif34good=0;
	/** Count of same-taxon pairs with valid 3/4-mer comparisons */
	long count3good=0;
	
	/** Accumulated 4-mer cosine differences for same-taxon pairs with 5-mers */
	double dif45good=0;
	/** Accumulated 5-mer cosine differences for same-taxon pairs */
	double dif5good=0;
	/** Count of same-taxon pairs with valid 5-mer comparisons */
	long count5good=0;

	/** Maximum GC content difference allowed between compared bins */
	float maxGCDif=1.0f;//0.15
	/** Maximum depth ratio allowed between compared bins */
	float maxDepthRatio=1000.0f;//2.4
	/** Maximum k-mer difference allowed between compared bins */
	float maxKmerDif=1.0f;//0.02
	/** Maximum product of k-mer difference and depth ratio */
	final float maxProduct;
	
	/** Minimum contig length in bases for inclusion */
	int minlen=100;
	/** Maximum contig length in bases for inclusion */
	int maxlen=2000000000;
	
	/** Complete list of loaded contigs */
	ArrayList<Contig> allContigs=null;
	ArrayList<ArrayList<Contig>> allSets=null;
	/** Reusable buffer for formatting output lines */
	private final ByteBuilder lineBuffer=new ByteBuilder();
	/** Reusable buffer for feature vector construction */
	private final FloatList vecBuffer=new FloatList();
	
	/*--------------------------------------------------------------*/
	
	/** Flag indicating whether errors occurred during processing */
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and results */
	private java.io.PrintStream outstream=System.err;
	/** Controls verbosity of status output */
	public static boolean verbose=false;
	
}
