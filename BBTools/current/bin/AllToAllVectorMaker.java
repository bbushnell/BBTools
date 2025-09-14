package bin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

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
	
	private FloatList makeVector(ArrayList<Contig> alist, ArrayList<Contig> blist, 
			int minSize, int maxSize, Oracle oracle) {
		int numClusters=randy.nextInt(3);
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
	
	private void decluster(Cluster clust) {
		for(Contig c : clust) {c.cluster=null; c.dest=-1;}
		clust.clear();
	}
	
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
	
	private Bin selectBin(ArrayList<Contig> list, int maxElements, int minSize, int maxSize, IntHashSet used) {
		if(maxElements==1) {return selectContig(list, minSize, maxSize, used);}
		return selectCluster(list, maxElements, minSize, maxSize, used, 1);
	}
	
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
	
	private int randomIndex(Random randy, int max, int rolls) {
		int idx=randy.nextInt(max);
		for(int i=randy.nextInt(rolls+1); i>0; i--) {
			idx=Math.min(idx, randy.nextInt(max));
		}
		return idx;
	}
	
	private static long toKey(int a, int b) {
		return (((long)Math.min(a, b))<<32)|(long)Math.max(a, b);
	}
	
	/*--------------------------------------------------------------*/
	
	FloatList getDifList(int size, int sameGenome) {
		int idx=KmerProb.quantizeLength(size);
		FloatList[] matrix=kmerDifMatrix[sameGenome];
		if(matrix[idx]==null) {matrix[idx]=new FloatList();}
		return matrix[idx];
	}
	
	FloatList[][] kmerDifMatrix=new FloatList[2][38];
	
	/*--------------------------------------------------------------*/
	
	private String out1=null;

	private String outKmerDif=null;
	private String outKmerFraction=null;
	
	private final FileFormat ffout1;
	
	DataLoader loader=null;
	long seed=-1;
	long lines=1000000;
	long linesOut=0;
	long posCount=0;
	long negCount=0;
	float positiveRate=0.5f;
	float edgeFraction=0.1f;
	int baseRolls=2;
	long positiveLines=0;
	long negativeLines=0;
	int maxClusterContigs=9;
	Random randy;

	double dif3good=0;
	double dif34good=0;
	long count3good=0;
	
	double dif45good=0;
	double dif5good=0;
	long count5good=0;

	float maxGCDif=1.0f;//0.15
	float maxDepthRatio=1000.0f;//2.4
	float maxKmerDif=1.0f;//0.02
	final float maxProduct;
	
	int minlen=100;
	int maxlen=2000000000;
	
	ArrayList<Contig> allContigs=null;
	ArrayList<ArrayList<Contig>> allSets=null;
	private final ByteBuilder lineBuffer=new ByteBuilder();
	private final FloatList vecBuffer=new FloatList();
	
	/*--------------------------------------------------------------*/
	
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
