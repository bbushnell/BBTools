package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import bin.SpectraCounter.LoadThread;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntHashMap;
import structures.IntHashSet;
import structures.IntLongHashMap;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

public class Binner extends BinObject implements Accumulator<Binner.CompareThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	Binner(PrintStream outstream_){outstream=outstream_;}

	boolean parse(String arg, String a, String b) {

		if(a.equalsIgnoreCase("productMult")){
			productMult=Float.parseFloat(b);
		}

		else if(a.equalsIgnoreCase("perfectoracle")){
			PERFECT_ORACLE=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("banbadmerges")){
			BAN_BAD_MERGES=Parse.parseBoolean(b);
		}

		else if(a.equalsIgnoreCase("basePasses")){
			basePasses=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("baseRange")){
			baseRange=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("residueRange")){
			residueRange=Integer.parseInt(b);
		}

		else if(a.equalsIgnoreCase("prepass") || a.equalsIgnoreCase("passAA")){
			runPassAA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassA") || a.equalsIgnoreCase("passA")){
			runPassA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassB") || a.equalsIgnoreCase("passB")){
			runPassB=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassC") || a.equalsIgnoreCase("passC")){
			runPassC=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassD") || a.equalsIgnoreCase("passD")){
			runPassD=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassE") || a.equalsIgnoreCase("passE")){
			runPassE=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassF") || a.equalsIgnoreCase("passF")){
			runPassF=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("runPassG") || a.equalsIgnoreCase("passG")){
			runPassG=Parse.parseBoolean(b);
		}

		else if(a.equalsIgnoreCase("overrideSetSamples")){
			overrideSetSamples=Integer.parseInt(b);
		}

		else if(a.equalsIgnoreCase("maxTrimerDif1") || a.equalsIgnoreCase("max3merDif1")){
			max3merDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxDif1") || a.equalsIgnoreCase("maxKmerDif1")){
			max4merDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxPentamerDif1") || a.equalsIgnoreCase("max5merDif1")){
			max5merDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxRatio1") || a.equalsIgnoreCase("maxDepthRatio1")){
			maxDepthRatio1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxGCDif1")){
			maxGCDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxCovariance1")){
			maxCovariance1=Float.parseFloat(b);
		}

		else if(a.equalsIgnoreCase("maxTrimerDif2") || a.equalsIgnoreCase("maxTrimerDif") || 
				a.equalsIgnoreCase("max3merDif")){
			max3merDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxKmerDif2") || a.equalsIgnoreCase("maxKmerDif")  || 
				a.equalsIgnoreCase("max4merDif")){
			max4merDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxPentamerDif2") || a.equalsIgnoreCase("maxPentamerDif") || 
				a.equalsIgnoreCase("max5merDif")){
			max5merDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxDepthRatio2") || a.equalsIgnoreCase("maxDepthRatio")){
			maxDepthRatio2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxGCDif2") || a.equalsIgnoreCase("maxGCDif")){
			maxGCDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxCovariance2") || a.equalsIgnoreCase("maxCovariance")){
			maxCovariance2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("minkmerprob") || a.equalsIgnoreCase("minkmerprob2")){
			minKmerProb2=Float.parseFloat(b);
		}

		else if(a.equalsIgnoreCase("minSimilarity")){
			minSimilarity=Float.parseFloat(b);
		}

		//		else if(a.equalsIgnoreCase("minSizeToCluster")){
		//			minSizeToCluster=Parse.parseIntKMG(b);
		//		}else if(a.equalsIgnoreCase("minSizeToRefine")){
		//			minSizeToRefine=Parse.parseIntKMG(b);
		//		}
		else if(a.equalsIgnoreCase("minseedsize") || a.equals("minsizeseed") || a.equals("minseed")){
			minSizeToCompare=minSizeToMerge=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeToCompare")){
			minSizeToCompare=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeToMerge")){
			minSizeToMerge=Parse.parseIntKMG(b);
			//		}else if(a.equalsIgnoreCase("minSizeToAdd")){
			//			minSizeToAdd=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeResidue") || a.equalsIgnoreCase("minResidue")){
			minSizeResidue=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minNetSize")){
			minNetSize=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("midNetSize")){
			midNetSize=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("largeNetSize")){
			largeNetSize=Parse.parseIntKMG(b);
		}

		else if(a.equalsIgnoreCase("residueStringency")){
			residueStringency=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("purifyStringency")){
			purifyStringency=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("smallthresh")){
			smallThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("smallmult")){
			smallMult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("bigthresh")){
			bigThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("bigmult")){
			bigMult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("hugethresh")){
			hugeThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("hugemult")){
			hugeMult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("tinythresh")){
			tinyThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("tinyMult") || a.equalsIgnoreCase("smallPenalty")){
			tinyMult=Float.parseFloat(b);
		}

		else if(a.equalsIgnoreCase("maxEdges")){
			maxEdges=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minEdgeWeight")){
			minEdgeWeight=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minEdgeRatio")){
			minEdgeRatio=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("reciprocalEdges")){
			reciprocalEdges=Parse.parseBoolean(b);
		}

		else if(a.equalsIgnoreCase("fuseLowerLimit")){
			fuseLowerLimit=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("fuseUpperLimit")){
			fuseUpperLimit=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("fuseStringency")){
			fuseStringency=Float.parseFloat(b);
		}

		else if(a.equalsIgnoreCase("lowDepthEdgeRatio")){
			lowDepthEdgeRatio=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("highDepthEdgeRatio")){
			highDepthEdgeRatio=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("goodEdgeMult")){
			goodEdgeMult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("transEdgeMult") || a.equalsIgnoreCase("goodTransEdgeMult") || 
				a.equalsIgnoreCase("goodTransitiveEdgeMult")){
			goodTransEdgeMult=Float.parseFloat(b);
		}

		else if(a.equals("netcutoff") || a.equals("cutoff")) {
			netCutoff_small=netCutoff_mid=netCutoff_large=Float.parseFloat(b);
		}else if(a.equals("netcutoff1") || a.equals("cutoff1")) {
			netCutoff1=Float.parseFloat(b);
		}else if(a.equals("cutoffsmall") || a.equals("smallcutoff")) {
			netCutoff_small=Float.parseFloat(b);
		}else if(a.equals("cutoffmid") || a.equals("midcutoff")) {
			netCutoff_mid=Float.parseFloat(b);
		}else if(a.equals("cutofflarge") || a.equals("largecutoff")) {
			netCutoff_large=Float.parseFloat(b);
		}else if(a.equals("netcutoffupper") || a.equals("cutoffupper") || a.equals("cutoff2")) {
			netCutoffUpper=Float.parseFloat(b);
		}else if(a.equals("netcutofflower") || a.equals("cutofflower")) {
			netCutoffLower=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("netMultUpper")) {
			netMultUpper=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("netMultLower")) {
			netMultLower=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("cutoffMultA")){
			cutoffMultA=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("cutoffMultB")){
			cutoffMultB=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("cutoffMultC")){
			cutoffMultC=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("cutoffMultD")){
			cutoffMultD=Float.parseFloat(b);
		}

		//		else if(a.equalsIgnoreCase("mtcompare") || a.equalsIgnoreCase("comparemt")){
		//			multiThreadedCompare=Parse.parseBoolean(b);
		//		}

		else {return false;}

		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Graph            ----------------*/
	/*--------------------------------------------------------------*/

	public long followEdges(ArrayList<Contig> contigs, ArrayList<? extends Bin> input, Oracle oracle){
		//		outstream.print("Following Edges:  \t");
		phaseTimer.start();
		long merges=0;
		merges=launchThreads(input, null, contigs, FOLLOW_MODE, 0, 0, oracle);
		//		phaseTimer.stopAndPrint();
		//		phaseTimer.start();

		//This phase can be much slower than edge-following, presumably when clusters get big.
		if(merges>0) {merges=mergeWithDest(contigs, input);}

		if(merges>0 && loud) {phaseTimer.stop("Merged "+merges+"/"+input.size()+" bins: \t");}
		return merges;
	}

	private int mergeWithDest(ArrayList<Contig> contigs, ArrayList<? extends Bin> input){
		int x=0, y=0, z=0;

		for(Bin a : input) {
			if(!a.isCluster() && a.cluster()!=null) {a.dest=-1;continue;}
			if(a.isCluster() && a.numContigs()==0) {a.dest=-1;continue;}
			assert(a.isValid());
			if(a.dest<0) {
				//ignore
				x++;
			}else {
				y++;
				Bin b=contigs.get(a.dest);
				if(b.cluster()!=null) {b=b.cluster();}
				if(!b.isEmpty() && !a.sameCluster(b)) {
					if(a.labelTaxid>=0 && b.labelTaxid>=0) {
						long minSize=Tools.min(a.size(), b.size());
						addMerge(minSize, a.labelTaxid==b.labelTaxid);
						if(a.labelTaxid==b.labelTaxid) {
							goodMergesFollow++;
							goodMergeSizeFollow+=minSize;
						}else {
							badMergesFollow++;
							badMergeSizeFollow+=minSize;
						}
					}
					assert(b.isValid());
					if(a.isCluster()) {
						try {
							((Cluster)a).add(b);
						} catch (Throwable e) {
							System.err.println(a.numContigs()+", "+((Cluster)a).contigs);
							// TODO Auto-generated catch block
							e.printStackTrace();
							throw new RuntimeException(e);
						}
					}else if(b.isCluster()) {
						((Cluster)b).add(a);
					}else if(a.cluster()!=null){
						assert(false);
						a.cluster().add(b);
					}else {
						Cluster c=a.toCluster();
						c.add(b);
					}
					z++;
				}
			}
			a.dest=-1;
		}
		return z;
	}

	private static final int countClustered(ArrayList<? extends Bin> list) {
		int clustered=0;
		for(Bin b : list) {
			clustered+=(b.cluster()!=null ? b.numContigs() : 0);
		}
		return clustered;
	}

	public static final ArrayList<Bin> toBinList(ArrayList<? extends Bin> list, int minSize){
		ArrayList<Bin> bins=new ArrayList<Bin>();
		IntHashSet clusterSet=new IntHashSet(255);
		for(Bin a : list) {
			Cluster c=a.cluster();
			long size=(c==null ? a.size() : c.size());
			if(size<minSize) {continue;}

			if(c==null) {//Contig
				bins.add(a);
			}else if(!clusterSet.contains(c.id())){
				bins.add(c);
				clusterSet.add(c.id());
			}
			//			Bin b=bins.get(bins.size()-1);
			//			assert(b.size()>=minSize);//123
		}
		return bins;
	}

	private int followEdges(Bin a, ArrayList<Contig> contigs, Oracle oracle) {
		oracle.clear();
		ArrayList<KeyValue> edges=KeyValue.toList(a.pairMap);//TODO: Slow
		float bestScore=0;
		Bin target=null;
		assert(a.isCluster() || a.cluster()==null);

		int max=(a.isCluster() ? maxEdges : maxEdges+Tools.min(2, maxEdges)*Tools.min(8, a.numContigs()-1));
		max=Tools.min(max, edges.size());
		int minWeight=(int)Math.ceil(minEdgeRatio*edges.get(0).value);
		minWeight=Tools.max(minWeight, minEdgeWeight);
		for(int i=0; i<max; i++) {
			KeyValue kv=edges.get(i);
			if(kv.value<minWeight) {break;}
			Contig c=contigs.get(kv.key);
			Bin b=(c.cluster()!=null ? c.cluster() : c);
			if(a==b || target==b) {continue;}
			int min=Tools.min((int)c.countEdgesTo(a), kv.value);
			//			if(min>minWeight) {min=c.countReciprocalEdges(a);}//Slow

			//			verbose=(a.id()==71 && b.id()==14499);
			if(verbose) {
				System.err.println("a="+a.id()+", b="+b.id()+", c="+c.id()+", v="+kv.value);
				System.err.println("ca="+c.countEdgesTo(a)+", cb="+c.countEdgesTo(b)+
						", ac="+a.countEdgesTo(c)+", bc="+b.countEdgesTo(c));
				System.err.println(b.numContigs()+", "+b.cluster().contigSet+
						"\naMap="+a.pairMap+"\nbMap="+b.pairMap+"\ncMap="+c.pairMap);
			}
			if(min>=minWeight) {
				//				float f=a.similarityTo(b, oracle.stringency0);
				float f=oracle.similarity(a, b, 1f);
				//				assert(f==f2) : f+", "+f2+", "+a.id()+", "+b.id()+", "+min+", "+minWeight;
				if(f>bestScore) {
					target=b;
					bestScore=f;
				}
			}
			assert(!verbose);
		}
		return target==null ? -1 : target.id();
	}

	public Cluster findLinkedCluster(Bin a, ArrayList<Contig> contigs) {
		if(a.pairMap==null) {return null;}
		assert(a.isCluster() || a.cluster()==null);
		int bestValue=0;
		Cluster bestCluster=null;
		int[] keys=a.pairMap.keys(), values=a.pairMap.values();
		//		final int max=Tools.max(values);
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			if(key!=a.pairMap.invalid()) {
				int v=values[i];
				Contig a2=contigs.get(key);
				assert(key==a2.id());
				Cluster c=contigs.get(key).cluster();
				if(c!=null && v>bestValue && c!=a) {
					assert(c.contigSet.contains(key));
					if(!a.isCluster() || a.id()<=c.id) {
						if(c.countEdgesTo(a)>0) {
							bestValue=v;
							bestCluster=c;
						}
					}
				}
			}
		}
		return bestCluster;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Binning           ----------------*/
	/*--------------------------------------------------------------*/

	public BinMap makeBinMap(ArrayList<Contig> contigList, ArrayList<? extends Bin> input) {
		outstream.print("Making BinMap:    \t");
		phaseTimer.start();
		makingBinMap=true;

		if(input==null) {
			input=contigList;
		}else {
			Collections.sort(input);
		}
		BinMap map=new BinMap(contigList);
		Key key=new Key();
		float stringency=1;
		long contigsAdded=0;
		long clustersCreated=0;

		for(int i=1; i<input.size(); i++) {
			assert(input.get(i).size()<=input.get(i-1).size());
		}

		float[] ret=new float[1];
		float maxProduct=max4merDif1*maxDepthRatio1*Binner.productMult;
		Oracle oracle=new Oracle(maxGCDif1, maxDepthRatio1, max3merDif1, max4merDif1, max5merDif1, maxProduct, 
				maxCovariance1, minKmerProb2, 0);
		for(Bin a : input) {
			int initialContigs=a.numContigs();
			//			Cluster c=map.addOrMerge(a, minSizeToCompare*8, minSizeToMerge*4, minSizeToCompare,	
			//					maxKmerDif1, maxDepthRatio1, maxGCDif1, maxCovariance1, stringency, 
			//					TaxTree.SPECIES, true, true, key, ret, 0);
			Cluster c=map.addOrMerge(a, minSizeToCompare*8, minSizeToMerge*4, minSizeToCompare,	
					oracle, key, 0);
			if(c==null) {
				//residual.add(b);//Automatic
			}else {
				contigsAdded+=initialContigs;
				if(c.numContigs()==1) {
					clustersCreated++;
				}else if(c.labelTaxid>0 && a.labelTaxid>0) {
					long minSize=Tools.min(a.size(), c.size());
					addMerge(minSize, a.labelTaxid==c.labelTaxid);
					goodMergesCreate+=(c.labelTaxid==a.labelTaxid) ? 1 : 0;
					badMergesCreate+=(c.labelTaxid==a.labelTaxid) ? 0 : 1;
					if(a.labelTaxid!=c.labelTaxid) {
						badMergeSizeCreate+=minSize;
					}else {
						goodMergeSizeCreate+=minSize;
					}
				}
			}
		}
		makingBinMap=false;
		fastComparisons.addAndGet(oracle.fastComparisons);
		trimerComparisons.addAndGet(oracle.trimerComparisons);
		tetramerComparisons.addAndGet(oracle.tetramerComparisons);
		slowComparisons.addAndGet(oracle.slowComparisons);
		netComparisons.addAndGet(oracle.netComparisons);
		phaseTimer.stopAndPrint();
		outstream.println("Made "+map.map.size()+" lists containing "+clustersCreated+
				" clusters and "+contigsAdded+" contigs from "+contigList.size()+" elements.");
		return map;
	}

	public int refineBinMapPass(BinMap map, float stringency, 
			int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID, int range, int minSize) {
		//		System.err.println("Merging clusters pass.");

		float maxTrimerDif=max3merDif2*stringency;
		float maxKmerDif=max4merDif2*stringency;
		float max5merDif=max5merDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;

		ArrayList<Cluster> clusters=map.toList(false);
		Collections.sort(clusters);
		Key key=new Key();

		for(int i=1; i<clusters.size(); i++) {
			assert(clusters.get(i).size()<=clusters.get(i-1).size());
		}

		//		System.err.println("maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+", maxGCDif="+maxGCDif
		//				+", maxProduct="+maxProduct+", allowNoTaxID="+allowNoTaxID+", allowHalfTaxID="+allowHalfTaxID
		//				+", range="+range);

		int merged=0;

		Oracle oracle=new Oracle(maxGCDif, maxDepthRatio, maxTrimerDif, maxKmerDif, max5merDif, 
				maxProduct, maxCovariance, minKmerProb2, 0);
		launchThreads(clusters, map, null, REFINE_MODE, range, minSize, oracle);
		for(Cluster c : clusters) {
			synchronized(c) {
				if(c.dest>=0) {
					merged++;
				}
			}
		}
		//		assert(false);
		//		assert(map.isValid());
		if(merged<1) {return 0;}

		for(int i=clusters.size()-1; i>=0; i--) {
			Cluster a=clusters.get(i);
			synchronized(a) {
				final int dest=a.dest;
				if(dest>=0 && dest!=a.id()) {
					Cluster b=map.contigList.get(dest).cluster;
					synchronized(b) {
						if(b!=a) {
							assert(!b.contigSet.contains(a.id));
							assert(!a.contigSet.contains(b.id));
							assert(a.id()!=dest && a.id()!=b.id()) : a.id+", "+dest+", "+b.id+", "+(a.id()==dest)+", "+(b.id()==dest)+", "+(a.id()!=b.id());

							if(a.labelTaxid>0 && b.labelTaxid>0) {
								long smaller=Tools.min(a.size(), b.size());
								addMerge(smaller, a.labelTaxid==b.labelTaxid);
								goodMergesRefine+=(a.labelTaxid==b.labelTaxid) ? 1 : 0;
								badMergesRefine+=(a.labelTaxid==b.labelTaxid) ? 0 : 1;
								if(a.labelTaxid!=b.labelTaxid) {
									badMergeSizeRefine+=smaller;
								}else {
									goodMergeSizeRefine+=smaller;
								}
							}
							b.add(a);
							clusters.set(i, null);
						}
					}
				}
			}
		}

		map.clear(false);
		//		assert(map.isValid());
		Tools.condense(clusters);
		Collections.sort(clusters);
		for(Cluster c : clusters) {
			map.add(c, key);
		}
		assert(map.isValid());
		return merged;
	}

	public int recluster(BinMap map, float stringency, int minSizeRecluster) {

		float maxTrimerDif=max3merDif2*stringency;
		float maxKmerDif=max4merDif2*stringency;
		float max5merDif=max5merDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;
		assert(isValid(map.toList(true), false));
		
		// Collect clusters worth reclustering
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		{
			ArrayList<Cluster> clusters0 = map.toList(false);
			Collections.sort(clusters0);
			for(Cluster c : clusters0) {
				if(c.size < minSizeRecluster) break;
				if(c.numContigs() > 3) { // Need minimum size for meaningful splits
					clusters.add(c);
				}
			}
		}
		System.err.println("Recluster Targets: "+clusters.size());
		if(clusters.isEmpty()){return 0;}

		assert(isValid(map.toList(true), false));
		Oracle oracle=new Oracle(maxGCDif, maxDepthRatio, maxTrimerDif, maxKmerDif, max5merDif, 
				maxProduct, maxCovariance, minKmerProb2*0.8f, 0);

		// Launch threads for parallel reclustering
		launchThreads(clusters, map, null, RECLUSTER_MODE, 0, 0, oracle);

		// Collect results and rebuild map if changes made
		int splits=0;
		for(Cluster clust : clusters) {
			synchronized(clust) {
				if(clust.wasReclustered) { // Flag set by thread
					splits++;
				}
			}
		}
		System.err.println("Reclustered: "+splits);

		if(splits>0) {
			map.clear(true);
			ArrayList<Bin> bins = toBinList(map.contigList, 0);
			Collections.sort(bins);
			Key key = new Key();
			for(Bin b : bins) {map.add(b, key);}
		}
		assert(isValid(map.toList(true), false));
		return splits;
	}

	public int purify(BinMap map, float stringency, int range, int minSizePurify, int minSizeCompare) {
		//		System.err.println("Merging clusters pass.");

		float maxTrimerDif=max3merDif2*stringency;
		float maxKmerDif=max4merDif2*stringency;
		float max5merDif=max5merDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;

		for(Contig c : map.contigList) {
			synchronized(c) {c.dest=-1;}
		}

		ArrayList<Cluster> clusters=new ArrayList<Cluster>();
		{
			ArrayList<Cluster> clusters0=map.toList(false);
			Collections.sort(clusters0);
			for(Cluster c : clusters0) {
				if(c.size<minSizePurify) {break;}
				if(c.numContigs()>1) {
					clusters.add(c);
				}
			}
		}
		//		System.err.println("Launching purify on "+clusters.size()+" clusters!");
		if(clusters.isEmpty()) {return 0;}

		int removed=0;

		Oracle oracle=new Oracle(maxGCDif, maxDepthRatio, maxTrimerDif, maxKmerDif, max5merDif, 
				maxProduct, maxCovariance, minKmerProb2*0.8f, 0);
		oracle.useEdges=false;
		launchThreads(clusters, map, null, PURIFY_MODE, range, minSizeCompare, oracle);
		for(Cluster clust : clusters) {
			synchronized(clust) {
				removed+=purifyCluster(clust);
			}
		}
		if(removed>0) {
			map.clear(true);
			ArrayList<Bin> bins=toBinList(map.contigList, 0);
			Collections.sort(bins);
			Key key=new Key();
			for(Bin b : bins) {
				map.add(b, key);
			}
		}
		return removed;
	}

	int purifyCluster(Cluster clust) {
		ArrayList<Contig> contigs=clust.contigs;
		int removed=0;
		for(int i=0; i<contigs.size(); i++) {
			Contig c=contigs.get(i);
			if(c.dest>=0 && c.dest!=clust.id) {
				c.cluster=null;
				removed++;
				contigs.set(i, null);
			}
		}
		if(removed>0) {
			@SuppressWarnings("unchecked")
			ArrayList<Contig> list=(ArrayList<Contig>) contigs.clone();
			clust.clear();
			for(Contig c : list) {
				if(c!=null) {
					c.cluster=null;
					clust.add(c);
				}
			}
			assert(clust.contigSet.contains(clust.id));
			if(clust.numDepths()>1) {clust.fillNormDepth();}
		}
		return removed;
	}

	public int processResidue(BinMap map, float stringency, 
			int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID, int range) {
		//		assert(map.isValid());
		System.err.println("Processing "+map.residual.size()+" residual contigs.");
		Timer t=new Timer(outstream, true);
		if(map.residual.isEmpty()) {return 0;}

		float maxTrimerDif=max3merDif2*stringency;
		float maxKmerDif=max4merDif2*stringency;
		float max5merDif=max5merDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;

		//		System.err.println("maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+", maxGCDif="+maxGCDif
		//				+", maxProduct="+maxProduct+", allowNoTaxID="+allowNoTaxID+", allowHalfTaxID="+allowHalfTaxID
		//				+", range="+range);
		int merged=0;
		map.residual=toBinList(map.residual, 0);
		//		assert(map.isValid());

		int minSize=Tools.max(minSizeToMerge, minClusterSize/5);
		Oracle oracle=new Oracle(maxGCDif, maxDepthRatio, maxTrimerDif, maxKmerDif, max5merDif, 
				maxProduct, maxCovariance, minKmerProb2, 0);
		launchThreads(map.residual, map, null, RESIDUE_MODE, range, minSize, oracle);
		for(Bin c : map.residual) {
			synchronized(c) {
				if(c.dest>=0) {merged++;}
			}
		}
		t.stopAndStart("Found "+merged+" merge targets.");

		if(merged<1) {return 0;}
		for(int i=0; i<map.residual.size(); i++) {
			Bin a=map.residual.get(i);
			if(a.dest>0) {
				Cluster b=map.contigList.get(a.dest).cluster;
				assert(a!=b) : "\n"+a.id()+"\n"+b.id()+"\n"+a.dest+"\n"+
				a.getClass()+"\n"+b.getClass()+"\n"+b.contigSet.contains(a.id());
				if(a.labelTaxid>0 && b.labelTaxid>0) {
					long smaller=Tools.min(a.size(), b.size());
					addMerge(smaller, a.labelTaxid==b.labelTaxid);
					goodMergesResidue+=(a.labelTaxid==b.labelTaxid) ? 1 : 0;
					badMergesResidue+=(a.labelTaxid==b.labelTaxid) ? 0 : 1;
					if(a.labelTaxid!=b.labelTaxid) {
						badMergeSizeResidue+=smaller;
					}else {
						goodMergeSizeResidue+=smaller;
					}
				}
				b.add(a);//Note at this point b could be a residue that becomes bigger than residue size
				map.residual.set(i, null);
			}
		}

		Tools.condenseStrict(map.residual);
		for(ArrayList<Cluster> list : map.map.values()) {
			Collections.sort(list);
		}
		assert(map.isValid());

		t.stop("Merged "+merged+" contigs into clusters.");
		return merged;
	}

	public Cluster findBestResidualCluster(Bin a, BinMap map, Key key, Oracle oracle, 
			int range, int minSize) {
		if(a==null || a.size()<minSizeResidue) {return null;}
		int minSize2=(int)Tools.max(minSizeToCompare, a.size(), minSize);
		Cluster b=map.findBestCluster(a, minSize2, key, range, oracle);
		return b;
	}

	public int refineBinMap(BinMap map) {
		System.err.println("Merging clusters.");
		if(sketchClusters) {sketcher.sketch(map.toList(false), false);}
		else {
			for(Cluster c : map) {
				if(c.sketchedSize()>=2*c.size()) {c.clearTax();}
			}
		}
		phaseTimer.start();

		int removedThisPhase=0;
		int removedTotal=0;
		//		net0small.cutoff+=0.02;
		//		net0mid.cutoff+=0.02;
		//		net0large.cutoff+=0.02;
		if(sketchContigs || sketchClusters) {
			removedThisPhase=refinePhase(map, "aa", 2.5f, TaxTree.SPECIES, false, false, baseRange, minSizeToMerge, 8);
			removedTotal+=removedThisPhase;

			removedThisPhase=refinePhase(map, "bb", 1.5f, TaxTree.SPECIES, true, false, baseRange, minSizeToMerge, 8);
			removedTotal+=removedThisPhase;

			removedThisPhase=refinePhase(map, "cc", 1.0f, TaxTree.GENUS, true, true, baseRange, minSizeToMerge, 8);
			removedTotal+=removedThisPhase;
		}else {
			if(runPassAA) {
				removedThisPhase=refinePhase(map, "aa", 0.8f, -1, true, true, baseRange+1, Tools.max(4*minSizeToMerge, 20000), basePasses+1);
				removedTotal+=removedThisPhase;
			}
			if(runPassA) {
				removedThisPhase=refinePhase(map, "a", 0.8f, -1, true, true, baseRange, minSizeToMerge, basePasses+1);
				removedTotal+=removedThisPhase;
			}
			if(runPassB) {
				removedThisPhase=refinePhase(map, "b", 0.85f, -1, true, true, baseRange+2, minSizeToMerge, basePasses+2);
				removedTotal+=removedThisPhase;
			}
		}
		//		net0small.cutoff-=0.02;
		//		net0mid.cutoff-=0.02;
		//		net0large.cutoff-=0.02;
		if(runPassD) {
			removedThisPhase=refinePhase(map, "d", 0.9f, -1, true, true, baseRange, minSizeToMerge, basePasses+0);
			removedTotal+=removedThisPhase;
		}
		if(runPassE) {
			removedThisPhase=refinePhase(map, "e", 1.0f, -1, true, true, baseRange+2, minSizeToMerge, basePasses+3);
			removedTotal+=removedThisPhase;
		}

		phaseTimer.stop("Refinement merged "+removedTotal+" clusters. ");
		return removedTotal;
	}

	int refinePhase(BinMap map, String phase,
			float stringency, int taxLevel, boolean noTax, boolean halfTax, int range, int initialMinSize, int passes) {
		Timer t=new Timer(outstream, true);
		int removedThisPhase=0;
		for(int pass=1; pass<=passes; pass++) {
			int initial=map.countClusters();
			int removed=refineBinMapPass(map, stringency, taxLevel, noTax, halfTax, range, initialMinSize*pass);
			removedThisPhase+=removed;
			if(loud) {
				System.err.print("Refinement Pass "+pass+phase+": Merged "+removed+"/"+initial+" clusters. ");
				t.stopAndStart("\t");
			}
			if(removed<2) {break;}
		}
		if(sketchClusters && removedThisPhase>0) {sketcher.sketch(map.toList(false), false);}
		return removedThisPhase;
	}

	public ArrayList<Bin> clusterByTaxid(ArrayList<? extends Bin> bins){
		outstream.print("Clustering by Taxid: \t");
		phaseTimer.start();
		Collections.sort(bins);
		HashMap<Integer, Bin> map=new HashMap<Integer, Bin>();
		int clustersMade=0;
		int contigsClustered=0;

		ArrayList<Bin> out=new ArrayList<Bin>();
		for(int i=0; i<bins.size(); i++) {
			Bin b=bins.get(i);
			if(b.taxid()>0) {
				Integer key=Integer.valueOf(b.taxid());
				Bin old=map.get(key);
				if(old==null) {
					map.put(key, b);
					b=null;
				}else if(old.getClass()==Cluster.class) {
					//Todo: Write "similar" function.
					((Cluster)old).add(b);
					contigsClustered++;
					b=null;
				}else {
					Cluster a=new Cluster(clustersMade);
					a.add(old);
					a.add(b);
					map.put(key, a);
					clustersMade++;
					contigsClustered++;
					b=null;
				}
			}
			if(b!=null) {out.add(b);}
		}

		out.addAll(map.values());
		Collections.sort(out);

		phaseTimer.stopAndPrint();
		outstream.println("Made "+clustersMade+" clusters containing "+contigsClustered+"/"+bins.size()+" elements.");
		return out;
	}

	public long fuse(ArrayList<Contig> contigs, ArrayList<? extends Bin> input, float stringency){
		if(loud) {outstream.print("Initiating Fusion:  \t");}
		phaseTimer.start();

		//		for(int i=1; i<input.size(); i++) {
		//			assert(input.get(i).size()<=input.get(i-1).size());
		//			assert(input.get(i).size()>=fuseLowerLimit);
		//		}
		//		assert(false) : fuseLowerLimit;

		float maxTrimerDif=max3merDif2*stringency;
		float maxKmerDif=max4merDif2*stringency;
		float max5merDif=max5merDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;
		Oracle oracle=new Oracle(maxGCDif, maxDepthRatio, maxTrimerDif, maxKmerDif, max5merDif, 
				maxProduct, maxCovariance, minKmerProb2, 0);

		//		System.err.println("maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+
		//				", maxGCDif="+maxGCDif+", maxProduct="+maxProduct+
		//				", maxCovariance="+maxCovariance+", minKmerProb2="+minKmerProb2);
		//		System.err.println("List size: "+input.size());

		long merges=0;
		merges=launchThreads(input, null, contigs, FUSE_MODE, 0, 0, oracle);
		//		phaseTimer.stopAndPrint();
		//		phaseTimer.start();

		//This phase can be slow.
		if(merges>0) {merges=mergeWithDest(contigs, input);}

		if(merges>0 && loud) {phaseTimer.stop("Merged "+merges+"/"+input.size()+" bins: \t");}
		return merges;
	}

	/*--------------------------------------------------------------*/

	void setSamples(int samples, float mult) {
		if(overrideSetSamples>0) {samples=overrideSetSamples;}
		System.err.println("Setting cutoffs for "+samples+" samples.");
		if(samples<2) {//Single mode //maxkmerdif=0.0060 maxgcdif=0.045
			max4merDif2=0.0048f;
			//			maxKmerDif2=0.008f;
			maxDepthRatio2=1.38f;
			maxGCDif2=0.032f;
			minKmerProb2=0.82f;
			max5merDif2=0.007f;

			fuseStringency=1.6f;
		}else if(samples<3){//Two mode
			max4merDif2=0.0065f;
			maxDepthRatio2=1.32f;
			maxGCDif2=0.032f;
			maxCovariance2=0.0040f;//Less than .35 greatly decreases CAMI contam, but less than 0.5 decreases synth3 completeness 
			minKmerProb2=0.80f;
			max5merDif2=0.008f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
		}else if(samples<4){//Three mode
			max4merDif2=0.0075f;
			maxDepthRatio2=1.40f;
			maxGCDif2=0.04f;
			maxCovariance2=0.0050f;
			minKmerProb2=0.75f;
			max5merDif2=0.009f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
		}else if(samples<5){//Four mode
			max4merDif2=0.0085f;
			maxDepthRatio2=1.4f;
			maxGCDif2=0.05f;
			maxCovariance2=0.0045f;
			minKmerProb2=0.7f;
			max5merDif2=0.010f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
			purifyStringency=3.0f;
			//			bigMult=0.70f;
			//			hugeMult=0.3f;
		}else if(samples<6){//5
			max4merDif2=0.0095f;
			maxDepthRatio2=1.45f;
			maxGCDif2=0.05f;
			maxCovariance2=0.0035f;
			minKmerProb2=0.7f;
			max5merDif2=0.011f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
			purifyStringency=2.5f;
			//			bigMult=0.70f;
			//			hugeMult=0.25f;
		}else if(samples<7){//6
			max4merDif2=0.0105f;
			maxDepthRatio2=1.45f;
			maxGCDif2=0.05f;
			maxCovariance2=0.0035f;
			minKmerProb2=0.7f;
			max5merDif2=0.012f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
			purifyStringency=2.5f;
			//			bigMult=0.70f;
			//			hugeMult=0.25f;
		}else{//7+
			max4merDif2=0.0115f;
			maxDepthRatio2=1.45f;
			maxGCDif2=0.05f;
			maxCovariance2=0.0035f;
			minKmerProb2=0.7f;
			max5merDif2=0.015f;

			netCutoffUpper=0.65f;
			netCutoff_small=netCutoff_mid=netCutoff_large=0.5f;
			fuseStringency=1.6f;
			purifyStringency=2.5f;
			//			bigMult=0.70f;
			//			hugeMult=0.20f;
		}

		if(mult!=1) {
			float lowMult=((mult-1)*0.5f)+1;
			float highMult=((mult-1)*2f)+1;
			max4merDif2*=lowMult;
			maxDepthRatio2=1+(maxDepthRatio2-1)*(mult>1 ? lowMult : mult);
			maxGCDif2*=mult;
			maxCovariance2*=(mult);
			purifyStringency*=mult;
			fuseStringency=1+(fuseStringency-1)*mult;
			smallMult=1+(smallMult-1)*mult;


			minKmerProb2=1-smallMult*(1-minKmerProb2);
			minKmerProb2=Tools.mid(0, minKmerProb2, 1);

			if(mult>1) {
				hugeMult=1/(1+(1/hugeMult-1)/highMult);
			}else {
				hugeMult=1/(1+(1/hugeMult-1)/lowMult);
			}
		}
		max3merDif2=max4merDif2*0.625f;
		max5merDif2=max4merDif2*2.000f;
	}

	void printThresholds(){
		System.err.println("maxTrimerDif:     "+max3merDif2);
		System.err.println("maxTetramerDif:   "+max4merDif2);
		System.err.println("minKmerProb:      "+minKmerProb2);
		System.err.println("maxDepthRatio:    "+maxDepthRatio2);
		System.err.println("maxGCDif:         "+maxGCDif2);
		System.err.println("maxCovariance:    "+maxCovariance2);
		System.err.println("netCutoff:        "+netCutoff_large);
		System.err.println("netCutoffUpper:   "+netCutoffUpper);
		System.err.println("fuseStringency:   "+fuseStringency);
		System.err.println("purifyStringency: "+purifyStringency);
		//		System.err.println("smallThresh:      "+smallThresh);
		System.err.println("smallMult:        "+smallMult);
		//		System.err.println("bigThresh:        "+bigThresh);
		System.err.println("bigMult:          "+bigMult);
		System.err.println("hugeMult:         "+hugeMult);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Classes and Threading    ----------------*/
	/*--------------------------------------------------------------*/

	private synchronized long launchThreads(ArrayList<? extends Bin> list, BinMap map, 
			ArrayList<Contig> contigs, int mode, int range, int minSize, Oracle oracle) {

		//Do anything necessary prior to processing

		//Determine how many threads may be used
		int threads=Tools.mid(1, list.size()/64, Shared.threads());
		if(compareThreadsOverride>0) {threads=compareThreadsOverride;}

		//		if(mode==FUSE_MODE) {threads=1;}//123
		threadMerges=0;

		//Fill a list with LoadThreads
		ArrayList<CompareThread> alpt=new ArrayList<CompareThread>(threads);
		for(int i=0; i<threads; i++){
			CompareThread ct=new CompareThread(list, map, contigs, i, 
					threads, mode, range, minSize, oracle.clone());
			alpt.add(ct);
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		assert(success);
		//Do anything necessary after processing

		return threadMerges;
	}

	@Override
	public synchronized void accumulate(CompareThread t) {
		synchronized(t){
			fastComparisons.addAndGet(t.oracle.fastComparisons);
			trimerComparisons.addAndGet(t.oracle.trimerComparisons);
			tetramerComparisons.addAndGet(t.oracle.tetramerComparisons);
			slowComparisons.addAndGet(t.oracle.slowComparisons);
			netComparisons.addAndGet(t.oracle.netComparisons);
			threadMerges+=t.mergesT;
			errorState|=(!t.success);
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public boolean success() {return !errorState;}

	class CompareThread extends Thread {

		public CompareThread(ArrayList<? extends Bin> list, BinMap map_, ArrayList<Contig> contigs_, 
				int tid_, int threads_, int mode_, int range_, int minSize_, Oracle oracle_) {
			input=list;
			map=map_;
			contigs=(contigs_!=null ? contigs_ : map.contigList);
			tid=tid_;
			threads=threads_;
			mode=mode_;
			range=range_;
			minSize=minSize_;
			oracle=oracle_;
		}


		public void run() {
			synchronized(this) {
				if(mode==REFINE_MODE) {
					refine();
				}else if(mode==RESIDUE_MODE) {
					residue();
				}else if(mode==PURIFY_MODE) {
					purify();
				}else if(mode==FOLLOW_MODE) {
					follow();
				}else if(mode==FUSE_MODE) {
					fuse();
				}else if(mode==RECLUSTER_MODE) {
					recluster();
				}else {
					throw new RuntimeException("Bad mode: "+mode);
				}
			}
			success=true;
		}

		private int follow() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				assert(a.dest<0);
				if(a.pairMap!=null) {
					int dest=followEdges(a, contigs, oracle);
					if(dest>=0) {
						mergesT++;
						synchronized(a) {a.dest=dest;}
					}
				}
			}
			return mergesT;
		}

		private int fuse() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				if(a.size()<fuseLowerLimit) {
					assert(false);
					break;
				}
				if(a.size()>fuseUpperLimit) {continue;}
				fuseSeeks++;
				Bin b=findBestFuseTarget(a);
				mergesT+=(b==null ? 0 : 1);
				fuseTargets+=(b==null ? 0 : 1);
			}
			//			System.err.println("fuseSeeks="+fuseSeeks);
			//			System.err.println("fuseCompares="+fuseCompares);
			//			System.err.println("fuseTargets="+fuseTargets);
			return mergesT;
		}

		private Bin findBestFuseTarget(Bin a) {
			oracle.best=null;
			oracle.topScore=-1;
			for(Bin b : input) {
				if(b.size()<fuseLowerLimit) {
					assert(false) : a.size()+", "+b.size();
					break;
				}
				if(b.size()<a.size()) {break;}
				if(b==a) {continue;}
				fuseCompares++;
				float f=oracle.similarity(a, b, 1);
				//				if(f<=0) {
				//					verbose=true;
				//					oracle.similarity(a, b, 1);
				//					assert(false);
				//				}
				if(f>oracle.topScore) {
					assert(f>0) : f;//Actually, could be -1; or clear should set to -1
					oracle.best=b;
					oracle.topScore=f;
				}
			}
			//			System.err.print('.');
			if(oracle.best!=null) {a.dest=oracle.best.id();}
			return oracle.best;
		}

		private int recluster() {
			AbstractRefiner refiner = AbstractRefiner.makeRefiner(oracle, AbstractRefiner.DEFAULT_TYPE);

			for(int i=tid; i<input.size(); i+=threads) {
				Cluster cluster = (Cluster)input.get(i);
				synchronized(cluster) {
					ArrayList<Bin> refined = refiner.refine(cluster);
					if(refined != null && refined.size() > 1) {
						// Set dest for each contig to its new cluster's representative (lowest ID)
						for(Bin newCluster : refined) {
							int representativeId = Integer.MAX_VALUE;
							// Find lowest contig ID in this refined cluster
							for(Contig c : newCluster) {
								representativeId = Math.min(representativeId, c.id());
							}
							// Set all contigs in this refined cluster to point to representative
							for(Contig c : newCluster) {
								c.dest = representativeId;
								cluster.wasReclustered=true;
							}
						}
						mergesT++;
					}
				}
			}
			return mergesT;
		}

		private int purify() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				assert(a.isCluster() && a.numContigs()>1);
				synchronized(a) {
					mergesT+=purifyCluster((Cluster)a);
				}
			}
			return mergesT;
		}

		private int purifyCluster(Cluster clust) {
			int changes=0;
			clust.dest=-1;
			@SuppressWarnings("unchecked")
			ArrayList<Contig> contigs=(ArrayList<Contig>) clust.contigs.clone();
			for(Contig a : contigs) {
				synchronized(a) {
					a.dest=-1;
					a.score=oracle.similarity(clust, a, 1);
					if(a.score<0) {
						a.dest=a.id();
						changes++;
						//						synchronized(Binner.class) {
						//							oracle.verbose2=true;
						//							oracle.similarity(clust, a, 1);
						//							assert(false);
						//						}
					}
				}
			}
			Collections.sort(contigs, ScoreComparator.comparator);
			for(Contig a : contigs) {
				if(a.id()!=clust.id()) {
					Cluster b=map.findBestCluster(a, minSize, key, range, oracle);
					if(b!=null && b!=clust && oracle.score>=1.125f*a.score) {
						synchronized(a) {
							if(a.dest<0) {changes++;}
							a.dest=b.id();
						}
						//						System.err.print(".");
					}else {
						break;
					}
				}
			}
			return changes;
		}

		private void residue() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				synchronized(a) {
					a.dest=-1;
					Cluster b=findBestResidualCluster(a, map, key, oracle, range, minSize);

					assert(a!=b);
					if(b==null) {
						b=findLinkedCluster(a, map.contigList);
						assert(a!=b);
					}

					if(b!=null) {
						a.dest=b.id;
						assert(a!=b);
						assert(a.cluster()!=b);
						assert(map.contigList.get(b.id).cluster==b);
					}
				}
			}
		}

		private void refine() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				synchronized(a) {
					a.dest=-1;
					if(a.size()>=minSizeToMerge) {
						refineBin(a);
					}
				}
			}
		}

		private void refineBin(Bin a) {
			a.dest=-1;
			int minSize2=(int)Tools.mid(minSize, a.size(), Integer.MAX_VALUE);
			Cluster b=map.findBestCluster(a, minSize2, key, range, oracle);

			if(b==null) {
				//do nothing
			}else {
				//				assert(a.size()<=b.size());
				assert(a!=b);
				assert(a.id()!=b.id);
				assert(!b.contigSet.contains(a.id()));
				assert(a.cluster()==null || !a.cluster().contigSet.contains(b.id));
				a.dest=b.id;
			}
		}

		final ArrayList<? extends Bin> input;
		final BinMap map;
		final ArrayList<Contig> contigs;
		final int tid;
		final int threads;
		final int mode;


		long fuseCompares=0;
		long fuseSeeks=0;
		long fuseTargets=0;

		//		final float maxKmerDif;
		//		final float maxDepthRatio;
		//		final float maxGCDif;
		//		final float maxProduct;
		//		final float maxCovariance;
		//		final int taxLevel;
		//		final boolean allowNoTaxID;
		//		final boolean allowHalfTaxID;
		final int range;
		final int minSize;
		final Key key=new Key();
		final float[] ret=new float[1];
		final Oracle oracle;
		boolean success=false;
		int mergesT=0;
	}

	static float sizeAdjustMult(long size) {
		float f=sizeAdjustMult2(size);
		if(size<tinyThresh){// && size>minSizeToCompare) {
			f*=tinyMult;//This is to correct the exemption from residueStringency
		}
		return f;
	}

	static float sizeAdjustMult2(long size) {
		if(size<smallThresh) {return 1f+smallMult*(smallThresh-size)/(float)smallThresh;}
		if(size>2*hugeThresh) {return hugeMult;}
		if(size>hugeThresh) {
			float range=1f-hugeMult;
			return Tools.min(bigThresh, 1f-(size-hugeThresh)*range/hugeThresh);
		}
		if(size>2*bigThresh) {return bigMult;}
		if(size>bigThresh) {
			float range=1f-bigMult;
			return 1f-(size-bigThresh)*range/bigThresh;
		}
		return 1f;
	}

	void addMerge(long size, boolean good) {
		int idx=(int)(Tools.log2(size));
		if(good) {goodMergeSize[idx]++;}
		else {badMergeSize[idx]++;}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	boolean multiThreadedCompare=true;

	//	long refinementComparisons=0;
	//	long refinementComparisonsSlow=0;
	BinSketcher sketcher;
	int baseRange=1;
	int basePasses=1;
	int residueRange=3;
	boolean errorState=false;
	boolean runPassAA=false;
	boolean runPassA=true;
	boolean runPassB=false;
	boolean runPassC=true;
	boolean runPassD=true;
	boolean runPassE=true;
	boolean runPassF=false;
	boolean runPassG=true;

	long goodMergesFollow=0;
	long badMergesFollow=0;
	long goodMergeSizeFollow=0;
	long badMergeSizeFollow=0;
	long goodMergesCreate=0;
	long badMergesCreate=0;
	long goodMergeSizeCreate=0;
	long badMergeSizeCreate=0;
	long goodMergesRefine=0;
	long badMergesRefine=0;
	long goodMergeSizeRefine=0;
	long badMergeSizeRefine=0;
	long goodMergesResidue=0;
	long badMergesResidue=0;
	long goodMergeSizeResidue=0;
	long badMergeSizeResidue=0;

	long[] goodMergeSize=new long[40];
	long[] badMergeSize=new long[40];

	BinMap binMap;

	/*--------------------------------------------------------------*/

	private int threadMerges=0;
	public AtomicLong fastComparisons=new AtomicLong(0);
	public AtomicLong trimerComparisons=new AtomicLong(0);
	public AtomicLong tetramerComparisons=new AtomicLong(0);
	public AtomicLong slowComparisons=new AtomicLong(0);
	public AtomicLong netComparisons=new AtomicLong(0);

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final Timer phaseTimer=new Timer();
	final PrintStream outstream;

	/*--------------------------------------------------------------*/
	/*----------------            Static Fields          ----------------*/
	/*--------------------------------------------------------------*/

	static boolean PERFECT_ORACLE=false;
	static boolean BAN_BAD_MERGES=false;

	static final int REFINE_MODE=0;
	static final int RESIDUE_MODE=1;
	static final int PURIFY_MODE=2;
	static final int FOLLOW_MODE=3;
	static final int FUSE_MODE=4;
	static final int RECLUSTER_MODE=5;

	static int fuseLowerLimit=5000;
	static int fuseUpperLimit=900000;
	static float fuseStringency=1.5f;
	static float purifyStringency=3f;

	static float residueStringency=0.65f;
	static float productMult=0.68f;
	static float minSimilarity=0.0f;

	static int maxEdges=2;
	static int minEdgeWeight=2;
	static boolean reciprocalEdges=true;
	static float minEdgeRatio=0.4f;

	static float lowDepthEdgeRatio=0.2f;
	static float highDepthEdgeRatio=2f;
	static float goodEdgeMult=1.35f;
	static float goodTransEdgeMult=1.25f;

	static int hugeThresh=1200000;
	static float hugeMult=0.375f;
	static int bigThresh=100000;
	static float bigMult=0.725f;
	static int smallThresh=8000;
	static float smallMult=2.0f;
	static int tinyThresh=1000;
	static float tinyMult=0.72f;

	//	static int minSizeToCluster=0;
	//	static int minSizeToRefine=500;
	/** Size of the bigger one */
	static int minSizeToCompare=2500;
	/** Size of the smaller one being compared */
	static int minSizeToMerge=2500;
	static int minSizeResidue=200;
	static int minNetSize=200;
	static int midNetSize=3000;
	static int largeNetSize=15000;
	static float netCutoffUpper=0.65f;
	static float netCutoffLower=0.547f;
	static float netCutoff_small=0.52f;
	static float netCutoff_mid=0.52f;
	static float netCutoff_large=0.52f;
	static float netCutoff1=0.65f;
	static float netMultUpper=1.4f;
	static float netMultLower=0.5f;
	static float cutoffMultA=2.7f;
	static float cutoffMultB=1.7f;
	static float cutoffMultC=1.6f;
	static float cutoffMultD=1.2f;

	static int overrideSetSamples=-1;
	static int compareThreadsOverride=-1;

	//Optimal selection when forming clusters
	static float max3merDif1=0.1f;
	static float max4merDif1=0.002f;
	static float max5merDif1=0.003f;
	static float maxDepthRatio1=1.05f;
	static float maxGCDif1=0.015f;
	static float maxCovariance1=0.0001f;

	//When merging clusters
	static float max3merDif2=0.1f;
	static float max4merDif2=0.005f; //.005 for k=4; .012 for k=5; .008 for Euclid
	static float max5merDif2=0.007f;
	static float maxDepthRatio2=1.36f;
	static float maxGCDif2=0.03f; //0.02f for 1 depth, 0.05 for 4 depths
	static float maxCovariance2=0.004f;
	static float minKmerProb2=0.9f;

}
