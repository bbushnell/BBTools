package bin;

import java.util.Arrays;

import clade.Clade;
import json.JsonObject;
import map.IntHashMap2;
import map.IntLongHashMap;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;

/**
 * Abstract base class for genomic bins containing contigs.
 * Provides depth analysis, taxonomic annotation, and similarity calculations
 * for metagenomic binning operations.
 * @author Brian Bushnell
 */
public abstract class Bin extends BinObject implements Sketchable, Iterable<Contig> {
	
	@Override
	public final int taxid() {return taxid;}
	
	/** Returns the name identifier for this bin */
	abstract String name();

	@Override
	public final float gc() {return gcSum/(float)size();}
	
	/** Clears depth measurements and normalized depth cache.
	 * Forces recalculation of depth-related metrics on next access. */
	public final void clearDepth() {
		depth.clear();
		normDepth=null;
	}
	
	/**
	 * Calculates depth ratio using specified method.
	 *
	 * @param b The bin to compare against
	 * @param method Calculation method (0-7): 1=max ratio, 2=harmonic mean,
	 * 3=adjusted max, 4=min ratio, 5-7=inverted versions of 2-4
	 * @param cutoff Early exit when ratio exceeds this
	 * @return Depth ratio metric based on method
	 */
	public float depthRatio(Bin b, int method, float cutoff) {
		if(method==0) {return 0;}
		if(method==1) {return depthRatio1(b, cutoff);}
		if(method==2) {return depthRatio2(b);}
		if(method==3) {return depthRatio3(b, cutoff);}
		if(method==4) {return depthRatio4(b);}
		if(method==5) {return 1-depthRatio2(b);}
		if(method==6) {return 1-depthRatio3(b, cutoff);}
		if(method==7) {return 1-depthRatio4(b);}
		throw new RuntimeException(""+method);
	}
	
	/**
	 * Calculates maximum depth ratio between bins across all samples.
	 * Uses depth boost to avoid division by zero.
	 * @param b The bin to compare against
	 * @param cutoff Early exit when ratio exceeds this
	 * @return Maximum ratio of depths across samples
	 */
	public float depthRatio1(Bin b, float cutoff) {
		float max=1;
		for(int i=0; i<depth.size && max<=cutoff; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=Tools.max(d1,d2)/Tools.min(d1,d2);
			max=Tools.max(max, ratio);
		}
		return max;
	}
	
	/**
	 * Calculates minimum harmonic mean of depth ratios across samples.
	 * Formula: 2*min(d1,d2)/(d1+d2)
	 * @param b The bin to compare against
	 * @return Minimum harmonic mean ratio
	 */
	private float depthRatio2(Bin b) {
		float min=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=2*Tools.min(d1,d2)/(d1+d2);
			min=Tools.min(min, ratio);
		}
		return min;
	}
	
	/**
	 * Calculates maximum adjusted depth ratio across samples.
	 * Formula: 2*(max(d1,d2)/(d1+d2))-1
	 * @param b The bin to compare against
	 * @return Maximum adjusted depth ratio
	 */
	private float depthRatio3(Bin b, float cutoff) {
		float max=0;
		for(int i=0; i<depth.size && max<=cutoff; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=2*(Tools.max(d1,d2)/(d1+d2))-1;
			max=Tools.max(max, ratio);
		}
		return max;
	}
	
	/**
	 * Calculates minimum simple depth ratio across samples.
	 * Formula: min(d1,d2)/max(d1,d2)
	 * @param b The bin to compare against
	 * @return Minimum simple depth ratio
	 */
	private float depthRatio4(Bin b) {
		float min=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=Tools.min(d1,d2)/Tools.max(d1,d2);
			min=Tools.min(min, ratio);
		}
		return min;
	}
	
	/**
	 * Returns whether this bin is taxonomically pure using 0.02 fraction threshold
	 */
	public boolean pure() {return pure(0.02f);}
	/**
	 * Returns whether this bin is taxonomically pure.
	 * Currently simplified to check if labelTaxid is positive.
	 * @param fraction Purity threshold fraction (not currently used)
	 * @return true if bin has a valid taxonomic label
	 */
	public boolean pure(float fraction) {return labelTaxid>0;}
	
	/**
	 * Sets depth value for a specific sample.
	 * @param d Depth value to set
	 * @param sample Sample index
	 */
	public final void setDepth(float d, int sample) {
		depth.set(sample, d);
	}
	
	/** Appends a depth value to the depth list.
	 * @param d Depth value to append */
	public final void appendDepth(float d) {
		depth.add(d);
	}
	
	/** Returns the number of depth measurements stored */
	public final int numDepths() {
		return depth.size;
	}
	
	/** Returns the number of pair-wise edges in the pair map */
	public final float numEdges() {
		return pairMap==null ? 0 : pairMap.size();
	}
	
	/**
	 * Returns depth for a specific sample.
	 * Special handling for single-sample case where sample 1 returns 0.
	 * @param sample Sample index
	 * @return Depth value for the sample
	 */
	public final float depth(int sample) {
		return (sample==1 && depth.size==1 ? 0 : depth.get(sample));
	}
	
	/** Returns minimum contig depth (delegates to depth()) */
	public float minContigDepth() {return depth();}
	/** Returns maximum contig depth (delegates to depth()) */
	public float maxContigDepth() {return depth();}
	
	/** Returns the maximum depth across all samples */
	public final float maxDepth() {
		return depth.max();
	}
	
	/** Returns the sum of depths across all samples */
	public final float depthTotal() {
		return (float)depth.sum();
	}
	
	/**
	 * Returns normalized depth array for multi-sample analysis.
	 * Calculates on first access and caches result.
	 * @return Normalized depth array, null for single-sample bins
	 */
	public float[] normDepth() {
		if(depth.size()<2) {return null;}
		if(normDepth==null) {fillNormDepth();}
		assert(normDepth.length==depth.size()) : normDepth.length+", "+depth.size;
		assert(normDepth.length>1);
		return normDepth;
	}
	
	/** Calculates normalized depth values using log transformation and library normalization.
	 * Normalizes by total library bases, applies log(x+1), then scales to 0-1 range. */
	synchronized void fillNormDepth() {
		if(!FILL_NORM_DEPTH) {return;}
		assert(normDepth==null || (normDepth.length>1 && normDepth.length==numDepths()));
		if(normDepth==null) {normDepth=new float[depth.size];}
		assert(normDepth.length>1) : normDepth.length;
		float sum=0;
		float max=0;
		for(int i=0; i<depth.size; i++) {
			float f=depth.get(i);
			f*=BinObject.invSampleDepthSum[i]; //Now normalized by total library bases
			f=(float)Math.log(f+1);
			sum+=f;
			max=Tools.max(max, f);
			normDepth[i]=f;
		}
		float inv=1/Tools.max(max, 0.1f);
		for(int i=0; i<normDepth.length; i++) {
			normDepth[i]=Tools.mid(0, normDepth[i]*inv, 1);
		}//Now the max should be 1
//		assert(false) : Arrays.toString(normDepth)+", "+Arrays.toString(BinObject.invSampleDepthSum);
	}
	
	/** Uses a weighted sum of linear and geometric means */
	public final float depth() {
		if(depthZeroProxy) {return depth.get(0);}
		if(avgDepthValid) {return avgDepth;}
		synchronized(this) {
			if(depth.size()==1) {avgDepth=depth.get(0);}
			else {
				double product=1;
				double sum=0;
				for(int i=0; i<depth.size; i++) {
					float d=depth.get(i);
					product*=(d+0.25f);
					sum+=d;
				}
				double inv=1.0/depth.size;
				float geo=(float)(Math.pow(product, inv)-0.25);
				float linear=(float)(sum*inv);
				avgDepth=geo*0.75f+linear*0.25f;
			}
			avgDepthValid=true;
		}
		return avgDepth;
	}
	
	@Override
	/** Biggest first */
	public final int compareTo(Sketchable o) {
		if(size()!=o.size()) {return size()>o.size() ? -1 : 1;}//Biggest first
		return id()-o.id();//lowest ID first
	}
	
	@Override
	public final void setFrom(JsonObject all) {
		assert(sketchedSize<size());
		clearTax();
		JsonObject top=null, second=null;
		if(all!=null && all.jmapSize()>0) {
			for(String key : all.jmap.keySet()){
				JsonObject hit=all.jmap.get(key);
				if(top==null) {top=hit;}
				else {
					if(hit.getLong("TaxID")!=1806490) {//Achromobacter sp. ATCC35328; messes with E.coli.
						second=hit;
						break;
					}
				}
			}
		}
		topHit=(top==null ? null : new SketchRecord(top));
		secondHit=(second==null ? null : new SketchRecord(second));
		taxid=(topHit==null ? -1 : topHit.taxid);
		genusTaxid=(topHit==null ? -1 : topHit.genusTaxid);
		sketchedSize=size();
	}
	
	@Override
	public final void clearTax() {
		taxid=genusTaxid=-1;
		topHit=secondHit=null;
		sketchedSize=0;
		clade=null;
		lineage=null;
	}
	
	@Override
	public final String toString() {
		return toBytes().toString();
	}
	
	public final ByteBuilder toSeq() {
		ByteBuilder bb=new ByteBuilder();
		for(Contig c : this) {
			bb.append(c.bases).nl();
		}
		return bb;
	}
	
	/**
	 * Converts bin information to formatted text representation.
	 * Includes size, contigs, GC content, depth, and taxonomic data.
	 * @return ByteBuilder with formatted bin information
	 */
	public final ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		bb.append(isCluster() ? "Cluster " : "Contig ").append(id()).append(":");
		bb.tab().append("Size ").append(size());
		bb.tab().append("Contigs ").append(numContigs());
		bb.tab().append("GC ").append(gc(), 3);
		bb.tab().append("Depth ").append(depth(), 1);
		if(depth.size()>1) {
			for(int i=0; i<depth.size; i++) {bb.comma().append(depth(i), 1);}
		}
		bb.tab().append("TaxID ").append(taxid);
		if(validation) {
			bb.tab().append("TaxID0 ").append(labelTaxid);
			if(completeness>=0) {
				bb.tab().append("Complt ").append(completeness*100, 2);
				bb.tab().append("Contam ").append(contam*100, 2);
			}
		}
//		if(labelTaxid>0) {bb.tab().append("TaxID0 ").append(labelTaxid);}
		if(lineage!=null) {bb.nl().append("Lineage ").append(lineage);}
		if(topHit!=null) {topHit.appendTo(bb.nl().tab().tab());}
		if(secondHit!=null) {secondHit.appendTo(bb.nl().tab().tab());}
		return bb;
	}
	
	/** Higher is more similar */
//	private final float similarityTo(Bin b) {
//		final float ratio=depthRatio(b);
//		final float gc=gc(), gc2=b.gc();
//		final float gcDif=Math.abs(gc-gc2)+1f;
//		final float simDif=SimilarityMeasures.calculateDifferenceAverage(counts, b.counts)*0.5f+1f;
//		final float covariance=1+covariance(b)*32;
//		float product=simDif*ratio*gcDif*covariance;
//		return 1f/product;
//	}
//	
//	/** Higher is more similar */
//	private static final float similarity(float ratio_, float gcDif_, 
//			float simDif_, float covariance_, long edges_) {
//		final float ratio=ratio_;
//		final float gcDif=gcDif_+1f;
//		final float simDif=simDif_*0.5f+1f;
//		final float covariance=1+covariance_*32;
//		float product=simDif*ratio*gcDif*covariance;
//		if(BinObject.verbose) {
//			System.err.println(product+"="+simDif+"*"+ratio+"*"+gcDif+"*"+covariance);
//		}
//		return 1f/product;
//	}
//	
//	public final float similarityTo(Bin b, float stringency) {
//		long size=Tools.min(size(), b.size());
//		float sizeMult=Binner.sizeAdjustMult(size);
//		stringency*=sizeMult;
//		
//		float maxKmerDif=Binner.maxKmerDif2*stringency;
//		float maxDepthRatio=1+((Binner.maxDepthRatio2-1)*stringency);
//		float maxGCDif=Binner.maxGCDif2*stringency;
//		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
//		float maxCovariance=Binner.maxCovariance2*stringency;
//		return similarityTo(b, maxGCDif, maxDepthRatio, maxKmerDif, maxProduct, maxCovariance);
//	}
//	
//	/** Higher is more similar */
//	public final float similarityTo(Bin b, float maxGCDif, float maxDepthRatio, 
//			float maxKmerDif, float maxProduct, float maxCovariance) {
//		long edges1=countEdgesTo(b);
//		long edges2=b.countEdgesTo(this);
//		float mult=(edges1>1 ? 1.4f : 1f)*(edges2>1 ? 1.4f : 1f);
//		
//		if(BinObject.verbose) {
//			System.err.println("Comparing to "+b.id()+": "+
//				"maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+
//				", maxProduct="+maxProduct+", maxGCDif="+maxGCDif+
//				", maxCovariance="+maxCovariance);
//		}
//		
//		float gcDif=Math.abs(gc()-b.gc());
////		float acDif=Math.abs(acRatio-b.acRatio);
////		float eDif=Math.abs(entropy-b.entropy);
////		if(gcDif>maxGCDif*mult || eDif>maxGCDif*mult*4.5f) {return -1;}
//		if(gcDif>maxGCDif*mult) {return -1;}
////		if(acDif>maxGCDif*mult*1f) {return -1;}
//		final float depthRatio=depthRatio(b);
//		final float covariance=covariance(b);
//		if(depthRatio>maxDepthRatio*mult || covariance>maxCovariance*mult) {return -1;}
//		final float kmerDif=SimilarityMeasures.calculateDifferenceAverage(counts, b.counts);
//		final float product=kmerDif*depthRatio;
//		if(kmerDif>maxKmerDif*mult || product>maxProduct*mult) {return -1;}
//		return similarity(depthRatio, gcDif, kmerDif, covariance, Tools.min(edges1, edges2));
//	}
	
	public float covariance(Bin b) {
		if(depth.size()<2) {return 0;}
		
		float[] nda=normDepth();
		float[] ndb=b.normDepth();
		assert(nda!=null) : getClass()+", "+b.getClass();
		assert(ndb!=null) : getClass()+", "+b.getClass();
		assert(nda.length==numDepths()) : getClass()+", "+b.getClass();
		assert(ndb.length==b.numDepths()) : getClass()+", "+b.getClass();
		assert(nda.length==ndb.length) : getClass()+", "+b.getClass();
		float f=SimilarityMeasures.cosineDifference(normDepth(), b.normDepth());
		return f>=0 && Float.isFinite(f) ? f : 0;
	}
	
	/** Returns the size when last sketched for taxonomic annotation */
	public final long sketchedSize() {return sketchedSize;}
	
	/**
	 * Calculates contamination and completeness metrics based on taxonomic labels.
	 * Identifies primary taxon by size, calculates completion vs expected genome size.
	 * @param sizeMap Map of taxon IDs to expected genome sizes
	 */
	public final void calcContam(IntLongHashMap sizeMap) {
		IntLongHashMap taxmap=new IntLongHashMap(7);
		long sum=0;
		for(Contig c : this) {
			int tid=c.labelTaxid;
			taxmap.increment(tid, c.size());
			sum+=c.size();
		}
		assert(sum==size());
		int[] keys=taxmap.keys();
		long[] values=taxmap.values();
		final int invalid=taxmap.invalid();
		int tid=-1;
		long maxSize=-1;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			long value=values[i];
			if(key!=invalid && value>maxSize) {
				tid=key;
				maxSize=value;
			}
		}
		taxid=tid;
		long targetSize=sizeMap.get(tid);
		if(targetSize==sizeMap.invalid()) {targetSize=sum;}//unknown...
		completeness=maxSize/(float)targetSize;
		contam=(sum-maxSize)/(float)sum;
		badContigs=0;
		for(Contig c : this) {
			if(c.labelTaxid>0 && c.labelTaxid!=taxid) {badContigs++;}
		}
	}
	
	/** Alternative contamination calculation focusing on contig-level metrics.
	 * @param sizeMap Map of taxon IDs to expected genome sizes */
	public final void calcContamContigs(IntLongHashMap sizeMap) {
//		taxid=primaryTaxID();
//		long targetSize=sizeMap.get(taxid);
		IntLongHashMap taxmap=new IntLongHashMap(7);
		long sum=0;
		for(Contig c : this) {
			int tid=c.labelTaxid;
			taxmap.increment(tid, c.size());
			sum+=c.size();
		}
		assert(sum==size());
		int[] keys=taxmap.keys();
		long[] values=taxmap.values();
		final int invalid=taxmap.invalid();
		int tid=-1;
		long maxSize=-1;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			long value=values[i];
			if(key!=invalid && value>maxSize) {
				tid=key;
				maxSize=value;
			}
		}
		taxid=tid;
		long targetSize=sizeMap.get(tid);
		if(targetSize==sizeMap.invalid()) {targetSize=sum;}//unknown...
		completeness=maxSize/(float)targetSize;
		contam=(sum-maxSize)/(float)sum;
	}
	
	/** Determines the primary taxonomic ID by finding the taxon with largest total size.
	 * @return Taxonomic ID of the dominant taxon in this bin */
	public final int primaryTaxid() {
		IntLongHashMap taxmap=new IntLongHashMap(7);
		long sum=0;
		for(Contig c : this) {
			int tid=c.labelTaxid;
			taxmap.increment(tid, c.size());
			sum+=c.size();
		}
		assert(sum==size());
		int[] keys=taxmap.keys();
		long[] values=taxmap.values();
		final int invalid=taxmap.invalid();
		int tid=-1;
		long maxSize=-1;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			long value=values[i];
			if(key!=invalid && value>maxSize) {
				tid=key;
				maxSize=value;
			}
		}
		return tid;
	}
	
	/**
	 * Abstract method to determine if this bin is in the same cluster as another.
	 * @param b Bin to compare clustering with
	 * @return true if bins are in the same cluster
	 */
	abstract boolean sameCluster(Bin b);
	
	/** Abstract method returning whether this bin represents a cluster */
	public abstract boolean isCluster();
	
	/** Abstract method returning whether this bin represents a cluster */
	public abstract boolean isContig();
	
	/** Abstract method to convert this bin to a cluster representation */
	public abstract Cluster toCluster();
	
	/** Abstract method returning the cluster containing this bin */
	public abstract Cluster cluster();
	
	/** Abstract method returning whether this bin is valid */
	public abstract boolean isValid();
	
	/** Returns true if this bin contains no contigs */
	public final boolean isEmpty() {return numContigs()<1;}
	
	/**
	 * Counts transitive edges to another bin.
	 * Dispatches to contig or cluster-specific method.
	 * @param b Target bin
	 * @return Number of transitive edges
	 */
	public int transEdgesTo(Bin b) {
		if(!b.isCluster()) {return transEdgesTo((Contig)b);}
		else {return transEdgesTo((Cluster)b);}
	}
	
	/**
	 * Counts transitive edges to a contig through intermediate contigs.
	 * @param b Target contig
	 * @return Maximum transitive edge count
	 */
	public int transEdgesTo(Contig b) {
		if(pairMap==null || Binner.goodTransEdgeMult==1) {return 0;}
		int[] keys=pairMap.keys();//, values=pairMap.values();
		int invalid=pairMap.invalid();
		int v=pairMap.get(b.id());
		int max=0;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			if(key!=invalid) {
				Contig c=DataLoader.allContigs.get(key);
				max=Tools.max(max, c.countEdgesTo(b));
			}
		}
		return Tools.max(v, max);
	}
	
	/**
	 * Counts transitive edges to a cluster through intermediate contigs.
	 * @param b Target cluster
	 * @return Maximum transitive edge count
	 */
	public int transEdgesTo(Cluster b) {
		if(pairMap==null || Binner.goodTransEdgeMult==1) {return 0;}
		int[] keys=pairMap.keys();//, values=pairMap.values();
		int invalid=pairMap.invalid();
		int v=pairMap.get(b.id());
		int max=0;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			if(key!=invalid) {
				Contig c=DataLoader.allContigs.get(key);
				max=Tools.max(max, c.countEdgesTo(b));
			}
		}
		return Tools.max(v, max);
	}
	
	/**
	 * Counts direct edges to another bin.
	 * Dispatches to contig or cluster-specific method.
	 * @param b Target bin
	 * @return Number of direct edges
	 */
	public int countEdgesTo(Bin b) {
		if(!b.isCluster()) {return countEdgesTo((Contig)b);}
		else {return countEdgesTo((Cluster)b);}
	}
	
	/**
	 * Counts direct edges to a specific contig.
	 * @param b Target contig
	 * @return Number of edges to the contig
	 */
	public int countEdgesTo(Contig b) {
		return pairMap==null ? 0 : Tools.max(0, pairMap.get(b.id()));
	}
	
	/**
	 * Counts maximum edges to any contig in a cluster.
	 * @param b Target cluster
	 * @return Maximum edge count to cluster members
	 */
	public int countEdgesTo(Cluster b) {
		if(pairMap==null) {return 0;}
		final int[] keys=pairMap.keys(), values=pairMap.values();
		final int invalid=pairMap.invalid();
		long sum=0;
		int max=0;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			if(key!=invalid && b.contigSet.contains(key)) {
				int v=values[i];
				sum+=v;
				max=Tools.max(max, v);
			}
		}
		return max;
	}
	
	/**
	 * Counts reciprocal edges between this bin and another.
	 * Dispatches to contig or cluster-specific method.
	 * @param b Target bin
	 * @return Number of reciprocal edges
	 */
	public int countReciprocalEdges(Bin b) {
		if(!b.isCluster()) {return countReciprocalEdges((Contig)b);}
		else {return countReciprocalEdges((Cluster)b);}
	}
	
	/**
	 * Counts reciprocal edges to a contig (minimum of both directions).
	 * @param b Target contig
	 * @return Minimum edge count in both directions
	 */
	public int countReciprocalEdges(Contig b) {
		if(pairMap==null || b.pairMap==null) {return 0;}
		return Tools.max(0, Tools.min(pairMap.get(b.id()), b.countEdgesTo(this)));
//		final int[] keysB=b.pairMap.keys();
//		final int invalidB=b.pairMap.invalid();
//		for(int i=0; i<keysB.length; i++) {
//			int key=keysB[i];
//			if(key!=invalid && b.contigSet.contains(key)) {
//				int v=values[i];
//				sum+=v;
//				max=Tools.max(max, v);
//			}
//		}
//		
//		return pairMap==null ? 0 : Tools.max(0, pairMap.get(b.id()));
	}
	
	/**
	 * Counts reciprocal edges to a cluster (slower operation).
	 * Optimizes by processing smaller cluster first.
	 * @param b Target cluster
	 * @return Maximum reciprocal edge count
	 */
	public int countReciprocalEdges(Cluster b) {//Slow
		if(pairMap==null || b.pairMap==null) {return 0;}
		if(!this.isCluster()) {return b.countReciprocalEdges((Contig)this);}
		if(size()>b.size()) {return b.countReciprocalEdges((Cluster)this);}
		//At this point they are both clusters and this one is smaller

		//		final int[] keys=pairMap.keys(), values=pairMap.values();
		//		final int invalid=pairMap.invalid();
		//		int max=0;
		//		for(int i=0; i<keys.length; i++) {
		//			int key=keys[i];
		//			if(key!=invalid && b.contigSet.contains(key)) {
		//				int v1=values[i], v2=b.pairMap.get(key)
		//				max=Tools.max(max, v);
		//			}
		//		}

		int max=0;
		for(Contig c : ((Cluster)this).contigs) {
			max=Tools.max(max, b.countReciprocalEdges(c));
		}

		return max;
	}
	
	/**
	 * Converts this bin to a Clade representation.
	 * Adds all contig sequences to the clade and finishes construction.
	 * @return Clade object containing this bin's sequences
	 */
	Clade toClade() {
		clade=new Clade(-1, -1, name());
		for(Contig c : this) {
			clade.add(c.bases, null);
		}
		clade.finish();
		return clade;
	}
	
	/** Returns the raw depth list for this bin */
	public FloatList depthList() {
		return depth;
	}
	
	/** Returns true if this bin contains 16S or 18S rRNA sequences */
	final boolean hasSSU() {return r16S!=null || r18S!=null;}

	/** Number of tetramers in k-mer spectrum */
	public int numTetramers;
	/** Number of pentamers in k-mer spectrum */
	public int numPentamers;
//	public float invKmers;
	
	/** Dimer frequency counts */
	public int[] dimers;
	/** Trimer frequency counts */
	public int[] trimers;
	/** Tetramer frequency counts */
	public int[] tetramers;
	/** Pentamer frequency counts */
	public int[] pentamers;
	/** Total GC bases across all contigs in this bin */
	public long gcSum;
	/** Size of bin when last sketched for taxonomic annotation */
	public long sketchedSize;
//	public float acRatio;
	
	/** Sum of depth values (initialized to 0) */
	public long depthSum=0;
	/** Cached average depth value */
	private float avgDepth=-1;
	/** Flag indicating if cached average depth is valid */
	boolean avgDepthValid=false;
	/** List of depth values across samples */
	private FloatList depth=new FloatList(1);
	/** Cached normalized depth array */
	private float[] normDepth;
	/** Map of contig pairs and their edge counts */
	public IntHashMap2 pairMap;
	/** Estimated contamination level (0-1) */
	/** Estimated genome completeness (0-1) */
	float completeness=0, contam=0;
	/** Number of contigs with incorrect taxonomic labels */
	int badContigs=0;
	/** Sequence entropy measure */
	float entropy;
	/** Measure of strand bias */
	float strandedness;
	/** H-H di-amino acid frequency metric */
	public float hh;
	/** CAGA motif frequency metric */
	public float caga;
	/** Overall quality score for this bin */
	float score;
	/** Flag indicating if bin was reclustered */
	public boolean wasReclustered=false;
	
	/** Destination cluster ID (-1 if unassigned) */
	int dest=-1;
	
	/** Primary taxonomic ID for this bin */
	public int taxid;
	/** Genus-level taxonomic ID */
	public int genusTaxid;
	/** True taxonomic ID for validation on labeled datasets */
	public int labelTaxid;//For validation on labeled data
	/** Best taxonomic match from sketching */
	SketchRecord topHit;
	/** Second best taxonomic match from sketching */
	SketchRecord secondHit;
	/** Full taxonomic lineage string */
	public String lineage;
	/** Clade object for taxonomic classification */
	public Clade clade;
	/** 16S rRNA sequence if found */
	public byte[] r16S;
	/** 18S rRNA sequence if found */
	public byte[] r18S;
	
}
