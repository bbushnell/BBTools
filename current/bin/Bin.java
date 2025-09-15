package bin;

import json.JsonObject;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntLongHashMap;

public abstract class Bin extends BinObject implements Sketchable, Iterable<Contig> {
	
	@Override
	public final int taxid() {return taxid;}
	
	abstract String name();

	@Override
	public final float gc() {return gcSum/(float)size();}
	
	public final void clearDepth() {
		depth.clear();
		normDepth=null;
	}
	
	public float depthRatio(Bin b, int method) {
		if(method==0) {return 0;}
		if(method==1) {return depthRatio(b);}
		if(method==2) {return depthRatio2(b);}
		if(method==3) {return depthRatio3(b);}
		if(method==4) {return depthRatio4(b);}
		if(method==5) {return 1-depthRatio2(b);}
		if(method==6) {return 1-depthRatio3(b);}
		if(method==7) {return 1-depthRatio4(b);}
		throw new RuntimeException(""+method);
	}
	
	public float depthRatio(Bin b) {
		float max=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=Tools.max(d1,d2)/Tools.min(d1,d2);
			max=Tools.max(max, ratio);
		}
		return max;
	}
	
	public float depthRatio2(Bin b) {
		float min=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=2*Tools.min(d1,d2)/(d1+d2);
			min=Tools.min(min, ratio);
		}
		return min;
	}
	
	public float depthRatio3(Bin b) {
		float max=0;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=2*(Tools.max(d1,d2)/(d1+d2))-1;
			max=Tools.max(max, ratio);
		}
		return max;
	}
	
	public float depthRatio4(Bin b) {
		float min=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+depthBoost;
			float d2=b.depth.get(i)+depthBoost;
			float ratio=Tools.min(d1,d2)/Tools.max(d1,d2);
			min=Tools.min(min, ratio);
		}
		return min;
	}
	
	public boolean pure() {return pure(0.02f);}
	public boolean pure(float fraction) {return labelTaxid>0;}
	
	public final void setDepth(float d, int sample) {
		depth.set(sample, d);
	}
	
	public final void appendDepth(float d) {
		depth.add(d);
	}
	
	public final int numDepths() {
		return depth.size;
	}
	
	public final float numEdges() {
		return pairMap==null ? 0 : pairMap.size();
	}
	
	public final float depth(int sample) {
		return (sample==1 && depth.size==1 ? 0 : depth.get(sample));
	}
	
	public float minContigDepth() {return depth();}
	public float maxContigDepth() {return depth();}
	
	public final float maxDepth() {
		return depth.max();
	}
	
	public final float depthTotal() {
		return (float)depth.sum();
	}
	
	public float[] normDepth() {
		if(depth.size()<2) {return null;}
		if(normDepth==null) {fillNormDepth();}
		assert(normDepth.length==depth.size()) : normDepth.length+", "+depth.size;
		assert(normDepth.length>1);
		return normDepth;
	}
	
	synchronized void fillNormDepth() {
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
//		assert(false) : Arrays.toString(normDepth);
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
		return o.id()-id();
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
	}
	
	@Override
	public final String toString() {
		return toBytes().toString();
	}
	
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
	
	public final long sketchedSize() {return sketchedSize;}
	
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
	
	abstract boolean sameCluster(Bin b);
	
	public abstract boolean isCluster();
	
	public abstract Cluster toCluster();
	
	public abstract Cluster cluster();
	
	public abstract boolean isValid();
	
	public final boolean isEmpty() {return numContigs()<1;}
	
	public int transEdgesTo(Bin b) {
		if(!b.isCluster()) {return transEdgesTo((Contig)b);}
		else {return transEdgesTo((Cluster)b);}
	}
	
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
	
	public int countEdgesTo(Bin b) {
		if(!b.isCluster()) {return countEdgesTo((Contig)b);}
		else {return countEdgesTo((Cluster)b);}
	}
	
	public int countEdgesTo(Contig b) {
		return pairMap==null ? 0 : Tools.max(0, pairMap.get(b.id()));
	}
	
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
	
	public int countReciprocalEdges(Bin b) {
		if(!b.isCluster()) {return countReciprocalEdges((Contig)b);}
		else {return countReciprocalEdges((Cluster)b);}
	}
	
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
	
	public FloatList depthList() {
		return depth;
	}
	
	final boolean hasSSU() {return r16S!=null || r18S!=null;}

	public int numTetramers;
	public int numPentamers;
//	public float invKmers;
	
	public int[] dimers;
	public int[] trimers;
	public int[] tetramers;
	public int[] pentamers;
	public long gcSum;
	public long sketchedSize;
//	public float acRatio;
	
	public long depthSum=0;
	private float avgDepth=-1;
	boolean avgDepthValid=false;
	private FloatList depth=new FloatList(1);
	private float[] normDepth;
	public IntHashMap pairMap;
	float completeness=0, contam=0;
	int badContigs=0;
	float entropy;
	float strandedness;
	float score;
	public boolean wasReclustered=false;
	
	int dest=-1;
	
	public int taxid;
	public int genusTaxid;
	public int labelTaxid;//For validation on labeled data
	SketchRecord topHit;
	SketchRecord secondHit;
	public byte[] r16S;
	public byte[] r18S;
	
}
