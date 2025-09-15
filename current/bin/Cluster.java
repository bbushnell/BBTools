package bin;

import java.util.ArrayList;
import java.util.Iterator;

import shared.Tools;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.IntHashMap;
import structures.IntHashSet;
import tracker.EntropyTracker;

public class Cluster extends Bin {

	public Cluster(int id_) {id=id_;}
	
	public Cluster(Contig c) {
		id=c.id();
		add(c);
	}
	
	@Override
	public String name() {return contigs.get(0).name;}
	
	@Override
	public boolean isCluster() {return true;}
	
	@Override
	public Cluster toCluster() {return this;}
	
	@Override
	public Cluster cluster() {return this;}
	
	@Override
	public final int id() {return id;}

	public boolean pure(float fraction) {return nonDominantSize<=size*fraction;}
	
	public float minContigDepth() {
		float d=contigs.size()==0 ? 0 : contigs.get(0).depth();
		for(Contig c : contigs) {d=Math.min(c.depth(), d);}
		return d;
	}
	
	public float maxContigDepth() {
		float d=0;
		for(Contig c : contigs) {d=Math.max(c.depth(), d);}
		return d;
	}
	
	@Override
	@Deprecated
	public void setID(int id_) {
		throw new RuntimeException("Not supported.");
	}
	
	public Cluster add(Bin b) {
		if(b.getClass()==Cluster.class) {return add((Cluster)b);}
		else {return add((Contig)b);}
	}
	
	public Cluster add(Cluster clust) {
		assert(clust!=this) : id();
		for(Contig c : clust.contigs) {add(c);}
		clust.contigs=null;
//		assert(isValid());
		return this;
	}
	
	boolean sameCluster(Bin b) {
		if(contigSet.contains(b.id())) {return true;}
		return b.isCluster() && ((Cluster)b).contigSet.contains(id());
	}
	
	public synchronized Cluster add(Contig c) {
		timestamp=globalTime;
//		verbose=(id()==52 || id()==3 || id()==3 || c.id()==52 || c.id()==3 || c.id()==3);
//		if(verbose) {
//			System.err.println("Adding "+c.id()+" to "+this.id()+(c.cluster==null ? "" : " from "+c.cluster.id));
//			assert(c.cluster==null || c.cluster.id()!=3);
//		}
//		assert(c.isValid());
		assert(c.cluster!=this);
		assert(!contigSet.contains(c.id()));
		c.cluster=this;
		final long size2=size()+c.size();
		final float invSize2=1f/size2;
		if(contigs.size()==0) {
			topHit=c.topHit;
			secondHit=c.secondHit;
			taxid=c.taxid;
			genusTaxid=c.genusTaxid;
			labelTaxid=c.labelTaxid;
			sketchedSize=c.sketchedSize;
			entropy=c.entropy;
//			acRatio=c.acRatio;
		}else {
			if(labelTaxid<1 || (c.labelTaxid!=labelTaxid && c.size()>=0.05f*size())) {nonDominantSize+=c.size();}
//			assert(isValid());
			entropy=(entropy*size()+c.entropy*c.size())*invSize2;
//			acRatio=(acRatio*size()+c.acRatio*c.size())*invSize2;
		}
		contigs.add(c);
		contigSet.add(c.id());
		numTetramers+=c.numTetramers;
		numPentamers+=c.numPentamers;
//		invKmers=1f/numTetramers;
		if(tetramers==null) {tetramers=c.tetramers.clone();}
		else {Tools.add(tetramers, c.tetramers);}
		if(c.dimers!=null) {
			if(dimers==null) {dimers=c.dimers.clone();}
			else {Tools.add(dimers, c.dimers);}
			strandedness=EntropyTracker.strandedness(dimers, 2);
		}
		if(c.trimers!=null) {
			if(trimers==null) {trimers=c.trimers.clone();}
			else {Tools.add(trimers, c.trimers);}
		}
		if(c.pentamers!=null) {
			if(pentamers==null) {pentamers=c.pentamers.clone();}
			else {Tools.add(pentamers, c.pentamers);}
		}
		
		if(c.pairMap!=null) {
			if(pairMap==null) {pairMap=new IntHashMap(5);}
//			pairMap.incrementAll(c.pairMap);//Remember, the targets are contig IDs.
			pairMap.setToMax(c.pairMap);
		}
		
		if(numDepths()==0) {//Empty cluster
			assert(size==0 || grading);
			for(int i=0, max=c.numDepths(); i<max; i++) {appendDepth(c.depth(i));}
		}else {//Nonempty cluster
			assert(c.numDepths()==numDepths()) : c.numDepths()+", "+numDepths();
			float cmult=invSize2*c.size();
			float bmult=invSize2*size;
			for(int i=0, max=c.numDepths(); i<max; i++) {
				float cdepth=c.depth(i)*cmult;
				float bdepth=depth(i)*bmult;
				setDepth(cdepth+bdepth, i);
			}
		}
		if(numDepths()>1) {
			fillNormDepth();
		}
		avgDepthValid=false;
		size=size2;
		gcSum+=c.gcSum;

		if(r16S==null) {r16S=c.r16S;}
		if(r18S==null) {r18S=c.r18S;}
		
		return this;
	}
	
	@Override
	public long size() {return size;}

	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		String name=Long.toString(id());
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		for(Contig c : contigs) {
			r.bases=c.bases;
			smm.processReadNucleotide(r);
		}
		return smm.toSketch(0);
	}
	
	@Override
	public boolean isValid() {
		if(id()<0) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(numDepths()<1) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(tetramers==null) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(contigs.isEmpty()) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		for(Contig c : contigs) {
			if(c.cluster!=this) {
				assert(c.cluster!=null);
				assert(false) : "Cluster "+id()+" contains a contig "+c.id()+" that points to "+c.cluster().id()+"\n"
						+ " set="+contigSet+", c.set="+c.cluster.contigSet+"\n"+contigs+"\n"+c.cluster.contigs+"\n";
				return false;
			}
			if(!contigSet.contains(c.id())) {
				assert(false) : id()+", "+contigSet;
				return false;
			}
//			if(c.pairMap!=null && !c.pairMap.isEmpty()) {
//				for(KeyValue kv : KeyValue.toList(c.pairMap)) {
//					int value2=pairMap.get(kv.key);
//					assert(value2>=kv.value) : c.pairMap+"\n"+pairMap;
//				}
//			}
		}
		return true;
	}
	
	public void clear() {
		size=0;
		contigs.clear();
		contigSet.clear();
		pairMap=null;
		
		numTetramers=0;
		numPentamers=0;
//		invKmers=0;
		dimers=null;
		trimers=null;
		tetramers=null;
		pentamers=null;
		gcSum=0;
		sketchedSize=0;
		clearDepth();
		completeness=contam=entropy=strandedness=score=0;
		dest=0;
		taxid=genusTaxid=labelTaxid=0;
		topHit=secondHit=null;
		nonDominantSize=0;
		r16S=r18S=null;
		wasReclustered=false;
	}
	
	@Override
	public int numContigs() {return contigs==null ? 0 : contigs.size();}
	
	@Override
	public Iterator<Contig> iterator() {
		return contigs.iterator();
	}
	
	final int id;
	public long size;
	public ArrayList<Contig> contigs=new ArrayList<Contig>(8);
	public IntHashSet contigSet=new IntHashSet(5);
	private int timestamp=-1;
	long nonDominantSize=0;
	
}
