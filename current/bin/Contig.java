package bin;

import java.util.ArrayList;
import java.util.Iterator;

import dna.AminoAcid;
import shared.Shared;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.ByteBuilder;

/**
 * Represents a genomic contig with sequence and metadata for bioinformatics analysis.
 * Manages contig sequence data, provides iterator, and supports sketch generation.
 * @author Brian Bushnell
 * @date 2013
 */
public class Contig extends Bin {

	/**
	 * Constructs a new Contig with specified name, sequence bases, and identifier.
	 * @param name_ Contig name/identifier
	 * @param bases_ Genomic sequence as byte array
	 * @param id_ Unique numeric identifier for this contig
	 */
	public Contig(String name_, byte[] bases_, int id_) {
		name=name_;
		shortName=ContigRenamer.toShortName(name);
		bases=bases_;
		id=id_;
	}

//	public Contig(String name_, byte[] bases_, int id_) {
//		name=name_;
//		shortName=ContigRenamer.toShortName(name);
//		bases=bases_;
//		setID(id_);
//	}
	
	@Override
	public String name() {return name;}
	
	@Override
	public boolean isCluster() {return false;}
	
	@Override
	public Cluster toCluster() {
		assert(cluster==null);
		return new Cluster(this);
	}
	
	/**
	 * Checks if this contig belongs to the same cluster as the specified bin.
	 * @param b Bin to compare cluster membership with
	 * @return true if both belong to the same cluster
	 */
	boolean sameCluster(Bin b) {
		return cluster!=null && cluster.contigSet.contains(b.id());
	}
	
//	@Override
//	public int clusterID() {return clusterID;}
	
	@Override
	public Cluster cluster() {return cluster;}
	
	@Override
	public void setID(int id_) {
//		assert(id==-1 && id_>-1);
		id=id_;
	}
	
	@Override
	public int id() {
//		assert(id>=0);
		return id;
	}
	
	@Override
	public boolean isValid() {
		if(numDepths()<1) {
			assert(false) : numDepths()+", "+name+"\n"+this;
			return false;
		}
		if(tetramers==null) {
			assert(false) : name;
			return false;
		}
//		assert(gcSum>0) : gcSum+", "+new String(bases);
		if(cluster!=null) {
			if(!cluster.contigSet.contains(id())) {
				assert(false) : id()+", "+cluster.id+", "+cluster.contigSet;
				return false;
			}
//			if(pairMap!=null && !pairMap.isEmpty()) {
//				for(KeyValue kv : KeyValue.toList(pairMap)) {
//					int value2=cluster.pairMap.get(kv.key);
//					assert(value2>=kv.value) : pairMap+"\n"+cluster.pairMap;
//				}
//			}
		}
		return true;
	}
	
	/**
	 * Loads k-mer frequency counts and GC content from the sequence.
	 * Calculates tetramers, optionally trimers and pentamers based on configuration.
	 * Also computes GC sum by counting G and C nucleotides in the sequence.
	 */
	public void loadCounts() {
		assert(numTetramers==0);
		tetramers=new int[canonicalKmers[4]];
		numTetramers=countKmers(bases, tetramers, 4);
//		invKmers=1f/Tools.max(1, numTetramers);
		for(byte b : bases) {
			int x=AminoAcid.baseToNumber[b];
			gcSum+=(x==1 || x==2) ? 1 : 0;
		}

		if(countTrimers) {
			trimers=new int[canonicalKmers[3]];
			countKmers(bases, trimers, 3);
		}

		if(countPentamers && size()>=minPentamerSizeCount) {
			pentamers=new int[canonicalKmers[5]];
			numPentamers=countKmers(bases, pentamers, 5);
		}
	}
	
	/** In fasta format */
	public void appendTo(ByteBuilder bb, int cluster) {
		bb.append('>').append(name);
		if(cluster>=0) {bb.tab().append("cluster_").append(cluster);}
		bb.nl();
		final int wrap=Shared.FASTA_WRAP;
		for(int i=0; i<bases.length; i+=wrap) {
			//Now with modified append I can just append(bases, wrap) ...???
			bb.append(bases, i, wrap).nl();
		}
	}
	
	@Override
	public long size() {return bases.length;}

	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		String name=Long.toString(id());
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		r.bases=bases;
		smm.processReadNucleotide(r);
		return smm.toSketch(0);
	}
	
	/**
	 * Generates coverage output format containing contig metadata and depth values.
	 * Includes name, ID, size, depth values for all samples, and edge information.
	 * @param bb ByteBuilder to use, created if null
	 * @return ByteBuilder containing formatted coverage data
	 */
	public final ByteBuilder toCov(ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		bb.append(shortName);
		bb.tab().append(id());
		bb.tab().append(size());
		for(int i=0, max=numDepths(); i<max; i++) {
			bb.tab().append(depth(i), 2);
		}
		ArrayList<KeyValue> list=KeyValue.toList(pairMap);
		if(list!=null) {
			for(int i=0; i<list.size() && i<DataLoader.MAX_EDGES_TO_PRINT; i++) {
				KeyValue ip=list.get(i);
				bb.tab().append(ip.key).tab().append(ip.value);
			}
		}
		return bb;
	}
	
	@Override
	public int numContigs() {return 1;}
	
	@Override
	public Iterator<Contig> iterator() {
		return new ContigIterator();
	}
	
	/** Iterator implementation that returns this single contig exactly once. */
	private class ContigIterator implements Iterator<Contig> {

		@Override
		public boolean hasNext() {
			return hasMore;
		}

		@Override
		public Contig next() {
			if(!hasMore) {return null;}
			hasMore=false;
			return Contig.this;
		}
		
		/** Flag tracking whether the iterator has more elements to return */
		boolean hasMore=true;
		
	}
	
	/** Numeric identifier for this contig */
	private int id=-1;
	/** The cluster this contig belongs to, or null if not clustered */
	public Cluster cluster=null;
	/** Full name/identifier of this contig */
	public final String name;
	/** Shortened version of the contig name for display purposes */
	public final String shortName;
	/** Genomic sequence data as byte array */
	public final byte[] bases;
	
}
