package bin;

import clade.Clade;
import structures.ByteBuilder;

/**
 * Represents statistical metadata and quality classification for genomic bins.
 * Captures bin characteristics like size, contamination, completeness, and RNA marker gene counts
 * for metagenomic assembly evaluation and quality assessment.
 * @author Brian Bushnell
 */
public class BinStats implements Comparable<BinStats> {
	
	/**
	 * Constructs BinStats from a Bin object and name.
	 * Extracts bin metadata: size, contamination, completeness, GC, and depth.
	 * @param b The bin to extract statistics from
	 * @param name_ Identifier for this bin
	 */
	BinStats(Bin b, String name_){
		name=name_;
		id=b.id();
		taxid=b.taxid;
		if(taxid<1) {taxid=b.labelTaxid;}
		size=b.size();
		contigs=b.numContigs();
		contam=b.contam;
		complt=b.completeness;
		badContigs=b.badContigs;
		gc=b.gc();
		depth=b.depth();
		minDepth=b.minContigDepth();
		maxDepth=b.maxContigDepth();
		lineage=b.lineage;
		clade=b.clade;

		assert(b.gc()!=0) : this;
	}
	
	@Override
	public int compareTo(BinStats b) {
		if(size!=b.size) {return size>b.size ? -1 : 1;}
		return name.compareTo(b.name);
	}
	
	/**
	 * Determines bin quality type using completeness, contamination, and RNA markers.
	 * @param useRNA Whether RNA markers are required for high-quality classification
	 * @return Quality type (UHQ, VHQ, HQ, MQ, LQ, VLQ, or HCN)
	 */
	String type(boolean useRNA) {
		return type(complt, contam, r16Scount, r23Scount, r5Scount, trnaCount, useRNA);
	}

	/**
	 * Checks if this bin qualifies as high quality (UHQ, VHQ, or HQ).
	 * @param useRNA Whether to require RNA marker genes for classification
	 * @return true if bin type ends with "HQ"
	 */
	public boolean hq(boolean useRNA) {
		String type=type(useRNA);
		return type.endsWith("HQ");
	}
	/**
	 * Checks if this bin qualifies as medium quality (MQ).
	 * @param useRNA Whether to require RNA marker genes for classification
	 * @return true if bin type equals "MQ"
	 */
	public boolean mq(boolean useRNA) {
		String type=type(useRNA);
		return type.equals("MQ");
	}
	/**
	 * Checks if this bin qualifies as medium quality (LQ).
	 * @param useRNA Whether to require RNA marker genes for classification
	 * @return true if bin type ends with "LQ"
	 */
	public boolean lq() {
		String type=type(false);
		return type.endsWith("LQ");
	}
	
	/**
	 * Classifies bin quality based on completeness, contamination, and RNA marker gene counts.
	 * Quality tiers: UHQ (>=99% complete, <=1% contamination), VHQ (>=95% complete, <=2% contamination),
	 * HQ (>90% complete, <5% contamination), MQ (>=50% complete, <10% contamination),
	 * LQ (<50%, <10%)
	 * VLQ (<20%, <10%)
	 * HCN (>=10% contamination).
	 * High quality tiers (UHQ, VHQ, HQ) require RNA markers if useRNA is true.
	 *
	 * @param complt Completeness fraction (0.0-1.0)
	 * @param contam Contamination fraction (0.0-1.0)
	 * @param r16S Count of 16S ribosomal RNA genes
	 * @param r23S Count of 23S ribosomal RNA genes
	 * @param r5S Count of 5S ribosomal RNA genes
	 * @param trna Count of tRNA genes
	 * @param useRNA Whether RNA markers are required for high quality classification
	 * @return Quality type string (UHQ, VHQ, HQ, MQ, LQ, VLQ, or HCN)
	 */
	static String type(float complt, float contam, int r16S, int r23S, int r5S, int trna, boolean useRNA) {
		boolean rnaOK=!useRNA || (r16S>0 && r23S>0 && r5S>0 && trna>=18);
		if(contam<0.05 && complt>0.90 && rnaOK) {//HQ
			if(contam<=0.01 && complt>=0.99) {return "UHQ";}
			if(contam<=0.02 && complt>=0.95) {return "VHQ";}
			return "HQ";
		}
		if(contam<0.10 && complt>=0.50) {return "MQ";}
		if(contam<0.10 && complt<0.50) {//LQ
			if(contam<0.10 && complt<0.20) {return "VLQ";}
			return "LQ";
		}
		assert(contam>=0.10);
		return "HCN";
	}
	
	public String toString() {return toBytes(null).toString();}
	
	public ByteBuilder toBytes(ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		bb.appendln(name);
		bb.appendln(id);
		bb.appendln(taxid);
		bb.appendln(size);
		bb.appendln(contigs);
		bb.appendln(badContigs);
		bb.append(contam, 4).nl();
		bb.append(complt, 4).nl();
		bb.append(gc, 4).nl();
		bb.append(depth, 4).nl();
		bb.append(minDepth, 4).nl();
		bb.append(maxDepth, 4).nl();
		bb.appendln(lineage);
		return bb;
	}
	
	/** Bin name identifier */
	final String name;
	/** Numeric bin identifier */
	int id;
	/** Taxonomic ID assigned to this bin */
	int taxid;
	/** Total size of bin in base pairs */
	long size;
	/** Number of contigs in this bin */
	int contigs;
	/** Number of contaminated or low-quality contigs in this bin */
	int badContigs;
	/** Contamination level as fraction (0.0-1.0) */
	float contam;
	/** Completeness level as fraction (0.0-1.0) */
	float complt;
	/** GC content as fraction (0.0-1.0) */
	float gc;
	/** Average sequencing depth across the bin */
	float depth;
	float minDepth, maxDepth;

	/** Count of 5S ribosomal RNA genes detected */
	int r5Scount=0;
	/** Count of 16S ribosomal RNA genes detected */
	int r16Scount=0;
	/** Count of 18S ribosomal RNA genes detected */
	int r18Scount=0;
	/** Count of 23S ribosomal RNA genes detected */
	int r23Scount=0;
	/** Count of transfer RNA genes detected */
	int trnaCount=0;
	/** Count of coding sequences detected */
	int cdsCount=0;
	/** Total length of coding sequences in base pairs */
	long cdsLength=0;
	/** Clade object associated with this bin's taxonomic classification */
	Clade clade;
	/** Taxonomic lineage string for this bin */
	String lineage;
	/** Short filename for this bin (e.g. "bin_0.fa"), set when loaded from files */
	String filename;
	/** Name of the first contig in this bin */
	String contigName;

}
