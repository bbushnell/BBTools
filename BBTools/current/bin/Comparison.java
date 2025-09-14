package bin;

import aligner.IDAligner;
import aligner.SingleStateAlignerFlat2;
import shared.Vector;
import structures.ByteBuilder;
import tax.TaxTree;

/**
 * Compares two Clade objects to determine their genomic similarity.
 * Implements multiple comparison methods including absolute difference, 
 * cosine similarity, Hellinger distance, Euclidean distance, and 
 * GC-compensated absolute difference.
 * 
 * @author Brian Bushnell
 * @date April 19, 2024
 */
public class Comparison extends BinObject implements Comparable<Comparison> {
	
	/**
	 * Creates an empty Comparison object.
	 */
	public Comparison() {}
	
	/**
	 * Creates a Comparison between query and reference Clades using default parameters.
	 * 
	 * @param query_ The query Clade
	 * @param ref_ The reference Clade
	 */
	public Comparison(Clade query_, Clade ref_) {
		compare(query_, ref_, 1, 1, 1);
	}
	
	/**
	 * Compares two Clades using a two-stage process.
	 * First performs a quick filtering based on GC content and strandedness,
	 * then performs a more detailed k-mer based comparison if the quick filter passes.
	 * 
	 * @param q The query Clade
	 * @param r The reference Clade
	 * @param gcLimit Maximum allowed GC content difference
	 * @param strLimit Maximum allowed strandedness difference
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated similarity measure (lower is more similar)
	 */
	float compare(Clade q, Clade r, float gcLimit, float strLimit, float k5Limit) {
		boolean pass=quickCompare(q, r, gcLimit, strLimit);
		if(!pass) {return 4+gcdif;}
		return slowCompare(q, r, k5Limit);
	}
	
	/**
	 * Performs a rapid preliminary comparison based on GC content and strandedness.
	 * Sets the query and ref fields and calculates basic difference metrics.
	 * 
	 * @param q The query Clade
	 * @param r The reference Clade
	 * @param gcLimit Maximum allowed GC content difference
	 * @param strLimit Maximum allowed strandedness difference
	 * @return true if the Clades pass the quick comparison filter
	 */
	boolean quickCompare(Clade q, Clade r, float gcLimit, float strLimit) {
		query=q;
		ref=r;
		assert(query.finished());
		assert(ref.finished());
		clearDif();

		gcdif=Math.abs(query.gc-ref.gc);
		strdif=Math.abs(query.strandedness-ref.strandedness);
		entdif=Math.abs(query.gcCompEntropy-ref.gcCompEntropy);
		return gcdif<=gcLimit && strdif<=strLimit;
	}
	
	/**
	 * Performs a detailed comparison using k-mer frequency profiles.
	 * The specific comparison method is determined by the static method field.
	 * 
	 * @param q The query Clade
	 * @param r The reference Clade
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated similarity measure (lower is more similar)
	 */
	float slowCompare(Clade q, Clade r, float k5Limit) {
		assert(query==q && ref==r);
		if(method==ABS) {return compareABS(k5Limit);}
		if(method==ABSCOMP) {return compareABSCOMP(k5Limit);}
		if(method==COS) {return compareCOS(k5Limit);}
		if(method==HEL) {return compareHEL(k5Limit);}
		if(method==EUC) {return compareEUC(k5Limit);}
		assert(false) : "Bad method: "+method;
		return 1;
	}
	
	/**
	 * Compares Clades using absolute difference between k-mer frequencies.
	 * Uses early exit optimizations to avoid unnecessary calculations when possible.
	 * 
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated absolute difference (lower is more similar)
	 */
	private float compareABS(float k5Limit) {
//		k3dif=SimilarityMeasures.absDif(query.counts[3], ref.counts[3]);
		k3dif=Vector.absDifFloat(query.trimers, ref.trimers);
		if(earlyExit && k3dif*comparisonCutoffMult2>k5Limit) {return k3dif*4;}
		k4dif=maxK<4 ? k3dif : SimilarityMeasures.absDif(query.counts[4], ref.counts[4]);
		if(earlyExit && k4dif*comparisonCutoffMult>k5Limit) {return k4dif*2;}
		k5dif=maxK<5 ? k4dif : SimilarityMeasures.absDif(query.counts[5], ref.counts[5]);
		return k5dif;
	}
	
	/**
	 * Compares Clades using GC-compensated absolute difference between k-mer frequencies.
	 * K-mers are grouped by GC content before comparison to reduce GC bias.
	 * Uses early exit optimizations to avoid unnecessary calculations when possible.
	 * 
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated GC-compensated absolute difference (lower is more similar)
	 */
	private float compareABSCOMP(float k5Limit) {
//		k3dif=SimilarityMeasures.absDifComp(query.counts[3], ref.counts[3], 3);
		k3dif=Vector.absDifFloat(query.trimers, ref.trimers);//Already compensated
		if(earlyExit && k3dif*comparisonCutoffMult2>k5Limit) {return k3dif*4;}
		k4dif=maxK<4 ? k3dif : SimilarityMeasures.absDifComp(query.counts[4], ref.counts[4], 4);
//		k4dif=maxK<4 ? k3dif : Vector.absDifComp(query.counts[4], ref.counts[4], 4, BinObject.gcmapMatrix[4]);
//		k4dif=Vector.absDifFloat(query.tetramers, ref.tetramers);//Already compensated.  This makes it 12% slower for some reason.
		if(earlyExit && k4dif*comparisonCutoffMult>k5Limit) {return k4dif*2;}
		k5dif=maxK<5 ? k4dif : SimilarityMeasures.absDifComp(query.counts[5], ref.counts[5], 5);
//		k5dif=maxK<5 ? k4dif : Vector.absDifComp(query.counts[5], ref.counts[5], 5, BinObject.gcmapMatrix[5]);
		return k5dif;
	}
	
	/**
	 * Compares Clades using cosine distance between k-mer frequency vectors.
	 * Uses early exit optimizations to avoid unnecessary calculations when possible.
	 * 
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated cosine distance (lower is more similar)
	 */
	private float compareCOS(float k5Limit) {
		k3dif=SimilarityMeasures.cosineDifference(query.counts[3], ref.counts[3]);
		if(earlyExit && k3dif*comparisonCutoffMult2>k5Limit) {return k3dif*4;}
		k4dif=maxK<4 ? k3dif : SimilarityMeasures.cosineDifference(query.counts[4], ref.counts[4]);
		if(earlyExit && k4dif*comparisonCutoffMult>k5Limit) {return k4dif*2;}
		k5dif=maxK<5 ? k4dif : SimilarityMeasures.cosineDifference(query.counts[5], ref.counts[5]);
		return k5dif;
	}
	
	/**
	 * Compares Clades using Hellinger distance between k-mer frequency distributions.
	 * Uses early exit optimizations to avoid unnecessary calculations when possible.
	 * 
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated Hellinger distance (lower is more similar)
	 */
	private float compareHEL(float k5Limit) {
		k3dif=SimilarityMeasures.hellingerDistance(query.counts[3], ref.counts[3]);
		if(earlyExit && k3dif*comparisonCutoffMult2>k5Limit) {return k3dif*4;}
		k4dif=maxK<4 ? k3dif : SimilarityMeasures.hellingerDistance(query.counts[4], ref.counts[4]);
		if(earlyExit && k4dif*comparisonCutoffMult>k5Limit) {return k4dif*2;}
		k5dif=maxK<5 ? k4dif : SimilarityMeasures.hellingerDistance(query.counts[5], ref.counts[5]);
		return k5dif;
	}
	
	/**
	 * Compares Clades using Euclidean distance between k-mer frequency vectors.
	 * Uses early exit optimizations to avoid unnecessary calculations when possible.
	 * 
	 * @param k5Limit Maximum allowed 5-mer profile difference
	 * @return The calculated Euclidean distance (lower is more similar)
	 */
	private float compareEUC(float k5Limit) {
		k3dif=SimilarityMeasures.euclideanDistance(query.counts[3], ref.counts[3]);
		if(earlyExit && k3dif*comparisonCutoffMult2>k5Limit) {return k3dif*4;}
		k4dif=maxK<4 ? k3dif : SimilarityMeasures.euclideanDistance(query.counts[4], ref.counts[4]);
		if(earlyExit && k4dif*comparisonCutoffMult>k5Limit) {return k4dif*2;}
		k5dif=maxK<5 ? k4dif : SimilarityMeasures.euclideanDistance(query.counts[5], ref.counts[5]);
		return k5dif;
	}
	
	/**
	 * Resets all difference measures to their default values.
	 */
	void clearDif() {
		gcdif=entdif=strdif=k3dif=k4dif=k5dif=ssudif=1;
	}
	
	/**
	 * Copies all values from another Comparison object to this one.
	 * 
	 * @param b The source Comparison to copy from
	 */
	public synchronized void setFrom(Comparison b) {
		query=b.query;
		ref=b.ref;
		gcdif=b.gcdif;
		entdif=b.entdif;
		strdif=b.strdif;
		k3dif=b.k3dif;
		k4dif=b.k4dif;
		k5dif=b.k5dif;
		ssudif=b.ssudif;
	}
	
	public final boolean align(IDAligner ssa){
		final byte[] q, r;
		ssudif=1;
		if(query.r16S!=null && ref.r16S!=null) {q=query.r16S; r=ref.r16S;}
		else if(query.r18S!=null && ref.r18S!=null) {q=query.r18S; r=ref.r18S;}
		else {return false;}
		
		float id=ssa.align(q, r);
		ssudif=1-id;
		return id>0;
	}
	
	/**
	 * Returns a string representation of this Comparison.
	 * 
	 * @return A descriptive string, or null if the reference Clade is null
	 */
	public String toString() {
		if(ref==null) {return null;}
		return toBytes(null).toString();
	}
	
	/**
	 * Creates a detailed text representation of this Comparison.
	 * Includes reference taxon ID, name, and all difference measures.
	 * If a taxonomy tree is available, also includes taxonomic information.
	 * 
	 * @param bb ByteBuilder to append to, or null to create a new one
	 * @return The ByteBuilder with appended Comparison information
	 */
	public synchronized ByteBuilder toBytes(ByteBuilder bb) {
		if(ref==null) {return bb;}
		if(bb==null) {bb=new ByteBuilder();}
		bb.append("tid=").append(ref.taxID).append("\tname=").append(ref.name).nl();
		bb.append("gcdif=").append(gcdif, 5).append("\tsdif=").append(strdif, 5);
		if(calcCladeEntropy || entdif<0.5f) {bb.append("\tedif=").append(entdif, 5);}
		bb.append("\tk3dif=").append(k3dif, 6).append("\tk4dif=").append(k4dif, 6);
		bb.append("\tk5dif=").append(k5dif, 6);
		if(ssudif<1) {bb.append("\tssu=").append(1-ssudif, 4);}
		if(ref!=null && (tree!=null || ref.lineage!=null)) {bb.nl().append(ref.lineage());}
		return bb;
	}
	
	/**
	 * Appends human-readable comparison results to a ByteBuilder.
	 * Includes both query and top hit information.
	 * 
	 * @param bb ByteBuilder to append to, or null to create a new one
	 * @return The ByteBuilder with appended human-readable results
	 */
	synchronized ByteBuilder appendResultHuman(ByteBuilder bb, int hitNum) {
		if(bb==null) {bb=new ByteBuilder();}
		if(hitNum==0) {bb.append("Query:\t").append(query.toString()).nl();}
		bb.append("Hit").append(hitNum+1).colon().tab();
		return toBytes(bb);
	}
	
	/**
	 * Creates a header line for machine-readable output format.
	 * 
	 * @param printQTID Whether to include the query taxon ID column
	 * @return ByteBuilder containing the header line
	 */
	public static ByteBuilder machineHeader(boolean printQTID) {
		ByteBuilder bb=new ByteBuilder();
		bb.appendt("#QueryName");
		if(printQTID) {bb.appendt("Q_TaxID");}
		bb.appendt("Q_GC");
		bb.appendt("Q_Bases");
		bb.appendt("Q_Contigs");
		bb.appendt("RefName");
		bb.appendt("R_TaxID");
		bb.appendt("R_GC");
		bb.appendt("R_Bases");
		bb.appendt("R_Contigs");
		bb.appendt("R_Level");
		bb.appendt("GCdif");
		bb.appendt("STRdif");
		if(calcCladeEntropy) {bb.appendt("ENTdif");}
		bb.appendt("k3dif");
		bb.appendt("k4dif");
		bb.append("k5dif");
		if(Clade.callSSU) {bb.append("ssuID");}
		bb.append("\tlineage");
		
		if(false) {bb.append("\tconfidence");}
		return bb;
	}
	
	/**
	 * Appends machine-readable comparison results to a ByteBuilder.
	 * Format matches the header created by machineHeader().
	 * 
	 * @param printQTID Whether to include the query taxon ID
	 * @param bb ByteBuilder to append to, or null to create a new one
	 * @return The ByteBuilder with appended machine-readable results
	 */
	synchronized ByteBuilder appendResultMachine(boolean printQTID, ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		Clade q=query;
		Clade r=ref;
		
		bb.appendt(q.name.replace('\t', ' '));
		if(printQTID) {bb.appendt(q.taxID);}
		bb.appendt(q.gc, 3);
		bb.appendt(q.bases);
		bb.appendt(q.contigs);
		bb.appendt(r.name.replace('\t', ' '));
		bb.appendt(r.taxID);
		bb.appendt(r.gc, 3);
		bb.appendt(r.bases);
		bb.appendt(r.contigs);
		bb.appendt(TaxTree.levelToString(r.level));
		bb.appendt(gcdif, 5);
		bb.appendt(strdif, 5);
		if(calcCladeEntropy) {bb.appendt(entdif, 5);}
		bb.appendt(k3dif, 5);
		bb.appendt(k4dif, 5);
		bb.append(k5dif, 5);
		if(Clade.callSSU) {bb.append(1-ssudif, 4);}
		bb.tab().append(lineage());
		
		if(false) {bb.append("\tNA");}
		return bb;
	}
	
	/**
	 * Gets the taxonomic lineage of the reference Clade.
	 * 
	 * @return Formatted taxonomic lineage string
	 */
	CharSequence lineage() {
		return ref.lineage();
	}
	
	/**
	 * Calculates a weighted composite measure for comparison ranking.
	 * Gives highest weight to 5-mer differences, with smaller weights for other measures.
	 * 
	 * @return The weighted comparison value (lower is more similar)
	 */
	float pivot() {//This could be a field, if it helped, which it doesn't.  A NN output might.
		return (k5dif+k4dif*0.05f+k3dif*0.03f+strdif*0.02f+8*ssudif);
	}
	
	/**
	 * Compares two Comparison objects to determine relative ordering.
	 * Uses a multi-tier comparison based on k-mer, GC, and other differences.
	 * 
	 * @param b The Comparison to compare against
	 * @return -1 if this is more similar than b, 1 if less similar, 0 if equal
	 */
	@Override
	public int compareTo(Comparison b) {
		if(k5dif<1 || b.k5dif<1 || ssudif<1 || b.ssudif<1) {
			float pa=pivot();
			float pb=b.pivot();
			return pa<pb ? -1 : pa>pb ? 1 : 0;
		}
		if(k5dif!=b.k5dif) {return k5dif<b.k5dif ? -1 : 1;}
		if(k4dif!=b.k4dif) {return k4dif<b.k4dif ? -1 : 1;}
		if(k3dif!=b.k3dif) {return k3dif<b.k3dif ? -1 : 1;}
		if(gcdif!=b.gcdif) {return gcdif<b.gcdif ? -1 : 1;}
		if(strdif!=b.strdif) {return strdif<b.strdif ? -1 : 1;}
		if(entdif!=b.entdif) {return entdif<b.entdif ? -1 : 1;}
		return 0;
	}
	
	/**
	 * Determines if this comparison represents a correct classification.
	 * A correct classification is when query and reference have the same taxonomic ID.
	 * 
	 * @return true if both query and reference have the same non-zero taxonomic ID
	 */
	boolean correct() {
		return ref!=null && query.taxID!=0 && query.taxID==ref.taxID;
	}
	
	/**
	 * Determines if this comparison represents an incorrect classification.
	 * An incorrect classification is when query and reference have different taxonomic IDs.
	 * 
	 * @return true if both query and reference have different non-zero taxonomic IDs
	 */
	boolean incorrect() {
		return ref!=null && query.taxID!=0 && ref.taxID!=0 && query.taxID!=ref.taxID;
	}
	
	/**
	 * Finds the taxonomic level of the lowest common ancestor between query and reference.
	 * Higher return values indicate more distant relationships.
	 * 
	 * @return Taxonomic level of the common ancestor, or TaxTree.LIFE if not determinable
	 */
	int correctLevel() {
		if(ref==null || query==null || query.taxID==0 || 
				ref.taxID==0 || BinObject.tree==null) {
			return TaxTree.LIFE;
		}
		return BinObject.tree.commonAncestorLevel(query.taxID, ref.taxID);
	}
	
	/**
	 * Gets the current comparison cutoff multiplier value.
	 * 
	 * @return The current comparisonCutoffMult value
	 */
	public static float ccm() {return comparisonCutoffMult;}
	
	/**
	 * Sets the comparison cutoff multiplier used in early exit optimizations.
	 * 
	 * @param f The new multiplier value (must be non-negative)
	 */
	public static void setComparisonCutoffMult(float f) {
		assert(f>=0);
		comparisonCutoffMult=f;
//		comparisonCutoffMult2=f*f;
	}
	
	/**
	 * Sets the secondary comparison cutoff multiplier used in early exit optimizations.
	 * 
	 * @param f The new multiplier value (must be non-negative)
	 */
	public static void setComparisonCutoffMult2(float f) {
		assert(f>=0);
		comparisonCutoffMult2=f;
	}
	
	/** The query Clade being compared */
	Clade query, ref;
	
	/** GC content difference between query and reference */
	float gcdif=1;
	/** Entropy difference between query and reference */
	float entdif=1;
	/** Strandedness difference between query and reference */
	float strdif=1;
	/** 3-mer profile difference between query and reference */
	float k3dif=1;
	/** 4-mer profile difference between query and reference */
	float k4dif=1;
	/** 5-mer profile difference between query and reference */
	float k5dif=1;
	
	float ssudif=1;

	/** Whether to use early exit optimization to avoid unnecessary calculations */
	static boolean earlyExit=true;
	/** Multiplier for 4-mer cutoff in early exit tests */
	private static float comparisonCutoffMult=1f;
	/** Multiplier for 3-mer cutoff in early exit tests; 1.0 is better for ABS, 1.5 for ABSCOMP */
	private static float comparisonCutoffMult2=1.5f;
	/** Constants for different comparison methods */
	static final int ABS=1, COS=2, HEL=3, EUC=4, ABSCOMP=5;
	/** The current comparison method (default is ABSCOMP) */
	static int method=ABSCOMP;
	/** Maximum k-mer size to use in comparisons */
	static int maxK=5;
}