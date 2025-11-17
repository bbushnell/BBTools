package consensus;

import structures.ByteBuilder;

/**
 * Superclass for consensus package classes.
 * 
 * @author Brian Bushnell
 * @date September 6, 2019
 *
 */
public abstract class ConsensusObject {

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Return the text representation of this object */
	public abstract ByteBuilder toText();
	
	@Override
	public final String toString(){return toText().toString();}
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum coverage depth required to call a consensus base */
	static int minDepth=2;
	/** Minimum allele frequency threshold for substitution variants */
	public static float MAF_sub=0.25f;
	/** Minimum allele frequency threshold for deletion variants */
	public static float MAF_del=0.5f;
	/** Minimum allele frequency threshold for insertion variants */
	public static float MAF_ins=0.5f;
	/** Minimum allele frequency threshold when no reference base is available */
	public static float MAF_noref=0.4f;
	/** If true, only convert ambiguous N bases to consensus calls */
	static boolean onlyConvertNs=false;
	/** If true, disable insertion and deletion variant calling */
	static boolean noIndels=false;
	/** Fraction of maximum depth below which bases are trimmed from ends */
	public static float trimDepthFraction=0.0f;
	/** If true, trim ambiguous N bases from sequence ends */
	public static boolean trimNs=false;
	
	/** If true, incorporate mapping quality scores into consensus calculations */
	public static boolean useMapq=false;
	/** If true, invert identity calculations for consensus scoring */
	public static boolean invertIdentity=false;
	/** Maximum identity score value for consensus calculations */
	public static int identityCeiling=150;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/* Possible types */
	/** Match/Sub, neutral-length node or edge to the next REF node */
	public static final int REF=2;
	/** Insertion node or edge to an insertion node */
	public static final int INS=1;
	/** Edge to a non-adjacent node */
	public static final int DEL=0;
	
	/** String names corresponding to variant type constants DEL, INS, REF */
	static final String[] TYPE_NAMES={"DEL", "INS", "REF"};
	
	/** If true, enable verbose output for debugging consensus operations */
	public static boolean verbose=false;
	
}
