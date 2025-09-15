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
	
	static int minDepth=2;
	public static float MAF_sub=0.25f;
	public static float MAF_del=0.5f;
	public static float MAF_ins=0.5f;
	public static float MAF_noref=0.4f;
	static boolean onlyConvertNs=false;
	static boolean noIndels=false;
	public static float trimDepthFraction=0.0f;
	public static boolean trimNs=false;
	
	public static boolean useMapq=false;
	public static boolean invertIdentity=false;
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
	
	static final String[] TYPE_NAMES={"DEL", "INS", "REF"};
	
	public static boolean verbose=false;
	
}
