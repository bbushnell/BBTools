package consensus;

import structures.ByteBuilder;

/**
 * Abstract superclass for consensus sequence analysis and variant calling.
 * Provides common framework for objects that represent consensus data,
 * including minimum allele frequency thresholds and depth requirements.
 *
 * @author Brian Bushnell
 * @date September 6, 2019
 */
public abstract class ConsensusObject {

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the text representation of this consensus object.
	 * Implementations should generate a formatted string representation
	 * suitable for output or display.
	 * @return ByteBuilder containing the text representation
	 */
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
	public static final int REF=2;
	public static final int INS=1;
	public static final int DEL=0;
	
	static final String[] TYPE_NAMES={"DEL", "INS", "REF"};
	
	public static boolean verbose=false;
	
}
