package consensus;

import structures.ByteBuilder;

/**
 * Abstract superclass for consensus sequence analysis and variant calling.
 * Supplies text representation, graph-element types and shared consensus policy.
 * Policy fields are mutable process-wide configuration. Set them before graph
 * processing and coordinate any temporary override; instances do not snapshot
 * the settings and these fields provide no synchronization.
 *
 * @author Brian Bushnell
 * @date September 6, 2019
 */
public abstract class ConsensusObject{

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
	
	/** Minimum supporting observations required by consensus selection. */
	static int minDepth=2;
	/** Minimum selected-allele fraction at a defined reference base. */
	public static float MAF_sub=0.25f;
	/** Minimum deletion-observation fraction. */
	public static float MAF_del=0.5f;
	/** Minimum insertion-node observation fraction against position depth. */
	public static float MAF_ins=0.5f;
	/** Minimum selected-allele fraction when the stored reference is uppercase N. */
	public static float MAF_noref=0.4f;
	/** Protects REF symbols other than uppercase N from substitution. */
	static boolean onlyConvertNs=false;
	/** Disables insertion/deletion output during traversal. */
	static boolean noIndels=false;
	/** Relative depth threshold for trimming, scaled by the maximum position depth. */
	public static float trimDepthFraction=0.0f;
	/** Requests terminal uppercase-N trimming. Enabling this flag also activates
	 * the legacy minimum-depth-one end filter, even when trimDepthFraction is zero. */
	public static boolean trimNs=false;
	
	/** Includes mapping-quality information in accumulated observation weights. */
	public static boolean useMapq=false;
	/** Adds a low-identity adjustment to observation weights when enabled. */
	public static boolean invertIdentity=false;
	/** Ceiling used by the optional max(0, ceiling-percentIdentity) weight adjustment. */
	public static int identityCeiling=150;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Reference-position node; it can select either the reference or a substitution. */
	public static final int REF=2;
	/** Insertion-chain node. */
	public static final int INS=1;
	/** Deletion-position node. */
	public static final int DEL=0;
	
	/** Type labels indexed by the stable DEL/INS/REF integer constants. */
	static final String[] TYPE_NAMES={"DEL", "INS", "REF"};
	
	/** Enables diagnostic output in graph operations that consult this flag. */
	public static boolean verbose=false;
	
}
