package consensus;

import java.io.Serializable;

/**
 * Abstract base class for graph elements in consensus sequence generation.
 * Serves as parent class for BaseEdge and BaseNode, providing fundamental
 * type classification and serialization support for graph-based consensus algorithms.
 * Supports three discrete element types: Reference (REF), Insertion (INS), and Deletion (DEL).
 *
 * @author Brian Bushnell
 * @date September 6, 2019
 */
public abstract class BaseGraphPart extends ConsensusObject implements Serializable {
	
	/** Serialization version identifier for class compatibility */
	private static final long serialVersionUID = 3854022870880887972L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a graph part with the specified type, validating it is REF, INS, or DEL.
	 * @param type_ Graph element type (REF=2, INS=1, DEL=0) */
	public BaseGraphPart(int type_){
		type=type_;
		assert(type==REF || type==INS || type==DEL) : type;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the string name for this graph element type */
	public final String typeString(){
		return TYPE_NAMES[type];
	}

	/** Returns a string representation of this specific graph part */
	public abstract String partString();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Graph element type: REF (2), INS (1), or DEL (0) */
	public final int type;
	
}
