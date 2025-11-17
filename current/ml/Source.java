package ml;

/**
 * Abstract base class defining fundamental properties for machine learning data sources.
 * Provides a lightweight template for data sources with basic value, ID, and bias
 * tracking capabilities.
 * Subclasses must implement validation and terminal status checking functionality.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public abstract class Source {
	
	public Source() {}
	
	public float value(){return value;}
	public int id() {return 0;}
	public float bias() {return 0;}
	
	/** Computational result value stored by the data source */
	protected float value=0f;
	
	/** Sets the computational result value for this data source.
	 * @param v The value to store */
	void setValue(float v) {
		value=v;
	}
	
	/**
	 * Validates the state or content of this data source, implementation-dependent
	 */
	public abstract boolean check();
	
	/** True only if this has no inputs */
	public abstract boolean terminal();

}
