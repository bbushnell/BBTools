package ml;

/** Base value slot for objects that expose a scalar to downstream cells.
 * Cell overrides id(), bias(), check() and terminal() for neurons with identity,
 * bias and topology; these base accessors only manage the stored value.
 *
 * @author Brian Bushnell
 * @date June 3, 2025 */
public abstract class Source {
	
	public Source() {}
	
	/** Returns the scalar value visible to downstream cells. */
	public float value(){return value;}
	/** Returns the source ID; the base implementation has none. */
	public int id(){return 0;}
	/** Returns the additive bias; the base implementation has none. */
	public float bias(){return 0;}
	
	/** Scalar value consumed by downstream cells. */
	protected float value=0f;
	
	/** Sets the scalar value visible to downstream cells.
	 * @param v value to store */
	void setValue(float v){
		value=v;
	}
	
	/** Validates implementation-specific cell state. */
	public abstract boolean check();
	
	/** True only if this source has no inputs. */
	public abstract boolean terminal();

}
