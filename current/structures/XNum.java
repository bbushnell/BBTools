package structures;

/**
 * Generic container pairing an object with a long value.
 * Used to attach a number (count, id, or score) to an arbitrary item.
 *
 * @param <X> Type of the stored object
 */
public class XNum<X>{

	/**
	 * Creates an XNum pairing an object with a numeric value.
	 * @param x_ Object to store
	 * @param num_ Numeric value to associate
	 */
	public XNum(X x_, long num_) {
		x=x_;
		num=num_;
	}

	/** The stored object */
	public X x;
	/** The associated numeric value */
	public long num;
	
}
