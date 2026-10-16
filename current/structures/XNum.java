package structures;

/**
 * Mutable pair of a public object reference and a public long value.
 *
 * @param <X> Type of the stored object reference
 */
public class XNum<X>{

	/** Stores an object reference and a long value.
	 * @param x_ Object reference
	 * @param num_ Long value */
	public XNum(X x_, long num_){
		x=x_;
		num=num_;
	}

	/** Stored object reference. */
	public X x;
	/** Associated long value. */
	public long num;

}
