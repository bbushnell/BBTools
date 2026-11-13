package barcode.stub;

/**
 * Public reflection target when the optional probability implementation is absent.
 * Its constructor has the expected factory signature, but the inherited validity
 * check rejects construction. This is not a usable probability matrix.
 *
 * @author Brian Bushnell
 * @date March 25, 2024
 */
public class PCRMatrixProb extends PCRMatrixProbAbstract{

	/**
	 * Matches the optional constructor signature. Arguments are not forwarded;
	 * the implicit superclass constructor rejects the missing implementation.
	 *
	 * @param a Requested first index length (unused by the stub)
	 * @param b Requested second index length (unused by the stub)
	 * @param c Requested delimiter (unused by the stub)
	 * @param d Requested distance-sum mode (unused by the stub)
	 */
	public PCRMatrixProb(int a, int b, int c, boolean d){}

}
