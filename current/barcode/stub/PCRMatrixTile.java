package barcode.stub;

/**
 * Public reflection target when the optional tile implementation is absent.
 * Construction is rejected by the inherited validity check; no usable tile
 * matrix or fallback computation is provided.
 *
 * @author Brian Bushnell
 * @date March 25, 2024
 */
public class PCRMatrixTile extends PCRMatrixProbAbstract{

	/**
	 * Matches the optional constructor signature; the superclass rejects construction.
	 * @param a Requested first index length (unused by the stub)
	 * @param b Requested second index length (unused by the stub)
	 * @param c Requested delimiter (unused by the stub)
	 * @param d Requested distance-sum mode (unused by the stub)
	 */
	public PCRMatrixTile(int a, int b, int c, boolean d){}

}
