package align2;

/**
 * Mutable row-index/length pair used to sort k-mer hit lists by their lengths.
 * Null rows have length zero. Sorting compares only value, without a key tie-break.
 * Keys retain original row indices after sorting; no matrix rows are retained.
 * Not thread-safe. The index callers use nonnegative array lengths as values.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Pointer implements Comparable<Pointer>{

	/**
	 * Creates an array of Pointers from a 2D matrix, mapping row indices to row lengths.
	 * Each Pointer has key=row_index and value=row_length.
	 * @param matrix The 2D matrix to analyze (null rows treated as length 0)
	 * @return Array of Pointers with keys as row indices and values as row lengths
	 */
	public static Pointer[] loadMatrix(final int[][] matrix){
		final Pointer[] out=new Pointer[matrix.length];
		for(int i=0; i<out.length; i++){
			final int len=(matrix[i]==null ? 0 : matrix[i].length);
			out[i]=new Pointer(i, len);
		}
		return out;
	}

	/**
	 * Reuses existing Pointer array to map matrix row indices to row lengths.
	 * Updates the provided array in-place rather than creating new objects.
	 *
	 * @param matrix The 2D matrix to analyze (null rows treated as length 0)
	 * @param out Same-length array of existing, nonnull Pointer objects; contents
	 * are overwritten in array order even if the pointers were previously sorted
	 * @return The updated Pointer array for method chaining
	 */
	public static Pointer[] loadMatrix(final int[][] matrix, final Pointer[] out){
		assert(out!=null) : "Reusable row pointers must already be allocated before loading matrix lengths";
		assert(out.length==matrix.length) : "One pointer is required per matrix row: pointers="+out.length+", rows="+matrix.length;
		for(int i=0; i<out.length; i++){
			final Pointer p=out[i];
			assert(p!=null) : "Reuse updates existing Pointer objects; missing pointer at row "+i;
			final int len=(matrix[i]==null ? 0 : matrix[i].length);
			p.key=i;
			p.value=len;
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a new Pointer with the specified key and value.
	 * @param key_ Key (often a matrix row/column index)
	 * @param value_ Value (often a length or count) used for comparisons
	 */
	public Pointer(final int key_, final int value_){
		key=key_;
		value=value_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Compares Pointers based on their values for sorting.
	 * Returns negative if this value is less, positive if greater, 0 if equal.
	 * @param o The Pointer to compare against
	 * @return Difference between this.value and o.value
	 */
	@Override
	public int compareTo(final Pointer o){
		//loadMatrix supplies nonnegative array lengths, so their difference fits int.
		//TODO: Probable bug for general use - unrestricted public value assignments
		// can make subtraction overflow. No such input was found in the index callers.
		return value-o.value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Original row index, independent of this pointer's current array position. */
	public int key;
	/** Row length for index callers; comparison assumes differences fit int. */
	public int value;
}
