package align2;

/**
 * Simple key-value pair data structure for matrix indexing operations.
 * Provides sortable pointers to matrix dimensions with value-based comparison.
 * Commonly used for sorting matrix rows/columns by their lengths.
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
	public static Pointer[] loadMatrix(int[][] matrix){
		Pointer[] out=new Pointer[matrix.length];
		for(int i=0; i<out.length; i++){
			int len=(matrix[i]==null ? 0 : matrix[i].length);
			out[i]=new Pointer(i, len);
		}
		return out;
	}
	
	/**
	 * Reuses existing Pointer array to map matrix row indices to row lengths.
	 * Updates the provided array in-place rather than creating new objects.
	 *
	 * @param matrix The 2D matrix to analyze (null rows treated as length 0)
	 * @param out Pre-allocated Pointer array to update (must match matrix length)
	 * @return The updated Pointer array for method chaining
	 */
	public static Pointer[] loadMatrix(int[][] matrix, Pointer[] out){
		assert(out!=null);
		assert(out.length==matrix.length);
		for(int i=0; i<out.length; i++){
			Pointer p=out[i];
			int len=(matrix[i]==null ? 0 : matrix[i].length);
			p.key=i;
			p.value=len;
		}
		return out;
	}
	
	/**
	 * Creates a new Pointer with the specified key-value pair.
	 * @param key_ The key (typically an index)
	 * @param value_ The value (typically a length or count)
	 */
	public Pointer(int key_, int value_){
		key=key_;
		value=value_;
	}
	
	@Override
	public int compareTo(Pointer o) {
		return value-o.value;
	}
	
	/** The key component, typically used as an index or identifier */
	public int key;
	/** The value component, typically used for lengths, counts, or priorities */
	public int value;
}