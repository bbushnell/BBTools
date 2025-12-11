package bloom;

/**
 * Matrix-based k-mer counting array implementation using integer arrays.
 * Stores k-mer counts in a two-dimensional matrix with bit-packed values.
 * Provides efficient storage and access for k-mer frequency counting operations.
 *
 * @author Brian Bushnell
 * @date Aug 17, 2012
 */
public class KCountArray3 extends KCountArray {
		
	/** Serial version UID for serialization compatibility */
	private static final long serialVersionUID = -5466091642729698944L;
	
	public KCountArray3(long cells_, int bits_){
		super(cells_, bits_);
		long words=cells/cellsPerWord;
		int wordsPerArray=(int)(words/numArrays);
		matrix=new int[numArrays][wordsPerArray];
	}
	
	@Override
	public int read(long key){
		if(verbose){System.err.println("Reading "+key);}
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
//		System.out.println("array="+arrayNum);
		key>>>=arrayBits;
//		System.out.println("key2="+key);
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
//		System.out.println("index="+index);
		int word=array[index];
//		System.out.println("word="+Integer.toHexString(word));
		int cellShift=(int)(cellBits*key);
//		System.out.println("cellShift="+cellShift);
		int value=(int)((word>>>cellShift)&valueMask);
		if(verbose){System.err.println("Read "+value);}
		return value;
	}
	
	/**
	 * Writes a count value for a given key to the matrix.
	 * Uses bit manipulation to pack the value into the integer array.
	 * @param key The key to write the count for
	 * @param value The count value to store
	 */
	@Override
	public void write(long key, int value){
		if(verbose){System.err.println("Writing "+key+", "+value);}
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
	}
	
//	static int count138=0;
	/**
	 * Increments the count for a given key by the specified amount.
	 * Updates cellsUsed counter when transitioning from/to zero count.
	 * Clamps the result to maxValue to prevent overflow.
	 *
	 * @param key The key to increment
	 * @param incr The amount to increment by (can be negative)
	 */
	@Override
	public void increment(long key, int incr){
		if(verbose){System.err.println("*** Incrementing "+key);}
//		if(key==138){
//			assert(count138==0) : count138;
//			count138++;
//		}
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		int value=((word>>>cellShift)&valueMask);
		if(value==0 && incr>0){cellsUsed++;}
		else if(incr<0 && value+incr==0){cellsUsed--;}
		value=min(value+incr, maxValue);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		if(verbose){System.err.println("Returning "+value);}
		//return (int)value;
	}
	
	/**
	 * Increments the count for a key and returns the original unincremented value.
	 * Useful for operations that need both the old and new values.
	 *
	 * @param key The key to increment
	 * @param incr The amount to increment by
	 * @return The count value before incrementing
	 */
	@Override
	public int incrementAndReturnUnincremented(long key, int incr){
		if(verbose){System.err.println("Incrementing2 "+key);}
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		final int value=((word>>>cellShift)&valueMask);
		final int value2=min(value+incr, maxValue);
		word=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		if(verbose){System.err.println("Returning "+value);}
		return value;
	}
	
	/**
	 * Transforms the count data into a frequency histogram.
	 * Delegates to the parent class implementation using the matrix data.
	 * @return Array where index represents count value and value represents frequency
	 */
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
	/**
	 * Converts the entire contents of the array to a string representation.
	 * Unpacks all bit-packed values and formats them as a comma-separated list.
	 * Primarily used for debugging and small array inspection.
	 * @return String representation of all stored values
	 */
	@Override
	public String toContentsString(){
		StringBuilder sb=new StringBuilder();
		sb.append("[");
		String comma="";
		for(int[] array : matrix){
			for(int i=0; i<array.length; i++){
				int word=array[i];
				for(int j=0; j<cellsPerWord; j++){
					int x=word&valueMask;
					sb.append(comma);
					sb.append(x);
					word>>>=cellBits;
					comma=", ";
				}
			}
		}
		sb.append("]");
		return sb.toString();
	}
	
	/** Returns the fraction of cells that contain non-zero values.
	 * @return Fraction of used cells (0.0 to 1.0) */
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/**
	 * Returns the fraction of cells with counts at or above the minimum depth.
	 * @param mindepth Minimum count value to consider as "used"
	 * @return Fraction of cells meeting the depth requirement
	 */
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/**
	 * Counts the number of cells with counts at or above the minimum depth.
	 * Iterates through all matrix elements and unpacks bit-packed values.
	 * @param mindepth Minimum count value to consider as "used"
	 * @return Number of cells meeting the depth requirement
	 */
	@Override
	public long cellsUsed(int mindepth){
		long count=0;
		for(int[] array : matrix){
			if(array!=null){
				for(int word : array){
					while(word>0){
						int x=word&valueMask;
						if(x>=mindepth){count++;}
						word>>>=cellBits;
					}
				}
			}
		}
		return count;
	}
	
	/**
	 * Hash function - not supported in this implementation.
	 * Always throws an assertion error when called.
	 *
	 * @param x Input value
	 * @param y Hash parameter
	 * @return Never returns - always throws assertion error
	 */
	@Override
	long hash(long x, int y) {
		assert(false) : "Unsupported.";
		return x;
	}
	
	private long cellsUsed;
	private final int[][] matrix;
	
}
