package bloom;

/**
 * Bit-packed exact count array backed by several ordinary {@code int[]} arrays.
 * This legacy variant stores fixed-width counters directly and is not thread-safe.
 *
 * @author Brian Bushnell
 * @date Aug 17, 2012
 */
public class KCountArray3 extends KCountArray{
		
	/** Serialization compatibility marker for the count-array family. */
	private static final long serialVersionUID=-5466091642729698944L;
	
	/** Constructs a count array with power-of-two cells and supported cell width. */
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
	
	/** Stores {@code value} at {@code key} without updating the cached nonzero-cell estimate. */
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
	/** Adds {@code incr} to {@code key}, saturating at {@link #maxValue}. */
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
		value=(int)min(value+(long)incr, maxValue);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		if(verbose){System.err.println("Returning "+value);}
		//return (int)value;
	}
	
	/** Adds {@code incr} to {@code key} and returns the old count, leaving cached occupancy unchanged. */
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
		final int value2=(int)min(value+(long)incr, maxValue);
		word=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		if(verbose){System.err.println("Returning "+value);}
		return value;
	}
	
	/** Returns a histogram where index is cell value and value is frequency. */
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
	/** Returns all packed cell values as a comma-separated debug string. */
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
	
	/** Returns the fraction of cells directly observed transitioning above zero. */
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/** Returns the fraction of cells meeting a positive count threshold. */
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/** Returns the number of cells meeting a positive count threshold. */
	@Override
	public long cellsUsed(int mindepth){
		long count=0;
		for(int[] array : matrix){
			if(array!=null){
				for(int word : array){
					long word2=word&0xFFFFFFFFL;
					while(word2!=0){
						int x=(int)(word2&valueMask);
						if(x>=mindepth){count++;}
						word2>>>=cellBits;
					}
				}
			}
		}
		return count;
	}
	
	/** Unsupported for direct-addressed arrays. */
	@Override
	long hash(long x, int y){
		assert(false) : "Unsupported.";
		return x;
	}
	
	private long cellsUsed;
	private final int[][] matrix;
	
}
