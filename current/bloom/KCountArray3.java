package bloom;

/**
 * @author Brian Bushnell
 * @date Aug 17, 2012
 *
 */
public class KCountArray3 extends KCountArray {
		
	/**
	 * 
	 */
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
	
	/** Returns unincremented value */
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
	
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
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
	
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
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
	
	@Override
	long hash(long x, int y) {
		assert(false) : "Unsupported.";
		return x;
	}
	
	private long cellsUsed;
	private final int[][] matrix;
	
}
