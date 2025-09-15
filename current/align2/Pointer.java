package align2;

public class Pointer implements Comparable<Pointer>{
	
	public static Pointer[] loadMatrix(int[][] matrix){
		Pointer[] out=new Pointer[matrix.length];
		for(int i=0; i<out.length; i++){
			int len=(matrix[i]==null ? 0 : matrix[i].length);
			out[i]=new Pointer(i, len);
		}
		return out;
	}
	
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
	
	public Pointer(int key_, int value_){
		key=key_;
		value=value_;
	}
	
	@Override
	public int compareTo(Pointer o) {
		return value-o.value;
	}
	
	public int key;
	public int value;
}