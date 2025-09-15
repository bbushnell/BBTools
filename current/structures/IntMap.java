package structures;
import java.util.Arrays;

import shared.Tools;


public class IntMap {
	
	public static void main(String[] args){
		
	}
	
	
	public IntMap(int from, int to){
		reset(from, to);
	}
	
	
	public int get(int key){
		assert(key>=min && key<=max);
		return array[key-min];
	}
	
	
	public boolean containsKey(int key){
		assert(key>=min && key<=max);
		return array[key-min]!=INVALID;
	}
	
	
	public int increment(int key){
		return increment(key, 1);
	}
	
	
	public int increment(int key, int value){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		int v2=(int)Tools.min(Integer.MAX_VALUE, old+(long)value);
		assert(array[index]!=INVALID);
		array[index]=v2;
		return v2;
	}
	
	
	public int put(int key, int value){
		assert(key>=min && key<=max);
		assert(value!=INVALID);
		int index=key-min;
		int old=array[index];
		array[index]=value;
		return old;
	}
	
	
	public int remove(int key){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		array[index]=INVALID;
		return old;
	}
	
	
	public int size(){
		int sum=0;
		for(int i=0; i<array.length; i++){
			if(array[i]!=INVALID){sum++;}
		}
		return sum;
	}
	
	
	public int[] keys(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=(min+i);
				j++;
			}
		}
		return r;
	}
	
	
	public int[] values(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=array[i];
				j++;
			}
		}
		return r;
	}
	
	
	public void clear(){
		Arrays.fill(array, INVALID);
	}
	
	
	public void reset(int from, int to){
		min=from;
		max=to;
		assert(max>=min);
		assert(((long)max)-((long)min)<Integer.MAX_VALUE);
		
		int size=max-min+1;
		if(array==null || array.length<size){
			array=new int[size];
		}
		clear();
	}
	
	
	public int min;
	public int max;
	public int[] array;
	
	private static final int INVALID=Integer.MIN_VALUE;
	
}
