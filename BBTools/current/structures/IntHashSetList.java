package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import shared.Shared;
import shared.Timer;

/**
 * Extends IntHashSet to track new numbers added for rapid clearing.
 * This is synchronized with clear, but not remove.
 * @author Brian Bushnell
 * @date September 13, 2017
 *
 */
public class IntHashSetList extends IntHashSet{
	
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		IntHashSetList set=new IntHashSetList(20, 0.7f);
		HashSet<Integer> set2=new HashSet<Integer>(20, 0.7f);
		ArrayList<Integer> list=new ArrayList<Integer>();
		ArrayList<Integer> list2=new ArrayList<Integer>();
		for(int i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.contains(i));
			list.add(Integer.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			list2.add(r);
		}
		
		for(int x : list){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(int x : list){
			assert(set.contains(x));
			assert(set2.contains(x));
			set.remove(x);
			set2.remove(x);
			assert(!set.contains(x));
			assert(!set2.contains(x));
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		
		for(int x : list2){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(int x : list2){
			assert(set.contains(x));
			assert(set2.contains(x));
			set.remove(x);
			set2.remove(x);
			assert(!set.contains(x));
			assert(!set2.contains(x));
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		
		int count=4000000;
		int runs=32;
		IntList ll=new IntList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextInt());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("IntHashSetList:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set.add(x);
//				}
				final int[] y=ll.array;
				for(int z=0; z<count; z++){
					final int value=y[z];
					set.add(value);
					set.contains(value);
					set.remove(value);
					set.add(value);
				}
//				for(int x : ll.array){
//					set.remove(x);
//				}
//				set.clear();
//				assert(set.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();
			
//			System.err.println("HashSet:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set2.add(x);
//				}
//				for(int x : ll.array){
//					set2.remove(x);
//				}
//				assert(set2.isEmpty());
////				System.err.println("Finished run "+i);
//			}
//			t.stop();
//			System.err.println(t);
//			Shared.printMemory();
		}
		t.stop();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public IntHashSetList(){
		super();
		list=new IntList(4);
	}
	
	public IntHashSetList(int initialSize){
		super(initialSize);
		list=new IntList(4);
	}
	
	public IntHashSetList(int initialSize, float loadFactor_){
		super(initialSize, loadFactor_);
		list=new IntList(4);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void clear(){
		assert(size()==list.size) : list.size+", "+size()+"\n"+list.toString()+"\n"+Arrays.toString(toArray())+"\n"+Arrays.toString(super.toArray())+"\n";
		if(size()<1){return;}
		if(size()>0.25*sizeLimit()){
			super.clear();
		}else{
			for(int i=0; i<list.size; i++){
				boolean b=remove(list.get(i));
				assert(b) : list.get(i)+", "+list.size+", "+size()+"\n"+list.toString()+"\n"+Arrays.toString(toArray())+"\n";
			}
		}
		list.clear();
	}
	
	@Override
	public boolean add(int value){
//		boolean b=super.contains(value);
		if(super.add(value)){
//			assert(!b);
//			assert(super.contains(value));
			list.add(value);
			return true;
		}
		return false;
	}
	
	@Override
	public int[] toArray(){
		int[] r=list.toArray();
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean verify(){
		if(size()!=list.size){return false;}
		if(list.containsDuplicates()){return false;}
		return super.verify();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private IntList list;
	
}
