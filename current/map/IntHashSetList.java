package map;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import shared.Random;

import shared.Shared;
import shared.Timer;
import structures.IntList;

/**
 * Hash set that tracks newly added integers for optimized bulk clearing.
 * Extends IntHashSet with an internal list to record added elements, enabling
 * faster clearing when the set is sparse relative to its capacity.
 * Thread-safe for clear() operations but not for remove() operations.
 *
 * Add/clear only: the tracking list is maintained by add() and reset by clear(), but remove() is
 * inherited unchanged from IntHashSet and does NOT update the list -- after a remove() the list
 * desyncs from the set, so clear() (assert size==list.size), the sparse-clear path, and verify()
 * would fail under -ea, and toArray() would return stale (removed) elements. The sole caller
 * (sketch/SketchIndex) uses only add/toArray/clear, so this is latent. See [map/IntHashSetList#001].
 *
 * @author Brian Bushnell
 * @date September 13, 2017
 */
public class IntHashSetList extends IntHashSet{
	
	/** Benchmark and validation entry point comparing IntHashSetList to standard HashSet.
	 * @param args Unused */
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
	
	/** Creates an empty set with default capacity and tracking list. */
	public IntHashSetList(){
		super();
		list=new IntList(4);
	}
	
	/** Creates an empty set with specified initial capacity.
	 * @param initialSize Initial hash table size */
	public IntHashSetList(int initialSize){
		super(initialSize);
		list=new IntList(4);
	}
	
	/**
	 * Creates an empty set with specified capacity and load factor.
	 * @param initialSize Initial hash table size
	 * @param loadFactor_ Load factor threshold for rehashing
	 */
	public IntHashSetList(int initialSize, float loadFactor_){
		super(initialSize, loadFactor_);
		list=new IntList(4);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void clear(){
		//assert the list-tracking invariant first (size==list.size); holds under add-only usage and would catch a
		//remove()-induced desync loudly HERE (-ea always on). Dense -> full super.clear(); sparse -> remove only the
		//tracked elements (cheaper than scanning the whole backing array). list.clear() resets tracking either way.
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
		//list-tracking invariant: list mirrors the set 1:1 under add-only usage. super.add returns true ONLY for a
		//genuinely new value, so list never gets a duplicate (verify() checks !containsDuplicates).
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
		//returns the tracking list (O(size)), NOT a scan of the backing array like super.toArray(). Same contents
		//under add-only usage; would include stale (removed) elements if remove() were ever used (see #001).
		int[] r=list.toArray();
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean verify(){
		//list must mirror the set exactly: same count, no duplicates, then the base set's own probe-integrity check.
		if(size()!=list.size){return false;}
		if(list.containsDuplicates()){return false;}
		return super.verify();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Possible bug [map/IntHashSetList#001] - the inherited remove() does NOT update this list, so any remove()
	//desyncs list from the set (clear()/verify() then fail under -ea; toArray() returns stale elements). Latent: the
	//sole caller sketch/SketchIndex uses add/toArray/clear only, never remove. Crash-loud fix = override remove() to
	//throw UnsupportedOperationException; deferred to Brian (API decision on his class).
	private IntList list;
	
}
