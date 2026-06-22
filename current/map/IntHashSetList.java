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
 * Add/clear only: the tracking list is maintained by add() and reset by clear(). remove() is NOT
 * supported -- it would desync the tracking list (which cannot be kept in sync without an O(n) scan),
 * so it is overridden to throw UnsupportedOperationException rather than silently corrupt
 * clear()/toArray()/verify(). The sole caller (sketch/SketchIndex) uses only add/toArray/clear.
 * See [map/IntHashSetList#001].
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
		}
		set.clear();//IntHashSetList is add/clear-only (remove() throws) -> validate via clear, not per-element remove
		set2.clear();
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
		}
		set.clear();
		set2.clear();
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
				}
				set.clear();//add/clear-only: clear each run instead of per-element remove (remove() throws)
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
	public boolean remove(int value){
		//FIXED [map/IntHashSetList#001]: the inherited remove() would silently desync the tracking list (clear()/
		//toArray()/verify() then corrupt). The list can't be kept in sync without an O(n) scan, and this subclass is
		//add/clear-only by design (sole caller sketch/SketchIndex), so remove() is UNSUPPORTED -> crash loud, never
		//corrupt (BBTools crash-loud-never-wrong). Brian-approved 2026-06-22.
		throw new UnsupportedOperationException("IntHashSetList does not support remove() (it would desync the "
			+"tracking list); use add/clear cycles, or a plain IntHashSet if you need remove().");
	}

	@Override
	public int[] toArray(){
		//returns the tracking list (O(size)), NOT a scan of the backing array like super.toArray(). Equals the set's
		//contents because the list is kept in sync by add/clear and remove() is unsupported (throws). See #001.
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
	
	//FIXED [map/IntHashSetList#001]: the inherited remove() would desync this list from the set; remove() is now
	//overridden to throw UnsupportedOperationException (crash-loud, not silent corruption). add/clear keep it in sync.
	private IntList list;
	
}
