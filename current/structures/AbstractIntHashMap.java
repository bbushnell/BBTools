package structures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 8, 2017
 *
 */
public abstract class AbstractIntHashMap{
	
	public static final void test(AbstractIntHashMap set){
		Random randy2=Shared.threadLocalRandom();
		HashMap<Integer, Integer> set2=new HashMap<Integer, Integer>(20, 0.7f);
		ArrayList<Integer> klist=new ArrayList<Integer>();
		ArrayList<Integer> klist2=new ArrayList<Integer>();
		ArrayList<Integer> vlist=new ArrayList<Integer>();
		ArrayList<Integer> vlist2=new ArrayList<Integer>();
		
		for(int i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.containsKey(i));
			klist.add(Integer.valueOf(i));
			vlist.add(Integer.valueOf(i*2+7));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			klist2.add(r);
			vlist2.add(randy2.nextInt()&Integer.MAX_VALUE);
		}
		
		
		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			set.put(k, v);
			set2.put(k, v);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			assert(set.size()==set2.size());
		}
		assert(!set.isEmpty());
		assert(!set2.isEmpty());
		assert(set.size()==set2.size());
		
		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
		}
		
		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			set.increment(k);
			assert(set.get(k)==v+1);
			set.increment(k, -1);
			assert(set.get(k)==v);
		}
		
		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			set.remove(k);
			set2.remove(k);
			assert(!set.containsKey(k));
			assert(!set2.containsKey(k));
			assert(set.size()==set2.size());
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		assert(set.size()==set2.size());
		
		
		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			set.put(k, v);
			set2.put(k, v);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			assert(set.size()==set2.size());
		}
		assert(!set.isEmpty());
		assert(!set2.isEmpty());
		assert(set.size()==set2.size());
		
		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
		}
		
		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			set.remove(k);
			set2.remove(k);
			assert(!set.containsKey(k));
			assert(!set2.containsKey(k));
			assert(set.size()==set2.size());
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		assert(set.size()==set2.size());
		
		
		int count=4000000;
		int runs=32;
		IntList kil=new IntList(count);
		IntList vil=new IntList(count);
		for(int i=0; i<count; i++){
			kil.add(randy2.nextInt());
			vil.add(randy2.nextInt()&Integer.MAX_VALUE);
		}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("IntHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set.add(x);
//				}
				final int[] kila=kil.array;
				final int[] vila=vil.array;
				for(int z=0; z<count; z++){
					final int key=kila[z];
					final int value=vila[z];
					set.set(key, value);
					set.contains(key);
					set.remove(key);
					set.set(key, value);
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
			
//			System.err.println("HashMap:");
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
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public abstract void clear();
	
	public final boolean contains(int key){
		return findCell(key)>=0;
	}
	
	public final boolean containsKey(int key){
		return findCell(key)>=0;
	}
	
	public abstract int get(int key);
	
	/**
	 * Set key to value.
	 * @param key
	 * @param value
	 * @return Old value.
	 */
	public abstract int put(int key, int value);
	
	/**
	 * Set key to value.
	 * @param key
	 * @param value
	 * @return Old value.
	 */
	public abstract int set(int key, int value);
	
	/**
	 * Increment key's value by 1.
	 * @param key
	 * @return New value.
	 */
	public abstract int increment(int key);
	
	/**
	 * Increment key's value.
	 * @param key
	 * @param incr
	 * @return New value.
	 */
	public abstract int increment(int key, int incr);

	public final void incrementAll(AbstractIntHashMap map) {
		final int[] keys=map.keys();
		final int[] values=map.values();
		final int invalid=map.invalid();
		for(int i=0; i<keys.length; i++){
			int a=keys[i], b=values[i];
			if(a!=invalid){
				increment(a, b);
			}
		}
	}
	
	/**
	 * Remove this key.
	 * @param key
	 * @return true if the key was removed, false if it was not present.
	 */
	public abstract boolean remove(int key);
	
	public final int cardinality(){return size();}
	public abstract int size();
	public abstract boolean isEmpty();
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final String toString(){
		return toStringListView();
	}
	
	public final String toStringSetView(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys().length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+"("+i+", "+keys[i]+", "+values[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public final String toStringListView(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys().length; i++){
			if(keys[i]!=invalid){
				sb.append(comma);
				sb.append('(');
				sb.append(keys[i]);
				sb.append(',');
				sb.append(values[i]);
				sb.append(')');
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public final int[] toKeyArray(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys();
		int[] x=new int[size];
		int i=0;
		for(int v : keys){
			if(v!=invalid){
				x[i]=v;
				i++;
			}
		}
		return x;
	}

	public final IntHashMapBinary toCountHistogram() {
		IntHashMapBinary counts=new IntHashMapBinary(Tools.mid(1, size(), 64));
		final int[] keys=keys();
		final int[] values=values();
		final int invalid=invalid();
		for(int i=0; i<keys.length; i++) {
			int a=keys[i];
			int b=values[i];
			if(a!=invalid){
				counts.increment(b);
			}
		}
		return counts;
	}
	
	public final boolean verify(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		int numValues=0;
		int numFound=0;
		if(keys.length!=values.length){return false;}
		for(int i=0; i<keys.length; i++){
			final int key=keys[i];
			if(key==invalid){
				if(values[i]!=0){return false;}
			}else{
				numValues++;
				final int cell=findCell(i);
				if(i==cell){
					numFound++;
				}else{
					return false;
				}
			}
		}
		return numValues==numFound && numValues==size;
	}
	
	abstract int findCell(final int key);
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public abstract int[] toArray();
	public abstract int[] keys();
	public abstract int[] values();
	public abstract int invalid();
	
	static final int MASK=Integer.MAX_VALUE;
	static final int MINMASK=Integer.MIN_VALUE;
	
	static final int extra=10;
	
}
