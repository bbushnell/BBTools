package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date December 12, 2017
 *
 */
public final class LongLongHashMap{
	
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongLongHashMap map=new LongLongHashMap(20, 0.7f);
		HashMap<Long, Integer> map2=new HashMap<Long, Integer>(20, 0.7f);
		ArrayList<Long> list=new ArrayList<Long>();
		ArrayList<Long> list2=new ArrayList<Long>();
//		ArrayList<Integer> vals=new ArrayList<Integer>();
		for(long i=0; i<1000; i++){
			assert(!map.contains(i));
			assert(!map2.containsKey(i));
			list.add(Long.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			long r=randy2.nextLong();
			list2.add(r);
		}
		
		for(long x : list){
			map.put(x, (int)(2*x));
			map2.put(x, (int)(2*x));
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}
		
		for(long x : list){
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());
		
		for(long x : list2){
			map.put(x, (int)(2*x));
			map2.put(x, (int)(2*x));
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}
		
		for(long x : list2){
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());
		
		int count=4000000;
		int runs=32;
		LongList ll=new LongList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextLong());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongLongHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				final long[] y=ll.array;
				for(int z=0; z<count; z++){
					final long key=y[z];
					map.add(key);
					map.contains(key);
					map.remove(key);
					map.add(key);
				}
//				for(long x : ll.array){
//					map.remove(x);
//				}
//				map.clear();
//				assert(map.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();
			
//			System.err.println("HashMap:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map2.add(x);
//				}
//				for(long x : ll.array){
//					map2.remove(x);
//				}
//				assert(map2.isEmpty());
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
	
	public LongLongHashMap(){
		this(256);
	}
	
	public LongLongHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	public LongLongHashMap(int initialSize, float loadFactor_){
		invalid=randy.nextLong()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
//		assert(verify()); //123
	}
	
	public boolean contains(long key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}
	
	public boolean containsKey(long key){
		return contains(key);
	}
	
	public long get(long key){
//		assert(verify()); //123
		long value=-1;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;//TODO: Should return invalid...
	}
	
	/**
	 * Increment this key's value by 1.
	 * @param key
	 * @return New value
	 */
	public long add(long key){
		return increment(key, 1);
	}
	
	/**
	 * Increment this key's value by incr.
	 * @param key
	 * @param incr
	 * @return New value
	 */
	public long increment(long key, long incr){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=incr;
			size++;
//			assert(verify()); //123
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return incr;
		}else{
			values[cell]+=incr;
//			assert(verify()); //123
			return values[cell];
		}
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the key was added, false if it was already contained.
	 */
	public boolean put(long key, long value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=value;
			size++;
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return true;
		}
		assert(keys[cell]==key);
//		assert(verify()); //123
		return false;
	}
	
	/**
	 * Remove this key from the map.
	 * @param key
	 * @return Old value.
	 */
	public long remove(long key){
//		assert(verify()); //123
		if(key==invalid){return -1;}
		final int cell=findCell(key);
		if(cell<0){return -1;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		final long value=values[cell];
		values[cell]=0;
		size--;
		
		rehashFrom(cell);
//		assert(verify()); //123
		return value;
	}
	
	public int size(){return size;}
	
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+"("+i+", "+keys[i]+", "+values[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+keys[i]);
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		int i=0;
		for(long key : keys){
			if(key!=invalid){
				x[i]=key;
				i++;
			}
		}
		return x;
	}
	
	public long[] toArray(long thresh){
		int len=0;
//		assert(verify());
		for(int i=0; i<values.length; i++){
			assert((values[i]==0)==(keys[i]==invalid)) : i+", "+values[i]+", "+keys[i]+", "+invalid+"\n"+toStringSetView();
			assert((keys[i]<0)==((keys[i]==invalid))) : toStringSetView();
			if(values[i]>=thresh){
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				len++;
			}
		}
		long[] x=KillSwitch.allocLong1D(len);
		for(int i=0, j=0; j<len; i++){
			if(values[i]>=thresh){
				x[j]=keys[i];
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				j++;
			}
		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean verify(){
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final long value=values[i];
			
			if(key==invalid){
				if(value!=0){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
			}else{
				numValues++;
				if(value<1){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
				final int cell=findCell(key);
				if(i==cell){
					numFound++;
				}else{
					assert(false) : i+", "+key+", "+value+", "+cell+"\n"+((cell>=0) ? keys[cell]+", "+values[cell]+"\n" : "");
					return false;
				}
			}
		}
		boolean pass=(numValues==numFound && numValues==size);
		assert(pass) : numValues+", "+numFound+", "+size;
		return pass;
	}
	
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}
	
	private boolean rehashCell(final int cell){
		final long key=keys[cell];
		final long value=values[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(keys[dest]==invalid);
		keys[cell]=invalid;
		values[cell]=0;
		keys[dest]=key;
		values[dest]=value;
		return true;
	}
	
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				assert(values[i]==0);
				keys[i]=invalid;
			}
		}
	}
	
	private int findCell(final long key){
		if(key==invalid){return -1;}
		
		final int limit=keys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}
	
	private int findCellOrEmpty(final long key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		
		final int limit=keys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
	private final void resize(final long size2){
//		assert(verify()); //123
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final long[] oldKeys=keys;
		final long[] oldValues=values;
		keys=KillSwitch.allocLong1D(size3);
		values=KillSwitch.allocLong1D(size3);
		Arrays.fill(keys, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(int i=0; i<oldKeys.length; i++){
			long key=oldKeys[i];
			if(key!=invalid){
				put(key, oldValues[i]);
			}
		}
//		assert(verify()); //123
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	public long[] keys() {return keys;}

	public long[] values() {return values;}

	public long invalid() {return invalid;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long[] keys;
	private long[] values;
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	private static final Random randy=new Random(1);
	private static final long MASK=Long.MAX_VALUE;
	private static final long MINMASK=Long.MIN_VALUE;
	
	private static final int extra=10;
	
}
