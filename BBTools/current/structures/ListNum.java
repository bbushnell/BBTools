package structures;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import shared.Shared;
import stream.Read;

public final class ListNum<K extends Serializable> implements Serializable, Iterable<K> {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7509242172010729386L;

	public ListNum(ArrayList<K> list_, long id_){
		list=list_;
		id=id_;
		if(GEN_RANDOM_NUMBERS && list!=null){
			for(K k : list){
				if(k!=null){
					((Read)k).rand=randy.nextDouble();
				}
			}
		}
	}
	
	public final int size(){
		return list==null ? 0 : list.size();
	}
	
	@Override
	public String toString(){return list==null ? "ln.list=null" : list.toString();}
	
	public final boolean isEmpty() {return list==null || list.isEmpty();}

	public final K get(int i){return list.get(i);}
	public final K set(int i, K k){return list.set(i, k);}
	public final K remove(int i){return list.remove(i);}
	public final void add(K k){list.add(k);}
	public final void clear(){list.clear();}
	
	@Override
	public Iterator<K> iterator() {return list==null ? null : list.iterator();}
	
	public final ArrayList<K> list;
	public final long id;
	
	public static synchronized void setDeterministicRandomSeed(long seed_){
		if(seed_>=0){seed=seed_;}
		else{seed=System.nanoTime()+(long)(Math.random()*10000000);}
	}
	
	public static synchronized void setDeterministicRandom(boolean b){
		GEN_RANDOM_NUMBERS=b;
		if(b){
			randy=Shared.threadLocalRandom(seed);
			seed++;
		}
	}
	public static boolean deterministicRandom(){
		return GEN_RANDOM_NUMBERS;
	}
	
	private static boolean GEN_RANDOM_NUMBERS=false;
	private static Random randy;
	private static long seed=0;
}
