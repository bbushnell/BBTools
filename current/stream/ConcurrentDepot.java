package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

public class ConcurrentDepot<K> {
	
	@SuppressWarnings("unchecked")
	public ConcurrentDepot(int bufSize, int numBufs){
		bufferSize=bufSize;
		bufferCount=numBufs;
		
		lists=new ArrayList[numBufs];
		empty=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1, fair);
		full=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1, fair);
		
		for(int i=0; i<lists.length; i++){
			lists[i]=new ArrayList<K>(bufSize);
			empty.add(lists[i]);
		}
		
	}
	
	
	public final ArrayBlockingQueue<ArrayList<K>> empty;
	public final ArrayBlockingQueue<ArrayList<K>> full;
	
	public final int bufferSize;
	public final int bufferCount;
	
	public static boolean fair=false;
	
	private final ArrayList<K>[] lists;
	
}
