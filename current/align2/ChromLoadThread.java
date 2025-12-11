package align2;

import dna.ChromosomeArray;
import shared.Shared;

/**
 * Thread for concurrent loading of chromosome data from files.
 * Manages thread pool limits and coordinates parallel loading of ChromosomeArray objects.
 * Used for efficient bulk loading of reference genome data.
 *
 * @author Brian Bushnell
 * @date Dec 31, 2012
 */
public class ChromLoadThread extends Thread {
	
	public static void main(String[] args){
		
	}
	
	/**
	 * Constructs a loader for one chromosome file and target array slot.
	 * @param fname_ Filename to load
	 * @param id_ Target index in the chromosome array
	 * @param r_ Array to populate
	 */
	public ChromLoadThread(String fname_, int id_, ChromosomeArray[] r_){
		fname=fname_;
		id=id_;
		array=r_;
	}
	
	/**
	 * Factory method to create and start a chromosome loading thread.
	 * Respects concurrency limits and only creates thread if array position is empty.
	 *
	 * @param fname Filename of chromosome data to load
	 * @param id Array index where loaded data will be stored
	 * @param r Array of ChromosomeArray objects to populate
	 * @return ChromLoadThread instance or null if position already filled
	 */
	public static ChromLoadThread load(String fname, int id, ChromosomeArray[] r){
		assert(r[id]==null);
		ChromLoadThread clt=null;
		if(r[id]==null){
			increment(1);
			clt=new ChromLoadThread(fname, id, r);
			clt.start();
		}
		return clt;
	}
	
	/**
	 * Loads multiple chromosome files in parallel using filename pattern.
	 * Pattern uses '#' as placeholder for chromosome number.
	 * Blocks until all chromosomes are loaded before returning.
	 *
	 * @param pattern Filename pattern with '#' placeholder for chromosome numbers
	 * @param min Minimum chromosome number to load (inclusive)
	 * @param max Maximum chromosome number to load (inclusive)
	 * @param r Array to store loaded chromosomes, created if null
	 * @return Array containing loaded ChromosomeArray objects
	 */
	public static ChromosomeArray[] loadAll(String pattern, int min, int max, ChromosomeArray[] r){
		if(r==null){r=new ChromosomeArray[max+1];}
		assert(r.length>=max+1);
		
		int pound=pattern.lastIndexOf('#');
		String a=pattern.substring(0, pound);
		String b=pattern.substring(pound+1);
		
		ChromLoadThread[] clta=new ChromLoadThread[max];
		for(int i=min; i<max; i++){
			String fname=(a+i+b);
			clta[i]=load(fname, i, r);
		}
		
		if(max>=min){ //Load last element in this thread instead of making a new thread.
			increment(1);
			r[max]=ChromosomeArray.read(a+max+b);
			increment(-1);
		}
		
		for(int i=min; i<max; i++){
			while(r[i]==null){
				synchronized(lock){
					while(lock[0]>0){
						try {
							lock.wait();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						lock.notify();
					}
				}
			}
		}
		
		return r;
	}
	
	/** Thread execution method that loads chromosome data from file.
	 * Decrements thread counter when complete or on exception. */
	@Override
	public void run(){
		try {
			array[id]=ChromosomeArray.read(fname);
		} catch (Exception e) {
			increment(-1);
			throw new RuntimeException(e);
		}
		increment(-1);
	}
	
	/**
	 * Thread-safe counter for managing concurrent loading operations.
	 * Blocks when maximum concurrent threads reached, notifies when threads complete.
	 * @param i Value to add to counter (negative values decrement)
	 * @return Current counter value after modification
	 */
	private static final int increment(int i){
		int r;
		synchronized(lock){
			if(i<=0){
				lock[0]+=i;
				lock.notify();
			}else{
				while(lock[0]>=MAX_CONCURRENT){
					try {
						lock.wait();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
			}
			r=lock[0];
		}
		return r;
	}
	
	private final int id;
	private final String fname;
	private final ChromosomeArray[] array;
	
	public static final int[] lock=new int[1];
	public static int MAX_CONCURRENT=Shared.threads();
	
}
