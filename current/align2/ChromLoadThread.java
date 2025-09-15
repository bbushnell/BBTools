package align2;

import dna.ChromosomeArray;
import shared.Shared;

/**
 * @author Brian Bushnell
 * @date Dec 31, 2012
 *
 */
public class ChromLoadThread extends Thread {
	
	public static void main(String[] args){
		
	}
	
	public ChromLoadThread(String fname_, int id_, ChromosomeArray[] r_){
		fname=fname_;
		id=id_;
		array=r_;
	}
	
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
