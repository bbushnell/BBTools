package sketch;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import shared.Tools;

public class AlignmentThreadPool {
	
	public AlignmentThreadPool(int maxThreads_) {
		maxThreads=maxThreads_;
		assert(maxThreads>0);
		tlist=new ArrayList<AlignmentThread>(maxThreads);
	}
	
	public void addJobs(ArrayList<Comparison> list, int maxRecords){
		if(list==null || list.isEmpty() || maxRecords<1){return;}
		final int limit=Tools.min(list.size(), maxRecords);
		ArrayBlockingQueue<Comparison> dest=new ArrayBlockingQueue<Comparison>(limit);
		int added=0;
		for(int i=0; i<limit; i++){
			Comparison c=list.get(i);
			if(c.needsAlignment()){
				addJob(c, dest);
				added++;
			}
		}
		for(int i=0; i<added; i++){
			take(dest);
		}
	}
	
	public void addJob(Comparison c, ArrayBlockingQueue<Comparison> dest){
		if(tlist.size()<maxThreads){spawnThread();}
		assert(!poisoned);
		AlignmentJob job=new AlignmentJob(c, dest);
		put(job);
	}
	
	private synchronized void spawnThread(){
		final int size=tlist.size();
		if(size<maxThreads && busy.get()>=size){
//			AlignmentThread alt=new AlignmentThread(source, busy);
			AlignmentThread alt=new AlignmentThread();
			tlist.add(alt);
			alt.start();
		}
	}
	
	synchronized void poison(){
		assert(!poisoned);
		if(poisoned){return;}
		put(poison);
		poisoned=true;
	}
	
	private void put(AlignmentJob job){
		busy.incrementAndGet();
		boolean success=false;
		while(!success){
			try {
				source.put(job);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private final <X> X take(ArrayBlockingQueue<X> queue){
		X x=null;
		while(x==null){
			try {
				x=queue.take();
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		return x;
	}
	
	private class AlignmentThread extends Thread {
		
		AlignmentThread(){}
		
		private final AlignmentJob next(){
			return take(source);
		}
		
		@Override
		public void run(){
			AlignmentJob job=null;
			for(job=next(); !job.isPoison(); job=next()){
				job.doWork();
				busy.decrementAndGet();
			}
			put(poison);
		}
		
//		private final ArrayBlockingQueue<AlignmentJob> source;
//		private final AtomicInteger busy;
//
//		private static final AlignmentJob poison=new AlignmentJob(null, null);
		
	}
	
	final ArrayList<AlignmentThread> tlist;
	final int maxThreads;
	final AtomicInteger busy=new AtomicInteger(0);
	private boolean poisoned=false;

	private static final AlignmentJob poison=new AlignmentJob(null, null);
	private static final ArrayBlockingQueue<AlignmentJob> source=new ArrayBlockingQueue<AlignmentJob>(4096);
	
}
