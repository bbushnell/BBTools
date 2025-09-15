package ukmer;

import java.util.concurrent.atomic.AtomicInteger;

import shared.Shared;
import shared.Tools;
import structures.SuperLongList;

public final class HistogramMakerU {
	
	public static long[] fillHistogram(final AbstractKmerTableU[] tables, final int histMax) {
		if(Shared.threads()>2){
			return fillHistogram_MT(tables, histMax);
		}else{
			return fillHistogram_ST(tables, histMax);
		}
	}
	
	private static long[] fillHistogram_ST(final AbstractKmerTableU[] tables, final int histMax) {
		long[] ca=new long[histMax+1];
		for(AbstractKmerTableU set : tables){
			set.fillHistogram(ca, histMax);
		}
		return ca;
	}
	
	private static long[] fillHistogram_MT(final AbstractKmerTableU[] tables, final int histMax) {
		boolean errorState=false;
		int threads=Shared.threads();
		threads=Tools.min((threads>20 ? threads/2 : threads), (tables.length+1)/2, 32);
		if(threads<2){return fillHistogram_ST(tables, histMax);}
		
		final FillThread[] array=new FillThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<threads; i++){array[i]=new FillThread(tables, histMax, next);}
		for(int i=0; i<threads; i++){array[i].start();}
		
		//Wait for completion of all threads
		final long[] ca=new long[histMax+1];
		boolean success=true;
		for(FillThread pt : array){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			//Accumulate per-thread statistics
			
			pt.sll.addTo(ca);
			pt.sll=null;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		return ca;
	}
	
	private static class FillThread extends Thread{
		
		FillThread(final AbstractKmerTableU[] tables_, int histMax_, AtomicInteger next_){
			tables=tables_;
			next=next_;
			sll=new SuperLongList(Tools.mid(5000, histMax_, 100000));
		}
		
		@Override
		public void run(){
			for(int tnum=next.getAndIncrement(); tnum<tables.length; tnum=next.getAndIncrement()){
				tables[tnum].fillHistogram(sll);
			}
		}
		
		final AbstractKmerTableU[] tables;
		final AtomicInteger next;
		SuperLongList sll;
		
	}
	
}
