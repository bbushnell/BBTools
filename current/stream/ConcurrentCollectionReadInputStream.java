package stream;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import dna.Data;
import structures.ListNum;

public class ConcurrentCollectionReadInputStream extends ConcurrentReadInputStream {
	
	public ConcurrentCollectionReadInputStream(List<Read> source1, List<Read> source2, long maxReadsToGenerate){
		super("list");
		assert(source1!=source2);
		producer1=source1;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		producer2=source2;
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
		
	}
	
	@Override
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("**************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("**************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try {
				list=depot.full.take();
				assert(list!=null);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("**************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisC:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisC:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	@Override
	public void run() {
//		producer.start();
		threads=new Thread[] {Thread.currentThread()};
		if(verbose){System.err.println("crisC started, thread="+threads[0]);}

//		readLists();
		readSingles();

		addPoison();
		
		//End thread

		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("Ending");
			if(verbose){System.err.println("B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
//		System.err.println("cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());
	}
	
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
//					System.err.println("Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("D: Adding list("+list.size()+") to full "+depot.full.size()+"/"+depot.bufferCount);}
				depot.full.add(list);
			}
		}
		//System.err.println("Added poison.");
	}
	
	private final void readSingles(){

		for(int i=0; !shutdown && i<producer1.size() && generated<maxReads; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
			}
			if(shutdown || list==null){break;}
			
			long bases=0;
			final long lim=producer1.size();
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA && generated<lim){
				Read a=producer1.get((int)generated);
				Read b=(producer2==null ? null : producer2.get((int)generated));
				if(a==null){break;}
				readsIn++;
				basesIn+=a.length();
				if(b!=null){
					readsIn++;
					basesIn+=b.length();
				}
				if(randy==null || randy.nextFloat()<samplerate){
					list.add(a);
					if(b!=null){
						assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
						a.mate=b;
						b.mate=a;

						assert(a.pairnum()==0);
						b.setPairnum(1);
						bases+=(b.bases==null ? 0 : b.length());
					}
					bases+=(a.bases==null ? 0 : a.length());
				}
				incrementGenerated(1);
			}

			if(verbose){System.err.println("E: Adding list("+list.size()+") to full "+depot.full.size()+"/"+depot.bufferCount);}
			depot.full.add(list);
		}
	}
	
	private boolean shutdown=false;
	
	@Override
	public void shutdown(){
		if(verbose){System.out.println("Called shutdown.");}
		shutdown=true;
		if(!shutdown){
			if(verbose){System.out.println("shutdown 2.");}
			for(Thread t : threads){
				if(verbose){System.out.println("shutdown 3.");}
				if(t!=null && t.isAlive()){
					if(verbose){System.out.println("shutdown 4.");}
					t.interrupt();
					if(verbose){System.out.println("shutdown 5.");}
				}
			}
		}
		if(verbose){System.out.println("shutdown 6.");}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
		nextProgress=PROGRESS_INCR;
	}
	
	@Override
	public synchronized void close(){
		if(verbose){System.out.println("Thread "+Thread.currentThread().getId()+" called close.");}
		shutdown();
//		producer1.close();
//		if(producer2!=null){producer2.close();}
//		System.out.println("A");
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			if(verbose){System.out.println("close 1.");}
			
			while(threads[0].isAlive()){
				if(verbose){System.out.println("close 2: Thread "+Thread.currentThread().getId()+" closing thread "+threads[0].getId()+" "+threads[0].getState());}
//				System.out.println("B");
				ArrayList<Read> list=null;
				for(int i=0; i<1 && list==null && threads[0].isAlive(); i++){
					if(verbose){System.out.println("close 3.");}
					try {
						if(verbose){System.out.println("close 4.");}
						list=depot.full.poll(100, TimeUnit.MILLISECONDS);
						if(verbose){System.out.println("close 5; list.size()="+depot.full.size()+", list="+(list==null ? "null" : list.size()+""));}
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						System.err.println("Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
				if(verbose){System.out.println("close 6.");}
				
//				System.out.println("isAlive? "+threads[0].isAlive());
			}
			if(verbose){System.out.println("close 7.");}
			
		}
		if(verbose){System.out.println("close 8.");}
		
		if(threads!=null){
			if(verbose){System.out.println("close 9.");}
			for(int i=1; i<threads.length; i++){
				if(verbose){System.out.println("close 10.");}
				while(threads[i]!=null && threads[i].isAlive()){
					if(verbose){System.out.println("close 11.");}
					try {
						if(verbose){System.out.println("close 12.");}
						threads[i].join();
						if(verbose){System.out.println("close 13.");}
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		if(verbose){System.out.println("close 14.");}
		
	}

	@Override
	public boolean paired() {
		return producer2!=null ? true : (producer1==null || producer1.isEmpty() ? false : producer1.get(0).mate!=null);
	}
	
	@Override
	public boolean verbose(){return verbose;}
	
	private void incrementGenerated(long amt){
		generated+=amt;
		if(SHOW_PROGRESS && generated>=nextProgress){
			Data.sysout.print('.');
			nextProgress+=PROGRESS_INCR;
		}
	}
	
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		if(rate>=1f){
			randy=null;
		}else if(seed>-1){
			randy=new java.util.Random(seed);
		}else{
			randy=new java.util.Random();
		}
	}
	
	@Override
	public long basesIn(){return basesIn;}
	@Override
	public long readsIn(){return readsIn;}
	
	@Override
	public boolean errorState(){return errorState;}
	/** TODO */
	private boolean errorState=false;
	
	private float samplerate=1f;
	private java.util.Random randy=null;
	
	private Thread[] threads;
	
	@Override
	public Object[] producers(){return new Object[] {producer1, producer2};}
	
	public final List<Read> producer1;
	public final List<Read> producer2;
	private ConcurrentDepot<Read> depot;
	
	private long basesIn=0;
	private long readsIn=0;
	
	private long maxReads;
	private long generated=0;
	private long listnum=0;
	private long nextProgress=PROGRESS_INCR;
	
	public static boolean verbose=false;
	
	private static final ArrayList<Read> poison=new ArrayList<Read>(0);
	
}
