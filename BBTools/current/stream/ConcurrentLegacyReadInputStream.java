package stream;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import structures.ListNum;

public class ConcurrentLegacyReadInputStream extends ConcurrentReadInputStream {
	
	public ConcurrentLegacyReadInputStream(ReadInputStream source, long maxReadsToGenerate){
		super(source.fname());
		producer=source;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
//		if(maxReads<Long.MAX_VALUE){System.err.println("maxReads="+maxReads);}
	}
	
	@Override
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		while(list==null){
			try {
				list=depot.full.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				if(shutdown){return null;}
			}
		}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("cris_:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("cris_:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	@Override
	public void run() {
//		producer.start();
		threads=new Thread[] {Thread.currentThread()};

		if(producer.preferLists()){
			readLists();
		}else{
			readSingles();
		}
		
		addPoison();
		
		//End thread
		
		while(!depot.empty.isEmpty()){
			depot.full.add(depot.empty.poll());
		}
//		System.err.println(depot.full.size()+", "+depot.empty.size());
	}
	
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
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
			if(list!=null){depot.full.add(list);}
		}
		//System.err.println("Added poison.");
	}
	
	private final void readSingles(){
		
		long bases=0;
		while(!shutdown && producer.hasMore() && generated<maxReads && bases<MAX_DATA){
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
			
			for(int i=0; i<depot.bufferSize && generated<maxReads && bases<MAX_DATA; i++){
				Read r=producer.next();
				if(r==null){break;}
				list.add(r);
				bases+=r.length();
				bases+=(r.mate==null || r.mate.bases==null ? 0 : r.mateLength());
				generated++;
			}
			depot.full.add(list);
		}
	}
	
	private final void readLists(){
		
		ArrayList<Read> buffer=null;
		ArrayList<Read> list=null;
		int next=0;
		while(buffer!=null || (!shutdown && producer.hasMore() && generated<maxReads)){
			while(list==null){
				//System.err.println("Fetching a list: generated="+generated+"/"+maxReads);
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
				//System.err.println("Fetched");
			}
			if(shutdown || list==null){
				//System.err.println("Shutdown triggered; breaking.");
				break;
			}
			
			long bases=0;
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
				if(buffer==null || next>=buffer.size()){
					buffer=producer.nextList();
					next=0;
				}
				if(buffer==null){break;}
				assert(buffer.size()<=BUF_LEN); //Although this is not really necessary.
				
				if(buffer.size()<=(BUF_LEN-list.size()) && (buffer.size()+generated)<maxReads && randy==null){
					//Then do a quicker bulk operation
					list.addAll(buffer);
					for(Read a : buffer){
						readsIn++;
						basesIn+=a.length();
						bases+=a.length();
						if(a.mate!=null){
							readsIn++;
							basesIn+=a.mateLength();
							bases+=a.mateLength();
						}
					}
					generated+=buffer.size();
					next=0;
					buffer=null;
				}else{
					while(next<buffer.size() && list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
						Read r=buffer.get(next);
						readsIn++;
						basesIn+=r.length();
						if(r.mate!=null){
							readsIn++;
							basesIn+=r.mateLength();
						}
						if(randy==null || randy.nextFloat()<samplerate){
							list.add(r);
							bases+=r.length();
							bases+=(r.mate==null || r.mate.bases==null ? 0 : r.mateLength());
						}
						generated++;
//						if(generated>1 && (generated%1000000)==0){System.err.println("Generated read #"+generated);}
						next++;
					}
					
					if(next>=buffer.size()){
						buffer=null;
						next=0;
					}
				}
			}
			//System.err.println("Adding list to full depot.");
			depot.full.add(list);
			//System.err.println("Added.");
			list=null;
		}

	}
	
	private boolean shutdown=false;
	
	@Override
	public void shutdown(){
		shutdown=true;
		if(threads[0]!=null && threads[0].isAlive()){
			threads[0].interrupt();
		}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		producer.restart();
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
	}
	
	@Override
	public synchronized void close(){
//		System.err.println("Closing cris: "+maxReads+", "+generated);
//		if(threads!=null){
//			for(int i=0; i<threads.length; i++){
//				if(threads[i]!=null){System.err.println(i+": "+threads[i].isAlive());}
//			}
//		}
		producer.close();
	}

	@Override
	public boolean paired() {
		return producer.paired();
	}
	
	@Override
	public boolean verbose(){return verbose;}
	
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
	public boolean errorState(){return errorState || (producer!=null && producer.errorState());}
	/** TODO */
	private boolean errorState=false;
	
	private float samplerate=1f;
	private java.util.Random randy=null;
	
	@Override
	public Object[] producers(){return new Object[] {producer};}

	private Thread[] threads;

	public final ReadInputStream producer;
	private ConcurrentDepot<Read> depot;
	
	public static boolean verbose=false;
	
	private long basesIn=0;
	private long readsIn=0;
	
	private long maxReads;
	private long generated=0;
	private long listnum=0;
	

}
