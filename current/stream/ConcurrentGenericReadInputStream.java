package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import dna.Data;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

public class ConcurrentGenericReadInputStream extends ConcurrentReadInputStream {
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		String in1=args[0];
		String in2=(args.length<2 || args[1].equalsIgnoreCase("null") || args[1].contains("=") ? null : args[1]);
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		long maxReads=-1;
		for(int i=1; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		Parser.processQuality();
		
		assert(FastaReadInputStream.settingsOK());
		Timer t=new Timer();
		
		ConcurrentReadInputStream cris=getReadInputStream(maxReads, false, true, in1, in2);
		System.out.println("Fetched "+cris.getClass().getName());
		{
			Object[] p=cris.producers();
//			while(p[0]==null){
//				p=cris.producers();
//			}
			System.out.print("Producers: ");
			String comma="";
			for(Object o : p){
				System.out.print(comma+(o==null ? "null" : o.getClass().getName()));
				comma=", ";
			}
			System.out.println();
		}
		boolean paired=cris.paired();
		System.out.println("paired="+paired);
		cris.start(); //4567
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((r.mate!=null)==paired);
		}
		
		long readCount=0;
		long baseCount=0;
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			
			for(Read r : reads){
				Read r2=r.mate;
				if(r!=null){
					readCount++;
					if(r.bases!=null){
						baseCount+=r.length();
					}
				}
				if(r2!=null){
					readCount++;
					if(r2.bases!=null){
						baseCount+=r2.length();
					}
				}
			}
			cris.returnList(ln);
//			System.err.println("fetching list");
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
//			System.out.println("reads: "+(reads==null ? "null" : reads.size()));
		}
		System.err.println("Finished reading");
		cris.returnList(ln);
		
		cris.close();
		t.stop();

		System.out.println("Reads:      \t"+readCount);
		System.out.println("Bases:      \t"+baseCount);
		System.out.println("Avg Length: \t"+Tools.format("%.2f",baseCount*1.0/readCount));
		System.out.println("Time:      \t"+t);
	}
	
	public ConcurrentGenericReadInputStream(ReadInputStream source1, ReadInputStream source2, long maxReadsToGenerate){
		super((source1==null ? "null" : source1.fname())+","+(source2==null ? "null" : source2.fname()));
		assert(source1!=source2);
		producer1=source1;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
//		assert(false) : BUF_LEN+", "+NUM_BUFFS;
		producer2=source2;
		assert(source2==null || !FASTQ.FORCE_INTERLEAVED) : "Please do not set 'interleaved=true' with dual input files.";
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("crisG:    Warning - created a read stream for 0 reads.");
			assert(false);
		}
//		if(maxReads<Long.MAX_VALUE){System.err.println("crisG:    maxReads="+maxReads);}

		if(producer1!=null){p1q=new ArrayBlockingQueue<ArrayList<Read>>(4);}
		if(producer2!=null){p2q=new ArrayBlockingQueue<ArrayList<Read>>(4);}
	}
	
	@Override
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("crisG:    **************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("crisG:    **************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
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

		if(verbose){System.err.println("crisG:    **************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	@Override
	public void returnList(long listNumber, boolean poison){
		if(poison){
			if(verbose){System.err.println("crisG:    A: Adding empty list to full.");}
			depot.full.add(new ArrayList<Read>(0));
		}else{
			if(verbose){System.err.println("crisG:    A: Adding empty list to empty.");}
			depot.empty.add(new ArrayList<Read>(BUF_LEN));
		}
	}
	
	@Override
	public void run() {
		try {
			run0();
		} catch (AssertionError e) {
			KillSwitch.assertionKill(e);
		}
	}
	
	private void run0() {
//		producer.start();
		synchronized(running){
			assert(!running[0]) : "This cris was started by multiple threads.";
			running[0]=true;
		}

		ReadThread rt1=null;
		ReadThread rt2=null;
		rt1=new ReadThread(producer1, p1q);
		rt2=(producer2==null ? null : new ReadThread(producer2, p2q));
		rt1.start();
		if(rt2!=null){rt2.start();}
		
		threads=(rt1==null ? new Thread[] {Thread.currentThread()} :
			rt2==null ? new Thread[] {Thread.currentThread(), rt1} :
				new Thread[] {Thread.currentThread(), rt1, rt2});

		try {
			readLists();
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
//		readSingles();

		addPoison();
		
		//End thread

		if(verbose){System.err.println("crisG:    cris finished addPoison.");}
		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("crisG:    Ending");
			if(verbose){System.err.println("crisG:    B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
		if(verbose){System.err.println("crisG:    cris thread syncing before shutdown.");}
		
		synchronized(running){//TODO Note: for some reason syncing on 'this' instead of 'running' causes a hang.  Something else must be syncing improperly on this.
			assert(running[0]);
			running[0]=false;
		}
		if(verbose){System.err.println("crisG:    cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());}
	}
	
	private final void addPoison(){
		//System.err.println("crisG:    Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("crisG:    C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
//					System.err.println("crisG:    Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("crisG:    D: Adding list("+list.size()+") to full.");}
				depot.full.add(list);
			}
		}
		if(verbose){System.err.println("crisG:    Added poison.");}
	}
	
//	@Deprecated
//	private final void readSingles(){
//
//		while(!shutdown && producer1.hasMore() && generated<maxReads){
//			ArrayList<Read> list=null;
//			while(list==null){
//				try {
//					list=depot.empty.take();
//				} catch (InterruptedException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//					if(shutdown){break;}
//				}
//			}
//			if(shutdown || list==null){break;}
//			
//			long bases=0;
//			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
//				Read a=producer1.next();
//				Read b=(producer2==null ? null : producer2.next());
//				if(a==null){break;}
//				readsIn++;
//				basesIn+=a.length();
//				if(b!=null){
//					readsIn++;
//					basesIn+=b.length();
//				}
//				if(randy==null || randy.nextFloat()<samplerate){
//					list.add(a);
//					if(b!=null){
//						assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
//						assert(a.mate==null) : "Please set interleaved=false when using dual input files.\n"+a.id+"\n"+a.mate.id+"\n"+producer1+"\n"+producer2;
//						assert(b.mate==null) : "Please set interleaved=false when using dual input files.";
//						a.mate=b;
//						b.mate=a;
//
//						assert(a.pairnum()==0);
//						b.setPairnum(1);
//						bases+=(b.bases==null ? 0 : b.length());
//					}
//					bases+=(a.bases==null ? 0 : a.length());
//				}
//				incrementGenerated(1);
//			}
//
//			if(verbose){System.err.println("crisG:    E: Adding list("+list.size()+") to full.");}
//			depot.full.add(list);
//		}
//	}
	
	private final void readLists(){
		ArrayList<Read> buffer1=null;
		ArrayList<Read> buffer2=null;
		ArrayList<Read> list=null;
		int next=0;
		
//		System.out.println("crisG:    a");
		if(verbose){System.err.println(getClass().getName()+" entering read lists loop.");}
		while(buffer1!=poison && (buffer1!=null || (!shutdown && generated<maxReads))){
//			System.out.println("crisG:    b");
			if(verbose){System.err.println("crisG:    looping: buffer1==null "+(buffer1==null)+", buffer1==poison "+(buffer1==poison)
					+", shutdown="+shutdown+", generated<maxReads="+(generated<maxReads));}
			while(list==null){
				if(verbose){System.err.println("crisG:    Fetching an empty list: generated="+generated+"/"+maxReads);}
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
				if(verbose){System.err.println("crisG:    Fetched "+(list==null ? "null" : ""+list.size()));}
			}
//			System.out.println("crisG:    c");
			if(verbose){System.err.println("crisG:    Left empty fetch loop.");}
			if(shutdown || list==null){
				//System.err.println("crisG:    Shutdown triggered; breaking.");
				break;
			}
//			System.out.println("crisG:    d");
			
			if(verbose){System.err.println("crisG:    Entering full fetch loop.");}
			long bases=0;
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
				if(verbose){System.err.println("crisG:    list.size()="+list.size()+", depot.bufferSize="+depot.bufferSize+", generated="+generated);}
				if(buffer1==null || next>=buffer1.size()){
					buffer1=null;
					while(!shutdown && buffer1==null){
						try {
							buffer1=p1q.take();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
//					System.out.println("crisG:    e");
					
					if(buffer1!=null && p2q!=null){
						buffer2=null;
						while(!shutdown && buffer2==null){
							try {
								buffer2=p2q.take();
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						if(buffer2!=null && FASTQ.PAIR_READS){pair(buffer1, buffer2);}
						if(REMOVE_DISCARDED_READS){removeDiscarded(buffer1, buffer2);}
					}
//					System.out.println("crisG:    f");
					next=0;
				}
//				System.out.println("crisG:    g");
				if(buffer1==null || buffer1==poison || shutdown){
//					if(list!=null && list.size()>0){
//						if(verbose){System.err.println("crisG:    G: Adding list("+list.size()+") to full.");}
//						depot.full.add(list);
//						list=null;
//					}
					if(verbose){System.err.println("crisG:    Breaking because buffer1==null: "+(buffer1==null)+" || buffer1==poison: "+(buffer1==poison)+" || shutdown: "+shutdown);}
					break;
				}
				assert(buffer1.size()<=BUF_LEN); //Although this is not really necessary.
				
//				assert(!set2.contains(buffer1)) : buffer1.hashCode();
//				set2.add(buffer1);
//				System.out.println(buffer1.hashCode());
				
				if(buffer2!=null){
//					System.out.println("crisG:    h");
					
					if(buffer2!=null && (buffer1==null || buffer2.size()!=buffer1.size()) && !ALLOW_UNEQUAL_LENGTHS){
						System.err.println("crisG:    Error: Misaligned read streams.");
						errorState=true;
						return;
					}
					assert(ALLOW_UNEQUAL_LENGTHS || buffer2==null || buffer2.size()==buffer1.size());
				}
				
				//Code disabled because it does not actually seem to make anything faster.
//				if(buffer1.size()<=(BUF_LEN-list.size()) && (buffer1.size()+generated)<maxReads && randy==null){
//					//System.out.println("crisG:    j");
//					//Then do a quicker bulk operation
//
//					for(Read a : buffer1){
//						list.add(a);
//						Read b=a.mate;
//						readsIn++;
//						basesIn+=a.length();
//						bases+=(a.bases==null ? 0 : a.length());
////						bases+=(b==null || b.bases==null ? 0 : b.length());
//						if(b!=null){
//							readsIn++;
//							basesIn+=b.length();
//							bases+=b.length();
//							assert(a.pairnum()==0 && b.pairnum()==1);
//						}
////						System.out.println(generated+", "+readsIn+", "+(b==null));
//					}
////					list.addAll(buffer1); //This is actually slower due to an array clone operation.
//					incrementGenerated(buffer1.size());
//
//					next=0;
//					buffer1=null;
//					buffer2=null;
//				}else
				{

					while(next<buffer1.size() && list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
						Read a=buffer1.get(next);
						Read b=a.mate;
						readsIn++;
						basesIn+=a.length();
						if(b!=null){
							readsIn++;
							basesIn+=b.length();
//							assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
//							a.mate=b;
//							b.mate=a;
//
//							assert(a.pairnum()==0);
//							b.setPairnum(1);
							assert(a.pairnum()==0 && b.pairnum()==1 && a.mate==b && b.mate==a && a.numericID==b.numericID) :
								"There is something wrong with the read pairing.\n"+
								a.pairnum()+", "+(b.pairnum())+", "+(a.mate==b)+", "+(b.mate==a)+", "+(a.numericID)+", "+(b.numericID);
						}
						if(randy==null || randy.nextFloat()<samplerate){
							list.add(a);
							bases+=a.length();
							if(a.mate!=null){
								bases+=a.mateLength();
//								assert(a.pairnum()==0 && a.mate.pairnum()==1);
							}
						}
						incrementGenerated(1);
						next++;
					}
//					System.out.println("crisG:    l");

					
					if(next>=buffer1.size()){
						buffer1=null;
						buffer2=null;
						next=0;
//						System.out.println("crisG:    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
					}else{
//						System.out.println("crisG:    ------------------------------------------------");
					}
//					System.out.println("crisG:    m");
				}
				if(verbose){System.err.println("crisG:    Loop end: list.size()="+(list.size()+", depot.bufferSize="+depot.bufferSize+", generated="+generated));}
//				System.out.println("crisG:    n");
				if(verbose){System.err.println(Thread.currentThread().getName());}
			}
			
//			System.out.println("crisG:    p");
//			System.err.println("crisG:    Adding list to full depot.  Shutdown="+shutdown);
			if(verbose){System.err.println("crisG:    F: Adding list("+list.size()+") to full.");}
			depot.full.add(list);
			list=null;
//			System.err.println("crisG:    Added.");
			
//			System.out.println("crisG:    o");
			if(buffer1==poison){
				if(verbose){System.err.println("crisG:    Detected poison from buffer1.");}
				break;
			}
			if(verbose){System.err.println("crisG:    Finished loop iteration.\n");}
			if(verbose){System.err.println("crisG:    loop end: buffer1==null "+(buffer1==null)+", buffer1==poison "+(buffer1==poison)
					+", shutdown="+shutdown+", generated<maxReads="+(generated<maxReads));}
//			System.out.println("crisG:    q");
		}
//		System.out.println("crisG:    r");
		
		
		p1q.clear();
		if(p2q!=null){p2q.clear();}
	}
	
	private final void pair(ArrayList<Read> buffer1, ArrayList<Read> buffer2){
		final int len1=buffer1.size(), len2=buffer2.size();
		assert(ALLOW_UNEQUAL_LENGTHS || len1==len2) : "\nThere appear to be different numbers of reads in the paired input files." +
				"\nThe pairing may have been corrupted by an upstream process.  It may be fixable by running repair.sh.";
		final int lim=Tools.min(len1, len2);
		
		if(FASTQ.FLIP_R2){
			for(Read r2 : buffer2) {r2.reverseComplement();}
		}
		
		for(int i=0; i<lim; i++){
			Read a=buffer1.get(i);
			Read b=buffer2.get(i);

			assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
			assert(a.mate==null) : "Please set interleaved=false when using dual input files.\n"+a.id+"\n"+a.mate.id+"\n"+b.id+"\n"+producer1+"\n"+producer2;
			assert(b.mate==null) : "Please set interleaved=false when using dual input files.";
			a.mate=b;
			b.mate=a;

			assert(a.pairnum()==0);
			b.setPairnum(1);
			//		assert(a.pairnum()!=b.pairnum());
		}
		
		if(len1>len2){
			//do nothing;
		}else if(len2>len1){
			for(int i=lim; i<len2; i++){
				Read b=buffer2.get(i);
				b.setPairnum(0);
				buffer1.add(b);
			}
		}
	}
	
	private static final int removeDiscarded(ArrayList<Read> buffer1, ArrayList<Read> buffer2){
		int removed=0;
		if(buffer2==null){
			for(int i=0; i<buffer1.size(); i++){
				Read a=buffer1.get(i);
				if(a.discarded()){
					buffer1.set(i, null);
					removed++;
				}
			}
		}else{
			for(int i=0; i<buffer1.size(); i++){
				Read a=buffer1.get(i);
				Read b=buffer2.get(i);
				if(a.discarded() || b.discarded()){
					buffer1.set(i, null);
					buffer2.set(i, null);
					removed++;
				}
			}
		}
		if(removed>0){
			Tools.condenseStrict(buffer1);
			if(buffer2!=null){Tools.condenseStrict(buffer2);}
		}
		return removed;
	}
	
	private boolean shutdown=false;
	
	@Override
	public void shutdown(){
//		System.err.println("crisG:    Called shutdown.");
		shutdown=true;
		if(!shutdown){//???
			for(Thread t : threads){
				if(t!=null && t.isAlive()){
					t.interrupt();
				}
			}
		}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		p1q.clear();
		if(p2q!=null){p2q.clear();}
		producer1.restart();
		if(producer2!=null){producer2.restart();}
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		basesIn=0;
		readsIn=0;
		listnum=0; //Added Oct 9, 2014
		nextProgress=PROGRESS_INCR;
		lastTime=System.nanoTime();
	}
	
	@Override
	public synchronized void close(){
		if(verbose){System.err.println("crisG:    Called shutdown for "+producer1+"; "+threads[0].getState());}
//		if(verbose){System.err.println(((FastqReadInputStream)producer1).tf.isOpen());}
		shutdown();
		errorState|=producer1.close();
		if(producer2!=null){errorState|=producer2.close();}
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			
			while(threads[0].isAlive()){
//				System.out.println("crisG:    B");
				ArrayList<Read> list=null;
				for(int i=0; i<1000 && list==null && threads[0].isAlive(); i++){
					try {
						list=depot.full.poll(200, TimeUnit.MILLISECONDS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						System.err.println("crisG:    Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
			}
			
		}
		
		if(threads!=null){
			for(int i=1; i<threads.length; i++){
				while(threads[i]!=null && threads[i].getState()!=Thread.State.TERMINATED){
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		//This assertion should be impossible but somehow it fired.
		//However, by the time the warning was issued, "isAlive" returned false anyway.
//		assert(threads==null || threads.length<2 || threads[1]==null || !threads[1].isAlive()) : ((ReadThread)threads[1]).generatedLocal+", "+threads[1].isAlive()+", "+threads[1].getState();
//		threads=null;
//		System.out.println("crisG:    C");

		if(verbose){System.err.println("crisG:    shutdown exited; errorState="+errorState);}
	}

	@Override
	public boolean paired() {
		return FASTQ.PAIR_READS && (producer1.paired() || producer2!=null);
	}
	
	@Override
	public boolean verbose(){return verbose;}
	
	private class ReadThread extends Thread{
		ReadThread(ReadInputStream producer_, ArrayBlockingQueue<ArrayList<Read>> pq_){
			producer=producer_;
			pq=pq_;
		}
		
		@Override
		public void run(){
			try {
				readLists();
			} catch (AssertionError e) {
				KillSwitch.assertionKill(e);
			}
		}
		
		private final void readLists(){
			
			ArrayList<Read> list=null;
			
			if(verbose){System.err.println(getClass().getName()+" entering read lists loop.");}
			while(list!=null || (!shutdown && producer.hasMore() && generatedLocal<maxReads)){

				if(verbose){System.err.println(getClass().getName()+" looping: buffer1==null "+(list==null)+", shutdown="+shutdown+
						", producer.hasMore()="+producer.hasMore()+", generated<maxReads="+(generatedLocal<maxReads));}

				
				
				if(verbose){System.err.println(getClass().getName()+" Entering full fetch loop.");}
				while(generatedLocal<maxReads){
//					System.out.println("crisG:    E");
					if(verbose){System.err.println(getClass().getName()+" depot.bufferSize="+depot.bufferSize+", generated="+generatedLocal);}
//					System.out.println("crisG:    F");
					try {
						list=producer.nextList();
					} catch (OutOfMemoryError e){
						KillSwitch.memKill(e);
					} catch (Throwable e1) {
						// TODO
//						System.err.print('*');
						e1.printStackTrace();
						list=null;
						shutdown=true;
						try {
							pq.put(new ArrayList<Read>(1));
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						errorState=true;
					}
					if(verbose){System.err.println(getClass().getName()+" grabbed a list of size "+(list==null ? "null" : list.size()+""));}
//					System.out.println("crisG:    G");
					if(list==null){
//						System.out.println("crisG:    H");
						if(verbose){System.err.println(getClass().getName()+" broke loop on null list.");}
						break;
					}
					assert(list.size()>0);
					assert(list.size()<=BUF_LEN); //Although this is not really necessary.
//					System.out.println("crisG:    I");
					if(list.size()+generatedLocal>maxReads){
//						System.out.println("crisG:    J");
						if(verbose){System.err.println("crisG:    Removing extra reads.");}
						while(list.size()+generatedLocal>maxReads){list.remove(list.size()-1);}
//						System.out.println("crisG:    K");
					}
//					System.out.println("crisG:    A");
					while(list!=null && !shutdown){
//						System.out.println("crisG:    B");
						try {
							if(verbose){System.err.println("crisG:    Trying to add list");}
							pq.put(list);
							generatedLocal+=list.size();
							list=null;
							if(verbose){
								System.out.println("crisG:    Added list; pq.size() = "+pq.size());
							}
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
//						System.out.println("crisG:    C");
					}
//					System.out.println("crisG:    D");
					if(verbose){System.err.println("crisG:    looping");}
				}

				if(verbose){System.err.println(getClass().getName()+" Finished inner loop iteration.\n");}
			}
			

			if(verbose){System.err.println(getClass().getName()+" attempting to poison output queue.");}
			boolean b=true;
			while(b){
				//TODO Note that this could cause a deadlock if there was a premature shutdown, so the consumer died while the queue was full.
				try {
//					pq.offer(poison, 10000, TimeUnit.SECONDS);
					pq.put(poison);
					b=false;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			

			if(verbose){System.err.println(getClass().getName()+" exited read lists loop: "+(list==null)+", "+shutdown+", "+producer.hasMore()+", "+generatedLocal+", "+maxReads);}

		}
		
		private final ArrayBlockingQueue<ArrayList<Read>> pq;
		private final ReadInputStream producer;
		private long generatedLocal=0;
	}
	
	private void incrementGenerated(long amt){
		generated+=amt;
		if(SHOW_PROGRESS && generated>=nextProgress){
			if(SHOW_PROGRESS2){
				nextProgress+=PROGRESS_INCR;
				long x=System.nanoTime();
				long duration=x-lastTime;
				lastTime=x;
				Data.sysout.println(Tools.format("%.1f", duration*0.000000001));
//				Data.sysout.println((long)(0.5+duration*0.000000001)+" ");
			}else{
				nextProgress+=PROGRESS_INCR;
				Data.sysout.print('.');
			}
		}
//		System.err.println("crisG:    generated="+generated+"\treadsIn="+readsIn);
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
	public boolean errorState(){
		return errorState ||
			(producer1==null ? false : producer1.errorState()) || (producer2==null ? false : producer2.errorState());}
	/** TODO */
	private boolean errorState=false;
	
	private boolean[] running=new boolean[] {false};
	
	private float samplerate=1f;
	private java.util.Random randy=null;
	
	private ArrayBlockingQueue<ArrayList<Read>> p1q;
	private ArrayBlockingQueue<ArrayList<Read>> p2q;
	
	
	@Override
	public Object[] producers(){return producer2==null ? new Object[] {producer1} : new Object[] {producer1, producer2};}

	private Thread[] threads;
	
	public final ReadInputStream producer1;
	public final ReadInputStream producer2;
	private ConcurrentDepot<Read> depot;
	
	private long basesIn=0;
	private long readsIn=0;
	
	private long maxReads;
	private long generated=0;
	private long listnum=0;
	private long nextProgress=PROGRESS_INCR;
	private long lastTime=System.nanoTime();
	
	public static boolean verbose=false;
	
	private static final ArrayList<Read> poison=new ArrayList<Read>(0);
	
}
