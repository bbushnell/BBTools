package icecream;

import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.SamReadStreamer;
import stream.SamStreamer;
import structures.ListNum;

/**
 * Wrapper for a ReadInputStream.
 * Produces one ZMW at a time for consumers.
 * Allows stopping after X reads or X ZMWs.
 * @author Brian Bushnell
 * @date June 5, 2020
 */
public class ZMWStreamer implements Runnable {
	
	public ZMWStreamer(FileFormat ff, int queuelen_, long maxReads_, long maxZMWs_){
		Data.USE_SAMBAMBA=false;//Sambamba changes PacBio headers.
		queuelen=Tools.mid(4, queuelen_, 64);
		maxReads=maxReads_;//(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		maxZMWs=maxZMWs_;
//		assert(false) : maxReads_+", "+maxReads;
		queue=new ArrayBlockingQueue<ZMW>(queuelen);
		if(ff.samOrBam() && useStreamer){
			cris=null;
			ss=makeStreamer(ff);
		}else{
			cris=makeCris(ff);
			ss=null;
		}
		assert((cris==null) != (ss==null)) : "Exactly one of cris or ss should exist.";
	}

	public ZMWStreamer(ConcurrentReadInputStream cris_, SamStreamer ss_, int queuelen_){
		cris=cris_;
		ss=ss_;
		queuelen=Tools.mid(4, queuelen_, 64);
		maxReads=-1;
		maxZMWs=-1;
		assert((cris==null) != (ss==null)) : "Exactly one of cris or ss should exist.";
		queue=new ArrayBlockingQueue<ZMW>(queuelen);
	}
	
	public Thread runStreamer(boolean makeThread){
		if(makeThread){
			Thread t=new Thread(this);
			t.start();
			return t;
		}else{
			run();
			return null;
		}
	}
	
	@Override
	public void run(){
		if(cris!=null){
			handleCris();
		}else{
			handleStreamer();
		}
	}
	
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
		cris.start(); //Start the stream
		if(verbose){System.err.println("Started cris");}
		return cris;
	}
	
	private SamReadStreamer makeStreamer(FileFormat ff){
		SamReadStreamer ss=new SamReadStreamer(ff, streamerThreads, true, maxReads);
		ss.start(); //Start the stream
		if(verbose){System.err.println("Started sam streamer");}
		return ss;
	}
	
	/** 
	 * Pull reads from the cris;
	 * organize them into lists of subreads from the same ZMW;
	 * put those lists into the shared queue.
	 */
	private void handleCris(){
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		ZMW buffer=new ZMW();buffer.id=ZMWs;
		long prevZmw=-1;

		long readsAdded=0;
//		long zmwsAdded=0;
		
		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){

			for(Read r : ln) {
				long zmw;
				try {
					zmw=Parse.parseZmw(r.id);
				} catch (Exception e) {
					zmw=r.numericID;//For testing only; disable for production
				}
				if(zmw<0){zmw=r.numericID;}//For testing only; disable for production
				if(verbose){System.err.println("Fetched read "+r.id+"; "+(zmw!=prevZmw)+", "+buffer.isEmpty()+", "+zmw+", "+prevZmw);}
				if(zmw!=prevZmw && !buffer.isEmpty()){
					ZMWs++;
					addToQueue(buffer);
					readsAdded+=buffer.size();
//					zmwsAdded++;
					buffer=new ZMW();buffer.id=ZMWs;
					if(maxZMWs>0 && ZMWs>=maxZMWs){break;}
				}
				buffer.add(r);
				prevZmw=zmw;
			}

			if(maxZMWs>0 && ZMWs>=maxZMWs){break;}
			cris.returnList(ln);
			
			//Fetch a new list
			ln=cris.nextList();
		}

		if(!buffer.isEmpty() && (maxZMWs<1 || ZMWs>=maxZMWs)){
			ZMWs++;
			readsAdded+=buffer.size();
			addToQueue(buffer);
		}
		
		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
//			cris.returnList(ln.id, true);
		}

		errorState|=ReadWrite.closeStreams(cris);
		addPoison();
	}
	
	/** 
	 * Pull reads from the streamer;
	 * organize them into lists of subreads from the same ZMW;
	 * put those lists into the shared queue.
	 */
	private void handleStreamer(){
		//Grab the first ListNum of reads
		ListNum<Read> ln=ss.nextList();

		ZMW buffer=new ZMW();buffer.id=ZMWs;
		long prevZmw=-1;
		
		long added=0;
		
		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){

			for(Read r : ln) {
				long zmw;
				try {
					zmw=Parse.parseZmw(r.id);
				} catch (Exception e) {
					zmw=r.numericID;//For testing only; disable for production
				}
				if(zmw<0){zmw=r.numericID;}//For testing only; disable for production
				if(verbose){System.err.println("Fetched read "+r.id+"; "+(zmw!=prevZmw)+", "+buffer.isEmpty()+", "+zmw+", "+prevZmw);}
				if(zmw!=prevZmw && !buffer.isEmpty()){
					ZMWs++;
					addToQueue(buffer);
					added+=buffer.size();
					buffer=new ZMW();buffer.id=ZMWs;
				}
				buffer.add(r);
				prevZmw=zmw;
			}
			
			//Fetch a new list
			ln=ss.nextList();
		}

		if(!buffer.isEmpty()){
			ZMWs++;
			added+=buffer.size();
			addToQueue(buffer);
		}
		
		addPoison();
	}
	
	private void addPoison(){
//		//Notify worker threads that there is no more data
//		for(int i=0; i<threads; i++){
//			addToQueue(POISON);
//		}
		addToQueue(POISON);
	}
	
	private void addToQueue(ZMW buffer){
		if(verbose) {System.err.println("Adding to queue "+(buffer==POISON ? "poison" : buffer.get(0).id));}
		while(buffer!=null) {
			try {
				queue.put(buffer);
				buffer=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public ZMW nextZMW(){
		ZMW buffer=null;
		while(buffer==null) {
			try {
				buffer=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose){System.err.println("Pulled from queue "+(buffer==POISON ? "poison" : buffer.get(0).id));}
		if(buffer==POISON){
			addToQueue(POISON);
			return null;
		}else{
			return buffer;
		}
	}
	
	private final ConcurrentReadInputStream cris;
	private final SamStreamer ss;
	private final int queuelen;
	public long ZMWs=0;
	private final long maxReads;
	private final long maxZMWs;
	public boolean errorState=false;
	
	private final ArrayBlockingQueue<ZMW> queue;
	private static final ZMW POISON=new ZMW(0);
	public static boolean verbose=false;
	
	//Streamer seems to give more highly variable timings... sometimes.  And it's not really needed.
	public static boolean useStreamer=false;
	//Only 1 thread for now to force ordered input
	public static final int streamerThreads=1;
	
}
