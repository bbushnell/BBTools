package icecream;

import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Single producer grouping consecutive reads by parsed numeric PacBio well ID.
 * The input must already keep each well contiguous. Parsed IDs distinguish movie
 * prefixes; unparsable or negative IDs retain the numericID fallback convention.
 * Publishes mutable Read references in ZMW lists through a bounded queue, then
 * a reusable terminal marker on normal completion. Multiple consumers may take
 * different wells; each owns the well it receives. Input errors and arbitrary
 * producer exceptions are not a general cancellation or terminal-delivery protocol.
 * @author Brian Bushnell, Sucrose
 * @date June 5, 2020
 */
public class ZMWStreamer implements Runnable{
	
	/** Creates and starts an input backend; grouping begins only when run is invoked.
	 * Also disables shared Sambamba input to preserve PacBio header names.
	 * @param ff Non-null input descriptor
	 * @param queuelen_ Requested well queue capacity, clamped to4-64
	 * @param maxReads_ Read limit forwarded to the selected backend
	 * @param maxZMWs_ Maximum published wells; nonpositive means unlimited
	 */
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

	/** Wraps exactly one caller-started backend, taking responsibility for normal closure.
	 * Does not start the supplied backend or impose new read/well limits.
	 * @param cris_ Started CRIS, or null when ss_ is used
	 * @param ss_ Started Streamer, or null when cris_ is used
	 * @param queuelen_ Requested well queue capacity, clamped to4-64
	 */
	public ZMWStreamer(ConcurrentReadInputStream cris_, Streamer ss_, int queuelen_){
		cris=cris_;
		ss=ss_;
		queuelen=Tools.mid(4, queuelen_, 64);
		maxReads=-1;
		maxZMWs=-1;
		assert((cris==null) != (ss==null)) : "Exactly one of cris or ss should exist.";
		queue=new ArrayBlockingQueue<ZMW>(queuelen);
	}
	
	/** Starts this single-use producer asynchronously or runs it on the calling thread.
	 * @param makeThread Whether to create/start a separate grouping thread
	 * @return Started thread for caller-managed joining, or null after synchronous execution
	 */
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
	
	/** Runs grouping once, closing the selected backend on its normal path.
	 * Does not reset counters or the queue and provides no general exception-to-terminal conversion.
	 */
	@Override
	public void run(){
		if(cris!=null){
			handleCris();
		}else{
			handleStreamer();
		}
	}
	
	/** Creates/starts the CRIS backend with the constructor's read limit. */
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
		cris.start(); //Start the stream
		if(verbose){System.err.println("Started cris");}
		return cris;
	}
	
	/** Creates/starts a SAM/BAM backend with header retention and current ordering request.
	 * The negative thread hint lets StreamerFactory choose its format defaults.
	 */
	private Streamer makeStreamer(FileFormat ff){
		Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreads, true, ordered, maxReads, true);
		ss.start(); //Start the stream
		if(verbose){System.err.println("Started sam streamer");}
		return ss;
	}
	
	/**
	 * Groups CRIS reads until EOF or the published-well cap. Returns each consumed
	 * backend list, closes input, records its close error and publishes the terminal.
	 * The backend read cap may end within a well; no whole-well read-limit guarantee.
	 */
	private void handleCris(){
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		ZMW buffer=new ZMW();buffer.id=ZMWs;
		long prevZmw=-1;
		String prevHeader=null;
		int prevMovieEnd=-1;

		long readsAdded=0;
//		long zmwsAdded=0;
		
		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){

			for(Read r : ln){
				long zmw;
				int movieEnd=-1;
				try{
					zmw=Parse.parseZmw(r.id);
					if(zmw>=0){movieEnd=r.id.indexOf('/');}
				}catch(Exception e){
					zmw=r.numericID;//For testing only; disable for production
				}
				if(zmw<0){zmw=r.numericID;}//For testing only; disable for production
				if(verbose){System.err.println("Fetched read "+r.id+"; "+(zmw!=prevZmw)+", "+buffer.isEmpty()+", "+zmw+", "+prevZmw);}
				final boolean newMovie=zmw==prevZmw && differentMovie(r.id, movieEnd, prevHeader, prevMovieEnd);
				if((zmw!=prevZmw || newMovie) && !buffer.isEmpty()){
					ZMWs++;
					addToQueue(buffer);
					readsAdded+=buffer.size();
//					zmwsAdded++;
					buffer=new ZMW();buffer.id=ZMWs;
					if(maxZMWs>0 && ZMWs>=maxZMWs){break;}
				}
				buffer.add(r);
				prevZmw=zmw;
				prevHeader=r.id;
				prevMovieEnd=movieEnd;
			}

			if(maxZMWs>0 && ZMWs>=maxZMWs){break;}
			cris.returnList(ln);
			
			//Fetch a new list
			ln=cris.nextList();
		}

		//PIC-100: EOF can arrive before the cap without another well boundary to flush this buffer.
		if(!buffer.isEmpty() && (maxZMWs<1 || ZMWs<maxZMWs)){
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
	 * Processes reads from a SAM/BAM Streamer. Groups subreads by ZMW ID,
	 * queues completed ZMW batches, forwards maxReads and enforces maxZMWs, then
	 * adds poison to signal completion after closing the streamer. At a well cap,
	 * drains remaining backend batches before closure: close alone does not guarantee
	 * that a producer blocked on its bounded output queue can finish. This limits
	 * published wells, without promising early physical input cancellation or a join.
	 */
	private void handleStreamer(){
		try{
			//Grab the first ListNum of reads
			ListNum<Read> ln=ss.nextList();

			ZMW buffer=new ZMW();buffer.id=ZMWs;
			long prevZmw=-1;
			String prevHeader=null;
			int prevMovieEnd=-1;

			long added=0;

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){

				for(Read r : ln){
					long zmw;
					int movieEnd=-1;
					try{
						zmw=Parse.parseZmw(r.id);
						if(zmw>=0){movieEnd=r.id.indexOf('/');}
					}catch(Exception e){
						zmw=r.numericID;//For testing only; disable for production
					}
					if(zmw<0){zmw=r.numericID;}//For testing only; disable for production
					if(verbose){System.err.println("Fetched read "+r.id+"; "+(zmw!=prevZmw)+", "+buffer.isEmpty()+", "+zmw+", "+prevZmw);}
					final boolean newMovie=zmw==prevZmw && differentMovie(r.id, movieEnd, prevHeader, prevMovieEnd);
					if((zmw!=prevZmw || newMovie) && !buffer.isEmpty()){
						ZMWs++;
						addToQueue(buffer);
						added+=buffer.size();
						buffer=new ZMW();buffer.id=ZMWs;
						if(maxZMWs>0 && ZMWs>=maxZMWs){break;}//PIC-101: Do not start another published well.
					}
					buffer.add(r);
					prevZmw=zmw;
					prevHeader=r.id;
					prevMovieEnd=movieEnd;
				}

				if(maxZMWs>0 && ZMWs>=maxZMWs){break;}
				//Fetch a new list
				ln=ss.nextList();
			}

			if(!buffer.isEmpty() && (maxZMWs<1 || ZMWs<maxZMWs)){
				ZMWs++;
				added+=buffer.size();
				addToQueue(buffer);
			}
			if(ln!=null && maxZMWs>0 && ZMWs>=maxZMWs){
				while(ss.nextList()!=null){
					//PIC-101: Release every queued backend batch through its normal terminal.
				}
			}
		}finally{
			ss.close();
			errorState|=ss.errorState();
		}
		
		addPoison();
	}
	
	/** Compares movie prefixes without extracting strings. An unknown movie preserves
	 * legacy numeric-ID fallback grouping rather than inventing an identity from malformed text.
	 */
	private static boolean differentMovie(String header, int movieEnd, String previous, int previousEnd){
		return movieEnd>=0 && previousEnd>=0 &&
				(movieEnd!=previousEnd || !header.regionMatches(0, previous, 0, movieEnd));
	}

	/** Publishes one shared terminal; consumers return it for the next consumer. */
	private void addPoison(){
//		//Notify worker threads that there is no more data
//		for(int i=0; i<threads; i++){
//			addToQueue(POISON);
//		}
		addToQueue(POISON);
	}
	
	/** Transfers a well reference, blocking for queue capacity and retrying interruptions. */
	private void addToQueue(ZMW buffer){
		if(verbose){System.err.println("Adding to queue "+(buffer==POISON ? "poison" : buffer.get(0).id));}
		while(buffer!=null){
			try{
				queue.put(buffer);
				buffer=null;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Blocks for a well, retrying interruptions; returned wells are not copied.
	 * @return Next owned well, or null after returning the shared terminal to the queue
	 */
	public ZMW nextZMW(){
		ZMW buffer=null;
		while(buffer==null){
			try{
				buffer=queue.take();
			}catch(InterruptedException e){
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
	
	/** Selected CRIS backend, mutually exclusive with ss. */
	private final ConcurrentReadInputStream cris;
	/** Streamer for SAM/BAM input (mutually exclusive with cris) */
	private final Streamer ss;
	/** Effective bounded-queue capacity. */
	private final int queuelen;
	/** Producer's published-well count; read final values after producer completion. */
	public long ZMWs=0;
	/** Backend read limit; may truncate a well. */
	private final long maxReads;
	/** Published-well cap, unlimited when nonpositive. */
	private final long maxZMWs;
	/** Collected normal close/backend errors; not a synchronized producer-completion flag. */
	public boolean errorState=false;
	/** Passed while constructing an owned Streamer; later writes do not reconfigure it. */
	public boolean ordered=true;
	
	/** Ownership handoff for completed well lists and the terminal marker. */
	private final ArrayBlockingQueue<ZMW> queue;
	/** Identity-tested terminal, returned to the queue by each consumer. */
	private static final ZMW POISON=new ZMW(0);
	/** Shared diagnostic output switch. */
	public static boolean verbose=false;
	
	//Streamer seems to give more highly variable timings... sometimes.  And it's not really needed.
	/** Selects the optional Streamer path for newly constructed SAM/BAM inputs. */
	public static boolean useStreamer=false;
	/** Negative hint selects backend defaults; ordering is requested separately. */
	public static final int streamerThreads=-1;
	
}
