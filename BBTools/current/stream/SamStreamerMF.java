package stream;

import java.io.PrintStream;
import java.util.ArrayDeque;

import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Loads multiple sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date March 6, 2019
 *
 */
public class SamStreamerMF {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static final void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		int threads=Shared.threads();
		if(args.length>1){threads=Integer.parseInt(args[1]);}
		SamStreamerMF x=new SamStreamerMF(args[0].split(","), threads, false, -1);
		
		//Run the object
		x.start();
		x.test();
		
		t.stop("Time: ");
	}
	
	/**
	 * Constructor.
	 */
	public SamStreamerMF(String[] fnames_, int threads_, boolean saveHeader_, long maxReads_){
		this(FileFormat.testInput(fnames_, FileFormat.SAM, null, true, false), threads_, saveHeader_, maxReads_);
	}
	
	/**
	 * Constructor.
	 */
	public SamStreamerMF(FileFormat[] ffin_, int threads_, boolean saveHeader_, long maxReads_){
		fname=ffin_[0].name();
		threads=threads_;
		ffin=ffin_;
		saveHeader=saveHeader_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	
	final void test(){
		for(ListNum<Read> list=nextReads(); list!=null; list=nextReads()){
			if(verbose){outstream.println("Got list of size "+list.size());}
		}
	}
	
	
	/** Create read streams and process all data */
	public final void start(){
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads();
		
		if(verbose){outstream.println("Finished; closing streams.");}
	}

	public final ListNum<Read> nextList(){return nextReads();}
	public final ListNum<Read> nextReads(){
		ListNum<Read> list=null;
		assert(activeStreamers!=null);
		synchronized(activeStreamers){
			if(activeStreamers.isEmpty()){return null;}
			while(list==null && !activeStreamers.isEmpty()){
				SamReadStreamer srs=activeStreamers.poll();
				list=srs.nextReads();
				if(list!=null){activeStreamers.add(srs);}
				else{
					readsProcessed+=srs.readsProcessed;
					basesProcessed+=srs.basesProcessed;
					if(srs.header!=null){
						SamReadInputStream.setSharedHeader(srs.header);
					}
					
					if(!streamerSource.isEmpty()){
						srs=streamerSource.poll();
						srs.start();
						activeStreamers.add(srs);
					}
				}
			}
		}
		return list;
	}
//	public abstract ListNum<SamLine> nextLines();
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	void spawnThreads(){
		final int maxActive=Tools.max(2, Tools.min((Shared.threads()+4)/5, ffin.length, MAX_FILES));
		streamerSource=new ArrayDeque<SamReadStreamer>(ffin.length);
		activeStreamers=new ArrayDeque<SamReadStreamer>(maxActive);
		for(int i=0; i<ffin.length; i++){
			SamReadStreamer srs=new SamReadStreamer(ffin[i], threads, saveHeader, maxReads);
			streamerSource.add(srs);
		}
		while(activeStreamers.size()<maxActive && !streamerSource.isEmpty()){
			SamReadStreamer srs=streamerSource.poll();
			srs.start();
			activeStreamers.add(srs);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	protected String fname;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	protected long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final boolean saveHeader;

	/** Primary input file */
	final FileFormat[] ffin;
	
	/** Readers */
//	final SamReadStreamer[] streamers;
	private ArrayDeque<SamReadStreamer> streamerSource;
	private ArrayDeque<SamReadStreamer> activeStreamers;
	
	final int threads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int DEFAULT_THREADS=6;
	public static int MAX_FILES=8;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
