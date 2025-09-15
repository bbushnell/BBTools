package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import structures.ListNum;

/**
 * Loads sam files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public abstract class SamStreamer {
	
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
		SamStreamer x=new SamReadStreamer(args[0], threads, false, -1);
		
		//Run the object
		x.start();
		x.test();
		
		t.stop("Time: ");
	}
	
	/**
	 * Constructor.
	 */
	public SamStreamer(String fname_, int threads_, boolean saveHeader_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_, saveHeader_, maxReads_);
	}
	
	/**
	 * Constructor.
	 */
	public SamStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, long maxReads_){
		fname=ffin_.name();
		threads=threads_;
		ffin=ffin_;
		saveHeader=saveHeader_;
		header=(saveHeader ? new ArrayList<byte[]>() : null);
		maxReads=(maxReads_<1 ? Long.MAX_VALUE : maxReads_);
		
		inq=new ArrayBlockingQueue<ListNum<byte[]>>(threads/2+2);
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
	public abstract ListNum<Read> nextReads();
	public abstract ListNum<SamLine> nextLines();
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Thread 0 reads the actual file, optionally processes the header,
	 * and produces lists of byte[] for other threads to format as SamLines.
	 * @param tid Thread number (should be 0).
	 */
	public final void processBytes0(int tid){
		if(verbose){outstream.println("tid "+tid+" started processBytes.");}

//		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=true;
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		long listNumber=0;
		long reads=0;
		
		final int limit=LIST_SIZE;
		ArrayList<byte[]> list=new ArrayList<byte[]>(limit);
		for(byte[] line=bf.nextLine(); line!=null && reads<maxReads; line=bf.nextLine()){
			assert(line!=null);
//			outstream.println("a");
			if(header!=null && line[0]=='@'){
				if(Shared.TRIM_RNAME){line=SamReadInputStream.trimHeaderSQ(line);}
				header.add(line);
			}else{
				if(header!=null){
					SamReadInputStream.setSharedHeader(header);
					header=null;
				}
				reads++;
				list.add(line);
				if(list.size()>=limit){
					//					outstream.println("b");
					//					outstream.println(inq.size()+", "+inq.remainingCapacity());
					putBytes(new ListNum<byte[]>(list, listNumber));
					listNumber++;
					//					outstream.println("c");
					list=new ArrayList<byte[]>(limit);
				}
			}
//			outstream.println("d");
		}
		if(header!=null){
			SamReadInputStream.setSharedHeader(header);
			header=null;
		}
		if(verbose){outstream.println("tid "+tid+" ran out of input.");}
		if(list.size()>0){
			putBytes(new ListNum<byte[]>(list, listNumber));
			listNumber++;
			list=null;
		}
		if(verbose || verbose2){outstream.println("tid "+tid+" done reading bytes.");}
		putBytes(POISON_BYTES);
		if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning.");}
		bf.close();
		if(verbose || verbose2){outstream.println("tid "+tid+" closed stream.");}
	}
	
	final void putBytes(ListNum<byte[]> list){
//		if(verbose){outstream.println("tid "+tid+" putting blist size "+list.size());}
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		if(verbose){outstream.println("tid "+tid+" done putting blist");}
	}
	
	final ListNum<byte[]> takeBytes(){
//		if(verbose){outstream.println("tid "+tid+" taking blist");}
		ListNum<byte[]> list=null;
		while(list==null){
			try {
				list=inq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		if(verbose){outstream.println("tid "+tid+" took blist size "+list.size());}
		return list;
	}
	
	/** Spawn process threads */
	abstract void spawnThreads();
	
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

//	public void setMaxReads(long x){
//		maxReads=(x<1 ? Long.MAX_VALUE : x);
//	}
	
	/** Quit after processing this many input reads */
	private final long maxReads;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final boolean saveHeader;
	
	/** Primary input file */
	final FileFormat ffin;
	
	final ArrayBlockingQueue<ListNum<byte[]>> inq;
	
	final int threads;
	
	public ArrayList<byte[]> header;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	static final ListNum<Read> POISON_READS=new ListNum<Read>(null, -1);
	static final ListNum<SamLine> POISON_LINES=new ListNum<SamLine>(null, -1);
	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, -1);
	public static int LIST_SIZE=200;
	public static int DEFAULT_THREADS=6;
	
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
