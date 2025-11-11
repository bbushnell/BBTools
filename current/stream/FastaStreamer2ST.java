package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile1F;
import fileIO.FileFormat;
import shared.Vector;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/**
 * Single-threaded FASTA file loader using record-based ByteFile1F.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 2026
 */
public class FastaStreamer2ST implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamer2ST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaStreamer2ST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		if(verbose){outstream.println("Made FastaStreamer2ST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamer2ST.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Open the file
		bf=new ByteFile1F(ffin);

		if(verbose){outstream.println("FastaStreamer2ST started.");}
	}

	@Override
	public void close(){
		if(bf!=null){
			bf.close();
			bf=null;
		}
	}

	@Override
	public String fname(){return fname;}

	@Override
	public boolean hasMore(){
		return !finished;
	}

	@Override
	public boolean errorState(){return errorState;}

	@Override
	public boolean paired(){return interleaved;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public long readsProcessed(){return readsProcessed;}

	@Override
	public long basesProcessed(){return basesProcessed;}

	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : new java.util.Random(seed));
	}

	@Override
	public ListNum<Read> nextList(){
		if(finished){return null;}

		ListNum<Read> list=interleaved ? nextListInterleaved() : nextListSingle();

		if(list==null || list.size()==0){
			finished=true;
			close();
			return null;
		}

		listNum++;
		return list;
	}

	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ListNum<Read> nextListSingle(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>(TARGET_LIST_SIZE);
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed;

		long readID=readsProcessed;
		int readsInList=0;
		int bytesInList=0;

		for(byte[] record=bf.nextLine(); record!=null && readsProcessed<maxReads; record=bf.nextLine()){
			if(samplerate>=1f || randy.nextFloat()<samplerate){
				Read r=Vector.fastaRecordToRead(record, readID, pairnum, bb, newlines);
				readList.add(r);
				readsProcessed++;
				basesProcessed+=r.length();
				readsInList++;
				bytesInList+=r.length();
			}
			readID++;

			// Check if we should ship current list
			if(readsInList>=TARGET_LIST_SIZE || bytesInList>=TARGET_LIST_BYTES){
				break;
			}
		}
		if(bb.array.length>8000000 && bytesInList<2000000){bb.shrinkTo(512000);} 
		return reads;
	}

	private ListNum<Read> nextListInterleaved(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>(TARGET_LIST_SIZE);
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed/2;

		long readID=readsProcessed/2;
		int readsInList=0;
		int bytesInList=0;

		Read pending=null;

		for(byte[] record=bf.nextLine(); record!=null && readsProcessed<maxReads; record=bf.nextLine()){
			Read r=Vector.fastaRecordToRead(record, readID, 0, bb, newlines);
			readsProcessed++;
			basesProcessed+=r.length();

			if(pending==null){
				pending=r;
				pending.setPairnum(0);
			}else{
				r.setPairnum(1);
				pending.mate=r;
				r.mate=pending;
				readID++;

				if(samplerate>=1f || randy.nextFloat()<samplerate){
					readList.add(pending);
					readsInList+=2;
					bytesInList+=pending.length()+r.length();
				}
				pending=null;

				// Check if we should ship current list
				if(readsInList>=TARGET_LIST_SIZE || bytesInList>=TARGET_LIST_BYTES){
					break;
				}
			}
		}

		// Pending should be null at batch boundaries, otherwise file has odd number of reads
		if(pending!=null && bf.nextLine()==null){
			throw new RuntimeException("Odd number of reads in interleaved FASTA file: "+fname);
		}
		if(bb.array.length>8000000 && bytesInList<2000000){bb.shrinkTo(512000);} 
		return reads;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;

	/** ByteFile for reading */
	private ByteFile1F bf;

	final int pairnum;
	final boolean interleaved;
	private final ByteBuilder bb=new ByteBuilder(4096);
	private final IntList newlines=new IntList(256);

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads */
	final long maxReads;

	/** Current list number */
	private long listNum=0;

	/** True when file is exhausted */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=200;
	public static int TARGET_LIST_BYTES=262144;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	private float samplerate=1f;
	private java.util.Random randy=null;

}