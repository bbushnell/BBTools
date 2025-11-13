package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile1Fc;
import fileIO.FileFormat;
import shared.KillSwitch;
import structures.IntList;
import structures.ListNum;

/**
 * Single-threaded FASTA file loader using ByteFile1Fc (pre-stripped records).
 * Not suitable for interleaved files.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 12, 2025
 */
public class FastaStreamer2ZT implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamer2ZT(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaStreamer2ZT(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(!interleaved) : "FastaStreamer2ZT does not support interleaved files";
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		if(verbose){outstream.println("Made FastaStreamer2ZT");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public synchronized void start(){
		if(verbose){outstream.println("FastaStreamer2ZT.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Open the file
		bf=new ByteFile1Fc(ffin);

		if(verbose){outstream.println("FastaStreamer2ZT started.");}
	}

	@Override
	public synchronized void close(){
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
	public boolean paired(){return false;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public synchronized long readsProcessed(){return readsProcessed;}

	@Override
	public synchronized long basesProcessed(){return basesProcessed;}

	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : new java.util.Random(seed));
	}

	@Override
	public synchronized ListNum<Read> nextList(){
		if(finished){return null;}

		ListNum<Read> list=nextListSingle();

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

		// Get block of records with newline positions
		byte[] block=bf.nextLine(newlines);
		if(block==null || block.length==0){return null;}

		ArrayList<Read> readList=new ArrayList<Read>();
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed;

		for(int i=0, nl0=-1; i<newlines.size() && readsProcessed<maxReads; i++){
			int nl1=newlines.get(i);
			int nl2=(newlines.size()>i+1 ? newlines.get(i+1) : nl1);
			assert(block[nl0+1]=='>') : nl0+", "+(char)block[nl0+1];
			final byte[] header=KillSwitch.copyOfRange(block, nl0+2, nl1);
			final byte[] bases=(nl2>nl1 ? KillSwitch.copyOfRange(block, nl1+1, nl2) : null);
			if(samplerate>=1f || randy.nextFloat()<samplerate){
				Read r=new Read(bases, null, new String(header, StandardCharsets.US_ASCII), readsProcessed);
				r.setPairnum(pairnum);
				readList.add(r);
				readsProcessed++;
				basesProcessed+=r.length();
			}
			if(bases!=null) {
				i++;
				nl0=nl2;
			}else {nl0=nl1;}
		}

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
	private ByteFile1Fc bf;

	final int pairnum;
	final boolean interleaved;
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