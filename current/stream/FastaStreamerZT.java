package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Single-threaded FASTA file loader.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 10, 2025
 */
public class FastaStreamerZT implements Streamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamerZT(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaStreamerZT(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		if(verbose){outstream.println("Made FastaStreamerZT");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public synchronized void start(){
		if(verbose){outstream.println("FastaStreamerZT.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Open the file
		bf=ByteFile.makeByteFile(ffin);

		if(verbose){outstream.println("FastaStreamerZT started.");}
	}

	@Override
	public synchronized void close(){
		if(bf!=null){
			bf.close();
			bf=null;
		}
	}

	@Override
	public String fname() {return fname;}

	@Override
	public boolean hasMore(){
		return !finished;
	}

	@Override
	public boolean errorState() {return errorState;}

	@Override
	public boolean paired(){return interleaved;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public synchronized long readsProcessed() {return readsProcessed;}

	@Override
	public synchronized long basesProcessed() {return basesProcessed;}

	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : new java.util.Random(seed));
	}

	@Override
	public synchronized ListNum<Read> nextList(){
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

		final ByteBuilder bb=new ByteBuilder(4096);

		int readsInList=0;
		int bytesInList=0;
		byte[] line=null;
		for(line=bf.nextLine(); line!=null && readsProcessed<maxReads; line=bf.nextLine()){

			if(line.length>0 && line[0]=='>'){
				if(header!=null) {
					Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readsProcessed);
					r.setPairnum(pairnum);
					readList.add(r);
					readsProcessed++;
					basesProcessed+=r.length();
					readsInList++;
					bytesInList+=r.length();
				}
				header=null;
				bb.clear();
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){header=line;}
				if(readsInList>=TARGET_LIST_SIZE || bytesInList>TARGET_LIST_BYTES) {break;}
			}else if(header!=null){
				bb.append(line);
			}
		}
		if(line==null && header!=null) {//EOF
			Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readsProcessed);
			r.setPairnum(pairnum);
			readList.add(r);
			readsProcessed++;
			basesProcessed+=r.length();
			readsInList++;
			bytesInList+=r.length();
			header=null;
			bb.clear();
		}
		return reads;
	}

	private ListNum<Read> nextListInterleaved(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>(TARGET_LIST_SIZE);
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed/2;

		final ByteBuilder bb=new ByteBuilder(4096);
		long readID=readsProcessed/2;

		int readsInList=0;
		int bytesInList=0;

		Read pending=null;
		byte[] line=null;
		for(line=bf.nextLine(); line!=null && readsProcessed<maxReads; line=bf.nextLine()){

			if(line.length>0 && line[0]=='>'){
				if(header!=null){
					// Finish current read
					Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0);
					readsProcessed++;
					basesProcessed+=r.length();
					bb.clear();

					if(pending==null){
						// This is read1
						pending=r;
						pending.setPairnum(0);
						header=line; // Start building read2
					}else{
						// This is read2
						r.setPairnum(1);
						pending.mate=r;
						r.mate=pending;
						pending.numericID=readID;
						r.numericID=readID++;

						readList.add(pending);
						readsInList+=2;
						bytesInList+=pending.length()+r.length();
						pending=null;

						// Decide whether to start building next pair
						header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;

						// Check if we should ship current list
						if(readsInList>=TARGET_LIST_SIZE || bytesInList>=TARGET_LIST_BYTES){
							break;
						}
					}
				}else{
					// Not currently building - decide whether to start
					header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;
					bb.clear();
				}
			}else if(header!=null){
				// Accumulate bases
				bb.append(line);
			}
		}
		
		if(line==null && header!=null) {//EOF
		    // Finish current read
		    Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0);
		    readsProcessed++;
		    basesProcessed+=r.length();
		    bb.clear();
		    
		    if(pending==null){
		        // File ends on read1 - incomplete pair
		        pending=r;
		        pending.setPairnum(0);
		        readList.add(pending);
		    }else{
		        // This is read2 - complete the pair
		        r.setPairnum(1);
		        pending.mate=r;
		        r.mate=pending;
		        pending.numericID=readID;
		        r.numericID=readID++;

		        readList.add(pending);
		        readsInList+=2;
		        bytesInList+=pending.length()+r.length();
		        pending=null;
		    }
		    header=null;
		}

		// Handle incomplete pair at end
		assert(pending==null) : "Odd number of reads in interleaved FASTA file: "+fname;

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
	private ByteFile bf;

	final int pairnum;
	final boolean interleaved;

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
	
	/** Header line pending for next list */
	private byte[] header=null;

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