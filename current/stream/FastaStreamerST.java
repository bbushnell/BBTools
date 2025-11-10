package stream;

import java.io.PrintStream;

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
public class FastaStreamerST implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastaStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}
	
	/** Constructor. */
	public FastaStreamerST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		if(verbose){outstream.println("Made FastaStreamerST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamerST.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Open the file
		bf=ByteFile.makeByteFile(ffin);
		
		if(verbose){outstream.println("FastaStreamerST started.");}
	}
	
	@Override
	public void close(){
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
	public long readsProcessed() {return readsProcessed;}
	
	@Override
	public long basesProcessed() {return basesProcessed;}
	
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
		
		ListNum<Read> reads=new ListNum<Read>(null, listNum);
		reads.firstRecordNum=readsProcessed;
		
		final ByteBuilder bb=new ByteBuilder(4096);
		byte[] header=null;
		long readID=readsProcessed;
		
		int readsInList=0;
		int bytesInList=0;
		
		for(byte[] line=bf.nextLine(); line!=null && readsProcessed<maxReads; line=bf.nextLine()){
			
			if(line.length>0 && line[0]=='>'){
				// Save previous record if exists
				if(header!=null){
					if(samplerate>=1f || randy.nextFloat()<samplerate){
						Read r=new Read(bb.toBytes(), null, new String(header), readID, true);
						r.setPairnum(pairnum);
						if(!r.validated()){r.validate(true);}
						reads.add(r);
						readsProcessed++;
						basesProcessed+=r.length();
						readsInList++;
						bytesInList+=r.length();
					}
					readID++;
				}
				
				header=line;
				bb.clear();
				
				// Check if we should ship current list
				if(readsInList>=TARGET_LIST_SIZE || bytesInList>=TARGET_LIST_BYTES){
					break;
				}
			}else{
				bb.append(line);
			}
		}
		
		// Save final record in this batch
		if(header!=null){
			if(samplerate>=1f || randy.nextFloat()<samplerate){
				Read r=new Read(bb.toBytes(), null, new String(header), readID, true);
				r.setPairnum(pairnum);
				if(!r.validated()){r.validate(true);}
				reads.add(r);
				readsProcessed++;
				basesProcessed+=r.length();
			}
		}
		
		return reads;
	}
	
	private ListNum<Read> nextListInterleaved(){
		if(bf==null){return null;}
		
		ListNum<Read> reads=new ListNum<Read>(null, listNum);
		reads.firstRecordNum=readsProcessed/2;
		
		final ByteBuilder bb=new ByteBuilder(4096);
		byte[] header=null;
		long readID=readsProcessed/2;
		
		int readsInList=0;
		int bytesInList=0;
		
		Read pending=null;
		
		for(byte[] line=bf.nextLine(); line!=null && readsProcessed<maxReads; line=bf.nextLine()){
			
			if(line.length>0 && line[0]=='>'){
				// Save previous record if exists
				if(header!=null){
					Read r=new Read(bb.toBytes(), null, new String(header), 0, true);
					if(!r.validated()){r.validate(true);}
					readsProcessed++;
					basesProcessed+=r.length();
					
					if(pending==null){
						pending=r;
						pending.setPairnum(0);
					}else{
						r.setPairnum(1);
						pending.mate=r;
						r.mate=pending;
						pending.numericID=readID;
						r.numericID=readID++;
						
						if(samplerate>=1f || randy.nextFloat()<samplerate){
							reads.add(pending);
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
				
				header=line;
				bb.clear();
			}else{
				bb.append(line);
			}
		}
		
		// Save final record in this batch
		if(header!=null){
			Read r=new Read(bb.toBytes(), null, new String(header), 0, true);
			if(!r.validated()){r.validate(true);}
			readsProcessed++;
			basesProcessed+=r.length();
			
			if(pending==null){
				pending=r;
				pending.setPairnum(0);
			}else{
				r.setPairnum(1);
				pending.mate=r;
				r.mate=pending;
				pending.numericID=readID;
				r.numericID=readID++;
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					reads.add(pending);
				}
				pending=null;
			}
		}
		
		// Pending should be null at batch boundaries, otherwise file has odd number of reads
		if(pending!=null && bf.nextLine()==null){
			throw new RuntimeException("Odd number of reads in interleaved FASTA file: "+fname);
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