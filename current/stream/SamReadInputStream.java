package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Multithreaded SAM/BAM input stream using SamStreamer.
 * 
 * Provides ReadInputStream interface for SAM and BAM files with automatic format detection
 * and multithreaded parsing. Delegates to SamStreamer for efficient parallel processing
 * while maintaining the familiar ReadInputStream API.
 * 
 * Supports both single-read and interleaved paired-read modes.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date Original, refactored October 23, 2025
 */
public class SamReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		SamReadInputStream sris=new SamReadInputStream(args[0], false, false, true, 
			SamStreamer.DEFAULT_THREADS);
		
		Read r=sris.next();
		System.out.println(r.toText(false));
		System.out.println();
		System.out.println(r.obj.toString());
		System.out.println();
	}
	
	/** Constructor with default thread count. */
	public SamReadInputStream(String fname, boolean loadHeader_, boolean interleaved_, 
			boolean allowSubprocess_, long maxReads_){
		this(fname, loadHeader_, interleaved_, allowSubprocess_, SamStreamer.DEFAULT_THREADS, maxReads_);
	}
	
	/** Constructor with explicit thread count. */
	public SamReadInputStream(String fname, boolean loadHeader_, boolean interleaved_, 
			boolean allowSubprocess_, int threads_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.SAM, null, allowSubprocess_, false), 
			loadHeader_, interleaved_, threads_, maxReads_);
	}
	
	/** Main constructor - creates and starts SamStreamer. */
	public SamReadInputStream(FileFormat ff, boolean loadHeader_, boolean interleaved_, 
			int threads_, long maxReads_){
		loadHeader=loadHeader_;
		interleaved=interleaved_;
		stdin=ff.stdio();
		
		if(!ff.samOrBam()){
			System.err.println("Warning: Did not find expected sam file extension for filename "+
				ff.name());
		}
		
		//Create streamer with appropriate thread count
		int threads=(threads_<=0 ? SamStreamer.DEFAULT_THREADS : threads_);
		streamer=SamStreamer.makeStreamer(ff, threads, loadHeader_, true, maxReads_, true);
		
		//Extract header if requested
		if(loadHeader){
			header=streamer.header;
			if(header!=null){setSharedHeader(header);}
		}
		streamer.start();
	}

	@Override
	public void start(){
//		streamer.start(); //Start streamer threads
	}
	
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			fillBuffer();
		}
		return (buffer!=null && next<buffer.size());
	}

	@Override
	public Read next(){
		if(!hasMore()){return null;}
		Read r=buffer.set(next, null);
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){
			throw new RuntimeException("'next' should not be used when doing blockwise access.");
		}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}
	
	/** Fill buffer from streamer. */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		//Get next list from streamer
		ListNum<Read> ln=streamer.nextReads();
		if(ln==null || ln.list==null){
			buffer=new ArrayList<Read>(); //Empty buffer signals EOF
			return;
		}
		
		buffer=ln.list;
		
		//Assign numeric IDs if not already set
		for(Read r : buffer){
			if(r.numericID<0){
				r.numericID=nextReadID++;
			}
		}
		
		generated+=buffer.size();
	}

	@Override
	public boolean close(){
		//Streamer cleanup handled automatically
		return errorState;
	}
	
	@Override
	public synchronized void restart(){
		throw new RuntimeException("SamReadInputStream does not support restart.");
	}
	
	/** Get shared header, optionally waiting for it to be read. */
	public static synchronized ArrayList<byte[]> getSharedHeader(boolean wait){
		if(!wait || SHARED_HEADER!=null){return SHARED_HEADER;}
		System.err.println("Waiting on header to be read from a sam file.");
		while(SHARED_HEADER==null){
			try{
				SamReadInputStream.class.wait(100);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		return SHARED_HEADER;
	}
	
	/** Set shared header for all SamReadInputStream instances. */
	public static synchronized void setSharedHeader(ArrayList<byte[]> list){
		SHARED_HEADER=list;
		SamReadInputStream.class.notifyAll();
	}
	
	/** Trim whitespace and annotations from SQ header reference names. */
	public static byte[] trimHeaderSQ(byte[] line){
		if(line==null || !Tools.startsWith(line, "@SQ")){return line;}
		
		final int idx=Tools.indexOfDelimited(line, "SN:", 2, (byte)'\t');
		if(idx<0){
			assert(false) : "Bad header: "+new String(line);
			return line;
		}
		
		int trimStart=-1;
		for(int i=idx; i<line.length; i++){
			final byte b=line[i];
			if(b=='\t'){return line;}
			if(Character.isWhitespace(b)){
				trimStart=i;
				break;
			}
		}
		if(trimStart<0){return line;}
		
		final int trimStop=Tools.indexOf(line, (byte)'\t', trimStart+1);
		final int bbLen=trimStart+(trimStop<0 ? 0 : line.length-trimStop);
		final ByteBuilder bb=new ByteBuilder(bbLen);
		for(int i=0; i<trimStart; i++){bb.append(line[i]);}
		if(trimStop>=0){
			for(int i=trimStop; i<line.length; i++){bb.append(line[i]);}
		}
		assert(bb.length==bbLen) : bbLen+", "+bb.length+", idx="+idx+", trimStart="+
			trimStart+", trimStop="+trimStop+"\n\n"+new String(line)+"\n\n"+bb+"\n\n";
		
		return bb.array;
	}
	
	@Override
	public String fname(){return streamer.fname;}
	
	@Override
	public boolean paired(){return interleaved;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared header across all SamReadInputStream instances */
	private static volatile ArrayList<byte[]> SHARED_HEADER;

	/** Current buffer of reads */
	private ArrayList<Read> buffer=null;
	/** Header lines from SAM/BAM file */
	private ArrayList<byte[]> header=null;
	/** Position in current buffer */
	private int next=0;
	
	/** Underlying multithreaded streamer */
	private final SamStreamer streamer;
	/** True if reads are interleaved paired-end */
	private final boolean interleaved;
	/** True if header should be loaded and shared */
	private final boolean loadHeader;
	
	/** Total reads generated */
	public long generated=0;
	/** Total reads consumed */
	public long consumed=0;
	/** Next numeric ID to assign */
	private long nextReadID=0;
	
	/** True if reading from stdin */
	public final boolean stdin;

}