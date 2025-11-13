package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Multithreaded SAM/BAM input stream using Streamer.
 * 
 * Provides ReadInputStream interface for SAM and BAM files with automatic format detection
 * and multithreaded parsing. Delegates to Streamer for efficient parallel processing
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
		SamReadInputStream sris=new SamReadInputStream(args[0], false, true, -1, -1);
		
		Timer t=new Timer();
		long reads=0, bases=0;
		for(ArrayList<Read> ln=sris.nextList(); ln!=null; ln=sris.nextList()) {
			for(Read r : ln) {bases+=r.pairLength();}
			reads+=ln.size();
		}
		t.stop();
		System.err.println();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/** Constructor with default thread count. */
	public SamReadInputStream(String fname, boolean loadHeader_, 
			boolean allowSubprocess_, long maxReads_){
		this(fname, loadHeader_, allowSubprocess_, -1, maxReads_);
	}
	
	/** Constructor with explicit thread count. */
	public SamReadInputStream(String fname, boolean loadHeader_, 
			boolean allowSubprocess_, int threads_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.SAM, null, allowSubprocess_, false), 
			loadHeader_, threads_, maxReads_);
	}
	
	/** Main constructor - creates and starts Streamer. */
	public SamReadInputStream(FileFormat ff, boolean loadHeader_, 
			int threads_, long maxReads_){
		loadHeader=loadHeader_;
		stdin=ff.stdio();
		
		if(!ff.samOrBam()){
			System.err.println("Warning: Did not find expected sam file extension for filename "+
				ff.name());
		}
		
		//Create streamer with appropriate thread count
		streamer=StreamerFactory.makeSamOrBamStreamer(ff, threads_, loadHeader_, true, maxReads_, true);
		
//		//Extract header if requested
//		if(loadHeader){
//			header=streamer.header;
//			if(header!=null){setSharedHeader(header);}
//		}
		streamer.start();
	}
	
	@Override
	public boolean hasMore(){
		return streamer.hasMore();
	}
	
	@Override
	public ArrayList<Read> nextList(){
		ListNum<Read> ln=streamer.nextList();
		return ln==null || ln.isEmpty() ? null : ln.list;
	}

	@Override
	public boolean close(){
		streamer.close();
		return errorState;
	}
	
	@Override
	public synchronized void restart(){
		throw new RuntimeException("SamReadInputStream does not support restart.");
	}
	
	/** Get shared header, optionally waiting for it to be read. */
	public static synchronized ArrayList<byte[]> getSharedHeader(boolean wait){
		if(!wait || SHARED_HEADER!=null){return SHARED_HEADER;}
		if(printHeaderWait) {System.err.println("Waiting on header to be read from a sam file.");}
		while(SHARED_HEADER==null){//TODO:  Test with headerless sam, should populate with an empty list
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
	public String fname(){return streamer.fname();}
	
	@Override
	public boolean paired(){return false;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared header across all SamReadInputStream instances */
	private static volatile ArrayList<byte[]> SHARED_HEADER;
	public static boolean printHeaderWait=false;
	
	/** Header lines from SAM/BAM file */
	private ArrayList<byte[]> header=null;
	
	/** Underlying multithreaded streamer */
	private final Streamer streamer;
	/** True if header should be loaded and shared */
	private final boolean loadHeader;
	
	/** True if reading from stdin */
	public final boolean stdin;

}