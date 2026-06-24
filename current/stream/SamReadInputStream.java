package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Multithreaded SAM/BAM reader built on Streamer.
 * Provides the ReadInputStream API for alignment files with automatic format detection,
 * optional header sharing, and streaming-friendly block iteration.
 * Supports single-read and interleaved paired-read inputs via the underlying streamer.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date Original, refactored October 23, 2025
 */
public class SamReadInputStream extends ReadInputStream {
	
	/** Demonstration entry point that iterates through a SAM/BAM file and reports throughput.
	 * @param args Command-line arguments; expects the input filename in args[0] */
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
	
	/**
	 * Creates a SamReadInputStream for the given filename with default thread count.
	 * Delegates to the threaded constructor and optionally loads the SAM/BAM header.
	 *
	 * @param fname Input SAM/BAM filename
	 * @param loadHeader_ Whether to parse and share the header
	 * @param allowSubprocess_ Allow use of subprocess for compressed input
	 * @param maxReads_ Maximum reads to stream (-1 for all)
	 */
	public SamReadInputStream(String fname, boolean loadHeader_, 
			boolean allowSubprocess_, long maxReads_){
		this(fname, loadHeader_, allowSubprocess_, -1, maxReads_);
	}
	
	/**
	 * Creates a SamReadInputStream with an explicit thread count.
	 * Initializes format detection, header loading, and streamer construction.
	 *
	 * @param fname Input SAM/BAM filename
	 * @param loadHeader_ Whether to parse and share the header
	 * @param allowSubprocess_ Allow use of subprocess for compressed input
	 * @param threads_ Number of threads for the Streamer (-1 for automatic)
	 * @param maxReads_ Maximum reads to stream (-1 for all)
	 */
	public SamReadInputStream(String fname, boolean loadHeader_, 
			boolean allowSubprocess_, int threads_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.SAM, null, allowSubprocess_, false), 
			loadHeader_, threads_, maxReads_);
	}
	
	/**
	 * Creates a SamReadInputStream from a FileFormat description.
	 * Sets stdin flag, warns on unexpected extensions, and starts a multithreaded streamer.
	 *
	 * @param ff FileFormat describing the input source
	 * @param loadHeader_ Whether to parse and share the header
	 * @param threads_ Number of threads for the Streamer (-1 for automatic)
	 * @param maxReads_ Maximum reads to stream (-1 for all)
	 */
	public SamReadInputStream(FileFormat ff, boolean loadHeader_,
			int threads_, long maxReads_){
		SAM_INPUT_PRESENT=true; //A SAM/BAM input now exists, so SAM output may legitimately wait for its shared header.
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
	
	/** Returns whether additional reads are available from the streamer. */
	@Override
	public boolean hasMore(){
		return streamer.hasMore();
	}
	
	/**
	 * Retrieves the next block of reads from the streamer.
	 * Returns null when the stream is exhausted.
	 * @return List of Read objects for the next chunk, or null if no data remains
	 */
	@Override
	public ArrayList<Read> nextList(){
		ListNum<Read> ln=streamer.nextList();
		return ln==null || ln.isEmpty() ? null : ln.list;
	}

	/** Closes the underlying streamer and returns the error state.
	 * @return true if any errors were detected, false otherwise */
	@Override
	public boolean close(){
		streamer.close();
		errorState|=streamer.errorState();//#001 fix: was returning the never-set LOCAL errorState, dropping the streamer's. SamStreamer/BamStreamer set errorState=true when a reader thread fails (e.g. a truncated/corrupt BAM), so without this a failed SAM/BAM read was reported as SUCCESS.
		return errorState;
	}

	/** @return true if this stream OR its underlying multithreaded Streamer detected an error.
	 * #001 fix: SamReadInputStream did not override errorState(), so it returned the always-false base field and masked streamer errors up the cris error chain (cris.errorState() -> producer.errorState()). */
	@Override
	public boolean errorState(){return errorState || streamer.errorState();}
	
	/** Unsupported operation for SamReadInputStream.
	 * Always throws RuntimeException because restarting is not implemented. */
	@Override
	public synchronized void restart(){
		throw new RuntimeException("SamReadInputStream does not support restart.");
	}
	
	/**
	 * Returns the globally shared SAM/BAM header, optionally waiting until it becomes available.
	 * @param wait Whether to block until the header is populated
	 * @return Shared header lines as byte arrays, or null if unavailable and wait is false
	 */
	public static synchronized ArrayList<byte[]> getSharedHeader(boolean wait){
		if(!wait || SHARED_HEADER!=null){return SHARED_HEADER;}
		//Crash-loud, never hang [SamReadInputStream#002]: if no SAM/BAM input stream was ever opened, no
		//reader will ever call setSharedHeader, so waiting here deadlocks (e.g. fastq/fasta -> sam, the
		//common case). Return null (unavailable) and let the caller generate a header instead. This also
		//resolves the original //TODO for the no-shared-header case.
		if(!SAM_INPUT_PRESENT){return SHARED_HEADER;}
		if(printHeaderWait) {System.err.println("Waiting on header to be read from a sam file.");}
		while(SHARED_HEADER==null){
			try{
				SamReadInputStream.class.wait(100);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		return SHARED_HEADER;
	}
	
	/** Sets the globally shared header for all SamReadInputStream instances and notifies waiters.
	 * @param list Header lines to share across instances */
	public static synchronized void setSharedHeader(ArrayList<byte[]> list){
		SHARED_HEADER=list;
		SamReadInputStream.class.notifyAll();
	}

	/** Marks that a SAM/BAM input source now exists in this JVM, so getSharedHeader(true) may legitimately
	 * block until that input parses and sets the shared header. Historically this was set in the
	 * SamReadInputStream constructor (L80), but the native BamStreamer/SamStreamer streamers bypass this
	 * class entirely, leaving the flag false; getSharedHeader(true) then hit the #002 no-hang gate and
	 * returned null instead of waiting, racing callers like var2/ScafMap.loadSamHeader (assert header!=null).
	 * StreamerFactory calls this synchronously, on the constructing thread, whenever it builds a native
	 * SAM/BAM streamer with saveHeader=true (which commits that streamer to calling setSharedHeader), so the
	 * gate invariant is visible before the streamer's worker thread starts. */
	public static void markSamInputPresent(){SAM_INPUT_PRESENT=true;}
	
	/**
	 * Normalizes an @SQ header line by stripping whitespace from the reference name.
	 * Leaves non-@SQ lines unchanged.
	 * @param line Header line to normalize
	 * @return New byte array with trimmed reference name, or the original line if unchanged
	 */
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
	
	/** Returns the input filename reported by the underlying streamer.
	 * @return Input filename */
	@Override
	public String fname(){return streamer.fname();}
	
	/**
	 * Indicates whether this stream reports paired reads.
	 * Always false; pairing is handled by the streamer rather than this interface.
	 * @return false
	 */
	@Override
	public boolean paired(){return false;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Globally shared header lines available to all SamReadInputStream instances.
	 */
	private static volatile ArrayList<byte[]> SHARED_HEADER;
	/** True once any SAM/BAM input stream has been opened in this JVM. Gates getSharedHeader's blocking
	 * wait so a non-SAM input (fastq/fasta), which never sets the shared header, does not deadlock a
	 * SAM-output writer that requested it. See stream/SamReadInputStream#002. */
	static volatile boolean SAM_INPUT_PRESENT=false;
	/** Controls diagnostic logging while waiting for a shared header. */
	public static boolean printHeaderWait=false;
	
	/** Header lines read from the current SAM/BAM stream, if loaded. */
	private ArrayList<byte[]> header=null;
	
	/** Multithreaded Streamer providing buffered SAM/BAM reads. */
	private final Streamer streamer;
	/** Whether this stream should parse and retain the SAM/BAM header. */
	private final boolean loadHeader;
	
	/** True if input is sourced from standard input rather than a file. */
	public final boolean stdin;

}