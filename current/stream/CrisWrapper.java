package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Wraps a cris to allow single-read next() capability, and the ability to go back.
 * @author Brian Bushnell
 * @date Jul 18, 2014
 *
 */
public class CrisWrapper {
	
	/**
	 * Creates a wrapper for reading from two FileFormats without quality files.
	 *
	 * @param maxReads Maximum number of reads to process
	 * @param keepSamHeader Whether to preserve SAM headers
	 * @param ff1 Primary input file format
	 * @param ff2 Secondary input file format (may be null)
	 */
	public CrisWrapper(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2){
		this(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null);
	}
	
	/**
	 * Creates a wrapper for reading from two FileFormats with optional quality files.
	 *
	 * @param maxReads Maximum number of reads to process
	 * @param keepSamHeader Whether to preserve SAM headers
	 * @param ff1 Primary input file format
	 * @param ff2 Secondary input file format (may be null)
	 * @param qf1 Primary quality file path (may be null)
	 * @param qf2 Secondary quality file path (may be null)
	 */
	public CrisWrapper(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		this(ConcurrentReadInputStream.getReadInputStream(maxReads, ff1.samOrBam(), ff1, ff2, qf1, qf2), true);
	}
		

	/**
	 * Creates a wrapper around an existing ConcurrentReadInputStream.
	 * @param cris_ The input stream to wrap
	 * @param start Whether to start the stream immediately
	 */
	public CrisWrapper(ConcurrentReadInputStream cris_, boolean start){
		initialize(cris_, start);
	}
	
	/**
	 * Initializes the wrapper with a ConcurrentReadInputStream and loads the first read list.
	 * Sets up internal state and handles empty streams by closing them immediately.
	 * @param cris_ The input stream to wrap
	 * @param start Whether to start the stream
	 */
	public void initialize(ConcurrentReadInputStream cris_, boolean start){
		cris=cris_;
		if(start){cris.start();}
		ln=cris.nextList();
		reads=(ln==null ? null : ln.list);
		if(reads==null || reads.size()==0){
			reads=null;
			//System.err.println("Empty.");
			cris.returnList(ln.id, true);
			errorState|=ReadWrite.closeStream(cris);
		}
		index=0;
		//System.err.println("Initialized.");
	}
	
	/**
	 * Returns the next read from the input stream.
	 * Automatically manages read list buffering and stream lifecycle.
	 * Returns null when no more reads are available.
	 * @return Next read or null if stream is exhausted
	 */
	public Read next(){
		//System.err.println("*******1");
		Read r=null;
		if(reads==null || index>=reads.size()){
			//System.err.println("*******2");
			if(reads==null){return null;}
			index=0;
			if(reads.size()==0){
				reads=null;
				cris.returnList(ln.id, true);
				errorState|=ReadWrite.closeStream(cris);
				return null;
			}
			cris.returnList(ln.id, false);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			if(reads==null){
				//System.err.println("*******3");
				cris.returnList(ln.id, true);
				errorState|=ReadWrite.closeStream(cris);
				//System.err.println("Returning null (2)");
				return null;
			}
		}
		//System.err.println("*******4");
		if(index<reads.size()){
			//System.err.println("*******5");
			r=reads.get(index);
			index++;
		}else{
			//System.err.println("*******6");
			//System.err.println("Recalling");
			return next();
		}
		//System.err.println("*******7");
		//System.err.println("Returning "+(r==null ? "null" : r.id));
		return r;
	}
	
	/** Moves the read position back by one to re-read the previous read.
	 * Requires that at least one read has been consumed from the current list. */
	public void goBack(){
		assert(index>0);
		index--;
	}
	
	/** Current read list container from the input stream */
	private ListNum<Read> ln;
	/** Current list of reads being iterated through */
	private ArrayList<Read> reads;
	/** Current position within the reads list */
	private int index;
	/** The underlying concurrent read input stream */
	public ConcurrentReadInputStream cris;
	/** Indicates whether an error occurred during stream operations */
	public boolean errorState=false;
	
}
