package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Abstract base for reading sequence data from various sources (single/batch).
 * Provides static helpers to load all reads and instance methods for streaming access.
 * @author Brian Bushnell
 */
public abstract class ReadInputStream {
	
	/**
	 * Loads all reads from a filename (auto-detect format) into an ArrayList.
	 * @param fname Input filename (nullable)
	 * @param defaultFormat Default format code if detection fails
	 * @param maxReads Max reads to load
	 * @return List of reads or null if fname is null
	 */
	public static final ArrayList<Read> toReads(String fname, int defaultFormat, long maxReads){
		if(fname==null){return null;}
		FileFormat ff=FileFormat.testInput(fname, defaultFormat, null, false, true);
		return toReads(ff, maxReads);
	}
	
	/**
	 * Loads all reads from a FileFormat into an array.
	 * @param ff FileFormat describing input
	 * @param maxReads Max reads to load
	 * @return Array of reads or null if none
	 */
	public static final Read[] toReadArray(FileFormat ff, long maxReads){
		ArrayList<Read> list=toReads(ff, maxReads);
		return list==null ? null : list.toArray(new Read[0]);
	}
	
	/**
	 * Loads all reads from a FileFormat via ConcurrentReadInputStream.
	 * @param ff FileFormat describing input
	 * @param maxReads Max reads to load
	 * @return List of reads
	 */
	public static final ArrayList<Read> toReads(FileFormat ff, long maxReads){
		ArrayList<Read> list=new ArrayList<Read>();

		/* Start an input stream */
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start(); //4567
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		/* Iterate through read lists from the input stream */
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			list.addAll(reads);
			
			/* Dispose of the old list and fetch a new one */
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		/* Cleanup */
		cris.returnList(ln);
		ReadWrite.closeStream(cris);
		return list;
	}
	
	/**
	 * Reads the next batch of sequences from the input stream.
	 * Batch processing is more efficient than individual reads.
	 * @return ArrayList of Read objects, or null if no more reads
	 */
	public abstract ArrayList<Read> nextList();
	
	/** Returns true if additional reads are available from the stream. */
	public abstract boolean hasMore();

	/** Resets the stream to the beginning to reread input. */
	public abstract void restart();
	
	/** Closes the input stream and releases associated resources.
	 * @return true if there was an error during closing, false otherwise */
	public abstract boolean close();

	/** Indicates whether the stream supplies paired-end reads.
	 * @return true if paired, false otherwise */
	public abstract boolean paired();

	protected static final ArrayList<Read> toList(Read[] array){
		if(array==null || array.length==0){return null;}
		ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}
	
	/** Returns true if this stream has detected an error */
	public boolean errorState(){return errorState;}
	/** Error state flag; true if stream encountered an error */
	protected boolean errorState=false;
	
	/** Returns the source identifier (e.g., filename) for this stream.
	 * @return Source name or null */
	public abstract String fname();
	
}
