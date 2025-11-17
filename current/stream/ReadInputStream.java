package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Abstract base class for reading biological sequence data from various input sources.
 * Provides unified interface for reading sequences with support for both single-read
 * and batch processing. Implementations handle different file formats and sources.
 * @author Brian Bushnell
 */
public abstract class ReadInputStream {
	
	/**
	 * Reads all sequences from a file into an ArrayList.
	 * Auto-detects file format and creates appropriate input stream.
	 *
	 * @param fname Input filename, or null to return null
	 * @param defaultFormat Default format code if detection fails
	 * @param maxReads Maximum number of reads to load
	 * @return ArrayList of Read objects, or null if fname is null
	 */
	public static final ArrayList<Read> toReads(String fname, int defaultFormat, long maxReads){
		if(fname==null){return null;}
		FileFormat ff=FileFormat.testInput(fname, defaultFormat, null, false, true);
		return toReads(ff, maxReads);
	}
	
	/**
	 * Reads all sequences from a FileFormat into an array.
	 * @param ff FileFormat object specifying input source and format
	 * @param maxReads Maximum number of reads to load
	 * @return Array of Read objects, or null if no reads loaded
	 */
	public static final Read[] toReadArray(FileFormat ff, long maxReads){
		ArrayList<Read> list=toReads(ff, maxReads);
		return list==null ? null : list.toArray(new Read[0]);
	}
	
	/**
	 * Reads all sequences from a FileFormat into an ArrayList.
	 * Uses ConcurrentReadInputStream for efficient batch processing.
	 *
	 * @param ff FileFormat object specifying input source and format
	 * @param maxReads Maximum number of reads to load
	 * @return ArrayList containing all loaded reads
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
	
	/** Checks if more sequences are available for reading.
	 * @return true if more reads can be obtained, false otherwise */
	public abstract boolean hasMore();

	/** Resets the input stream to the beginning.
	 * Allows re-reading the same input source from the start. */
	public abstract void restart();
	
	/** Returns true if there was an error, false otherwise */
	public abstract boolean close();

	/** Indicates whether this stream contains paired-end reads.
	 * @return true if reads are paired-end, false for single-end */
	public abstract boolean paired();

	/**
	 * Converts a Read array to an ArrayList.
	 * @param array Array of Read objects to convert
	 * @return ArrayList containing the same reads, or null for empty/null input
	 */
	protected static final ArrayList<Read> toList(Read[] array){
		if(array==null || array.length==0){return null;}
		ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){return errorState;}
	/** TODO */
	protected boolean errorState=false;
	
	/** Returns the filename or source identifier for this input stream.
	 * @return Source filename/identifier, or null if not applicable */
	public abstract String fname();
	
}
