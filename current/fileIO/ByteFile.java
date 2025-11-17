package fileIO;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;

import shared.Shared;
import structures.ListNum;


/**
 * Abstract base class for efficient byte-level file reading and processing
 * across different file formats and input strategies.
 * Provides a flexible, extensible framework for reading files as byte arrays
 * with multiple implementation strategies and performance optimizations.
 * Factory methods automatically select optimal implementation based on
 * system resources and file characteristics.
 *
 * @author Brian Bushnell
 */
public abstract class ByteFile {
	
//	public static final ByteFile makeByteFile(String fname){
//		return makeByteFile(fname, false, true);
//	}
	
	/**
	 * Creates a ByteFile1 instance for reading the specified file.
	 * Forces usage of ByteFile1 implementation regardless of system settings.
	 *
	 * @param fname The filename to read
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return ByteFile1 instance for the specified file
	 */
	public static final ByteFile makeByteFile1(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return new ByteFile1(ff);
	}
	
	/**
	 * Creates an appropriate ByteFile instance for reading the specified file.
	 * Automatically selects optimal implementation based on system resources.
	 *
	 * @param fname The filename to read
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return ByteFile instance (ByteFile1 or ByteFile2) for the specified file
	 */
	public static final ByteFile makeByteFile(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return makeByteFile(ff);
	}
	
	/**
	 * Creates an appropriate ByteFile instance for the specified FileFormat.
	 * Uses default type selection logic.
	 * @param ff The FileFormat describing the input file
	 * @return ByteFile instance for the specified format
	 */
	public static final ByteFile makeByteFile(FileFormat ff){
		return makeByteFile(ff, 0);
	}
	
	/**
	 * Creates a ByteFile instance with explicit type selection.
	 * Type 1 forces ByteFile1, type 2 forces ByteFile2.
	 * Type 0 uses automatic selection based on system resources and settings.
	 * Prefers ByteFile2 for multi-threaded systems unless forced otherwise.
	 *
	 * @param ff The FileFormat describing the input file
	 * @param type Implementation type (0=auto, 1=ByteFile1, 2=ByteFile2)
	 * @return ByteFile instance of the specified or selected type
	 */
	public static final ByteFile makeByteFile(FileFormat ff, int type){
		type=pickType(type);
		if(type==4){return new ByteFile4(ff);}
		if(type==3){return new ByteFile3(ff);}
		if(type==2){return new ByteFile2(ff);}
		return new ByteFile1(ff);
	}
	
	/**
	 * Protected constructor for ByteFile instances.
	 * Validates that the FileFormat is configured for reading.
	 * @param ff_ The FileFormat for this ByteFile instance
	 */
	protected ByteFile(FileFormat ff_){
		ff=ff_;
		assert(ff.read()) : ff;
	}
	
	/**
	 * Reads the entire file and returns all lines as byte arrays.
	 * Convenient method for loading entire files into memory.
	 * Uses initial capacity of 4096 lines for efficiency.
	 * @return ArrayList containing all lines as byte arrays
	 */
	public final ArrayList<byte[]> toByteLines(){
		
		byte[] s=null;
		ArrayList<byte[]> list=new ArrayList<byte[]>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list;
	}
	
	/**
	 * Static utility method to read all lines from a FileFormat.
	 * Creates temporary ByteFile instance and automatically closes it.
	 * @param ff The FileFormat to read from
	 * @return ArrayList containing all lines as byte arrays
	 */
	public static final ArrayList<byte[]> toLines(FileFormat ff){
		ByteFile bf=makeByteFile(ff);
		ArrayList<byte[]> lines=bf.toByteLines();
		bf.close();
		return lines;
	}
	
	/**
	 * Static utility method to read all lines from a filename.
	 * Creates FileFormat and ByteFile instances automatically.
	 * @param fname The filename to read
	 * @return ArrayList containing all lines as byte arrays
	 */
	public static final ArrayList<byte[]> toLines(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, false);
		return toLines(ff);
	}
	
	/**
	 * Counts the total number of lines in the file.
	 * Resets the file position after counting for subsequent reading.
	 * Efficient method that doesn't store line content in memory.
	 * @return Total number of lines in the file
	 */
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		
		return count;
	}

	/** Resets the file reader to the beginning of the file.
	 * Implementation varies by subclass. */
	public abstract void reset();
	/** Resets the internal list ID counter to zero.
	 * Called by subclass reset() implementations. */
	final void superReset(){
		nextID=0;
	}
	
	/**
	 * Reads the next batch of lines (up to 200) as a numbered list.
	 * Provides efficient batch processing for large files.
	 * Each ListNum contains up to 200 lines with a unique sequential ID.
	 * @return ListNum containing up to 200 byte array lines, or null if EOF
	 */
	public synchronized ListNum<byte[]> nextList(){
		byte[] line=nextLine();
		if(line==null){return null;}
		final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
		ArrayList<byte[]> list=new ArrayList<byte[]>(slimit);
		list.add(line);
		int bytes=line.length;
		
		for(int i=1; i<slimit && bytes<blimit; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
			bytes+=line.length;
		}
		ListNum<byte[]> ln=new ListNum<byte[]>(list, nextID);
		nextID++;
		return ln;
	}
	
	/**
	 * Checks if the file exists or is a special input source.
	 * Returns true for stdin, stdin variants, JAR resources, or existing files.
	 * @return true if the file exists or is accessible
	 */
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists(); //TODO Ugly and unsafe hack for files in jars
	}

	/**
	 * Returns the underlying InputStream for this ByteFile.
	 * Implementation varies by subclass.
	 * @return InputStream for direct access to file data
	 */
	public abstract InputStream is();
	/**
	 * Returns the current line number being read.
	 * Implementation varies by subclass.
	 * @return Current line number (1-based)
	 */
	public abstract long lineNum();
	
	/** Returns true if there was an error */
	public abstract boolean close();
	
	/**
	 * Reads the next line from the file as a byte array.
	 * Core method for line-by-line file reading.
	 * Implementation varies by subclass.
	 * @return Next line as byte array, or null if EOF
	 */
	public abstract byte[] nextLine();
	
//	public final void pushBack(byte[] line){
//		assert(pushBack==null);
//		pushBack=line;
//	}
	
	/**
	 * Pushes a line back to be returned by the next nextLine() call.
	 * Allows single-line lookahead functionality.
	 * Implementation varies by subclass.
	 * @param line The line to push back for re-reading
	 */
	public abstract void pushBack(byte[] line);
	
	/**
	 * Returns whether the file is currently open for reading.
	 * Implementation varies by subclass.
	 * @return true if the file is open
	 */
	public abstract boolean isOpen();
	
	/** Returns the filename from the associated FileFormat */
	public final String name(){return ff.name();}
	/** Returns whether subprocess decompression is allowed */
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	private static final int pickType(int type) {
		if(type==4 && !ALLOW_BF4) {type=0;}
		else if(type==3 && !ALLOW_BF3) {type=0;}
		else if(type==2 && !ALLOW_BF2) {type=0;}
		else if(type==1 && !ALLOW_BF1) {type=0;}
		if(type>0) {return type;}
		assert(type==0) : type;
		
		final int threads=Shared.threads();
		if(FORCE_MODE_BF1) {return 1;}
		if(FORCE_MODE_BF4) {return 4;}
		if(FORCE_MODE_BF3) {return 3;}
		if(FORCE_MODE_BF2) {return 2;}
		
		if(Shared.LOW_MEMORY || threads<12) {return 1;}
		if(ALLOW_BF4) {return 4;}
		if(ALLOW_BF3) {return 3;}
		if(ALLOW_BF2) {return 2;}
		return 1;
	}
	
	/** The FileFormat describing this file's characteristics and location */
	public final FileFormat ff;
	
	/** Force usage of ByteFile1 */
	public static boolean FORCE_MODE_BF1=false;
	/** Force usage of ByteFile2 implementation regardless of system settings */
	public static boolean FORCE_MODE_BF2=false;
	/** Unused legacy flag for ByteFile3 implementation */
	public static boolean FORCE_MODE_BF3=false;
	public static boolean FORCE_MODE_BF4=false;

	public static boolean ALLOW_BF1=true;
	public static boolean ALLOW_BF2=true;
	public static boolean ALLOW_BF3=false;//Hung on fasta input
	public static boolean ALLOW_BF4=true;//Hangs if there are limited reads and it is never closed. 
	
	/** Carriage return character constant */
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';
	protected static final byte[] plusLine=new byte[] {plus};
	
	public static int TARGET_LIST_SIZE=800;
	public static int TARGET_LIST_BYTES=262144;
	
//	byte[] pushBack=null;
	protected long nextID=0;
	
}
