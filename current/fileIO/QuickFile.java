package fileIO;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;


/**
 * Written for testing a NERSC slowdown in multithreaded file reading.
 * The problem was sometimes sidestepped by eliminating "if(pushBack!=null){" and reimplementing pushback.
 * However, that does not address the cause, so is not an overall solution; the cause remains a mystery.
 * This class may safely be deleted.
 * 
 * @author Brian Bushnell
 *
 */
public class QuickFile {
	
	
	/**
	 * Main method for testing file reading performance and functionality.
	 * Supports reading from stdin or specified file with optional line range limits.
	 * @param args Command line arguments: [filename] [start_line|"speedtest"] [end_line]
	 */
	public static void main(String[] args){
		QuickFile tf=new QuickFile(args.length>0 ? args[0] : "stdin", true);
		long first=0, last=100;
		boolean speedtest=false;
		if(args.length>1){
			if(args[1].equalsIgnoreCase("speedtest")){
				speedtest=true;
				first=0;
				last=Long.MAX_VALUE;
			}else{
				first=Integer.parseInt(args[1]);
				last=first+100;
			}
		}
		if(args.length>2){
			last=Integer.parseInt(args[2]);
		}
		speedtest(tf, first, last, !speedtest);
		
		tf.close();
		tf.reset();
		tf.close();
	}
	
	/**
	 * Performs speed test or line reading within specified range.
	 * @param tf QuickFile instance to read from
	 * @param first Starting line number (0-based)
	 * @param last Ending line number (exclusive)
	 * @param reprint If true, prints lines to stdout; if false, only counts performance
	 */
	private static void speedtest(QuickFile tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length;
				System.out.println(new String(s));
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
			}
		}
		t.stop();
		
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}
	
	/**
	 * Constructs QuickFile from filename with subprocess permission.
	 * @param fname Input filename or "stdin" for standard input
	 * @param allowSubprocess_ Whether to allow subprocess for compressed files
	 */
	public QuickFile(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}
	
	/** Constructs QuickFile from FileFormat specification.
	 * @param ff_ FileFormat object defining input source and properties */
	public QuickFile(FileFormat ff_){
		ff=ff_;
		assert(ff.read()) : ff;
		if(verbose){System.err.println("ByteFile1("+ff+")");}
		is=open();
	}
	
	/** Resets file reader to beginning, clearing pushback buffer and line counter.
	 * Closes and reopens the input stream to start from the beginning. */
	public final void reset(){
		close();
		is=open();

		pushBack=null;
		nextID=0;
	}
	
	/**
	 * Closes the input stream and cleans up resources.
	 * Thread-safe method that handles proper stream closure and error state tracking.
	 * @return Error state - true if errors occurred during reading/closing
	 */
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
//		assert(false) : name()+","+allowSubprocess();
		errorState|=ReadWrite.finishReading(is, name(), (allowSubprocess() || FileFormat.isBamFile(name())));
		
		is=null;
		lineNum=-1;
		pushBack=null;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}
	
	/**
	 * Reads next line from input as byte array, handling various line endings.
	 * Strips carriage returns and manages internal buffer for efficient reading.
	 * @return Next line as byte array, or null if end of file reached
	 */
	public byte[] nextLine(){
//		if(pushBack!=null){
//			byte[] temp=pushBack;
//			pushBack=null;
//			return temp;
//		}
		
		if(verbose){System.err.println("Reading line "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

//		System.out.println("\nCalled nextLine() for line "+lineNum);
//		System.out.println("A: bstart="+bstart+", bstop="+bstop);
		
		//if(bstart<bstop && lasteol==slasher && buffer[bstart]==slashn){bstart++;}
//		assert(bstart>=bstop || (buffer[bstart]!=slashn)/*buffer[bstart]>slasher || buffer[bstart]==slashn*/);
		int nlpos=bstart;
		
//		System.out.println("B: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
//		while(nlpos<bstop && (buffer[nlpos]>slasher || buffer[nlpos]==tab)){nlpos++;}
		while(nlpos<bstop && buffer[nlpos]!=slashn){nlpos++;}
//		System.out.println("C: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		if(nlpos>=bstop){
			nlpos=fillBuffer();
//			System.out.println("Filled buffer.");
		}
//		System.out.println("D: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		
		if(nlpos<0 || bstop<1){
			close();
			return null;
		}

		lineNum++;
		//Limit is the position after the last position to copy.
		//Limit equals nlpos unless there was a \r before the \n.
		final int limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		if(bstart==limit){//Empty line.
			bstart=nlpos+1;
//			System.out.println("E: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos+", returning='"+printNL(blankLine)+"'");
			return blankLine;
		}
		
		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);
		
		assert(line.length>0) : bstart+", "+nlpos+", "+limit;
		bstart=nlpos+1;
//		System.out.println("F: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos+", returning='"+printNL(line)+"'");
		return line;
	}
	
	/** Dummy byte array for testing purposes */
	final byte[] dummy=new byte[100];
	
	/**
	 * Fills internal buffer with data from input stream.
	 * Shifts remaining bytes to buffer start and reads new data to fill remainder.
	 * Automatically expands buffer size if needed to accommodate long lines.
	 * @return Position of next newline character, or -1 if end of stream
	 */
	private int fillBuffer(){
		if(bstart<bstop){ //Shift end bytes to beginning
//			System.err.println("Shift: "+bstart+", "+bstop);
			assert(bstart>0);
//			assert(bstop==buffer.length);
			int extra=bstop-bstart;
			for(int i=0; i<extra; i++, bstart++){
//				System.err.print((char)buffer[bstart]);
				//System.err.print('.');
				buffer[i]=buffer[bstart];
//				assert(buffer[i]>=slasher || buffer[i]==tab);
				assert(buffer[i]!=slashn);
			}
			bstop=extra;
//			System.err.println();

//			{//for debugging only
//				buffer=new byte[bufferlen];
//				bstop=0;
//				bstart=0;
//			}
		}else{
			bstop=0;
		}

		bstart=0;
		int len=bstop;
		int r=-1;
		while(len==bstop){//hit end of input without encountering a newline
			if(bstop==buffer.length){
//				assert(false) : len+", "+bstop;
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			try {
				r=is.read(buffer, bstop, buffer.length-bstop);
//				byte[] x=new byte[buffer.length-bstop];
//				r=is.read(x);
//				if(r>0){
//					for(int i=0, j=bstop; i<r; i++, j++){
//						buffer[j]=x[i];
//					}
//				}
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("open="+open);
			}
			if(r>0){
				bstop=bstop+r;
//				while(len<bstop && (buffer[len]>slasher || buffer[len]==tab)){len++;}
				while(len<bstop && buffer[len]!=slashn){len++;}
			}else{
				len=bstop;
				break;
			}
		}
		
//		System.err.println("After Fill: ");
//		printBuffer();
//		System.err.println();
		
//		System.out.println("Filled buffer; r="+r+", returning "+len);
		assert(r==-1 || buffer[len]==slashn);
		
//		System.err.println("lasteol="+(lasteol=='\n' ? "\\n" : lasteol==slashr ? "\\r" : ""+(int)lasteol));
//		System.err.println("First="+(int)buffer[0]+"\nLastEOL="+(int)lasteol);
		
		return len;
	}
	
	/**
	 * Opens input stream for reading, preventing double-opening.
	 * Thread-safe method that initializes buffer positions and stream state.
	 * @return InputStream for the file
	 * @throws RuntimeException if file is already open
	 */
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess());
		bstart=-1;
		bstop=-1;
		return is;
	}
	
	/**
	 * Reads all remaining lines into ArrayList of byte arrays.
	 * Continues reading until end of file is reached.
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
	 * Counts total number of lines in file and resets to beginning.
	 * Reads through entire file without storing line contents.
	 * @return Total number of lines in the file
	 */
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		
		return count;
	}
	
	/**
	 * Reads up to 200 lines into a numbered list structure.
	 * Thread-safe method for batch reading with automatic ID assignment.
	 * @return ListNum containing up to 200 lines, or null if no more data
	 */
	public synchronized final ListNum<byte[]> nextList(){
		byte[] line=nextLine();
		if(line==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>(200);
		list.add(line);
		for(int i=1; i<200; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
		}
		ListNum<byte[]> ln=new ListNum<byte[]>(list, nextID);
		nextID++;
		return ln;
	}
	
	/**
	 * Checks if the input file exists or is a special stream.
	 * Handles stdin, jar resources, and regular files.
	 * @return True if file exists or is accessible special stream
	 */
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists(); //TODO Ugly and unsafe hack for files in jars
	}
	
	/**
	 * Pushes line back to be returned by next nextLine() call.
	 * Can only hold one line at a time for simple lookahead functionality.
	 * @param line Line to push back for next read
	 */
	public final void pushBack(byte[] line){
		assert(pushBack==null);
		pushBack=line;
	}
	
	/** Returns the name of the input file or stream */
	public final String name(){return ff.name();}
	/** Returns whether subprocess usage is allowed for this file */
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	/** FileFormat object defining input source properties */
	public final FileFormat ff;

	/** Force mode flag for ByteFile1 compatibility testing */
	public static boolean FORCE_MODE_BF1=false;//!(Shared.GENEPOOL || Shared.DENOVO || Shared.CORI || Shared.WINDOWS);
	/** Force mode flag for ByteFile2 compatibility testing */
	public static boolean FORCE_MODE_BF2=false;
	/** Force mode flag for ByteFile3 compatibility testing */
	public static boolean FORCE_MODE_BF3=false;
	
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';

	/** Test variable for debugging purposes */
	long a=1;
	/** Test variable for debugging purposes */
	long b=2;
	/** Test variable for debugging purposes */
	long c=3;
	/** Test variable for debugging purposes */
	long d=4;
	/** Test pointer variable for debugging purposes */
	byte[] p0=null;
	/** Test pointer variable for debugging purposes */
	byte[] p1=null;
	/** Test pointer variable for debugging purposes */
	byte[] p2=null;
	/** Test pointer variable for debugging purposes */
	byte[] p3=null;
	/** Single line buffer for pushback functionality */
	private byte[] pushBack=null;
	/** Counter for assigning IDs to ListNum objects */
	private long nextID=0;
	
	/** Flag indicating whether input stream is currently open */
	private boolean open=false;
	/** Internal buffer for reading data from input stream */
	private byte[] buffer=new byte[bufferlen];
	/** Shared empty byte array for representing blank lines */
	private static final byte[] blankLine=new byte[0];
	private int bstart=0, bstop=0;
	/** Input stream for reading file data */
	public InputStream is;
	/** Current line number being read (0-based) */
	public long lineNum=-1;
	
	/** Enable verbose debug output for file operations */
	public static boolean verbose=false;
	/** Whether to use buffered input streams */
	public static boolean BUFFERED=false;
	/** Size of internal buffer for reading data */
	public static int bufferlen=16384;

	/** Flag tracking whether errors occurred during file operations */
	private boolean errorState=false;
	
}
