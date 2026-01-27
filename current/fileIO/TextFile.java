package fileIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * Robust utility class for reading and processing text files with flexible
 * input handling and advanced line processing capabilities.
 * Supports multiple input sources including files, standard input, JAR resources,
 * and compressed file formats with optional subprocess decompression.
 * Provides efficient memory management using BufferedReader with configurable
 * blank line handling and line-by-line processing.
 *
 * @author Brian Bushnell
 */
public class TextFile {
	
	
	/**
	 * Program entry point for testing TextFile functionality.
	 * Supports reading from files or stdin with optional speed testing.
	 * @param args Command-line arguments: [filename] [start_line|"speedtest"] [end_line]
	 */
	public static void main(String[] args){
		TextFile tf=new TextFile(args.length>0 ? args[0] : "stdin", false);
		int first=0;
		long last=100;
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
		
//		long lines=0;
//		long bytes=0;
//		if(args.length>1){
//			first=Integer.parseInt(args[1]);
//			last=first+100;
//		}
//		if(args.length>2){
//			last=Integer.parseInt(args[2]);
//		}
//		
//		for(int i=0; i<first; i++){tf.readLine();}
//		for(int i=first; i<last; i++){
//			String s=tf.readLine();
//			if(s==null){break;}
//
//			lines++;
//			bytes+=s.length();
//			System.out.println(s);
////			System.out.println(Arrays.toString(s.getBytes()));
//		}
//		
//		System.err.println("\n");
//		System.err.println("Lines: "+lines);
//		System.err.println("Bytes: "+bytes);
//		tf.close();
//		tf.reset();
//		tf.close();
//		
////		for(int i=first; i<last; i++){
////			String s=tf.readLine();
////			if(s==null){break;}
////
////			lines++;
////			bytes+=s.length();
////			System.out.println(s);
////		}
	}
	
	/**
	 * Performs speed testing of file reading operations with timing and throughput metrics.
	 * Optionally reprints lines or runs silent performance measurement.
	 *
	 * @param tf TextFile instance to test
	 * @param first First line number to start reading
	 * @param last Last line number to stop reading
	 * @param reprint Whether to print lines to output or run silently
	 */
	private static void speedtest(TextFile tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				String s=tf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length();
				System.out.println(s);
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				String s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length();
			}
		}
		t.stop();
		
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}

	/** Creates a TextFile with subprocess support disabled by default.
	 * @param name File path or "stdin" for standard input */
	public TextFile(String name){this(name, false);}
	
	/**
	 * Creates a TextFile from a FileFormat object with format-specific settings.
	 * Uses FileFormat's subprocess allowance and file path information.
	 * @param ff FileFormat containing file metadata and processing options
	 */
	public TextFile(FileFormat ff){
		file=new File(ff.name());
		allowSubprocess=ff.allowSubprocess();
		name=ff.name();
		
		br=open();
	}
	
	/**
	 * Creates a TextFile with configurable subprocess decompression support.
	 * Normalizes path separators and opens the file for reading.
	 * @param fname File path or "stdin" for standard input
	 * @param allowSubprocess_ Whether to allow subprocess decompression for compressed files
	 */
	public TextFile(String fname, boolean allowSubprocess_){
		fname=fname.replace('\\', '/');
		file=new File(fname);
		allowSubprocess=allowSubprocess_;
		name=fname;
		
		br=open();
	}
	
	/**
	 * Convenience method to read entire file into string array from FileFormat.
	 * Automatically closes the TextFile after reading.
	 * @param ff FileFormat specifying the file to read
	 * @return Array of all lines in the file as strings
	 */
	public static final String[] toStringLines(FileFormat ff){return toStringLines(ff, -1);}
	
	/**
	 * Convenience method to read entire file into string array from filename.
	 * Automatically closes the TextFile after reading.
	 * @param fname File path to read
	 * @return Array of all lines in the file as strings
	 */
	public static final String[] toStringLines(String fname){return toStringLines(fname, -1);}
	
	/**
	 * Convenience method to read entire file into string array from FileFormat.
	 * Automatically closes the TextFile after reading.
	 * @param ff FileFormat specifying the file to read
	 * @return Array of all lines in the file as strings
	 */
	public static final String[] toStringLines(FileFormat ff, int max){
		TextFile tf=new TextFile(ff);
		String[] lines=tf.toStringLines(max);
		tf.close();
		return lines;
	}
	
	/**
	 * Convenience method to read entire file into string array from filename.
	 * Automatically closes the TextFile after reading.
	 * @param fname File path to read
	 * @return Array of all lines in the file as strings
	 */
	public static final String[] toStringLines(String fname, int max){
		TextFile tf=new TextFile(fname);
		String[] lines=tf.toStringLines(max);
		tf.close();
		return lines;
	}
	
	/** Generate an array of the lines in this TextFile */
	public final String[] toStringLines(){
		return toStringLines(-1);
	}
	
	/** Generate an array of the lines in this TextFile */
	public final String[] toStringLines(long max){
		if(max<0) {max=Long.MAX_VALUE;}
		
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		
		for(s=nextLine(); s!=null && list.size()<max; s=nextLine()){
			list.add(s);
		}
		
		return list.toArray(new String[list.size()]);
		
	}
	
	/**
	 * Counts total number of lines in the file and resets to beginning.
	 * Performs full file traversal then automatically resets file position.
	 * @return Total number of lines in the file
	 */
	public final long countLines(){
		
		String s=null;
		long count=0;
		
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		
		reset();
		
		return count;
		
	}
	
	/**
	 * Splits an array of strings on tab characters to create 2D string array.
	 * Optionally trims whitespace before splitting each line.
	 *
	 * @param lines Array of strings to split
	 * @param trim Whether to trim whitespace before splitting
	 * @return 2D array where each row contains tab-separated fields from input line
	 */
	public static String[][] doublesplitTab(String[] lines, boolean trim){
		String[][] lines2=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			if(trim){
				lines2[i]=lines[i].trim().split("\t", -1);
			}else{
				lines2[i]=lines[i].split("\t", -1);
			}
		}
		return lines2;
	}
	
	
	/**
	 * Splits an array of strings on whitespace to create 2D string array.
	 * Uses Java whitespace pattern for splitting, optionally trims before splitting.
	 *
	 * @param lines Array of strings to split
	 * @param trim Whether to trim whitespace before splitting
	 * @return 2D array where each row contains whitespace-separated fields from input line
	 */
	public static String[][] doublesplitWhitespace(String[] lines, boolean trim){
		String[][] lines2=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			if(trim){
				lines2[i]=lines[i].trim().split("\\p{javaWhitespace}+");
			}else{
				lines2[i]=lines[i].split("\\p{javaWhitespace}+");
			}
		}
		return lines2;
	}
	
	/** Closes current file handle and reopens from the beginning.
	 * Allows reading the same file multiple times from the start. */
	public final void reset(){
		close();
		br=open();
	}
	
	/**
	 * Checks if the file or input source exists and is accessible.
	 * Handles special cases for stdin, JAR resources, and regular files.
	 * @return true if file exists or is a valid input source, false otherwise
	 */
	public boolean exists(){
		return name.equals("stdin") || name.startsWith("stdin.") || name.startsWith("jar:") || file.exists(); //TODO Ugly and unsafe hack for files in jars
	}
	
	/**
	 * Closes all file handles and streams, cleaning up resources.
	 * Sets open flag to false and nullifies all stream references.
	 * Updates error state based on stream closure success.
	 * @return Always returns false after closing
	 */
	public final boolean close(){
		if(!open){return false;}
		open=false;
		assert(br!=null);
		
		errorState|=ReadWrite.finishReading(is, name, allowSubprocess, br, isr);
		
		br=null;
		is=null;
		isr=null;
		lineNum=-1;
		return false;
	}
	
	/** Reads the next line from the file, skipping blank lines by default.
	 * @return Next non-blank line as string, or null if end of file */
	public String nextLine(){
		return readLine(true);
	}
	
	/** Reads the next line from the file, skipping blank lines by default.
	 * @return Next non-blank line as string, or null if end of file */
	public final String readLine(){
		return readLine(true);
	}
	
	/**
	 * Reads the next line from the file with configurable blank line handling.
	 * Increments line number counter and handles various error conditions.
	 * Blank line detection includes whitespace-only lines when skipBlank is true.
	 *
	 * @param skipBlank Whether to skip blank and whitespace-only lines
	 * @return Next line as string, or null if end of file or error
	 */
	public final String readLine(boolean skipBlank){
		String currentLine=null;
		
		
		//Note:  Disabling this block seems to speed things up maybe 5%.
//		boolean ready=false;
//		try {
//			ready=br.ready();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		if(!ready){return null;}
		
		if(!open || br==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name);}
			return null;
		}
		try{
			lineNum++;
			currentLine=br.readLine();
//			System.out.println(lineNum+":\t"+currentLine);
		}catch(Exception e){
			System.err.println("Oops! Bad read in file "+name+" at line "+lineNum);
			System.err.println(""+open+", "+(br==null));
			try {
				File f=new File(name);
				System.err.println("path and length: \t"+f.getAbsolutePath()+"\t"+f.length());
			} catch (Exception e1) {
				//e1.printStackTrace();
			}
			throw new RuntimeException(e);
		}
		if(currentLine==null){return null;}
//		System.out.println("Read "+line);
		
//		currentLine=currentLine.trim();
		
		//Note! This may generate a new String for every line and thus be slow.
//		if(currentLine.trim().length()==0){return readLine();} //Skips blank lines
		if(skipBlank && (currentLine.length()==0 ||
				(Character.isWhitespace(currentLine.charAt(0)) &&
						(Character.isWhitespace(currentLine.charAt(currentLine.length()-1)))) &&
						currentLine.trim().length()==0)){
			return readLine(skipBlank); //Skips blank lines
		}
		
		return currentLine;
	}
	
	/**
	 * Opens input stream and creates BufferedReader with 32KB buffer.
	 * Uses ReadWrite utility to handle various input sources and compression.
	 * Sets open flag and initializes all stream objects.
	 * @return Configured BufferedReader for the file
	 */
	private final BufferedReader open(){
		
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name);
		}
		open=true;
		
		is=ReadWrite.getInputStream(name, true, allowSubprocess);
		isr=new InputStreamReader(is);
		
		BufferedReader b=new BufferedReader(isr, 32768);
		
		return b;
	}
	
	/** Returns whether the file is currently open for reading */
	public boolean isOpen(){return open;}

	/** Flag indicating whether file is currently open for reading */
	private boolean open=false;
	/** Flag indicating if any errors occurred during file operations */
	public boolean errorState=false;
	
	/** File path or input source name */
	public final String name;
	/** File object representing the input file */
	public File file;
	/** Whether subprocess decompression is allowed for compressed files */
	private final boolean allowSubprocess;
	
	/** Underlying input stream for reading file data */
	public InputStream is;
	/** InputStreamReader wrapping the InputStream for character conversion */
	public InputStreamReader isr;
	/** BufferedReader providing efficient line-by-line reading with 32KB buffer */
	public BufferedReader br;
	
	/** Current line number being read, starts at -1 before first line */
	public long lineNum=-1;

	/** Global flag for verbose output during file operations */
	public static boolean verbose=false;
	
}
