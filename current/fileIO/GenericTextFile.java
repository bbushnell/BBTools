package fileIO;

import java.util.ArrayList;

/**
 * Generic text file reader that extends TextFile with specialized line processing.
 * Provides convenient methods for reading all lines into an array and customized
 * line-by-line processing with optional filtering capabilities.
 * @author Brian Bushnell
 */
public class GenericTextFile extends TextFile {

	/**
	 * Constructs a GenericTextFile for reading from the specified file.
	 * Opens the file in read-only mode without write capabilities.
	 * @param name Path to the text file to read
	 */
	public GenericTextFile(String name) {
		super(name, false);
	}
	
	

	
	/**
	 * Reads all remaining lines from the file into a String array.
	 * Uses an ArrayList with initial capacity of 4096 for efficient growth.
	 * Each line is added sequentially until end of file is reached.
	 * @return Array containing all lines from the current position to end of file
	 */
	public String[] toLines(){
		
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list.toArray(new String[list.size()]);
		
	}
	
	@Override
	public String nextLine(){
		String line=readLine();
		while(line!=null && false){
			line=readLine();
		}
		return line;
	}
	

}
