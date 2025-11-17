package driver;

import fileIO.TextFile;

/**
 * Simple text search utility that mimics basic grep functionality.
 * Reads a text file line by line and prints lines containing a specified search string.
 * @author Brian Bushnell
 */
public class Grep {
	
	/**
	 * Program entry point that performs text search on a file.
	 * Expects two command-line arguments: filename and search string.
	 * Prints all lines from the file that contain the search string.
	 * @param args Command-line arguments: args[0] = filename, args[1] = search string
	 */
	public static void main(String[] args){
		
		TextFile tf=new TextFile(args[0], true);
		
		String s=null;
		
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.contains(args[1])){System.out.println(s);}
		}
		tf.close();
		
	}
	
}
