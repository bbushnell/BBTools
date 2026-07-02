package driver;

import fileIO.TextFile;

/**
 * Utility for counting non-blank lines in text files.
 * Simple command-line tool that takes a filename as argument and reports the line count.
 * @author Brian Bushnell
 */
public class LineCount {
	
	/**
	 * Program entry point that counts lines in the specified file.
	 * Opens the file, counts non-blank lines, and prints the result.
	 * @param args Command-line arguments where args[0] should be the input filename
	 */
	public static void main(String[] args){
		
		//n [LineCount] CLEAN dev-trivial (no .sh, no callers): only risk is args[0] unguarded → AIOOBE with 0 args. tf is
		//n properly closed (L21); countLines returns non-blank count. Nothing else to flag.
		TextFile tf=new TextFile(args[0], false);
		long lines=tf.countLines();
		tf.close();
		System.out.println(args[0]+" has "+lines+" non-blank lines.");
		
	}
	
}
