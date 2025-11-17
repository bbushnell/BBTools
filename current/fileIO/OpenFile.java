package fileIO;

import java.io.IOException;
import java.io.InputStream;

/**
 * Simple file reading utility for testing raw input stream functionality.
 * Reads the first 100 bytes from a specified file and prints them to stderr.
 * @author Brian Bushnell
 */
public class OpenFile {
	
	/**
	 * Reads and displays first 100 bytes from specified file.
	 * Opens file using ReadWrite.getRawInputStream, reads up to 100 bytes,
	 * prints content to stderr, then closes stream.
	 * @param args Command line arguments, args[0] should be file path
	 */
	public static void main(String[] args){
		InputStream is=ReadWrite.getRawInputStream(args[0], false);
		byte[] line=new byte[100];
		try {
			int r=is.read(line, 0, 100);
			if(r>0){
				System.err.println("'"+new String(line, 0, r)+"'");
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
