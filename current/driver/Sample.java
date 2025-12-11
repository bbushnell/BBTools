package driver;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/**
 * Simple file I/O example class that reads a text file line by line and writes
 * each line to an output file. Demonstrates basic BufferedReader and PrintWriter usage
 * for text file processing in BBTools.
 * @author Brian Bushnell
 */
public class Sample {
	
	/**
	 * Primary entry point for the Sample program. Reads command line arguments for input
	 * and output file names, creates readers/writers, and processes the data.
	 * @param args Command line arguments: args[0] = input filename, args[1] = output filename
	 */
	public static void main(String[] args){

		String fnameIn=args[0];
		String fnameOut=args[1];
		
		BufferedReader br=getReader(fnameIn);
		PrintWriter pw=getWriter(fnameOut);
		
		try {
			processData(br, pw);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		
	}
	
	/**
	 * Processes data by reading lines from input BufferedReader and writing each line
	 * to output PrintWriter. Currently performs a simple line-by-line copy operation.
	 * @param br BufferedReader for input file
	 * @param pw PrintWriter for output file
	 * @throws IOException if file reading or writing fails
	 */
	static void processData(BufferedReader br, PrintWriter pw) throws IOException{
		for(String s=br.readLine(); s!=null; s=br.readLine()){
			//Parsing goes here
			pw.println(s);
		}
	}
	
	/**
	 * Creates a BufferedReader for line-by-line text file reading. Handles FileInputStream
	 * creation and wraps with InputStreamReader for character encoding.
	 * @param fname Input filename
	 * @return BufferedReader for the specified file
	 */
	static BufferedReader getReader(String fname){
		FileInputStream fis=null;
		try {
			fis=new FileInputStream(fname);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		InputStreamReader isr=new InputStreamReader(fis);
		BufferedReader br=new BufferedReader(isr);
		return br;
	}
	
	/**
	 * Creates a PrintWriter for text output with buffering. Uses FileOutputStream wrapped
	 * with BufferedOutputStream for efficient writing.
	 * @param fname Output filename
	 * @return PrintWriter for the specified file
	 */
	static PrintWriter getWriter(String fname){
		FileOutputStream fos=null;
		try {
			fos=new FileOutputStream(fname);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		BufferedOutputStream bos=new BufferedOutputStream(fos);
		PrintWriter pw=new PrintWriter(bos);
		return pw;
	}
	
}
