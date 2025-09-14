package driver;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/**
 * @author Brian Bushnell
 * @date Oct 13, 2015
 *
 * This class will read a file and write it to another file.
 *
 */
public class Sample {
	
	/** Primary method, called by java */
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
	
	/** Do stuff */
	static void processData(BufferedReader br, PrintWriter pw) throws IOException{
		for(String s=br.readLine(); s!=null; s=br.readLine()){
			//Parsing goes here
			pw.println(s);
		}
	}
	
	/** Fetches a BufferedReader, which allows line-by-line String iteration over text files */
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
	
	/** Fetches a PrintWriter, which transforms Strings into a byte stream. */
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
