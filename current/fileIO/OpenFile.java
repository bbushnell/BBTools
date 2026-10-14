package fileIO;

import java.io.IOException;
import java.io.InputStream;

/**
 * Diagnostic that performs one raw-input read of at most 100 bytes.
 * Prints any returned bytes in single quotes to stderr using the default charset.
 * A short read is not retried, and empty input produces no output.
 * @author Brian Bushnell
 */
public class OpenFile{
	
	/**
	 * Opens raw input without buffering, reads once, and closes the acquired stream.
	 * Raw-input names supported by {@link ReadWrite#getRawInputStream(String, boolean)}
	 * are accepted, including standard input; this method closes standard input too.
	 * A read error remains the cause when close also fails, with the close error
	 * suppressed on it. Already printed bytes are not retracted after a close error.
	 * @param args Input name at index zero; additional arguments are ignored
	 * @throws RuntimeException If reading or closing fails, retaining the IOException
	 * as its cause; acquisition errors from the raw-input factory also propagate
	 */
	public static void main(String[] args){
		InputStream is=ReadWrite.getRawInputStream(args[0], false);
		byte[] line=new byte[100];
		try(InputStream input=is){
			int r=input.read(line, 0, 100);
			if(r>0){
				System.err.println("'"+new String(line, 0, r)+"'");
			}
		}catch(IOException e){
			throw new RuntimeException("Could not read or close input: "+args[0], e);
		}
	}
	
}
