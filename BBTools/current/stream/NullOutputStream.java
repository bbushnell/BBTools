package stream;

import java.io.IOException;
import java.io.OutputStream;

/** Writes to nowhere.
 * Courtesy of https://stackoverflow.com/a/692580 and https://stackoverflow.com/a/691835 */
public class NullOutputStream extends OutputStream {
	
	@Override
	public void write(int b) throws IOException {}
	
	@Override
	public void write(byte[] b) throws IOException {}
	
	@Override
	public void write(byte[] b, int off, int len) throws IOException {}
	
}
