package stream;

import java.io.IOException;
import java.io.OutputStream;

/**
 * An OutputStream that discards all data written to it.
 * Useful for suppressing output without performing I/O.
 * These write methods do not inspect their arguments or store data: null arrays
 * and invalid ranges are ignored rather than validated.
 * @author Brian Bushnell
 */
public class NullOutputStream extends OutputStream{
	
	/**
	 * Writes a single byte to the null stream (discards the data).
	 * @param b The byte to write (ignored)
	 * @throws IOException Never thrown by this implementation
	 */
	@Override
	public void write(final int b) throws IOException{}
	
	/**
	 * Writes an array of bytes to the null stream (discards the data).
	 * @param b The byte array to write (ignored)
	 * @throws IOException Never thrown by this implementation
	 */
	@Override
	public void write(final byte[] b) throws IOException{}
	
	/**
	 * Writes a portion of a byte array to the null stream (discards the data).
	 *
	 * @param b The byte array containing data (ignored)
	 * @param off Starting offset in the array (ignored)
	 * @param len Number of bytes to write (ignored)
	 * @throws IOException Never thrown by this implementation
	 */
	@Override
	public void write(final byte[] b, final int off, final int len) throws IOException{}
	
}
