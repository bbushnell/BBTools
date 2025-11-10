package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import shared.LineParser1;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Helper class for writing BAM binary structures with little-endian byte order.
 * All multi-byte integers in BAM format are little-endian.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamWriterHelper {

	public BamWriterHelper(OutputStream out) {
		this.out = out;
	}
	
	/**
	 * Write BAM header from SAM header lines.
	 * Uses efficient binary operations with ByteBuilder and LineParser.
	 * 
	 * @param headerLines SAM header lines (as byte arrays)
	 * @param supressHeader If true, skip writing header
	 * @param supressSequences If true, skip writing reference sequences
	 * @throws IOException
	 */
	public void writeHeaderFromLines(ArrayList<byte[]> headerLines, 
		boolean supressHeader, boolean supressSequences) throws IOException {

		if(supressHeader) {return;}

		// Extract reference names and lengths using LineParser
		ArrayList<byte[]> refNames = new ArrayList<byte[]>();
		IntList refLengths = new IntList();
		LineParser1 lp = new LineParser1('\t');

		ByteBuilder textBuilder = new ByteBuilder();

		for(byte[] line : headerLines) {
			// Append to header text
			textBuilder.append(line).nl();

			// Parse @SQ lines for reference info
			if(line.length > 3 && line[0] == '@' && line[1] == 'S' && line[2] == 'Q') {
				lp.set(line);
				byte[] sn = null;
				int ln = 0;

				for(int i = 1; i < lp.terms(); i++) {
					if(lp.termStartsWith("SN:", i)) {
						sn = lp.parseByteArray(i, 3);  // Skip "SN:" prefix
					} else if(lp.termStartsWith("LN:", i)) {
						ln = lp.parseInt(i, 3);  // Skip "LN:" prefix
					}
				}

				if(sn != null && ln > 0) {
					refNames.add(sn);
					refLengths.add(ln);
				}
			}
		}

		// Write magic "BAM\1"
		writeBytes(new byte[]{'B', 'A', 'M', 1});

		// Write header text section
		byte[] textBytes = textBuilder.toBytes();
		writeUint32(textBytes.length);
		writeBytes(textBytes);

		// Write reference dictionary section
		if(supressSequences) {
			writeInt32(0);
		} else {
			writeInt32(refNames.size());
			for(int i = 0; i < refNames.size(); i++) {
				byte[] nameBytes = refNames.get(i);
				writeUint32(nameBytes.length + 1);
				writeBytes(nameBytes);
				writeUint8(0);
				writeUint32(refLengths.get(i));
			}
		}
	}

	/**
	 * Write a 32-bit signed integer in little-endian format.
	 */
	public void writeInt32(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
		out.write((val >> 16) & 0xFF);
		out.write((val >> 24) & 0xFF);
	}

	/**
	 * Write a 32-bit unsigned integer (stored as long) in little-endian format.
	 */
	public void writeUint32(long val) throws IOException {
		writeInt32((int)val);
	}

	/**
	 * Write a 64-bit unsigned integer in little-endian format.
	 */
	public void writeUint64(long val) throws IOException {
		out.write((int)(val & 0xFF));
		out.write((int)((val >> 8) & 0xFF));
		out.write((int)((val >> 16) & 0xFF));
		out.write((int)((val >> 24) & 0xFF));
		out.write((int)((val >> 32) & 0xFF));
		out.write((int)((val >> 40) & 0xFF));
		out.write((int)((val >> 48) & 0xFF));
		out.write((int)((val >> 56) & 0xFF));
	}

	/**
	 * Write a 16-bit signed integer in little-endian format.
	 */
	public void writeInt16(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
	}

	/**
	 * Write a 16-bit unsigned integer in little-endian format.
	 */
	public void writeUint16(int val) throws IOException {
		writeInt16(val);
	}

	/**
	 * Write an 8-bit unsigned integer.
	 */
	public void writeUint8(int val) throws IOException {
		out.write(val & 0xFF);
	}

	/**
	 * Write a byte array.
	 */
	public void writeBytes(byte[] data) throws IOException {
		out.write(data);
	}

	/**
	 * Write a byte array subset.
	 */
	public void writeBytes(byte[] data, int off, int len) throws IOException {
		out.write(data, off, len);
	}

	/**
	 * Write a string as US-ASCII bytes (no null terminator).
	 */
	public void writeString(String s) throws IOException {
		out.write(s.getBytes("US-ASCII"));
	}

	/**
	 * Write a 32-bit IEEE float in little-endian format.
	 */
	public void writeFloat(float val) throws IOException {
		int bits = Float.floatToIntBits(val);
		writeInt32(bits);
	}

	private final OutputStream out;
}
