package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Helper class for writing BAM binary structures with little-endian byte order.
 * All multi-byte integers in BAM format are little-endian.
 * Retains a caller-owned OutputStream without closing or flushing it. Calls are not
 * synchronized here; callers must coordinate the order of complete structures.
 * Primitive methods emit low-width bit patterns without numeric range validation.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamWriterHelper{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the supplied stream without wrapping, writing, closing or flushing it.
	 * @param out Output destination used by all subsequent writes
	 */
	public BamWriterHelper(OutputStream out){this.out=out;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes BAM magic, SAM header text and the selected binary reference dictionary.
	 * Appends every supplied line plus a newline to header text. Lines beginning with @SQ
	 * and longer than three bytes are parsed for SN/LN fields; repeated matches overwrite
	 * earlier values. A dictionary entry requires nonnull SN and positive LN. This is not
	 * general header validation, and text lines are retained even when no entry is selected.
	 * Sequence suppression retains header text but emits a zero-entry binary dictionary.
	 * All-header suppression returns before inspecting input or writing any bytes.
	 * @param headerLines Nonnull list of stable, nonnull byte arrays when output is enabled
	 * @param supressHeader Skip the entire operation when true
	 * @param supressSequences Emit no binary dictionary entries when true; text is unchanged
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeHeaderFromLines(ArrayList<byte[]> headerLines, 
		boolean supressHeader, boolean supressSequences) throws IOException{

		if(supressHeader){return;}

		// Extract reference names and lengths using LineParser
		ArrayList<byte[]> refNames=new ArrayList<byte[]>();
		IntList refLengths=new IntList();
		LineParser1 lp=new LineParser1('\t');

		ByteBuilder textBuilder=new ByteBuilder();

		for(byte[] line : headerLines){
			// Append to header text
			textBuilder.append(line).nl();

			// Parse @SQ lines for reference info
			if(line.length>3 && line[0]=='@' && line[1]=='S' && line[2]=='Q'){
				lp.set(line);
				byte[] sn=null;
				int ln=0;

				for(int i=1; i<lp.terms(); i++){
					if(lp.termStartsWith("SN:", i)){
						sn=lp.parseByteArray(i, 3);  // Skip "SN:" prefix
					}else if(lp.termStartsWith("LN:", i)){
						ln=lp.parseInt(i, 3);  // Skip "LN:" prefix
					}
				}

				//@SQ → ref dict: a ref is added iff it has BOTH a non-null SN AND a positive LN. A @SQ
				//missing SN or with LN<=0 is silently skipped (invalid ref). Matches the inline @SQ parse in
				//BamOutputStream.writeHeader (which uses ln!=null; here ln>0 is slightly stricter — a malformed
				//LN<=0 ref is dropped rather than emitted). refNames/refLengths stay index-aligned (both add or
				//neither). textBuilder always gets the raw line regardless (full @-header text is preserved).
				if(sn!=null && ln>0){
					refNames.add(sn);
					refLengths.add(ln);
				}
			}
		}

		// Write magic "BAM\1"
		writeBytes(new byte[]{'B', 'A', 'M', 1});

		// Write header text section
		byte[] textBytes=textBuilder.toBytes();
		writeUint32(textBytes.length);
		writeBytes(textBytes);

		// Write reference dictionary section
		if(supressSequences){
			writeInt32(0);
		}else{
			writeInt32(refNames.size());
			for(int i=0; i<refNames.size(); i++){
				byte[] nameBytes=refNames.get(i);
				writeUint32(nameBytes.length+1);
				writeBytes(nameBytes);
				writeUint8(0);
				writeUint32(refLengths.get(i));
			}
		}
	}

	/** Writes all 32 bits in little-endian order.
	 * @param val Integer bit pattern to emit
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeInt32(int val) throws IOException{
		out.write(val&0xFF);
		out.write((val>>8)&0xFF);
		out.write((val>>16)&0xFF);
		out.write((val>>24)&0xFF);
	}

	/** Narrows to the low 32 bits and writes them in little-endian order.
	 * No unsigned-range check is performed.
	 * @param val Value whose low 32 bits are emitted
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeUint32(long val) throws IOException{writeInt32((int)val);}

	/** Writes all 64 bits in little-endian order, including negative long bit patterns.
	 * @param val Long bit pattern representing the unsigned field
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeUint64(long val) throws IOException{
		out.write((int)(val&0xFF));
		out.write((int)((val>>8)&0xFF));
		out.write((int)((val>>16)&0xFF));
		out.write((int)((val>>24)&0xFF));
		out.write((int)((val>>32)&0xFF));
		out.write((int)((val>>40)&0xFF));
		out.write((int)((val>>48)&0xFF));
		out.write((int)((val>>56)&0xFF));
	}

	/** Writes the low 16 bits in little-endian order without a range check.
	 * @param val Integer supplying the two-byte bit pattern
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeInt16(int val) throws IOException{
		out.write(val&0xFF);
		out.write((val>>8)&0xFF);
	}

	/** Delegates unsigned-field emission to the same low-16-bit encoding as writeInt16.
	 * @param val Value whose low 16 bits are emitted
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeUint16(int val) throws IOException{writeInt16(val);}

	/** Writes the low eight bits without a range check.
	 * @param val Value whose low byte is emitted
	 * @throws IOException If the stream write fails
	 */
	public void writeUint8(int val) throws IOException{out.write(val&0xFF);}

	/** Writes the complete array without a length prefix or terminator.
	 * @param data Byte array passed directly to the stream
	 * @throws IOException If the stream write fails
	 */
	public void writeBytes(byte[] data) throws IOException{out.write(data);}

	/** Writes an array range without a length prefix or terminator.
	 * @param data Source byte array
	 * @param off Starting offset passed to the stream
	 * @param len Byte count passed to the stream
	 * @throws IOException If the stream write fails
	 */
	public void writeBytes(byte[] data, int off, int len) throws IOException{out.write(data, off, len);}

	/** Writes String.getBytes("US-ASCII") output without a length prefix or terminator.
	 * @param s String to encode
	 * @throws IOException If encoding or the delegated stream write fails
	 */
	public void writeString(String s) throws IOException{out.write(s.getBytes("US-ASCII"));}

	/** Writes Float.floatToIntBits output as a little-endian 32-bit field.
	 * Uses that conversion's NaN handling rather than preserving every raw NaN payload.
	 * @param val Floating-point value to encode
	 * @throws IOException If a delegated stream write fails
	 */
	public void writeFloat(float val) throws IOException{
		int bits=Float.floatToIntBits(val);
		writeInt32(bits);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained caller-owned destination; this helper never closes or flushes it. */
	private final OutputStream out;
}
