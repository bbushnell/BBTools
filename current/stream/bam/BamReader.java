package stream.bam;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import fileIO.ReadWrite;
import structures.BinaryByteWrapperLE;

/**
 * Reads little-endian scalar values and byte ranges from a caller-supplied stream.
 * Instances borrow the stream and do not close it. The 16-bit and 32-bit methods share an eight-byte
 * temporary array and a BinaryByteWrapperLE cursor, so one caller must control access.
 * Exact-length reads throw EOFException for either boundary EOF or a partial value;
 * they do not distinguish record boundaries. This helper does not validate BAM records.
 * Static sort-order helpers separately open and close their own files and inspect
 * only a bounded portion of the header, as documented on those methods.
 *
 * @author Brian Bushnell, Chloe, Isla
 * @date November 5, 2025
 */
public class BamReader{

	/**
	 * Borrows an input stream and creates reusable scalar-decoding scratch storage.
	 * @param in Non-null stream, already positioned at the next value to read
	 */
	//Historical caller examples below do not guarantee exclusive access for every caller;
	//the owner must serialize use of this instance and its borrowed stream.
	//NOT thread-safe by design: the readInt/Short methods share one 8-byte `temp` + `wrapper`, reused per
	//call (readFully into temp, wrapper.position(0), getInt/getShort). Safe because every caller drives one
	//BamReader on a single thread (BamIndexWriter, getBamSortOrder) - a per-stream reader, not a shared one.
	public BamReader(InputStream in){
		this.in=in;
		this.temp=new byte[8];
		this.wrapper=new BinaryByteWrapperLE(temp);
	}

	/**
	 * Consumes four bytes and interprets them as a signed little-endian integer.
	 * @return Signed 32-bit value
	 * @throws EOFException If fewer than four bytes remain
	 * @throws IOException If the underlying stream reports another read error
	 */
	public int readInt32() throws IOException{
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt();
	}

	/**
	 * Consumes four bytes and widens the unsigned little-endian value to a long.
	 * @return Value in the range 0 through 4294967295
	 * @throws EOFException If fewer than four bytes remain
	 * @throws IOException If the underlying stream reports another read error
	 */
	public long readUint32() throws IOException{
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt()&0xFFFFFFFFL;
	}

	/**
	 * Consumes two bytes and interprets them as a signed little-endian short.
	 * @return Signed 16-bit value
	 * @throws EOFException If fewer than two bytes remain
	 * @throws IOException If the underlying stream reports another read error
	 */
	public short readInt16() throws IOException{
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort();
	}

	/**
	 * Consumes two bytes and widens the unsigned little-endian value to an int.
	 * @return Value in the range 0 through 65535
	 * @throws EOFException If fewer than two bytes remain
	 * @throws IOException If the underlying stream reports another read error
	 */
	public int readUint16() throws IOException{
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort()&0xFFFF;
	}

	/**
	 * Consumes one byte without using the scalar scratch buffer.
	 * @return Value in the range 0 through 255
	 * @throws EOFException If the underlying stream returns end-of-stream
	 * @throws IOException If the underlying stream reports another read error
	 */
	public int readUint8() throws IOException{
		int b=in.read();
		if(b<0){throw new EOFException();}
		return b;
	}

	/**
	 * Allocates a byte array and fills it from the current stream position.
	 * @param n Nonnegative byte count; zero produces an empty array without reading
	 * @return Newly allocated array containing exactly n bytes
	 * @throws EOFException If the stream ends before n bytes are read
	 * @throws IOException If the underlying stream reports another read error
	 */
	public byte[] readBytes(int n) throws IOException{
		byte[] result=new byte[n];
		readFully(result, 0, n);
		return result;
	}

	/**
	 * Reads n bytes and decodes them as US-ASCII, without removing a terminator.
	 * Callers reading a NUL-terminated BAM name must exclude and consume that byte separately.
	 * @param n Nonnegative number of bytes to decode
	 * @return Newly decoded string
	 * @throws EOFException If the stream ends before n bytes are read
	 * @throws IOException If the underlying stream reports another read error
	 */
	public String readString(int n) throws IOException{
		byte[] bytes=readBytes(n);
		return new String(bytes, 0, n, java.nio.charset.StandardCharsets.US_ASCII);
	}

	/**
	 * Repeats reads until the requested destination range is full.
	 * A partial read followed by EOF leaves the bytes already obtained in the array
	 * and throws; neither the input position nor the destination is rolled back.
	 * @param array Destination array
	 * @param offset First destination index
	 * @param n Nonnegative count fitting the destination range; zero performs no read
	 * @throws EOFException If the stream ends before the range is filled
	 * @throws IOException If the underlying stream reports another read error
	 */
	private void readFully(byte[] array, int offset, int n) throws IOException{
		int total=0;
		while(total<n){
			int bytesRead=in.read(array, offset+total, n-total);
			if(bytesRead<0){throw new EOFException("Expected "+n+" bytes, got "+total);}
			total+=bytesRead;
		}
	}

	/**
	 * Infers a SAM file's sort order from its first decompressed line using the platform charset.
	 * Uses ReadWrite's compression-aware input selection (STR378), including gzip SAM.
	 * A line starting with "@HD" is searched for SO:coordinate, SO:queryname, then
	 * SO:unsorted substrings in that order; this is not field-by-field SAM validation.
	 * The internally opened reader is closed before returning.
	 * @param samPath Filename to inspect
	 * @return "coordinate", "queryname", or "unsorted" for a recognized substring;
	 * "unknown" for an @HD-prefixed line without one; otherwise "headerless"
	 * @throws RuntimeException Wrapping an exception from opening, reading, or closing
	 */
	public static String getSamSortOrder(String samPath){
		try(BufferedReader br=new BufferedReader(new InputStreamReader(ReadWrite.getInputStream(samPath, false, false)))){
			String line=br.readLine();
			if(line==null || !line.startsWith("@HD")){return "headerless";}
			//Parse @HD line for SO: tag
			if(line.contains("SO:coordinate")){return "coordinate";}
			if(line.contains("SO:queryname")){return "queryname";}
			if(line.contains("SO:unsorted")){return "unsorted";}
			return "unknown"; //Has @HD but no SO tag
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/**
	 * Infers sort order after checking BAM magic and reading at most 100 header-text bytes.
	 * Inspects only the first line or the available prefix of that line, decoded as
	 * US-ASCII. An @HD-prefixed fragment is searched for recognized SO substrings in
	 * coordinate, queryname, unsorted order. This does not validate the full header,
	 * reference dictionary, or alignment records. Internally opened streams are closed.
	 * @param bamPath BAM filename to inspect through the BGZF input wrapper
	 * @return A recognized order, "unknown" for an unrecognized @HD-prefixed fragment,
	 * or "headerless" for absent header text or a fragment without the @HD prefix
	 * @throws RuntimeException Wrapping a bad-magic, opening, reading, or closing exception
	 */
	public static String getBamSortOrder(String bamPath){
		try(FileInputStream fis=new FileInputStream(bamPath);
			BufferedInputStream bis=new BufferedInputStream(fis, 1024);
			BgzfInputStream bgzf=new BgzfInputStream(bis)){

			BamReader reader=new BamReader(bgzf);

			//Validate BAM magic
			byte[] magic=reader.readBytes(4);
			if(magic[0]!='B' || magic[1]!='A' || magic[2]!='M' || magic[3]!=1){throw new IOException("Not a BAM file: "+bamPath);}

			//Read header text
			long lText=reader.readUint32();
			if(lText==0){return "headerless";}

			//Inspects only the first min(lText,100) header bytes to find the @HD line - same window as
			//BamIndexWriter. A @HD line >100 bytes (many optional fields before SO:) could push SO: past the
			//window -> misreported as "unknown"/sort-order missed. LOW/edge (most @HD lines are <100 bytes).
			int checkLen=(int)Math.min(lText, 100);
			byte[] headerStart=reader.readBytes(checkLen);

			//Find first line
			int newline=0;
			while(newline<headerStart.length && headerStart[newline]!='\n'){newline++;}
			String firstLine=new String(headerStart, 0, newline, java.nio.charset.StandardCharsets.US_ASCII);

			if(!firstLine.startsWith("@HD")){return "headerless";}
			if(firstLine.contains("SO:coordinate")){return "coordinate";}
			if(firstLine.contains("SO:queryname")){return "queryname";}
			if(firstLine.contains("SO:unsorted")){return "unsorted";}
			return "unknown";
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed stream; instance methods never close it. */
	private final InputStream in;
	/** Shared eight-byte scratch storage used for 16-bit and 32-bit reads. */
	private final byte[] temp;
	/** Little-endian cursor over temp, reset before decoding each 16-bit or 32-bit value. */
	private final BinaryByteWrapperLE wrapper;
}
