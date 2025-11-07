package stream.bam;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;

import structures.BinaryByteWrapperLE;

/**
 * Reads BAM binary structures with proper little-endian handling.
 * All multi-byte integers in BAM format are little-endian.
 * Optimized version using BinaryByteWrapperLE instead of ByteBuffer.
 *
 * @author Brian Bushnell, Chloe, Isla
 * @date November 5, 2025
 */
public class BamReader {

	public BamReader(InputStream in){
		this.in=in;
		this.temp=new byte[8];
		this.wrapper=new BinaryByteWrapperLE(temp);
	}

	/**
	 * Read a 32-bit signed integer (little-endian).
	 */
	public int readInt32() throws IOException {
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt();
	}

	/**
	 * Read a 32-bit unsigned integer as long (little-endian).
	 */
	public long readUint32() throws IOException {
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt()&0xFFFFFFFFL;
	}

	/**
	 * Read a 16-bit signed integer (little-endian).
	 */
	public short readInt16() throws IOException {
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort();
	}

	/**
	 * Read a 16-bit unsigned integer as int (little-endian).
	 */
	public int readUint16() throws IOException {
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort()&0xFFFF;
	}

	/**
	 * Read an 8-bit unsigned integer.
	 */
	public int readUint8() throws IOException {
		int b=in.read();
		if(b<0){
			throw new EOFException();
		}
		return b;
	}

	/**
	 * Read exactly n bytes.
	 */
	public byte[] readBytes(int n) throws IOException {
		byte[] result=new byte[n];
		readFully(result, 0, n);
		return result;
	}

	/**
	 * Read n bytes and interpret as ASCII string (no NUL terminator included).
	 */
	public String readString(int n) throws IOException {
		byte[] bytes=readBytes(n);
		return new String(bytes, 0, n, java.nio.charset.StandardCharsets.US_ASCII);
	}

	/**
	 * Read exactly n bytes into array at offset.
	 */
	private void readFully(byte[] array, int offset, int n) throws IOException {
		int total=0;
		while(total<n){
			int bytesRead=in.read(array, offset+total, n-total);
			if(bytesRead<0){
				throw new EOFException("Expected "+n+" bytes, got "+total);
			}
			total+=bytesRead;
		}
	}
	
	/**
	 * Determine the sort order of a SAM file from its header.
	 * @param samPath Path to SAM file
	 * @return "coordinate", "queryname", "unsorted", or "headerless" (if no @HD line)
	 * @throws IOException if file cannot be read
	 */
	public static String getSamSortOrder(String samPath){
		try(BufferedReader br=new BufferedReader(new FileReader(samPath))){
			String line=br.readLine();
			if(line==null || !line.startsWith("@HD")){
				return "headerless";
			}
			//Parse @HD line for SO: tag
			if(line.contains("SO:coordinate")){return "coordinate";}
			if(line.contains("SO:queryname")){return "queryname";}
			if(line.contains("SO:unsorted")){return "unsorted";}
			return "unknown"; //Has @HD but no SO tag
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	/**
	 * Determine the sort order of a BAM file from its header.
	 * @param bamPath Path to BAM file
	 * @return "coordinate", "queryname", "unsorted", or "headerless" (if no @HD line)
	 * @throws IOException if file cannot be read
	 */
	public static String getBamSortOrder(String bamPath){
		try(FileInputStream fis=new FileInputStream(bamPath);
			BufferedInputStream bis=new BufferedInputStream(fis, 1024);
			BgzfInputStream bgzf=new BgzfInputStream(bis)){
			
			BamReader reader=new BamReader(bgzf);
			
			//Validate BAM magic
			byte[] magic=reader.readBytes(4);
			if(magic[0]!='B' || magic[1]!='A' || magic[2]!='M' || magic[3]!=1){
				throw new IOException("Not a BAM file: "+bamPath);
			}
			
			//Read header text
			long lText=reader.readUint32();
			if(lText==0){return "headerless";}
			
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
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	private final InputStream in;
	private final byte[] temp;
	private final BinaryByteWrapperLE wrapper;
}