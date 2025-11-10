package illumina;

import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Quick inspector to examine CBCL binary file structure.
 * Helps understand the format before building the full parser.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclInspector {

	public static void main(String[] args) throws Exception {
		if(args.length<1){
			System.err.println("Usage: CbclInspector <cbcl_file>");
			System.exit(1);
		}

		String fname=args[0];
		inspectCbcl(fname);
	}

	static void inspectCbcl(String fname) throws Exception {
		System.out.println("Inspecting: " + fname);

		FileInputStream fis=new FileInputStream(fname);
		byte[] headerBytes=new byte[100]; //Read first 100 bytes
		int read=fis.read(headerBytes);
		fis.close();

		System.out.println("Read " + read + " bytes");
		System.out.println("\nFirst 100 bytes as hex:");
		for(int i=0; i<Math.min(read, 100); i++){
			System.out.printf("%02X ", headerBytes[i]);
			if((i+1)%16==0){System.out.println();}
		}
		System.out.println("\n");

		//Try parsing as CBCL header
		ByteBuffer bb=ByteBuffer.wrap(headerBytes);
		bb.order(ByteOrder.LITTLE_ENDIAN);

		short version=bb.getShort();
		int headerSize=bb.getInt();
		byte bitsPerBasecall=bb.get();
		byte bitsPerQscore=bb.get();

		System.out.println("Parsed header:");
		System.out.println("  Version: " + version);
		System.out.println("  Header size: " + headerSize);
		System.out.println("  Bits per basecall: " + bitsPerBasecall);
		System.out.println("  Bits per qscore: " + bitsPerQscore);
	}
}
