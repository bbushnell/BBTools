package stream.bam;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Copies decoded input through BgzfOutputStreamMT, then reports a decoded byte-count comparison.
 * Always opens debug_output.bam for replacement in the working directory.
 * Does not check whether input names that same file; opening output can truncate the input.
 * Does not compare contents: equal decoded lengths do not establish data equality.
 */
public class DebugBgzfWrite{
	/**
	 * Copies input with an 8192-byte buffer, then reads the output back and counts bytes.
	 * The historical "Block" messages label nonnegative read calls, not necessarily
	 * compressed blocks. "Bytes written" reports the decoded read-back total, not file size.
	 * A count mismatch is printed without explicitly throwing or setting failure status.
	 * After copying, closes the writer explicitly; resource cleanup also invokes close.
	 * File streams and readers are owned by their try-with-resources statements.
	 * @param args Optional input path, defaulting to mapped.bam; later arguments ignored
	 * @throws IOException If opening, reading, writing or closing a stream fails
	 */
	public static void main(String[] args) throws IOException{
		String input=args.length>0 ? args[0] : "mapped.bam";
		String output="debug_output.bam";

		long bytesRead=0;
		long bytesWritten=0;

		//TODO: Probable bug STR381 - input can name debug_output.bam. BgzfInputStream's
		//constructor reads no bytes, so opening this output truncates input before its first read.
		//Check that the files differ before opening output in a separately verified repair.
		try(FileInputStream fis=new FileInputStream(input);
			BgzfInputStream in=new BgzfInputStream(fis);
			FileOutputStream fos=new FileOutputStream(output);
			BgzfOutputStreamMT out=new BgzfOutputStreamMT(fos, 1, 6)){

			byte[] buffer=new byte[8192];
			int n;
			int blockNum=0;

			while((n=in.read(buffer))>=0){
				bytesRead+=n;
				out.write(buffer, 0, n);
				System.out.println("Block "+blockNum+": read "+n+" bytes (total read: "+bytesRead+")");
				blockNum++;
			}

			System.out.println("\nBefore close");
			out.close();
			System.out.println("After close");
		}

		// Now read back and count
		try(FileInputStream fis=new FileInputStream(output);
			BgzfInputStreamMT in=new BgzfInputStreamMT(fis, 1)){

			byte[] buffer=new byte[8192];
			int n;
			while((n=in.read(buffer))>=0){
				bytesWritten+=n;
			}
		}

		System.out.println("\n===================");
		System.out.println("Bytes read:    "+bytesRead);
		System.out.println("Bytes written: "+bytesWritten);
		System.out.println("Match: "+(bytesRead==bytesWritten ? "YES" : "NO - LOST "+(bytesRead-bytesWritten)));
	}
}
