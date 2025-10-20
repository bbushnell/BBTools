package stream.bam;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Debug tool to trace BGZF writing in detail.
 */
public class DebugBgzfWrite {
	public static void main(String[] args) throws IOException {
		String input = args.length > 0 ? args[0] : "mapped.bam";
		String output = "debug_output.bam";

		long bytesRead = 0;
		long bytesWritten = 0;

		try (FileInputStream fis = new FileInputStream(input);
		     BgzfInputStream in = new BgzfInputStream(fis);
		     FileOutputStream fos = new FileOutputStream(output);
		     BgzfOutputStreamMT out = new BgzfOutputStreamMT(fos, 1, 6)) {

			byte[] buffer = new byte[8192];
			int n;
			int blockNum = 0;

			while ((n = in.read(buffer)) >= 0) {
				bytesRead += n;
				out.write(buffer, 0, n);
				System.out.println("Block " + blockNum + ": read " + n + " bytes (total read: " + bytesRead + ")");
				blockNum++;
			}

			System.out.println("\nBefore close");
			out.close();
			System.out.println("After close");
		}

		// Now read back and count
		try (FileInputStream fis = new FileInputStream(output);
		     BgzfInputStreamMT in = new BgzfInputStreamMT(fis, 1)) {

			byte[] buffer = new byte[8192];
			int n;
			while ((n = in.read(buffer)) >= 0) {
				bytesWritten += n;
			}
		}

		System.out.println("\n===================");
		System.out.println("Bytes read:    " + bytesRead);
		System.out.println("Bytes written: " + bytesWritten);
		System.out.println("Match: " + (bytesRead == bytesWritten ? "YES" : "NO - LOST " + (bytesRead - bytesWritten)));
	}
}
