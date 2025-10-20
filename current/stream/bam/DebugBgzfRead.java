package stream.bam;

import java.io.FileInputStream;
import java.io.IOException;

/**
 * Debug tool to see exactly what happens when reading BGZF files.
 */
public class DebugBgzfRead {
	public static void main(String[] args) throws IOException {
		String file = args.length > 0 ? args[0] : "test_mt_output.bam";

		try (FileInputStream fis = new FileInputStream(file);
		     BgzfInputStreamMT in = new BgzfInputStreamMT(fis, 1)) {

			byte[] buf = new byte[8192];
			long total = 0;
			int blockNum = 0;

			while (true) {
				int n = in.read(buf);
				if (n < 0) {
					System.out.println("Block " + blockNum + ": EOF");
					break;
				}
				System.out.println("Block " + blockNum + ": " + n + " bytes (total: " + (total + n) + ")");
				total += n;
				blockNum++;

				// Stop after a few blocks near the end
				if (blockNum >= 374 && blockNum <= 380) {
					// Keep going
				} else if (blockNum > 374) {
					// We're past the interesting part
				}
			}

			System.out.println("\nTotal: " + total + " bytes in " + blockNum + " blocks");
		}
	}
}
