package stream.bam;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Creates a repeating fixture, performs a BGZF round trip and compares recovered bytes.
 * This standalone diagnostic replaces fixed files in the working directory and leaves
 * them in place. The small repeating fixture is not exhaustive compression coverage.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfSimple{

	/**
	 * Creates test_simple.txt, compresses it to test_simple.txt.bgz, and decodes it to
	 * test_simple_decompressed.txt. Exits with status one if the paired-read comparison
	 * reports a mismatch. Prints suggested external gzip commands without executing them.
	 * Compression sets the JVM property bgzf.debug to true without restoring its old value.
	 * @param args Ignored
	 * @throws IOException If a file or codec operation throws an I/O exception
	 */
	//Original small-fixture rationale: the usual ASCII-compatible encoding is 21000
	//bytes, below a full-sized block. The actual fixture size depends on the default
	//charset; this harness does not establish general block-boundary coverage.
	public static void main(String[] args) throws IOException{
		System.out.println("BGZF Simple Data Test");
		System.out.println("=====================\n");

		// Create simple test data
		String testFile="test_simple.txt";
		String bgzfFile="test_simple.txt.bgz";

		createTestData(testFile);

		// Test 1: Compress with our compressor
		System.out.println("Test 1: Compressing with BgzfOutputStreamMT...");
		compressFile(testFile, bgzfFile);
		System.out.println("  Output: "+bgzfFile);
		System.out.println("  ✅ Compression complete\n");

		// Test 2: Decompress with our decompressor
		System.out.println("Test 2: Decompressing with BgzfInputStreamMT...");
		String decompressed="test_simple_decompressed.txt";
		decompressFile(bgzfFile, decompressed);

		// Verify
		if(filesMatch(testFile, decompressed)){System.out.println("  ✅ Round-trip successful! Files match perfectly.");}else{
			System.out.println("  ❌ MISMATCH! Files differ.");
			System.exit(1);
		}

		System.out.println("\n✅ All tests passed!");
		System.out.println("\nNow test with command-line gzip:");
		System.out.println("  gunzip -c "+bgzfFile+" > test_gzip_decompressed.txt");
		System.out.println("  diff "+testFile+" test_gzip_decompressed.txt");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Replaces a file with 1000 copies of twenty A characters followed by a newline.
	 * Encodes the line once with the default charset. The printed size is the fixed
	 * nominal 21000-byte size, not a measurement of the encoded output.
	 * @param filename Output path to replace
	 * @throws IOException If opening, writing or closing the file fails
	 */
	private static void createTestData(String filename) throws IOException{
		System.out.println("Creating test data: "+filename);

		try(FileOutputStream fos=new FileOutputStream(filename)){
			// Write "AAAAAAAAAAAAAAAAAAAA\n" 1000 times
			String line="AAAAAAAAAAAAAAAAAAAA\n";
			byte[] lineBytes=line.getBytes();

			for(int i=0; i<1000; i++){
				fos.write(lineBytes);
			}
		}

		System.out.println("  Created: 1000 lines of 'AAAAAAAAAAAAAAAAAAAA'");
		System.out.println("  Size: "+(21*1000)+" bytes\n");
	}

	/**
	 * Copies input through BgzfOutputStreamMT using an 8192-byte buffer.
	 * Sets bgzf.debug to true and leaves that JVM property set. The writer receives
	 * constructor arguments one and six; this method owns and closes its streams.
	 * @param input Input file to read
	 * @param output Compressed output path to replace
	 * @throws IOException If a file or codec operation throws an I/O exception
	 */
	private static void compressFile(String input, String output) throws IOException{
		// Set the debug property for the codec; retain this process-wide setting.
		System.setProperty("bgzf.debug", "true");

		try(FileInputStream fis=new FileInputStream(input);
			FileOutputStream fos=new FileOutputStream(output);
			BgzfOutputStreamMT bgzf=new BgzfOutputStreamMT(fos, 1, 6)){

			byte[] buffer=new byte[8192];
			int n;
			long total=0;

			while((n=fis.read(buffer))>=0){
				bgzf.write(buffer, 0, n);
				total+=n;
			}

			System.out.println("  Compressed "+total+" bytes");
			// Resource cleanup invokes writer close; finalization belongs to the codec.
		}
	}

	/**
	 * Copies decoded bytes from BgzfInputStreamMT to a replacement output file.
	 * Uses an 8192-byte buffer, reports the returned byte total and closes its streams.
	 * @param input Compressed file to read
	 * @param output Decoded output path to replace
	 * @throws IOException If a file or codec operation throws an I/O exception
	 */
	private static void decompressFile(String input, String output) throws IOException{
		try(FileInputStream fis=new FileInputStream(input);
			BgzfInputStreamMT bgzf=new BgzfInputStreamMT(fis, 1);
			FileOutputStream fos=new FileOutputStream(output)){

			byte[] buffer=new byte[8192];
			int n;
			long total=0;

			while((n=bgzf.read(buffer))>=0){
				fos.write(buffer, 0, n);
				total+=n;
			}

			System.out.println("  Decompressed "+total+" bytes");
		}
	}

	/**
	 * Compares paired reads using two 8192-byte buffers and closes both file streams.
	 * Rejects unequal read lengths immediately; it does not reconcile different read
	 * segmentation. For equal lengths, compares each returned byte before continuing.
	 * @param file1 First file to compare
	 * @param file2 Second file to compare
	 * @return true when all paired counts and bytes match through simultaneous EOF
	 * @throws IOException If opening, reading or closing either file fails
	 */
	private static boolean filesMatch(String file1, String file2) throws IOException{
		try(FileInputStream fis1=new FileInputStream(file1);
			FileInputStream fis2=new FileInputStream(file2)){

			byte[] buf1=new byte[8192];
			byte[] buf2=new byte[8192];
			long pos=0;

			while(true){
				int n1=fis1.read(buf1);
				int n2=fis2.read(buf2);

				if(n1!=n2){
					System.err.println("  Read size mismatch at position "+pos+": "+n1+" vs "+n2);
					return false;
				}

				if(n1<0){break;} // EOF

				for(int i=0; i<n1; i++){
					if(buf1[i]!=buf2[i]){
						System.err.println("  Byte mismatch at position "+(pos+i));
						return false;
					}
				}

				pos+=n1;
			}

			return true;
		}
	}
}
