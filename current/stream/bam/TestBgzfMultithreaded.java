package stream.bam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Round-trip diagnostic using single-threaded and multithreaded BGZF streams.
 * First rewrites the original compressed input using an ST reader and MT writer,
 * then compares original and rewritten decoded streams using MT readers on both sides.
 * Next rewrites that first output with MT input/output and compares the two outputs,
 * again using the same MT decoder. This is not an independent decoded-byte oracle.
 * Fixed output filenames in the current directory are overwritten; deletion is attempted
 * after normal completion and its results are ignored. A failed run can leave artifacts.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfMultithreaded{

	//The first comparison includes the original compressed input, decoded by the MT reader.
	//Later comparison uses both outputs. Shared decoder behavior can mask correlated errors;
	//this remains a main-only diagnostic, not a broad correctness certification.

	/** Runs the three fixed stages with a positional filename and optional workers (default one).
	 * Missing input exits 1. Cleanup ignores File.delete() results.
	 * @throws IOException If a stream operation or decoded comparison fails */
	public static void main(String[] args) throws IOException{
		if(args.length<1){
			System.err.println("Usage: java stream.bam.TestBgzfMultithreaded <input.bam> [threads]");
			System.err.println("Example: java stream.bam.TestBgzfMultithreaded mapped.bam 1");
			System.exit(1);
		}

		String inputFile=args[0];
		int threads=args.length>1 ? Integer.parseInt(args[1]) : 1;
		String outputFile="test_mt_output.bam";

		System.out.println("BGZF Multithreaded Round-Trip Test");
		System.out.println("===================================");
		System.out.println("Input file:   "+inputFile);
		System.out.println("Output file:  "+outputFile);
		System.out.println("Threads:      "+threads);
		System.out.println();

		// Test 1: Read with single-threaded, write with multithreaded
		System.out.println("Test 1: Single-threaded read -> Multithreaded write");
		testReadWrite(inputFile, outputFile, false, true, threads);

		// Test 2: Read with multithreaded, compare byte-by-byte
		System.out.println("\nTest 2: Multithreaded read -> Byte comparison");
		testRoundTrip(inputFile, outputFile, threads);

		// Test 3: Full multithreaded round-trip
		System.out.println("\nTest 3: Full multithreaded round-trip");
		String output2="test_mt_output2.bam";
		testReadWrite(outputFile, output2, true, true, threads);
		testRoundTrip(outputFile, output2, threads);

		System.out.println("\n✅ All tests passed!");
		System.out.println("Multithreaded BGZF streams working correctly with "+threads+" thread(s)");

		// Cleanup
		new File(outputFile).delete();
		new File(output2).delete();
	}

	/** Copies decoded bytes through the selected ST/MT stream types into overwritten output.
	 * Explicitly closes the selected wrappers on normal completion; the raw file streams
	 * are also managed by try-with-resources. Reported throughput divides by 1024 twice. */
	private static void testReadWrite(String input, String output,
			boolean mtRead, boolean mtWrite, int threads) throws IOException{
		long startTime=System.currentTimeMillis();
		long bytesRead=0;

		try(FileInputStream fis=new FileInputStream(input);
			FileOutputStream fos=new FileOutputStream(output)){

			// Create input stream
			Object inStream=mtRead ?
				new BgzfInputStreamMT(fis, threads) :
				new BgzfInputStream(fis);

			// Create output stream
			Object outStream=mtWrite ?
				new BgzfOutputStreamMT(fos, threads, 6) :
				new BgzfOutputStream(fos, 6);

			// Copy data
			byte[] buffer=new byte[8192];
			int n;

			if(mtRead){
				BgzfInputStreamMT in=(BgzfInputStreamMT)inStream;
				if(mtWrite){
					BgzfOutputStreamMT out=(BgzfOutputStreamMT)outStream;
					while((n=in.read(buffer))>=0){
						if(n>0){
							out.write(buffer, 0, n);
							bytesRead+=n;
						}
					}
					out.close();
				}else{
					BgzfOutputStream out=(BgzfOutputStream)outStream;
					while((n=in.read(buffer))>=0){
						if(n>0){
							out.write(buffer, 0, n);
							bytesRead+=n;
						}
					}
					out.writeEOF();
					out.close();
				}
				in.close();
			}else{
				BgzfInputStream in=(BgzfInputStream)inStream;
				if(mtWrite){
					BgzfOutputStreamMT out=(BgzfOutputStreamMT)outStream;
					while((n=in.read(buffer))>=0){
						if(n>0){
							out.write(buffer, 0, n);
							bytesRead+=n;
						}
					}
					out.close();
				}else{
					BgzfOutputStream out=(BgzfOutputStream)outStream;
					while((n=in.read(buffer))>=0){
						if(n>0){
							out.write(buffer, 0, n);
							bytesRead+=n;
						}
					}
					out.writeEOF();
					out.close();
				}
				in.close();
			}
		}

		long elapsed=System.currentTimeMillis()-startTime;
		double mbps=(bytesRead/1024.0/1024.0)/(elapsed/1000.0);

		System.out.println("  Bytes transferred: "+bytesRead);
		System.out.println("  Time:             "+elapsed+" ms");
		System.out.println("  Throughput:       "+String.format("%.2f", mbps)+" MB/s");
	}

	/** Compares the decoded streams using two instances of the same MT reader.
	 * Checks read lengths and byte values; the printed block count counts bulk reads.
	 * @throws IOException If reading fails or any compared length/byte differs */
	//TODO: Validation limitation [stream/bam/TestBgzfMultithreaded#146] - both sides use
	//the same MT decoder. The first call does include the original compressed input,
	//but an independent decoder or known raw-byte oracle is needed for broader validation.
	private static void testRoundTrip(String file1, String file2, int threads) throws IOException{
		long startTime=System.currentTimeMillis();

		try(FileInputStream fis1=new FileInputStream(file1);
			FileInputStream fis2=new FileInputStream(file2);
			BgzfInputStreamMT in1=new BgzfInputStreamMT(fis1, threads);
			BgzfInputStreamMT in2=new BgzfInputStreamMT(fis2, threads)){

			byte[] buf1=new byte[8192];
			byte[] buf2=new byte[8192];
			long totalBytes=0;
			int blockNum=0;

			while(true){
				int n1=in1.read(buf1);
				int n2=in2.read(buf2);

				if(n1!=n2){
					throw new IOException("Read size mismatch at block "+blockNum+
						": file1="+n1+", file2="+n2);
				}

				if(n1<0){break;}// EOF

				// Compare bytes
				for(int i=0; i<n1; i++){
					if(buf1[i]!=buf2[i]){
						throw new IOException("Byte mismatch at block "+blockNum+
							", offset "+i+": file1="+(buf1[i]&0xFF)+
							", file2="+(buf2[i]&0xFF));
					}
				}

				totalBytes+=n1;
				blockNum++;
			}

			long elapsed=System.currentTimeMillis()-startTime;
			System.out.println("  Bytes compared:   "+totalBytes);
			System.out.println("  Blocks compared:  "+blockNum);
			System.out.println("  Time:             "+elapsed+" ms");
			System.out.println("  ✅ Files match perfectly!");
		}
	}
}
