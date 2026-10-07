package stream.bam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Repeated file round-trip diagnostic using BgzfOutputStreamMT and BgzfInputStreamMT.
 * Creates the input only when absent, using 40000 repetitions of a default-charset
 * encoded line containing twenty A characters and a newline. Existing input is reused
 * without checking its contents. Compressed and decoded outputs are overwritten and
 * retained; Files.mismatch compares the input and decoded files after each round-trip.
 *
 * Positional arguments select iterations (default 20), maximum workers (default 4),
 * and block size (default 65536). Lower bounds are 1, 1 and 256 respectively. The
 * default block size conflicts with the current writer assertion, as noted below.
 * This diagnostic also retains post-Java8 API calls; it is not Java8-compatible.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfA80k{

	//Historical stress intent: sweep workers and repetitions over a repetitive fixture.
	//The familiar 840 KB size assumes a one-byte encoding for the generated ASCII text;
	//reused files can contain different data. This is not incompressible-boundary coverage.

	/** Input filename, reused unchanged if already present. */
	private static final String INPUT_FILE="A80k.txt";
	/** Fixed compressed-output filename, overwritten each iteration. */
	private static final String COMPRESSED_FILE="A40k.txt.gz";
	/** Fixed decoded-output filename, overwritten each iteration. */
	private static final String OUTPUT_FILE="A80k.roundtrip.txt";
	/** Generated line encoded using the platform default charset. */
	private static final byte[] LINE_BYTES="AAAAAAAAAAAAAAAAAAAA\n".getBytes();
	/** Number of lines generated when the input is absent. */
	private static final int LINE_REPETITIONS=40_000;

	/** Runs each worker count from one through maxThreads and retains all three artifacts.
	 * Arguments are iterations, maximum workers and block size; invalid integers propagate.
	 * @throws IOException If fixture creation, a round-trip, comparison or size reporting fails */
	public static void main(String[] args) throws IOException{
		int iterations=args.length>0 ? Integer.parseInt(args[0]) : 20;
		int maxThreads=args.length>1 ? Integer.parseInt(args[1]) : 4;
		int blockSize=args.length>2 ? Integer.parseInt(args[2]) : 65536;
		//TODO: Probable bug [TestBgzfA80k#001] - default 65536 exceeds BgzfOutputStreamMT's
		//asserted DEFAULT_BLOCK_SIZE=65280. The default is rejected when that constructor
		//is reached under -ea; this diagnostic has no upper-bound adjustment.
		if(iterations<1){iterations=1;}
		if(maxThreads<1){maxThreads=1;}
		if(blockSize<256){blockSize=256;}

		System.out.println("BGZF A80k Stress Test");
		System.out.println("=====================");
		System.out.println("Iterations per thread count: "+iterations);
		System.out.println("Thread counts: 1.."+maxThreads);
		System.out.println("Block size:   "+blockSize+" bytes");
		System.out.println();

		generateInputFile();

		for(int threads=1; threads<=maxThreads; threads++){
			System.out.println("Thread count: "+threads);
			for(int iter=1; iter<=iterations; iter++){
				System.out.print("  Iteration "+iter+"/"+iterations+" ... ");
				runRoundTrip(threads, blockSize);
				System.out.println("OK");
			}
			System.out.println();
		}

		System.out.println("✅ Completed BGZF round-trip stress test.");

		// Leave artifacts for inspection, but ensure they exist.
		printFileSizes();
	}

	/** Writes the repetitive fixture only if INPUT_FILE does not already exist. */
	private static void generateInputFile() throws IOException{
		File file=new File(INPUT_FILE);
		if(file.exists()){return;}

		try(FileOutputStream fos=new FileOutputStream(file)){
			for(int i=0; i<LINE_REPETITIONS; i++){fos.write(LINE_BYTES);}
		}
	}

	/** Compresses, decodes, and reports the first differing byte through Files.mismatch. */
	private static void runRoundTrip(int threads, int blockSize) throws IOException{
		compressFile(INPUT_FILE, COMPRESSED_FILE, threads, blockSize);
		decompressFile(COMPRESSED_FILE, OUTPUT_FILE, threads);

		//TODO: Compatibility [TestBgzfA80k#002] - Path.of and Files.mismatch are unavailable
		//in Java8. Preserve exact comparison semantics if replacing these APIs.
		long mismatch=Files.mismatch(Path.of(INPUT_FILE), Path.of(OUTPUT_FILE));
		if(mismatch!=-1){throw new IOException("Round-trip mismatch at byte position "+mismatch);}
	}

	/** Copies positive input reads to an MT writer using the supplied level-6 block configuration. */
	private static void compressFile(String input, String output, int threads, int blockSize) throws IOException{
		try(FileInputStream fis=new FileInputStream(input);
			FileOutputStream fos=new FileOutputStream(output);
			BgzfOutputStreamMT bgzf=new BgzfOutputStreamMT(fos, threads, 6, blockSize)){

			byte[] buffer=new byte[32*1024];
			int n;
			while((n=fis.read(buffer))>=0){
				if(n>0){bgzf.write(buffer, 0, n);}
			}
			//Explicit finalization retained; try-with-resources closes the writer again.
			bgzf.close();
		}
	}

	/** Decodes through the MT reader into the overwritten output, copying positive reads. */
	private static void decompressFile(String input, String output, int threads) throws IOException{
		try(FileInputStream fis=new FileInputStream(input);
			BgzfInputStreamMT bgzf=new BgzfInputStreamMT(fis, threads);
			FileOutputStream fos=new FileOutputStream(output)){

			byte[] buffer=new byte[32*1024];
			int n;
			while((n=bgzf.read(buffer))>=0){
				if(n>0){fos.write(buffer, 0, n);}
			}
		}
	}

	/** Prints current sizes of the input, compressed output and decoded output. */
	private static void printFileSizes() throws IOException{
		Path input=Path.of(INPUT_FILE);
		Path compressed=Path.of(COMPRESSED_FILE);
		Path output=Path.of(OUTPUT_FILE);

		System.out.println("Artifacts:");
		System.out.println("  "+INPUT_FILE+"       "+Files.size(input)+" bytes");
		System.out.println("  "+COMPRESSED_FILE+"  "+Files.size(compressed)+" bytes");
		System.out.println("  "+OUTPUT_FILE+"      "+Files.size(output)+" bytes");
	}
}
