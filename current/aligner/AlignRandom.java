package aligner;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Benchmarking tool for testing sequence alignment performance on random DNA sequences.
 * Generates pairs of random sequences of varying lengths and aligns them using GlocalPlusAligner5
 * to measure alignment identity distributions and performance characteristics.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class AlignRandom {

	/**
	 * Program entry point for alignment benchmarking.
	 * Runs alignment tests across multiple sequence lengths with specified parameters.
	 * @param args Command-line arguments: [min_length] [step_multiplier] [intervals] [iterations] [buckets] [max_loops] [output_file]
	 */
	public static void main(String[] args) {
		args=new PreParser(args, System.err, null, false, true, false).args;
		Shared.SIMD=true;
		int min=args.length<1 ? 10 : Integer.parseInt(args[0]);
		int step=args.length<2 ? 10 : Integer.parseInt(args[1]);
		int intervals=args.length<3 ? 4 : Integer.parseInt(args[2]);
		int iters=args.length<4 ? 200 : Integer.parseInt(args[3]);
		int buckets=args.length<5 ? 100 : Integer.parseInt(args[4]);
		long maxLoops=args.length<6 ? Long.MAX_VALUE : Parse.parseKMG(args[5]);
		String out=args.length<7 ? "stdout.txt" : args[6];
		
		final boolean mt=true;
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(out, true, false, false);
		bsw.println(header(buckets));
		for(int i=0, len=min; i<intervals; i++, len*=step) {
			long iters2=(maxLoops/len)/len;
			int iters3=(int)Tools.min(iters, iters2);
			System.err.println(len+", "+iters+", "+iters2+", "+iters3);
			int[] hist=(mt ? runMT(len, iters3, buckets) : runInterval(len, iters3, buckets));
			printInterval(len, iters3, hist, bsw);
		}
		bsw.poisonAndWait();
	}

	/**
	 * Runs single-threaded alignment benchmark for sequences of specified length.
	 * Generates random sequence pairs, aligns them, and records identity scores in histogram.
	 *
	 * @param len Length of random sequences to generate
	 * @param iters Number of alignment iterations to perform
	 * @param buckets Number of histogram buckets for identity scores
	 * @return Histogram array of alignment identity frequencies
	 */
	private static int[] runInterval(int len, int iters, int buckets) {
		int[] hist=new int[buckets+1];
		Random randy=Shared.threadLocalRandom(-1);
		for(int i=0; i<iters; i++) {
			byte[] a=randomSequence(len, randy);
			byte[] b=randomSequence(len, randy);
			float id=GlocalPlusAligner5.alignStatic(a, b, null);
			hist[Math.round(id*buckets)]++;
			
		}
		return hist;
	}

	/**
	 * Generates a random DNA sequence of specified length.
	 * Uses equal probability for all four standard bases (A, C, G, T).
	 *
	 * @param len Desired sequence length in bases
	 * @param randy Random number generator instance
	 * @return Random DNA sequence as byte array
	 */
	public static byte[] randomSequence(int len, Random randy) {
		byte[] array=new byte[len];
		for(int i=0; i<len; i++) {
			int x=randy.nextInt(4);
			byte b=AminoAcid.numberToBase[x];
			array[i]=b;
		}
		return array;
	}

	/**
	 * Prints alignment benchmark results for one sequence length interval.
	 * Outputs sequence length followed by normalized frequency values for each identity bucket.
	 *
	 * @param len Sequence length tested
	 * @param iters Total iterations performed
	 * @param hist Histogram of alignment identity frequencies
	 * @param bsw Output stream writer
	 */
	private static void printInterval(int len, int iters, int[] hist, ByteStreamWriter bsw) {
		assert(iters==Tools.sum(hist));
		float inv=1f/iters;
		ByteBuilder bb=new ByteBuilder();
		bb.append(len);
		for(int i=0; i<hist.length; i++) {
			bb.tab().append(hist[i]*inv, 5);
		}
		bsw.print(bb.nl());
	}

	/**
	 * Creates header row for benchmark output table.
	 * Contains "ANI" label followed by normalized identity bucket values.
	 * @param buckets Number of identity histogram buckets
	 * @return Header row as ByteBuilder
	 */
	private static ByteBuilder header(int buckets) {
		float inv=1f/buckets;
		ByteBuilder bb=new ByteBuilder();
		bb.append("ANI");
		for(int i=0; i<=buckets; i++) {bb.tab().append(i*inv, 4);}
		return bb;
	}
	
	/**
	 * Runs multi-threaded alignment benchmark for sequences of specified length.
	 * Distributes alignment tasks across thread pool using atomic counters for thread safety.
	 * Each thread generates its own random sequences and updates shared histogram.
	 *
	 * @param len Length of random sequences to generate
	 * @param iters Number of alignment iterations to perform
	 * @param buckets Number of histogram buckets for identity scores
	 * @return Histogram array of alignment identity frequencies
	 */
	static int[] runMT(final int len, final int iters, final int buckets) {
		ExecutorService executor = Executors.newFixedThreadPool(Shared.threads());
		final AtomicIntegerArray atomicHist = new AtomicIntegerArray(buckets + 1);

		List<Future<?>> futures = new ArrayList<>();
		for (int i = 0; i < iters; i++) {
		    futures.add(executor.submit(new Runnable() {
		        public void run() {
		            Random randy = new Random(); // Thread-local
		            byte[] a = randomSequence(len, randy);
		            byte[] b = randomSequence(len, randy);
		            float id = GlocalPlusAligner5.alignStatic(a, b, null);
		            int bucket = Math.round(id * buckets);
		            atomicHist.incrementAndGet(bucket);
		        }
		    }));
		}

		// Wait for completion
		for (Future<?> f : futures) { try {
			f.get();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} }
		executor.shutdown();

		// Convert to regular array
		int[] hist = new int[buckets + 1];
		for (int i = 0; i <= buckets; i++) {
		    hist[i] = atomicHist.get(i);
		}
		return hist;
	}
	
}
