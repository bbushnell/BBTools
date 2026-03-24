package idaligner;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import aligner.AlignRandom;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Benchmark search space vs sequence length at fixed ANI.
 * For each length, generates random sequence pairs, runs one aligner at a time
 * with multithreaded work-stealing. Loops are tracked per-swarm:
 * clear before launch, read after join, since only one aligner class runs at a time.
 *
 * Output columns:
 *   Name  Length  ANI  rStart  rStop  Loops  Space%  Time  ANI_err  rStart_err  rStop_err
 *
 * @author Brian Bushnell
 * @contributor Neptune (Opus)
 * @date March 23, 2026
 */
public class TestAlignerLength {

	public static void main(String[] args) {
		int samples = 100;
		int threads = Shared.threads();
		float targetANI = 0.75f;
		boolean subsOnly = false;
		boolean equalRates = false;
		long seed = 54321;
		String lengthList = null;

		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			String[] split = arg.split("=");
			String a = split[0].toLowerCase();
			String b = split.length>1 ? split[1] : null;

			if(a.equals("samples") || a.equals("iterations") || a.equals("iters")) {
				samples = Integer.parseInt(b);
			} else if(a.equals("threads") || a.equals("t")) {
				threads = Integer.parseInt(b);
			} else if(a.equals("ani")) {
				targetANI = Float.parseFloat(b);
				if(targetANI > 1) { targetANI /= 100f; }
			} else if(a.equals("subsonly") || a.equals("subs")) {
				subsOnly = Parse.parseBoolean(b);
			} else if(a.equals("equalrates") || a.equals("equal")) {
				equalRates = Parse.parseBoolean(b);
			} else if(a.equals("seed")) {
				seed = Long.parseLong(b);
			} else if(a.equals("lengths") || a.equals("lens") || a.equals("length") || a.equals("len")) {
				lengthList = b;
			}
		}

		int[] lengths;
		if(lengthList != null) {
			String[] parts = lengthList.split(",");
			lengths = new int[parts.length];
			for(int i=0; i<parts.length; i++) {
				lengths[i] = Parse.parseIntKMG(parts[i].trim());
			}
		} else {
			lengths = new int[]{64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536,
					2048, 3072, 4096, 6144, 8192, 12288, 16384, 24576, 32768, 49152, 65536};
		}

		int mutMode = subsOnly ? 1 : (equalRates ? 2 : 0);

		System.err.println("TestAlignerLength: ani=" + targetANI + " samples=" + samples
				+ " lengths=" + lengths.length + " threads=" + threads
				+ " mutMode=" + mutMode + " seed=" + seed);

		IDAligner[] templates = {
			new GlocalPlusAligner5(),
			new BandedAligner(), new DriftingAligner(), new WobbleAligner(),
			new ScrabbleAligner2(), new QuantumAligner(), new QuabbleAligner(),
			new XDropHAligner(), new WaveFrontAligner2()
		};

		System.err.println(header());

		for(int len : lengths) {
			System.err.println("--- len=" + len + " ---");

			// Generate sequences for this length
			shared.Random masterRandom = new rand.FastRandomXoshiro(seed + len);
			byte[] ref = AlignRandom.randomSequence(len, masterRandom);

			byte[][] queries = new byte[samples][];
			for(int s=0; s<samples; s++) {
				shared.Random randy = new rand.FastRandomXoshiro(seed + len * 100000L + s);
				queries[s] = TestAlignerSuite.mutateSequence(ref, targetANI, randy, 0, mutMode);
			}

			// Run Glocal first as truth
			float[] trueANI = new float[samples];
			int[] trueRStart = new int[samples];
			int[] trueRStop = new int[samples];

			IDAligner glocal = Test.createNewInstance(templates[0]);
			glocal.setLoops(0);
			runSwarm(templates[0], ref, queries, trueANI, trueRStart, trueRStop, threads);
			long glocalLoops = glocal.loops();
			printRow("Glocal+5", len, samples, trueANI, trueRStart, trueRStop,
					glocalLoops, ref.length, queries,
					trueANI, trueRStart, trueRStop);

			// Run each other aligner
			for(int a=1; a<templates.length; a++) {
				float[] ani = new float[samples];
				int[] rStart = new int[samples];
				int[] rStop = new int[samples];

				IDAligner instance = Test.createNewInstance(templates[a]);
				instance.setLoops(0);
				runSwarm(templates[a], ref, queries, ani, rStart, rStop, threads);
				long totalLoops = instance.loops();
				printRow(templates[a].name(), len, samples, ani, rStart, rStop,
						totalLoops, ref.length, queries,
						trueANI, trueRStart, trueRStop);
			}
		}
	}

	/**
	 * Run one aligner on all samples using multithreaded work-stealing.
	 * Loop counting works via static AtomicLong per aligner class —
	 * caller must clear loops before calling and read via any instance after.
	 */
	static void runSwarm(IDAligner template, byte[] ref, byte[][] queries,
			float[] aniOut, int[] rStartOut, int[] rStopOut, int threads) {

		AtomicLong counter = new AtomicLong(0);
		int n = queries.length;

		Thread[] workers = new Thread[threads];
		for(int t=0; t<threads; t++) {
			final IDAligner ida = Test.createNewInstance(template);
			workers[t] = new Thread(() -> {
				int[] pos = new int[2];
				long idx;
				while((idx = counter.getAndIncrement()) < n) {
					int i = (int)idx;
					float id = ida.align(queries[i], ref, pos);
					aniOut[i] = id;
					rStartOut[i] = pos[0];
					rStopOut[i] = pos[1];
				}
			});
			workers[t].start();
		}

		for(Thread w : workers) {
			try { w.join(); } catch(InterruptedException e) { e.printStackTrace(); }
		}
	}

	static void printRow(String name, int len, int samples,
			float[] ani, int[] rStart, int[] rStop,
			long totalLoops, int refLen, byte[][] queries,
			float[] trueANI, int[] trueRStart, int[] trueRStop) {

		double sumANI=0, sumRStart=0, sumRStop=0;
		double sumANIerr=0, sumRStartErr=0, sumRStopErr=0;
		double sumStateSpace=0;

		for(int s=0; s<samples; s++) {
			sumANI += ani[s];
			sumRStart += rStart[s];
			sumRStop += rStop[s];
			sumStateSpace += (double)queries[s].length * refLen;
			sumANIerr += Math.abs(ani[s] - trueANI[s]);
			sumRStartErr += Math.abs(rStart[s] - trueRStart[s]);
			sumRStopErr += Math.abs(rStop[s] - trueRStop[s]);
		}

		double n = samples;
		double avgLoops = (double)totalLoops / n;
		double avgSpace = sumStateSpace / n;

		ByteBuilder bb = new ByteBuilder();
		bb.append(name);
		while(bb.length() < 9) { bb.append(' '); }
		bb.tab();
		bb.appendt(len);
		bb.appendt(sumANI/n, 4);
		bb.appendt((int)(sumRStart/n));
		bb.appendt((int)(sumRStop/n));
		bb.appendt((long)avgLoops);
		bb.appendt(avgLoops/avgSpace*100, 3);
		bb.appendt(0, 9); // no per-sample timing in swarm mode
		bb.appendt(sumANIerr/n, 9);
		bb.appendt(sumRStartErr/n, 9);
		bb.append(sumRStopErr/n, 9);
		System.err.println(bb);
	}

	static String header() {
		return "Name     \tLength\tANI\trStart\trStop\tLoops\tSpace%\tTime\tANI_err\trStart_err\trStop_err";
	}
}
