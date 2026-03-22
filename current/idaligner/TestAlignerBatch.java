package idaligner;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import aligner.AlignRandom;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Batch benchmark for long sequences (e.g. 40kbp).
 * Pre-generates all sequence pairs, runs Glocal first to establish truth,
 * then runs each other aligner against the same pairs.
 * Uses atomic work-stealing for efficient multithreading.
 *
 * Output columns match TestAlignerSuite:
 *   Name  ANI  rStart  rStop  Loops  Space%  Time  ANI_err  rStart_err  rStop_err
 *
 * @author Brian Bushnell
 * @contributor Neptune (Opus)
 * @date March 22, 2026
 */
public class TestAlignerBatch {

	public static void main(String[] args) {
		int samplesPerANI = 10;
		int length = 40000;
		int threads = Shared.threads();
		int sinewaves = 0;
		boolean subsOnly = false;
		boolean equalRates = false;
		long seed = 12345;
		String aniList = null; // comma-separated ANI values

		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			String[] split = arg.split("=");
			String a = split[0].toLowerCase();
			String b = split.length>1 ? split[1] : null;

			if(a.equals("samples") || a.equals("samplesperani") || a.equals("iterations") || a.equals("iters")) {
				samplesPerANI = Integer.parseInt(b);
			} else if(a.equals("length") || a.equals("len")) {
				length = Parse.parseIntKMG(b);
			} else if(a.equals("threads") || a.equals("t")) {
				threads = Integer.parseInt(b);
			} else if(a.equals("sinewaves") || a.equals("waves")) {
				sinewaves = Integer.parseInt(b);
			} else if(a.equals("subsonly") || a.equals("subs")) {
				subsOnly = Parse.parseBoolean(b);
			} else if(a.equals("equalrates") || a.equals("equal")) {
				equalRates = Parse.parseBoolean(b);
			} else if(a.equals("seed")) {
				seed = Long.parseLong(b);
			} else if(a.equals("ani") || a.equals("anis") || a.equals("anilist")) {
				aniList = b;
			}
		}

		// Default ANI levels: 100, 99.99, 99.98, 99.95, 99.9, 99.8, 99.5, 99, 98, 96...4
		float[] targetANIs;
		if(aniList != null) {
			String[] parts = aniList.split(",");
			targetANIs = new float[parts.length];
			for(int i=0; i<parts.length; i++) {
				targetANIs[i] = Float.parseFloat(parts[i].trim());
				if(targetANIs[i] > 1) { targetANIs[i] /= 100f; }
			}
		} else {
			// Build default list
			ArrayList<Float> list = new ArrayList<>();
			list.add(1.0f);
			for(float v : new float[]{99.99f, 99.98f, 99.95f, 99.9f, 99.8f, 99.5f}) {
				list.add(v / 100f);
			}
			for(float v = 99; v >= 4; v -= 2) {
				list.add(v / 100f);
			}
			targetANIs = new float[list.size()];
			for(int i=0; i<list.size(); i++) { targetANIs[i] = list.get(i); }
		}

		int mutMode = subsOnly ? 1 : (equalRates ? 2 : 0);
		int totalPairs = targetANIs.length * samplesPerANI;

		System.err.println("TestAlignerBatch: length=" + length + " samples=" + samplesPerANI
				+ " aniLevels=" + targetANIs.length + " totalPairs=" + totalPairs
				+ " threads=" + threads + " mutMode=" + mutMode + " seed=" + seed);

		// ============================
		// Phase 0: Generate all sequences
		// ============================

		Timer tGen = new Timer();
		shared.Random masterRandom = new rand.FastRandomXoshiro(seed);
		byte[] ref = AlignRandom.randomSequence(length, masterRandom);

		byte[][] queries = new byte[totalPairs][];
		int[] pairToANI = new int[totalPairs]; // maps pair index to ANI level index

		for(int a=0; a<targetANIs.length; a++) {
			for(int s=0; s<samplesPerANI; s++) {
				int idx = a * samplesPerANI + s;
				pairToANI[idx] = a;
				// Deterministic seed per pair for reproducibility
				shared.Random randy = new rand.FastRandomXoshiro(seed + a * 100000L + s);
				queries[idx] = TestAlignerSuite.mutateSequence(ref, targetANIs[a], randy, sinewaves, mutMode);
			}
		}
		tGen.stop();
		System.err.println("Generated " + totalPairs + " query sequences in " + tGen);

		// ============================
		// Phase 1: Run Glocal (truth)
		// ============================

		float[] trueANI = new float[totalPairs];
		int[] trueRStart = new int[totalPairs];
		int[] trueRStop = new int[totalPairs];
		long[] trueLoops = new long[totalPairs];
		long[] trueTime = new long[totalPairs]; // nanos

		IDAligner glocalTemplate = new GlocalPlusAligner5();
		System.err.println("Phase 1: Running Glocal+5 on all " + totalPairs + " pairs...");
		Timer tGlocal = new Timer();
		runAligner(glocalTemplate, ref, queries, trueANI, trueRStart, trueRStop,
				trueLoops, trueTime, threads);
		tGlocal.stop();
		System.err.println("Glocal+5 done in " + tGlocal);

		// ============================
		// Phase 2: Run each other aligner and compare to truth
		// ============================

		System.err.println(header());

		// Print Glocal results first
		printResults("Glocal+5", targetANIs, samplesPerANI, trueANI, trueRStart, trueRStop,
				trueLoops, trueTime, ref.length, queries,
				trueANI, trueRStart, trueRStop); // compare to self = zero error

		IDAligner[] templates = {
			new BandedAligner(), new DriftingAligner(), new WobbleAligner(),
			new ScrabbleAligner2(), new QuantumAligner(), new QuabbleAligner(),
			new XDropHAligner(), new WaveFrontAligner2()
		};

		for(IDAligner template : templates) {
			float[] ani = new float[totalPairs];
			int[] rStart = new int[totalPairs];
			int[] rStop = new int[totalPairs];
			long[] loops = new long[totalPairs];
			long[] time = new long[totalPairs];

			Timer tAligner = new Timer();
			runAligner(template, ref, queries, ani, rStart, rStop, loops, time, threads);
			tAligner.stop();
			System.err.println("  " + template.name() + " done in " + tAligner);

			printResults(template.name(), targetANIs, samplesPerANI, ani, rStart, rStop,
					loops, time, ref.length, queries,
					trueANI, trueRStart, trueRStop);
		}
	}

	/** Run one aligner on all pairs using atomic work-stealing. */
	static void runAligner(IDAligner template, byte[] ref, byte[][] queries,
			float[] aniOut, int[] rStartOut, int[] rStopOut,
			long[] loopsOut, long[] timeOut, int threads) {

		AtomicLong counter = new AtomicLong(0);
		int totalPairs = queries.length;

		Thread[] workers = new Thread[threads];
		for(int t=0; t<threads; t++) {
			final IDAligner ida = Test.createNewInstance(template);
			workers[t] = new Thread(() -> {
				int[] pos = new int[2];
				long idx;
				while((idx = counter.getAndIncrement()) < totalPairs) {
					int i = (int)idx;
					// Don't track loops — static counter races across threads
					ida.setLoops(-1);
					long startNanos = System.nanoTime();
					float id = ida.align(queries[i], ref, pos);
					long elapsed = System.nanoTime() - startNanos;

					// No synchronization needed — each index is unique
					aniOut[i] = id;
					rStartOut[i] = pos[0];
					rStopOut[i] = pos[1];
					loopsOut[i] = 0; // loops disabled in MT mode
					timeOut[i] = elapsed;
				}
			});
			workers[t].start();
		}

		for(Thread w : workers) {
			try { w.join(); } catch(InterruptedException e) { e.printStackTrace(); }
		}
	}

	/** Print per-ANI-level averaged results with per-pair error vs truth. */
	static void printResults(String name, float[] targetANIs, int samplesPerANI,
			float[] ani, int[] rStart, int[] rStop, long[] loops, long[] time,
			int refLen, byte[][] queries,
			float[] trueANI, int[] trueRStart, int[] trueRStop) {

		for(int a=0; a<targetANIs.length; a++) {
			double sumANI=0, sumRStart=0, sumRStop=0, sumLoops=0, sumTime=0;
			double sumANIerr=0, sumRStartErr=0, sumRStopErr=0;
			double sumStateSpace=0;

			for(int s=0; s<samplesPerANI; s++) {
				int idx = a * samplesPerANI + s;
				sumANI += ani[idx];
				sumRStart += rStart[idx];
				sumRStop += rStop[idx];
				sumLoops += loops[idx];
				sumTime += time[idx];
				sumStateSpace += (double)queries[idx].length * refLen;

				// Per-pair error vs truth
				sumANIerr += Math.abs(ani[idx] - trueANI[idx]);
				sumRStartErr += Math.abs(rStart[idx] - trueRStart[idx]);
				sumRStopErr += Math.abs(rStop[idx] - trueRStop[idx]);
			}

			double n = samplesPerANI;
			ByteBuilder bb = new ByteBuilder();
			bb.append(name);
			while(bb.length() < 9) { bb.append(' '); }
			bb.tab();
			bb.appendt(sumANI/n, 4);                          // avg ANI
			bb.appendt((int)(sumRStart/n));                    // avg rStart
			bb.appendt((int)(sumRStop/n));                     // avg rStop
			bb.appendt((long)(sumLoops/n));                    // avg loops
			double avgSpace = sumStateSpace / n;
			bb.appendt((sumLoops/n)/avgSpace*100, 3);          // avg space%
			bb.appendt(sumTime/n/1e9d, 9);                     // avg time (seconds)
			bb.appendt(sumANIerr/n, 9);                        // avg ANI error
			bb.appendt(sumRStartErr/n, 9);                     // avg rStart error
			bb.append(sumRStopErr/n, 9);                       // avg rStop error
			System.err.println(bb);
		}
	}

	static String header() {
		return "Name     \tANI\trStart\trStop\tLoops\tSpace%\tTime\tANI_err\trStart_err\trStop_err";
	}
}
