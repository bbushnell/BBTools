package idaligner;

import java.util.ArrayList;
import shared.Random;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;

import aligner.AlignRandom;
import bin.ConservationModel;
import dna.AminoAcid;
import parse.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Tests multiple aligners using random sequences.
 * The sequences have variable pairwise ANI, and each
 * ANI level is tested multiple times for average accuracy
 * and loop count.
 * @author Brian Bushnell
 * @contributor Opus
 * @date May 31, 2025
 */
public class TestAlignerSuite {

	public static void main(String[] args) {
		// Parse arguments
		int iterations = 32;
		float maxANI = 80;
		float minANI = 30;
		float step = 2;
		int length = 40000;
		int threads = Shared.threads();
		int sinewaves = 0; // 0 for uniform, >0 for conservation model
		boolean subsOnly = false; // true = subs only, no indels
		boolean equalRates = false; // true = 33/33/33 S/D/I

		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			String[] split = arg.split("=");
			String a = split[0].toLowerCase();
			String b = split.length>1 ? split[1] : null;

			if(a.equals("iterations") || a.equals("iters")) {
				iterations = Integer.parseInt(b);
			} else if(a.equals("maxani")) {
				maxANI = Float.parseFloat(b);
				if(maxANI<1) {maxANI*=100;}
			} else if(a.equals("minani")) {
				minANI = Float.parseFloat(b);
				if(minANI<1) {minANI*=100;}
			} else if(a.equals("step")) {
				step = Float.parseFloat(b);
				if(step<=0.05f) {step*=100;}
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
			}
		}

		// Create aligners to test
		ArrayList<IDAligner> aligners = new ArrayList<>();
//		aligners.add(new GlocalAligner());
		aligners.add(new GlocalPlusAligner5());
		aligners.add(new BandedAligner());
		aligners.add(new DriftingAligner());
		aligners.add(new WobbleAligner());
		aligners.add(new ScrabbleAligner());
		aligners.add(new QuantumAligner());
		aligners.add(new QuabbleAligner());
		aligners.add(new XDropHAligner());
		aligners.add(new WaveFrontAligner2());

		// Mutation mode: 0=default (75/12.5/12.5), 1=subsOnly, 2=equal (33/33/33)
		int mutMode = subsOnly ? 1 : (equalRates ? 2 : 0);

		// Run tests for each ANI value
		for(float ani = maxANI; ani >= (minANI-0.00001f); ani -= step) {
			System.err.println(Test.header());// Print header
			runTestsForANI(aligners, ani, iterations, length, threads, sinewaves, mutMode);
		}
	}

	private static void runTestsForANI(ArrayList<IDAligner> aligners, float targetANI,
	        int iterations, int length, int threads, int sinewaves, int mutMode){

	    // Create thread-safe accumulator for results
		final double[][] results=new double[aligners.size()][7]; // ani, rStart, rStop, loops, stateSpace, time, (unused)
	    
	    // Capture loop counts BEFORE running all jobs
	    long[] loopsBefore=new long[aligners.size()];
	    for(int a=0; a<aligners.size(); a++){
	        loopsBefore[a]=aligners.get(a).loops();
	    }
	    
	    // Create job queue
	    ConcurrentLinkedQueue<Job> jobs=new ConcurrentLinkedQueue<>();
	    for(int iter=0; iter<iterations; iter++){
	        jobs.add(new Job(iter, aligners, targetANI, length, sinewaves, mutMode, results));
	    }
	    
	    // Create and start worker threads
	    AtomicLong jobsProcessed=new AtomicLong(0);
	    Thread[] workers=new Thread[threads];
	    for(int t=0; t<threads; t++){
	        workers[t]=new Thread(new Worker(jobs, jobsProcessed));
	        workers[t].start();
	    }
	    
	    // Wait for all threads to complete
	    for(Thread t : workers){
	        while(t.isAlive()){
	            try{
	                t.join();
	            }catch(InterruptedException e){
	                e.printStackTrace();
	            }
	        }
	    }
	    
	    // Capture loop counts AFTER all jobs complete
	    for(int a=0; a<aligners.size(); a++){
	        long totalLoops=aligners.get(a).loops()-loopsBefore[a];
	        results[a][3]=totalLoops; // Store total, not per-iteration
	    }

	    // Print averaged results
	    for(int a=0; a<aligners.size(); a++){
	        IDAligner ida=aligners.get(a);
	        ByteBuilder bb=new ByteBuilder();

	        bb.append(ida.name());
	        while(bb.length()<9){bb.append(' ');}
	        bb.tab();
	        bb.appendt(results[a][0]/iterations, 4); // avg ANI
	        bb.appendt((int)(results[a][1]/iterations)); // avg rStart
	        bb.appendt((int)(results[a][2]/iterations)); // avg rStop
	        bb.appendt((long)(results[a][3]/iterations)); // avg loops (total/iterations)
	        float avgStateSpace=(float)(results[a][4]/iterations);
	        bb.appendt((results[a][3]/iterations)/avgStateSpace*100, 3); // avg space%
	        bb.append(results[a][5]/iterations/1e9d, 6); // avg time
	        
	        System.err.println(bb);
	    }
	}

	private static class Worker implements Runnable {

		Worker(ConcurrentLinkedQueue<Job> jobs_, AtomicLong jobsProcessed_) {
			jobs = jobs_;
			jobsProcessed = jobsProcessed_;
		}

		@Override
		public void run() {
			Job job;
			while((job = jobs.poll()) != null) {
				job.run();
				jobsProcessed.incrementAndGet();
			}
		}

		private final ConcurrentLinkedQueue<Job> jobs;
		private final AtomicLong jobsProcessed;
	}

	private static class Job {

		Job(int iteration_, ArrayList<IDAligner> aligners_, float targetANI_,
				int length_, int sinewaves_, int mutMode_, double[][] results_) {
			iteration = iteration_;
			// Create local copies of aligners to avoid thread conflicts
			aligners = new ArrayList<>();
			for(IDAligner ida : aligners_) {
				aligners.add(Test.createNewInstance(ida));
			}
			targetANI = targetANI_;
			length = length_;
			sinewaves = sinewaves_;
			mutMode = mutMode_;
			results = results_;
		}

		void run() {
			Random randy = shared.Shared.random();

			// Generate random reference sequence
			byte[] ref = AlignRandom.randomSequence(length, randy);

			// Mutate to create query at target ANI
			byte[] query = mutateSequence(ref, targetANI/100.0f, randy, sinewaves, mutMode);

			// Test each aligner
			for(int a = 0; a < aligners.size(); a++) {
				IDAligner ida = aligners.get(a);

				int[] pos = new int[2];
				long startTime = System.nanoTime();
				float id = ida.align(query, ref, pos);
				long time = System.nanoTime() - startTime;

				// Thread-safe accumulation
				float stateSpace = query.length * ref.length;
				synchronized(results[a]){
				    results[a][0]+=id; // ANI
				    results[a][1]+=pos[0]; // rStart
				    results[a][2]+=pos[1]; // rStop
				    // loops handled separately
				    results[a][4]+=stateSpace; // actual state space
				    results[a][5]+=time;
				}
			}
		}

		private final int iteration;
		private final ArrayList<IDAligner> aligners;
		private final float targetANI;
		private final int length;
		private final int sinewaves;
		private final int mutMode;
		private final double[][] results;
	}

	/**
	 * Mutate a sequence to achieve target identity
	 */
	public static byte[] mutateSequence(byte[] bases, float targetIdentity, Random randy, int sinewaves) {
		return mutateSequence(bases, targetIdentity, randy, sinewaves, 0);
	}

	public static byte[] mutateSequence(byte[] bases, float targetIdentity, Random randy, int sinewaves, int mutMode) {
		float errorRate = 1.0f - targetIdentity;
		float subRate, delRate, insRate;
		if(mutMode==1) { // subs only
			subRate = errorRate;
			delRate = 0;
			insRate = 0;
		} else if(mutMode==2) { // equal S/D/I
			subRate = errorRate / 3.0f;
			delRate = errorRate / 3.0f;
			insRate = errorRate - subRate - delRate;
		} else { // default: 75/12.5/12.5
			subRate = 0.75f * errorRate;
			delRate = 0.125f * errorRate;
			insRate = errorRate - subRate - delRate;
		}

		int maxIndel = 9;
		int indelSpacing = 0;
		boolean banHomopolymers = false;

		ByteBuilder bb = new ByteBuilder((int)(bases.length * 1.1f));

		ConservationModel conservator = (sinewaves < 1 ? null : 
			new ConservationModel(0.0f, sinewaves, randy));

		char prevChar = 'N';
		int lastIndel = -1;

		for(int i = 0; i < bases.length; ) {
			final byte b0 = bases[i];
			final boolean defined = AminoAcid.isFullyDefined(b0);
			final boolean spaceOK = (i - lastIndel > indelSpacing);
			byte b = b0;

			if(conservator != null && !conservator.shouldMutatePosition(i, randy)) {
				// Skip mutation, this is a conserved region
				bb.append(b);
				i++;
				continue;
			}

			float x = randy.nextFloat();
			boolean addSub = x <= subRate;
			boolean addDel = x > subRate && x <= subRate + delRate && spaceOK;
			boolean addIns = x > subRate + delRate && x <= errorRate && spaceOK;
			boolean addRef = (i == 0 || i == bases.length - 1 || 
					x > errorRate || !defined || !(addDel || addIns || addSub));
			boolean success = false;

			if(bb.length() > 0) {prevChar = (char)bb.get(bb.length() - 1);}

			if(addRef) {
				bb.append(b);
				success = true;
				i++;
			} else if(addSub) {
				b = AminoAcid.numberToBase[((AminoAcid.baseToNumber[b] + randy.nextInt3() + 1) & 3)];
				bb.append(b);
				success = true;
				i++;
			} else if(addIns) {
				int lim = Tools.min(maxIndel, bases.length - i - 2);
				int len = 1;
				if(lim >= 1) {
					len = 1 + Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim));
				}

				if(len > 0) {
					for(int j = 0; j < len; j++) {
						b = AminoAcid.numberToBase[randy.nextInt()&3];
						while(banHomopolymers && ((j == 0 && b == prevChar) || (j == len - 1 && b == b0))) {
							b = AminoAcid.numberToBase[randy.nextInt()&3];
						}
						bb.append(b);
					}
					success = true;
					lastIndel = i;
					// No i++ for insertions
				}
			} else if(addDel) {
				int lim = Tools.min(maxIndel, bases.length - i - 2);
				int len = 1;
				if(lim >= 1) {
					len = 1 + Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim));
				}

				if(len > 0) {
					// Skip the deletion
					i = i + len;
					lastIndel = i;
					success = true;
				}
			}

			if(!success) {
				// Problem encountered; advance
				bb.append(b);
				success = true;
				i++;
			}
		}

		return bb.toBytes();
	}
}