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

public class AlignRandom {

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

	public static byte[] randomSequence(int len, Random randy) {
		byte[] array=new byte[len];
		for(int i=0; i<len; i++) {
			int x=randy.nextInt(4);
			byte b=AminoAcid.numberToBase[x];
			array[i]=b;
		}
		return array;
	}

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

	private static ByteBuilder header(int buckets) {
		float inv=1f/buckets;
		ByteBuilder bb=new ByteBuilder();
		bb.append("ANI");
		for(int i=0; i<=buckets; i++) {bb.tab().append(i*inv, 4);}
		return bb;
	}
	
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
