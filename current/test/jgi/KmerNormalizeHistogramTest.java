package test.jgi;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLongArray;

import bloom.KCountArray;
import stream.ConcurrentReadOutputStream;
import jgi.KmerNormalize;
import stream.Streamer;

/**
 * Literal private-helper oracle for {@code KmerNormalize.incrementHistogramSorted}.
 * The fixture constructs an unstarted worker, invokes only the histogram helper,
 * and verifies local bins, shared bins, conservation, and unchanged input arrays.
 * @author Jean
 */
public final class KmerNormalizeHistogramTest{

	public static void main(String[] args) throws Exception{
		final int originalThreadHistLen=KmerNormalize.THREAD_HIST_LEN;
		final int originalHistLen=KmerNormalize.HIST_LEN;
		final long originalHistLenPrint=KmerNormalize.HIST_LEN_PRINT;
		final AtomicLongArray originalKHistogram=KmerNormalize.khistogram;
		try{
			final Harness harness=new Harness();
			int failed=0;
			for(final Case c : CASES){
				try{
					if(!harness.run(c)){failed++;}
				}catch(Throwable t){
					failed++;
					System.err.println("FAIL\t"+c.name+"\t"+format(c.input)+"\t"+t);
					t.printStackTrace(System.err);
				}
			}
			if(failed>0){
				throw new AssertionError("KmerNormalize histogram fixture failed "+failed+
					" of "+CASES.length+" cases");
			}
			System.out.println("KMER_NORMALIZE_HISTOGRAM_PASS cases="+CASES.length);
		}finally{
			KmerNormalize.THREAD_HIST_LEN=originalThreadHistLen;
			KmerNormalize.HIST_LEN=originalHistLen;
			KmerNormalize.HIST_LEN_PRINT=originalHistLenPrint;
			KmerNormalize.khistogram=originalKHistogram;
		}
	}

	private static final class Harness{

		Harness() throws Exception{
			processThread=Class.forName("jgi.KmerNormalize$ProcessThread");
			constructor=processThread.getDeclaredConstructor(Streamer.class, KCountArray.class, KCountArray.class,
				int.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
				ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
				ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ArrayList.class);
			constructor.setAccessible(true);
			increment=processThread.getDeclaredMethod("incrementHistogramSorted", int[].class);
			increment.setAccessible(true);
			hist=processThread.getDeclaredField("hist");
			hist.setAccessible(true);
		}

		boolean run(final Case c) throws Exception{
			KmerNormalize.THREAD_HIST_LEN=LOCAL_LENGTH;
			KmerNormalize.HIST_LEN=HIST_LENGTH;
			KmerNormalize.HIST_LEN_PRINT=HIST_LENGTH;
			KmerNormalize.khistogram=new AtomicLongArray(HIST_LENGTH);

			final Object worker=constructor.newInstance(null, null, null, Integer.valueOf(1),
				null, null, null, null, null, null, null);
			final int[] observed=(c.input==null ? null : c.input.clone());
			try{
				increment.invoke(worker, (Object)observed);
			}catch(InvocationTargetException e){
				final Throwable cause=e.getCause();
				if(cause instanceof Exception){throw (Exception)cause;}
				if(cause instanceof Error){throw (Error)cause;}
				throw e;
			}

			final long[] local=(long[])hist.get(worker);
			final long[] shared=sharedBins();
			boolean ok=true;
			ok&=check(c, Arrays.equals(c.expectedLocal, local),
				"local expected="+Arrays.toString(c.expectedLocal)+" actual="+Arrays.toString(local));
			ok&=check(c, Arrays.equals(c.expectedShared, shared),
				"shared expected="+Arrays.toString(c.expectedShared)+" actual="+Arrays.toString(shared));
			ok&=check(c, sum(local)+sum(shared)==c.eligible,
				"conservation expected="+c.eligible+" actual="+(sum(local)+sum(shared)));
			ok&=check(c, Arrays.equals(c.input, observed),
				"input changed expected="+format(c.input)+" actual="+format(observed));
			System.out.println((ok ? "PASS" : "FAIL")+"\t"+c.name+"\tinput="+format(c.input)+
				"\tlocal="+Arrays.toString(local)+"\tshared="+Arrays.toString(shared)+
				"\teligible="+c.eligible);
			return ok;
		}

		private final Class<?> processThread;
		private final Constructor<?> constructor;
		private final Method increment;
		private final Field hist;
	}

	private static boolean check(final Case c, final boolean ok, final String message){
		if(!ok){System.err.println("CHECK_FAILED "+c.name+": "+message);}
		return ok;
	}

	private static long[] sharedBins(){
		final long[] bins=new long[HIST_LENGTH];
		for(int i=0; i<bins.length; i++){bins[i]=KmerNormalize.khistogram.get(i);}
		return bins;
	}

	private static long sum(final long[] array){
		long total=0;
		for(final long value : array){total+=value;}
		return total;
	}

	private static String format(final int[] array){
		return array==null ? "null" : Arrays.toString(array);
	}

	private static final class Case{

		Case(String name_, int[] input_, long[] local_, long[] shared_, long eligible_){
			name=name_;
			input=input_;
			expectedLocal=local_;
			expectedShared=shared_;
			eligible=eligible_;
		}

		final String name;
		final int[] input;
		final long[] expectedLocal;
		final long[] expectedShared;
		final long eligible;
	}

	private static final int LOCAL_LENGTH=3;
	private static final int HIST_LENGTH=5;
	private static final long[] ZERO_SHARED={0, 0, 0, 0, 0};
	private static final Case[] CASES={
		new Case("leading-negative-ascending-2", new int[]{-1, 2, 2},
			new long[]{0, 0, 2}, ZERO_SHARED, 2),
		new Case("trailing-negative-descending-2", new int[]{2, 2, -1},
			new long[]{0, 0, 2}, ZERO_SHARED, 2),
		new Case("leading-negative-ascending-capped", new int[]{-1, 0, 1, 2, 3, 4, 10},
			new long[]{1, 1, 1}, new long[]{0, 0, 0, 1, 2}, 6),
		new Case("trailing-negative-descending-capped", new int[]{10, 4, 3, 2, 1, 0, -1},
			new long[]{1, 1, 1}, new long[]{0, 0, 0, 1, 2}, 6),
		new Case("ascending-shared", new int[]{1, 2, 3, 4},
			new long[]{0, 1, 1}, new long[]{0, 0, 0, 1, 1}, 4),
		new Case("descending-shared", new int[]{4, 3, 2, 1},
			new long[]{0, 1, 1}, new long[]{0, 0, 0, 1, 1}, 4),
		new Case("null", null, new long[]{0, 0, 0}, ZERO_SHARED, 0),
		new Case("empty", new int[0], new long[]{0, 0, 0}, ZERO_SHARED, 0),
		new Case("all-negative", new int[]{-2, -1}, new long[]{0, 0, 0}, ZERO_SHARED, 0)
	};
}
