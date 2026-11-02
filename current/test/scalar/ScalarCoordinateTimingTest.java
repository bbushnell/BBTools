package test.scalar;

import java.lang.reflect.Constructor;
import java.nio.charset.StandardCharsets;

import scalar.ScalarData;
import stream.Read;
import tracker.KmerTracker;

/** Matched single-thread timing of interval generation, with metric/weight checksum.
 * Reuses a one-million-base input and tracker; no files or external services.
 * @author Shinobu
 */
public final class ScalarCoordinateTimingTest{

	/** Runs the requested number of complete sequence scans after ten warmup scans. */
	public static void main(String[] args) throws Exception{
		int rounds=100, gap=0;
		for(String arg : args){
			if(arg.startsWith("rounds=")){rounds=Integer.parseInt(arg.substring(7));}
			else if(arg.startsWith("gap=")){gap=Integer.parseInt(arg.substring(4));}
			else{throw new IllegalArgumentException(arg);}
		}
		if(rounds<1 || gap<0){throw new IllegalArgumentException("rounds must be positive and gap nonnegative");}
		ScalarData.makeClade=ScalarData.makeSketch=ScalarData.parseTID=false;
		final byte[] bases=new byte[1000000], motif="ACGTTGCAAACTCGAG".getBytes(StandardCharsets.US_ASCII);
		for(int i=0; i<bases.length; i++){bases[i]=(gap>0 && i%gap==0 ? (byte)'N' : motif[i%motif.length]);}
		final Read read=new Read(bases, null, "timing", 0);
		final KmerTracker tracker=new KmerTracker(2, 50000);
		final Constructor<ScalarData> constructor=ScalarData.class.getDeclaredConstructor(boolean.class, long.class);
		constructor.setAccessible(true);
		for(int i=0; i<10; i++){measure(constructor, read, tracker);}
		long checksum=0;
		final long start=System.nanoTime();
		for(int i=0; i<rounds; i++){checksum+=measure(constructor, read, tracker);}
		final double seconds=(System.nanoTime()-start)*1e-9;
		System.out.println("rounds="+rounds+" gap="+gap+" bases="+((long)rounds*bases.length)+
				" seconds="+seconds+" checksum="+checksum);
	}

	/** Includes every original numeric metric and Length, excluding intentionally changed coordinates. */
	private static long measure(Constructor<ScalarData> constructor, Read read, KmerTracker tracker) throws Exception{
		final ScalarData data=constructor.newInstance(false, -1L);
		data.add(read, tracker, null, 10000, 500, -1, true);
		long sum=data.gc.size();
		for(int i=0; i<data.gc.size(); i++){
			sum=31*sum+Float.floatToIntBits(data.gc.get(i));
			sum=31*sum+Float.floatToIntBits(data.hh.get(i));
			sum=31*sum+Float.floatToIntBits(data.caga.get(i));
			sum=31*sum+Float.floatToIntBits(data.depth.get(i));
			sum=31*sum+Float.floatToIntBits(data.length.get(i));
		}
		assert(data.gc.size()>0) : "Million-base timing input must produce interval rows";
		return sum;
	}
}
