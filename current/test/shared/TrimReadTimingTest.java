package test.shared;

import java.util.Arrays;

import shared.TrimRead;
import stream.Read;
import stream.SamLine;

/** Matched timing of ordinary unmapped and perfect-SAM trimming.
 * Inputs and read/SAM containers are prepared once. Each call trims one base
 * from each end, retaining production allocation costs while excluding input
 * construction. Boundary failures belong to TrimReadMinimumTest, not this
 * equal-output performance comparison.
 * @author Sucrose
 */
public final class TrimReadTimingTest{

	/** Runs mode=unmapped|mapped with positive rounds and warmup counts. */
	public static void main(String[] args){
		long rounds=100000, warmup=10000;
		String mode="unmapped";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("rounds")){rounds=Long.parseLong(split[1]);}
			else if(split[0].equals("warmup")){warmup=Long.parseLong(split[1]);}
			else if(split[0].equals("mode")){mode=split[1];}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(rounds<1 || warmup<1 || (!mode.equals("unmapped") && !mode.equals("mapped"))){
			throw new IllegalArgumentException("Positive rounds/warmup and mode=unmapped|mapped required");
		}
		final TrimReadTimingTest test=new TrimReadTimingTest(mode.equals("mapped"));
		if(rounds>Long.MAX_VALUE/test.expected || warmup>Long.MAX_VALUE/test.expected){
			throw new IllegalArgumentException("Repetition count would overflow checksum");
		}
		require(test.loop(warmup)==test.expected*warmup, "warmup checksum");
		final long start=System.nanoTime();
		final long sum=test.loop(rounds);
		final long elapsed=System.nanoTime()-start;
		require(sum==test.expected*rounds, "measured checksum");
		System.out.println("TRIM_READ_TIMING mode="+mode+" rounds="+rounds+" warmup="+warmup+
				" calls="+(rounds*test.reads.length)+" nanos="+elapsed+" checksum="+sum);
	}

	/** Creates 64 valid reads of 150..213 bases, with literal common-path expectations. */
	private TrimReadTimingTest(boolean mapped){
		long total=0;
		for(int i=0; i<reads.length; i++){
			final int length=150+i;
			bases[i]=new byte[length];
			quals[i]=new byte[length];
			for(int j=0; j<length; j++){
				bases[i][j]=(byte)"ACGT".charAt((i+j)&3);
				quals[i][j]=(byte)(20+(j%20));
			}
			reads[i]=new Read(bases[i], quals[i], "timing", i, false);
			cigars[i]=length+"M";
			if(mapped){sams[i]=new SamLine();sams[i].flag=0;}
			// Two bases removed, length-2 retained, first retained base/quality at index1.
			total+=2+(length-2)+bases[i][1]+quals[i][1]+(mapped ? 101 : 0);
		}
		expected=total;
		// Validate complete outputs once outside either timed loop.
		loop(1);
		for(int i=0; i<reads.length; i++){
			final Read read=reads[i];
			require(Arrays.equals(read.bases, Arrays.copyOfRange(bases[i], 1, bases[i].length-1)), "base oracle");
			require(Arrays.equals(read.quality, Arrays.copyOfRange(quals[i], 1, quals[i].length-1)), "quality oracle");
			if(mapped){
				final SamLine sam=sams[i];
				final int length=bases[i].length-2;
				require(sam.pos==101 && (sam.cigar.equals(length+"M") || sam.cigar.equals(length+"=")), "SAM coordinate/CIGAR oracle");
				require(Arrays.equals(sam.seq, read.bases) && Arrays.equals(sam.qual, read.quality), "SAM arrays oracle");
			}
		}
	}

	/** Resets borrowed inputs and times the public trimming dispatcher, including its allocations. */
	private long loop(long rounds){
		long sum=0;
		for(long round=0; round<rounds; round++){
			for(int i=0; i<reads.length; i++){
				final Read read=reads[i];
				read.bases=bases[i];
				read.quality=quals[i];
				read.samline=sams[i];
				final SamLine sam=sams[i];
				if(sam!=null){
					sam.pos=100;
					sam.setCigar(cigars[i]);
					sam.setSeq(bases[i]);
					sam.setQual(quals[i]);
				}
				final int trimmed=TrimRead.trimByAmount(read, 1, 1, 1, false);
				sum+=trimmed+read.length()+read.bases[0]+read.quality[0]+(sam==null ? 0 : sam.pos);
			}
		}
		return sum;
	}

	/** Checks expectations even when JVM assertions are disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reused mutable read containers. */
	private final Read[] reads=new Read[64];
	/** Immutable input arrays retained for resetting each call. */
	private final byte[][] bases=new byte[64][], quals=new byte[64][];
	/** Optional reused forward-strand SAM containers. */
	private final SamLine[] sams=new SamLine[64];
	/** Prebuilt input CIGAR strings. */
	private final String[] cigars=new String[64];
	/** Independent sum for one complete corpus traversal. */
	private final long expected;
}
