package cardinality;

import rand.FastRandomXoshiro;
import parse.Parse;

/**
 * Measures register read and write rates for all estimator types
 * at various cardinalities. Reports what fraction of adds result in
 * a register read (passed early exit) and a register write (actual update).
 *
 * @author Brian Bushnell, Nahida
 */
public class BranchRateTest {

	public static void main(String[] args){
		int memBytes=2048;
		long[] cards={100_000, 1_000_000, 5_000_000, 10_000_000, 20_000_000, 40_000_000, 80_000_000};

		for(int i=0; i<args.length; i++){
			final String arg=args[i].toLowerCase();
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String key=arg.substring(0, eq);
			final String val=arg.substring(eq+1);
			if(key.equals("mem")){memBytes=(int)Parse.parseKMGBinary(val);}
			else{throw new RuntimeException("Unknown parameter '"+args[i]+"'");}
		}

		final String[] types={"avll", "exa", "udll6m", "ull", "dll4", "ll6", "hll4", "hlll", "htc4", "htb"};
		final String[] labels={"AVLL", "EXA", "UDLL", "ULL", "DLL4", "LL6", "HLL4", "HLLL", "HTC4", "HTB"};

		System.out.println("Type\tBuckets\tCard\tReadRate%\tWriteRate%\tWritePerRead%\tWritesPerBucket");
		for(int i=0; i<types.length; i++){
			final int requestedBuckets=CardinalityTracker.memToBuckets(types[i], memBytes);
			for(long target : cards){
				final CardinalityTracker ct=DDLCalibrationDriver.makeInstance(types[i], requestedBuckets, 31, 0, 0);
				final int actualBuckets=ct.actualBuckets();
				final FastRandomXoshiro rng=new FastRandomXoshiro(99);
				for(long c=0; c<target; c++){
					ct.add(rng.nextLong());
				}
				final long reads=ct.registerReads();
				final long writes=ct.registerWrites();
				final double readRate=100.0*reads/target;
				final double writeRate=100.0*writes/target;
				final double writePerRead=(reads>0) ? 100.0*writes/reads : 0;
				final double writesPerBucket=(double)writes/actualBuckets;
				System.out.printf("%s\t%d\t%d\t%.6f\t%.6f\t%.4f\t%.1f%n",
					labels[i], actualBuckets, target, readRate, writeRate, writePerRead, writesPerBucket);
			}
		}
	}

}
