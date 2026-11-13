package test.barcode;

import java.util.Locale;

import barcode.Barcode;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Measures left-distance calls over 1024 deterministic, valid string pairs.
 * Families i%4 are plain8, plain16, dual8+8 and dual16+16; only a known prefix
 * of the left component differs. Pair validation, construction and output are
 * outside timing. This is a helper benchmark, not file or demultiplexer throughput.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeLeftDistanceBenchmark{

	/** Parses loops and optional warmup; default warmup is max(64, loops/8). */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BarcodeLeftDistanceBenchmark.class, false);
		try{
			int loops=0, warmup=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("loops")){loops=Integer.parseInt(value);}
				else if(key.equals("warmup")){warmup=Integer.parseInt(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check(loops>0 && loops<=MAX_LOOPS, "Use 1.."+MAX_LOOPS+" measured passes.");
			if(warmup<0){warmup=Math.max(64, loops/8);}
			check(warmup>0 && warmup<=MAX_LOOPS, "Use 1.."+MAX_LOOPS+" warmup passes.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Left-distance benchmark requires -ea for production preconditions.");
			Shared.setThreads(1);
			final Inputs inputs=new Inputs();
			inputs.verifyPairs();
			final long warmupChecksum=inputs.run(warmup);
			check(warmupChecksum==6108L*warmup, "Warmup checksum differs from the independent literal.");
			final long start=System.nanoTime();
			final long checksum=inputs.run(loops);
			final long elapsed=System.nanoTime()-start;
			check(checksum==6108L*loops && elapsed>0, "Timed checksum differs or clock did not advance: "+checksum+", "+elapsed);
			blackhole=checksum;
			System.out.println("BARCODE_LEFT_DISTANCE_BENCHMARK loops="+loops+" warmup="+warmup+
					" pairs="+PAIRS+" checksum="+checksum+" expected="+(6108L*loops)+" elapsed_nanos="+elapsed);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Owns immutable strings and reused arrays; each family has 256 pairs. */
	private static final class Inputs{

		/** Base-4 digits of i vary references; rotating a prefix base guarantees each intended mismatch. */
		Inputs(){
			final int[] sums=new int[4];
			for(int i=0; i<PAIRS; i++){
				final int family=i%4, leftLength=(family&1)==0 ? 8 : 16;
				final boolean dual=family>=2;
				final int mutations=(i/4)%(leftLength+1);
				final int length=dual ? leftLength*2+1 : leftLength;
				final char[] query=new char[length], reference=new char[length];
				for(int pos=0; pos<length; pos++){
					if(dual && pos==leftLength){query[pos]=reference[pos]='+';}
					else{
						final int base=((i>>>(2*(pos%5)))+pos)&3;
						reference[pos]=BASES[base];
						query[pos]=BASES[(base+(pos<mutations ? 1 : 0))&3];
					}
				}
				queries[i]=new String(query);
				references[i]=new String(reference);
				expected[i]=mutations;
				sums[family]+=mutations;
			}
			check(sums[0]==1014 && sums[1]==2040 && sums[2]==1014 && sums[3]==2040 &&
					sums[0]+sums[1]+sums[2]+sums[3]==6108,
					"Independent cohort sums must be 1014/2040/1014/2040, totaling 6108.");
		}

		/** Checks each literal mutation count before warmup, without deriving an oracle from production results. */
		void verifyPairs(){
			for(int i=0; i<PAIRS; i++){
				final int actual=Barcode.hdistL(queries[i], references[i]);
				check(actual==expected[i], "Pair "+i+" expected distance="+expected[i]+", observed="+actual);
			}
		}

		/** No timed allocation or output; every pass consumes all pairs into the returned long checksum. */
		long run(final int loops){
			assert(loops>0) : "A benchmark checksum must include at least one complete input pass.";
			long checksum=0;
			for(int pass=0; pass<loops; pass++){
				for(int i=0; i<PAIRS; i++){checksum+=Barcode.hdistL(queries[i], references[i]);}
			}
			return checksum;
		}

		private final String[] queries=new String[PAIRS], references=new String[PAIRS];
		private final int[] expected=new int[PAIRS];
	}

	/** Keeps independent checksum and pair-oracle failures visible under any assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private BarcodeLeftDistanceBenchmark(){}
	private static final int PAIRS=1024, MAX_LOOPS=10000000;
	private static final char[] BASES={'A', 'C', 'G', 'T'};
	private static volatile long blackhole;
}
