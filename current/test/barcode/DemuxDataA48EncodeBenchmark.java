package test.barcode;

import java.util.Locale;

import barcode.DemuxData;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Ordinary-path timing harness for DemuxData A48 base-5 barcode encoding.
 * Generates padded legal strings from known numeric codes so the expected
 * checksum is independent of the production encoder.
 * @author Jean
 */
public final class DemuxDataA48EncodeBenchmark{

	/** Parses length=16|27 plus loop counts, verifies known codes, then times encodeACGTN. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DemuxDataA48EncodeBenchmark.class, false);
		try{
			int length=-1;
			long loops=-1, warmup=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("length")){length=Integer.parseInt(value);}
				else if(key.equals("loops")){loops=Long.parseLong(value);}
				else if(key.equals("warmup")){warmup=Long.parseLong(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check(length==16 || length==27, "Set length=16 or length=27 for the ordinary A48 paths; length="+length);
			check(loops>0, "Set loops to a positive full-pool pass count; loops="+loops);
			check(loops<=Long.MAX_VALUE/POOL, "Loop count would overflow the reported call count; loops="+loops);
			if(warmup<0){warmup=loops/16>1 ? loops/16 : 1;}
			check(warmup>=0, "Warmup loops must be nonnegative; warmup="+warmup);
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "DemuxData A48 encode benchmark requires -ea to exercise production assertions.");

			final String[] barcodes=makeBarcodes(length);
			verifyKnownCodes(barcodes, length);
			final long expected=expectedChecksum(length, loops);
			final long warmupChecksum=run(barcodes, warmup);
			blackhole=warmupChecksum;
			final long start=System.nanoTime();
			final long observed=run(barcodes, loops);
			final long elapsed=System.nanoTime()-start;
			blackhole=observed;
			check(observed==expected, "Measured checksum mismatch for known A48 codes; observed="+
					observed+", expected="+expected+", length="+length+", loops="+loops);
			final long calls=loops*POOL;
			final long step=codeStep(length);
			System.out.println("DEMUXDATA_A48_ENCODE length="+length+" loops="+loops+
					" warmup="+warmup+" pool="+POOL+" stride="+STRIDE+" offset="+OFFSET+
					" step="+step+" calls="+calls+" checksum="+observed+" expected="+expected+
					" warmup_checksum="+warmupChecksum+" elapsed_nanos="+elapsed+
					" calls_per_second="+(calls*1000000000.0/elapsed));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each fresh JVM owns one length's deterministic timing loop. */
	private DemuxDataA48EncodeBenchmark(){}

	/** Builds fixed-length A/C/G/N/T strings whose base-5 numeric values are known. */
	private static String[] makeBarcodes(final int length){
		final String[] barcodes=new String[POOL];
		for(int i=0; i<POOL; i++){
			final long code=knownCode(i, length);
			barcodes[i]=toPaddedBase5(code, length);
		}
		return barcodes;
	}

	/** Verifies the generated inputs once before the timed loop protects encodeACGTN timing. */
	private static void verifyKnownCodes(final String[] barcodes, final int length){
		check(barcodes.length==POOL, "Generated barcode pool size changed; expected="+POOL+
				", observed="+barcodes.length);
		for(int i=0; i<barcodes.length; i++){
			final String barcode=barcodes[i];
			final long expected=knownCode(i, length);
			check(barcode.length()==length, "Generated barcode length changed at index "+i+
					"; expected="+length+", observed="+barcode.length()+", barcode="+barcode);
			final long observed=DemuxData.encodeACGTN(barcode, 0, length);
			check(observed==expected, "Generated barcode code mismatch at index "+i+
					"; observed="+observed+", expected="+expected+", barcode="+barcode);
		}
	}

	/** Times full passes over the fixed barcode pool and returns a rolling checksum. */
	private static long run(final String[] barcodes, final long loops){
		long checksum=0;
		for(long loop=0; loop<loops; loop++){
			for(int i=0; i<barcodes.length; i++){
				checksum=roll(checksum, DemuxData.encodeACGTN(barcodes[i], 0, barcodes[i].length()));
			}
		}
		return checksum;
	}

	/** Computes the expected rolling checksum from the numeric codes, without calling production. */
	private static long expectedChecksum(final int length, final long loops){
		Transform cycle=Transform.identity();
		for(int i=0; i<POOL; i++){
			cycle=Transform.rollKnown(knownCode(i, length)).after(cycle);
		}
		return cycle.pow(loops).apply(0);
	}

	/** Returns one unique known base-5 code in a deterministic odd-stride order. */
	private static long knownCode(final int index, final int length){
		return permutedRank(index)*codeStep(length);
	}

	/** Returns one unique full-pool rank in a deterministic odd-stride order. */
	private static long permutedRank(final int index){
		return (((long)index)*STRIDE+OFFSET)&(POOL-1);
	}

	/** Returns a spacing that spreads the pool across the selected base-5 domain. */
	private static long codeStep(final int length){
		return pow5(length)/POOL;
	}

	/** Returns 5 to the selected A48 length. */
	private static long pow5(final int length){
		long x=1;
		for(int i=0; i<length; i++){
			check(x<=Long.MAX_VALUE/5, "Base-5 domain exceeds signed-long range; length="+length);
			x*=5;
		}
		return x;
	}

	/** Renders a known code as a left-A-padded A/C/G/N/T string. */
	private static String toPaddedBase5(final long code0, final int length){
		final char[] out=new char[length];
		long code=code0;
		for(int i=length-1; i>=0; i--){
			out[i]=BASES[(int)(code%5)];
			code/=5;
		}
		check(code==0, "Known code does not fit selected A48 length; remaining="+code+
				", original="+code0+", length="+length);
		return new String(out);
	}

	/** Updates the loop-carried checksum. */
	private static long roll(final long checksum, final long value){
		return checksum*MULT+value;
	}

	/** Always enforces benchmark oracles, independent of JVM assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Affine checksum transform x -> multiplier*x + addend, in signed-long arithmetic. */
	private static final class Transform{
		Transform(final long multiplier_, final long addend_){
			multiplier=multiplier_;
			addend=addend_;
		}

		/** Returns the identity transform. */
		static Transform identity(){return new Transform(1, 0);}

		/** Returns the transform for appending one independently known code. */
		static Transform rollKnown(final long code){return new Transform(MULT, code);}

		/** Returns this transform applied after {@code inner}. */
		Transform after(final Transform inner){
			return new Transform(multiplier*inner.multiplier, multiplier*inner.addend+addend);
		}

		/** Returns this transform repeated {@code times} times by exponentiation. */
		Transform pow(long times){
			check(times>=0, "Transform exponent must be nonnegative; times="+times);
			Transform result=identity(), base=this;
			while(times>0){
				if((times&1L)!=0){result=base.after(result);}
				base=base.after(base);
				times>>=1;
			}
			return result;
		}

		/** Applies this transform to one checksum value. */
		long apply(final long checksum){return multiplier*checksum+addend;}

		private final long multiplier;
		private final long addend;
	}

	private static final int POOL=65536;
	private static final long MULT=1000003L;
	private static final int STRIDE=40503;
	private static final int OFFSET=17;
	private static final char[] BASES={'A', 'C', 'G', 'N', 'T'};
	private static volatile long blackhole;
}
