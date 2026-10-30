package test.simd;

import shared.Shared;
import simd.Vector;

/** Measures byte predicate throughput with retained, per-call result checks.
 * Input buffers are preallocated and mutated in-place with valid read bases so
 * repeated calls cannot read a compile-time-stable array.
 * @author Sucrose */
public final class VectorBytePredicateBenchmark{

	/** Runs warmup and a fixed amount of timed single-threaded work.
	 * @param args predicate=protein|nucleotide, variant=valid|late, loops=positive
	 * integer, length=integer at least 2, and simd=t|f
	 */
	public static void main(final String[] args){
		String predicate=PREDICATE_NUCLEOTIDE, variant=VARIANT_VALID;
		int loops=65536, length=150;
		boolean simd=false;
		for(final String arg : args){
			if(arg.startsWith("predicate=")){predicate=arg.substring(10);}
			else if(arg.startsWith("variant=")){variant=arg.substring(8);}
			else if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.startsWith("length=")){length=Integer.parseInt(arg.substring(7));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!validPredicate(predicate)){
			throw new IllegalArgumentException("Unknown predicate: "+predicate);
		}
		if(!validVariant(variant)){throw new IllegalArgumentException("Unknown variant: "+variant);}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(length<2){throw new IllegalArgumentException("length must be at least 2: "+length);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD benchmark requires 256-bit Vector support");
		}
		final byte[][] reads=makeReads(length, VARIANT_LATE.equals(variant));
		final boolean oldSimd=Shared.SIMD;
		try{
			Shared.SIMD=simd;
			final boolean expected=expected(predicate, variant);
			final long warmup=repeat(predicate, reads, expected, WARMUP_LOOPS);
			final long start=System.nanoTime();
			final long checksum=repeat(predicate, reads, expected, loops);
			final long elapsed=System.nanoTime()-start;
			System.out.println("BYTE_PREDICATE_BENCH predicate="+predicate+
					" variant="+variant+" simd="+simd+" length="+length+
					" loops="+loops+" buffers="+BANKS+" predicate_calls="+loops+
					" tested_bytes="+((long)length*loops)+" wall_ns="+elapsed+
					" checksum="+checksum+" warmup_checksum="+warmup);
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Repeats one predicate over preallocated mutable read buffers.
	 * @param predicate Predicate name
	 * @param reads Mutable read buffers
	 * @param expected Expected predicate result
	 * @param loops Number of predicate calls
	 * @return Retained checksum of call results and sampled bases
	 */
	private static long repeat(final String predicate, final byte[][] reads,
			final boolean expected, final int loops){
		long checksum=0;
		final int length=reads[0].length;
		final int mutableLength=length-(reads[0][length-1]==LATE_AMINO ? 1 : 0);
		final int[] cursors=new int[BANKS];
		for(int i=0; i<loops; i++){
			final int bank=i&(BANKS-1);
			final byte[] read=reads[bank];
			final int pos=cursors[bank]++;
			if(cursors[bank]>=mutableLength){cursors[bank]=0;}
			read[pos]=(read[pos]=='A' ? (byte)'C' : (byte)'A');
			final boolean result=apply(predicate, read);
			if(result!=expected){
				throw new AssertionError("Mismatch predicate="+predicate+" bank="+bank+
						" pos="+pos+" expected="+expected+" actual="+result);
			}
			checksum=mix(mix(checksum, result ? 1 : 0), read[pos]);
		}
		return checksum;
	}

	/** Applies the selected predicate.
	 * @param predicate Predicate name
	 * @param read Read buffer
	 * @return Predicate result
	 */
	private static boolean apply(final String predicate, final byte[] read){
		if(PREDICATE_NUCLEOTIDE.equals(predicate)){return Vector.isNucleotide(read);}
		if(PREDICATE_PROTEIN.equals(predicate)){return Vector.isProtein(read);}
		throw new IllegalArgumentException("Unknown predicate: "+predicate);
	}

	/** Makes realistic valid read bases with an optional final amino-only symbol.
	 * @param length Read length
	 * @param lateAmino Whether to place a late invalid/amino symbol
	 * @return Bank of mutable reads
	 */
	private static byte[][] makeReads(final int length, final boolean lateAmino){
		final byte[][] bank=new byte[BANKS][length];
		for(int b=0; b<bank.length; b++){
			final byte[] read=bank[b];
			for(int i=0; i<read.length; i++){read[i]=VALID_BASES[(i+b)&VALID_BASE_MASK];}
			if(lateAmino){read[read.length-1]=LATE_AMINO;}
		}
		return bank;
	}

	/** Returns the expected result for a predicate and variant.
	 * @param predicate Predicate name
	 * @param variant Variant name
	 * @return Expected predicate result
	 */
	private static boolean expected(final String predicate, final String variant){
		final boolean late=VARIANT_LATE.equals(variant);
		if(PREDICATE_NUCLEOTIDE.equals(predicate)){return !late;}
		if(PREDICATE_PROTEIN.equals(predicate)){return late;}
		throw new IllegalArgumentException("Unknown predicate: "+predicate);
	}

	/** Checks whether a predicate name is supported.
	 * @param predicate Predicate name
	 * @return True if supported
	 */
	private static boolean validPredicate(final String predicate){
		return PREDICATE_NUCLEOTIDE.equals(predicate) || PREDICATE_PROTEIN.equals(predicate);
	}

	/** Checks whether a variant name is supported.
	 * @param variant Variant name
	 * @return True if supported
	 */
	private static boolean validVariant(final String variant){
		return VARIANT_VALID.equals(variant) || VARIANT_LATE.equals(variant);
	}

	/** Mixes an int value into the checksum.
	 * @param checksum Previous checksum
	 * @param value Retained value
	 * @return Updated checksum
	 */
	private static long mix(long checksum, final int value){
		checksum^=(value&0xFFL)+0x9E3779B97F4A7C15L;
		checksum*=0xBF58476D1CE4E5B9L;
		return Long.rotateLeft(checksum, 11);
	}

	/** Prevents construction of this standalone benchmark. */
	private VectorBytePredicateBenchmark(){}

	/** Number of mutable read buffers. */
	private static final int BANKS=128;
	/** Warmup loops run before timed work. */
	private static final int WARMUP_LOOPS=8192;
	/** Realistic read bases. */
	private static final byte[] VALID_BASES={'A', 'C', 'G', 'T'};
	/** Mask for VALID_BASES lookup. */
	private static final int VALID_BASE_MASK=VALID_BASES.length-1;
	/** Late amino-only and non-nucleotide symbol recognized by old and expanded masks. */
	private static final byte LATE_AMINO='E';
	/** Operation selecting Vector.isNucleotide. */
	private static final String PREDICATE_NUCLEOTIDE="nucleotide";
	/** Operation selecting Vector.isProtein. */
	private static final String PREDICATE_PROTEIN="protein";
	/** All bases are valid nucleotide bases. */
	private static final String VARIANT_VALID="valid";
	/** The last base is an amino-only letter and invalid nucleotide. */
	private static final String VARIANT_LATE="late";
}
