package test.simd;

import shared.Shared;
import simd.Vector;

/** Measures one repeated byte-vector read transform with retained scalar-checked results.
 * The timed loop copies an ordinary input buffer before every transform; arraycopy is
 * part of the measured workload for both scalar and SIMD modes.
 * @author Sucrose */
public final class VectorByteTransformBenchmark{

	/** Runs warmup and a fixed amount of timed single-threaded work.
	 * @param args op=dotdash|uppercase|addandcapmin, loops=positive integer,
	 * length=positive integer, and simd=t|f
	 */
	public static void main(final String[] args){
		String op=OP_DOTDASH;
		int loops=65536, length=150;
		boolean simd=false;
		for(final String arg : args){
			if(arg.startsWith("op=")){op=arg.substring(3);}
			else if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.startsWith("length=")){length=Integer.parseInt(arg.substring(7));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!validOp(op)){throw new IllegalArgumentException("Unknown op: "+op);}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(length<1){throw new IllegalArgumentException("length must be positive: "+length);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD benchmark requires 256-bit Vector support");
		}
		final byte[][] inputs=makeInputs(op, length), outputs=makeOutputs(length);
		final boolean oldSimd=Shared.SIMD;
		try{
			Shared.SIMD=false;
			final Expected expected=makeExpected(op, inputs);
			Shared.SIMD=simd;
			final long warmup=repeat(op, inputs, outputs, expected, WARMUP_LOOPS);
			final long start=System.nanoTime();
			final long checksum=repeat(op, inputs, outputs, expected, loops);
			final long elapsed=System.nanoTime()-start;
			System.out.println("BYTE_TRANSFORM_BENCH op="+op+" simd="+simd+
					" length="+length+" loops="+loops+" buffers="+BANKS+
					" transformed_bytes="+((long)length*loops)+
					" copy_included=t wall_ns="+elapsed+" checksum="+checksum+
					" warmup_checksum="+warmup);
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Repeats one transform over copied bank entries.
	 * @param op Operation name
	 * @param inputs Immutable input templates
	 * @param outputs Mutable preallocated destination buffers
	 * @param expected Expected scalar results for each input template
	 * @param loops Number of transform calls
	 * @return Retained checksum of sampled outputs and return values
	 */
	private static long repeat(final String op, final byte[][] inputs, final byte[][] outputs,
			final Expected expected, final int loops){
		long checksum=0;
		final int length=inputs[0].length;
		for(int i=0; i<loops; i++){
			final int bank=i&(BANKS-1), sample=i%length;
			final byte[] output=outputs[bank];
			System.arraycopy(inputs[bank], 0, output, 0, length);
			final int result=apply(op, output);
			final int expectedValue=expected.outputs[bank][sample], observedValue=output[sample];
			if(result!=expected.results[bank] || observedValue!=expectedValue){
				throw new AssertionError("Mismatch op="+op+" bank="+bank+
						" sample="+sample+" expectedResult="+expected.results[bank]+
						" actualResult="+result+" expectedValue="+expectedValue+
						" actualValue="+observedValue);
			}
			checksum=mix(mix(mix(checksum, result), sample), observedValue);
		}
		return checksum;
	}

	/** Applies the selected transform.
	 * @param op Operation name
	 * @param array Mutable input/output array
	 * @return Transform return value normalized to int
	 */
	private static int apply(final String op, final byte[] array){
		if(OP_DOTDASH.equals(op)){return Vector.dotDashXToN(array) ? 1 : 0;}
		if(OP_UPPERCASE.equals(op)){return Vector.toUpperCase(array) ? 1 : 0;}
		if(OP_ADD_AND_CAP_MIN.equals(op)){
			return Vector.addAndCapMin(array, QUAL_DELTA, QUAL_FLOOR);
		}
		throw new IllegalArgumentException("Unknown op: "+op);
	}

	/** Precomputes scalar results for every input bank.
	 * @param op Operation name
	 * @param inputs Immutable input templates
	 * @return Expected outputs and return values
	 */
	private static Expected makeExpected(final String op, final byte[][] inputs){
		final byte[][] outputs=makeOutputs(inputs[0].length);
		final int[] results=new int[BANKS];
		for(int bank=0; bank<BANKS; bank++){
			System.arraycopy(inputs[bank], 0, outputs[bank], 0, inputs[bank].length);
			results[bank]=apply(op, outputs[bank]);
		}
		return new Expected(outputs, results);
	}

	/** Makes ordinary inputs for the selected operation.
	 * @param op Operation name
	 * @param length Array length
	 * @return Bank of input templates
	 */
	private static byte[][] makeInputs(final String op, final int length){
		return OP_ADD_AND_CAP_MIN.equals(op) ? makeQualities(length) : makeBases(length);
	}

	/** Makes uppercase A/C/G/T/N read bases that dotDashXToN and toUpperCase leave unchanged.
	 * @param length Array length
	 * @return Bank of base templates
	 */
	private static byte[][] makeBases(final int length){
		final byte[][] bank=new byte[BANKS][length];
		for(int b=0; b<bank.length; b++){
			final byte[] array=bank[b];
			for(int i=0; i<array.length; i++){
				array[i]=VALID_BASES[(i+b)%VALID_BASES.length];
			}
		}
		return bank;
	}

	/** Makes ordinary positive quality scores.
	 * @param length Array length
	 * @return Bank of quality templates
	 */
	private static byte[][] makeQualities(final int length){
		final byte[][] bank=new byte[BANKS][length];
		for(int b=0; b<bank.length; b++){
			final byte[] array=bank[b];
			for(int i=0; i<array.length; i++){
				array[i]=(byte)(2+((i*7+b*3)&31));
			}
		}
		return bank;
	}

	/** Makes mutable destination buffers.
	 * @param length Array length
	 * @return Bank of destination buffers
	 */
	private static byte[][] makeOutputs(final int length){
		return new byte[BANKS][length];
	}

	/** Checks whether an operation name is supported.
	 * @param op Operation name
	 * @return True if supported
	 */
	private static boolean validOp(final String op){
		return OP_DOTDASH.equals(op) || OP_UPPERCASE.equals(op) ||
				OP_ADD_AND_CAP_MIN.equals(op);
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
	private VectorByteTransformBenchmark(){}

	/** Scalar outputs and return values for each input bank. */
	private static final class Expected{
		/** Stores expected outputs and return values.
		 * @param outputs_ Expected transformed outputs
		 * @param results_ Expected return values
		 */
		private Expected(final byte[][] outputs_, final int[] results_){
			outputs=outputs_;
			results=results_;
		}

		/** Expected transformed outputs. */
		final byte[][] outputs;
		/** Expected return values. */
		final int[] results;
	}

	/** Number of immutable input buffers. */
	private static final int BANKS=128;
	/** Warmup loops run before timed work. */
	private static final int WARMUP_LOOPS=8192;
	/** Ordinary quality delta for addAndCapMin. */
	private static final byte QUAL_DELTA=5;
	/** Ordinary minimum quality for addAndCapMin. */
	private static final int QUAL_FLOOR=2;
	/** Valid read bases that exercise unchanged dot/dash and uppercase outputs. */
	private static final byte[] VALID_BASES={'A', 'C', 'G', 'T', 'N'};
	/** Operation selecting Vector.dotDashXToN. */
	private static final String OP_DOTDASH="dotdash";
	/** Operation selecting Vector.toUpperCase. */
	private static final String OP_UPPERCASE="uppercase";
	/** Operation selecting Vector.addAndCapMin. */
	private static final String OP_ADD_AND_CAP_MIN="addandcapmin";
}
