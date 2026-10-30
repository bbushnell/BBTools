package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks primitive reductions and arithmetic at fixed-width and tail boundaries.
 * Exact small binary fractions avoid imposing an unsupported floating reduction
 * order; long expectations retain Java overflow semantics.
 * @author Shinobu */
public final class VectorPrimitiveWidthTest{

	/** Runs finite primitive controls in one isolated JVM.
	 * @param args Optional simd=t or simd=f */
	public static void main(final String[] args){
		boolean simd=true;
		for(String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd128){throw new IllegalStateException("128-bit Vector support is required");}
		final boolean oldSimd=Shared.SIMD, oldFma=Shared.SIMD_FMA;
		final boolean oldCopy=Vector.SIMDCOPY, oldSparse=Vector.SIMD_MULT_SPARSE;
		try{
			Shared.SIMD=simd;
			Shared.SIMD_FMA=Vector.SIMDCOPY=Vector.SIMD_MULT_SPARSE=true;
			for(int length : LENGTHS){
				for(int pattern=0; pattern<3; pattern++){checkReductions(length, pattern);}
				checkArithmetic(length);
			}
			checkSparseBlocks();
			System.out.println("PRIMITIVE_WIDTH_RESULT simd="+simd+" simd_bits="+
					Vector.SIMD_BITS+" checks="+checks+" failures=0");
		}finally{
			Shared.SIMD=oldSimd;
			Shared.SIMD_FMA=oldFma;
			Vector.SIMDCOPY=oldCopy;
			Vector.SIMD_MULT_SPARSE=oldSparse;
		}
	}

	/** Checks widening, wrapping, prefixes, and exact binary-fraction reductions.
	 * @param length Array length
	 * @param pattern Mixed, all-positive, or all-negative byte pattern */
	private static void checkReductions(final int length, final int pattern){
		final byte[] bytes=new byte[length];
		final int[] ints=new int[length];
		final long[] longs=new long[length];
		final double[] doubles=new double[length];
		long byteSum=0, intSum=0, longSum=0, longMax=Long.MIN_VALUE;
		double doubleSum=0;
		for(int i=0; i<length; i++){
			bytes[i]=(byte)(pattern==1 ? 127 : pattern==2 ? -128 : i*37-121);
			ints[i]=(i%3==0 ? Integer.MAX_VALUE : i%3==1 ? Integer.MIN_VALUE : i);
			longs[i]=(i%3==0 ? Long.MAX_VALUE : i%3==1 ? Long.MIN_VALUE : i);
			doubles[i]=(i%17-8)*0.5;
			byteSum+=bytes[i]; intSum+=ints[i]; longSum+=longs[i]; doubleSum+=doubles[i];
			longMax=Math.max(longMax, longs[i]);
		}
		require(Vector.sum(bytes)==byteSum, "byte sum length="+length+" pattern="+pattern);
		require(Vector.sum(ints)==intSum, "int sum length="+length);
		require(Vector.sum(longs)==longSum, "long sum length="+length);
		require(Vector.sum(doubles)==doubleSum, "double sum length="+length);
		if(length>0){require(Vector.max(longs)==longMax, "long max length="+length);}
		final int from=Math.min(3, length), to=length-1;
		long intRange=0, longRange=0;
		for(int i=from; i<=to; i++){intRange+=ints[i]; longRange+=longs[i];}
		require(Vector.sum(ints, from, to)==intRange, "int range length="+length);
		require(Vector.sum(longs, from, to)==longRange, "long range length="+length);
	}

	/** Checks dense arithmetic and copying with an untouched destination suffix.
	 * @param length Array length */
	private static void checkArithmetic(final int length){
		final float[] a=new float[length], b=new float[length];
		float expectedDot=0;
		for(int i=0; i<length; i++){
			a[i]=(i%7-3)*0.5f;
			b[i]=(i%5-2)*0.25f;
			expectedDot+=a[i]*b[i];
		}
		require(Vector.fma(a, b)==expectedDot, "dense dot length="+length);
		float[] actual=a.clone(), expected=a.clone();
		for(int i=0; i<length; i++){expected[i]+=b[i];}
		Vector.add(actual, b);
		require(Arrays.equals(actual, expected), "dense add length="+length);
		actual=a.clone(); expected=a.clone();
		for(int i=0; i<length; i++){expected[i]+=b[i]*0.5f;}
		Vector.addProduct(actual, b, 0.5f);
		require(Arrays.equals(actual, expected), "dense addProduct length="+length);
		actual=new float[length+5]; expected=new float[length+5];
		Arrays.fill(actual, -77); Arrays.fill(expected, -77);
		System.arraycopy(b, 0, expected, 0, b.length);
		Vector.copy(actual, b);
		require(Arrays.equals(actual, expected), "copy suffix length="+length);
	}

	/** Tests four-entry contiguous blocks that are not eight-entry contiguous blocks. */
	private static void checkSparseBlocks(){
		final int[] starts={8, 0, 20, 12};
		final int[] indexes=new int[16];
		final float[] a=new float[16], b=new float[32];
		for(int i=0; i<b.length; i++){b[i]=(i+1)*0.25f;}
		float dot=0;
		final float[] expected=new float[a.length];
		for(int i=0; i<a.length; i++){
			indexes[i]=starts[i/4]+i%4;
			a[i]=(i%3+1)*0.5f;
			dot+=a[i]*b[indexes[i]];
			expected[i]=a[i]+b[indexes[i]]*0.5f;
		}
		require(Vector.fma(a, b, indexes, 4, true)==dot, "four-entry sparse dot");
		Vector.addProduct(a, b, indexes, 0.5f, 4);
		require(Arrays.equals(a, expected), "four-entry sparse addProduct");
	}

	/** Records a checked expectation or fails loudly.
	 * @param condition Expected result
	 * @param message Case identifier */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Prevents construction of this standalone fixture. */
	private VectorPrimitiveWidthTest(){}

	/** Counts assertions completed by this isolated process. */
	private static int checks=0;
	/** Boundaries for two-, four-, eight-, sixteen- and thirty-two-lane operations. */
	private static final int[] LENGTHS={0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31, 32, 33, 65};
}
