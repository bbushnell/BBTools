package test.prot;

import prot.ProteinSearcher;

import java.util.Arrays;
import java.util.Random;

/** Exact arithmetic/alias regression against the pre-interleaving scalar loop. */
public final class ProteinRowDotTest {

	public static void main(final String[] args){
		check(args.length==0, "Row-dot regression accepts no arguments");
		final Random random=new Random(8537400);
		for(int rows=0; rows<=9; rows++){
			for(int dims : new int[]{0, 1, 3, 399, 400, 401}){
				for(int pad : new int[]{0, 3}){
					final int stride=Math.max(1, dims+pad);
					compare(values(random, dims), values(random, rows*stride+5), stride,
						rows, values(random, rows+3));
				}
			}
		}
		for(int pad : new int[]{0, 3}){
			compare(values(random, 400), values(random, 8537*(400+pad)+5), 400+pad,
				8537, values(random, 8540));
		}
		arithmetic(); aliases(); errors(); hotCalls();
		System.err.println("PROTEIN_ROW_DOT_TEST_PASS shapes=120 large_rows=8537 aliases=true "+
			"arithmetic=true errors=true hot_calls=12000");
	}

	private static void arithmetic(){
		final float[] ones=new float[36]; Arrays.fill(ones, 1);
		final float[] out=new float[9];
		compare(new float[]{1e20f, 1, -1e20f, 1}, ones, 4, 9, out);
		ProteinSearcher.dotRowsScalar(new float[]{1e20f, 1, -1e20f, 1}, ones, 4, 9, out);
		for(float x : out){check(x==1, "Ascending additions must retain the final unit after cancellation");}
		final float[] matrix=new float[18];
		for(int i=0; i<matrix.length; i+=2){matrix[i]=1; matrix[i+1]=0x1.fffffcp-1f;}
		final float[] q={-1, 0x1.000002p0f};
		compare(q, matrix, 2, 9, out);
		ProteinSearcher.dotRowsScalar(q, matrix, 2, 9, out);
		for(float x : out){check(Float.floatToIntBits(x)==0, "Rounded product then add must not become FMA");}
		final float[] special={0, -0f, 1, -1, Float.MIN_VALUE, -Float.MIN_VALUE,
			Float.MIN_NORMAL, Float.MAX_VALUE, -Float.MAX_VALUE,
			Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY, Float.intBitsToFloat(0x7fc12345)};
		for(float x : special){compare(new float[]{x}, special, 1, special.length, new float[special.length+2]);}
		compare(new float[0], new float[4], 1, 4, special);
		final Random random=new Random(173);
		for(int trial=0; trial<100; trial++){
			final float[] a=new float[7], b=new float[9*10+3];
			for(int i=0; i<a.length; i++){a[i]=Float.intBitsToFloat(random.nextInt());}
			for(int i=0; i<b.length; i++){b[i]=Float.intBitsToFloat(random.nextInt());}
			compare(a, b, 10, 9, new float[12]);
		}
	}

	private static void aliases(){
		final Random random=new Random(71);
		for(int rows=0; rows<=9; rows++){
			final float[] q=values(random, 9), m=values(random, rows*12+5);
			compare(q, m, 12, rows, q);
			compare(new float[]{1, 2, 3}, m, 12, rows, m);
		}
		// q==matrix permits at most one row under the stride/capacity validation.
		final float[] empty=new float[0];
		compare(empty, empty, 1, 0, empty);
		final float[] both={1, 2, 3, 4};
		compare(both, both, 4, 0, new float[3]);
		compare(both, both, 4, 1, new float[3]);
		compare(both, both, 4, 1, both);
		compare(both, both, 4, 0, both);
		final float[] q={1, 2, 3, 4}, ones=new float[16]; Arrays.fill(ones, 1);
		ProteinSearcher.dotRowsScalar(q, ones, 4, 4, q);
		equal(new float[]{10, 19, 36, 69}, q, "Earlier aliased row writes must feed later rows", false);
	}

	/** Clone the alias graph, not just each argument, so the oracle exercises the same writes. */
	private static void compare(final float[] q, final float[] matrix, final int stride,
			final int rows, final float[] out){
		final float[] eq=q.clone(), em=(matrix==q ? eq : matrix.clone());
		final float[] eo=(out==q ? eq : out==matrix ? em : out.clone());
		final float[] aq=q.clone(), am=(matrix==q ? aq : matrix.clone());
		final float[] ao=(out==q ? aq : out==matrix ? am : out.clone());
		original(eq, em, stride, rows, eo);
		ProteinSearcher.dotRowsScalar(aq, am, stride, rows, ao);
		equal(eo, ao, "Original row-dot outputs including untouched suffix", false);
		equal(eq, aq, "Query state including output aliases", false);
		equal(em, am, "Matrix state including output aliases", false);
		if(q!=out){equal(q, aq, "Non-output query must be bitwise unchanged", true);}
		if(matrix!=out){equal(matrix, am, "Non-output matrix and padding must be bitwise unchanged", true);}
		for(int i=rows; i<out.length; i++){
			check(Float.floatToRawIntBits(out[i])==Float.floatToRawIntBits(ao[i]),
				"Output suffix must not be touched at index "+i);
		}
	}

	private static void hotCalls(){
		final Random random=new Random(31);
		final float[] q=values(random, 7), matrix=values(random, 9*10);
		final float[] expected=new float[12], actual=new float[12];
		for(int trial=0; trial<12000; trial++){
			q[trial%q.length]=(trial%17-8)*0.125f;
			final int rows=trial%10;
			original(q, matrix, 10, rows, expected);
			ProteinSearcher.dotRowsScalar(q, matrix, 10, rows, actual);
			equal(expected, actual, "Repeated calls must preserve bits and untouched rows", false);
		}
	}

	private static void errors(){
		final float[] q={1, 2, 3}, matrix=new float[20], out={7, 8, 9, 10};
		reject(null, matrix, 3, 4, out); reject(q, null, 3, 4, out); reject(q, matrix, 3, 4, null);
		reject(q, matrix, 2, 4, out); reject(new float[0], matrix, 0, 4, out);
		reject(q, matrix, -1, 4, out); reject(q, matrix, 3, -1, out);
		reject(q, matrix, 3, 5, out); reject(q, matrix, 6, 4, out);
		reject(q, matrix, Integer.MAX_VALUE, 4, out);
		// Zero rows still validate stride, but require no matrix capacity.
		reject(q, matrix, 2, 0, out);
		compare(q, new float[0], Integer.MAX_VALUE, 0, out);
	}

	private static void reject(final float[] q, final float[] matrix, final int stride,
			final int rows, final float[] out){
		final float[] saved=(out==null ? null : out.clone());
		String expected=null, actual=null;
		try{original(q, matrix, stride, rows, out);}catch(IllegalArgumentException e){expected=e.getMessage();}
		try{ProteinSearcher.dotRowsScalar(q, matrix, stride, rows, out);}catch(IllegalArgumentException e){actual=e.getMessage();}
		check(expected!=null && expected.equals(actual), "Invalid arguments must retain original exception: "+expected+" / "+actual);
		if(saved!=null){equal(saved, out, "Rejected arguments must not write output", true);}
	}

	/** Frozen from b16ca995 ProteinSearcher.dotRowsScalar; independent of the production helper. */
	private static void original(final float[] q, final float[] matrix, final int stride,
			final int nRows, final float[] out){
		if(q==null || matrix==null || out==null || stride<q.length || stride<1 ||
				nRows<0 || out.length<nRows || (nRows>0 && (long)nRows*stride>matrix.length)){
			throw new IllegalArgumentException("Invalid scalar row-dot arguments: q="+
				(q==null ? -1 : q.length)+" matrix="+(matrix==null ? -1 : matrix.length)+
				" stride="+stride+" rows="+nRows+" out="+(out==null ? -1 : out.length));
		}
		for(int row=0; row<nRows; row++){
			final int base=row*stride;
			float sum=0;
			for(int d=0; d<q.length; d++){sum+=q[d]*matrix[base+d];}
			out[row]=sum;
		}
	}

	private static float[] values(final Random random, final int length){
		final float[] a=new float[length];
		for(int i=0; i<length; i++){a[i]=(random.nextFloat()-0.5f)*8;}
		return a;
	}

	/** Canonicalize NaNs only for computed values; finite values and signed zero remain exact. */
	private static void equal(final float[] expected, final float[] actual, final String why, final boolean raw){
		check(expected.length==actual.length, why+": array lengths differ");
		for(int i=0; i<expected.length; i++){
			final int e=raw ? Float.floatToRawIntBits(expected[i]) : Float.floatToIntBits(expected[i]);
			final int a=raw ? Float.floatToRawIntBits(actual[i]) : Float.floatToIntBits(actual[i]);
			if(e!=a){throw new AssertionError(why+": index="+i+" expected="+Integer.toHexString(e)+" actual="+Integer.toHexString(a));}
		}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
