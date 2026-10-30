package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks batched profile dot products and deterministic top-N selection.
 * Requires the optional Vector module and support for the fixed256-bit backend.
 * Temporarily changes process-wide SIMD flags and restores them on every exit;
 * run in an isolated JVM, not alongside other work using those flags.
 * The selected finite inputs are correctness controls, not a throughput benchmark.
 * @contributor Shinobu */
public final class SIMDProfileDotTest{

	/** Runs normalized and unnormalized dot-product cases plus literal ranking controls.
	 * @param args Unused command-line arguments */
	public static void main(final String[] args){
		if(!Vector.simd256){
			throw new IllegalStateException("Profile SIMD checks require the Vector module and 256-bit support");
		}
		final boolean oldSimd=Shared.SIMD, oldFma=Shared.SIMD_FMA;
		int cases=0;
		try{
			Shared.SIMD_FMA=true;
			for(int dim : DIMS){
				for(int rows : ROWS){checkDotRows(dim, rows, true); cases++;}
			}
			for(int dim : new int[]{401, 1000}){
				for(int rows : ROWS){checkDotRows(dim, rows, false); cases++;}
			}
			checkTopNTies();
			System.out.println("SIMDProfileDotTest PASS cases="+(cases+1)+
				" dims="+DIMS.length+" row_shapes="+ROWS.length);
		}finally{
			Shared.SIMD=oldSimd; Shared.SIMD_FMA=oldFma;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this standalone test. */
	private SIMDProfileDotTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares scalar and vector scores within a norm-scaled error bound, then their top indices.
	 * @param dim Positive query dimension
	 * @param rows Positive matrix row count
	 * @param normalize Whether query and rows are normalized to unit length */
	private static void checkDotRows(final int dim, final int rows, final boolean normalize){
		final int stride=((dim+7)/8)*8;
		final float[] query=new float[dim];
		final float[] matrix=new float[rows*stride];
		final java.util.Random random=new java.util.Random(0xD07A5EEDL+dim*31L+rows);
		for(int i=0; i<query.length; i++){query[i]=random.nextFloat()*2f-1f;}
		final double queryNorm=norm(query);
		if(normalize){scale(query, queryNorm);}
		for(int row=0; row<rows; row++){
			final int base=row*stride;
			for(int d=0; d<dim; d++){matrix[base+d]=random.nextFloat()*2f-1f;}
			if(normalize){scale(matrix, base, dim, norm(matrix, base, dim));}
		}
		final float[] scalar=new float[rows], vector=new float[rows];
		Shared.SIMD=false;
		Vector.dotRows(query, matrix, stride, rows, scalar);
		Shared.SIMD=true;
		Vector.dotRows(query, matrix, stride, rows, vector);
		for(int row=0; row<rows; row++){
			final float delta=Math.abs(scalar[row]-vector[row]);
			final float rowNorm=(normalize ? 1f : (float)norm(matrix, row*stride, dim));
			final float bound=ABSOLUTE_TOLERANCE*(normalize ? 1f : (float)(queryNorm*rowNorm));
			if(!(delta<=bound)){
				throw new AssertionError("dotRows mismatch dim="+dim+" rows="+rows+
					" normalize="+normalize+" row="+row+" scalar="+scalar[row]+" vector="+vector[row]+" delta="+delta+" bound="+bound);
			}
		}
		final int k=Math.min(100, rows);
		final int[] scalarTop=new int[k], vectorTop=new int[k];
		Vector.topN(scalar, k, scalarTop);
		Vector.topN(vector, k, vectorTop);
		if(!Arrays.equals(scalarTop, vectorTop)){
			throw new AssertionError("topN mismatch dim="+dim+" rows="+rows+
				" scalar="+Arrays.toString(scalarTop)+" vector="+Arrays.toString(vectorTop));
		}
	}

	/** Computes the Euclidean length with double products and accumulation.
	 * @param values Finite input vector
	 * @return Its Euclidean length */
	private static double norm(final float[] values){
		double sum=0;
		for(float value : values){sum+=(double)value*value;}
		return Math.sqrt(sum);
	}

	/** Computes the Euclidean length of one matrix row.
	 * @param values Matrix storage
	 * @param from Inclusive row start
	 * @param length Number of row elements
	 * @return Euclidean row length */
	private static double norm(final float[] values, final int from, final int length){
		double sum=0;
		for(int i=from; i<from+length; i++){sum+=(double)values[i]*values[i];}
		return Math.sqrt(sum);
	}

	/** Normalizes a vector using its previously computed positive length.
	 * @param values Vector to modify
	 * @param length Finite positive Euclidean length */
	private static void scale(final float[] values, final double length){
		final float inverse=(float)(1.0/length);
		for(int i=0; i<values.length; i++){values[i]*=inverse;}
	}

	/** Normalizes one row while preserving its padding and neighboring rows.
	 * @param values Matrix storage to modify
	 * @param from Inclusive row start
	 * @param length Number of row elements
	 * @param norm Finite positive row length */
	private static void scale(final float[] values, final int from, final int length, final double norm){
		final float inverse=(float)(1.0/norm);
		for(int i=from; i<from+length; i++){values[i]*=inverse;}
	}

	/** Checks descending score order, ascending index ties, NaN exclusion and zero selection. */
	private static void checkTopNTies(){
		final float[] scores={1f, 2f, 2f, Float.NaN, 2f, -1f};
		final int[] top=new int[5];
		Vector.topN(scores, top.length, top);
		final int[] expected={1, 2, 4, 0, 5};
		if(!Arrays.equals(top, expected)){
			throw new AssertionError("topN tie/NaN order mismatch: "+Arrays.toString(top));
		}
		Vector.topN(scores, 0, new int[0]);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Dimensions spanning vector widths, padding and scalar tails. */
	private static final int[] DIMS={14, 20, 196, 200, 400, 401, 1000};
	/** Matrix shapes covering one row, a partial tile and many complete tiles. */
	private static final int[] ROWS={1, 7, 4432};
	/** Absolute score error allowance, multiplied by the two vector norms. */
	private static final float ABSOLUTE_TOLERANCE=2e-5f;
}
