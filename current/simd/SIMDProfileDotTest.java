package simd;

import java.util.Arrays;

import shared.Shared;

/** Correctness selftest for batched profile dot products and deterministic top-N selection. */
public final class SIMDProfileDotTest {

	private static final int[] DIMS={14,20,196,200,400,401,1000};
	private static final int[] ROWS={1,7,4432};
	private static final float ABSOLUTE_TOLERANCE=2e-5f;

	private SIMDProfileDotTest(){}

	public static void main(final String[] args){
		final boolean oldSimd=Shared.SIMD, oldFma=Shared.SIMD_FMA;
		int cases=0;
		try{
			Shared.SIMD_FMA=true;
			for(int dim : DIMS){
				for(int rows : ROWS){checkDotRows(dim, rows, true); cases++;}
			}
			for(int dim : new int[]{401,1000}){
				for(int rows : ROWS){checkDotRows(dim, rows, false); cases++;}
			}
			checkTopNTies();
			System.out.println("SIMDProfileDotTest PASS cases="+(cases+1)+
				" dims="+DIMS.length+" row_shapes="+ROWS.length);
		}finally{
			Shared.SIMD=oldSimd; Shared.SIMD_FMA=oldFma;
		}
	}

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

	private static double norm(final float[] values){
		double sum=0;
		for(float value : values){sum+=(double)value*value;}
		return Math.sqrt(sum);
	}

	private static double norm(final float[] values, final int from, final int length){
		double sum=0;
		for(int i=from; i<from+length; i++){sum+=(double)values[i]*values[i];}
		return Math.sqrt(sum);
	}

	private static void scale(final float[] values, final double length){
		final float inverse=(float)(1.0/length);
		for(int i=0; i<values.length; i++){values[i]*=inverse;}
	}

	private static void scale(final float[] values, final int from, final int length, final double norm){
		final float inverse=(float)(1.0/norm);
		for(int i=from; i<from+length; i++){values[i]*=inverse;}
	}

	private static void checkTopNTies(){
		final float[] scores={1f,2f,2f,Float.NaN,2f,-1f};
		final int[] top=new int[5];
		Vector.topN(scores, top.length, top);
		final int[] expected={1,2,4,0,5};
		if(!Arrays.equals(top, expected)){
			throw new AssertionError("topN tie/NaN order mismatch: "+Arrays.toString(top));
		}
		Vector.topN(scores, 0, new int[0]);
	}
}
