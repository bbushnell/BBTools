package test.simd;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;
import simd.Vector;

/** Exercises legacy literal-I sparse CellNet rows whose edge count is full-width,
 * and malformed literal sparse rows the parser must reject.
 * @author Sucrose */
public final class VectorSparseCellNetTest{

	/** Loads the identity sparse model, then checks malformed literal sparse models are rejected.
	 * @param args Optional simd=t or simd=f selector */
	public static void main(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd128){
			throw new IllegalStateException("Sparse CellNet checks require 128/256-bit Vector support");
		}

		final boolean oldSimd=Shared.SIMD, oldSimdFma=Shared.SIMD_FMA;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			Shared.SIMD_FMA=true;
			failures+=checkModel("identity", IDENTITY_MODEL, EXPECTED); cases++;
			failures+=checkRejected("reordered", REORDERED_MODEL); cases++;
			failures+=checkRejected("duplicate", DUPLICATE_MODEL); cases++;
			failures+=checkRejected("out-of-range", OUT_OF_RANGE_MODEL); cases++;
			failures+=checkRejected("negative", NEGATIVE_MODEL); cases++;
			System.out.println("SPARSE_CELLNET_RESULT simd="+simd+" cases="+cases+
					" failures="+failures);
			if(failures>0){throw new AssertionError("Sparse CellNet failures="+failures);}
		}finally{
			Shared.SIMD_FMA=oldSimdFma;
			Shared.SIMD=oldSimd;
		}
	}

	/** Loads a literal sparse model and compares its first output to an independent value.
	 * @param label Diagnostic model label
	 * @param rows Literal BBNet rows
	 * @param expected Expected first output
	 * @return One on mismatch, otherwise zero */
	private static int checkModel(final String label, final String[] rows, final float expected){
		assertModelWidth(label, rows);
		final CellNet net=CellNetParser.loadInferenceFromLines(lines(rows));
		final float actual=net.applyInput(INPUT).feedForward();
		if(Float.floatToIntBits(actual)==Float.floatToIntBits(expected)){return 0;}
		System.out.println("FAIL CellNet "+label+" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Checks that a malformed literal sparse model crashes loudly during parsing.
	 * @param label Diagnostic model label
	 * @param rows Literal BBNet rows
	 * @return One on failed rejection, otherwise zero */
	private static int checkRejected(final String label, final String[] rows){
		assertModelWidth(label, rows);
		try{
			CellNetParser.loadInferenceFromLines(lines(rows));
			System.out.println("FAIL CellNet "+label+" accepted malformed literal inputs");
			return 1;
		}catch(IllegalArgumentException expected){
			final String message=expected.getMessage();
			if(message!=null && message.startsWith(REJECTION_PREFIX)){return 0;}
			System.out.println("FAIL CellNet "+label+" wrong rejection message: "+message);
			return 1;
		}
	}

	/** Converts ASCII fixture rows into the production parser's borrowed row format.
	 * @param rows Text rows
	 * @return New byte-row list */
	private static ArrayList<byte[]> lines(final String[] rows){
		final ArrayList<byte[]> list=new ArrayList<byte[]>(rows.length);
		for(final String row : rows){list.add(row.getBytes(StandardCharsets.US_ASCII));}
		return list;
	}

	/** Checks that a fixture row has sixteen literal inputs and sixteen weights.
	 * @param label Diagnostic model label
	 * @param rows Literal BBNet rows */
	private static void assertModelWidth(final String label, final String[] rows){
		int inputs=-1, weights=-1;
		for(final String row : rows){
			if(row.startsWith("I"+CELL_ID+" ")){inputs=countFields(row)-1;}
			else if(row.startsWith("W"+CELL_ID+" ")){weights=countFields(row)-3;}
		}
		if(inputs!=WIDTH || weights!=WIDTH){
			throw new AssertionError(label+" literal widths: inputs="+inputs+", weights="+weights);
		}
	}

	/** Counts single-space-delimited fields in a literal fixture row.
	 * @param row Literal row
	 * @return Number of fields */
	private static int countFields(final String row){
		int count=1;
		for(int i=0; i<row.length(); i++){
			if(row.charAt(i)==' '){count++;}
		}
		return count;
	}

	/** Prevents construction of this standalone test. */
	private VectorSparseCellNetTest(){}

	/** Input width that is the smallest width to trigger the equal-length sparse shortcut. */
	private static final int WIDTH=16;

	/** Global output-cell ID in a 16-input, one-output network. */
	private static final int CELL_ID=17;

	/** Input row whose first two values expose whether the literal I row is honored. */
	private static final float[] INPUT={2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/** Independently computed prediction: bias 3 plus 100*2 plus 10*5. */
	private static final float EXPECTED=253;

	/** Parser diagnostic expected for unordered, duplicate and out-of-range literal indices. */
	private static final String REJECTION_PREFIX=
		"Legacy sparse inputs must be strictly increasing within previous layer:";

	/** Legacy literal sparse model whose row order matches weight order. */
	private static final String[] IDENTITY_MODEL={
		"##bbnet",
		"#version 1",
		"#concise",
		"#sparse",
		"#density 1",
		"#blocksize 1",
		"#seed 1",
		"#layers 2",
		"#dims 16 1",
		"#coding decimal",
		"I17 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15",
		"W17 LINEAR 3 100 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
	};

	/** Legacy literal sparse model whose first two input indices are reordered. */
	private static final String[] REORDERED_MODEL={
		"##bbnet",
		"#version 1",
		"#concise",
		"#sparse",
		"#density 1",
		"#blocksize 1",
		"#seed 1",
		"#layers 2",
		"#dims 16 1",
		"#coding decimal",
		"I17 1 0 2 3 4 5 6 7 8 9 10 11 12 13 14 15",
		"W17 LINEAR 3 10 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
	};

	/** Literal sparse model with duplicate input 1 and missing input 2. */
	private static final String[] DUPLICATE_MODEL={
		"##bbnet",
		"#version 1",
		"#concise",
		"#sparse",
		"#density 1",
		"#blocksize 1",
		"#seed 1",
		"#layers 2",
		"#dims 16 1",
		"#coding decimal",
		"I17 0 1 1 3 4 5 6 7 8 9 10 11 12 13 14 15",
		"W17 LINEAR 3 100 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
	};

	/** Literal sparse model with input 16 outside the 0..15 input range. */
	private static final String[] OUT_OF_RANGE_MODEL={
		"##bbnet",
		"#version 1",
		"#concise",
		"#sparse",
		"#density 1",
		"#blocksize 1",
		"#seed 1",
		"#layers 2",
		"#dims 16 1",
		"#coding decimal",
		"I17 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16",
		"W17 LINEAR 3 100 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
	};

	/** Literal sparse model with negative input -1. */
	private static final String[] NEGATIVE_MODEL={
		"##bbnet",
		"#version 1",
		"#concise",
		"#sparse",
		"#density 1",
		"#blocksize 1",
		"#seed 1",
		"#layers 2",
		"#dims 16 1",
		"#coding decimal",
		"I17 -1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15",
		"W17 LINEAR 3 100 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
	};
}
