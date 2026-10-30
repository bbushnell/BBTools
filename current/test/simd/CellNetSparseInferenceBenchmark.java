package test.simd;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;
import simd.Vector;

/** Measures sparse CellNet inference on generated linear one-output networks.
 * The ordered mode has one sorted full-width sparse row, exercising the dense
 * fast path for identity sparse inputs at blockSize 1.  The pruned mode keeps
 * every second input, exercising indexed sparse FMA with the same parser and
 * feed-forward call shape.
 * @author Sucrose */
public final class CellNetSparseInferenceBenchmark{

	/** Runs warmup and timed sparse CellNet inference over mutable input banks.
	 * @param args mode=ordered|pruned, loops=positive integer,
	 * width=integer from 16 to 32768, and simd=t|f
	 */
	public static void main(final String[] args){
		String mode=MODE_ORDERED;
		int loops=65536, width=256;
		boolean simd=false;
		for(final String arg : args){
			if(arg.startsWith("mode=")){mode=arg.substring(5);}
			else if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.startsWith("width=")){width=Integer.parseInt(arg.substring(6));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!validMode(mode)){throw new IllegalArgumentException("Unknown mode: "+mode);}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(width<MIN_WIDTH){throw new IllegalArgumentException("width must be at least 16: "+width);}
		if(width>MAX_WIDTH){throw new IllegalArgumentException(
				"width must be at most 32768 for an exact float oracle: "+width);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("Sparse CellNet benchmark requires 256-bit Vector support");
		}

		final Workload workload=makeWorkload(mode, width);
		final CellNet net=CellNetParser.loadInferenceFromLines(makeModel(mode, width));
		final boolean oldSimd=Shared.SIMD, oldSimdFma=Shared.SIMD_FMA;
		try{
			Shared.SIMD=simd;
			Shared.SIMD_FMA=true;
			final long warmup=repeat(net, workload, WARMUP_LOOPS);
			final long start=System.nanoTime();
			final long checksum=repeat(net, workload, loops);
			final long elapsed=System.nanoTime()-start;
			System.out.println("CELLNET_SPARSE_INFERENCE_BENCH mode="+mode+
					" simd="+simd+" width="+width+" edges="+workload.indices.length+
					" loops="+loops+" buffers="+BANKS+" blocksize=1"+
					" fma_elements="+((long)workload.indices.length*loops)+
					" wall_ns="+elapsed+" checksum="+checksum+
					" warmup_checksum="+warmup);
		}finally{
			Shared.SIMD_FMA=oldSimdFma;
			Shared.SIMD=oldSimd;
		}
	}

	/** Repeats sparse inference without allocating in the timed path.
	 * @param net Parsed sparse inference network
	 * @param workload Mutable inputs and exact expected outputs
	 * @param loops Number of predictions
	 * @return Retained checksum of exact outputs
	 */
	private static long repeat(final CellNet net, final Workload workload, final int loops){
		long checksum=0;
		for(int i=0; i<loops; i++){
			final int bank=i&(BANKS-1);
			final float expected=mutate(workload, bank);
			final float actual=net.applyInput(workload.inputs[bank]).feedForward();
			if(Float.floatToIntBits(actual)!=Float.floatToIntBits(expected)){
				throw new AssertionError("feedForward bank="+bank+" expected="+
						expected+" actual="+actual);
			}
			checksum=mix(mix(checksum, Float.floatToIntBits(actual)), bank);
		}
		return checksum;
	}

	/** Mutates one active input and updates that bank's independent exact output.
	 * @param workload Mutable benchmark state
	 * @param bank Input bank to mutate
	 * @return Expected output after the mutation
	 */
	private static float mutate(final Workload workload, final int bank){
		final int edge=workload.cursors[bank];
		int nextEdge=edge+1;
		if(nextEdge>=workload.indices.length){nextEdge=0;}
		workload.cursors[bank]=nextEdge;

		final int input=workload.indices[edge];
		final float[] values=workload.inputs[bank];
		final int old=(int)values[input];
		final int value=(old>=MAX_INPUT ? MIN_INPUT : old+1);
		values[input]=value;
		workload.expected[bank]+=(value-old)*workload.weights[edge];
		return workload.expected[bank];
	}

	/** Builds mutable inputs and independent exact expected outputs.
	 * @param mode Benchmark mode
	 * @param width Input width
	 * @return Mutable benchmark state
	 */
	private static Workload makeWorkload(final String mode, final int width){
		final int[] indices=makeIndices(mode, width);
		final float[] weights=makeWeights(indices.length);
		final Workload workload=new Workload(width, indices, weights);
		for(int bank=0; bank<BANKS; bank++){
			final float[] values=workload.inputs[bank];
			for(int i=0; i<values.length; i++){
				values[i]=((i*17+bank*5)&63)-31;
			}
			workload.expected[bank]=expected(values, indices, weights);
			workload.cursors[bank]=bank%indices.length;
		}
		return workload;
	}

	/** Independently evaluates the generated linear sparse row.
	 * @param values Input values
	 * @param indices Sparse input indices
	 * @param weights Sparse input weights
	 * @return Exact generated model output
	 */
	private static float expected(final float[] values, final int[] indices,
			final float[] weights){
		float sum=BIAS;
		for(int i=0; i<indices.length; i++){sum+=values[indices[i]]*weights[i];}
		return sum;
	}

	/** Builds a sparse BBNet model for the selected row shape.
	 * @param mode Benchmark mode
	 * @param width Input width
	 * @return ASCII BBNet rows for CellNetParser
	 */
	private static ArrayList<byte[]> makeModel(final String mode, final int width){
		final int[] indices=makeIndices(mode, width);
		final float[] weights=makeWeights(indices.length);
		final ArrayList<byte[]> rows=new ArrayList<byte[]>(MODEL_HEADER_ROWS+2);
		add(rows, "##bbnet");
		add(rows, "#version 1");
		add(rows, "#concise");
		add(rows, "#sparse");
		add(rows, "#density 1");
		add(rows, "#blocksize 1");
		add(rows, "#seed 1");
		add(rows, "#layers 2");
		add(rows, "#dims "+width+" 1");
		add(rows, "#coding decimal");
		final StringBuilder inputs=new StringBuilder(width*4+8);
		inputs.append('I').append(width+1);
		for(int i=0; i<indices.length; i++){inputs.append(' ').append(indices[i]);}
		add(rows, inputs.toString());

		final StringBuilder weightRow=new StringBuilder(width*6+16);
		weightRow.append('W').append(width+1).append(" LINEAR ").append(BIAS);
		for(int i=0; i<weights.length; i++){weightRow.append(' ').append(weights[i]);}
		add(rows, weightRow.toString());
		return rows;
	}

	/** Appends one ASCII parser row.
	 * @param rows Destination rows
	 * @param row Text row
	 */
	private static void add(final ArrayList<byte[]> rows, final String row){
		rows.add(row.getBytes(StandardCharsets.US_ASCII));
	}

	/** Creates the sparse input indices for one benchmark mode.
	 * @param mode Benchmark mode
	 * @param width Input width
	 * @return Sorted sparse input indices
	 */
	private static int[] makeIndices(final String mode, final int width){
		final int edges=MODE_ORDERED.equals(mode) ? width : width/2;
		final int[] indices=new int[edges];
		for(int i=0; i<indices.length; i++){
			indices[i]=MODE_ORDERED.equals(mode) ? i : i*2;
		}
		return indices;
	}

	/** Creates exact small integer weights for one sparse row.
	 * @param edges Number of weights
	 * @return Generated weights
	 */
	private static float[] makeWeights(final int edges){
		final float[] weights=new float[edges];
		for(int i=0; i<weights.length; i++){
			int weight=(i*5)%7-3;
			if(weight==0){weight=4;}
			weights[i]=weight;
		}
		return weights;
	}

	/** Checks whether a benchmark mode is supported.
	 * @param mode Mode name
	 * @return True if supported
	 */
	private static boolean validMode(final String mode){
		return MODE_ORDERED.equals(mode) || MODE_PRUNED.equals(mode);
	}

	/** Mixes a retained integer into the checksum.
	 * @param checksum Previous checksum
	 * @param value Retained value
	 * @return Updated checksum
	 */
	private static long mix(long checksum, final int value){
		checksum^=(value&0xFFFFFFFFL)+0x9E3779B97F4A7C15L;
		checksum*=0xBF58476D1CE4E5B9L;
		return Long.rotateLeft(checksum, 11);
	}

	/** Prevents construction of this standalone benchmark. */
	private CellNetSparseInferenceBenchmark(){}

	/** Mutable inference inputs and independently retained exact outputs. */
	private static final class Workload{
		/** Preallocates benchmark state.
		 * @param width Input width
		 * @param indices_ Sparse input indices
		 * @param weights_ Sparse input weights
		 */
		private Workload(final int width, final int[] indices_, final float[] weights_){
			inputs=new float[BANKS][width];
			expected=new float[BANKS];
			cursors=new int[BANKS];
			indices=indices_;
			weights=weights_;
		}

		/** Mutable input banks passed directly to CellNet.applyInput. */
		final float[][] inputs;
		/** Exact linear output expected for each bank. */
		final float[] expected;
		/** Per-bank mutation cursor into indices and weights. */
		final int[] cursors;
		/** Generated sparse input indices. */
		final int[] indices;
		/** Generated sparse input weights. */
		final float[] weights;
	}

	/** Number of mutable input buffers. */
	private static final int BANKS=128;
	/** Warmup predictions run before timing starts. */
	private static final int WARMUP_LOOPS=8192;
	/** Minimum width that exercises the 16-wide sparse-vector threshold. */
	private static final int MIN_WIDTH=16;
	/** Maximum width preserving exact integer products, sums and oracle updates. */
	private static final int MAX_WIDTH=32768;
	/** Minimum generated input value. */
	private static final int MIN_INPUT=-127;
	/** Maximum generated input value. */
	private static final int MAX_INPUT=127;
	/** Exact small integer linear bias. */
	private static final float BIAS=3;
	/** Number of rows before the sparse I and W rows. */
	private static final int MODEL_HEADER_ROWS=10;
	/** Ordered full-width identity sparse-row benchmark. */
	private static final String MODE_ORDERED="ordered";
	/** Genuinely pruned indexed sparse-row control. */
	private static final String MODE_PRUNED="pruned";
}
