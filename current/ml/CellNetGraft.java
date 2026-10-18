package ml;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Evaluates new private nodes and heads from a completed, unchanged base forward pass.
 * The three borrowed networks are mutable inference workers used serially by one
 * owner of this graft object. The base's topology, parameters and activations are
 * read-only here. Serving preserves FP32 activations; the trainer's reduced-precision
 * cache is not part of inference.
 * @author Brian Bushnell, Yoimiya */
public final class CellNetGraft{

	/** Validates geometry and per-head private edge partitions before allocating reusable scratch. */
	public CellNetGraft(final CellNet base_, final CellNet private_, final CellNet heads_,
			final int upstreamLayer_, final int shared_, final int privatePerHead_, final String activation){
		require(base_!=null && private_!=null && heads_!=null, "Graft requires three private inference workers");
		base=base_; privateNodes=private_; heads=heads_; upstreamLayer=upstreamLayer_; shared=shared_;
		require(upstreamLayer>0 && upstreamLayer==base.layers-3,
			"Graft upstream must be the input to the base's semifinal layer");
		require(shared>0 && shared<=base.dims[upstreamLayer+1] && privatePerHead_>0,
			"Graft shared/private widths exceed the declared base geometry");
		require(privateNodes.layers==2 && heads.layers==2, "Graft components must be single affine layers");
		requireIdentity(privateNodes); requireIdentity(heads);
		require(privateNodes.numInputs()==base.dims[upstreamLayer] &&
			privateNodes.numOutputs()==heads.numOutputs()*privatePerHead_ &&
			heads.numInputs()==shared+privateNodes.numOutputs(), "Graft component and bypass dimensions disagree");
		if("tanh".equals(activation)){hidden=Function.getFunction(Function.TANH);}
		else if("swish".equals(activation)){hidden=Function.getFunction(Function.SWISH);}
		else{throw new IllegalArgumentException("Unsupported graft hidden activation: "+activation);}
		for(Cell cell:privateNodes.finalLayer){
			require(cell.function.type()==Function.LINEAR, "Private affine component must defer its hidden activation");
		}
		for(int row=0; row<heads.numOutputs(); row++){
			final Cell cell=heads.finalLayer[row];
			require(cell.function.type()==Function.SIG, "Graft endpoint heads require their trained sigmoid activation");
			for(int edge=0; edge<cell.weights.length; edge++){
				final int column=cell.inputs==null ? edge : cell.inputs[edge];
				final boolean own=column<shared || (column>=shared+row*privatePerHead_ && column<shared+(row+1)*privatePerHead_);
				require(own || cell.weights[edge]==0, "Graft head has a connection to another head's private nodes");
			}
		}
		upstream=new float[privateNodes.numInputs()]; joined=new float[heads.numInputs()];
	}

	/** The shared exporter stores explicit zero/one normalization even on affine components. */
	private static void requireIdentity(final CellNet model){
		if(!model.hasInputNormalization()){return;}
		final float[] mean=model.inputMeanCopy(), inverse=model.inputInverseStdCopy();
		for(int i=0; i<mean.length; i++){
			require(Float.floatToRawIntBits(mean[i])==0 && inverse[i]==1,
				"Graft component normalization must be positive zero and unit inverse deviation");
		}
	}

	/** Requires the caller to complete the base pass and copy base outputs first; writes only suffix heads. */
	public void appendAfterBase(final float[] output, final int offset){
		require(output!=null && offset==base.numOutputs() && output.length>=offset+heads.numOutputs(),
			"Graft outputs must follow all original base heads");
		final float[] source=base.values[upstreamLayer], semifinal=base.values[upstreamLayer+1];
		System.arraycopy(source, 0, upstream, 0, upstream.length);
		System.arraycopy(semifinal, 0, joined, 0, shared);
		privateNodes.applyInput(upstream); privateNodes.feedForward();
		for(int i=0; i<privateNodes.numOutputs(); i++){
			joined[shared+i]=(float)hidden.activate(privateNodes.getOutput(i));
		}
		heads.applyInput(joined); heads.feedForward();
		for(int i=0; i<heads.numOutputs(); i++){
			final float value=heads.getOutput(i);
			if(!Float.isFinite(value)){throw new IllegalStateException("Nonfinite graft output at "+i);}
			output[offset+i]=value;
		}
	}

	/** Writes the last completed inference snapshot for the matching caller-supplied raw input.
	 * The destination must be new, so repeated bin IDs or accidental reruns fail loudly.
	 * Freshness is the caller's responsibility; no inference is repeated or verified here.
	 * Each float has a round-trip decimal value and its exact signed integer bit pattern. */
	public void writeActivationTrace(final String destination, final float[] raw){
		require(destination!=null && raw!=null && raw.length==base.numInputs(),
			"Activation tracing requires a fresh destination and the actual composite input");
		final ByteBuilder out=new ByteBuilder();
		out.append("#columns\tstage\tindex\tvalue\tfloat_bits\n");
		appendTrace(out, "raw_input", raw);
		for(int layer=0; layer<base.values.length; layer++){
			appendTrace(out, "base_layer_"+layer, base.values[layer]);
		}
		appendTrace(out, "upstream", upstream);
		appendTrace(out, "private_affine", privateNodes.values[1]);
		appendTrace(out, "joined", joined);
		appendTrace(out, "new_heads", heads.values[1]);
		try{Files.createFile(Paths.get(destination));}
		catch(IOException e){throw new IllegalArgumentException("Could not create new activation trace: "+destination, e);}
		final ByteStreamWriter writer=new ByteStreamWriter(destination, true, false, true);
		writer.start(); writer.print(out);
		if(writer.poisonAndWait()){throw new IllegalStateException("Could not write activation trace: "+destination);}
	}

	/** Appends one traced activation vector with both decimal and exact FP32 bit text. */
	private static void appendTrace(final ByteBuilder out, final String stage, final float[] values){
		assert(values!=null) : "Completed graft inference owns every traced activation array";
		for(int i=0; i<values.length; i++){
			out.append(stage).tab().append(i).tab().append(Float.toString(values[i])).tab()
				.append(Float.floatToRawIntBits(values[i])).nl();
		}
	}

	/** Rejects invalid external geometry even with assertions disabled. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Borrowed mutable inference networks; this graft serially reuses their activation arrays. */
	private final CellNet base, privateNodes, heads;
	/** Layer and prefix widths defining the base-to-private bypass geometry. */
	private final int upstreamLayer, shared;
	/** Activation applied between private affine nodes and endpoint heads. */
	private final Function hidden;
	/** Reusable per-call scratch; this object is not thread-safe. */
	private final float[] upstream, joined;
}
