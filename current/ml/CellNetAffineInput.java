package ml;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedHashMap;

import fileIO.ByteStreamWriter;

/** Converts the declared two-layer explicit-affine normalization to float32 input headers.
 * Trained weights and biases are copied without folding or quantization.
 * @author Ganyu
 */
public final class CellNetAffineInput {

	public static void main(String[] args) throws IOException{
		String in=null,out=null,x=null;
		int rows=1000;
		for(String arg:args){
			if(arg.startsWith("in=") && in==null){in=arg.substring(3);}
			else if(arg.startsWith("out=") && out==null){out=arg.substring(4);}
			else if(arg.startsWith("x=") && x==null){x=arg.substring(2);}
			else if(arg.startsWith("rows=")){rows=Integer.parseInt(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown/duplicate argument: "+arg);}
		}
		if(in==null || out==null || x==null || rows<1){throw new IllegalArgumentException("Require in= out=fresh x=NPY rows=1000");}
		final Path destination=Paths.get(out);
		if(Files.exists(destination)){throw new IllegalArgumentException("Fresh output required: "+out);}
		final CellNet source=CellNetParser.load(in,false), converted=convert(source);
		CellNet.codingA48Out=true; CellNet.OUT_DENSE=false; CellNet.OUT_SPARSE=false;
		final Path temporary=Files.createTempFile(destination.toAbsolutePath().getParent(),".affine-input-",".partial");
		double maximum=0;
		try{
			final ByteStreamWriter writer=new ByteStreamWriter(temporary.toString(),true,false,false);
			writer.start();
			try{writer.print(converted.toBytes());}finally{if(writer.poisonAndWait()){throw new IOException("Affine conversion output failed");}}
			final CellNet reloaded=CellNetParser.load(temporary.toString(),false);
			try(CellNetWeightBitsEval.Npy data=new CellNetWeightBitsEval.Npy(x)){
				if(data.cols!=source.numInputs() || data.rows<rows){throw new IllegalArgumentException("Parity matrix shape/row count mismatch");}
				final float[] input=new float[data.cols];
				for(int row=0; row<rows; row++){
					data.read(input); source.applyInput(input); reloaded.applyInput(input);
					source.feedForward(); reloaded.feedForward();
					for(int k=0; k<source.numOutputs(); k++){
						final float a=source.getOutput(k), b=reloaded.getOutput(k);
						if(!Float.isFinite(a) || !Float.isFinite(b)){throw new IllegalStateException("Nonfinite affine parity output");}
						maximum=Math.max(maximum,Math.abs((double)a-b));
						if(Float.floatToRawIntBits(a)!=Float.floatToRawIntBits(b)){
							throw new IllegalStateException("Affine conversion must preserve exact float32 output: row="+row+" head="+k+" original="+a+" converted="+b);
						}
					}
				}
			}
			Files.createLink(destination,temporary);
		}finally{Files.deleteIfExists(temporary);}
		System.out.println("AFFINE_INPUT_PASS rows="+rows+" outputs="+source.numOutputs()+" maxAbs="+maximum+" trainedEdges="+converted.countEdges());
	}

	/** Only the explicitly documented diagonal subtract-then-scale representation is accepted. */
	static CellNet convert(CellNet source){
		if(!"explicit-affine".equals(source.tags.get("normalization")) || !"2".equals(source.tags.get("affine_prefix_layers"))
				|| source.hasInputNormalization() || source.dims.length<4 || source.dims[0]!=source.dims[1] || source.dims[0]!=source.dims[2]){
			throw new IllegalArgumentException("Require declared two-layer explicit-affine input normalization and no existing input headers");
		}
		final int width=source.numInputs();
		final float[] mean=new float[width], inv=new float[width];
		for(int i=0; i<width; i++){
			final Cell subtract=source.net[1][i], scale=source.net[2][i];
			final float first=diagonal(subtract,i,width), second=diagonal(scale,i,width);
			if(first!=1f || Float.floatToRawIntBits(scale.bias())!=0 || !Float.isFinite(subtract.bias()) || !Float.isFinite(second) || second<=0){
				throw new IllegalArgumentException("Unexpected affine bias/scale at input="+i);
			}
			mean[i]=-subtract.bias(); inv[i]=second;
		}
		final int[] dims=new int[source.dims.length-2]; dims[0]=width;
		System.arraycopy(source.dims,3,dims,1,dims.length-1);
		final CellNet target=new CellNet(dims,source.seed,source.density,source.density1,source.edgeBlockSize,source.commands);
		target.tags=new LinkedHashMap<String,String>(source.tags);
		target.tags.put("normalization","explicit-input");
		target.tags.remove("affine_prefix_layers"); target.tags.remove("affine_prefix_semantics");
		target.tags.put("native_affine_conversion","two diagonal LINEAR layers to exact float32 input headers");
		target.epochsTrained=source.epochsTrained;
		target.setInputNormalization(mean,inv);
		for(int layer=1; layer<target.net.length; layer++){
			for(int i=0; i<target.net[layer].length; i++){
				final Cell from=source.net[layer+2][i], to=target.net[layer][i];
				to.weights=from.weights.clone(); to.inputs=from.inputs==null ? null : from.inputs.clone();
				to.bias=from.bias; to.function=from.function;
				assert(Arrays.equals(to.weights,from.weights)) : "Normalization conversion must not fold or change any trained edge";
			}
		}
		// CellNet.toBytes/countEdges require matrix references; sparse copies also need reverse connectivity.
		if(!CellNet.DENSE){CellNet.makeOutputSets(target.net);}
		target.makeWeightMatrices();
		return target;
	}

	/** Validates a diagonal LINEAR cell without assuming dense/sparse serialization. */
	private static float diagonal(Cell cell,int index,int width){
		if(!"LINEAR".equalsIgnoreCase(cell.function.name())){throw new IllegalArgumentException("Affine prefix must be LINEAR");}
		if(cell.inputs!=null){
			if(cell.inputs.length!=1 || cell.weights.length!=1 || cell.inputs[0]!=index){throw new IllegalArgumentException("Non-diagonal sparse affine prefix");}
			return cell.weights[0];
		}
		if(cell.weights.length!=width){throw new IllegalArgumentException("Affine prefix width mismatch");}
		for(int j=0; j<width; j++){if(j!=index && cell.weights[j]!=0f){throw new IllegalArgumentException("Non-diagonal dense affine prefix");}}
		return cell.weights[index];
	}
}
