package test.ml;

import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/** Temporary diagnostic only: preserves the failed parity result and measures its cause.
 * Compile on Dori against the pinned E97 classes; never install into the production tree.
 * @author Yoimiya */
public final class SubnetCopyDiagnostic {

	public static void main(String[] args){
		final PreParser pre=new PreParser(args, SubnetCopyDiagnostic.class, false);
		try{run(pre.args);}finally{Shared.closeStream(pre.outstream);}
	}

	private static void run(String[] args){
		String model=null,vectors=null,predictions=null;
		for(String arg:args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("net")){model=value;}
			else if(key.equals("vectors")){vectors=value;}
			else if(key.equals("predictions")){predictions=value;}
			else{throw new IllegalArgumentException("Unknown argument: "+key);}
		}
		if(model==null || vectors==null || predictions==null){
			throw new IllegalArgumentException("Required net= vectors= predictions=");
		}
		Parser.parseCommonStatic("forceJavaParseDouble=t", "forcejavaparsedouble", "t");
		final CellNet inference=CellNetParser.loadInferenceFromLines(ByteFile.toLines(model));
		SubnetGpuRoundTripTest.compare(inference, vectors, predictions, 70, .00002, null);
		final CellNet loaded=CellNetParser.load(model, false);
		SubnetGpuRoundTripTest.compare(loaded, vectors, predictions, 70, .00002, inference);
		final CellNet copy=loaded.copy(false);
		parameters(loaded, copy);
		System.out.println("flags\tSIMD="+Shared.SIMD+"\tFF="+Shared.SIMD_FEED_FORWARD+
			"\tFMA="+Shared.SIMD_FMA+"\tSPECIAL_FMA="+CellNet.SPECIAL_FMA+
			"\tsourceFF="+loaded.simdFF+"\tcopyFF="+copy.simdFF+"\tedgeBlockSize="+loaded.edgeBlockSize);
		final ByteFile bf=ByteFile.makeByteFile(vectors, true);
		final LineParser1 parser=new LineParser1((byte)'\t');
		final float[][] rows=new float[70][loaded.numInputs()];
		try{
			parser.set(bf.nextLine());
			if(!parser.termEquals("#dims", 0) || parser.parseInt(1)!=loaded.numInputs()){
				throw new IllegalArgumentException("Diagnostic requires the original input width");
			}
			for(int r=0; r<rows.length; r++){
				byte[] line=bf.nextLine();
				while(line!=null && (line.length==0 || line[0]=='#')){line=bf.nextLine();}
				if(line==null){throw new IllegalArgumentException("Missing row "+r);}
				parser.set(line);
				for(int i=0; i<rows[r].length; i++){
					rows[r][i]=parser.parseFloat(i);
					if(!Float.isFinite(rows[r][i])){throw new IllegalArgumentException("Nonfinite input");}
				}
			}
		}finally{if(bf.close()){throw new IllegalStateException("Vector close failed");}}
		// Fixed diagnostic passes, not retries seeking a passing acceptance result.
		inspect(loaded, copy, rows, 0, false);
		for(int pass=1; pass<=3; pass++){inspect(loaded, copy, rows, pass, true);}
		parameters(loaded, copy);
		System.out.println("DIAGNOSTIC_COMPLETE_NOT_PARITY_ACCEPTANCE");
	}

	/** Compare every persisted forward parameter before attributing a discrepancy to arithmetic. */
	private static void parameters(CellNet a, CellNet b){
		assert(a.net.length==b.net.length) : "copy(false) must retain layer topology for worker inference";
		if(a.dense!=b.dense || a.edgeBlockSize!=b.edgeBlockSize ||
			!Arrays.equals(a.dims,b.dims) || a.simdFF!=b.simdFF ||
			!Arrays.equals(a.inputMeanCopy(), b.inputMeanCopy()) ||
			!Arrays.equals(a.inputInverseStdCopy(), b.inputInverseStdCopy()) ||
			Float.floatToRawIntBits(a.inputClip())!=Float.floatToRawIntBits(b.inputClip())){
			throw new IllegalStateException("Network shape, mode or preprocessing differs");
		}
		for(int l=0; l<a.net.length; l++){
			for(int i=0; i<a.net[l].length; i++){
				final Cell x=a.net[l][i],y=b.net[l][i];
				if(x.function!=y.function || Float.floatToRawIntBits(x.bias)!=Float.floatToRawIntBits(y.bias)
						|| !Arrays.equals(x.weights,y.weights) || !Arrays.equals(x.inputs,y.inputs)){
					throw new IllegalStateException("Forward parameter mismatch layer="+l+" cell="+i);
				}
			}
		}
		System.out.println("forward_parameters\tEXACT");
	}

	/** Inspect first differing layers, output error, and repeat stability on identical inputs. */
	private static void inspect(CellNet loaded, CellNet copy, float[][] rows, int pass, boolean repeat){
		assert(rows.length==70) : "Diagnostic population must remain the original 64 held-out and six edge rows";
		long different=0,unstableLoaded=0,unstableCopy=0,argmaxChanges=0;
		double max=0;
		for(int r=0; r<rows.length; r++){
			copy.applyInput(rows[r]); copy.feedForward();
			loaded.applyInput(rows[r]); loaded.feedForward();
			boolean reported=false;
			for(int l=0; l<loaded.values.length; l++){
				for(int i=0; i<loaded.values[l].length; i++){
					final float x=loaded.values[l][i],y=copy.values[l][i];
					if(!Float.isFinite(x) || !Float.isFinite(y)){throw new IllegalStateException("Nonfinite activation");}
					if(Float.floatToRawIntBits(x)!=Float.floatToRawIntBits(y) && !reported){
						if(different<12){System.out.println("first_layer_difference\tpass="+pass+"\trow="+r+
							"\tlayer="+l+"\tcell="+i+"\tsource="+Float.toHexString(x)+
							"\tcopy="+Float.toHexString(y)+"\tsourceSum="+Double.toHexString(loaded.net[l][i].sum)+
							"\tcopySum="+Double.toHexString(copy.net[l][i].sum));}
						reported=true;
					}
				}
			}
			final float[] a=loaded.values[loaded.values.length-1].clone();
			final float[] b=copy.values[copy.values.length-1].clone();
			int ai=0,bi=0;
			for(int i=0; i<a.length; i++){
				if(Float.floatToRawIntBits(a[i])!=Float.floatToRawIntBits(b[i])){different++;}
				max=Math.max(max,Math.abs((double)a[i]-b[i]));
				if(a[i]>a[ai]){ai=i;} if(b[i]>b[bi]){bi=i;}
			}
			if(ai!=bi){argmaxChanges++;}
			if(repeat){
				loaded.applyInput(rows[r]); loaded.feedForward();
				copy.applyInput(rows[r]); copy.feedForward();
				for(int i=0; i<a.length; i++){
					if(Float.floatToRawIntBits(a[i])!=Float.floatToRawIntBits(loaded.getOutput(i))){unstableLoaded++;}
					if(Float.floatToRawIntBits(b[i])!=Float.floatToRawIntBits(copy.getOutput(i))){unstableCopy++;}
				}
			}
		}
		System.out.println("pass\t"+pass+"\trepeat_evaluations="+repeat+"\tdifferent_output_values="+different+"\tmax_absolute="+max+
			"\targmax_changes="+argmaxChanges+"\tsource_repeat_changes="+unstableLoaded+
			"\tcopy_repeat_changes="+unstableCopy);
	}
}
