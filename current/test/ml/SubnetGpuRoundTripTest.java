package test.ml;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Compares exported Java predictions with saved GPU predictions or a second
 * native network. Reference-network mode requires exact output bits and also
 * checks that compression/metadata changes preserve all declared output heads.
 * Raw validation inputs are read directly, including their declared target columns;
 * Only input columns enter CellNet. Declared training targets may differ from the
 * network's output count because live error heads need no extra stored labels.
 * No normalization is recomputed here.
 * Adapted from Yoimiya's fixed donor checker; ordinary=t additionally exercises
 * the ordinary file loader and the caller's worker-owned copy(false) path.
 * @author Yoimiya, Raiden
 */
public final class SubnetGpuRoundTripTest {

	/** Requires dimensions and either saved predictions or an exact native reference. */
	public static void main(String[] args){
		final PreParser pre=new PreParser(args, SubnetGpuRoundTripTest.class, false);
		try{run(pre.args);}finally{Shared.closeStream(pre.outstream);}
	}

	/** Keeps the original comparison/report behavior after standard argument handling. */
	private static void run(String[] args){
		String netFile=null,vectors=null,predictions=null,referenceFile=null,outBits=null;
		int rows=1000;
		boolean ordinary=false,argmax=false;
		double tolerance=0.00002;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String a=arg.substring(0,eq).toLowerCase(Locale.ROOT),b=arg.substring(eq+1);
			if(a.equals("net")){netFile=b;}
			else if(a.equals("vectors")){vectors=b;}
			else if(a.equals("predictions")){predictions=b;}
			else if(a.equals("reference")){referenceFile=b;}
			else if(a.equals("outbits")){outBits=b;}
			else if(a.equals("rows")){rows=Parse.parseIntKMG(b);}
			else if(a.equals("tolerance")){tolerance=Double.parseDouble(b);}
			else if(a.equals("ordinary")){ordinary=Parse.parseBoolean(b);}
			else if(a.equals("argmax")){argmax=Parse.parseBoolean(b);}
			// Exported TSV values may use exponents; PreParser does not consume this common flag.
			else if(a.equals("forcejavaparsedouble") && Parser.parseCommonStatic(arg,a,b)){}
			else{throw new IllegalArgumentException("Unknown argument: "+a);}
		}
		if(netFile==null || vectors==null || (predictions==null)==(referenceFile==null) || rows<1 ||
				!Double.isFinite(tolerance) || tolerance<0 || (predictions!=null && tolerance==0)){
			throw new IllegalArgumentException("Required: net= vectors= and exactly one of predictions= or reference=; positive rows and valid tolerance");
		}
		final CellNet net=CellNetParser.loadInferenceFromLines(ByteFile.toLines(netFile));
		if(net==null){throw new IllegalArgumentException("Could not load network: "+netFile);}
		assert(net.numInputs()>0 && net.numOutputs()>0) : "CellNet inference requires nonempty input and output layers";
		final CellNet reference=referenceFile==null ? null : CellNetParser.loadInferenceFromLines(ByteFile.toLines(referenceFile));
		if(referenceFile!=null && reference==null){throw new IllegalArgumentException("Could not load reference network: "+referenceFile);}
		final long[] changes=argmax ? new long[1] : null;
		try(PredictionBits captured=outBits==null ? null : new PredictionBits(outBits,net.numOutputs())){
			final Comparison first=comparePaths(net,vectors,predictions,rows,tolerance,reference,true,captured);
			final double maximum=first.maximum;
			if(changes!=null){changes[0]=first.oracleArgmaxChanges;}
			Comparison ordinaryResult=null,copyResult=null;
			if(ordinary){
				final CellNet loaded=CellNetParser.load(netFile, false);
				if(predictions==null){
					// Keep the historical reference= exact gate unchanged.
					compare(loaded, vectors, null, rows, tolerance, net);
					compare(loaded.copy(false), vectors, null, rows, tolerance, loaded);
				}else{
					final CellNet copy=loaded.copy(false);
					requireForwardCopy(loaded, copy);
					ordinaryResult=compareExportPaths(loaded, vectors, predictions, rows, tolerance, net);
					copyResult=compareExportPaths(copy, vectors, predictions, rows, tolerance, loaded);
					requireForwardCopy(loaded, copy);
					// The scalar control isolates copy correctness from the Vector API's
					// permitted time-varying ADD reduction order. Restore global policy.
					final boolean simd=Shared.SIMD;
					try{
						Shared.SIMD=false;
						compare(copy, vectors, null, rows, 0, loaded);
					}finally{Shared.SIMD=simd;}
				}
			}
			if(captured!=null){captured.publish(rows);}
			System.out.println("status\t"+(referenceFile==null ? "JAVA_GPU_PREDICTION_PASS" : "JAVA_JAVA_EXACT_PREDICTION_PASS"));
			System.out.println("rows\t"+rows);
			System.out.println("inputs\t"+net.numInputs());
			System.out.println("outputs\t"+net.numOutputs());
			System.out.println("max_absolute_difference\t"+maximum);
			System.out.println("tolerance\t"+tolerance);
			if(ordinary){System.out.println("ordinary_load_and_worker_copy\tPASS");}
			if(ordinaryResult!=null){
				System.out.println("forward_copy_parameters\tEXACT");
				System.out.println("scalar_worker_copy\tEXACT");
				System.out.println("accelerated_native_comparison\tBOUNDED_NOT_BIT_EXACT");
				ordinaryResult.print("ordinary_vs_inference");
				copyResult.print("copy_vs_ordinary");
			}
			if(argmax){System.out.println("first_maximum_changes\t"+changes[0]);}
		}
	}

	/** Streams aligned raw/GPU rows and fails on truncation, extra predictions, or mismatch. */
	private static double compare(CellNet net,String vectors,String predictions,int rows,double tolerance){
		return compare(net,vectors,predictions,rows,tolerance,null);
	}

	/** Checks exact native equality with an optional additional saved-prediction bound. */
	static double compare(CellNet net,String vectors,String predictions,int rows,double tolerance,CellNet loose){
		return compare(net, vectors, predictions, rows, tolerance, loose, null);
	}

	/** Optional first-maximum differences are reported independently of numeric tolerance. */
	static double compare(CellNet net,String vectors,String predictions,int rows,double tolerance,CellNet loose,long[] changes){
		final Comparison result=comparePaths(net, vectors, predictions, rows, tolerance, loose, true);
		if(changes!=null){changes[0]+=result.oracleArgmaxChanges;}
		return result.maximum;
	}

	/** New-export mode alone permits bounded native differences; reference= stays exact.
	 * JDK25 FloatVector.reduceLanes(ADD) allows its floating-point reduction order
	 * to vary over time. This is not permission to change copied model parameters. */
	static Comparison compareExportPaths(CellNet net,String vectors,String predictions,int rows,double tolerance,CellNet loose){
		if(predictions==null || loose==null || !Double.isFinite(tolerance) || tolerance<=0){
			throw new IllegalArgumentException("Bounded native comparison requires both networks and the exported oracle");
		}
		return comparePaths(net, vectors, predictions, rows, tolerance, loose, false);
	}

	/** Stream both native paths against the oracle, recording every actual output. */
	private static Comparison comparePaths(CellNet net,String vectors,String predictions,int rows,double tolerance,CellNet loose,boolean exactNative){
		return comparePaths(net,vectors,predictions,rows,tolerance,loose,exactNative,null);
	}

	/** Captures this very inference pass; no recomputation may replace the checked bits. */
	private static Comparison comparePaths(CellNet net,String vectors,String predictions,int rows,double tolerance,
			CellNet loose,boolean exactNative,PredictionBits captured){
		if(predictions==null && loose==null){throw new IllegalArgumentException("A prediction file or reference network is required");}
		if(loose!=null && (loose.numInputs()!=net.numInputs() || loose.numOutputs()!=net.numOutputs())){
			throw new IllegalArgumentException("Loose and bundled network dimensions differ");
		}
		final ByteFile vf=ByteFile.makeByteFile(vectors,true);
		ByteFile pf=null;
		try{
			pf=predictions==null ? null : ByteFile.makeByteFile(predictions,true);
			final LineParser1 v=new LineParser1((byte)'\t'),p=new LineParser1((byte)'\t');
			final byte[] header=vf.nextLine();
			if(header==null){throw new IllegalArgumentException("Empty vector file");}
			v.set(header);
			if(v.terms()!=4 || !v.termEquals("#dims",0) || v.parseInt(1)!=net.numInputs()
					|| v.parseInt(2)<0 || v.parseInt(3)!=0){
				throw new IllegalArgumentException("Vector dimensions do not match exported network inputs");
			}
			final int targetColumns=v.parseInt(2);
			// Live residual and classification heads derive targets from the two stored value labels.
			if(targetColumns!=net.numOutputs() && !(targetColumns==2 &&
					(net.numOutputs()==4 || net.numOutputs()==6))){
				throw new IllegalArgumentException("Stored targets must match network outputs, or use two value labels for four or six live outputs");
			}
			final float[] input=new float[net.numInputs()];
			final Comparison result=new Comparison();
			double maximum=0;
			for(int row=0; row<rows; row++){
				final byte[] vl=nextData(vf),pl=pf==null ? null : nextData(pf);
				if(vl==null || (pf!=null && pl==null)){throw new IllegalArgumentException("Missing raw/reference row "+row);}
				v.set(vl); if(pf!=null){p.set(pl);}
				if(v.terms()!=net.numInputs()+targetColumns || (pf!=null && (p.terms()!=net.numOutputs()+1
						|| p.parseInt(0)!=row))){
					throw new IllegalArgumentException("Input/prediction shape or row identity mismatch at "+row);
				}
				for(int i=0; i<input.length; i++){
					input[i]=v.parseFloat(i);
					if(!Float.isFinite(input[i])){throw new IllegalArgumentException("Nonfinite input at row "+row+" column "+i);}
				}
				net.applyInput(input);
				net.feedForward();
				if(loose!=null){loose.applyInput(input); loose.feedForward();}
				int actualBest=0,expectedBest=0,looseBest=0;
				double actualMaximum=Double.NEGATIVE_INFINITY,expectedMaximum=Double.NEGATIVE_INFINITY,
					looseMaximum=Double.NEGATIVE_INFINITY;
				for(int i=0; i<net.numOutputs(); i++){
					if(exactNative && loose!=null && Float.floatToIntBits(net.getOutput(i))!=Float.floatToIntBits(loose.getOutput(i))){
						throw new IllegalStateException("Loose/bundled output bits differ at row="+row+" output="+i);
					}
					final double expected=pf==null ? loose.getOutput(i) : p.parseDouble(i+1),actual=net.getOutput(i);
					final double error=Math.abs(actual-expected);
					if(!Double.isFinite(error) || error>tolerance){
						throw new IllegalStateException("GPU/Java prediction mismatch row="+row+" output="+i+
							" expected="+expected+" actual="+actual+" absolute_error="+error+" tolerance="+tolerance);
					}
					maximum=Math.max(maximum,error);
					if(loose!=null){
						final double other=loose.getOutput(i),nativeError=Math.abs(actual-other),oracleError=Math.abs(other-expected);
						if(!Double.isFinite(nativeError) || !Double.isFinite(oracleError) ||
								(!exactNative && (nativeError>tolerance || oracleError>tolerance))){
							throw new IllegalStateException("Native/reference prediction mismatch row="+row+" output="+i+
								" native_error="+nativeError+" reference_oracle_error="+oracleError+" tolerance="+tolerance);
						}
						result.nativeMaximum=Math.max(result.nativeMaximum,nativeError);
						result.referenceMaximum=Math.max(result.referenceMaximum,oracleError);
						if(other>looseMaximum){looseMaximum=other; looseBest=i;}
					}
					if(actual>actualMaximum){actualMaximum=actual; actualBest=i;}
					if(expected>expectedMaximum){expectedMaximum=expected; expectedBest=i;}
				}
				if(actualBest!=expectedBest){result.oracleArgmaxChanges++;}
				if(loose!=null){
					if(actualBest!=looseBest){result.nativeArgmaxChanges++;}
					if(looseBest!=expectedBest){result.referenceOracleArgmaxChanges++;}
				}
				if(captured!=null){captured.add(row,net);}
			}
			if(pf!=null && nextData(pf)!=null){throw new IllegalArgumentException("GPU predictions exceed requested row count "+rows);}
			assert(maximum<=tolerance) : "Every compared output must satisfy the absolute-error bound";
			result.maximum=maximum;
			return result;
		}finally{
			//TODO: Probable bug inherited from the donor: a close failure can mask the
			//original comparison failure. This scoped adaptation leaves comparison cleanup unchanged.
			final boolean vError=vf.close(),pError=pf!=null && pf.close();
			if(vError || pError){throw new IllegalStateException("I/O error closing prediction comparison inputs");}
		}
	}

	/** Lossless first-pass output. A failed run retains only the clearly partial file. */
	static final class PredictionBits implements AutoCloseable{
		PredictionBits(String name,int outputs_){
			if(outputs_<1){throw new IllegalArgumentException("Prediction capture requires output heads");}
			outputs=outputs_;
			target=Paths.get(name).toAbsolutePath().normalize();
			partial=target.resolveSibling(target.getFileName()+".partial");
			if(Files.exists(target,LinkOption.NOFOLLOW_LINKS)){
				throw new IllegalArgumentException("Refusing to replace prediction output: "+target);
			}
			try{Files.createFile(partial);}
			catch(IOException e){throw new IllegalStateException("Could not create fresh partial prediction output: "+partial,e);}
			try{
				writer=new ByteStreamWriter(partial.toString(),true,false,false);
				writer.start();
			}catch(RuntimeException | Error e){
				//TODO: Probable bug - failed Thread.start leaves the constructor-opened stream unclosed.
				//NcrnaTextIO.Output.close uses closeBeforeStart for this case; repair separately under review.
				System.err.println("unaccepted_prediction_bits\t"+partial);
				throw e;
			}
			buffer.append("#prediction_bits\t1\t").append(outputs).nl();
		}

		void add(int row,CellNet net){
			assert(!closed && row==written && net.numOutputs()==outputs) :
				"Checked prediction rows must retain their original order and output width";
			buffer.append(row);
			for(int i=0; i<outputs; i++){buffer.tab().append(Float.floatToRawIntBits(net.getOutput(i)));}
			buffer.nl();
			written++;
			if(buffer.length()>=65536){writer.print(buffer); buffer.clear();}
		}

		void publish(int expectedRows){
			if(written!=expectedRows || published){throw new IllegalStateException("Prediction output is incomplete or already published");}
			finish();
			try{
				// Same-directory hard-link publication is atomic and fails if the target
				// appeared meanwhile. It never replaces somebody else's accepted output.
				Files.createLink(target,partial);
				published=true;
			}catch(IOException e){throw new IllegalStateException("Could not publish checked prediction output: "+target,e);}
			try{Files.delete(partial);}
			catch(IOException e){System.err.println("Checked output published; redundant partial retained: "+partial+"; "+e);}
		}

		private void finish(){
			if(closed){return;}
			closed=true;
			Throwable failure=null;
			try{
				if(buffer.length()>0){writer.print(buffer); buffer.clear();}
			}catch(RuntimeException | Error e){failure=e; throw e;}
			finally{
				try{
					if(writer.poisonAndWait()){throw new IllegalStateException("Could not close prediction output: "+partial);}
				}catch(RuntimeException | Error e){
					if(failure==null){throw e;}
					if(e!=failure){failure.addSuppressed(e);}
				}
			}
		}

		@Override
		public void close(){
			try{finish();}
			finally{if(!published){System.err.println("unaccepted_prediction_bits\t"+partial);}}
		}

		private final Path target,partial;
		private final int outputs;
		private final ByteStreamWriter writer;
		private final ByteBuilder buffer=new ByteBuilder(65536);
		private int written;
		private boolean closed,published;
	}

	/** A numerical allowance must never conceal a changed forward model in copy(false). */
	static void requireForwardCopy(CellNet source,CellNet copy){
		if(source==null || copy==null || source.dense!=copy.dense || source.simdFF!=copy.simdFF ||
				source.edgeBlockSize!=copy.edgeBlockSize || !Arrays.equals(source.dims,copy.dims) ||
				!sameBits(source.inputMeanCopy(),copy.inputMeanCopy()) ||
				!sameBits(source.inputInverseStdCopy(),copy.inputInverseStdCopy()) ||
				Float.floatToRawIntBits(source.inputClip())!=Float.floatToRawIntBits(copy.inputClip())){
			throw new IllegalStateException("Worker copy changed forward shape, mode or normalization");
		}
		assert(source.net.length==copy.net.length) : "Matched dimensions must retain every layer for copy inference";
		for(int layer=0; layer<source.net.length; layer++){
			for(int i=0; i<source.net[layer].length; i++){
				final Cell a=source.net[layer][i],b=copy.net[layer][i];
				if(a.function!=b.function || Float.floatToRawIntBits(a.bias)!=Float.floatToRawIntBits(b.bias) ||
						!sameBits(a.weights,b.weights) || !Arrays.equals(a.inputs,b.inputs)){
					throw new IllegalStateException("Worker copy changed forward parameters at layer="+layer+" cell="+i);
				}
			}
		}
	}

	/** Raw bits include signed zero and NaN payloads; numeric equality is insufficient here. */
	private static boolean sameBits(float[] a,float[] b){
		if(a==b){return true;}
		if(a==null || b==null || a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(Float.floatToRawIntBits(a[i])!=Float.floatToRawIntBits(b[i])){return false;}
		}
		return true;
	}

	/** Separate oracle and native-pair errors; no warmed-up run replaces a failed run. */
	static final class Comparison{
		void print(String label){
			System.out.println(label+"_max_absolute_difference\t"+maximum);
			System.out.println(label+"_reference_oracle_max_absolute_difference\t"+referenceMaximum);
			System.out.println(label+"_native_max_absolute_difference\t"+nativeMaximum);
			System.out.println(label+"_first_maximum_changes\t"+oracleArgmaxChanges);
			System.out.println(label+"_reference_oracle_first_maximum_changes\t"+referenceOracleArgmaxChanges);
			System.out.println(label+"_native_first_maximum_changes\t"+nativeArgmaxChanges);
		}
		double maximum,referenceMaximum,nativeMaximum;
		long oracleArgmaxChanges,referenceOracleArgmaxChanges,nativeArgmaxChanges;
	}

	/** Skips only empty/comment lines; data order must remain unchanged. */
	private static byte[] nextData(ByteFile bf){
		assert(bf!=null) : "Prediction comparison must open its input before reading";
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]!='#'){return line;}
		}
		return null;
	}
}
