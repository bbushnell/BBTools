package ml;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fileIO.ByteStreamWriter;
import parse.Parser;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Paired native inference on pinned C-order float32 NPY validation rows.
 * Measures the first two value outputs, never error heads or training loss.
 * @author Ganyu
 */
public final class CellNetWeightBitsEval {

	public static void main(String[] args) throws IOException{
		// PreParser opens outstream before validation; this read-only gate rejects it instead.
		args=Parser.parseConfig(args);
		String reference=null, candidate=null, x=null, y=null, out=null, storage="bfloat16";
		int rows=-1, limit=Integer.MAX_VALUE, threads=1;
		double originalMse=Double.NaN, absoluteLimit=Double.NaN;
		boolean explicitLimit=false;
		for(String arg:args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("reference")){reference=value;}
			else if(key.equals("candidate")){candidate=value;}
			else if(key.equals("x")){x=value;}
			else if(key.equals("y")){y=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("store")){storage=value;}
			else if(key.equals("rows")){rows=Integer.parseInt(value);}
			else if(key.equals("limit")){limit=Integer.parseInt(value);}
			else if(key.equals("t") || key.equals("threads")){threads=Integer.parseInt(value);}
			else if(key.equals("originalmse")){originalMse=Double.parseDouble(value);}
			else if(key.equals("limitmse")){absoluteLimit=Double.parseDouble(value); explicitLimit=true;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(reference==null || candidate==null || x==null || y==null || out==null || rows<1 || limit<1
				|| threads<1 || threads>128 || !Double.isFinite(originalMse) || originalMse<=0
				|| (explicitLimit && (!Double.isFinite(absoluteLimit) || absoluteLimit<0))
				|| !(storage.equals("bfloat16") || storage.equals("float32"))){
			throw new IllegalArgumentException("Require reference= candidate= x= y= out= rows= originalmse=; optional limitmse= t= limit= store=bfloat16|float32");
		}
		final double limitMse=explicitLimit ? absoluteLimit : originalMse*1.005;
		if(!Double.isFinite(limitMse)){throw new IllegalArgumentException("MSE limit overflow");}
		if(Files.exists(Paths.get(out))){throw new IllegalArgumentException("Fresh output required: "+out);}
		final long started=System.nanoTime();
		final CellNet ref=CellNetParser.load(reference,false);
		final boolean dense=CellNet.DENSE;
		final CellNet cand=CellNetParser.load(candidate,false);
		if(CellNet.DENSE!=dense){throw new IllegalArgumentException("Paired nets must use the same native dense/sparse representation");}
		checkPair(ref,cand);
		try(Npy rx=new Npy(x); Npy ry=new Npy(y)){
			if(rx.rows!=rows || ry.rows!=rows || rx.cols!=ref.numInputs() || ry.cols!=2){
				throw new IllegalArgumentException("Pinned row/input/label shape mismatch: X="+rx.rows+"x"+rx.cols+", Y="+ry.rows+"x"+ry.cols);
			}
		}
		final int count=Math.min(rows,limit), jobs=Math.min(threads,count);
		final ArrayList<Worker> workers=new ArrayList<Worker>(jobs);
		for(int i=0; i<jobs; i++){
			workers.add(new Worker(ref,cand,x,y,(int)((long)count*i/jobs),(int)((long)count*(i+1)/jobs),storage.equals("bfloat16")));
		}
		ThreadWaiter.startAndWait(workers);
		final double[] refSq=new double[2], candSq=new double[2], deltaSq=new double[2];
		double maxDelta=0;
		long seen=0;
		for(Worker worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Paired evaluation worker failed; no result published",worker.failure);}
			seen+=worker.seen; maxDelta=Math.max(maxDelta,worker.maxDelta);
			for(int k=0; k<2; k++){refSq[k]+=worker.refSq[k]; candSq[k]+=worker.candSq[k]; deltaSq[k]+=worker.deltaSq[k];}
		}
		if(seen!=count){throw new IllegalStateException("Workers must conserve every requested row: "+seen+" != "+count);}
		final double refMse=valueMse(refSq[0]+refSq[1],seen), candMse=valueMse(candSq[0]+candSq[1],seen);
		final ByteBuilder bb=new ByteBuilder();
		bb.append("rows\tsource_rows\tstore\treference_mse\tcandidate_mse\tcandidate_over_reference\toriginal_mse\tlimit_mse\t")
			.append(explicitLimit ? "within_absolute_bar" : "within_original_bar")
			.append("\tprediction_delta_mse\tmax_abs_value_delta\telapsed_seconds\n");
		bb.append(seen).tab().append(rows).tab().append(storage).tab().append(Double.toString(refMse)).tab()
			.append(Double.toString(candMse)).tab().append(Double.toString(candMse/refMse)).tab()
			.append(Double.toString(originalMse)).tab().append(Double.toString(limitMse)).tab()
			.append(count==rows ? Boolean.toString(candMse<=limitMse) : "CANARY_ONLY").tab()
			.append(Double.toString(valueMse(deltaSq[0]+deltaSq[1],seen))).tab().append(Double.toString(maxDelta)).tab()
			.append(Double.toString((System.nanoTime()-started)*1e-9)).nl();
		final ByteStreamWriter writer=new ByteStreamWriter(out,false,false,false);
		writer.start(); writer.print(bb);
		if(writer.poisonAndWait()){throw new IOException("Could not write paired evaluation: "+out);}
		System.err.println("WEIGHT_BITS_EVAL_PASS rows="+seen+" sourceRows="+rows+" store="+storage);
	}

	/** Frozen train_mlp.evaluate_validation_files uses sum of two head SSEs / rows, not /2rows. */
	static double valueMse(double sumSquared,long rows){
		if(rows<1 || !Double.isFinite(sumSquared) || sumSquared<0){throw new IllegalArgumentException("Invalid value-loss accumulator");}
		return sumSquared/rows;
	}

	/** Independent bit oracle: truncate omitted bits while retaining a nonzero edge at the smallest stored magnitude. */
	static void checkPair(CellNet a, CellNet b){
		final int bits=b.weightBits();
		if(bits!=18 && bits!=24){throw new IllegalArgumentException("Candidate must declare18-bit or24-bit trained edges");}
		final int mask=bits==18 ? 0xffffc000 : 0xffffff00;
		if(a.tags.containsKey("affine_prefix_layers") || b.tags.containsKey("affine_prefix_layers")){
			throw new IllegalArgumentException("Convert explicit-affine prefixes to input headers before paired storage evaluation");
		}
		if(!Arrays.equals(a.dims,b.dims) || a.numOutputs()<2){throw new IllegalArgumentException("Paired topology/value-output mismatch");}
		for(int layer=1; layer<a.net.length; layer++){
			for(int i=0; i<a.net[layer].length; i++){
				final Cell ca=a.net[layer][i], cb=b.net[layer][i];
				if(!Arrays.equals(ca.inputs,cb.inputs) || ca.weights.length!=cb.weights.length
						|| Float.floatToRawIntBits(ca.bias())!=Float.floatToRawIntBits(cb.bias())
						|| !ca.function.name().equals(cb.function.name())){
					throw new IllegalArgumentException("Topology/bias/activation changed at layer="+layer+" cell="+i);
				}
				for(int j=0; j<ca.weights.length; j++){
					final int raw=Float.floatToRawIntBits(ca.weights[j]);
					int expected=raw&mask;
					if((raw&0x7fffffff)!=0 && (expected&0x7fffffff)==0){expected|=1<<(32-bits);}
					if(expected!=Float.floatToRawIntBits(cb.weights[j])){
						throw new IllegalArgumentException("Reduced edge oracle mismatch at layer="+layer+" cell="+i+" edge="+j);
					}
				}
			}
		}
	}

	/** Same RNE bfloat16 conversion as trainer storage, after float32 normalization. */
	static float stored(float x, boolean bf16){
		if(!Float.isFinite(x)){throw new IllegalArgumentException("Nonfinite normalized input");}
		if(!bf16){return x;}
		x=Math.max(-10000f,Math.min(10000f,x));
		final int bits=Float.floatToRawIntBits(x);
		return Float.intBitsToFloat((bits+0x7fff+((bits>>>16)&1))&0xffff0000);
	}

	/** Disjoint row ranges with private readers, nets and buffers; merge once after join. */
	private static final class Worker extends Thread {
		Worker(CellNet a,CellNet b,String xp,String yp,int lo,int hi,boolean bf){
			assert(lo>=0 && hi>lo) : "Worker ranges partition the positive requested row count";
			ref=a; cand=b; x=xp; y=yp; start=lo; end=hi; bf16=bf;
		}
		@Override public void run(){
			try(Npy rx=new Npy(x); Npy ry=new Npy(y)){
				final CellNet a=ref.copy(false), b=cand.copy(false);
				final float[] input=new float[rx.cols], labels=new float[2];
				rx.seek(start); ry.seek(start);
				for(int row=start; row<end; row++){
					rx.read(input); ry.read(labels);
					a.applyInput(input); b.applyInput(input);
					for(int j=0; j<input.length; j++){
						if(Float.floatToRawIntBits(a.values[0][j])!=Float.floatToRawIntBits(b.values[0][j])){
							throw new IllegalStateException("Paired input normalization differs at row="+row+" column="+j);
						}
						a.values[0][j]=stored(a.values[0][j],bf16); b.values[0][j]=a.values[0][j];
					}
					a.feedForward(); b.feedForward();
					for(int k=0; k<a.numOutputs(); k++){
						if(!Float.isFinite(a.getOutput(k)) || !Float.isFinite(b.getOutput(k))){throw new IllegalStateException("Nonfinite output at row="+row);}
						if(k<2){
							final double av=a.getOutput(k), bv=b.getOutput(k), label=labels[k];
							refSq[k]+=(av-label)*(av-label); candSq[k]+=(bv-label)*(bv-label);
							deltaSq[k]+=(av-bv)*(av-bv); maxDelta=Math.max(maxDelta,Math.abs(av-bv));
						}
					}
					seen++;
					if(seen%10000==0){System.err.println("PROGRESS start="+start+" rows="+seen+" end="+end);}
				}
			}catch(Throwable t){failure=t;}
		}
		final CellNet ref,cand;
		final String x,y;
		final int start,end;
		final boolean bf16;
		final double[] refSq=new double[2],candSq=new double[2],deltaSq=new double[2];
		long seen;
		double maxDelta;
		Throwable failure;
	}

	/** Strict NPY v1/v2 reader for only the pinned uncompressed C-order <f4 matrices. */
	static final class Npy implements AutoCloseable {
		Npy(String path) throws IOException{
			file=new RandomAccessFile(path,"r");
			try{
				final byte[] magic=new byte[6]; file.readFully(magic);
				if(!Arrays.equals(magic,new byte[]{(byte)0x93,'N','U','M','P','Y'})){throw new IOException("Not NPY: "+path);}
				final int major=file.readUnsignedByte(), minor=file.readUnsignedByte();
				if((major!=1 && major!=2) || minor!=0){throw new IOException("Unsupported NPY version");}
				final int size=major==1 ? Short.toUnsignedInt(Short.reverseBytes(file.readShort())) : Integer.reverseBytes(file.readInt());
				if(size<1 || size>1048576){throw new IOException("Invalid NPY header length");}
				final byte[] raw=new byte[size]; file.readFully(raw);
				final String header=new String(raw,StandardCharsets.US_ASCII);
				if(!header.matches("(?s).*'descr'\\s*:\\s*'<f4'.*") || !header.matches("(?s).*'fortran_order'\\s*:\\s*False.*")){
					throw new IOException("Require little-endian C-order float32 NPY: "+header);
				}
				final Matcher shape=Pattern.compile("'shape'\\s*:\\s*\\(\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,?\\s*\\)").matcher(header);
				if(!shape.find()){throw new IOException("Require two-dimensional NPY");}
				rows=Integer.parseInt(shape.group(1)); cols=Integer.parseInt(shape.group(2)); offset=file.getFilePointer();
				if(rows<1 || cols<1 || cols>1000000 || file.length()!=offset+4L*rows*cols){throw new IOException("NPY dimensions/payload length mismatch");}
				bytes=new byte[cols*4]; buffer=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
			}catch(Throwable t){file.close(); throw t;}
		}
		void seek(int row) throws IOException{
			if(row<0 || row>=rows){throw new IllegalArgumentException("NPY row outside matrix");}
			file.seek(offset+4L*row*cols);
		}
		void read(float[] values) throws IOException{
			assert(values.length==cols) : "NPY row buffer must match the declared matrix width";
			file.readFully(bytes); buffer.position(0);
			for(int i=0; i<cols; i++){
				values[i]=buffer.getFloat();
				if(!Float.isFinite(values[i])){throw new IOException("Nonfinite NPY value at column="+i);}
			}
		}
		/** Reads consecutive rows with one bounded bulk read, avoiding one system read per row. */
		void readBatch(float[] values, int rowCount) throws IOException{
			assert(cols==2) : "Batch path is currently for two-column truth matrices";
			if(batchBytes==null){
				final int batchRows=Math.min(rows, 65536/(cols*4));
				batchBytes=new byte[Math.max(1, batchRows)*cols*4];
				batchBuffer=ByteBuffer.wrap(batchBytes).order(ByteOrder.LITTLE_ENDIAN);
			}
			assert(rowCount>0 && rowCount<=batchBytes.length/(cols*4)) : "Invalid NPY batch row count";
			assert(values.length>=rowCount*cols) : "NPY batch buffer is too small";
			final int byteCount=rowCount*cols*4;
			file.readFully(batchBytes, 0, byteCount); batchBuffer.position(0).limit(byteCount);
			for(int i=0; i<rowCount*cols; i++){
				values[i]=batchBuffer.getFloat();
				if(!Float.isFinite(values[i])){throw new IOException("Nonfinite NPY value at batch index="+i);}
			}
		}
		@Override public void close() throws IOException{file.close();}
		final RandomAccessFile file;
		final int rows,cols;
		final long offset;
		final byte[] bytes;
		final ByteBuffer buffer;
		byte[] batchBytes;
		ByteBuffer batchBuffer;
	}
}
