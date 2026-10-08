package ml;

import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import parse.LineParser1;
import parse.PreParser;
import prot.MagQCNetBundle;
import shared.Shared;

/**
 * Measures native sparse and dense inference on one real subnet and input row.
 * The model and row are loaded before timing; each mode warms up separately and
 * runs for at least ten seconds. Reports numerical differences without treating
 * a speed measurement as model-quality or full-pipeline acceptance.
 * @author Yoimiya
 */
public final class CellNetInferenceBenchmark {

	/** Requires a bundle/subnet, its original subnet vectors, and an optional duration. */
	public static void main(String[] args) throws Exception{
		String bundle=null,subnet=null,input=null;
		int seconds=10;
		boolean selftest=false;
		for(String arg:new PreParser(args,CellNetInferenceBenchmark.class,false).args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("bundle")){bundle=value;}
			else if(key.equals("subnet")){subnet=value;}
			else if(key.equals("in")){input=value;}
			else if(key.equals("seconds")){seconds=Integer.parseInt(value);}
			else if(key.equals("selftest")){selftest=parse.Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+key);}
		}
		if(selftest){testDenseCopy(); return;}
		if(bundle==null || subnet==null || input==null || seconds<10){
			throw new IllegalArgumentException("Required: bundle= subnet= in=; seconds>=10");
		}
		final MagQCNetBundle loaded=MagQCNetBundle.loadMultiOutput(Paths.get(bundle));
		final MagQCNetBundle.Subnet selected=loaded.subnet(subnet);
		if(selected==null){throw new IllegalArgumentException("Subnet is absent: "+subnet);}
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(String line:new String(selected.netBytes,StandardCharsets.UTF_8).split("\\r?\\n")){
			lines.add(line.getBytes(StandardCharsets.UTF_8));
		}
		final CellNet source=CellNetParser.loadFromLines(lines);
		final boolean sourceDense=CellNet.DENSE;
		final CellNet dense=source.copyDenseForInference();
		final float[] row=firstInput(input,source.numInputs());
		System.out.println("subnet\t"+subnet+"\ndims\t"+Arrays.toString(source.dims));
		System.out.println("source_dense\t"+sourceDense+"\ndensity\t"+source.density+
			"\ndensity1\t"+source.density1+"\nblock_size\t"+source.edgeBlockSize);
		System.out.println("simd\t"+Shared.SIMD+"\nsimd_fma\t"+Shared.SIMD_FMA+
			"\nsimd_feed_forward\t"+Shared.SIMD_FEED_FORWARD+"\nsparse_simd_fma\t"+simd.Vector.SIMD_FMA_SPARSE);
		System.out.println("source_simd_ff\t"+source.simdFF+"\ndense_simd_ff\t"+dense.simdFF);
		source.applyInput(row); explicitForward(source,sourceDense);
		dense.applyInput(row); dense.feedForward();
		reportOutputs("cold",source,dense);
		time("source",source,row,sourceDense,seconds);
		time("dense",dense,row,true,seconds);
		// Timed loops leave each network's final outputs populated. Report these too:
		// JIT warmup can change SIMD reduction behavior relative to its first invocation.
		reportOutputs("warm",source,dense);
	}

	/** Reports all finite outputs before/after warmup without inventing an acceptance tolerance. */
	private static void reportOutputs(String phase,CellNet source,CellNet dense){
		assert(source.numOutputs()==dense.numOutputs()) : "Dense inference copies preserve output width";
		double maxDifference=0;
		for(int i=0; i<source.numOutputs(); i++){
			final float a=source.getOutput(i),b=dense.getOutput(i);
			if(!Float.isFinite(a) || !Float.isFinite(b)){throw new IllegalStateException("Nonfinite inference output");}
			maxDifference=Math.max(maxDifference,Math.abs((double)a-b));
			System.out.println(phase+"_output_"+i+"\t"+a+"\t"+b);
		}
		System.out.println(phase+"_max_abs_difference\t"+maxDifference);
	}

	/** Reads only the first vector, checking its declared input width before scoring. */
	private static float[] firstInput(String path,int width){
		assert(width>0) : "A network input must have a positive width";
		final ByteFile file=ByteFile.makeByteFile1(path,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		boolean dimensions=false;
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(line.length==0){continue;}
				parser.set(line);
				if(line[0]=='#'){
					if(parser.termEquals("#dims",0)){
						if(parser.terms()!=4 || parser.parseInt(1)!=width){throw new IllegalArgumentException("Subnet vector width mismatch");}
						dimensions=true;
					}
					continue;
				}
				if(!dimensions || parser.terms()<width){throw new IllegalArgumentException("Missing dimensions or truncated vector row");}
				final float[] row=new float[width];
				for(int i=0; i<width; i++){
					row[i]=parser.parseFloat(i);
					if(!Float.isFinite(row[i])){throw new IllegalArgumentException("Nonfinite input feature "+i);}
				}
				return row;
			}
		}finally{if(file.close()){throw new IllegalStateException("Could not close subnet vector input");}}
		throw new IllegalArgumentException("No vector rows: "+path);
	}

	/** Uses explicit native dispatch so another parser's global flag cannot affect the comparison. */
	private static void explicitForward(CellNet net,boolean dense){
		assert(net!=null) : "Benchmark dispatch requires a fully loaded model";
		if(dense){net.feedForwardDense();}else{net.feedForwardSparse();}
	}

	/** Times repeated scoring, with a consumed checksum and a minimum ten-second compute window. */
	private static void time(String label,CellNet net,float[] row,boolean dense,int seconds){
		assert(seconds>=10 && row.length==net.numInputs()) : "Timing excludes short startup-dominated runs and wrong-width input";
		for(int i=0; i<2000; i++){net.applyInput(row); explicitForward(net,dense);}
		long operations=0;
		for(int layer=1; layer<net.layers; layer++){
			for(Cell cell:net.net[layer]){operations+=cell.weights.length;}
		}
		long count=0;
		double sum=0;
		final long start=System.nanoTime(),minimum=seconds*1000000000L;
		do{
			for(int i=0; i<256; i++){
				net.applyInput(row); explicitForward(net,dense); sum+=net.getOutput(0);
			}
			count+=256;
		}while(System.nanoTime()-start<minimum);
		final double elapsed=(System.nanoTime()-start)/1e9;
		checksum=sum;
		System.out.println(label+"\trows="+count+"\tseconds="+elapsed+"\trows_per_second="+(count/elapsed)+
			"\tmacs_per_row="+operations+"\tmacs_per_second="+(operations*(double)count/elapsed)+"\tchecksum="+checksum);
	}

	/** Checks sparse index expansion, dense input, private arrays and dispatch without timing. */
	private static void testDenseCopy(){
		final boolean previous=CellNet.DENSE;
		try{
			final CellNet source=new CellNet(new int[]{3,2,1},1,0.5f,0,1,new ArrayList<String>());
			set(source.net[1][0],new int[]{2,0},new float[]{2,3},1);
			set(source.net[1][1],new int[]{1},new float[]{-1},2);
			set(source.net[2][0],new int[]{0,1},new float[]{2,3},-1);
			final CellNet dense=source.copyDenseForInference(),second=dense.copyDenseForInference();
			if(CellNet.DENSE!=previous){throw new AssertionError("Dense copy changed global dispatch");}
			final float[] row={1,2,3}; // Hidden: 10,0. Output: 19.
			for(boolean global:new boolean[]{false,true}){
				CellNet.DENSE=global;
				dense.applyInput(row); second.applyInput(row);
				if(dense.feedForward()!=19 || second.feedForward()!=19){throw new AssertionError("Dense copy changed indexed-edge arithmetic");}
			}
			dense.net[1][0].weights[0]=100;
			if(source.net[1][0].weights[1]!=3 || second.net[1][0].weights[0]!=3){throw new AssertionError("Inference weights are shared between copies");}
			System.out.println("CellNetInferenceBenchmark SELFTEST PASS");
		}finally{CellNet.DENSE=previous;}
	}

	/** Constructs one explicit sparse fixture cell with a linear activation. */
	private static void set(Cell cell,int[] indices,float[] weights,float bias){
		assert(indices.length==weights.length) : "Fixture edges must pair each input with its trained weight";
		cell.inputs=indices; cell.weights=weights; cell.bias=bias;
		cell.function=Function.getFunction(Function.LINEAR);
	}

	private static volatile double checksum;
}
