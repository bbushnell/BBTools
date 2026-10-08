package ml;

import java.io.DataOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import structures.ByteBuilder;

/** Independent NPY/endian, BF16 tie and paired-loss fixtures.
 * @author Ganyu
 */
public final class CellNetWeightBitsEvalTest {

	public static void main(String[] args) throws Exception{
		boolean assertions=false; assert(assertions=true);
		if(!assertions || args.length!=0){throw new IllegalArgumentException("Fixture requires -ea and no arguments");}
		assert(CellNetWeightBitsEval.valueMse(2,1)==2) : "Frozen trainer convention: two unit value residuals yield MSE2, not1";
		affineFixture(false); affineFixture(true);
		assert(Float.floatToRawIntBits(CellNetWeightBitsEval.stored(Float.intBitsToFloat(0x3f808000),true))==0x3f800000) : "BF16 halfway-even must round down to even";
		assert(Float.floatToRawIntBits(CellNetWeightBitsEval.stored(Float.intBitsToFloat(0x3f818000),true))==0x3f820000) : "BF16 halfway-odd must round up to even";
		assert(Float.floatToRawIntBits(CellNetWeightBitsEval.stored(-0f,true))==0x80000000) : "Storage must preserve negative zero";
		assert(CellNetWeightBitsEval.stored(20000,true)==9984f) : "Trainer clamps to10000 before BF16 conversion";
		final Path root=Files.createTempDirectory("weightbits-eval-fixture-");
		floorFixture(root);
		final Path sentinel=root.resolve("sentinel.txt"), config=root.resolve("bad.config");
		Files.write(sentinel,"preserve".getBytes(StandardCharsets.US_ASCII));
		rejectArgs(new String[]{"outstream="+sentinel});
		Files.write(config,("outstream="+sentinel+"\n").getBytes(StandardCharsets.US_ASCII));
		rejectArgs(new String[]{"config="+config});
		assert(new String(Files.readAllBytes(sentinel),StandardCharsets.US_ASCII).equals("preserve")) : "Rejected CLI/config must not truncate a pre-existing outstream";
		final float[][] x={{1,2},{-1,3},{2,-2},{0,1}}, y={{0,0},{0,0},{0,0},{0,0}};
		final Path xp=root.resolve("x.npy"), yp=root.resolve("y.npy");
		write(xp,x); write(yp,y);
		try(CellNetWeightBitsEval.Npy input=new CellNetWeightBitsEval.Npy(xp.toString())){
			final float[] row=new float[2]; input.seek(2); input.read(row);
			assert(row[0]==2 && row[1]==-2) : "NPY reader must preserve little-endian floats and random row offsets";
		}
		final float weight=Float.intBitsToFloat(0x3f812345);
		final ByteBuilder text=new ByteBuilder();
		text.append("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 2 2\n#coding A48\n");
		text.append("C3 LINEAR ").appendFloatA48(0f).space().appendFloatA48(weight).space().appendFloatA48(0f).nl();
		text.append("C4 LINEAR ").appendFloatA48(0f).space().appendFloatA48(0f).space().appendFloatA48(weight).nl();
		final Path ref=root.resolve("reference.bbnet");
		Files.write(ref,text.toBytes());
		for(int precision:new int[]{24,18}){
			final Path cand=root.resolve("candidate"+precision+".bbnet");
			CellNetWeightBits.main(new String[]{"in="+ref,"out="+cand,"bits="+precision});
			final float truncated=Float.intBitsToFloat(Float.floatToRawIntBits(weight)&(precision==18 ? 0xffffc000 : 0xffffff00));
			double expected=0;
			for(float[] row:x){for(float value:row){final float prediction=value*truncated; final double p=prediction; expected+=p*p;}}
			expected/=4;
			for(int threads:new int[]{1,2}){
				final Path out=root.resolve("result"+precision+"_"+threads+".tsv");
				CellNetWeightBitsEval.main(new String[]{"reference="+ref,"candidate="+cand,"x="+xp,"y="+yp,"out="+out,
					"rows=4","t="+threads,"store=float32","originalmse=8"});
				final String[] fields=Files.readAllLines(out).get(1).split("\t");
				assert(fields[0].equals("4")) : "Parallel ranges must conserve all four rows";
				assert(Math.abs(Double.parseDouble(fields[4])-expected)<1e-12) : "MSE sums BOTH value-output squared errors and divides only by rows, matching the fixed trainer bar";
				assert(fields[8].equals("true")) : "Original fixed error bar must govern the result";
				assert(Double.parseDouble(fields[6])==8 && Double.parseDouble(fields[7])==8*1.005
					&& Files.readAllLines(out).get(0).contains("within_original_bar")) : "Default bar and schema remain unchanged";
			}
			absoluteLimitFixture(root, ref, cand, xp, yp, expected, precision);
			final CellNet original=CellNetParser.load(ref.toString(),false), altered=CellNetParser.load(cand.toString(),false);
			altered.net[1][0].weights[0]=Float.intBitsToFloat(Float.floatToRawIntBits(altered.net[1][0].weights[0])|1);
			boolean badEdge=false;
			try{CellNetWeightBitsEval.checkPair(original,altered);}catch(IllegalArgumentException expectedFailure){badEdge=true;}
			assert(badEdge) : "A candidate claiming reduced precision must not retain an omitted low bit";
		}
		final byte[] truncatedFile=Files.readAllBytes(xp);
		final Path broken=root.resolve("broken.npy"); Files.write(broken,java.util.Arrays.copyOf(truncatedFile,truncatedFile.length-1));
		boolean rejected=false;
		try(CellNetWeightBitsEval.Npy ignored=new CellNetWeightBitsEval.Npy(broken.toString())){throw new AssertionError("Short payload accepted: rows="+ignored.rows);}
		catch(java.io.IOException expectedFailure){rejected=true;}
		assert(rejected) : "NPY file-size gate must reject truncated matrices";
		try(java.nio.file.DirectoryStream<Path> files=Files.newDirectoryStream(root)){for(Path file:files){Files.delete(file);}}
		Files.delete(root);
		System.out.println("CellNetWeightBitsEvalTest PASS BF16_RNE_clamp NPY_seek_endian truncated_payload sum2_SSE_per_row threads1_2 CLI_config_preservation precision18_24 absolute_limit");
	}

	/** Positive and negative subnormal edges survive reduced storage; the evaluator must accept the writer's floor. */
	private static void floorFixture(final Path root) throws Exception{
		for(int precision:new int[]{18, 24}){
			final int quantum=1<<(32-precision);
			final int[] raw={0, 0x80000000, 1, 0x80000001, quantum-1, 0x80000000|(quantum-1), quantum, 0x80000000|quantum};
			final ByteBuilder text=new ByteBuilder();
			text.append("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 8 2\n#coding A48\n");
			for(int output=0; output<2; output++){
				text.append('C').append(9+output).append(" LINEAR ").appendFloatA48(0);
				for(int bits:raw){text.space().appendFloatA48(Float.intBitsToFloat(bits));}
				text.nl();
			}
			final Path input=root.resolve("floor"+precision+".bbnet"), output=root.resolve("floor"+precision+"_reduced.bbnet");
			Files.write(input, text.toBytes());
			CellNetWeightBits.main(new String[]{"in="+input, "out="+output, "bits="+precision});
			final CellNet original=CellNetParser.load(input.toString(), false), reduced=CellNetParser.load(output.toString(), false);
			for(int i=0; i<raw.length; i++){
				final int expected=(raw[i]&0x7fffffff)==0 ? raw[i] : (raw[i]&0x80000000)|quantum;
				assert(Float.floatToRawIntBits(reduced.net[1][0].weights[i])==expected) :
					"The smallest stored nonzero magnitude preserves sign at18/24-bit boundaries";
			}
			CellNetWeightBitsEval.checkPair(original, reduced);
			reduced.net[1][0].weights[2]=0;
			boolean rejected=false;
			try{CellNetWeightBitsEval.checkPair(original, reduced);}catch(IllegalArgumentException expected){rejected=true;}
			assert(rejected) : "The paired oracle must reject a nonzero reference edge erased by truncation";
		}
	}

	/** CLI absolute limits change acceptance without changing measured losses or original provenance. */
	private static void absoluteLimitFixture(Path root, Path ref, Path cand, Path x, Path y,
			double expectedMse, int precision) throws Exception{
		final double[] bars={1, expectedMse, 10};
		for(int i=0; i<bars.length; i++){
			final Path out=root.resolve("absolute"+precision+"_"+i+".tsv");
			CellNetWeightBitsEval.main(new String[]{"reference="+ref,"candidate="+cand,"x="+x,"y="+y,"out="+out,
				"rows=4","t=1","store=float32","originalmse=8","limitmse="+bars[i]});
			final String[] fields=Files.readAllLines(out).get(1).split("\t");
			assert(Double.parseDouble(fields[4])==expectedMse && Double.parseDouble(fields[6])==8
				&& Double.parseDouble(fields[7])==bars[i]) : "Absolute bar must preserve actual loss and original MSE";
			assert(fields[8].equals(i==0 ? "false" : "true")) : "Strict bar fails; inclusive equality and larger bar pass";
			assert(Files.readAllLines(out).get(0).contains("within_absolute_bar")) : "Do not mislabel an absolute verdict as the original relative bar";
		}
		final Path partial=root.resolve("absolute_partial"+precision+".tsv");
		CellNetWeightBitsEval.main(new String[]{"reference="+ref,"candidate="+cand,"x="+x,"y="+y,"out="+partial,
			"rows=4","limit=1","store=float32","originalmse=8","limitmse=100"});
		assert(Files.readAllLines(partial).get(1).split("\t")[8].equals("CANARY_ONLY")) : "A partial validation cannot qualify under any limit";
		for(String bad:new String[]{"NaN","Infinity","-1"}){
			final Path out=root.resolve("bad_limit"+precision+bad+".tsv");
			rejectArgs(new String[]{"reference="+ref,"candidate="+cand,"x="+x,"y="+y,"out="+out,
				"rows=4","originalmse=8","limitmse="+bad});
			assert(!Files.exists(out)) : "Invalid limits must be rejected before output creation";
		}
	}

	/** Malformed options must fail before opening any output. */
	private static void rejectArgs(String[] args) throws Exception{
		try{CellNetWeightBitsEval.main(args);}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid evaluation arguments accepted");
	}

	/** Diagonal-prefix conversion must preserve actual float32 predictions and learned bits. */
	private static void affineFixture(boolean sparse){
		final String text="##bbnet\n##normalization explicit-affine\n##affine_prefix_layers 2\n#version 1\n#concise\n"+(sparse ? "#sparse\n" : "#dense\n")
			+"#density 1\n#blocksize 1\n#seed 1\n#layers 4\n#dims 2 2 2 2\n#coding decimal\n"
			+(sparse ? "I3 0\nW3 LINEAR -2 1\nI4 1\nW4 LINEAR 3 1\nI5 0\nW5 LINEAR 0 2\nI6 1\nW6 LINEAR 0 0.5\nI7 0 1\nW7 LINEAR 0.25 0.125 0.375\nI8 0 1\nW8 LINEAR -0.5 0.5 -0.25\n"
			: "C3 LINEAR -2 1 0\nC4 LINEAR 3 0 1\nC5 LINEAR 0 2 0\nC6 LINEAR 0 0 0.5\nC7 LINEAR 0.25 0.125 0.375\nC8 LINEAR -0.5 0.5 -0.25\n");
		final java.util.ArrayList<byte[]> lines=new java.util.ArrayList<byte[]>();
		for(String line:text.split("\n")){lines.add(line.getBytes(StandardCharsets.US_ASCII));}
		final CellNet original=CellNetParser.loadFromLines(lines), converted=CellNetAffineInput.convert(original);
		CellNet.codingA48Out=true; CellNet.OUT_DENSE=false; CellNet.OUT_SPARSE=false;
		lines.clear();
		for(String line:converted.toBytes().toString().split("\n")){lines.add(line.getBytes(StandardCharsets.US_ASCII));}
		final CellNet reloaded=CellNetParser.loadFromLines(lines);
		assert(converted.hasInputNormalization() && converted.dims.length==2) : "Exactly two affine layers become input headers";
		for(float[] row:new float[][]{{0,0},{1,2},{-3,4},{Float.MIN_VALUE,-0f}}){
			original.applyInput(row); original.feedForward(); converted.applyInput(row); converted.feedForward(); reloaded.applyInput(row); reloaded.feedForward();
			for(int k=0; k<2; k++){
				assert(Float.floatToRawIntBits(original.getOutput(k))==Float.floatToRawIntBits(converted.getOutput(k))) : "Affine normalization conversion must be exact before truncation";
				assert(Float.floatToRawIntBits(original.getOutput(k))==Float.floatToRawIntBits(reloaded.getOutput(k))) : "Serialized affine conversion must preserve all output bits in dense and sparse modes";
			}
		}
		for(int k=0; k<2; k++){
			assert(java.util.Arrays.equals(original.net[3][k].weights,converted.net[1][k].weights)) : "Trained weights remain untouched";
		}
		original.tags.put("affine_prefix_layers","3");
		boolean rejected=false;
		try{CellNetAffineInput.convert(original);}catch(IllegalArgumentException expected){rejected=true;}
		assert(rejected) : "Undocumented affine structures must not be guessed";
	}

	/** Literal little-endian float payload with a standard NPYv1 header. */
	private static void write(Path path,float[][] rows) throws Exception{
		assert(rows.length>0 && rows[0].length>0) : "Fixture matrices must be nonempty";
		String header="{'descr': '<f4', 'fortran_order': False, 'shape': ("+rows.length+", "+rows[0].length+"), }";
		while((10+header.length()+1)%64!=0){header+=" ";} header+="\n";
		try(DataOutputStream out=new DataOutputStream(Files.newOutputStream(path))){
			out.write(new byte[]{(byte)0x93,'N','U','M','P','Y',1,0}); out.writeShort(Short.reverseBytes((short)header.length()));
			out.write(header.getBytes(StandardCharsets.US_ASCII));
			for(float[] row:rows){for(float value:row){out.writeInt(Integer.reverseBytes(Float.floatToRawIntBits(value)));}}
		}
	}
}
