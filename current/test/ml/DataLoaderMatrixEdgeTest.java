package test.ml;

import ml.DataLoader;
import ml.Matrix;
import ml.SampleSet;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Stream;

/**
 * Selectable -ea cases for DataLoader and Matrix edge contracts. Cases mutate
 * Matrix globals, so execute them serially; run(name) lets an external harness
 * isolate cases, while no-arg main stops at the first failure.
 * @author Brian Bushnell, Sucrose
 */
public class DataLoaderMatrixEdgeTest{
	
	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("shuffle_truncates");
		run("shuffle_no_truncate");
		run("inclusive_all_rows_validation");
		run("inclusive_partial_validation");
		run("convert_to_01_finite_mean");
		run("detect_range_finite_mean");
		run("target_2_3_equal_width");
		run("target_neg1_0_equal_width");
		run("target_0_1_unchanged");
		run("target_2_4_rescaled");
		run("balanced_target_0_2");
		run("balanced_target_2_3");
		run("balanced_target_0_1");
		run("unbalanced_target_0_2");
		run("balanced_binarize");
		run("unbalanced_binarize");
	}
	
	/** Runs one case in a private temporary directory and restores Matrix globals. */
	private static void run(String name) throws Exception{
		final boolean oldConvertTo01=Matrix.convertTo01;
		final float oldTargetMin=Matrix.targetOutputRangeMin;
		final float oldTargetMax=Matrix.targetOutputRangeMax;
		final boolean oldSetMin=Matrix.setTargetOutputRangeMin;
		final boolean oldSetMax=Matrix.setTargetOutputRangeMax;
		final Path dir=Files.createTempDirectory("dataloader-matrix-edge-");
		try{
			if(name.equals("shuffle_truncates")){
				checkShuffle(load(writeSingleOutput(dir), 2, true, 0f, 0, true), 2);
			}else if(name.equals("shuffle_no_truncate")){
				checkShuffle(load(writeSingleOutput(dir), 3, true, 0f, 0, true), 3);
			}else if(name.equals("inclusive_all_rows_validation")){
				checkSplit(load(writeSingleOutput(dir), 3, false, 0.1f, 3, false), 3, 3);
			}else if(name.equals("inclusive_partial_validation")){
				checkSplit(load(writeSingleOutput(dir), 3, false, 0.1f, 2, false), 3, 2);
			}else if(name.equals("convert_to_01_finite_mean")){
				Matrix.convertTo01=true;
				checkMatrix(load(writeMultiOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{0f, 1f}, {0f, 1f}}, "mean=0.5\n");
			}else if(name.equals("detect_range_finite_mean")){
				checkMatrix(load(writeMultiOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{0.25f, 0.75f}, {0.125f, 0.875f}}, "mean=0.5\n");
			}else if(name.equals("target_2_3_equal_width")){
				setTargetOutputRange(2f, 3f);
				checkMatrix(load(writeZeroOneOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{2f, 3f}, {3f, 2f}}, "mean=2.5\n",
					"midpoint=2.5\n");
			}else if(name.equals("target_neg1_0_equal_width")){
				setTargetOutputRange(-1f, 0f);
				checkMatrix(load(writeZeroOneOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{-1f, 0f}, {0f, -1f}}, "mean=-0.5\n",
					"midpoint=-0.5\n");
			}else if(name.equals("target_0_1_unchanged")){
				setTargetOutputRange(0f, 1f);
				checkMatrix(load(writeZeroOneOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{0f, 1f}, {1f, 0f}}, "mean=0.5\n");
			}else if(name.equals("target_2_4_rescaled")){
				setTargetOutputRange(2f, 4f);
				checkMatrix(load(writeZeroOneOutput(dir), 2, false, 0f, 0, true)[0].matrix,
					new float[][]{{2f, 4f}, {4f, 2f}}, "mean=3.0\n",
					"midpoint=3.0\n");
			}else if(name.equals("balanced_target_0_2")){
				checkBalancedRange(dir, 0f, 2f, 1f, 6, 1f);
			}else if(name.equals("balanced_target_2_3")){
				checkBalancedRange(dir, 2f, 3f, 1f, 6, 2.5f);
			}else if(name.equals("balanced_target_0_1")){
				checkBalancedRange(dir, 0f, 1f, 1f, 6, 0.5f);
			}else if(name.equals("unbalanced_target_0_2")){
				checkBalancedRange(dir, 0f, 2f, 0f, 4, 0.5f);
			}else if(name.equals("balanced_binarize")){
				checkBalancedBinarize(dir, 1f, 6, 0.5f);
			}else if(name.equals("unbalanced_binarize")){
				checkBalancedBinarize(dir, 0f, 4, 0.75f);
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("DataLoaderMatrixEdgeTest PASS "+name);
		}finally{
			Matrix.convertTo01=oldConvertTo01;
			Matrix.targetOutputRangeMin=oldTargetMin;
			Matrix.targetOutputRangeMax=oldTargetMax;
			Matrix.setTargetOutputRangeMin=oldSetMin;
			Matrix.setTargetOutputRangeMax=oldSetMax;
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}
	
	/** Calls the public loader route used by Trainer. */
	private static SampleSet[] load(Path path, int maxLines0, boolean shuffleRaw,
			float splitFraction, int maxLines1, boolean exclusive){
		return DataLoader.load(path.toString(), maxLines0, shuffleRaw, splitFraction,
			maxLines1, exclusive, 0f);
	}
	
	/** Three distinct one-output rows expose raw shuffle truncation and inclusive splits. */
	private static Path writeSingleOutput(Path dir) throws Exception{
		final Path path=dir.resolve("single-output.tsv");
		Files.write(path, ("#dims\t1\t1\n"+
			"##input\ttarget\n"+
			"0\t0.1\n"+
			"1\t0.4\n"+
			"2\t0.9\n").getBytes(StandardCharsets.US_ASCII));
		return path;
	}
	
	/** Two rows and two outputs make the binarized mean require four counted values. */
	private static Path writeMultiOutput(Path dir) throws Exception{
		final Path path=dir.resolve("multi-output.tsv");
		Files.write(path, ("#dims\t1\t2\n"+
			"##input\ttarget0\ttarget1\n"+
			"0\t0.25\t0.75\n"+
			"1\t0.125\t0.875\n").getBytes(StandardCharsets.US_ASCII));
		return path;
	}

	/** Two mirrored rows span exactly 0..1 in each output column. */
	private static Path writeZeroOneOutput(Path dir) throws Exception{
		final Path path=dir.resolve("zero-one-output.tsv");
		Files.write(path, ("#dims\t1\t2\n"+
			"##input\ttarget0\ttarget1\n"+
			"0\t0\t1\n"+
			"1\t1\t0\n").getBytes(StandardCharsets.US_ASCII));
		return path;
	}

	/** Requests a concrete target output range through Matrix's global flags. */
	private static void setTargetOutputRange(float min, float max){
		Matrix.targetOutputRangeMin=min;
		Matrix.targetOutputRangeMax=max;
		Matrix.setTargetOutputRangeMin=true;
		Matrix.setTargetOutputRangeMax=true;
	}

	/** Checks one affine transform per logical target even when balancing repeats its row. */
	private static void checkBalancedRange(Path dir, float min, float max, float balance,
			int expectedRows, float expectedMean) throws Exception{
		final Path path=dir.resolve("balanced-range.tsv");
		Files.write(path, ("#dims\t1\t1\n##input\ttarget\n"
			+"0\t0\n1\t0\n2\t0\n3\t1\n").getBytes(StandardCharsets.US_ASCII));
		Matrix.convertTo01=false;
		setTargetOutputRange(min, max);
		final Matrix m=DataLoader.load(path.toString(), 4, false, 0f, 0, true, balance)[0].matrix;
		assert(m.inputs.length==expectedRows) : "Balancing must retain all original rows and add only minority entries";
		final int[] counts=new int[4];
		for(int i=0; i<m.inputs.length; i++){
			final float input=m.inputs[i][0];
			final int id=(int)input;
			assert(input==id && id>=0 && id<counts.length) : "Unexpected balanced input ID "+input;
			counts[id]++;
			final float expected=id==3 ? max : min;
			assert(m.outputs[i].length==1 && m.outputs[i][0]==expected) :
				"Balanced range row "+id+": expected "+expected+", observed "+m.outputs[i][0]
				+"; each original target needs exactly one affine transform";
			assert(m.weights[i].length==1 && m.weights[i][0]==1f) : "Balancing must preserve unit row weights";
		}
		assert(counts[0]==1 && counts[1]==1 && counts[2]==1 && counts[3]==expectedRows-3) :
			"Only the positive minority row should be repeated";
		final String summary=m.toString();
		assert(summary.contains("mean="+expectedMean+"\n")) : "Cached mean must match literal transformed targets: "+summary;
		assert(m.outputMidpoint()==(min+max)*0.5f) : "Cached midpoint must describe the requested range";
	}
	
	/** Uses a negative midpoint so revisiting a binarized minority row would change 0 to 1. */
	private static void checkBalancedBinarize(Path dir, float balance, int expectedRows,
			float expectedMean) throws Exception{
		final Path path=dir.resolve("balanced-binarize.tsv");
		Files.write(path, ("#dims\t1\t1\n##input\ttarget\n"
			+"0\t-4\n1\t1\n2\t1\n3\t1\n").getBytes(StandardCharsets.US_ASCII));
		Matrix.convertTo01=true;
		Matrix.setTargetOutputRangeMin=Matrix.setTargetOutputRangeMax=false;
		final Matrix m=DataLoader.load(path.toString(), 4, false, 0f, 0, true, balance)[0].matrix;
		assert(m.inputs.length==expectedRows) : "Balancing must add only minority entries before binarization";
		final int[] counts=new int[4];
		for(int i=0; i<m.inputs.length; i++){
			final float input=m.inputs[i][0];
			final int id=(int)input;
			assert(input==id && id>=0 && id<counts.length) : "Unexpected balanced input ID "+input;
			counts[id]++;
			final float expected=id==0 ? 0f : 1f;
			assert(m.outputs[i].length==1 && m.outputs[i][0]==expected) :
				"Binarized balanced row "+id+": expected "+expected+", observed "+m.outputs[i][0]
				+"; compare the original target once against midpoint -1.5";
			assert(m.weights[i].length==1 && m.weights[i][0]==1f) : "Balancing must preserve unit row weights";
		}
		assert(counts[0]==expectedRows-3 && counts[1]==1 && counts[2]==1 && counts[3]==1) :
			"Only the negative minority row should be repeated";
		final String summary=m.toString();
		assert(summary.contains("mean="+expectedMean+"\n")) : "Cached mean must match literal binary targets: "+summary;
		assert(m.outputMidpoint()==0.5f) : "Binarization must refresh its nominal midpoint";
	}

	/** Verifies the requested row count without depending on deterministic shuffle order. */
	private static void checkShuffle(SampleSet[] sets, int rows){
		assert(sets!=null && sets.length==1) : "Raw shuffle without split must return one set";
		final Matrix m=sets[0].matrix;
		assert(m.inputs.length==rows) : "Loaded "+m.inputs.length+" rows, expected "+rows;
		final boolean[] seen=new boolean[3];
		for(int i=0; i<m.inputs.length; i++){
			final float in=m.inputs[i][0];
			final int index=(int)in;
			assert(in==index && index>=0 && index<seen.length) : "Unexpected retained input "+in;
			assert(!seen[index]) : "Duplicate retained input "+in;
			seen[index]=true;
			assert(m.outputs[i][0]==expectedSingleTarget(index)) :
				"Input/target pair drifted for input "+index+": "+m.outputs[i][0];
			assert(m.weights[i].length==1 && m.weights[i][0]==1f) :
				"Unweighted TSV should retain unit weight for row "+index;
		}
	}
	
	/** Inclusive splitting keeps the full original matrix as training and adds a validation subset. */
	private static void checkSplit(SampleSet[] sets, int trainRows, int validationRows){
		assert(sets!=null && sets.length==2) : "Validation split must return two sets";
		assert(sets[0].matrix.inputs.length==trainRows) :
			"Training rows "+sets[0].matrix.inputs.length+" != "+trainRows;
		assert(sets[1].matrix.inputs.length==validationRows) :
			"Validation rows "+sets[1].matrix.inputs.length+" != "+validationRows;
		checkShuffle(new SampleSet[]{sets[0]}, trainRows);
		checkShuffle(new SampleSet[]{sets[1]}, validationRows);
	}
	
	/** Checks output cells and Matrix.toString's private range summary with default midpoint. */
	private static void checkMatrix(Matrix matrix, float[][] expectedOutputs, String expectedMean){
		checkMatrix(matrix, expectedOutputs, expectedMean, "midpoint=0.5\n");
	}

	/** Checks output cells and Matrix.toString's private range summary. */
	private static void checkMatrix(Matrix matrix, float[][] expectedOutputs,
			String expectedMean, String expectedMidpoint){
		assert(matrix.inputs.length==expectedOutputs.length) :
			"Matrix rows "+matrix.inputs.length+" != "+expectedOutputs.length;
		for(int i=0; i<expectedOutputs.length; i++){
			final float[] expected=expectedOutputs[i];
			assert(matrix.outputs[i].length==expected.length) :
				"Output width "+matrix.outputs[i].length+" != "+expected.length+" at row "+i;
			for(int j=0; j<expected.length; j++){
				assert(matrix.outputs[i][j]==expected[j]) :
					"Output["+i+"]["+j+"] "+matrix.outputs[i][j]+" != "+expected[j];
			}
		}
		final String text=matrix.toString();
		assert(text.startsWith("[input, target0, target1]\n")) :
			"## columns did not reach Matrix.toString: "+text;
		assert(text.contains("outputs=2\n")) : text;
		assert(text.contains(expectedMidpoint)) : text;
		assert(text.contains(expectedMean)) : text;
		assert(!text.contains("Infinity") && !text.contains("NaN")) : text;
	}
	
	/** Returns the literal target paired with each single-output fixture input. */
	private static float expectedSingleTarget(int input){
		return input==0 ? 0.1f : input==1 ? 0.4f : input==2 ? 0.9f :
			Float.NaN;
	}
}
