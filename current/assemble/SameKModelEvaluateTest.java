package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.LongList;

/** Analytic score/curve tests independent of trained weights or biological labels. @author Fischl */
public final class SameKModelEvaluateTest {

	/** Exercises exact inference, tied thresholds, population routing and loud malformed-input failures. */
	public static void main(final String[] args){
		if(args.length!=1 || !new File(args[0]).mkdir()){throw new IllegalArgumentException("Supply a new fixture directory.");}
		final String root=args[0];
		fixture(root, "valid");
		invoke(root, root+"/good", false);
		checkPredictions(root+"/good/scores.tsv");
		checkCurve(root+"/good/candidate.tsv", 0, 6);
		checkCurve(root+"/good/control.tsv", 1, 3);
		checkCurve(root+"/good/candidate_points.tsv", 0, 20);
		checkCurve(root+"/good/control_points.tsv", 1, 20);
		invoke(root, root+"/simd", true);
		checkPredictions(root+"/simd/scores.tsv");
		for(String variant : new String[]{"badtarget", "badordinal", "badwidth", "nonfinite", "badscore"}){
			final String dir=root+"/"+variant;
			if(!new File(dir).mkdir()){throw new AssertionError("Cannot create negative fixture.");}
			fixture(dir, variant);
			try{invoke(dir, dir+"/out", false); throw new AssertionError("Accepted malformed "+variant);}
			catch(IllegalArgumentException expected){System.err.println("EXPECTED_REJECTION "+variant+": "+expected.getMessage());}
			if(new File(dir+"/out/complete.txt").exists()){throw new AssertionError("Bad input published a completion marker.");}
		}
		try{invoke(root, root, false); throw new AssertionError("Accepted existing input directory as output.");}
		catch(IllegalArgumentException expected){System.err.println("EXPECTED_REJECTION alias: "+expected.getMessage());}
		for(float score : new float[]{Float.NaN, Float.POSITIVE_INFINITY, -1, 2}){
			try{SameKModelEvaluate.pack(score, 1); throw new AssertionError("Accepted invalid score.");}
			catch(IllegalArgumentException expected){/* Expected probability-domain failure. */}
		}
		if(SameKModelEvaluate.pack(-0.0f, 1)!=SameKModelEvaluate.pack(0.0f, 1)){throw new AssertionError("Signed zero splits one attainable threshold.");}
		final ByteStreamWriter empty=new ByteStreamWriter(root+"/empty.tsv", false, false, false);
		empty.start();
		final long[] counts=SameKModelEvaluate.writeCurve(new LongList(), empty);
		if(empty.poisonAndWait() || !Arrays.equals(counts, new long[3])){throw new AssertionError("Empty curve accounting failed.");}
		System.err.println("SAME_K_EVALUATE_TEST_PASS analytic_rows=8 malformed=6 simd_parity=true");
	}

	/** Calls the public CLI parser, keeping external heldout rows out of the selected feature stream. */
	private static void invoke(final String root, final String out, final boolean simd){
		SameKModelEvaluate.main(new String[]{"in="+root+"/validation.tsv", "rows="+root+"/rows.tsv",
				"net="+root+"/model.bbnet", "fold=validation", "out="+out, "simdff="+simd, "t=1"});
	}

	/** A27-to1 linear model returns the first feature exactly, including both probability endpoints. */
	private static void fixture(final String root, final String variant){
		final ByteBuilder model=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 27 1\n#edges 27\n#coding decimal\n##layer 1\nC28 LINEAR 0 1");
		for(int i=1; i<27; i++){model.append(" 0");}
		write(root+"/model.bbnet", model.nl());
		final ByteBuilder data=new ByteBuilder("#dims\t27\t1\n#schema=samek27-v1; unknown labels excluded; no class balancing\n");
		final ByteBuilder meta=new ByteBuilder("fold\trow\tgenome\tkind\tevent\ttarget\ntrain\t1\tg0001\tcandidate\t1\t1\n");
		for(int i=0; i<SCORES.length; i++){
			if(i==0 && variant.equals("nonfinite")){data.append("NaN");}
			else{data.appendSlow(i==0 && variant.equals("badscore") ? 2 : SCORES[i]);}
			for(int j=1; j<27; j++){data.append("\t0");}
			data.tab().append(LABELS[i]);
			if(i==0 && variant.equals("badwidth")){data.append("\t0");}
			data.nl();
			meta.append("validation\t").append(i==0 && variant.equals("badordinal") ? 2 : i+1)
					.append("\tg0002\t").append(i<6 ? "candidate" : "control").tab().append(i<6 ? i+1 : i-5)
					.tab().append(i==0 && variant.equals("badtarget") ? 0 : LABELS[i]).nl();
			if(i==3){meta.append("test\t1\tg0003\tcandidate\t1\t0\n");}
		}
		write(root+"/validation.tsv", data); write(root+"/rows.tsv", meta);
	}

	/** Checks exact float bits and row provenance, without calling the evaluator's packing helper. */
	private static void checkPredictions(final String path){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 fields=new LineParser1('\t');
		try{
			file.nextLine();
			for(int i=0; i<SCORES.length; i++){
				final byte[] line=file.nextLine();
				if(line==null){throw new AssertionError("Missing prediction.");}
				fields.set(line);
				if(fields.terms()!=8 || !fields.termEquals("validation", 0) || fields.parseInt(1)!=i+1 ||
						!fields.termEquals("g0002", 2) || !fields.termEquals(i<6 ? "candidate" : "control", 3) ||
						fields.parseInt(4)!=(i<6 ? i+1 : i-5) || fields.parseInt(5)!=LABELS[i] ||
						fields.parseInt(6)!=Float.floatToRawIntBits(SCORES[i]) || fields.parseFloat(7)!=SCORES[i]){
					throw new AssertionError("Prediction or source identity differs from analytic fixture.");
				}
			}
			if(file.nextLine()!=null){throw new AssertionError("Unselected fold leaked into prediction output.");}
		}finally{if(file.close()){throw new AssertionError("Cannot read prediction fixture.");}}
	}

	/** Brute-forces every threshold against unsorted scalar truth, including tied positive/negative scores. */
	private static void checkCurve(final String path, final int kind, final int expectedRows){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 fields=new LineParser1('\t');
		float previous=Float.POSITIVE_INFINITY;
		int rows=0;
		try{
			file.nextLine();
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				fields.set(line);
				final float cutoff=Float.intBitsToFloat(fields.parseInt(0));
				long tp=0, fp=0, tn=0, fn=0, gp=0, gn=0;
				for(int i=kind==0 ? 0 : 6; i<(kind==0 ? 6 : 8); i++){
					if(SCORES[i]>=cutoff){if(LABELS[i]==1){tp++;}else{fp++;}}
					else{if(LABELS[i]==1){fn++;}else{tn++;}}
					if(SCORES[i]==cutoff){if(LABELS[i]==1){gp++;}else{gn++;}}
				}
				if(fields.terms()!=12 || fields.parseFloat(1)!=cutoff || !(cutoff<previous) ||
						fields.parseLong(2)!=gp || fields.parseLong(3)!=gn || fields.parseLong(4)!=tp ||
						fields.parseLong(5)!=fp || fields.parseLong(6)!=tn || fields.parseLong(7)!=fn){
					throw new AssertionError("Curve threshold/group/confusion counts differ from brute force.");
				}
				final long[] num={tp, tp, fp, fn}, den={tp+fp, tp+fn, fp+tn, tp+fn};
				for(int i=0; i<4; i++){
					if(den[i]==0){if(!fields.termEquals("NA", 8+i)){throw new AssertionError("Undefined rate was not NA.");}}
					else if(Math.abs(fields.parseDouble(8+i)-num[i]/(double)den[i])>1e-9){throw new AssertionError("Rate denominator incorrect.");}
				}
				previous=cutoff; rows++;
			}
			if(rows!=expectedRows){throw new AssertionError("Lost a declared threshold or split a tied group.");}
		}finally{if(file.close()){throw new AssertionError("Cannot read curve fixture.");}}
	}

	/** Writes tiny fixture text through native BBTools output, observing close errors. */
	private static void write(final String path, final ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new IllegalStateException("Fixture write failed.");}
	}

	private static final float[] SCORES={1, 0.75f, 0.75f, 0.5f, 0.25f, 0, 0.75f, 0.25f};
	private static final int[] LABELS={1, 1, 0, 1, 0, 1, 1, 1};
}
