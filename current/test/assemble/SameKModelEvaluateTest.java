package test.assemble;

import assemble.SameKModelEvaluate;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Runs an eight-row analytic model fixture through the real evaluator in one JVM.
 * Checks scalar predictions and curves, then repeats prediction checks with simdff=true;
 * this tests the per-model switch, without forcing a particular hardware execution path.
 * Score/label expectations are literal and independent of trained weights or biological truth.
 * Explicit checker failures remain active when Java assertions are disabled.
 * @author Fischl
 */
public final class SameKModelEvaluateTest{

	/**
	 * Requires one new fixture directory. Runs positive inference and malformed-input cases;
	 * malformed corpus cases require IllegalArgumentException and no completion marker. Also checks
	 * invalid score packing, signed-zero grouping and empty-curve totals. Leaves artifacts
	 * for inspection; an aborted run can leave partial outputs.
	 * Optional ratecheck=t runs only the three ASM027 rate-checker contracts.
	 */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKModelEvaluateTest.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs fixture modes after removing JVM/global launcher arguments. */
	private static void run(final String[] args){
		final boolean ratecheck=args.length==2 && args[1].equals("ratecheck=t");
		if((args.length!=1 && !ratecheck) || !new File(args[0]).mkdir()){
			throw new IllegalArgumentException("Supply a new fixture directory, optionally followed by ratecheck=t.");
		}
		final String root=args[0];
		if(ratecheck){checkRateRejection(root); return;}
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

	/** Captures normal/NA controls and the unchanged checker's rejection of a single NaN recall field. */
	private static void checkRateRejection(final String root){
		System.err.println("RATE_CHECK_ASSERTIONS enabled="+SameKModelEvaluateTest.class.desiredAssertionStatus());
		fixture(root, "valid");
		invoke(root, root+"/good", false);
		checkCurve(root+"/good/candidate.tsv", 0, 6);
		checkCurve(root+"/good/control.tsv", 1, 3);
		writeNanRateCurve(root+"/good/candidate.tsv", root+"/nan.tsv");
		boolean rejected=false;
		try{checkCurve(root+"/nan.tsv", 0, 6);}
		catch(AssertionError failure){
			if(!"Rate denominator incorrect.".equals(failure.getMessage())){throw failure;}
			rejected=true;
		}
		recordRateChecks(root, rejected);
	}

	/**
	 * Copies the six-row analytic candidate/grey curve, changing only the first row's recall
	 * token to NaN. Its reject-all row has four positives, two negatives and defined recall0.
	 * Callers must first validate the original curve with their unmodified checker.
	 */
	static void writeNanRateCurve(final String input, final String output){
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 fields=new LineParser1('\t');
		final ByteBuilder changed=new ByteBuilder();
		int rows=0;
		try{
			final byte[] header=file.nextLine();
			if(header==null){throw new IllegalStateException("Missing source curve header.");}
			changed.append(header).nl();
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				fields.set(line);
				if(fields.terms()!=12 || rows>=6){throw new IllegalStateException("Unexpected analytic curve shape.");}
				if(rows==0 && (fields.parseInt(0)!=Float.floatToRawIntBits(Math.nextUp(1.0f)) ||
						fields.parseLong(4)!=0 || fields.parseLong(5)!=0 || fields.parseLong(6)!=2 ||
						fields.parseLong(7)!=4 || !fields.termEquals("NA", 8) || fields.parseDouble(9)!=0)){
					throw new IllegalStateException("Expected reject-all row with defined zero recall and undefined precision.");
				}
				if(rows==0){
					for(int column=0; column<12; column++){
						if(column>0){changed.tab();}
						if(column==9){changed.append("NaN");}else{fields.appendTerm(changed, column);}
					}
				}else{changed.append(line);}
				changed.nl(); rows++;
			}
		}finally{if(file.close()){throw new IllegalStateException("Source curve read failed.");}}
		if(rows!=6){throw new IllegalStateException("Expected six analytic curve rows.");}
		write(output, changed);
		verifyNanRateCurve(input, output);
	}

	/** Reads both saved files and requires exactly the selected recall-token replacement. */
	private static void verifyNanRateCurve(final String input, final String output){
		final ByteFile original=ByteFile.makeByteFile(input, true);
		try{
			final ByteFile mutated=ByteFile.makeByteFile(output, true);
			try{
				final LineParser1 before=new LineParser1('\t'), after=new LineParser1('\t');
				final ByteBuilder restored=new ByteBuilder();
				int row=0;
				for(byte[] expected=original.nextLine(); expected!=null; expected=original.nextLine()){
					byte[] actual=mutated.nextLine();
					if(actual==null){throw new IllegalStateException("Mutated curve lost a row.");}
					if(row==1){
						before.set(expected); after.set(actual);
						if(before.terms()!=12 || after.terms()!=12 || before.parseDouble(9)!=0 || !after.termEquals("NaN", 9)){
							throw new IllegalStateException("Expected only the defined recall0 token to become NaN.");
						}
						for(int column=0; column<12; column++){
							if(column>0){restored.tab();}
							if(column==9){before.appendTerm(restored, column);}else{after.appendTerm(restored, column);}
						}
						actual=restored.toBytes();
					}
					if(!Arrays.equals(expected, actual)){throw new IllegalStateException("Unexpected mutation outside the selected rate field.");}
					row++;
				}
				if(row!=7 || mutated.nextLine()!=null){throw new IllegalStateException("Mutated curve row count changed.");}
			}finally{if(mutated.close()){throw new IllegalStateException("Mutated curve verification read failed.");}}
		}finally{if(original.close()){throw new IllegalStateException("Original curve verification read failed.");}}
	}

	/** Saves all three contracts after both unchanged controls passed; baseline failure remains visible. */
	static void recordRateChecks(final String root, final boolean rejected){
		final ByteBuilder checks=new ByteBuilder("case\tpassed\tdetail\n");
		checks.append("normal-curve\t1\tliteral group/confusion counts and rate tolerances\n");
		checks.append("na-denominator\t1\tall-positive controls retain undefined negative-class rates as NA\n");
		checks.append("nan-defined-rate\t").append(rejected ? 1 : 0)
				.append("\trecall denominator=4; required rate-diagnostic AssertionError\n");
		write(root+"/rate_checks.tsv", checks);
		System.err.println("RATE_CHECK_CONTRACTS checks=3 failed="+(rejected ? 0 : 1));
		if(!rejected){throw new AssertionError("NaN defined rate was accepted; inspect "+root+"/rate_checks.tsv");}
	}

	/** Calls the public CLI parser, keeping external heldout rows out of the selected feature stream. */
	private static void invoke(final String root, final String out, final boolean simd){
		SameKModelEvaluate.main(new String[]{"in="+root+"/validation.tsv", "rows="+root+"/rows.tsv",
				"net="+root+"/model.bbnet", "fold=validation", "out="+out, "simdff="+simd, "t=1"});
	}

	/**
	 * Writes a27-to1 linear model returning the first feature, six candidates and two controls.
	 * All eight selected rows belong to validation; interleaved train/test provenance rows
	 * have no entries in its feature table. A named mutation alters the first selected row.
	 */
	private static void fixture(final String root, final String variant){
		final ByteBuilder model=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 27 1\n#edges 27\n#coding decimal\n##layer 1\nC28 LINEAR 0 1");
		for(int i=1; i<27; i++){model.append(" 0");}
		write(root+"/model.bbnet", model.nl());
		final ByteBuilder data=new ByteBuilder("#dims\t27\t1\n#schema=samek27-v1; unknown labels excluded; no class balancing\n");
		final ByteBuilder meta=new ByteBuilder("fold\trow\tgenome\tkind\tevent\ttarget\ntrain\t1\tg0001\tcandidate\t1\t1\n");
		for(int i=0; i<SCORES.length; i++){
			if(i==0 && variant.equals("nonfinite")){data.append("NaN");}else{data.appendSlow(i==0 && variant.equals("badscore") ? 2 : SCORES[i]);}
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

	/**
	 * Checks eight ordered predictions, literal float bits and row provenance without calling
	 * the evaluator's packing helper. Skips the header; requires no extra prediction rows.
	 */
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

	/**
	 * Brute-forces each emitted cutoff against the unsorted literal population. Checks descending
	 * cutoffs, tied group/confusion counts, rate tolerances and NA for zero denominators.
	 * Defined-rate comparisons reject unordered NaN errors as well as excessive errors (ASM027).
	 * Skips the header and checks the expected row count; does not compare a literal cutoff list.
	 */
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
					if(SCORES[i]>=cutoff){if(LABELS[i]==1){tp++;}else{fp++;}}else{if(LABELS[i]==1){fn++;}else{tn++;}}
					if(SCORES[i]==cutoff){if(LABELS[i]==1){gp++;}else{gn++;}}
				}
				if(fields.terms()!=12 || fields.parseFloat(1)!=cutoff || !(cutoff<previous) ||
						fields.parseLong(2)!=gp || fields.parseLong(3)!=gn || fields.parseLong(4)!=tp ||
						fields.parseLong(5)!=fp || fields.parseLong(6)!=tn || fields.parseLong(7)!=fn){
					throw new AssertionError("Curve threshold/group/confusion counts differ from brute force.");
				}
				final long[] num={tp, tp, fp, fn}, den={tp+fp, tp+fn, fp+tn, tp+fn};
				//ASM027: Negate the ordered tolerance test so NaN cannot pass as a small error.
				//Finite tolerance decisions and rejection of infinite errors are unchanged.
				for(int i=0; i<4; i++){
					if(den[i]==0){if(!fields.termEquals("NA", 8+i)){throw new AssertionError("Undefined rate was not NA.");}}else if(!(Math.abs(fields.parseDouble(8+i)-num[i]/(double)den[i])<=1e-9)){throw new AssertionError("Rate denominator incorrect.");}
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

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** First six scores are candidates (including a mixed-label tie); last two are controls. */
	private static final float[] SCORES={1, 0.75f, 0.75f, 0.5f, 0.25f, 0, 0.75f, 0.25f};
	/** Literal binary labels in SCORES order; both controls are positive. */
	private static final int[] LABELS={1, 1, 0, 1, 0, 1, 1, 1};
}
