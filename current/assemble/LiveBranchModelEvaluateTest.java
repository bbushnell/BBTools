package assemble;

import java.io.File;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Analytic tests of the12-input adapter, independent of fitted models. @author Fischl */
public final class LiveBranchModelEvaluateTest {

	/** Checks all attainable thresholds against literal scores and refuses malformed pairs. */
	public static void main(final String[] args){
		if(args.length!=1 || !new File(args[0]).mkdir()){throw new IllegalArgumentException("Supply a new fixture directory.");}
		final String root=args[0];
		fixture(root, "valid");
		for(boolean simd : new boolean[]{false, true}){
			for(boolean controls : new boolean[]{false, true}){
				final String out=root+"/result-"+simd+"-"+controls;
				invoke(root, out, simd, controls);
				checkPredictions(out+"/scores.tsv", controls);
				checkCurve(out+"/"+(controls ? "control" : "grey")+".tsv", controls);
			}
		}
		for(String variant : VARIANTS){
			final String dir=root+"/"+variant;
			if(!new File(dir).mkdir()){throw new AssertionError("Cannot create negative fixture.");}
			fixture(dir, variant);
			try{invoke(dir, dir+"/out", false, false); throw new AssertionError("Accepted malformed "+variant);}
			catch(IllegalArgumentException expected){System.err.println("EXPECTED_REJECTION "+variant+": "+expected.getMessage());}
			if(new File(dir+"/out/complete.txt").exists()){throw new AssertionError("Bad input published completion: "+variant);}
		}
		try{invoke(root, root, false, false); throw new AssertionError("Accepted existing output.");}
		catch(IllegalArgumentException expected){/* Existing files must remain untouched. */}
		System.err.println("LIVE_BRANCH_EVALUATE_TEST_PASS analytic_rows=8 malformed="+(VARIANTS.length+1)+" scalar_simd_fixture_parity=true");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fixture Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Calls the actual CLI, explicitly selecting both SIMD switches. */
	private static void invoke(final String root, final String out, final boolean simd, final boolean controls){
		final String partition=controls ? "validation_controls" : "validation";
		LiveBranchModelEvaluate.main(new String[]{"in="+root+"/"+partition+".tsv", "rows="+root+"/rows.tsv",
				"net="+root+"/model.bbnet", "partition="+partition, "out="+out,
				"simdff="+simd, "simdfma="+simd, "t=1"});
	}

	/** A12-to1 linear network returns the first feature exactly, including probability endpoints. */
	private static void fixture(final String root, final String variant){
		assert(new File(root).isDirectory()) : "Fixture writers require an existing isolated scratch directory.";
		final int width=variant.equals("modelwidth") ? 11 : 12;
		final ByteBuilder model=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims ");
		model.append(width).append(" 1\n#edges ").append(width).append("\n#coding decimal\n##layer 1\nC")
				.append(width+1).append(" LINEAR 0 1");
		for(int i=1; i<width; i++){model.append(" 0");}
		write(root+"/model.bbnet", model.nl());
		final String header="#dims\t12\t1\n"+LiveBranchTrainingCorpus.SCHEMA+"\n";
		final ByteBuilder data=new ByteBuilder(variant.equals("header") ? "#dims\t27\t1\nwrong\n" : header);
		final ByteBuilder controls=new ByteBuilder(header), meta=new ByteBuilder(LiveBranchModelEvaluate.ROW_HEADER+"\n");
		appendMeta(meta, "train", 1, "g0001", 24, 8, "grey", 1);
		for(int i=0; i<SCORES.length; i++){
			final ByteBuilder row=new ByteBuilder();
			if(i==0 && variant.equals("nonfinite")){row.append("NaN");}
			else{row.appendSlow(i==0 && variant.equals("badscore") ? 2 : SCORES[i]);}
			for(int j=1; j<12; j++){row.tab().appendSlow(j==9 ? (float)(0.125*Math.log(24)/Math.log(2)) : 0);}
			row.tab().append(LABELS[i]);
			if(i==0 && variant.equals("width")){row.append("\t0");}
			if(!(i==5 && variant.equals("truncated"))){(i<6 ? data : controls).append(row).nl();}
			appendMeta(meta, i<6 ? "validation" : "validation_controls",
					i==0 && variant.equals("ordinal") ? 2 : i<6 ? i+1 : i-5,
					i==0 && variant.equals("genome") ? "bad" : "g0002",
					i==0 && variant.equals("k") ? 25 : 24,
					i==0 && variant.equals("word") ? -1 : i,
					i<6 && !(i==0 && variant.equals("role")) ? "grey" : "control",
					i==0 && variant.equals("target") ? 0 : LABELS[i]);
			if(i==3){appendMeta(meta, "train", variant.equals("otherordinal") ? 3 : 2, "g0001", 24, 9, "control", 0);}
		}
		if(variant.equals("extra")){data.append("0\n");}
		write(root+"/validation.tsv", data);
		write(root+"/validation_controls.tsv", controls);
		write(root+"/rows.tsv", meta);
	}

	/** Emits unambiguous canonical words and explicit partition-local ordinals. */
	private static void appendMeta(final ByteBuilder out, final String partition, final int row,
			final String genome, final int k, final int id, final String role, final int label){
		assert(row>0) : "The producer sidecar uses one-based row numbers within each partition.";
		out.append(partition).tab().append(row).tab().append(genome).tab().append(k).tab();
		for(int i=0; i<22; i++){out.append('A');}
		out.append(id<0 ? 'N' : "ACGT".charAt((id>>>2)&3)).append("ACGT".charAt(id&3));
		out.tab().append(role).tab().append(label).nl();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Check Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks source identity and exact score bits without using the packing implementation. */
	private static void checkPredictions(final String path, final boolean controls){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 p=new LineParser1('\t');
		try{
			file.nextLine();
			for(int i=controls ? 6 : 0; i<(controls ? 8 : 6); i++){
				final byte[] line=file.nextLine();
				if(line==null){throw new AssertionError("Missing analytic prediction.");}
				p.set(line);
				if(p.terms()!=9 || p.parseInt(1)!=(controls ? i-5 : i+1) || !p.termEquals("g0002", 2) ||
						!p.termEquals(controls ? "validation_controls" : "validation", 0) ||
						!p.termEquals(controls ? "control" : "grey", 5) || p.parseInt(6)!=LABELS[i] ||
						p.parseInt(7)!=Float.floatToRawIntBits(SCORES[i]) || p.parseFloat(8)!=SCORES[i]){
					throw new AssertionError("Analytic prediction or provenance mismatch.");
				}
			}
			if(file.nextLine()!=null){throw new AssertionError("Unselected partition leaked into scores.");}
		}finally{if(file.close()){throw new AssertionError("Prediction read failed.");}}
	}

	/** Brute-forces every tied group and class-conditional denominator against unsorted truth. */
	private static void checkCurve(final String path, final boolean controls){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 p=new LineParser1('\t');
		float previous=Float.POSITIVE_INFINITY;
		int rows=0;
		try{
			file.nextLine();
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				p.set(line);
				final float cutoff=Float.intBitsToFloat(p.parseInt(0));
				long tp=0, fp=0, tn=0, fn=0, gp=0, gn=0;
				for(int i=controls ? 6 : 0; i<(controls ? 8 : 6); i++){
					if(SCORES[i]>=cutoff){if(LABELS[i]==1){tp++;}else{fp++;}}
					else{if(LABELS[i]==1){fn++;}else{tn++;}}
					if(SCORES[i]==cutoff){if(LABELS[i]==1){gp++;}else{gn++;}}
				}
				if(p.terms()!=12 || p.parseFloat(1)!=cutoff || !(cutoff<previous) || p.parseLong(2)!=gp ||
						p.parseLong(3)!=gn || p.parseLong(4)!=tp || p.parseLong(5)!=fp || p.parseLong(6)!=tn || p.parseLong(7)!=fn){
					throw new AssertionError("Curve differs from brute-force threshold decisions.");
				}
				final long[] num={tp, tp, fp, fn}, den={tp+fp, tp+fn, fp+tn, tp+fn};
				for(int j=0; j<4; j++){
					if(den[j]==0){if(!p.termEquals("NA", j+8)){throw new AssertionError("Undefined rate must be NA.");}}
					else if(Math.abs(p.parseDouble(j+8)-num[j]/(double)den[j])>1e-9){throw new AssertionError("Wrong rate denominator.");}
				}
				previous=cutoff; rows++;
			}
			if(rows!=(controls ? 3 : 6)){throw new AssertionError("Missing reject-all row, lost threshold, or split tie.");}
		}finally{if(file.close()){throw new AssertionError("Curve read failed.");}}
	}

	/** Writes isolated fixture files through the native writer, checking termination. */
	private static void write(final String path, final ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new IllegalStateException("Fixture write failed.");}
	}

	private static final float[] SCORES={1, 0.75f, 0.75f, 0.5f, 0.25f, 0, 0.75f, 0.25f};
	private static final int[] LABELS={1, 1, 0, 1, 0, 1, 1, 1};
	private static final String[] VARIANTS={"target", "ordinal", "width", "nonfinite", "badscore", "header",
			"modelwidth", "genome", "k", "word", "role", "otherordinal", "extra", "truncated"};
}
