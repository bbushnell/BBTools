package test.assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;

import assemble.LiveBranchCurveSummary;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Exercises curve summarization through its real CLI using literal six-observation curves.
 * Covers normal and exponent-form thresholds, malformed curves and negative-only TP ties.
 * Uses explicit failures for checks, independent of whether JVM assertions are enabled.
 * @author Fischl
 */
public final class LiveBranchCurveSummaryTest{

	/**
	 * Requires a new fixture directory and retains inputs/outputs for inspection.
	 * Removes recognized JVM flags before validating the positional directory argument.
	 * Seven malformed variants require IllegalArgumentException and no published summary.
	 * The negative-only control checks reject-all selection in each observed output row;
	 * that control does not separately validate its header or number of rows.
	 */
	public static void main(String[] args){
		args=SameKJoinTest.removeJavaFlags(args);
		if(args.length!=1 || !new File(args[0]).mkdir()){throw new IllegalArgumentException("Supply a new fixture directory.");}
		final String root=args[0];
		write(root+"/curve.tsv", CURVE);
		invoke(root+"/curve.tsv", root+"/summary.tsv");
		check(root+"/summary.tsv", 0);
		final String[] exponent=CURVE.clone();
		exponent[6]=Float.floatToRawIntBits(0.0001f)+"\t1.0E-4\t1\t0\t4\t2\t0\t0\t0.666666666667\t1\t1\t0";
		write(root+"/exponent.tsv", exponent);
		invoke(root+"/exponent.tsv", root+"/exponent.out");
		check(root+"/exponent.out", 0.0001f);
		for(String variant : new String[]{"header", "width", "population", "order", "group", "rate", "truncated"}){
			final String[] changed=CURVE.clone();
			if(variant.equals("header")){changed[0]="wrong";}
			else if(variant.equals("width")){changed[1]+="\t0";}
			else if(variant.equals("population")){changed[1]=changed[1].replace("\t2\t4\t", "\t2\t5\t");}
			else if(variant.equals("order")){changed[3]=changed[2];}
			else if(variant.equals("group")){changed[3]=changed[3].replace("\t1\t1\t2\t1\t", "\t0\t1\t2\t1\t");}
			else if(variant.equals("rate")){changed[2]=changed[2].replace("\t0.25\t", "\t0.5\t");}
			else{changed[changed.length-1]=null;}
			write(root+"/"+variant+".tsv", changed);
			try{invoke(root+"/"+variant+".tsv", root+"/"+variant+".out"); throw new AssertionError("Accepted malformed "+variant);}
			catch(IllegalArgumentException expected){System.err.println("EXPECTED_REJECTION "+variant+": "+expected.getMessage());}
			if(new File(root+"/"+variant+".out").exists()){throw new AssertionError("Published malformed curve summary.");}
		}
		write(root+"/negative.tsv", new String[]{CURVE[0],
				"1065353217\t1.0000001\t0\t0\t0\t0\t6\t0\tNA\tNA\t0\tNA",
				"1065353216\t1.0\t0\t6\t0\t6\t0\t0\t0\tNA\t1\tNA"});
		invoke(root+"/negative.tsv", root+"/negative.out");
		final ByteFile negative=ByteFile.makeByteFile(root+"/negative.out", true);
		final LineParser1 p=new LineParser1('\t');
		try{
			negative.nextLine();
			for(byte[] row=negative.nextLine(); row!=null; row=negative.nextLine()){
				p.set(row);
				if(p.parseInt(1)!=1065353217 || p.parseLong(5)!=0 || p.parseLong(6)!=0){
					throw new AssertionError("A TP tie must retain reject-all, not spend false-rescue budget.");
				}
			}
		}finally{if(negative.close()){throw new AssertionError("Negative-only fixture read failed.");}}
		System.err.println("LIVE_BRANCH_CURVE_SUMMARY_TEST_PASS budgets=15 malformed=7 conservative_ties=true");
	}

	/** Exercises the public CLI with the literal six-observation population. */
	private static void invoke(final String input, final String output){
		LiveBranchCurveSummary.main(new String[]{"in="+input, "out="+output, "expected=6", "t=1"});
	}

	/**
	 * Checks literal optima: TP 1/FP 0 at 1, TP 3/FP 1 at 0.5, TP 4/FP 2 at minimum.
	 * Minimum is 0 for the ordinary fixture or 0.0001 for exponent-form input. Requires the
	 * exact header and 15 rows, but does not compare the complete budget list or copied rate fields.
	 */
	private static void check(final String path, final float minimum){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try{
			if(!new String(file.nextLine(), StandardCharsets.UTF_8).equals("fp_budget\t"+CURVE[0])){
				throw new AssertionError("Summary changed the source curve columns.");
			}
			for(byte[] row=file.nextLine(); row!=null; row=file.nextLine()){
				p.set(row);
				final long budget=p.parseLong(0), tp=budget==0 ? 1 : budget==1 ? 3 : 4, fp=budget==0 ? 0 : budget==1 ? 1 : 2;
				final float cutoff=budget==0 ? 1 : budget==1 ? 0.5f : minimum;
				if(p.terms()!=13 || p.parseLong(5)!=tp || p.parseLong(6)!=fp || Float.parseFloat(p.parseString(2))!=cutoff){
					throw new AssertionError("Summary did not choose the highest-TP attainable row under its FP budget.");
				}
				rows++;
			}
			if(rows!=15){throw new AssertionError("Lost a declared false-rescue budget.");}
		}finally{if(file.close()){throw new AssertionError("Summary fixture read failed.");}}
	}

	/** Writes only tiny literal fixtures; production parsing remains byte-oriented. */
	private static void write(final String path, final String[] lines){
		assert(lines!=null && lines.length>0) : "A curve fixture includes its schema header and observations.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		final ByteBuilder out=new ByteBuilder();
		for(String line : lines){if(line!=null){out.append(line).nl();}}
		writer.print(out);
		if(writer.poisonAndWait()){throw new IllegalStateException("Fixture write failed.");}
	}

	/** Four positives and two negatives, including a mixed-label tie at 0.75. */
	private static final String[] CURVE={
			"threshold_bits\tthreshold\tgroup_positive\tgroup_negative\ttp\tfp\ttn\tfn\tprecision\trecall\tfpr\tfnr",
			"1065353217\t1.0000001\t0\t0\t0\t0\t2\t4\tNA\t0\t0\t1",
			"1065353216\t1.0\t1\t0\t1\t0\t2\t3\t1\t0.25\t0\t0.75",
			"1061158912\t0.75\t1\t1\t2\t1\t1\t2\t0.666666666667\t0.5\t0.5\t0.5",
			"1056964608\t0.5\t1\t0\t3\t1\t1\t1\t0.75\t0.75\t0.5\t0.25",
			"1048576000\t0.25\t0\t1\t3\t2\t0\t1\t0.6\t0.75\t1\t0.25",
			"0\t0.0\t1\t0\t4\t2\t0\t0\t0.666666666667\t1\t1\t0"};
}
