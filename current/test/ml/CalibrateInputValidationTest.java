package test.ml;

import ml.Calibrate;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import structures.ByteBuilder;

/**
 * Exercises Calibrate argument and raw-score validation through the real CLI.
 * Rejecting cases retain inputs, stderr, and the observed rejection; finite
 * controls retain ordinary #cal, reliability, and LUT outputs.
 * @author Jean
 */
public class CalibrateInputValidationTest{

	/** Runs one validation case in a new artifact directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: CalibrateInputValidationTest <case> <artifact-dir>";
		final String name=args[0];
		final Path dir=Paths.get(args[1]);
		Files.createDirectory(dir);
		if(name.equals("restarts-zero")){
			rejectPairs(name, dir, writePairs(dir), "restarts", "restarts=0");
		}else if(name.equals("restarts-negative")){
			rejectPairs(name, dir, writePairs(dir), "restarts", "restarts=-1");
		}else if(name.equals("scorecol-negative")){
			rejectPairs(name, dir, writePairs(dir), "scorecol", "restarts=1", "scorecol=-1");
		}else if(name.equals("labelcol-negative")){
			rejectPairs(name, dir, writePairs(dir), "labelcol", "restarts=1", "labelcol=-1");
		}else if(name.equals("pair-nan")){
			rejectPairs(name, dir, writePairs(dir, "NaN"), "finite", "restarts=1");
		}else if(name.equals("pair-posinf")){
			rejectPairs(name, dir, writePairs(dir, "Infinity"), "finite", "restarts=1");
		}else if(name.equals("pair-neginf")){
			rejectPairs(name, dir, writePairs(dir, "-Infinity"), "finite", "restarts=1");
		}else if(name.equals("net-overflow")){
			rejectNet(name, dir, "finite");
		}else if(name.equals("finite-pair")){
			pass(name, dir, writePairs(dir), null, "restarts=1");
		}else if(name.equals("finite-net")){
			pass(name, dir, writeVectors(dir), writeNet(dir, false), "restarts=1");
		}else if(name.equals("net-unused-columns")){
			pass(name, dir, writeVectors(dir), writeNet(dir, false), "restarts=1",
				"scorecol=-1", "labelcol=-1");
		}else if(name.equals("column-override")){
			pass(name, dir, writeColumnPairs(dir), null, "restarts=1",
				"scorecol=2", "labelcol=0");
		}else{
			throw new IllegalArgumentException("Unknown Calibrate validation case: "+name);
		}
		System.err.println("CalibrateInputValidationTest PASS "+name);
	}

	/** Runs one direct-pair rejection case and asserts a clear validation failure. */
	private static void rejectPairs(final String label, final Path dir, final Path in,
			final String expected, final String... args) throws Exception{
		reject(label, dir, in, null, expected, args);
	}

	/** Runs the finite-input overflow net case and asserts a clear validation failure. */
	private static void rejectNet(final String label, final Path dir, final String expected)
			throws Exception{
		reject(label, dir, writeOverflowVectors(dir), writeNet(dir, true),
			expected, "restarts=1");
	}

	/** Runs one rejecting case and verifies that no optional output file was created. */
	private static void reject(final String label, final Path dir, final Path in, final Path net,
			final String expected, final String... args) throws Exception{
		final Path cal=dir.resolve("cal.txt");
		final Path rel=dir.resolve("reliability.tsv");
		final Path lut=dir.resolve("lut.tsv");
		final RunResult result=run(label, in, net, cal, rel, lut, args);
		writeResult(dir, result);
		if(!(result.rejection instanceof IllegalArgumentException)){
			throw new AssertionError("Expected IllegalArgumentException for "+label+
				"; observed "+result.diagnostic(), result.rejection);
		}
		final String message=String.valueOf(result.rejection.getMessage()).toLowerCase(Locale.ROOT);
		if(!message.contains(expected)){
			throw new AssertionError("Expected rejection for "+label+
				" to mention '"+expected+"': "+result.diagnostic());
		}
		check(!Files.exists(cal), "Rejected case created "+cal);
		check(!Files.exists(rel), "Rejected case created "+rel);
		check(!Files.exists(lut), "Rejected case created "+lut);
	}

	/** Runs one finite control and checks all ordinary outputs are finite. */
	private static void pass(final String label, final Path dir, final Path in, final Path net,
			final String... args) throws Exception{
		final Path cal=dir.resolve("cal.txt");
		final Path rel=dir.resolve("reliability.tsv");
		final Path lut=dir.resolve("lut.tsv");
		final RunResult result=run(label, in, net, cal, rel, lut, args);
		writeResult(dir, result);
		if(result.rejection!=null){
			throw new AssertionError("Ordinary case rejected: "+label+"; "+
				result.diagnostic(), result.rejection);
		}
		checkCalFile(cal);
		checkReliability(rel, label, 8);
		checkLut(lut, 8);
		checkReport(result.report, label, 8, cal, rel, lut);
	}

	/** Runs Calibrate with actual CLI flags and deterministic bounded fitting settings. */
	private static RunResult run(final String label, final Path in, final Path net,
			final Path cal, final Path rel, final Path lut, final String... args)
			throws Exception{
		final String[] base={
			"in="+in, "out="+cal, "rel="+rel, "lutout="+lut,
			"label="+label, "bins=4", "fitbins=16", "holdout=0", "maxiters=25"
		};
		final String[] all=Arrays.copyOf(base, base.length+args.length+(net==null ? 0 : 1));
		int idx=base.length;
		if(net!=null){all[idx++]="net="+net;}
		for(String arg : args){all[idx++]=arg;}

		final PrintStream oldErr=System.err;
		final Locale oldLocale=Locale.getDefault();
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(buffer, true, "UTF-8");
		RuntimeException rejection=null;
		try{
			System.setErr(ps);
			Locale.setDefault(Locale.ROOT);
			Calibrate.main(all);
		}catch(RuntimeException e){
			rejection=e;
		}finally{
			System.setErr(oldErr);
			Locale.setDefault(oldLocale);
			ps.close();
		}
		return new RunResult(new String(buffer.toByteArray(), StandardCharsets.UTF_8),
			rejection);
	}

	/** Retains the exact stderr report and observed rejection. */
	private static void writeResult(final Path dir, final RunResult result) throws Exception{
		Files.write(dir.resolve("report.stderr"), result.report.getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("rejection.txt"),
			(result.diagnostic()+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Writes eight finite score/label pairs in distinct fit bins. */
	private static Path writePairs(final Path dir) throws Exception{
		return writePairs(dir, "0.10");
	}

	/** Writes one caller-supplied raw score plus seven finite raw scores. */
	private static Path writePairs(final Path dir, final String firstScore) throws Exception{
		final Path in=dir.resolve("pairs.tsv");
		Files.write(in, Arrays.asList(
			firstScore+"\t0",
			"0.20\t0",
			"0.30\t0",
			"0.40\t0",
			"0.60\t1",
			"0.70\t1",
			"0.80\t1",
			"0.90\t1"
		), StandardCharsets.US_ASCII);
		return in;
	}

	/** Writes finite score/label pairs with score in column 2 and label in column 0. */
	private static Path writeColumnPairs(final Path dir) throws Exception{
		final Path in=dir.resolve("pairs.tsv");
		Files.write(in, Arrays.asList(
			"0\tignored\t0.10",
			"0\tignored\t0.20",
			"0\tignored\t0.30",
			"0\tignored\t0.40",
			"1\tignored\t0.60",
			"1\tignored\t0.70",
			"1\tignored\t0.80",
			"1\tignored\t0.90"
		), StandardCharsets.US_ASCII);
		return in;
	}

	/** Writes eight one-input vector rows scored by the finite retained network. */
	private static Path writeVectors(final Path dir) throws Exception{
		final Path in=dir.resolve("vectors.tsv");
		Files.write(in, Arrays.asList(
			"#dims\t1\t1",
			"0\t0",
			"1\t0",
			"2\t0",
			"3\t0",
			"4\t1",
			"5\t1",
			"6\t1",
			"7\t1"
		), StandardCharsets.US_ASCII);
		return in;
	}

	/** Writes finite vector rows, one of which overflows the output cell. */
	private static Path writeOverflowVectors(final Path dir) throws Exception{
		final Path in=dir.resolve("vectors.tsv");
		Files.write(in, Arrays.asList(
			"#dims\t1\t1",
			"2\t0",
			"0\t0",
			"0\t0",
			"0\t0",
			"0\t1",
			"0\t1",
			"0\t1",
			"0\t1"
		), StandardCharsets.US_ASCII);
		return in;
	}

	/** Writes a deterministic one-input linear net. */
	private static Path writeNet(final Path dir, final boolean overflow) throws Exception{
		final Path net=dir.resolve("net.bbnet");
		if(overflow){
			final ByteBuilder model=new ByteBuilder();
			model.append("##bbnet\n#version 1\n#concise\n#dense\n"+
				"#density 1\n#blocksize 1\n#seed 1\n#layers 2\n"+
				"#dims 1 1\n#edges 1\n#coding A48\n##layer 1\nC2 LINEAR ");
			model.appendFloatA48(0.25f).space().appendFloatA48(Float.MAX_VALUE).nl();
			Files.write(net, model.toBytes());
		}else{
			Files.write(net, Arrays.asList(
				"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
				"#seed 1", "#layers 2", "#dims 1 1", "#edges 1", "#coding decimal",
				"##layer 1", "C2 LINEAR 0.25 0.10"
			), StandardCharsets.US_ASCII);
		}
		return net;
	}

	/** Checks the standalone #cal file. */
	private static void checkCalFile(final Path path) throws Exception{
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(lines.size()==1, "Expected one #cal line in "+path+": "+lines);
		checkCalLine(lines.get(0));
	}

	/** Checks reliability table structure, row counts, and finite numeric columns. */
	private static void checkReliability(final Path path, final String label, final int total)
			throws Exception{
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(lines.size()>=5, "Expected dataset, #cal, headers, and rows in "+path+": "+lines);
		check(lines.get(0).equals("#dataset\t"+label),
			"Wrong reliability dataset line: "+lines.get(0));
		checkCalLine(lines.get(1));
		check(lines.get(2).equals("#reliability (binned by calibrated probability)"),
			"Missing reliability title: "+lines.get(2));
		check(lines.get(3).equals("#bin\tlo\thi\tn\tfrac\tmeanRaw\tmeanPred\tempirical\tgap"),
			"Missing reliability header: "+lines.get(3));
		int sum=0;
		for(int i=4; i<lines.size(); i++){
			final String[] split=lines.get(i).split("\t", -1);
			check(split.length==9, "Wrong reliability column count: "+lines.get(i));
			sum+=Integer.parseInt(split[3]);
			for(int j=1; j<split.length; j++){checkFinite(split[j], lines.get(i));}
		}
		check(sum==total, "Reliability rows cover "+sum+" observations, expected "+total);
	}

	/** Checks LUT structure, row count, and finite numeric columns. */
	private static void checkLut(final Path path, final int expectedRows) throws Exception{
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(lines.size()==expectedRows+2, "Unexpected LUT row count: "+lines);
		check(lines.get(0).equals("#lut\t"+expectedRows), "Wrong LUT count: "+lines.get(0));
		check(lines.get(1).equals("#rawScore\tprobability"), "Wrong LUT header: "+lines.get(1));
		for(int i=2; i<lines.size(); i++){
			final String[] split=lines.get(i).split("\t", -1);
			check(split.length==2, "Wrong LUT column count: "+lines.get(i));
			checkFinite(split[0], lines.get(i));
			checkFinite(split[1], lines.get(i));
		}
	}

	/** Checks the process report and every optional writer line. */
	private static void checkReport(final String report, final String label, final int rows,
			final Path cal, final Path rel, final Path lut){
		check(report.contains("Dataset:   "+label+"\n"),
			"Missing dataset label in report:\n"+report);
		check(report.contains("Rows="+rows+"  mean_label=0.50000  occupied_bins="),
			"Missing row summary in report:\n"+report);
		check(report.contains("log-loss:  raw="), "Missing log-loss line:\n"+report);
		check(report.contains("ECE(4):     raw="), "Missing ECE line:\n"+report);
		check(report.contains("params:    K="), "Missing parameter line:\n"+report);
		check(report.contains("#reliability (binned by calibrated probability)\n"),
			"Missing displayed reliability table:\n"+report);
		check(report.contains("Wrote "+cal+"\n"), "Missing cal write report:\n"+report);
		check(report.contains("Wrote "+rel+"\n"), "Missing reliability write report:\n"+report);
		check(report.contains("Wrote "+lut+" (8 knots)\n"), "Missing LUT write report:\n"+report);
		check(report.contains("\nTime: "), "Missing elapsed-time line:\n"+report);
		check(!report.contains("NaN"), "Unexpected NaN in report:\n"+report);
		check(!report.contains("Infinity"), "Unexpected infinity in report:\n"+report);
		for(String line : report.split("\n")){
			if(!line.startsWith("Time: ")){
				check(!line.toLowerCase(Locale.ROOT).contains("exception"),
					"Unexpected error-looking report line: "+line);
			}
		}
	}

	/** Checks the #cal prefix and four finite constants. */
	private static void checkCalLine(final String line){
		final String[] split=line.split(" ", -1);
		check(split.length==5 && split[0].equals("#cal"), "Malformed #cal line: "+line);
		for(int i=1; i<split.length; i++){checkFinite(split[i], line);}
	}

	/** Parses one finite double literal from a retained output. */
	private static void checkFinite(final String value, final String source){
		final double d=Double.parseDouble(value);
		check(Double.isFinite(d), "Expected finite value in "+source);
	}

	/** Throws a diagnostic assertion failure for a mismatched literal. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}

	/** Captured result of one direct Calibrate invocation. */
	private static final class RunResult{
		RunResult(final String report_, final RuntimeException rejection_){
			report=report_;
			rejection=rejection_;
		}

		/** Returns a stable one-line rejection diagnostic. */
		String diagnostic(){
			return rejection==null ? "NONE" :
				rejection.getClass().getName()+": "+rejection.getMessage();
		}

		final String report;
		final RuntimeException rejection;
	}

	private CalibrateInputValidationTest(){}
}
