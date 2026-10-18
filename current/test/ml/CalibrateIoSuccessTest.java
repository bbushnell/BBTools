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

/**
 * Preserves Calibrate's ordinary score/label and net-driven output paths while
 * reader and writer completion handling is tightened.
 * @author Jean
 */
public class CalibrateIoSuccessTest{

	/** Runs one ordinary success case. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: CalibrateIoSuccessTest <case> <artifact-dir>";
		final String name=args[0];
		final Path dir=Paths.get(args[1]);
		Files.createDirectory(dir);
		if(name.equals("pair-out")){
			pairOut(dir);
		}else if(name.equals("pair-rel")){
			pairReliability(dir);
		}else if(name.equals("pair-lutout")){
			pairLut(dir);
		}else if(name.equals("net")){
			net(dir);
		}else{
			throw new IllegalArgumentException("Unknown Calibrate I/O case: "+name);
		}
		System.err.println("CalibrateIoSuccessTest PASS "+name);
	}

	/** Runs direct score/label mode and writes the paste-ready #cal line. */
	private static void pairOut(final Path dir) throws Exception{
		final Path in=writePairs(dir);
		final Path out=dir.resolve("cal.txt");
		final String report=run("pair-out", "in="+in, "out="+out);
		writeReports(dir, report);
		checkCalFile(out);
		checkReport(report, "pair-out", "Wrote "+out, 8);
	}

	/** Runs direct score/label mode and writes the reliability table. */
	private static void pairReliability(final Path dir) throws Exception{
		final Path in=writePairs(dir);
		final Path out=dir.resolve("reliability.tsv");
		final String report=run("pair-rel", "in="+in, "rel="+out);
		writeReports(dir, report);
		checkReliability(out, "pair-rel", 8);
		checkReport(report, "pair-rel", "Wrote "+out, 8);
	}

	/** Runs direct score/label mode and writes the isotonic lookup table. */
	private static void pairLut(final Path dir) throws Exception{
		final Path in=writePairs(dir);
		final Path out=dir.resolve("lut.tsv");
		final String report=run("pair-lutout", "in="+in, "lutout="+out);
		writeReports(dir, report);
		checkLut(out, 8);
		checkReport(report, "pair-lutout", "Wrote "+out+" (8 knots)", 8);
	}

	/** Runs net-scored vector mode without optional output files. */
	private static void net(final Path dir) throws Exception{
		final Path in=writeVectors(dir);
		final Path net=writeNet(dir);
		final String report=run("net", "in="+in, "net="+net);
		writeReports(dir, report);
		checkReport(report, "net", null, 4);
	}

	/** Runs Calibrate with deterministic bounded optimizer settings. */
	private static String run(final String label, final String... args) throws Exception{
		final String[] all=Arrays.copyOf(args, args.length+6);
		all[args.length]="label="+label;
		all[args.length+1]="bins=4";
		all[args.length+2]="fitbins=16";
		all[args.length+3]="holdout=0";
		all[args.length+4]="restarts=1";
		all[args.length+5]="maxiters=25";
		final PrintStream oldErr=System.err;
		final Locale oldLocale=Locale.getDefault();
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(buffer, true, "UTF-8");
		try{
			System.setErr(ps);
			Locale.setDefault(Locale.ROOT);
			Calibrate.main(all);
		}finally{
			System.setErr(oldErr);
			Locale.setDefault(oldLocale);
			ps.close();
		}
		return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Writes the exact report and a path-normalized variant for baseline/fixed parity. */
	private static void writeReports(final Path dir, final String report) throws Exception{
		Files.write(dir.resolve("report.stderr"), report.getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("report-parity.stderr"),
			parityReport(report, dir).getBytes(StandardCharsets.US_ASCII));
	}

	/** Removes elapsed time and the exact artifact directory from a report. */
	private static String parityReport(final String report, final Path dir){
		return withoutTime(report).replace(dir.toString(), "<CASE_DIR>");
	}

	/** Removes only the elapsed-time line from an otherwise exact report. */
	private static String withoutTime(final String report){
		final StringBuilder sb=new StringBuilder();
		for(String line : report.split("\n", -1)){
			if(line.startsWith("Time: ")){continue;}
			if(sb.length()>0){sb.append('\n');}
			sb.append(line);
		}
		return sb.toString();
	}

	/** Writes eight finite score/label pairs in distinct fit bins. */
	private static Path writePairs(final Path dir) throws Exception{
		final Path in=dir.resolve("pairs.tsv");
		Files.write(in, Arrays.asList(
			"0.10\t0",
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

	/** Writes four one-input vector rows scored by the retained literal network. */
	private static Path writeVectors(final Path dir) throws Exception{
		final Path in=dir.resolve("vectors.tsv");
		Files.write(in, Arrays.asList(
			"#dims\t1\t1",
			"0\t0",
			"1\t0",
			"2\t1",
			"3\t1"
		), StandardCharsets.US_ASCII);
		return in;
	}

	/** Writes a deterministic one-input linear net whose predictions stay finite. */
	private static Path writeNet(final Path dir) throws Exception{
		final Path net=dir.resolve("net.bbnet");
		Files.write(net, Arrays.asList(
			"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
			"#seed 1", "#layers 2", "#dims 1 1", "#edges 1", "#coding decimal",
			"##layer 1", "C2 LINEAR 0.25 0.125"
		), StandardCharsets.US_ASCII);
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

	/** Checks the process report while ignoring only elapsed-time variability. */
	private static void checkReport(final String report, final String label,
			final String writeLine, final int rows){
		check(report.contains("Dataset:   "+label+"\n"),
			"Missing dataset label in report:\n"+report);
		check(report.contains("Rows="+rows+"  mean_label=0.50000  occupied_bins="),
			"Missing row summary in report:\n"+report);
		check(report.contains("log-loss:  raw="), "Missing log-loss line:\n"+report);
		check(report.contains("ECE(4):     raw="), "Missing ECE line:\n"+report);
		check(report.contains("params:    K="), "Missing parameter line:\n"+report);
		check(report.contains("#reliability (binned by calibrated probability)\n"),
			"Missing displayed reliability table:\n"+report);
		check(report.contains("\nTime: "), "Missing elapsed-time line:\n"+report);
		check(!report.contains("NaN"), "Unexpected NaN in report:\n"+report);
		check(!report.contains("Infinity"), "Unexpected infinity in report:\n"+report);
		if(writeLine!=null){
			check(report.contains(writeLine+"\n"), "Missing write report:\n"+report);
		}
		for(String line : withoutTime(report).split("\n")){
			check(!line.toLowerCase(Locale.ROOT).contains("exception"),
				"Unexpected error-looking report line: "+line);
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

	private CalibrateIoSuccessTest(){}
}
