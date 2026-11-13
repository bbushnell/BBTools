package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

import barcode.CompareLabels;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks empty normalized reports, ordinary reports and runs with no report requested.
 * Baseline decimal zeros do not establish finite arithmetic: the formatter can hide NaN.
 * Each combination requires a fresh JVM and a new scratch directory.
 * @author Brian Bushnell, Sucrose
 */
public final class CompareLabelsEmptyInputTest{

	/** Requires mode=baseline|fixed, report=summary|labels|both|none, scenario=empty|nonempty, dir=NEW. */
	public static void main(final String[] args){
		int exit=1;
		try{
			run(args);
			exit=0;
		}catch(final Throwable failure){
			failure.printStackTrace(System.err);
		}
		//Contain unexpected production failures that might leave I/O threads alive.
		System.exit(exit);
	}

	/** Parses fixture selectors and verifies the assertion-enabled execution contract. */
	private static void run(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CompareLabelsEmptyInputTest.class, false);
		try{
			String mode=null, report=null, scenario=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("report")){report=value;}
				else if(key.equals("scenario")){scenario=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("summary".equals(report) || "labels".equals(report) || "both".equals(report) || "none".equals(report),
					"Set report=summary, labels, both or none.");
			check("empty".equals(scenario) || "nonempty".equals(scenario), "Set scenario=empty or scenario=nonempty.");
			check(directory!=null, "Set dir to a new directory whose parent exists.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Run the production and fixture checks with -ea.");
			final Path dir=Files.createDirectory(Paths.get(directory));
			checkReports(dir, mode.equals("baseline"), report, scenario.equals("empty"));
			System.out.println("COMPARE_LABELS_EMPTY mode="+mode+" report="+report+" scenario="+scenario+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Runs the real CLI and checks either exact reports or a qualified pre-writer rejection. */
	private static void checkReports(final Path dir, final boolean baseline, final String report, final boolean empty) throws IOException{
		final boolean summary=report.equals("summary") || report.equals("both");
		final boolean labels=report.equals("labels") || report.equals("both");
		final Path input=dir.resolve("input.fq"), summaryFile=dir.resolve("summary.tsv"), labelsFile=dir.resolve("labels.tsv");
		Files.write(input, (empty ? "" : INPUT).getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final String[] args={"in="+input,
				"out="+(summary ? summaryFile.toString() : "null"), "labelstats="+(labels ? labelsFile.toString() : "null"),
				"ow=f", "append=f", "delimiter=tab", "interleaved=f", "qin=33", "qout=33",
				"trimreadcomments=f", "t=1", "simd=f", "pigz=f", "unpigz=f"};
		Throwable observed=null;
		try{CompareLabels.main(args);}
		catch(final Throwable failure){observed=failure;}
		final boolean reject=!baseline && empty && (summary || labels);
		if(reject){
			check(observed!=null && observed.getClass()==IllegalStateException.class &&
					EMPTY_DIAGNOSTIC.equals(observed.getMessage()) && hasProcessFrame(observed),
					"Expected the normalized-report empty-input guard, not an unrelated failure: "+observed);
			check(!Files.exists(summaryFile) && !Files.exists(labelsFile),
					"The empty-input guard must run before either report writer opens a destination.");
			System.out.println("EMPTY_REPORT rejected=true report_files_created=false exception="+observed.getClass().getName());
		}else{
			if(observed!=null){throw new AssertionError("Expected successful report/no-report processing", observed);}
			checkReport(summaryFile, summary, empty ? EMPTY_SUMMARY : SUMMARY);
			checkReport(labelsFile, labels, empty ? EMPTY_LABELS : LABELS);
			System.out.println("EMPTY_REPORT rejected=false input_records="+(empty ? 0 : 1)+" requested_reports="+((summary ? 1 : 0)+(labels ? 1 : 0)));
		}
	}

	/** Requires a literal complete report, or absence when that report was not requested. */
	private static void checkReport(final Path path, final boolean requested, final String expected) throws IOException{
		if(!requested){check(!Files.exists(path), "Unrequested report was created: "+path); return;}
		final String actual=new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
		check(expected.equals(actual), "Report "+path+" expected:\n"+expected+"observed:\n"+actual);
	}
	private static boolean hasProcessFrame(final Throwable failure){
		for(final StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals("barcode.CompareLabels") && frame.getMethodName().equals("process")){return true;}
		}
		return false;
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private CompareLabelsEmptyInputTest(){}

	private static final String EMPTY_DIAGNOSTIC="CompareLabels cannot write normalized reports without observed labels.";
	private static final String INPUT="@one\tA\tA\nACGT\n+\nIIII\n";
	private static final String SUMMARY="#Labels\t1\n#Valid\t1\t1\n#TermsPerRead\t3\n"+
			"#RelYield1\t1\n#RelYield2\t1\n#AbsYield1\t1\n#AbsYield2\t1\n"+
			"#Contam1_PPM\t0\n#Contam2_PPM\t0\n#Metric\tCount\tRate\tPPM\n"+
			"AACount\t1\t1\t1000000\nUUCount\t0\t0\t0\n"+
			"AUCount\t0\t0\t0\nUACount\t0\t0\t0\nABCount\t0\t0\t0\n\n";
	//These baseline zero strings expose formatter output, not a claim that the underlying ratios are finite.
	private static final String EMPTY_SUMMARY="#Labels\t0\n#Valid\t0\t0\n#TermsPerRead\t0\n"+
			"#RelYield1\t0\n#RelYield2\t0\n#AbsYield1\t0\n#AbsYield2\t0\n"+
			"#Contam1_PPM\t0\n#Contam2_PPM\t0\n#Metric\tCount\tRate\tPPM\n"+
			"AACount\t0\t0\t0.00\nUUCount\t0\t0\t0.00\n"+
			"AUCount\t0\t0\t0.00\nUACount\t0\t0\t0.00\nABCount\t0\t0\t0.00\n\n";
	private static final String LABEL_HEADER="#Name\tCount\tCount1\tCount2\tAA\tAU\tUA\tAB\tBA\tYield1\tYield2\tContam1\tContam2\n";
	private static final String LABELS="#Labels\t1\n#Valid\t1\t1\n#TermsPerRead\t3\n"+LABEL_HEADER+
			"A\t1\t1\t1\t1\t0\t0\t0\t0\t1\t1\t0\t0\n";
	private static final String EMPTY_LABELS="#Labels\t0\n#Valid\t0\t0.00000\n#TermsPerRead\t0\n"+LABEL_HEADER;
}
