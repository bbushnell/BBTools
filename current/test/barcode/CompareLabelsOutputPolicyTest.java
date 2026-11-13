package test.barcode;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import barcode.CompareLabels;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks each CompareLabels report's output policy against literal text and a sentinel.
 * One report/scenario per fresh JVM; baseline records the ignored overwrite/append flags.
 * Does not exercise writer shutdown failures or assertion-disabled behavior.
 * @author Brian Bushnell, Sucrose
 */
public final class CompareLabelsOutputPolicyTest{

	/** Requires mode=baseline|fixed, report=summary|labels, scenario=fresh|reject|replace|append, dir=NEW. */
	public static void main(final String[] args){
		int exit=1;
		try{
			run(args);
			exit=0;
		}catch(final Throwable failure){
			failure.printStackTrace(System.err);
		}
		//A production failure need not clean up all I/O threads; contain each case in its own JVM.
		System.exit(exit);
	}

	/** Parses fixture selectors separately from the production CLI and requires assertions. */
	private static void run(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CompareLabelsOutputPolicyTest.class, false);
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
			check("summary".equals(report) || "labels".equals(report), "Set report=summary or report=labels.");
			check("fresh".equals(scenario) || "reject".equals(scenario) || "replace".equals(scenario) || "append".equals(scenario),
					"Set scenario=fresh, reject, replace or append.");
			check(directory!=null, "Set dir to a new directory whose parent exists.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Output-policy checks require the production canWrite assertion (-ea).");
			final Path dir=Files.createDirectory(Paths.get(directory));
			checkPolicy(dir, mode.equals("baseline"), report, scenario);
			System.out.println("COMPARE_LABELS_OUTPUT mode="+mode+" report="+report+" scenario="+scenario+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Checks exact replacement/append bytes or a qualified rejection with the sentinel intact. */
	private static void checkPolicy(final Path dir, final boolean baseline, final String report, final String scenario) throws Exception{
		final Path input=dir.resolve("input.fq"), output=dir.resolve("report.tsv");
		Files.write(input, INPUT.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(!scenario.equals("fresh")){
			Files.write(output, SENTINEL.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		}
		final boolean labels=report.equals("labels"), append=scenario.equals("append"), overwrite=scenario.equals("replace");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input,
				"out="+(labels ? "null" : output.toString()), "ow="+overwrite, "append="+append,
				"delimiter=tab", "interleaved=f", "qin=33", "qout=33", "trimreadcomments=f",
				"t=1", "simd=f", "pigz=f", "unpigz=f"));
		if(labels){args.add("labelstats="+output);}
		Throwable observed=null;
		try{CompareLabels.main(args.toArray(new String[args.size()]));}
		catch(final Throwable failure){observed=failure;}
		final boolean reject=!baseline && scenario.equals("reject");
		if(reject){
			check(observed!=null, "ow=false must reject the existing report without replacing its bytes.");
			final String message=observed.getMessage();
			//Require the standard output-policy diagnostic, not an unrelated input/parse failure.
			check((observed instanceof AssertionError || observed instanceof RuntimeException) &&
					hasCallerFrame(observed) && message!=null && message.contains(output.toString()) && message.contains("overwrite=false"),
					"Expected an output-policy failure naming the report and overwrite=false; observed "+observed);
			checkText(SENTINEL, read(output));
			System.out.println("POLICY rejected=true sentinel_preserved=true exception="+observed.getClass().getName());
		}else{
			if(observed!=null){throw new AssertionError("Expected successful report generation for "+scenario, observed);}
			final String normal=labels ? LABEL_REPORT : SUMMARY_REPORT;
			checkText(!baseline && append ? SENTINEL+normal : normal, read(output));
			System.out.println("POLICY rejected=false appended="+(!baseline && append)+" report_bytes="+Files.size(output));
		}
	}

	/** Attributes a caught policy failure to this production invocation without pinning line numbers. */
	private static boolean hasCallerFrame(final Throwable failure){
		for(final StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals("barcode.CompareLabels")){return true;}
		}
		return false;
	}

	/** Reads only the fixture's bounded report. */
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	/** Includes all headers, column ordering and final newlines in the independent oracle. */
	private static void checkText(final String expected, final String actual){
		check(expected.equals(actual), "Report expected:\n"+expected+"observed:\n"+actual);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private CompareLabelsOutputPolicyTest(){}

	private static final String INPUT="@one\tA\tA\nACGT\n+\nIIII\n";
	private static final String SENTINEL="KEEP\n";
	private static final String SUMMARY_REPORT="#Labels\t1\n#Valid\t1\t1\n#TermsPerRead\t3\n"+
			"#RelYield1\t1\n#RelYield2\t1\n#AbsYield1\t1\n#AbsYield2\t1\n"+
			"#Contam1_PPM\t0\n#Contam2_PPM\t0\n#Metric\tCount\tRate\tPPM\n"+
			"AACount\t1\t1\t1000000\nUUCount\t0\t0\t0\n"+
			"AUCount\t0\t0\t0\nUACount\t0\t0\t0\nABCount\t0\t0\t0\n\n";
	private static final String LABEL_REPORT="#Labels\t1\n#Valid\t1\t1\n#TermsPerRead\t3\n"+
			"#Name\tCount\tCount1\tCount2\tAA\tAU\tUA\tAB\tBA\tYield1\tYield2\tContam1\tContam2\n"+
			"A\t1\t1\t1\t1\t0\t0\t0\t0\t1\t1\t0\t0\n";
}
