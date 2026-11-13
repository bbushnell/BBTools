package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

import barcode.AnalyzeBarcodes;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite expected-count validation regression for AnalyzeBarcodes.
 * Numeric mode proves expected=&lt;number&gt; can select the same top observed pairs
 * as an explicit expected file. File mode preserves the existing accepted path.
 * @author Brian Bushnell, Jean
 */
public final class AnalyzeBarcodesExpectedTest{

	/** Requires mode=baseline|fixed, case=numeric|file and a new dir=path. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, AnalyzeBarcodesExpectedTest.class, false);
		try{
			String mode=null, testCase=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			if(!"numeric".equals(testCase) && !"file".equals(testCase)){
				throw new IllegalArgumentException("Set case=numeric or case=file.");
			}
			if(directory==null){throw new IllegalArgumentException("Set dir to a nonexistent output directory.");}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("AnalyzeBarcodes expected regression requires -ea.");}

			final boolean baseline=mode.equals("baseline");
			final Path dir=Files.createDirectory(Paths.get(directory));
			writeInput(dir);
			final int checks="numeric".equals(testCase) ? checkNumeric(baseline, dir) : checkExpectedFile(dir);
			System.out.println("ANALYZE_BARCODES_EXPECTED mode="+mode+" case="+testCase+" checks="+checks+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each case owns one fresh process and filesystem root. */
	private AnalyzeBarcodesExpectedTest(){}

	/** Writes two high-count expected pairs and two swapped unexpected-but-valid pairs. */
	private static void writeInput(final Path dir) throws IOException{
		Files.write(dir.resolve("counts.tsv"), ("#code\tcount\n"+
				"ACGT+TGCA\t8\n"+
				"TGCA+ACGT\t7\n"+
				"ACGT+ACGT\t3\n"+
				"TGCA+TGCA\t2\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		Files.write(dir.resolve("expected.tsv"), ("ACGT+TGCA\n"+
				"TGCA+ACGT\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Baseline must reject expected=2 in validation; fixed mode must match the file-selected control. */
	private static int checkNumeric(final boolean baseline, final Path dir) throws IOException{
		boolean failedValidation=false;
		try{
			AnalyzeBarcodes.main(command(dir, "2"));
		}catch(final AssertionError e){
			if(!baseline){throw e;}
			failedValidation=fromValidateParams(e);
			if(!failedValidation){throw e;}
		}
		if(baseline){
			check(failedValidation, "Numeric baseline must fail in AnalyzeBarcodes.validateParams.");
			check(!Files.exists(dir.resolve("left.tsv")), "Baseline must fail before writing left.tsv.");
			check(!Files.exists(dir.resolve("right.tsv")), "Baseline must fail before writing right.tsv.");
			check(!Files.exists(dir.resolve("pair.tsv")), "Baseline must fail before writing pair.tsv.");
			System.out.println("numeric_validate_params_assertion="+failedValidation);
			return 4;
		}
		checkOutputs(dir);
		return 3;
	}

	/** Explicitly lists the two top pairs selected by the equivalent numeric case. */
	private static int checkExpectedFile(final Path dir) throws IOException{
		AnalyzeBarcodes.main(command(dir, dir.resolve("expected.tsv").toString()));
		checkOutputs(dir);
		return 3;
	}

	/** Builds a complete command with explicit reports to avoid default-path collisions. */
	private static String[] command(final Path dir, final String expected){
		return new String[] {
				"in="+dir.resolve("counts.tsv"),
				"expected="+expected,
				"left="+dir.resolve("left.tsv"),
				"right="+dir.resolve("right.tsv"),
				"pair="+dir.resolve("pair.tsv"),
				"ow=f",
				"t=1",
				"pigz=f",
				"unpigz=f"
		};
	}

	/** Checks the literal reports from the chosen two expected barcode pairs. */
	private static void checkOutputs(final Path dir) throws IOException{
		checkText("left", "#Name\tTotal\tBadPair\n"+
				"ACGT\t11\t3\n"+
				"TGCA\t9\t2\n", read(dir.resolve("left.tsv")));
		checkText("right", "#Name\tTotal\tBadPair\n"+
				"ACGT\t10\t3\n"+
				"TGCA\t10\t2\n", read(dir.resolve("right.tsv")));
		checkText("pair", "#Name\tExpectedCount\tActualCount\n"+
				"ACGT+ACGT\t0\t3\n"+
				"TGCA+TGCA\t0\t2\n", read(dir.resolve("pair.tsv")));
	}

	/** Confirms the failing assertion is the over-strict expected-file validation. */
	private static boolean fromValidateParams(final AssertionError e){
		for(final StackTraceElement frame : e.getStackTrace()){
			if(frame.getClassName().equals("barcode.AnalyzeBarcodes") && frame.getMethodName().equals("validateParams")){
				return true;
			}
		}
		return false;
	}

	/** Reads only bounded fixture output. */
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Requires exact independent text, including headers and final newlines. */
	private static void checkText(final String label, final String expected, final String observed){
		if(!expected.equals(observed)){
			throw new AssertionError(label+" expected:\n"+expected+"observed:\n"+observed);
		}
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
