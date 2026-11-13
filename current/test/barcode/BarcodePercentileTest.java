package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import barcode.Barcode;
import barcode.BarcodeCounter;
import barcode.BarcodeStats;
import barcode.CountBarcodes2;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;

/**
 * Finite read-weighted cutoff and opt-in matrix-filter checks.
 * Literal count/map oracles distinguish read mass from distinct-name percentiles,
 * keep cutoff ties, and require the original count/statistics collection to survive.
 * Each CLI case runs once in a fresh fixture JVM. Reflection exposes only the
 * existing CountBarcodes2 process method and its original BarcodeStats instance.
 * @author Brian Bushnell, Shinobu
 */
public final class BarcodePercentileTest{

	/** Requires mode=baseline|fixed, case=helper|off|on|ties|nan, and dir for CLI cases. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, BarcodePercentileTest.class, false);
		try{
			String mode=null, testCase=null, directory=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
				if(a.equals("mode")){mode=b;}
				else if(a.equals("case")){testCase=b;}
				else if(a.equals("dir")){directory=b;}
				else if(!Parser.parseStatic(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or fixed.");
			check(Arrays.asList("helper", "off", "on", "ties", "nan").contains(testCase), "Set a supported case.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "The baseline qualifies the unfinished helper's actual assertion; -ea is required.");
			Shared.setThreads(1);
			final boolean baseline="baseline".equals(mode);
			if(testCase.equals("helper")){helper(baseline);}
			else{
				check(directory!=null, "CLI checks require a fresh dir.");
				cli(baseline, testCase, Files.createDirectory(Paths.get(directory)));
			}
			System.out.println("BARCODE_PERCENTILE mode="+mode+" case="+testCase+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Literal cutoffs exercise weighted bins, ties, histogram overflow and endpoint precision. */
	private static void helper(final boolean baseline){
		cutoff("off", new long[]{1, 1, 8, 90}, 0f, 0, baseline, false);
		cutoff("lower-tail", new long[]{1, 1, 8, 90}, 0.1f, 8, baseline, false);
		cutoff("weighted-median", new long[]{1, 1, 8, 90}, 0.5f, 90, baseline, false);
		cutoff("ties", new long[]{1, 1, 8, 8, 82}, 0.1f, 8, baseline, false);
		cutoff("tiny-fraction", new long[]{1, 1, 8, 90}, Float.MIN_VALUE, 0, baseline, false);
		cutoff("maximum", new long[]{1, 1, 8, 90}, 1f, 90, baseline, false);
		cutoff("dense-overflow-boundary", new long[]{99999, 100000, 200000}, 0.5f, 100000, baseline, false);
		cutoff("long-endpoint", new long[]{9223372036854775806L}, 1f,
				9223372036854775806L, baseline, true);
		if(!baseline){
			invalid("negative-fraction", counts(1), -0.1f, IllegalArgumentException.class, null);
			invalid("large-fraction", counts(1), 1.1f, IllegalArgumentException.class, null);
			invalid("nan-fraction", counts(1), Float.NaN, IllegalArgumentException.class, null);
			invalid("infinite-fraction", counts(1), Float.POSITIVE_INFINITY, IllegalArgumentException.class, null);
			invalid("null-collection", null, 0.5f, IllegalArgumentException.class, null);
			final ArrayList<Barcode> nullEntry=counts(1); nullEntry.add(null);
			invalid("null-entry", nullEntry, 0.5f, IllegalArgumentException.class, null);
			invalid("negative-count", counts(-1, 2), 0.5f, IllegalArgumentException.class, null);
			invalid("empty", counts(), 0.5f, IllegalStateException.class, null);
			invalid("zero-mass", counts(0, 0), 0.5f, IllegalStateException.class, null);
			invalid("overflow", counts(Long.MAX_VALUE, 1), 0.5f, IllegalArgumentException.class, "overflow");
			System.out.println("domain_cases=10");
		}
	}

	/** Baseline requires the prototype assertion, with its computed value or the known long endpoint failure. */
	private static void cutoff(final String label, final long[] values, final float fraction,
			final long expected, final boolean baseline, final boolean endpoint){
		final ArrayList<Barcode> list=counts(values);
		AssertionError failure=null;
		long actual=-1;
		try{actual=BarcodeCounter.barcodeCountPercentile(list, fraction);}
		catch(AssertionError e){failure=e;}
		if(baseline){
			check(failure!=null && frame(failure, "barcode.BarcodeCounter", "barcodeCountPercentile"),
					label+": expected an assertion on the actual unfinished helper path.");
			if(endpoint){
				check(frame(failure, "structures.SuperLongList", "percentileValueBySum")
						&& "1.0, 1, 9223372036854775806".equals(failure.getMessage()),
						label+": expected the rounded cumulative target to exceed the represented total.");
			}else{check(Long.toString(expected).equals(failure.getMessage()), label+": wrong baseline cutoff: "+failure);}
		}else{check(failure==null && actual==expected, label+": expected cutoff "+expected+", observed "+actual+", "+failure);}
		for(int i=0; i<values.length; i++){check(list.get(i).count()==values[i], label+": helper mutated a count.");}
		System.out.println(label+" cutoff="+expected+" prototype_assertion="+baseline);
	}

	/** Checks the new helper's declared invalid-input contract, without qualifying unrelated exceptions. */
	private static void invalid(final String label, final ArrayList<Barcode> values, final float fraction,
			final Class<? extends RuntimeException> type, final String term){
		RuntimeException failure=null;
		try{BarcodeCounter.barcodeCountPercentile(values, fraction);}
		catch(RuntimeException e){failure=e;}
		check(type.isInstance(failure) && frame(failure, "barcode.BarcodeCounter", "barcodeCountPercentile"),
				label+": expected "+type.getSimpleName()+" from the cutoff helper, observed "+failure);
		check(failure.getMessage()!=null && !failure.getMessage().isEmpty(), label+": missing diagnostic.");
		if(term!=null){check(failure.getMessage().toLowerCase(Locale.ROOT).contains(term), label+": missing "+term);}
	}

	/** Invokes the real argument parser and processing path while retaining access to original statistics. */
	private static void cli(final boolean baseline, final String testCase, final Path dir) throws Exception{
		final boolean ties=testCase.equals("ties"), nan=testCase.equals("nan"), off=testCase.equals("off");
		final String fraction=nan ? "NaN" : off ? "0" : "0.1";
		final long[] values=ties ? new long[]{1, 1, 8, 8, 82} : new long[]{1, 1, 8, 90};
		final String input=ties ? "ACGT\t1\nTCGA\t1\nGATC\t8\nCTAG\t8\nTAGC\t82\n" :
				"ACGT\t1\nTCGA\t1\nGATC\t8\nCTAG\t90\n";
		write(dir.resolve("input.tsv"), input);
		write(dir.resolve("expected.txt"), ties ? "ACGT\nTCGA\nGATC\nCTAG\nTAGC\n" : "ACGT\nTCGA\nGATC\nCTAG\n");
		final ArrayList<String> arguments=new ArrayList<String>(Arrays.asList("countsin="+dir.resolve("input.tsv"), "expected="+dir.resolve("expected.txt"),
				"out="+dir.resolve("raw.tsv"), "mapout="+dir.resolve("map.tsv"), "useserver=f", "matrixtype=hdist",
				"mincounta=1", "mincountr=1", "t=1", "matrixthreads=1",
				"ow=f", "pigz=f", "unpigz=f", "simd=f", "printstats=t"));
		if(!off){arguments.add("mincountpercentile="+fraction);}
		final String[] command=arguments.toArray(new String[arguments.size()]);
		write(dir.resolve("cli.args"), String.join("\n", command)+"\n");
		final PrintStream previous=System.err;
		final ByteArrayOutputStream captured=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(captured, true, "UTF-8");
		CountBarcodes2 tool=null;
		Throwable failure=null;
		try{
			System.setErr(capture);
			tool=new CountBarcodes2(command);
			final Method process=CountBarcodes2.class.getDeclaredMethod("process", Timer.class);
			process.setAccessible(true);
			try{process.invoke(tool, new Timer());}
			catch(InvocationTargetException e){failure=e.getCause();}
		}catch(RuntimeException e){failure=e;}
		finally{System.setErr(previous); capture.close();}
		write(dir.resolve("cli.stderr"), new String(captured.toByteArray(), StandardCharsets.UTF_8));
		if(failure!=null){
			try(PrintStream error=new PrintStream(Files.newOutputStream(dir.resolve("failure.txt")), false, "UTF-8")){
				failure.printStackTrace(error);
			}
		}
		if(!baseline && nan){
			check(failure instanceof IllegalArgumentException && failure.getMessage()!=null
					&& failure.getMessage().toLowerCase(Locale.ROOT).contains("mincountpercentile"),
					"NaN CLI input must fail with an explicit option diagnostic, observed "+failure);
			check(!Files.exists(dir.resolve("raw.tsv")) && !Files.exists(dir.resolve("map.tsv")),
					"Invalid percentile must be rejected before output writing.");
			return;
		}
		if(baseline && !off && !nan){
			check(failure instanceof AssertionError && "8".equals(failure.getMessage())
					&& frame(failure, "barcode.BarcodeCounter", "barcodeCountPercentile"),
					"Enabled baseline must fail at the unfinished cutoff8 assertion, observed "+failure);
			check(!Files.exists(dir.resolve("map.tsv")), "Baseline must not write an assignment map after its assertion.");
		}else{check(failure==null, "Unexpected CLI failure: "+failure);}
		final String raw=ties ? "#Barcodes\t100\n#Unique\t5\nTAGC\t82\nCTAG\t8\nGATC\t8\nACGT\t1\nTCGA\t1\n" :
				"#Barcodes\t100\n#Unique\t4\nCTAG\t90\nGATC\t8\nACGT\t1\nTCGA\t1\n";
		check(raw.equals(read(dir.resolve("raw.tsv"))), "Raw counts changed or were filtered.");
		final Field statsField=CountBarcodes2.class.getDeclaredField("bs");
		statsField.setAccessible(true);
		final BarcodeStats stats=(BarcodeStats)statsField.get(tool);
		if(off){
			final Field setting=CountBarcodes2.class.getDeclaredField("minCountPercentile");
			setting.setAccessible(true);
			check(((Number)setting.get(tool)).doubleValue()==0, "The omitted option must remain disabled by default.");
		}
		check(stats.codeMap.size()==values.length, "Percentile filtering mutated the original observation map.");
		for(int i=0; i<values.length; i++){
			check(stats.codeMap.get(NAMES[i]).count()==values[i], "Original observation count changed: "+NAMES[i]);
		}
		final String statistics=stats.toStats("Barcodes");
		check(statistics.startsWith("Barcodes\t100       \t("+values.length+" unique)"), "Raw statistics lost observations.");
		write(dir.resolve("raw.stats"), statistics+"\n");
		if(failure==null){
			final String map=off || nan ? "ACGT\tACGT\t1\t0\nCTAG\tCTAG\t90\t0\nGATC\tGATC\t8\t0\nTCGA\tTCGA\t1\t0\n" :
					ties ? "CTAG\tCTAG\t8\t0\nGATC\tGATC\t8\t0\nTAGC\tTAGC\t82\t0\n" :
					"CTAG\tCTAG\t90\t0\nGATC\tGATC\t8\t0\n";
			check(map.equals(read(dir.resolve("map.tsv"))), "Wrong retained barcode map or count.");
		}
		System.out.println("raw_reads=100 raw_unique="+values.length+" percentile="+fraction+" prototype_assertion="+(failure!=null));
	}

	private static ArrayList<Barcode> counts(final long... values){
		final ArrayList<Barcode> list=new ArrayList<Barcode>(values.length);
		for(int i=0; i<values.length; i++){list.add(new Barcode(NAMES[i%NAMES.length], values[i]));}
		return list;
	}
	private static boolean frame(final Throwable failure, final String type, final String method){
		if(failure==null){return false;}
		for(StackTraceElement element : failure.getStackTrace()){
			if(element.getClassName().equals(type) && element.getMethodName().equals(method)){return true;}
		}
		return false;
	}
	private static void write(final Path file, final String text) throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}
	private static String read(final Path file) throws Exception{return new String(Files.readAllBytes(file), StandardCharsets.UTF_8);}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static final String[] NAMES={"ACGT", "TCGA", "GATC", "CTAG", "TAGC"};
	private BarcodePercentileTest(){}
}
