package test.barcode;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

import barcode.BarcodeMappingStats;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite public-source checks for UNKNOWN label ordering and empty barcode reports.
 * Reflection reaches only the nested Label class in public CompareLabels source.
 * No optional or private probability implementation is referenced.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeReportRegressionTest{

	/** Requires mode=baseline|fixed and case=order|mapping; mapping also requires a new dir=path. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, BarcodeReportRegressionTest.class, false);
		try{
			String mode=null, name=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT);
				final String value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){name=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Barcode report checks require -ea.");
			final boolean baseline=mode.equals("baseline");
			final int checks;
			if("order".equals(name)){checks=checkOrder(baseline);}
			else if("mapping".equals(name)){
				if(directory==null){throw new IllegalArgumentException("Mapping case requires a new dir=path.");}
				checks=checkMapping(baseline, Files.createDirectory(Paths.get(directory)));
			}else{throw new IllegalArgumentException("Set case=order or case=mapping.");}
			System.out.println("BARCODE_REPORT case="+name+" mode="+mode+" checks="+checks+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; root runs each case in a fresh JVM. */
	private BarcodeReportRegressionTest(){}

	/** Checks reciprocal signs directly; it does not rely on TimSort throwing for inconsistent input. */
	private static int checkOrder(final boolean baseline) throws ReflectiveOperationException{
		final Labels labels=new Labels();
		final Object unknown=labels.create("UNKNOWN", 1), normal=labels.create("ZETA", 2);
		check(labels.compare(unknown, normal)<0, "UNKNOWN must sort before a higher-count ordinary label");
		final int reverse=labels.compare(normal, unknown);
		check(baseline ? reverse<0 : reverse>0, "Reverse UNKNOWN comparison violates the requested mode: "+reverse);
		final Object tiedUnknown=labels.create("unknown", 2), alpha=labels.create("ALPHA", 2);
		check(labels.compare(tiedUnknown, alpha)<0, "UNKNOWN recognition must remain case-insensitive at equal counts");
		final int tiedReverse=labels.compare(alpha, tiedUnknown);
		check(baseline ? tiedReverse<0 : tiedReverse>0, "Equal-count UNKNOWN reverse order is wrong: "+tiedReverse);
		final Object lower=labels.create("BETA", 1), tied=labels.create("BETA", 2);
		check(labels.compare(alpha, lower)<0, "Higher ordinary count must sort first");
		check(labels.compare(lower, alpha)>0, "Lower ordinary count must sort last");
		check(labels.compare(alpha, tied)<0, "Ordinary equal-count names must sort lexically");
		check(labels.compare(tied, alpha)>0, "Ordinary lexical order must be reciprocal");
		check(labels.compare(alpha, alpha)==0, "Ordinary self comparison must be zero");
		check(labels.compare(unknown, unknown)==0, "UNKNOWN self comparison must be zero");
		System.out.println("unknown_reverse_sign="+Integer.signum(reverse)+" tied_reverse_sign="+Integer.signum(tiedReverse));
		return 10;
	}

	/** Checks a two-barcode control and the empty report's literal headers and absence of data rows. */
	private static int checkMapping(final boolean baseline, final Path dir) throws IOException{
		final BarcodeMappingStats normal=new BarcodeMappingStats();
		normal.incrementCodeMap("ACGT", 4);
		normal.incrementSourceMap("ACGT", "refA", 3);
		normal.incrementSourceMap("ACGT", "refB", 1);
		normal.incrementCodeMap("TGCA", 2);
		normal.incrementSourceMap("TGCA", "refC", 2);
		final String columns="#Barcode\tSource\tCount\tFraction\n";
		final String expected="#Reads\t6\n"+columns+
				"ACGT\trefA\t3\t0.750000\nACGT\trefB\t1\t0.250000\nTGCA\trefC\t2\t1\n";
		final Path control=dir.resolve("normal.tsv");
		normal.writeStats(control.toString(), false);
		check(!normal.errorState, "Normal mapping report recorded a writer error");
		checkText("Normal mapping report", expected, read(control));
		final BarcodeMappingStats empty=new BarcodeMappingStats();
		final Path target=dir.resolve("empty.tsv");
		boolean threw=false;
		try{empty.writeStats(target.toString(), false);}
		catch(final NullPointerException e){threw=true;}
		check(threw==baseline, "Empty report NullPointerException: expected="+baseline+", observed="+threw);
		if(baseline){check(!Files.exists(target), "Empty baseline must fail before creating output");}
		else{
			check(!empty.errorState, "Empty mapping report recorded a writer error");
			checkText("Empty mapping report", "#Reads\t0\n"+columns, read(target));
		}
		System.out.println("empty_report_null_exception="+threw);
		return baseline ? 4 : 5;
	}

	/** Reflective view of the public-source nested label, using its own comparison method. */
	private static final class Labels{
		Labels() throws ReflectiveOperationException{
			final Class<?> type=Class.forName("barcode.CompareLabels$Label");
			constructor=type.getDeclaredConstructor(String.class);
			agreement=type.getDeclaredField("aa");
			compare=type.getDeclaredMethod("compareTo", type);
			constructor.setAccessible(true);
			agreement.setAccessible(true);
			compare.setAccessible(true);
		}
		/** Uses the agreement counter to supply an exact total without changing production visibility. */
		Object create(final String name, final long count) throws ReflectiveOperationException{
			final Object label=constructor.newInstance(name);
			agreement.setLong(label, count);
			return label;
		}
		/** Returns the production comparator's sign-bearing result. */
		int compare(final Object left, final Object right) throws ReflectiveOperationException{
			return (Integer)compare.invoke(left, right);
		}
		final Constructor<?> constructor;
		final Field agreement;
		final Method compare;
	}

	/** Reads the fixture's small completed output files. */
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	/** Requires exact independent text, including integer and fractional number formatting. */
	private static void checkText(final String label, final String expected, final String actual){
		if(!expected.equals(actual)){throw new AssertionError(label+" expected:\n"+expected+"observed:\n"+actual);}
	}
	/** A failed oracle always throws, independently of the JVM assertion switch. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
