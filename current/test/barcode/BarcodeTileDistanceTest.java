package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

import barcode.Barcode;
import barcode.NovaDemux;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks equal-core tile suffixes in left-distance calculation and legacy map reports.
 * The first argument continues to define the comparison boundary. Each legacy case
 * owns a fresh JVM and scratch directory; no sequence output or server is requested.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeTileDistanceTest{

	/** Requires mode=baseline|fixed, case=api|legacy-extra|legacy-control; legacy cases need a new dir. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, BarcodeTileDistanceTest.class, false);
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
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("api".equals(testCase) || "legacy-extra".equals(testCase) || "legacy-control".equals(testCase),
					"Set case=api, legacy-extra or legacy-control.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Tile-distance checks require -ea for the suffix-domain assertion.");
			Shared.setThreads(1);
			final boolean baseline="baseline".equals(mode);
			if("api".equals(testCase)){checkApi(baseline);}
			else{
				check(directory!=null, "Legacy cases require a new dir=path.");
				checkLegacy(Files.createDirectory(Paths.get(directory)), baseline, "legacy-extra".equals(testCase));
			}
			System.out.println("BARCODE_TILE_DISTANCE mode="+mode+" case="+testCase+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Literal distances preserve ordinary, reversed and first-argument separator behavior. */
	private static void checkApi(final boolean baseline){
		checkDistance("suffix-first", "ACGT1101", "ACGT", 0, baseline);
		checkDistance("suffix-first-mismatch", "ACGA1101", "ACGT", 1, baseline);
		checkDistance("suffix-second", "ACGT", "ACGT1101", 0, false);
		checkDistance("suffix-second-mismatch", "ACGT", "ACGA1101", 1, false);
		checkDistance("plain", "ACGT", "ACGT", 0, false);
		checkDistance("plain-mismatch", "ACGA", "ACGT", 1, false);
		checkDistance("dual-left", "AAAA+CCCC", "AAAT+GGGG", 1, false);
		//These intentionally differ: only the first argument determines where comparison stops.
		checkDistance("later-first-boundary", "AA+AA", "A+AAA", 1, false);
		checkDistance("earlier-first-boundary", "A+AAA", "AA+AA", 0, false);
		System.out.println("api_cases=9");
	}

	/** A caught baseline failure must originate in the actual two-string left-distance helper. */
	private static void checkDistance(final String label, final String a, final String b,
			final int expected, final boolean expectedFailure){
		StringIndexOutOfBoundsException failure=null;
		int actual=-1;
		try{actual=Barcode.hdistL(a, b);}
		catch(final StringIndexOutOfBoundsException e){failure=e;}
		check((failure!=null)==expectedFailure, label+": expected bounds failure="+expectedFailure+", observed="+failure);
		if(failure!=null){check(stackContains(failure, "barcode.Barcode", "hdistL"), label+": wrong failure path.");}
		else{check(actual==expected, label+": expected distance="+expected+", observed="+actual);}
		System.out.println(label+" bounds_failure="+(failure!=null)+" distance="+actual);
	}

	/** An unused raw suffix row must not break reports for the same matched/unknown reads as the control. */
	private static void checkLegacy(final Path dir, final boolean baseline, final boolean extra) throws IOException{
		write(dir.resolve("input.fq"), "@INST:1:FC:1:1101:100:100 1:N:0:ACGT\nACGT\n+\nIIII\n"+
				"@INST:1:FC:1:1101:100:101 1:N:0:TGCA\nTGCA\n+\nIIII\n");
		final String map="ACGT\tACGT\n"+(extra ? "ACGT1101\tACGT\n" : "");
		write(dir.resolve("map.tsv"), map);
		final Path legacy=Files.createDirectory(dir.resolve("legacy"));
		StringIndexOutOfBoundsException failure=null;
		try{
			NovaDemux.main(new String[] {"in="+dir.resolve("input.fq"), "expected=ACGT",
					"mapin="+dir.resolve("map.tsv"), "legacypath="+legacy, "matrixtype=hdist",
					"server=f", "interleaved=f", "qin=33", "ow=f", "t=1", "pigz=f", "unpigz=f",
					"loglog=f", "trimreadcomments=f"});
		}catch(final StringIndexOutOfBoundsException e){failure=e;}
		finally{waitForLegacyWriters(legacy);}
		final boolean expectedFailure=baseline && extra;
		check((failure!=null)==expectedFailure, "Legacy unused row: expected bounds failure="+expectedFailure+", observed="+failure);
		if(failure!=null){
			check(stackContains(failure, "barcode.Barcode", "hdistL") &&
					stackContains(failure, "barcode.LegacyFileWriter", "writeDemultiplexStats"), "Wrong legacy failure path.");
			check(!Files.exists(legacy.resolve("Demultiplex_Stats.csv")), "Baseline must fail before creating demultiplex stats.");
		}else{
			checkText(legacy.resolve("Demultiplex_Stats.csv"),
					"Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads\n"+
					"0,,ACGT,1,1,0,0,0.5000,1,0,0\n0,Undetermined,,1,1,0,0,0.5000,1,0,0\n");
		}
		checkText(legacy.resolve("Top_Unknown_Barcodes.csv"),
				"Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads\n0,TGCA,,1,1,0.500000\n");
		checkText(legacy.resolve("Quality_Metrics.csv"),
				"Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30\n"+
				"0,Undetermined,ACGT,,1,4,4,160,40,1\n0,Undetermined,,,1,4,4,160,40,1\n");
		checkText(dir.resolve("map.tsv"), map);
		System.out.println("legacy_unused_row="+extra+" qualified_bounds_failure="+(failure!=null));
	}

	/** Joins report writers even after the expected baseline exception; no new writer is created afterward. */
	private static void waitForLegacyWriters(final Path legacy){
		final String prefix="BSW for "+legacy;
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.toString().startsWith(prefix)){
				final ByteStreamWriter writer=(ByteStreamWriter)thread;
				writer.waitForFinish();
				check(!writer.errorState, "Legacy report writer failed: "+writer);
			}
		}
	}

	/** Keeps fixture files exclusive; a previous result is never overwritten. */
	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Checks literal bytes, including headers and final newlines. */
	private static void checkText(final Path path, final String expected) throws IOException{
		final String actual=new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
		check(expected.equals(actual), path+": expected:\n"+expected+"observed:\n"+actual);
	}

	/** Qualifies the production origin rather than accepting any caught bounds error. */
	private static boolean stackContains(final Throwable failure, final String type, final String method){
		for(final StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals(type) && frame.getMethodName().equals(method)){return true;}
		}
		return false;
	}

	/** Keeps literal oracle failures visible independently of downstream assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private BarcodeTileDistanceTest(){}
}
