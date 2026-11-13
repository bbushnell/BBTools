package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Locale;

import barcode.BarcodeStats;
import jgi.DemuxByName2;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks barcode-specific expected-name normalization and exact two-read routing.
 * Generic loader and genuine tab-delimited dual names remain unaffected controls.
 * Each case/source uses one fresh JVM and a nonexistent filesystem root.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeExpectedNamesTest{

	/** Requires mode=baseline|fixed, case=annotated|plain|rc|dual|tab|generic, source=file|literal, dir=NEW. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, BarcodeExpectedNamesTest.class, false);
		try{
			String mode=null, testCase=null, source=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("source")){source=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or fixed.");
			check(Arrays.asList("annotated", "plain", "rc", "dual", "tab", "generic").contains(testCase), "Unknown case: "+testCase);
			check("file".equals(source) || "literal".equals(source), "Set source=file or literal.");
			check(directory!=null, "Set dir to a nonexistent output directory.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Expected-name cases require -ea for production invariants.");
			Shared.setThreads(1);
			final Path dir=Files.createDirectory(Paths.get(directory));
			run(mode, testCase, source.equals("file"), dir);
			System.out.println("BARCODE_EXPECTED_NAMES mode="+mode+" case="+testCase+" source="+source+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Literal independent names distinguish count annotations, transforms and real index separators. */
	private static void run(final String mode, final String testCase, final boolean file, final Path dir) throws IOException{
		final boolean annotated=testCase.equals("annotated"), rc=testCase.equals("rc"), dual=testCase.equals("dual");
		final String[] names=testCase.equals("tab") ? new String[] {"AAGC\tTTGA", "TTGA\tAAGC"} :
				testCase.equals("generic") ? new String[] {"sample_one\tmetadata", "sample two\t12"} :
				dual ? new String[] {"AAGC+TTGA\t7", "TTGA+AAGC\t12"} :
				annotated || rc ? new String[] {"AAGC\t7", "TTGA\t12"} : new String[] {"AAGC", "TTGA"};
		final Path namesFile=dir.resolve("names.tsv");
		write(namesFile, names[0]+"\n"+names[1]+"\n");
		final byte delimiter=(byte)(testCase.equals("tab") ? '\t' : dual ? '+' : 0);
		final String[] transformed=rc ? new String[] {"GCTT", "TCAA"} :
				dual ? new String[] {"AAGC+TTGA", "TTGA+AAGC"} : names;
		//The legacy overload also serves generic names. Its existing no-transform semantics must survive.
		final LinkedHashSet<String> input=new LinkedHashSet<String>();
		if(file){input.add(namesFile.toString());}else{input.addAll(Arrays.asList(names));}
		final LinkedHashSet<String> result=BarcodeStats.loadBarcodeSet(input, delimiter, rc, false, false);
		check(result==input && Arrays.equals(result.toArray(new String[0]), transformed),
				"Generic expansion, source order or existing delimiter/RC behavior changed: "+result);
		if(testCase.equals("tab")){
			check("GCTT\tTCAA".equals(BarcodeStats.fixBarcode(names[0], '\t', true, true)),
					"A genuine tab separator must survive independent reverse complementation of both indices.");
			return;//Illumina barcode headers use '+', so do not claim tab-header CLI coverage.
		}
		if(testCase.equals("generic")){return;}
		final String[] expected=rc ? transformed : dual ? transformed : new String[] {"AAGC", "TTGA"};
		//The ordinary barcode-list API already strips annotations; keep that unaffected control.
		final BarcodeStats stats=new BarcodeStats(delimiter, dual ? 2 : 1, null);
		stats.loadBarcodeList(file ? namesFile.toString() : names[0]+","+names[1], delimiter, rc, false);
		check(stats.expectedCodeList.size()==2 && stats.expectedCodeList.get(0).name.equals(expected[0]) &&
				stats.expectedCodeList.get(1).name.equals(expected[1]), "Barcode-list normalization changed.");
		final String first=record(1, expected[0]), second=record(2, expected[1]);
		//Extensionless input forces a real format/header probe despite explicit single-end selection.
		write(dir.resolve("reads"), first+second);
		DemuxByName2.main(new String[] {"in="+dir.resolve("reads"), "names="+
				(file ? namesFile.toString() : names[0]+","+names[1]), "barcode=t", "pcr=f", "hdist=0",
				"rc1="+rc, "rc2=f", "out="+dir.resolve("matched_%.fq"), "outu="+dir.resolve("unmatched.fq"),
				"mapout="+dir.resolve("assignment.tsv"), "interleaved=f", "qin=33", "qout=33",
				"minreads=1", "ow=f", "append=f", "t=1", "pigz=f", "unpigz=f"});
		final boolean misses=mode.equals("baseline") && annotated;
		final String[] mapNames=misses ? names : expected;
		final String line1=mapNames[0]+"\t"+mapNames[0]+"\n", line2=mapNames[1]+"\t"+mapNames[1]+"\n";
		final String map=read(dir.resolve("assignment.tsv"));
		check(map.equals(line1+line2) || map.equals(line2+line1), "Assignment keys/values differ from the literal two-name oracle: "+map);
		check(readOrEmpty(dir.resolve("matched_"+expected[0]+".fq")).equals(misses ? "" : first), "First read routed incorrectly.");
		check(readOrEmpty(dir.resolve("matched_"+expected[1]+".fq")).equals(misses ? "" : second), "Second read routed incorrectly.");
		check(readOrEmpty(dir.resolve("unmatched.fq")).equals(misses ? first+second : ""), "Unexpected unmatched reads or changed FASTQ records.");
		System.out.println("matched_reads="+(misses ? 0 : 2)+" unmatched_reads="+(misses ? 2 : 0)+" assignment_names=2");
	}

	/** Two single reads with distinct coordinates and exactly known sequence/quality bytes. */
	private static String record(final int number, final String barcode){
		return "@INST:1:FLOW:1:1101:"+number+":1 1:N:0:"+barcode+"\nACGT\n+\nIIII\n";
	}
	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	/** No reads may produce no matched file or an empty file; any payload must be exactly accounted for. */
	private static String readOrEmpty(final Path path) throws IOException{return Files.exists(path) ? read(path) : "";}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private BarcodeExpectedNamesTest(){}
}
