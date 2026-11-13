package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Locale;

import barcode.Barcode;
import barcode.BarcodeCounter;
import barcode.CompareLabels;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite public regressions for count-writer delegation, requested homopolymer
 * bases and label-summary PPM units. Each case runs in a fresh JVM with assertions.
 * Baseline mode requires the original defect; fixed mode requires its repair.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeCoreRegressionTest{

	/** Requires mode=baseline|fixed and case=counts|homopolymer|labels; I/O cases require a new dir=path. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, BarcodeCoreRegressionTest.class, false);
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
			if(!assertions){throw new IllegalStateException("Barcode regression checks require -ea.");}
			final boolean baseline=mode.equals("baseline");
			final int checks;
			if("homopolymer".equals(name)){checks=checkHomopolymers(baseline);}
			else if("counts".equals(name) || "labels".equals(name)){
				if(directory==null){throw new IllegalArgumentException("I/O cases require a new dir=path.");}
				final Path dir=Files.createDirectory(Paths.get(directory));
				checks=name.equals("counts") ? checkCounts(baseline, dir) : checkLabels(baseline, dir);
			}else{throw new IllegalArgumentException("Set case=counts, homopolymer or labels.");}
			System.out.println("BARCODE_CORE case="+name+" mode="+mode+" checks="+checks+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each process exercises one bounded contract. */
	private BarcodeCoreRegressionTest(){}

	/** Compares the String writer control and the FileFormat overload on the same two literal records. */
	private static int checkCounts(final boolean baseline, final Path dir) throws IOException{
		final ArrayList<Barcode> counts=new ArrayList<Barcode>();
		counts.add(new Barcode("TGCA", 1));
		counts.add(new Barcode("ACGT", 3));
		final String expected="#Barcodes\t4\n#Unique\t2\nACGT\t3\nTGCA\t1\n";
		final Path control=dir.resolve("string-counts.tsv");
		check(!BarcodeCounter.writeCounts(counts, 1, control.toString(), true, false, false),
				"String writer reported an output error");
		checkText("String writer", expected, read(control));
		final Path target=dir.resolve("format-counts.tsv");
		final FileFormat ff=FileFormat.testOutput(target.toString(), FileFormat.TXT, null, true, false, false, false);
		boolean overflow=false;
		try{
			check(!BarcodeCounter.writeCounts(counts, 1, ff, true, false, false),
					"FileFormat writer reported an output error");
		}catch(final StackOverflowError e){overflow=true;}
		check(overflow==baseline, "FileFormat recursion: expected StackOverflowError="+baseline+", observed="+overflow);
		if(baseline){check(!Files.exists(target), "Recursive baseline must not have written its target");}
		else{checkText("FileFormat writer", expected, read(target));}
		System.out.println("counts\tFileFormat_StackOverflowError="+overflow);
		return baseline ? 4 : 5;
	}

	/** Uses literal base-specific truth tables plus mixed and delimited controls. */
	private static int checkHomopolymers(final boolean baseline){
		final String[] names={"AAAA", "CCCC", "GGGG", "TTTT"};
		final char[] requested={'A', 'C', 'G', 'T'};
		final boolean[][] expected={{true, false, false, false}, {false, true, false, false},
				{false, false, true, false}, {false, false, false, true}};
		int checks=0;
		for(int i=0; i<names.length; i++){
			final Barcode barcode=new Barcode(names[i]);
			for(int j=0; j<requested.length; j++){
				check(barcode.isHomopolymer(requested[j])==(baseline || expected[i][j]),
						"Instance homopolymer mismatch: "+names[i]+", requested="+requested[j]);
				check(Barcode.isHomopolymer(names[i], requested[j])==expected[i][j],
						"Static homopolymer control mismatch: "+names[i]+", requested="+requested[j]);
				checks+=2;
			}
		}
		check(!new Barcode("AGGA").isHomopolymer('G'), "Mixed bases must not be a G homopolymer");
		check(new Barcode("GG-GG").isHomopolymer('G'), "Existing delimiter tolerance must remain");
		check(new Barcode("AA+AA").isHomopolymer('G')==baseline, "Delimited A must follow requested-base policy");
		check(new Barcode("CCCC").isHomopolymer(), "Unqualified homopolymer control must remain true");
		return checks+4;
	}

	/** Runs the public CLI on five mixed labels and an all-agreeing control, checking complete summary text. */
	private static int checkLabels(final boolean baseline, final Path dir) throws IOException{
		final String mixed="@r1\tA\tA\nACGT\n+\nIIII\n"+
				"@r2\tB\tB\nACGT\n+\nIIII\n"+
				"@r3\tC\tC\nACGT\n+\nIIII\n"+
				"@r4\tA\tB\nACGT\n+\nIIII\n"+
				"@r5\tUNKNOWN\tUNKNOWN\nACGT\n+\nIIII\n";
		// ByteBuilder.append(double,int) emits exact integers without decimals; other values retain the requested places.
		// summarize ends with a newline and ByteStreamWriter.println(ByteBuilder) adds the second final newline.
		final String ppm=baseline ? "0.25" : "250000";
		final String expected="#Labels\t5\n#Valid\t5\t1\n#TermsPerRead\t3\n"+
				"#RelYield1\t0.75000\n#RelYield2\t0.75000\n#AbsYield1\t0.80000\n#AbsYield2\t0.80000\n"+
				"#Contam1_PPM\t"+ppm+"\n#Contam2_PPM\t"+ppm+"\n"+
				"#Metric\tCount\tRate\tPPM\nAACount\t3\t0.60000\t600000\n"+
				"UUCount\t1\t0.20000\t200000\nAUCount\t0\t0\t0\n"+
				"UACount\t0\t0\t0\nABCount\t1\t0.20000\t200000\n\n";
		checkText("Mixed label summary", expected, runLabels(dir, "mixed", mixed));
		final String control="#Labels\t1\n#Valid\t1\t1\n#TermsPerRead\t3\n"+
				"#RelYield1\t1\n#RelYield2\t1\n#AbsYield1\t1\n#AbsYield2\t1\n"+
				"#Contam1_PPM\t0\n#Contam2_PPM\t0\n#Metric\tCount\tRate\tPPM\n"+
				"AACount\t1\t1\t1000000\nUUCount\t0\t0\t0\n"+
				"AUCount\t0\t0\t0\nUACount\t0\t0\t0\nABCount\t0\t0\t0\n\n";
		checkText("All-agreeing label summary", control, runLabels(dir, "agree", "@r1\tA\tA\nACGT\n+\nIIII\n"));
		return 2;
	}

	/** Creates one new input and output pair and invokes CompareLabels with explicit read/header settings. */
	private static String runLabels(final Path dir, final String stem, final String input) throws IOException{
		final Path in=dir.resolve(stem+".fq"), out=dir.resolve(stem+".tsv");
		Files.write(in, input.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		CompareLabels.main(new String[] {"in="+in, "out="+out, "interleaved=f", "qin=33",
				"trimreadcomments=f", "t=1", "simd=f", "pigz=f", "unpigz=f"});
		return read(out);
	}

	/** Reads only the fixture's bounded output files. */
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Requires an independent literal text oracle, including headers, units and final newlines. */
	private static void checkText(final String label, final String expected, final String actual){
		if(!expected.equals(actual)){
			throw new AssertionError(label+" expected:\n"+expected+"observed:\n"+actual);
		}
	}

	/** A false predicate always fails, independently of the JVM assertion switch. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}
}
