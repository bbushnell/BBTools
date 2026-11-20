package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import jgi.CallPeaks;
import parse.Parse;

/**
 * Checks the static peak-report helper's existing nonpositive-k fallback.
 * Literal counts and GC totals describe a 31-mer histogram with 50% GC.
 * Baseline mode records the inconsistent raw-k reports; fixed mode requires
 * all three calls to report the effective constructor value, 31.
 * @author Shinobu
 */
public final class CallPeaksKNormalizationTest{

	/** Requires phase=baseline|fixed and a new outdir; creates three small reports. */
	public static void main(final String[] args) throws IOException{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
		}
		if((!"baseline".equals(phase) && !"fixed".equals(phase)) || directory==null){
			throw new IllegalArgumentException("Set phase=baseline|fixed and a new outdir.");
		}
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final boolean baseline=phase.equals("baseline");
		int mismatches=0;
		mismatches+=checkCase(root, 31, baseline);
		mismatches+=checkCase(root, 0, baseline);
		mismatches+=checkCase(root, -1, baseline);
		assert(mismatches==(baseline ? 2 : 0)) :
			"The two nonpositive arguments must share the constructor's 31-mer report after repair; mismatches="+mismatches;
		System.out.println("PASS phase="+phase+" cases=3 mismatches="+mismatches);
	}

	/** Calls the production API once and compares its complete report to literal expectations. */
	private static int checkCase(final Path root, final int k, final boolean baseline) throws IOException{
		final long[] counts={0, 2, 6, 2, 0}, gc={0, 31, 93, 31, 0};
		final long[] countsBefore=counts.clone(), gcBefore=gc.clone();
		final Path report=root.resolve("k"+k+".txt");
		final boolean error=CallPeaks.printPeaks(counts, gc, report.toString(), false,
			2, 5, 3, 2, Integer.MAX_VALUE, 10, k, -1, false, 0.1, null);
		if(error){throw new AssertionError("Writer reported an error for k="+k);}
		assert(Arrays.equals(counts, countsBefore) && Arrays.equals(gc, gcBefore)) :
			"The static report API borrows histogram arrays and must leave both inputs unchanged; k="+k;
		final String actual=new String(Files.readAllBytes(report), StandardCharsets.UTF_8);
		final String expected=expectedReport(k, baseline);
		if(!actual.equals(expected)){
			throw new AssertionError("Unexpected complete report for k="+k+" at "+report+"\n"+actual);
		}
		final boolean mismatch=!actual.equals(expectedReport(31, false));
		System.out.println("k="+k+" writer_error=false mismatch="+mismatch);
		return mismatch ? 1 : 0;
	}

	/** Independent 50%-GC control, with only the known raw-k baseline differences. */
	private static String expectedReport(final int k, final boolean baseline){
		assert(k==31 || k==0 || k==-1) : "Literal oracle covers only the 31-mer control and its two fallback arguments: "+k;
		final String gc=baseline && k<=0 ? (k==0 ? "1.000" : "0.000") : "0.500";
		final String rowGc=baseline && k<=0 ? (k==0 ? "1" : "0") : "0.500";//ByteBuilder prints integral doubles as integers.
		return (baseline && k<=0 ? "" : "#k\t31\n")+
			"#unique_kmers\t10\n#error_kmers\t0\n#genomic_kmers\t10\n#main_peak\t2\n"+
			"#genome_size_in_peaks\t10\n#genome_size\t12\n#main_peak_gc\t"+gc+
			"\n#gc_content_in_peaks\t"+gc+"\n#gc_content\t"+gc+
			"\n#haploid_genome_size\t12\n#fold_coverage\t2\n#haploid_fold_coverage\t2\n"+
			"#ploidy\t1\n#percent_repeat_in_peaks\t0.000\n#percent_repeat\t0.000\n"+
			"#start\tcenter\tstop\tmax\tvolume\tgc\n1\t2\t4\t6\t10\t"+rowGc+"\n";
	}
}
