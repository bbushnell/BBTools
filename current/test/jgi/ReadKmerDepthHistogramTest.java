package test.jgi;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLongArray;

import jgi.ReadKmerDepthDistribution;
import parse.Parse;
import parse.PreParser;
import shared.Shared;

/**
 * Checks histogram conservation through the real depth-distribution entry point.
 * Each invocation requires a fresh JVM because the target maintains shared static state.
 * Sets the target's public histogram capacities before invoking main; the histlen option
 * alone controls the printed range and does not shrink its default allocation.
 * Input and histogram files must not already exist in the supplied directory.
 * @author Shinobu
 */
public final class ReadKmerDepthHistogramTest{

	/** Runs one literal clipped or unclipped case in a fresh JVM.
	 * @param args Required case=clipped|unclipped and dir=existing-directory options
	 * @throws IOException If creating the owned input or inspecting output fails */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, ReadKmerDepthHistogramTest.class, false);
		try{
			run(pp.args, pp.outstream);
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** Creates a literal input, runs the target, and checks independent bin/count oracles. */
	private static void run(final String[] args, final PrintStream out) throws IOException{
		String caseName=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("case")){caseName=pair[1];}
			else if(key.equals("dir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown histogram-test option: "+arg);}
		}
		require("clipped".equals(caseName) || "unclipped".equals(caseName),
				"Choose case=clipped or case=unclipped.");
		require(directory!=null && !directory.isEmpty(), "Supply dir=existing-directory.");
		final boolean clipped=caseName.equals("clipped");
		final String bases=clipped ? "AAAAAAAAAAAAAAAAAAAAAAAA" : "AAAAAAAAAAAA";
		final int expectedPositions=clipped ? 20 : 8;
		final int expectedBin=clipped ? 15 : 8;
		verify(bases.length()-K+1==expectedPositions,
				"The literal input must contain the independently specified number of k-mer windows.");
		final Path dir=Paths.get(directory).toAbsolutePath().normalize();
		require(Files.isDirectory(dir), "dir must name an existing directory: "+dir);
		final Path input=dir.resolve("input.fa"), histogramFile=dir.resolve("histogram.tsv");
		require(!Files.exists(histogramFile), "Histogram output already exists: "+histogramFile);
		Files.write(input, (">"+caseName+"\n"+bases+"\n").getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);

		ReadKmerDepthDistribution.HIST_LEN=HISTOGRAM_LENGTH;
		ReadKmerDepthDistribution.HIST_LEN_PRINT=HISTOGRAM_LENGTH;
		ReadKmerDepthDistribution.main(new String[]{"in="+input, "hist="+histogramFile,
				"k=5", "cells=4096", "hashes=1", "threads=1", "minq=0", "minprob=0",
				"fixspikes=f", "histlen=15", "canonical=t", "rdk=f", "ow=f", "simd=f"});

		final AtomicLongArray histogram=ReadKmerDepthDistribution.histogram_total;
		verify(histogram!=null && histogram.length()==HISTOGRAM_LENGTH,
				"The fixture must exercise a 16-bin internal histogram, not merely a short printed range.");
		verify(ReadKmerDepthDistribution.HIST_LEN_PRINT==HISTOGRAM_LENGTH
				&& ReadKmerDepthDistribution.THREAD_HIST_LEN==HISTOGRAM_LENGTH,
				"Printed and thread-local histogram limits must agree with the fixture's 16 bins.");
		long positions=0;
		for(int bin=0; bin<histogram.length(); bin++){positions+=histogram.get(bin);}
		out.println("Observed histogram: case="+caseName+", positions="+positions
				+", bin="+expectedBin+", binCount="+histogram.get(expectedBin));
		verify(positions==expectedPositions, "Histogram must conserve all "+expectedPositions
				+" valid k-mer positions; case="+caseName+", observed="+positions);
		for(int bin=0; bin<histogram.length(); bin++){
			final long expected=(bin==expectedBin ? expectedPositions : 0);
			verify(histogram.get(bin)==expected, "Identical all-A windows must share the expected depth bin; case="
					+caseName+", bin="+bin+", expected="+expected+", observed="+histogram.get(bin));
		}
		final String expectedText="#Depth\tRaw_Count\tUnique_Kmers\n"+expectedBin+"\t"+expectedPositions+"\t1\n";
		verify(Arrays.equals(Files.readAllBytes(histogramFile), expectedText.getBytes(StandardCharsets.US_ASCII)),
				"The written histogram must match the literal bin/count oracle, including its single unique key.");
		out.println("READKMERDEPTH_HISTOGRAM_PASS case="+caseName+" positions="+positions+" bin="+expectedBin);
	}

	/** Rejects invalid fixture configuration before invoking the target. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Enforces regression oracles even if a caller disables JVM assertions. */
	private static void verify(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Canonical short keys require k greater than three. */
	private static final int K=5;
	/** Small public capacity makes clipping reachable with a finite literal input. */
	private static final int HISTOGRAM_LENGTH=16;
}
