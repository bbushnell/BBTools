package align2;

import java.io.File;
import java.nio.file.Files;

import dna.Data;

/**
 * Runs BBMapS on two reference sets using BBSplitter's merged-reference arguments.
 * Checks only that both named output files exist and are nonempty. A header-only
 * SAM can satisfy that check; record counts, routing and mapping accuracy require
 * separate tests. Assertions must be enabled for the output checks.
 * Run from a private working/index directory: reference preparation and mapping
 * write files, overwrite is enabled, and the mapper uses four threads. The read
 * fixture must exercise both reference sets. Use a single-output pattern with %,
 * without a paired-file # placeholder; this helper checks one file per set.
 * Mapper and scaffold-prefix settings remain process-global after the run.
 */
public final class BBMapSBBSplitterSmokeTest{

	/** Arguments: first reference, second reference, reads, output pattern, build ID. */
	public static void main(String[] args) throws Exception{
		if(args.length!=5){
			throw new IllegalArgumentException("usage: ref1 ref2 reads output_pattern build");
		}
		final String pattern=args[3];
		final String[] splitterArgs={
			"ref_a="+args[0], "ref_b="+args[1], "in="+args[2],
			"basename="+pattern, "build="+args[4], "overwrite=t", "ordered=t", "t=4"
		};
		final String[] mapperArgs=BBSplitter.processArgs(splitterArgs);
		Data.scaffoldPrefixes=true;
		BBMapS.main(mapperArgs);

		final File outA=new File(pattern.replaceFirst("%", "a"));
		final File outB=new File(pattern.replaceFirst("%", "b"));
		assert(outA.isFile() && outA.length()>0) : outA;
		assert(outB.isFile() && outB.length()>0) : outB;
		System.out.println("PASS BBMapSBBSplitterSmokeTest: out_a="+outA+" ("+Files.size(outA.toPath())
			+" bytes), out_b="+outB+" ("+Files.size(outB.toPath())+" bytes).");
	}
}
