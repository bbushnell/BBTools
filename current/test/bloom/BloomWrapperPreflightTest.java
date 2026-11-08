package test.bloom;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilterWrapper;

/** Checks wrapper path validation with tiny real reference filters.
 * Constructs wrappers but never calls process or requests serialization, so no
 * output writer is started. Run in a fresh JVM because construction changes globals.
 * Does not test serialized-filter settings, stream cleanup or filtering accuracy.
 * @author Brian Bushnell, Shinobu */
public final class BloomWrapperPreflightTest{

	/** Takes one new scratch-directory path; reports all cases before failing. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Files.createDirectory(Paths.get(args[0]));
		int failures=0;
		failures+=exercise(root, "distinct", false);
		failures+=exercise(root, "count-input", true);
		failures+=exercise(root, "count-unmatched", true);
		failures+=exercise(root, "count-matched", true);
		failures+=exercise(root, "count-existing", true);
		failures+=exercise(root, "orphan-mate", true);
		failures+=exercise(root, "matched-pair", false);
		failures+=exercise(root, "primary-input", true);
		failures+=exercise(root, "qfout-distinct", false);
		failures+=exercise(root, "qfout-input", true);
		failures+=exercise(root, "qfout-unmatched", true);
		failures+=exercise(root, "qfout-matched", true);
		failures+=exercise(root, "qfout-count", true);
		failures+=exercise(root, "qfout-reference", true);
		failures+=exercise(root, "qfout-existing", true);
		if(failures!=0){throw new AssertionError("Wrapper preflight failed "+failures+" of 15 cases.");}
		System.out.println("BLOOM_WRAPPER_PREFLIGHT_TEST_PASS cases=15");
	}

	/** Exercises constructor validation only, preserving each original file for byte comparison. */
	private static int exercise(final Path root, final String name, final boolean reject) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path input=dir.resolve("input.fq"), ref=dir.resolve("ref.fa");
		final Path existing=dir.resolve("existing.fa"), matched=dir.resolve("matched.fq");
		final Path unmatched=(name.equals("primary-input") ? input : dir.resolve("unmatched.fq"));
		final Path counts=(name.equals("count-input") ? input : name.equals("count-unmatched") ? unmatched :
				name.equals("count-matched") ? matched : name.equals("count-existing") ? existing : dir.resolve("counts.fa"));
		final Path qfout=(name.equals("qfout-input") ? input :
				name.equals("qfout-unmatched") ? unmatched :
				name.equals("qfout-matched") ? matched :
				name.equals("qfout-count") ? counts :
				name.equals("qfout-reference") ? ref :
				name.equals("qfout-existing") ? existing : dir.resolve("qualities.fq"));
		Files.write(input, INPUT);
		Files.write(ref, REFERENCE);
		Files.write(existing, EXISTING);
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag : new String[]{"in="+input, "ref="+ref, "out="+unmatched, "outc="+counts,
				"ow="+!(name.equals("count-existing") || name.equals("qfout-existing")),
				"k=3", "minhits=1", "bits=4", "hashes=1",
				"cells=10007", "seed=1", "minprob=0", "t=1", "interleaved=f", "pigz=f", "unpigz=f", "qin=33"}){
			flags.add(flag);
		}
		if(!name.equals("orphan-mate")){flags.add("outm="+matched);}
		if(name.equals("orphan-mate") || name.equals("matched-pair")){flags.add("outm2="+dir.resolve("mate2.fq"));}
		if(name.startsWith("qfout-")){flags.add("qfout="+qfout);}
		boolean rejected=false;
		String error="none";
		try{new BloomFilterWrapper(flags.toArray(new String[flags.size()]));}
		catch(RuntimeException | AssertionError e){rejected=true; error=e.getClass().getSimpleName();}
		final boolean intact=Arrays.equals(INPUT, Files.readAllBytes(input)) &&
				Arrays.equals(REFERENCE, Files.readAllBytes(ref)) && Arrays.equals(EXISTING, Files.readAllBytes(existing));
		final boolean pass=rejected==reject && intact;
		System.out.println("case="+name+" expectedRejected="+reject+" rejected="+rejected+
				" exception="+error+" originalBytesIntact="+intact+" pass="+pass);
		return pass ? 0 : 1;
	}

	private static final byte[] INPUT="@input\nCCCCCC\n+\nIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REFERENCE=">reference\nCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] EXISTING="keep existing output\n".getBytes(StandardCharsets.US_ASCII);
}
