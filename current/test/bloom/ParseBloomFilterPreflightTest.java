package test.bloom;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import bloom.ParseBloomFilter;

/** Checks output-path validation without starting readers or writers.
 * Each case uses a fresh directory; constructor rejection must leave input and
 * existing output bytes intact. Does not exercise processing or stream cleanup.
 * Run in a fresh JVM: parser construction changes process-wide I/O settings.
 * @author Brian Bushnell, Shinobu */
public final class ParseBloomFilterPreflightTest{

	/** Takes one scratch-directory path and reports every case before failing. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Paths.get(args[0]);
		Files.createDirectory(root);
		int failures=0;
		failures+=exercise(root, "distinct", "new", true, false, false);
		failures+=exercise(root, "invalid-input", "input", true, false, true);
		failures+=exercise(root, "invalid-summary", "summary", true, false, true);
		failures+=exercise(root, "existing-protected", "existing", false, false, true);
		failures+=exercise(root, "existing-overwrite", "existing", true, false, false);
		failures+=exercise(root, "invalid-disabled", "null", true, false, false);
		failures+=exercise(root, "summary-input", "new", true, true, true);
		failures+=exercise(root, "equals=in-path", "new", true, false, false);
		if(failures!=0){throw new AssertionError("Output preflight failed "+failures+" of 8 cases.");}
		System.out.println("PARSE_BLOOM_PREFLIGHT_TEST_PASS cases=8");
	}

	/** Constructs only the parser; fixtures are never passed to its processing loop. */
	private static int exercise(final Path root, final String name, final String invalidKind,
			final boolean overwrite, final boolean summaryIsInput, final boolean reject) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path input=dir.resolve("input.txt"), existing=dir.resolve("existing.txt");
		final Path summary=(summaryIsInput ? input : dir.resolve("summary.txt"));
		Files.write(input, INPUT);
		Files.write(existing, EXISTING);
		final String invalid=(invalidKind.equals("input") ? input.toString() :
				invalidKind.equals("summary") ? summary.toString() :
				invalidKind.equals("existing") ? existing.toString() :
				invalidKind.equals("null") ? "null" : dir.resolve("invalid.txt").toString());
		boolean rejected=false;
		String error="none";
		try{
			new ParseBloomFilter(new String[]{"in="+input, "out="+summary,
					"invalid="+invalid, "ow="+overwrite, "t=1", "pigz=f", "unpigz=f"});
		}catch(RuntimeException | AssertionError e){
			rejected=true;
			error=e.getClass().getSimpleName();
		}
		final boolean intact=Arrays.equals(INPUT, Files.readAllBytes(input)) &&
				Arrays.equals(EXISTING, Files.readAllBytes(existing));
		final boolean pass=rejected==reject && intact;
		System.out.println("case="+name+" expectedRejected="+reject+" rejected="+rejected+
				" exception="+error+" originalBytesIntact="+intact+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Literal independent fixtures; neither is a generated expected parser result. */
	private static final byte[] INPUT="# source\nother\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] EXISTING="keep existing output\n".getBytes(StandardCharsets.US_ASCII);
}
