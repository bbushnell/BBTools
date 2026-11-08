package test.bloom;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilterCorrectorWrapper;

/** Constructor-only output-path contracts for the correction wrapper (BLM215).
 * Builds tiny real reference tables but never starts processing/output writers.
 * Existing bytes must survive every case. Run in a fresh JVM because constructors
 * alter shared parser/counting settings. Does not test exceptional stream cleanup.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorPreflightTest{

	/** Captures eight rejection/control cases under the supplied fresh directory. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Files.createDirectory(Paths.get(args[0].substring(7)));
		int failures=0;
		failures+=exercise(root, "distinct", false);
		failures+=exercise(root, "orphan-rejected-mate", true);
		failures+=exercise(root, "rejected-pair", false);
		failures+=exercise(root, "quality-input", true);
		failures+=exercise(root, "quality-primary", true);
		failures+=exercise(root, "quality-rejected", true);
		failures+=exercise(root, "quality-existing", true);
		failures+=exercise(root, "primary-input", true);
		if(failures>0){throw new AssertionError("Corrector preflight: "+failures+" failures of 8 cases.");}
		System.out.println("BLOOM_CORRECTOR_PREFLIGHT_TEST_PASS cases=8");
	}

	/** Constructs the real wrapper and checks rejection plus all original bytes. */
	private static int exercise(final Path root, final String name, final boolean expectReject) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path input=dir.resolve("input.fq"), mate=dir.resolve("mate.fq"), reference=dir.resolve("reference.fa");
		final Path existing=dir.resolve("existing.qual"), rejected=dir.resolve("rejected.fq");
		final Path primary=(name.equals("primary-input") ? input : dir.resolve("output.fq"));
		final Path quality=(name.equals("quality-input") ? input : name.equals("quality-primary") ? primary :
				name.equals("quality-rejected") ? rejected : name.equals("quality-existing") ? existing : dir.resolve("output.qual"));
		Files.write(input, INPUT);
		Files.write(mate, MATE);
		Files.write(reference, REFERENCE);
		Files.write(existing, EXISTING);
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag : new String[]{"in="+input, "ref="+reference, "out="+primary, "qfout1="+quality,
				"ow="+!name.equals("quality-existing"), "k=15", "ksmall=15", "bits=4", "hashes=1",
				"cells=10007", "seed=1", "minprob=0", "t=1", "interleaved=f", "pigz=f", "unpigz=f", "qin=33",
				"ecc=f", "merge=f", "ecco=f", "markerrors=f", "mincount=0", "tossjunk=f", "loglog=f"}){
			flags.add(flag);
		}
		if(!name.equals("orphan-rejected-mate")){flags.add("outbad1="+rejected);}
		if(name.equals("orphan-rejected-mate") || name.equals("rejected-pair")){
			flags.add("in2="+mate);
			flags.add("outbad2="+dir.resolve("rejected-mate.fq"));
		}
		boolean rejectedInput=false;
		String error="none";
		try{new BloomFilterCorrectorWrapper(flags.toArray(new String[flags.size()]));}
		catch(RuntimeException | AssertionError failure){
			rejectedInput=true;
			error=failure.getClass().getSimpleName();
			System.err.println("EXCEPTION_CAPTURE case="+name);
			failure.printStackTrace(System.err);
		}
		final boolean intact=Arrays.equals(INPUT, Files.readAllBytes(input)) && Arrays.equals(MATE, Files.readAllBytes(mate)) &&
				Arrays.equals(REFERENCE, Files.readAllBytes(reference)) && Arrays.equals(EXISTING, Files.readAllBytes(existing));
		final boolean pass=rejectedInput==expectReject && intact;
		System.out.println("case="+name+" expectedRejected="+expectReject+" rejected="+rejectedInput+
				" exception="+error+" originalBytesIntact="+intact+" pass="+pass);
		return pass ? 0 : 1;
	}

	private static final byte[] INPUT="@pair/1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] MATE="@pair/2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REFERENCE=">reference\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] EXISTING="retain this preexisting quality output\n".getBytes(StandardCharsets.US_ASCII);
}
