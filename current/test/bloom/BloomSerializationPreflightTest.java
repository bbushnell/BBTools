package test.bloom;

import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterWrapper;
import bloom.KCountArray7MTA;
import fileIO.ReadWrite;
import shared.Shared;

/** Checks serialization/reference path protection before wrapper processing.
 * Uses only newly created scratch files. Constructor-started object writes are
 * drained before inspecting bytes; no read-processing output streams are started.
 * Run in a fresh JVM because construction changes shared settings.
 * @author Brian Bushnell, Shinobu */
public final class BloomSerializationPreflightTest{

	/** Takes a new scratch directory and retains all nine case outcomes. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path root=Files.createDirectory(Paths.get(args[0]));
		Shared.setThreads(1);
		BloomFilter.OVERRIDE_CELLS=10007;
		KCountArray7MTA.setSeed(1);
		int failures=0;
		for(String name : new String[]{"fresh-write", "write-input", "write-reference", "write-output",
				"write-protected", "fresh-load", "load-output", "load-counts", "counts-reference"}){
			failures+=exercise(root, name);
		}
		if(failures!=0){throw new AssertionError("Serialization preflight failed "+failures+" of 9 cases.");}
		System.out.println("BLOOM_SERIALIZATION_PREFLIGHT_TEST_PASS cases=9");
	}

	/** Keeps reference/input/snapshot/protected bytes independent of the rejection result. */
	private static int exercise(final Path root, final String name) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final Path input=dir.resolve("input.fq"), ref=dir.resolve("ref.fa");
		final Path snapshot=dir.resolve("loaded.ser"), protectedFile=dir.resolve("protected.ser");
		Files.write(input, INPUT);
		Files.write(ref, REFERENCE);
		Files.write(protectedFile, PROTECTED);
		final BloomFilter initial=new BloomFilter(null, null, null, 3, 3, 4, 1, 1, true, false, false, 0.01f, false);
		initial.filter.increment(42L, 4);//Canonical CCC is GGG, packed key42.
		try(ObjectOutputStream out=new ObjectOutputStream(Files.newOutputStream(snapshot))){out.writeObject(initial);}
		final byte[] snapshotBytes=Files.readAllBytes(snapshot);
		final boolean load=name.equals("fresh-load") || name.startsWith("load-");
		final boolean reject=!name.startsWith("fresh-");
		final Path output=(name.equals("load-output") ? snapshot : dir.resolve("output.fq"));
		final Path counts=(name.equals("load-counts") ? snapshot : name.equals("counts-reference") ? ref : dir.resolve("counts.fa"));
		final Path serialOut=(name.equals("write-input") ? input : name.equals("write-reference") ? ref :
				name.equals("write-output") ? output : name.equals("write-protected") ? protectedFile : dir.resolve("written.ser"));
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag : new String[]{"in="+input, "ref="+ref, "out="+output, "outc="+counts,
				"ow="+!name.equals("write-protected"), "k=3", "minhits=1", "bits=4", "hashes=1",
				"cells=10007", "seed=1", "minprob=0", "t=1", "interleaved=f", "pigz=f", "unpigz=f", "qin=33"}){
			flags.add(flag);
		}
		if(load){flags.add("serialin="+snapshot);}
		else if(!name.equals("counts-reference")){flags.add("serialout="+serialOut);}
		boolean rejected=false;
		String error="none";
		try{new BloomFilterWrapper(flags.toArray(new String[flags.size()]));}
		catch(RuntimeException | AssertionError e){rejected=true; error=e.getClass().getSimpleName();}
		finally{ReadWrite.waitForWritingToFinish();}
		final boolean intact=Arrays.equals(INPUT, Files.readAllBytes(input)) &&
				Arrays.equals(REFERENCE, Files.readAllBytes(ref)) &&
				Arrays.equals(PROTECTED, Files.readAllBytes(protectedFile)) &&
				Arrays.equals(snapshotBytes, Files.readAllBytes(snapshot));
		boolean wroteValid=true;
		if(name.equals("fresh-write")){
			wroteValid=Files.exists(serialOut);
			if(wroteValid){
				final BloomFilter saved=ReadWrite.read(BloomFilter.class, serialOut.toString(), false);
				wroteValid=saved.getCount(42L)>0;
			}
		}
		final boolean pass=rejected==reject && intact && wroteValid;
		System.out.println("case="+name+" expectedRejected="+reject+" rejected="+rejected+
				" exception="+error+" originalBytesIntact="+intact+" validWrite="+wroteValid+" pass="+pass);
		return pass ? 0 : 1;
	}

	private static final byte[] INPUT="@read\nCCCCCC\n+\nIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REFERENCE=">reference\nCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] PROTECTED="preserve this existing file\n".getBytes(StandardCharsets.US_ASCII);
}
