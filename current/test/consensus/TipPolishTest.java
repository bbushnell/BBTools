package test.consensus;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Tests standalone publication, input preservation and fail-closed orchestration.
 * @author Brian Bushnell, Fischl */
public final class TipPolishTest {

	/** Runs real child launchers; the test directory is retained for inspection. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Pass the BBTools installation directory.");
		final Path code=Paths.get(args[0]).toAbsolutePath();
		final Path root=Files.createTempDirectory("tippolish-test-");
		final Path assembly=root.resolve("input with spaces.fa"), reads=root.resolve("reads.fq");
		write(assembly, ">short\nACGTAGCTAGCATCGATCGA\n>ambiguous\nNNNNACGTNNNN\n");
		write(reads, "@one\nACGTAGCTAGCATCGATCGA\n+\nIIIIIIIIIIIIIIIIIIII\n");
		final byte[] original=Files.readAllBytes(assembly), originalReads=Files.readAllBytes(reads);
		final Path output=root.resolve("output with spaces.fa"), report=root.resolve("edits.tsv");
		run(code, root, "no_context", true, "assembly="+assembly, "in="+reads, "out="+output, "report="+report);
		check(Arrays.equals(original, Files.readAllBytes(output)), "No-context output must preserve all records and bases.");
		check(Files.readAllLines(report).size()==1, "No-context edit ledger must contain only its header.");
		run(code, root, "existing", false, "assembly="+assembly, "in="+reads, "out="+output);
		run(code, root, "overwrite", true, "assembly="+assembly, "input1="+reads, "out1="+output,
				"report="+report, "ow=t", "threads=auto", "threads=1", "maxoverlap=0.2k");
		check(Arrays.equals(original, Files.readAllBytes(output)), "Explicit overwrite must publish the completed assembly.");
		run(code, root, "alias", false, "assembly="+assembly, "in="+reads, "out="+assembly);
		run(code, root, "overwrite_alias", false, "assembly="+assembly, "in="+reads, "out="+assembly, "ow=t");
		run(code, root, "append", false, "assembly="+assembly, "in="+reads, "out="+root.resolve("append.fa"), "app=t");
		final Path same=root.resolve("same");
		run(code, root, "same_outputs", false, "assembly="+assembly, "in="+reads, "out="+same, "report="+same);
		check(!Files.exists(same), "Aliasing outputs must fail before publication.");
		run(code, root, "pair_conflict", false, "assembly="+assembly, "in="+reads, "in2="+reads,
				"interleaved=t", "out="+root.resolve("conflict.fa"));
		run(code, root, "unknown", false, "assembly="+assembly, "in="+reads, "out="+root.resolve("unknown.fa"), "relax=t");
		run(code, root, "bounds", false, "assembly="+assembly, "in="+reads, "out="+root.resolve("bounds.fa"), "minoverlap=32");
		final Path bad=root.resolve("duplicate.fa");
		write(bad, ">same\nACGT\n>same\nAGCT\n");
		run(code, root, "duplicate", false, "assembly="+bad, "in="+reads, "out="+root.resolve("duplicate.out.fa"));
		write(root.resolve("whitespace.fa"), ">name description\nACGT\n");
		run(code, root, "whitespace", false, "assembly="+root.resolve("whitespace.fa"), "in="+reads, "out="+root.resolve("whitespace.out.fa"));
		final Path config=root.resolve("recipe.config"), compressed=root.resolve("compressed.fa.gz");
		write(config, "assembly="+assembly+"\nin="+reads+"\nout="+compressed+"\ninterleaved=f\nminoverlap=33\nmaxoverlap=200\n");
		run(code, root, "config_gzip", true, "config="+config);
		final ArrayList<byte[]> lines=ByteFile.toLines(compressed.toString());
		check(lines.size()==4 && new String(lines.get(1), StandardCharsets.US_ASCII).equals("ACGTAGCTAGCATCGATCGA"),
				"Compressed output must contain FASTA, not uncompressed bytes with a gz suffix.");
		final Path fake=root.resolve("failing tools");
		Files.createDirectories(fake.resolve("dev"));
		for(String tool : new String[]{"dev/tipcontextpanel.sh", "dev/tipcontextguard.sh", "dev/kmerwindowdepth.sh",
				"bbmap.sh", "samstreamer.sh", "consensus.sh"}){write(fake.resolve(tool), "#!/bin/bash\nexit 17\n");}
		final Path failed=root.resolve("failed.fa");
		run(code, root, "child_failure", false, "bbtools="+fake, "assembly="+assembly, "in="+reads, "out="+failed);
		check(!Files.exists(failed) && !Files.exists(Paths.get(failed+".edits.tsv")), "Child failure must publish neither output.");
		check(Arrays.equals(original, Files.readAllBytes(assembly)) && Arrays.equals(originalReads, Files.readAllBytes(reads)),
				"Every success/failure must leave inputs byte-identical.");
		final Path linked=root.resolve("launcher links");
		Files.createDirectory(linked);
		Files.createSymbolicLink(linked.resolve("tippolish.sh"), code.resolve("tippolish.sh"));
		run(linked, root, "symlink", true, "assembly="+assembly, "in="+reads, "out="+root.resolve("linked.fa"));
		check(Arrays.equals(original, Files.readAllBytes(root.resolve("linked.fa"))), "Symlink launch must resolve the real installation.");
		System.out.println("TIP_POLISH_TEST_PASS checks="+checks+" dir="+root);
	}

	/** Each invocation has independent static state and a retained failure/success log. */
	private static void run(final Path code, final Path root, final String name, final boolean success, final String... args) throws Exception{
		assert(code!=null && root!=null) : "Tests must execute a concrete installation and retain diagnostics.";
		final List<String> command=new ArrayList<String>(Arrays.asList("bash", code.resolve("tippolish.sh").toString(),
				"-Xmx256m", "t=1", "temp="+root));
		command.addAll(Arrays.asList(args));
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(root.resolve(name+".log").toFile()).start();
		child.getOutputStream().close();
		check((child.waitFor()==0)==success, "Unexpected status for "+name+"; see "+root.resolve(name+".log"));
	}

	/** Fixture data uses the native buffered writer and checks asynchronous completion. */
	private static void write(final Path path, final String text){
		assert(path!=null && text!=null) : "Fixture paths and literal content must be explicit.";
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start(); writer.print(new ByteBuilder(text));
		check(!writer.poisonAndWait(), "Could not write fixture: "+path);
	}

	/** Test expectations remain fatal with assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}
	private static int checks;
}
