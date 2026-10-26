package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/** Native CLI tests for scaffold defaults, publication and input preservation.
 * @author Brian Bushnell, Fischl */
public final class TadpoleScaffoldTest {

	/** All generated data are tiny deterministic fixtures, not benchmark libraries. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Supply the tested installation path.");
		final TadpoleScaffoldTest test=new TadpoleScaffoldTest(Paths.get(args[0]).toAbsolutePath());
		test.run();
		System.out.println("TADPOLE_SCAFFOLD_PASS checks="+checks+" dir="+test.root);
	}

	/** Creates one repeat-free fixture and three independent supporting records. */
	private TadpoleScaffoldTest(final Path code_) throws Exception{
		code=code_;
		root=Files.createTempDirectory("tadpole-scaffold-test-");
		final Random random=new Random(971);
		final StringBuilder sequence=new StringBuilder();
		for(int i=0; i<750; i++){sequence.append("ACGT".charAt(random.nextInt(4)));}
		reference=sequence.toString();
		reads=root.resolve("reads.fa");
		Files.write(reads, Arrays.asList(">r0", reference, ">r1", reference, ">r2", reference), StandardCharsets.US_ASCII);
		contigs=root.resolve("contigs.fa");
		Files.write(contigs, Arrays.asList(">left", reference.substring(0, 300), ">right", reference.substring(340, 640)), StandardCharsets.US_ASCII);
	}

	/** Exercises real standalone spelling and native single/multi-K dispatch, not a fake assembler. */
	private void run() throws Exception{
		for(boolean fill : new boolean[]{false, true}){
			final String label=fill ? "filled" : "gapped";
			final ArrayList<String> options=new ArrayList<String>(Arrays.asList("in="+reads, "contigs="+contigs,
					"out="+root.resolve(label+".tsv"), "outjoins="+root.resolve(label+".joins.tsv"),
					"outassembly="+root.resolve(label+".fa"), "k=128", "tipbases=128", "join=t", "t=1"));
			if(fill){options.add("fillgaps=t");}
			launch(label, "dev/hollowbridge.sh", true, options.toArray(new String[0]));
			final String expected=fill ? reference.substring(0, 640)
					: reference.substring(0, 300)+"NNNNNNNNNN"+reference.substring(340, 640);
			check(sequence(root.resolve(label+".fa")).equals(expected), "Native Hollow CLI failed "+label+" spelling.");
		}
		for(String k : new String[]{"31", "31,64"}){
			final String label=k.equals("31") ? "single" : "multi";
			final ArrayList<String> options=new ArrayList<String>(Arrays.asList("in="+reads, "k="+k,
					"t=1", "bt=1", "prealloc=f", "prefilter=f", "branchrescue=f",
					"scaffold=t", "out="+root.resolve(label+".fa")));
			// Even fusenet=null selects multi-K dispatch, so it must be absent from the single-K fixture.
			if(!k.equals("31")){options.add("fusenet=null");}
			launch(label, "tadpole.sh", true, options.toArray(new String[0]));
			check(sequence(root.resolve(label+".fa")).length()>=500, "Native scaffolding must publish the completed assembly.");
			check(text(root.resolve(label+".log")).contains("Unscaffolded Tadpole contigs and link reports retained:"),
					"Native dispatch must report its preserved intermediate evidence.");
		}
		final String before=text(reads);
		launch("alias", "tadpole.sh", false, "in="+reads, "out="+reads, "scaffold=t", "ow=t", "t=1");
		check(text(reads).equals(before), "Rejected output aliases must not overwrite evidence.");
		launch("extension", "tadpole.sh", false, "in="+reads, "out="+root.resolve("invalid.fa"), "scaffold=t", "mode=extend", "t=1");
		check(!Files.exists(root.resolve("invalid.fa")), "Scaffolding must reject incompatible modes before writing output.");
	}

	/** Executes the shipped shell launcher in a fresh process and retains diagnostics. */
	private void launch(final String name, final String script, final boolean success, final String... options) throws Exception{
		assert(name!=null && script.endsWith(".sh")) : "Each test names a real launcher and a unique diagnostic log.";
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", code.resolve(script).toString(), "-Xmx256m", "-Xms32m"));
		command.addAll(Arrays.asList(options));
		final Path log=root.resolve(name+".log");
		final Process process=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		process.getOutputStream().close();
		check((process.waitFor()==0)==success, "Unexpected exit for "+name+"; inspect "+log);
	}

	/** Concatenates only these tiny fixture FASTA records; production readers remain native. */
	private static String sequence(final Path path) throws Exception{
		final StringBuilder bases=new StringBuilder();
		for(String line : Files.readAllLines(path, StandardCharsets.US_ASCII)){if(!line.startsWith(">")){bases.append(line);}}
		return bases.toString();
	}
	private static String text(final Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}
	private static void check(final boolean value, final String message){checks++; if(!value){throw new AssertionError(message);}}
	private final Path code, root, reads, contigs;
	private final String reference;
	private static int checks;
}
