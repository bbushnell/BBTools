package test.tax;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

import tax.ShrinkAccession;

/** Tests ShrinkAccession assembly_summary header/data alignment. */
public final class ShrinkAccessionAssemblyHeaderTest{

	/** Runs one bounded assembly-mode header/data check.
	 * @param args Optional expectation flag, expect=baseline or expect=fixed
	 * @throws Exception If a temporary fixture cannot be written or read */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			runChild(parseExpectation(new String[] {args[1]}), parseRoot(args));
			return;
		}
		runParent(parseExpectation(args));
	}

	/** Prevents construction of this standalone fixture. */
	private ShrinkAccessionAssemblyHeaderTest(){}

	/** Maximum time to wait for one ShrinkAccession child. */
	private static final long CHILD_TIMEOUT_MS=15000;
	/** Maximum time to wait for a destroyed child. */
	private static final long DESTROY_TIMEOUT_MS=1000;

	/** Assembly-summary input with one header, one skipped body comment and two rows. */
	private static final String INPUT=headerLine()+
			assemblyLine("GCF_000001215.4", "7227", "haploid", "143706478",
					"142553500", "7", "1869", "2441", "17872", "13962", "3543")+
			// Retain a body comment to exercise the existing mid-stream skip control.
			"#assembly_accession\tbody_comment_should_be_skipped\n"+
			assemblyLine("GCF_000005845.2", "562", "haploid", "4641652",
					"4641652", "1", "1", "1", "4494", "4140", "143");

	/** Current output: the header omits the blank compatibility column. */
	private static final String BASELINE_OUTPUT=
			"#assembly_accession\ttaxid\tassembly_type\tgenome_size\tgenome_size_ungapped\t"+
			"replicon_count\tscaffold_count\tcontig_count\ttotal_gene_count\t"+
			"protein_coding_gene_count\tnon_coding_gene_count\n"+
			"GCF_000001215.4\t\t7227\thaploid\t143706478\t142553500\t7\t"+
			"1869\t2441\t17872\t13962\t3543\n"+
			"GCF_000005845.2\t\t562\thaploid\t4641652\t4641652\t1\t1\t1\t"+
			"4494\t4140\t143\n";

	/** Fixed output: the header keeps the same blank compatibility column as data rows. */
	private static final String FIXED_OUTPUT=
			"#assembly_accession\t\ttaxid\tassembly_type\tgenome_size\t"+
			"genome_size_ungapped\treplicon_count\tscaffold_count\tcontig_count\t"+
			"total_gene_count\tprotein_coding_gene_count\tnon_coding_gene_count\n"+
			"GCF_000001215.4\t\t7227\thaploid\t143706478\t142553500\t7\t"+
			"1869\t2441\t17872\t13962\t3543\n"+
			"GCF_000005845.2\t\t562\thaploid\t4641652\t4641652\t1\t1\t1\t"+
			"4494\t4140\t143\n";

	/** Runs the actual ShrinkAccession case in one bounded owned child JVM.
	 * @param expectation Expected behavior: baseline or fixed
	 * @throws Exception If the child fails, hangs or cannot be cleaned up */
	private static void runParent(final String expectation) throws Exception{
		final Path root=Files.createTempDirectory("shrink_assembly_header_");
		final Path childRoot=root.resolve("child");
		final Path stdout=root.resolve("child.stdout");
		final Path stderr=root.resolve("child.stderr");
		Files.createDirectories(childRoot);
		Process child=null;
		final long started=System.nanoTime();
		try{
			child=new ProcessBuilder(childCommand(expectation, childRoot))
					.redirectOutput(stdout.toFile())
					.redirectError(stderr.toFile())
					.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				child.destroy();
				if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
					child.destroyForcibly().waitFor();
				}
				throw new AssertionError("ShrinkAccession child timed out\nstdout\n"+
						readAscii(stdout)+"\nstderr\n"+readAscii(stderr));
			}
			final String childOut=readAscii(stdout);
			final String childErr=readAscii(stderr);
			System.out.println("shrink-assembly-header child exit="+child.exitValue()+
					" elapsed_ms="+
					TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started)+
					" stdout='"+compact(childOut)+"' stderr='"+compact(childErr)+"'");
			check(child.exitValue()==0, "child exit "+child.exitValue()+"\nstdout\n"+
					childOut+"\nstderr\n"+childErr+"\nactual\n"+
					readAscii(childRoot.resolve("assembly.out")));
		}finally{
			if(child!=null && child.isAlive()){
				child.destroy();
				if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
					child.destroyForcibly().waitFor();
				}
			}
			deleteTree(root);
		}
		System.out.println("shrink-assembly-header expect="+expectation);
	}

	/** Builds a child JVM command with the parent assertion mode.
	 * @param expectation Expected behavior: baseline or fixed
	 * @param root Child artifact directory
	 * @return Child command */
	private static List<String> childCommand(final String expectation, final Path root){
		final List<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+
				File.separator+"java");
		command.add(assertionsEnabled() ? "-ea" : "-da");
		command.add("-Xmx512m");
		command.add("-Xms512m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ShrinkAccessionAssemblyHeaderTest.class.getName());
		command.add("child");
		command.add("expect="+expectation);
		command.add("root="+root.toString());
		return command;
	}

	/** Runs the real ShrinkAccession invocation inside a bounded child JVM.
	 * @param expectation Expected behavior: baseline or fixed
	 * @param root Child artifact directory
	 * @throws Exception If the output diverges from the oracle */
	private static void runChild(final String expectation, final Path root) throws Exception{
		Files.createDirectories(root);
		final Path in=root.resolve("assembly.tsv");
		final Path out=root.resolve("assembly.out");
		Files.write(in, INPUT.getBytes(StandardCharsets.US_ASCII));
		ShrinkAccession.main(new String[] {"in="+in.toString(),
				"out="+out.toString(), "asm=t", "ow=t"});
		final String actual=readAscii(out);
		final String expected="baseline".equals(expectation) ? BASELINE_OUTPUT : FIXED_OUTPUT;
		check(expected.equals(actual), "unexpected assembly output\nexpected:\n"+
				expected+"\nactual:\n"+actual);
		checkDataValues(actual, expectation);
		System.out.println("case=asm expect="+expectation+" header="+
				assemblyFields(actual, 0).length+" data="+assemblyFields(actual, 1).length);
	}

	/** Parses the BBTools-style expectation flag.
	 * @param args Command-line arguments
	 * @return baseline or fixed */
	private static String parseExpectation(final String[] args){
		if(args.length<1){return "fixed";}
		check(args.length==1, "unknown extra argument "+(args.length>1 ? args[1] : ""));
		final String arg=args[0];
		final String prefix="expect=";
		check(arg.regionMatches(true, 0, prefix, 0, prefix.length()), "unknown argument "+arg);
		final String expectation=arg.substring(prefix.length());
		check("baseline".equals(expectation) || "fixed".equals(expectation),
				"expectation must be baseline or fixed, not "+expectation);
		return expectation;
	}

	/** Parses the child artifact directory.
	 * @param args Child command-line arguments
	 * @return Child artifact directory */
	private static Path parseRoot(final String[] args){
		check(args.length==3, "child requires expect=... and root=...");
		final String prefix="root=";
		check(args[2].startsWith(prefix), "unknown child argument "+args[2]);
		return Paths.get(args[2].substring(prefix.length()));
	}

	/** Builds a literal assembly_summary header line.
	 * @return Header line */
	private static String headerLine(){
		final String[] fields=filledLine();
		fields[0]="#assembly_accession";
		fields[5]="taxid";
		fields[23]="assembly_type";
		fields[25]="genome_size";
		fields[26]="genome_size_ungapped";
		fields[28]="replicon_count";
		fields[29]="scaffold_count";
		fields[30]="contig_count";
		fields[34]="total_gene_count";
		fields[35]="protein_coding_gene_count";
		fields[36]="non_coding_gene_count";
		return join(fields);
	}

	/** Builds one literal assembly_summary data line.
	 * @param accession Assembly accession
	 * @param taxid NCBI TaxID
	 * @param assemblyType Assembly type
	 * @param genomeSize Genome size
	 * @param ungappedSize Ungapped genome size
	 * @param replicons Replicon count
	 * @param scaffolds Scaffold count
	 * @param contigs Contig count
	 * @param genes Total gene count
	 * @param coding Protein-coding gene count
	 * @param noncoding Non-coding gene count
	 * @return Data line */
	private static String assemblyLine(final String accession, final String taxid,
			final String assemblyType, final String genomeSize, final String ungappedSize,
			final String replicons, final String scaffolds, final String contigs,
			final String genes, final String coding, final String noncoding){
		final String[] fields=filledLine();
		fields[0]=accession;
		fields[5]=taxid;
		fields[23]=assemblyType;
		fields[25]=genomeSize;
		fields[26]=ungappedSize;
		fields[28]=replicons;
		fields[29]=scaffolds;
		fields[30]=contigs;
		fields[34]=genes;
		fields[35]=coding;
		fields[36]=noncoding;
		return join(fields);
	}

	/** Builds the ignored fields in a literal 37-column assembly_summary row.
	 * @return Filled field array */
	private static String[] filledLine(){
		final String[] fields=new String[37];
		Arrays.fill(fields, "na");
		return fields;
	}

	/** Joins one assembly_summary row.
	 * @param fields Fields to join
	 * @return Tab-delimited row with trailing newline */
	private static String join(final String[] fields){
		final StringBuilder sb=new StringBuilder(256);
		for(int i=0; i<fields.length; i++){
			if(i>0){sb.append('\t');}
			sb.append(fields[i]);
		}
		return sb.append('\n').toString();
	}

	/** Checks that data fields are unchanged by the header repair.
	 * @param actual Actual ShrinkAccession output
	 * @param expectation Expected behavior: baseline or fixed */
	private static void checkDataValues(final String actual, final String expectation){
		final String[] lines=actual.split("\n", -1);
		check(lines.length==4 && lines[3].length()==0,
				"expected header, two data lines and one final newline");
		final String[] header=assemblyFields(actual, 0);
		check(header.length==("baseline".equals(expectation) ? 11 : 12),
				"header field count did not match "+expectation);
		check(assemblyFields(actual, 1).length==12, "first data row lost a field");
		check(assemblyFields(actual, 2).length==12, "second data row lost a field");
		checkDataLine(lines[1], "GCF_000001215.4", "7227", "haploid", "143706478",
				"142553500", "7", "1869", "2441", "17872", "13962", "3543");
		checkDataLine(lines[2], "GCF_000005845.2", "562", "haploid", "4641652",
				"4641652", "1", "1", "1", "4494", "4140", "143");
	}

	/** Splits an output row into fields.
	 * @param text Output text
	 * @param lineNumber Zero-based line index
	 * @return Fields in the requested row */
	private static String[] assemblyFields(final String text, final int lineNumber){
		return text.split("\n", -1)[lineNumber].split("\t", -1);
	}

	/** Checks all retained assembly output values.
	 * @param line Output line
	 * @param accession Expected accession
	 * @param taxid Expected taxid
	 * @param assemblyType Expected assembly type
	 * @param genomeSize Expected genome size
	 * @param ungappedSize Expected ungapped genome size
	 * @param replicons Expected replicon count
	 * @param scaffolds Expected scaffold count
	 * @param contigs Expected contig count
	 * @param genes Expected gene count
	 * @param coding Expected protein-coding gene count
	 * @param noncoding Expected non-coding gene count */
	private static void checkDataLine(final String line, final String accession,
			final String taxid, final String assemblyType, final String genomeSize,
			final String ungappedSize, final String replicons, final String scaffolds,
			final String contigs, final String genes, final String coding,
			final String noncoding){
		final String[] fields=line.split("\t", -1);
		check(fields[0].equals(accession), "accession changed");
		check(fields[1].length()==0, "blank compatibility column changed");
		check(fields[2].equals(taxid), "taxid changed");
		check(fields[3].equals(assemblyType), "assembly type changed");
		check(fields[4].equals(genomeSize), "genome size changed");
		check(fields[5].equals(ungappedSize), "ungapped genome size changed");
		check(fields[6].equals(replicons), "replicon count changed");
		check(fields[7].equals(scaffolds), "scaffold count changed");
		check(fields[8].equals(contigs), "contig count changed");
		check(fields[9].equals(genes), "gene count changed");
		check(fields[10].equals(coding), "coding gene count changed");
		check(fields[11].equals(noncoding), "noncoding gene count changed");
	}

	/** Checks whether Java assertions are enabled.
	 * @return True if assertions are enabled */
	private static boolean assertionsEnabled(){
		boolean enabled=false;
		assert enabled=true;
		return enabled;
	}

	/** Adds the incubator vector module on JDKs where it is present.
	 * @param command Command under construction */
	private static void addVectorModule(final List<String> command){
		try{
			Class.forName("jdk.incubator.vector.ByteVector");
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}catch(final ClassNotFoundException e){
			// Java 8: no module flag needed or accepted.
		}
	}

	/** Reads one US-ASCII file, returning empty text when absent.
	 * @param path File to read
	 * @return Decoded file contents
	 * @throws Exception If the file exists but cannot be read */
	private static String readAscii(final Path path) throws Exception{
		return Files.exists(path) ?
				new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}

	/** Collapses process output for a parent-visible one-line child report.
	 * @param text Text to compact
	 * @return Single-line compact text */
	private static String compact(final String text){
		return text.trim().replaceAll("\\s+", " ");
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		try(final Stream<Path> paths=Files.walk(root)){
			for(final Path path : (Iterable<Path>)paths.sorted(
					Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under both assertion modes.
	 * @param b Required condition
	 * @param message Failure explanation */
	private static void check(final boolean b, final String message){
		if(!b){throw new AssertionError(message);}
	}
}
