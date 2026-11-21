package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parse;

/** Finite CLI oracles for per-input, nonempty FASTA-record limits in Stats3.
 * @author Shinobu */
public final class AssemblyStats3LimitTest{

	/** Runs fresh CLI children against literal source and corrected expectations. */
	public static void main(String[] args) throws Exception{
		String phase=null, directory=null, root=null;
		for(String arg : args){
			String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){root=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled || !("baseline".equals(phase) || "fixed".equals(phase)) || directory==null || root==null){
			throw new IllegalArgumentException("Use phase=baseline|fixed outdir=<new directory> bbtools=<root> with -ea");
		}
		new AssemblyStats3LimitTest(Paths.get(directory), Paths.get(root), phase).run();
	}

	private AssemblyStats3LimitTest(Path out_, Path root_, String phase_){
		out=out_.toAbsolutePath().normalize();
		root=root_.toAbsolutePath().normalize();
		phase=phase_;
		baseline=phase.equals("baseline");
	}

	private void run() throws Exception{
		Files.createDirectory(out);
		Path inputs=out.resolve("inputs");
		Files.createDirectory(inputs);
		String[] names={"first.fa", "second.fa", "empty.fa", "empty-records.fa", "prefix.fa"};
		String[] text={">a\n"+repeat('A', 5000)+"\n>b\n"+repeat('C', 10000)+"\n>c\nGGGGGGGGGG\n",
			">d\nGGG\n>e\nTTTTT\n", "", ">empty\n>first\nAC\n>empty2\n>second\nGGG\n", "AC\n>later\nGGG\n"};
		for(int i=0; i<names.length; i++){
			files[i]=inputs.resolve(names[i]);
			bytes[i]=text[i].getBytes(StandardCharsets.US_ASCII);
			Files.write(files[i], bytes[i], StandardOpenOption.CREATE_NEW);
		}
		Files.write(out.resolve("cases.tsv"), "case\tphase\texit\tsourceMatches\tcorrectedMatches\tstdoutEmpty\tstderrClean\tinputsUnchanged\n".getBytes(StandardCharsets.UTF_8));
		for(Case c : CASES){runCase(c);}
		check(children==10, "Expected ten fresh children, found "+children);
		check(sourceMatches==(baseline ? 10 : 4), "Unexpected source matches: "+sourceMatches);
		check(correctedMatches==(baseline ? 4 : 10), "Unexpected corrected matches: "+correctedMatches);
		String summary="ASSEMBLY_STATS3_LIMIT_"+(failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" cases=10 children="+children+" sourceMatches="+sourceMatches+
			" correctedMatches="+correctedMatches+" checks="+checks+" failures="+failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8));
		if(failures>0){throw new AssertionError(summary+"See "+out);}
		System.out.print(summary);
	}

	private void runCase(Case c) throws Exception{
		Path dir=out.resolve(c.name);
		Files.createDirectory(dir);
		Path result=dir.resolve("stats.tsv"), stdout=dir.resolve("stdout.txt"), stderr=dir.resolve("stderr.txt");
		ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(root.resolve("stats3.sh").toString());
		command.add("-Xmx128m"); command.add("-ea"); command.add("t=1"); command.add("simd=f");
		command.add("ow=t"); command.add("app=f"); command.add("out="+result);
		StringBuilder in=new StringBuilder("in=");
		for(int i=0; i<c.inputs.length; i++){
			if(i>0){in.append(',');}
			in.append(files[c.inputs[i]]);
		}
		command.add(in.toString());
		if(c.limit!=null){command.add(c.limit);}
		StringBuilder argv=new StringBuilder();
		for(String arg : command){argv.append(arg).append('\n');}
		Files.write(dir.resolve("argv.txt"), argv.toString().getBytes(StandardCharsets.UTF_8));
		byte[] source=expected(c, c.source), corrected=expected(c, c.corrected);
		Files.write(dir.resolve("source.tsv"), source);
		Files.write(dir.resolve("corrected.tsv"), corrected);
		ProcessBuilder pb=new ProcessBuilder(command).directory(root.toFile()).redirectOutput(stdout.toFile()).redirectError(stderr.toFile());
		pb.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		int exit=pb.start().waitFor();
		children++;
		Files.write(dir.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII));
		byte[] actual=Files.isRegularFile(result) ? Files.readAllBytes(result) : new byte[0];
		boolean sourceMatch=Arrays.equals(source, actual), correctedMatch=Arrays.equals(corrected, actual);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		boolean quiet=Files.size(stdout)==0;
		String err=new String(Files.readAllBytes(stderr), StandardCharsets.UTF_8);
		boolean clean=err.contains("Picked up JAVA_TOOL_OPTIONS") && !err.contains("Exception") && !err.contains("Error") && !err.contains("Unknown parameter");
		boolean unchanged=true;
		for(int i=0; i<files.length; i++){unchanged&=Arrays.equals(bytes[i], Files.readAllBytes(files[i]));}
		Files.write(out.resolve("cases.tsv"), (c.name+"\t"+phase+"\t"+exit+"\t"+sourceMatch+"\t"+correctedMatch+"\t"+quiet+"\t"+clean+"\t"+unchanged+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(exit==0, c.name+": child exit "+exit);
		check(baseline ? sourceMatch : correctedMatch, c.name+": wrong statistics table");
		check(quiet, c.name+": unexpected stdout");
		check(clean, c.name+": unexpected diagnostic: "+err);
		check(unchanged, c.name+": input changed");
	}

	/** Combines literal numeric rows with the actual input names, without deriving statistics. */
	private byte[] expected(Case c, String[] rows){
		StringBuilder result=new StringBuilder(HEADER);
		for(int i=0; i<rows.length; i++){result.append(files[c.inputs[i]]).append('\t').append(rows[i]).append('\n');}
		return result.toString().getBytes(StandardCharsets.UTF_8);
	}

	private void check(boolean ok, String message){checks++; if(!ok){failures++; System.err.println(message);}}
	private static String repeat(char base, int length){char[] a=new char[length]; Arrays.fill(a, base); return new String(a);}

	private static final class Case{
		Case(String name_, int[] inputs_, String limit_, String[] source_, String[] corrected_){
			name=name_; inputs=inputs_; limit=limit_; source=source_; corrected=corrected_;
		}
		final String name, limit;
		final int[] inputs;
		final String[] source, corrected;
	}

	private final Path out, root;
	private final String phase;
	private final boolean baseline;
	private final Path[] files=new Path[5];
	private final byte[][] bytes=new byte[5][];
	private int children, sourceMatches, correctedMatches, checks, failures;
	private static final String HEADER="fname\tsize\tcontigs\tgc\tmaxContig\t5kplus\t10kplus\t25kplus\t50kplus\n";
	private static final String FULL="15010\t3\t0.667\t10000\t15000\t10000\t0\t0";
	private static final String FIRST="5000\t1\t0\t5000\t5000\t0\t0\t0";
	private static final String TWO="15000\t2\t0.667\t10000\t15000\t10000\t0\t0";
	private static final String ZERO="0\t0\t0\t0\t0\t0\t0\t0";
	private static final String SECOND="8\t2\t0.375\t5\t0\t0\t0\t0";
	private static final String SECOND_FIRST="3\t1\t1\t3\t0\t0\t0\t0";
	private static final String SHORT_FULL="5\t2\t0.800\t3\t0\t0\t0\t0";
	private static final String SHORT_FIRST="2\t1\t0.500\t2\t0\t0\t0\t0";
	private static final Case[] CASES={
		new Case("first-record", new int[]{0}, "reads=1", new String[]{FULL}, new String[]{FIRST}),
		new Case("two-record-alias", new int[]{0}, "maxreads=2", new String[]{FULL}, new String[]{TWO}),
		new Case("zero", new int[]{0}, "reads=0", new String[]{FULL}, new String[]{ZERO}),
		new Case("unlimited", new int[]{0}, "reads=-1", new String[]{FULL}, new String[]{FULL}),
		new Case("beyond-eof", new int[]{0}, "reads=99", new String[]{FULL}, new String[]{FULL}),
		new Case("default", new int[]{0}, null, new String[]{FULL}, new String[]{FULL}),
		new Case("per-input-reset", new int[]{0, 1}, "reads=1", new String[]{FULL, SECOND}, new String[]{FIRST, SECOND_FIRST}),
		new Case("empty-input", new int[]{2}, "reads=1", new String[]{ZERO}, new String[]{ZERO}),
		new Case("empty-records", new int[]{3}, "reads=1", new String[]{SHORT_FULL}, new String[]{SHORT_FIRST}),
		new Case("preheader-sequence", new int[]{4}, "reads=1", new String[]{SHORT_FULL}, new String[]{SHORT_FIRST})
	};
}
