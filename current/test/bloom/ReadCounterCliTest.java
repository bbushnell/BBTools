package test.bloom;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import bloom.ReadCounter;

/**
 * Checks ReadCounter CLI handling of the literal null mate placeholder.
 * Each case invokes bloom.ReadCounter in a fresh JVM against one exact-k FASTQ
 * record per file and records the raw command, input, stdout, stderr and outcome.
 *
 * @author Brian Bushnell, Jean
 */
public final class ReadCounterCliTest{

	/** Runs the null-mate regression and nearby real-mate control. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		require(!Files.exists(root), "Result directory must be new.");
		require(FIRST_BASES.length()==128 && SECOND_BASES.length()==128,
				"ReadCounter CLI fixture records must stay exactly one 128-base kmer.");
		Files.createDirectories(root);

		int failures=0;
		for(CaseSpec spec : CASES){failures+=runCase(root, spec);}
		if(failures>0){
			throw new AssertionError("ReadCounter CLI null-mate contract failed "+failures+
					" of "+CASES.length+" cases.");
		}
		System.out.println("READ_COUNTER_CLI_TEST_PASS cases="+CASES.length);
	}

	/** Writes one case's literal inputs, starts ReadCounter, and checks natural completion. */
	private static int runCase(final Path root, final CaseSpec spec) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(spec.name));
		final Path input1=dir.resolve("input1.fq");
		final Path input2=spec.paired ? dir.resolve("input2.fq") : null;
		final byte[] bytes1=fastq(spec.name+"/1", FIRST_BASES);
		final byte[] bytes2=spec.paired ? fastq(spec.name+"/2", SECOND_BASES) : null;
		Files.write(input1, bytes1);
		if(input2!=null){Files.write(input2, bytes2);}

		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(ReadCounterCliTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx128m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		if(hasVectorModule()){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ReadCounter.class.getName());
		command.add(input1.toString());
		command.add(spec.paired ? input2.toString() : "null");
		command.addAll(Arrays.asList("k=128", "cbits=8", "matrixbits=12", "hashes=1"));
		Files.write(dir.resolve("command.txt"), command, StandardCharsets.UTF_8);

		final Process process=new ProcessBuilder(command)
				.redirectInput(Paths.get("/dev/null").toFile())
				.redirectOutput(dir.resolve("stdout.txt").toFile())
				.redirectError(dir.resolve("stderr.txt").toFile()).start();
		boolean completed=false;
		boolean terminated=true;
		try{
			completed=process.waitFor(15, TimeUnit.SECONDS);
		}finally{
			if(process.isAlive()){
				process.destroy();
				if(!process.waitFor(2, TimeUnit.SECONDS)){
					process.destroyForcibly();
					terminated=process.waitFor(5, TimeUnit.SECONDS);
				}
			}
		}

		final int exit=terminated && completed ? process.exitValue() : Integer.MIN_VALUE;
		final String stdout=new String(Files.readAllBytes(dir.resolve("stdout.txt")), StandardCharsets.UTF_8);
		final boolean inputsPreserved=Arrays.equals(bytes1, Files.readAllBytes(input1)) &&
				(input2==null || Arrays.equals(bytes2, Files.readAllBytes(input2)));
		final boolean pass=completed && terminated && exit==0 && inputsPreserved &&
				hasLine(stdout, "Keys Counted:  \t         \t"+spec.expectedKeys);
		final List<String> outcome=Arrays.asList(
				"case="+spec.name,
				"completed="+completed,
				"terminated="+terminated,
				"exit="+exit,
				"inputsPreserved="+inputsPreserved,
				"expectedKeys="+spec.expectedKeys,
				"pass="+pass);
		Files.write(dir.resolve("outcome.txt"), outcome, StandardCharsets.UTF_8);
		System.out.println(String.join(" ", outcome));
		if(!terminated){throw new AssertionError("Owned child survived forced destruction: "+spec.name);}
		return pass ? 0 : 1;
	}

	/** Creates one exact-128-base FASTQ record with fixed Phred-40 qualities. */
	private static byte[] fastq(final String name, final String bases){
		final char[] quality=new char[bases.length()];
		Arrays.fill(quality, 'I');
		return ("@"+name+"\n"+bases+"\n+\n"+new String(quality)+"\n")
				.getBytes(StandardCharsets.US_ASCII);
	}

	/** Matches complete output lines so a longer incorrect count cannot pass by prefix. */
	private static boolean hasLine(final String report, final String expected){
		for(String line : report.split("\\r?\\n")){
			if(trimRight(line).equals(trimRight(expected))){return true;}
		}
		return false;
	}

	/** Removes only the formatter's trailing padding. */
	private static String trimRight(final String line){
		int end=line.length();
		while(end>0 && Character.isWhitespace(line.charAt(end-1))){end--;}
		return line.substring(0, end);
	}

	/** Distinguishes invalid fixture setup from the targeted CLI failure. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Returns true only when this JVM can resolve the incubator Vector module. */
	private static boolean hasVectorModule(){
		try{
			Class.forName("jdk.incubator.vector.ByteVector", false, ReadCounterCliTest.class.getClassLoader());
			return true;
		}catch(ClassNotFoundException e){
			return false;
		}catch(LinkageError e){
			return false;
		}
	}

	/** Literal CLI case. */
	private static final class CaseSpec{
		CaseSpec(final String name_, final boolean paired_, final int expectedKeys_){
			name=name_;
			paired=paired_;
			expectedKeys=expectedKeys_;
		}

		final String name;
		final boolean paired;
		final int expectedKeys;
	}

	private static final CaseSpec[] CASES={
			new CaseSpec("single-null-placeholder", false, 1),
			new CaseSpec("paired-real-mate", true, 2)
	};

	private static final String FIRST_BASES="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"+
			"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	private static final String SECOND_BASES="CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"+
			"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
}
