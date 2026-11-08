package test.bloom;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import bloom.KCountArray2;
import bloom.KmerCount3;
import stream.FASTQ;

/**
 * Finite KmerCount3 cases over matching FASTQ and FASTA inputs.
 * Spawns one child JVM per case so a baseline null-quality failure cannot strand
 * the parent in ConcurrentReadInputStream.
 *
 * @author Brian Bushnell, Jean
 */
public final class KmerCount3FastaNullQualityTest{

	/** Runs bounded child JVMs against real FASTQ and FASTA paths. */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			child(args);
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path dir=Paths.get(args[0].substring(7)).toAbsolutePath();
		require(!Files.exists(dir), "Result directory must be new.");
		Files.createDirectories(dir);

		final List<CaseSpec> cases=Arrays.asList(
				new CaseSpec("fastq-single", false, "ACGTAC", null, null, new int[][]{
					{ACG, 1}, {CGT, 1}, {GTA, 1}, {TAC, 1}}),
				new CaseSpec("fasta-single", true, "ACGTAC", null, null, new int[][]{
					{ACG, 1}, {CGT, 1}, {GTA, 1}, {TAC, 1}}),
				new CaseSpec("fastq-paired", false, "ACGTAC", "TTTT", null, new int[][]{
					{ACG, 1}, {CGT, 1}, {GTA, 1}, {TAC, 1}, {TTT, 2}}),
				new CaseSpec("fasta-paired", true, "ACGTAC", "TTTT", null, new int[][]{
					{ACG, 1}, {CGT, 1}, {GTA, 1}, {TAC, 1}, {TTT, 2}}),
				new CaseSpec("fastq-low-quality", false, "ACGTAC", null, "???!!!", new int[][]{
					{ACG, 1}}));

		int failures=0;
		for(CaseSpec spec : cases){
			final Path first=dir.resolve(spec.name+"-1."+(spec.fasta ? "fa" : "fq"));
			final Path second=(spec.second==null ? null : dir.resolve(spec.name+"-2."+(spec.fasta ? "fa" : "fq")));
			if(spec.fasta){
				writeFasta(first, spec.first, spec.name+"/1");
				if(second!=null){writeFasta(second, spec.second, spec.name+"/2");}
			}else{
				writeFastq(first, spec.first, spec.firstQuality, spec.name+"/1");
				if(second!=null){writeFastq(second, spec.second, null, spec.name+"/2");}
			}
			final int code=runChild(dir, spec, first, second);
			if(code!=0){failures++;}
		}
		if(failures>0){throw new AssertionError("KmerCount3 FASTA/null-quality contract failed "+failures+" of "+cases.size()+" cases.");}
		System.out.println("KMER_COUNT3_FASTA_NULL_QUALITY_TEST_PASS cases="+cases.size());
	}

	/** Calls the legacy path-based counter and checks exact three-base cells. */
	private static void child(final String[] args){
		if(args.length!=4){throw new IllegalArgumentException("Usage: child <in1> <in2|-> <expectedCells>");}
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=false;
		FASTQ.ASCII_OFFSET=33;
		KmerCount3.minQuality=6;
		KmerCount3.maxReads=-1;
		KmerCount3.readsProcessed=0;

		final String in1=args[1];
		final String in2=("-".equals(args[2]) ? null : args[2]);
		final int[] expected=parseExpected(args[3]);
		final KCountArray2 count=KmerCount3.countFastq(in1, in2, K, 4);
		for(int key=0; key<expected.length; key++){
			final int observed=count.read(key);
			require(observed==expected[key],
					"key "+key+" expected "+expected[key]+" observed "+observed);
		}
		System.out.println("child pass in1="+in1+" in2="+in2);
	}

	/** Starts a bounded child and stores combined stdout/stderr beside its inputs. */
	private static int runChild(final Path dir, final CaseSpec spec, final Path first,
			final Path second) throws Exception{
		final List<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(KmerCount3FastaNullQualityTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx128m");
		if(!System.getProperty("java.specification.version").startsWith("1.")){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(KmerCount3FastaNullQualityTest.class.getName());
		command.add("child");
		command.add(first.toString());
		command.add(second==null ? "-" : second.toString());
		command.add(spec.expectedArgument());

		final Path log=dir.resolve(spec.name+".log");
		final Path commandFile=dir.resolve(spec.name+".command.txt");
		final Path outcome=dir.resolve(spec.name+".outcome.txt");
		Files.write(commandFile, command, StandardCharsets.UTF_8);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.redirectInput(Paths.get("/dev/null").toFile());
		pb.redirectErrorStream(true);
		pb.redirectOutput(log.toFile());
		final Process process=pb.start();
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
		if(!completed || !terminated){
			final String summary=spec.name+" completed="+completed+
					" terminated="+terminated+" exit=timeout log="+log;
			Files.write(outcome, Arrays.asList(summary), StandardCharsets.UTF_8);
			System.out.println(summary);
			if(!terminated){throw new AssertionError("Owned child survived forced destruction: "+spec.name);}
			return 1;
		}
		final int code=process.exitValue();
		final String summary=spec.name+" completed="+completed+
				" terminated="+terminated+" exit="+code+" log="+log;
		Files.write(outcome, Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return code;
	}

	/** Writes one ASCII-33 read. Null qualities mean every base gets Phred thirty. */
	private static void writeFastq(final Path path, final String bases,
			final String qualities, final String id) throws IOException{
		final char[] qualityChars=new char[bases.length()];
		if(qualities==null){
			Arrays.fill(qualityChars, '?');
		}else{
			require(qualities.length()==bases.length(), "Quality length must match bases.");
			qualities.getChars(0, qualities.length(), qualityChars, 0);
		}
		Files.write(path, ("@"+id+"\n"+bases+"\n+\n"+
				new String(qualityChars)+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Writes one FASTA read; production should receive a null quality array. */
	private static void writeFasta(final Path path, final String bases,
			final String id) throws IOException{
		Files.write(path, (">"+id+"\n"+bases+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Expands key=count pairs into a complete direct-address expectation. */
	private static int[] parseExpected(final String text){
		final int[] expected=new int[1<<(2*K)];
		if(text.length()<1){return expected;}
		for(String pair : text.split(",")){
			final int equals=pair.indexOf('=');
			require(equals>0, "Expected key=count entry: "+pair);
			final int key=Integer.parseInt(pair.substring(0, equals));
			final int value=Integer.parseInt(pair.substring(equals+1));
			expected[key]=value;
		}
		return expected;
	}

	/** Distinguishes invalid fixture setup from the baseline's targeted failure. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Literal input case and expected direct-address cells. */
	private static final class CaseSpec{
		CaseSpec(final String name_, final boolean fasta_, final String first_,
				final String second_, final String firstQuality_, final int[][] expected_){
			name=name_;
			fasta=fasta_;
			first=first_;
			second=second_;
			firstQuality=firstQuality_;
			expected=expected_;
		}

		/** Encodes only nonzero expected cells for the child JVM. */
		String expectedArgument(){
			final StringBuilder sb=new StringBuilder();
			for(int i=0; i<expected.length; i++){
				if(i>0){sb.append(',');}
				sb.append(expected[i][0]).append('=').append(expected[i][1]);
			}
			return sb.toString();
		}

		final String name, first, second, firstQuality;
		final boolean fasta;
		final int[][] expected;
	}

	private static final int K=3;
	private static final int ACG=6;
	private static final int CGT=27;
	private static final int GTA=44;
	private static final int TAC=49;
	private static final int TTT=63;
}
