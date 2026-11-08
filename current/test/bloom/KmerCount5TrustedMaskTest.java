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

import bloom.ErrorCorrect;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import bloom.KmerCount5;
import shared.Shared;
import stream.FASTQ;

/**
 * Finite KmerCount5 trusted-mask cases over real file inputs.
 * Spawns one child JVM per case so a baseline short-read or null-quality failure
 * cannot strand the parent in ConcurrentReadInputStream.
 *
 * @author Brian Bushnell, Jean
 */
public final class KmerCount5TrustedMaskTest{

	/** Runs bounded child JVMs against real FASTQ and FASTA paths. */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			child(args);
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){throw new IllegalArgumentException("Usage: outdir=<new directory>");}
		final Path dir=Paths.get(args[0].substring(7)).toAbsolutePath();
		require(!Files.exists(dir), "Result directory must be new.");
		Files.createDirectories(dir);

		final List<CaseSpec> cases=Arrays.asList(
				new CaseSpec("fq-short-single", false, "CC", null, 0),
				new CaseSpec("fq-long-single", false, "CCCCCA", null, 1),
				new CaseSpec("fq-short-long", false, "CC", "CCCCCA", 1),
				new CaseSpec("fq-long-short", false, "CCCCCA", "CC", 1),
				new CaseSpec("fq-short-short", false, "CC", "CC", 0),
				new CaseSpec("fq-long-long", false, "CCCCCA", "CCCCCA", 2),
				new CaseSpec("fa-long-null-quality", true, "CCCCCA", null, 1));

		int failures=0;
		for(CaseSpec spec : cases){
			final Path first=dir.resolve(spec.name+"-1."+(spec.fasta ? "fa" : "fq"));
			final Path second=(spec.second==null ? null : dir.resolve(spec.name+"-2."+(spec.fasta ? "fa" : "fq")));
			if(spec.fasta){
				writeFasta(first, spec.first, spec.name+"/1");
				if(second!=null){writeFasta(second, spec.second, spec.name+"/2");}
			}else{
				writeFastq(first, spec.first, spec.name+"/1");
				if(second!=null){writeFastq(second, spec.second, spec.name+"/2");}
			}
			final int code=runChild(dir, spec, first, second);
			if(code!=0){failures++;}
		}
		if(failures>0){throw new AssertionError("KmerCount5 trusted-mask contract failed "+failures+" of "+cases.size()+" cases.");}
		System.out.println("KMER_COUNT5_TRUSTED_MASK_TEST_PASS cases="+cases.size());
	}

	/** Calls KmerCount5 exactly through its path-based trusted-mask overload. */
	private static void child(final String[] args){
		if(args.length!=4){throw new IllegalArgumentException("Usage: child <in1> <in2|-> <expectedLongReads>");}
		Shared.setThreads(1);
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=false;
		FASTQ.ASCII_OFFSET=33;
		ErrorCorrect.MIN_ADVANCE=1;
		ErrorCorrect.MAX_ERROR_BURST=3;
		ErrorCorrect.OUTPUT_INFO=false;
		ErrorCorrect.verbose=false;

		final String in1=args[1];
		final String in2=("-".equals(args[2]) ? null : args[2]);
		final int expectedLongReads=Integer.parseInt(args[3]);
		KCountArray7MTA.setSeed(20261105);
		final KCountArray trusted=KCountArray.makeNew(2000003, 4, 3);
		final KCountArray output=KCountArray.makeNew(2000003, 4, 3);
		trusted.increment(CCCCC, 5);
		trusted.increment(GGGGG, 5);
		for(int key=0; key<1<<(2*K); key++){
			require(trusted.read(key)==(key==CCCCC || key==GGGGG ? 5 : 0) && output.read(key)==0,
					"Literal five-base tables must be collision-free before the path case.");
		}

		KmerCount5.count(in1, in2, K, 4, false, output, trusted, -1, 2, 1, true);
		require(output.read(CCCCC)==expectedLongReads,
				"Expected "+expectedLongReads+" trusted CCCCC windows, observed "+output.read(CCCCC));
		require(output.read(CCCCA)==0,
				"Unsupported terminal base should be masked before CCCCA is counted.");
		System.out.println("child pass in1="+in1+" in2="+in2+" ccccc="+output.read(CCCCC));
	}

	/** Starts a bounded child and stores combined stdout/stderr beside its inputs. */
	private static int runChild(final Path dir, final CaseSpec spec, final Path first, final Path second) throws Exception{
		final List<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(KmerCount5TrustedMaskTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx128m");
		if(!System.getProperty("java.specification.version").startsWith("1.")){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(KmerCount5TrustedMaskTest.class.getName());
		command.add("child");
		command.add(first.toString());
		command.add(second==null ? "-" : second.toString());
		command.add(Integer.toString(spec.expectedLongReads));

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
			final String summary=spec.name+" completed="+completed+" terminated="+terminated+" exit=timeout log="+log;
			Files.write(outcome, Arrays.asList(summary), StandardCharsets.UTF_8);
			System.out.println(summary);
			if(!terminated){throw new AssertionError("Owned child survived forced destruction: "+spec.name);}
			return 1;
		}
		final int code=process.exitValue();
		final String summary=spec.name+" completed="+completed+" terminated="+terminated+" exit="+code+" log="+log;
		Files.write(outcome, Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return code;
	}

	/** Writes one unambiguous ASCII-33 read with quality thirty. */
	private static void writeFastq(final Path path, final String bases, final String id) throws IOException{
		final char[] qualities=new char[bases.length()];
		Arrays.fill(qualities, '?');
		Files.write(path, ("@"+id+"\n"+bases+"\n+\n"+new String(qualities)+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Writes one FASTA read; production should receive a null quality array. */
	private static void writeFasta(final Path path, final String bases, final String id) throws IOException{
		Files.write(path, (">"+id+"\n"+bases+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Distinguishes invalid fixture setup from the baseline's targeted process failure. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Literal input case and expected number of long reads that survive to key CCCCC. */
	private static final class CaseSpec{
		CaseSpec(final String name_, final boolean fasta_, final String first_, final String second_, final int expectedLongReads_){
			name=name_;
			fasta=fasta_;
			first=first_;
			second=second_;
			expectedLongReads=expectedLongReads_;
		}
		final String name, first, second;
		final boolean fasta;
		final int expectedLongReads;
	}

	private static final int K=5;
	private static final int CCCCC=341;
	private static final int CCCCA=340;
	private static final int GGGGG=682;
}
