package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import jgi.KmerCountExact;
import shared.Timer;

/**
 * Checks KmerCountExact propagation of late k-mer dump output failures.
 * Successful controls write the same tiny valid input to normal files.  The
 * /dev/full cases are Linux-only and must report a process-level error through
 * the serial and two-output-thread dump paths.  Run each case in a fresh JVM.
 * @author Jean
 */
public final class KmerCountExactOutputTest{

	/** Accepts case=serial-control|parallel-control|serial-full|parallel-full and a new dir. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(dir==null || !("serial-control".equals(mode) || "parallel-control".equals(mode) ||
				"serial-full".equals(mode) || "parallel-full".equals(mode))){
			throw new IllegalArgumentException("Expected a serial/parallel control/full case and dir=<new directory>");
		}

		final boolean parallel=mode.startsWith("parallel"), full=mode.endsWith("full");
		if(full){requireDevFull();}
		Files.createDirectories(dir);
		final Path input=dir.resolve("input.fa"), output=full ? Paths.get("/dev/full") : dir.resolve("kmers.fa"),
			hist=parallel ? dir.resolve("hist.txt") : null;
		Files.write(input, ">input\nACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(!full){check(!Files.exists(output), "Output must be new so a previous dump cannot satisfy the control");}
		if(hist!=null){check(!Files.exists(hist), "Histogram must be new so a previous run cannot satisfy the control");}

		final ArrayList<String> options=new ArrayList<String>();
		options.add("in="+input);
		options.add("out="+output);
		options.add("k=3");
		options.add("ways=3");
		options.add("threads="+(parallel ? 2 : 1));
		options.add("prefilter=0");
		options.add("prealloc=f");
		options.add("onepass=f");
		options.add("rcomp=t");
		options.add("minprob=0");
		options.add("interleaved=f");
		options.add("minavgquality=0");
		options.add("qtrim=f");
		options.add("ecco=f");
		options.add("merge=f");
		options.add("overwrite=t");
		if(hist!=null){options.add("hist="+hist);}
		RuntimeException failure=null;
		try{
			new KmerCountExact(options.toArray(new String[0])).process(new Timer());
		}catch(RuntimeException re){
			failure=re;
		}
		if(full){
			check(failure!=null && failure.getMessage()!=null &&
				failure.getMessage().contains("KmerCountExact terminated in an error state"),
				"/dev/full must report KmerCountExact's terminal error state, not "+failure);
		}else{
			check(failure==null, "Successful control must not fail: "+failure);
		}
		if(!full){check(Files.size(output)>0, "Successful control must emit dumped k-mers");}
		if(hist!=null){check(Files.size(hist)>0, "Parallel control must emit the histogram side output");}
		System.out.println("PASS case="+mode+" full="+full+" parallel="+parallel+" failure="+(failure!=null));
	}

	/** Restricts /dev/full oracles to Linux hosts with the expected special file. */
	private static void requireDevFull(){
		check(System.getProperty("os.name").toLowerCase().contains("linux"),
			"/dev/full failure cases are Linux-only");
		check(Files.exists(Paths.get("/dev/full")) && !Files.isRegularFile(Paths.get("/dev/full")),
			"/dev/full must exist as a special write-failing device");
	}

	/** Keeps error-propagation oracles active with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
