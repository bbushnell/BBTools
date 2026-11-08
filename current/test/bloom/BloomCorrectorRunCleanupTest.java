package test.bloom;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;

import bloom.BloomFilterCorrectorWrapper;

/** Linux integration checks for wrapper failure propagation and diagnostic closure.
 * Uses /dev/full for a real output failure and /proc/self/fd for open-file evidence.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorRunCleanupTest{

	/** Retains tiny inputs, outputs and diagnostics beneath the supplied fresh directory. */
	public static void main(final String[] args) throws IOException{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new result directory>");
		}
		check(Files.isDirectory(Paths.get("/proc/self/fd")) && Files.exists(Paths.get("/dev/full")),
				"This Linux-only integration test requires /proc/self/fd and /dev/full.");
		final Path dir=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(dir), "The result directory must be new so prior artifacts cannot satisfy checks.");
		Files.createDirectories(dir);
		final Path input=dir.resolve("input.fq");
		Files.write(input, ("@fixture\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\n"+
				"55555555555555555555555555555555\n").getBytes(StandardCharsets.US_ASCII));
		//Prove that the descriptor oracle detects an open file, then its closure.
		final Path probe=dir.resolve("descriptor-probe");
		try(FileOutputStream held=new FileOutputStream(probe.toFile())){
			held.flush();
			check(openDescriptors(probe)>0, "The descriptor oracle must detect the deliberately open control file.");
		}
		check(openDescriptors(probe)==0, "The descriptor control must disappear after close.");
		exercise(dir, input, false);
		exercise(dir, input, true);
		System.out.println("BLOOM_CORRECTOR_RUN_CLEANUP_TEST_PASS cases=2");
	}

	/** Checks successful processing and a failure after construction, in the same JVM. */
	private static void exercise(final Path dir, final Path input, final boolean fail) throws IOException{
		final String label=(fail ? "output-failure" : "success");
		final Path diagnostic=dir.resolve(label+".log");
		final Path output=(fail ? Paths.get("/dev/full") : dir.resolve(label+".fq"));
		Throwable failure=null;
		try{
			BloomFilterCorrectorWrapper.main(new String[]{"in="+input, "out="+output, "outstream="+diagnostic,
					"k=31", "ksmall=31", "bits=4", "cells=2000003", "t=1", "seed=0", "qin=33", "qout=33",
					"ecc=f", "merge=f", "ecco=f", "markerrors=f", "mincount=0", "tossjunk=f", "loglog=f",
					"interleaved=f", "pigz=f", "unpigz=f", "ow=t"});
		}catch(Throwable caught){failure=caught;}
		check(Files.size(diagnostic)>0, label+": wrapper must actually open and write its diagnostic file.");
		final int open=openDescriptors(diagnostic);
		check(open==0, label+": diagnostic descriptor remained open after wrapper return/throw.");
		if(fail){
			check(failure instanceof RuntimeException && failure.getMessage()!=null &&
					failure.getMessage().contains("terminated in an error state"),
					"A real output failure must escape wrapper main/run as its processing failure, not success: "+failure);
		}else{
			check(failure==null, "The success control unexpectedly failed: "+failure);
			check(java.util.Arrays.equals(Files.readAllBytes(input), Files.readAllBytes(output)),
					"The successful correction-disabled control must preserve the literal input bytes.");
		}
		System.out.println(label+" failure="+(failure==null ? "none" : failure.getClass().getName())+
				" diagnosticOpenDescriptors="+open);
	}

	/** Counts this JVM's descriptors pointing to the exact retained file, without reading other files. */
	private static int openDescriptors(final Path target) throws IOException{
		final Path wanted=target.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path descriptor : descriptors){
				try{
					if(Files.readSymbolicLink(descriptor).equals(wanted)){count++;}
				}catch(NoSuchFileException vanished){/* A transient JVM descriptor closed during enumeration. */}
			}
		}
		return count;
	}

	/** Keeps every integration assertion active under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
