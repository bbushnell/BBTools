package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import kmer.AbstractKmerTableSet;
import kmer.DumpThread;
import kmer.KmerTableSet;
import shared.Shared;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Checks dump-worker failure reporting through both real table-set callers.
 * A null shard deliberately injects an internal worker failure before traversal;
 * it is not a claim that null shards are valid application input. Nearby controls
 * dump the same loaded tables without injection. Run each case in a fresh JVM.
 * @author Shinobu
 */
public final class KmerDumpFailureTest{

	/** Accepts case=short-control|short-failure|long-control|long-failure and a new dir. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(dir==null || !("short-control".equals(mode) || "short-failure".equals(mode) ||
				"long-control".equals(mode) || "long-failure".equals(mode))){
			throw new IllegalArgumentException("Expected a short/long control/failure case and dir=<new directory>");
		}
		final boolean longKeys=mode.startsWith("long"), failure=mode.endsWith("failure");
		Files.createDirectories(dir);
		final Path input=dir.resolve("input.fa"), output=dir.resolve("output.fa");
		final StringBuilder record=new StringBuilder(">input\n");
		for(int i=0; i<63; i++){record.append('A');}
		record.append('\n');
		Files.write(input, record.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		check(!Files.exists(output), "Output must be new so a previous dump cannot satisfy the control");

		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		Kmer.PACKED=true;
		DumpThread.NUM_THREADS=3;
		final String[] options={"in="+input, "k="+(longKeys ? 63 : 3), "ways=3", "threads=3",
			"prefilter=0", "prealloc=f", "onepass=f", "rcomp=t", "minprob=0", "interleaved=f",
			"minavgquality=0", "qtrim=f", "ecco=f", "merge=f", "overwrite=f"};
		final AbstractKmerTableSet tables=longKeys ? new KmerTableSetU(options, 12) : new KmerTableSet(options, 12);
		tables.processInput();
		check(!tables.errorState && tables.readsIn==1 && tables.basesIn==63,
			"Fixture input must finish successfully before the injected dump failure");
		check(Shared.threads()>=3 && tables.tables().length==3,
			"Three shards and workers are required to enter the table-set parallel dump path");
		if(failure){
			if(longKeys){((KmerTableSetU)tables).tables()[0]=null;}
			else{((KmerTableSet)tables).tables()[0]=null;}
		}
		final boolean reported=tables.dumpKmersAsBytes_MT(output.toString(), 1, Integer.MAX_VALUE, false, null);
		check(reported==failure, "Dump worker failure must reach the caller error result: injected="+failure+", reported="+reported);
		if(!failure){check(Files.size(output)>0, "Successful control must emit the loaded k-mer");}
		System.out.println("PASS case="+mode+" reportedError="+reported);
	}

	/** Keeps failure-reporting oracles active with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
