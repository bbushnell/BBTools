package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import kmer.AbstractKmerTable;
import kmer.KmerTable;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Finite lookup, singleton-buffer and collision-chain dump contracts for KmerTable.
 * Run each case in a fresh JVM under both assertion modes. Dump oracles are literal
 * bytes; positive chain nodes isolate the table's head gate from link pruning.
 * @author Sucrose
 */
public final class KmerTableContractTest{

	/** Accepts a case name and new output directory through case= and dir= arguments. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(mode==null || dir==null){throw new IllegalArgumentException("Expected case=<case name> and dir=<new directory>");}
		switch(mode){
			case "scalar-absent": scalar(false); break;
			case "scalar-present": scalar(true); break;
			case "singleton-present": singleton(true); break;
			case "singleton-absent": singleton(false); break;
			case "dump-tail": dump(dir, 1, 5, 2, 10, ">5\nAAT\n", false); break;
			case "dump-head": dump(dir, 5, 1, 2, 10, ">5\nAAA\n", false); break;
			case "dump-none": dump(dir, 1, 1, 2, 10, "", false); break;
			case "dump-max": dump(dir, 11, 5, 2, 10, ">5\nAAT\n", false); break;
			case "dump-two": dump(dir, 2, 5, 2, 5, ">2\nAAA\n>5\nAAT\n", false); break;
			case "dump-empty": dump(dir, 0, 0, 2, 10, "", true); break;
			default: throw new IllegalArgumentException("Unknown case: "+mode);
		}
		System.out.println("PASS case="+mode);
	}

	/** Checks absence separately from stored zero and positive values (KMR006). */
	private static void scalar(final boolean present){
		final KmerTable table=table();
		check(table.set(0, 0)==1 && table.set(3, 5)==1, "Both distinct fixture keys must be new");
		check(table.size()==2 && table.contains(0) && table.contains(3) && !table.contains(6),
			"Fixture must retain two colliding keys and exclude the third");
		if(present){
			check(table.getValue(0)==0 && table.getValue(3)==5, "Stored zero and positive counts must remain distinguishable from absence");
		}else{
			final int value=table.getValue(6);
			check(value==-1, "KMR006: AbstractKmerTable requires -1 for an absent key; got "+value);
		}
	}

	/** Checks the concrete scalar-copy API with a valid one-element buffer (KMR005/006). */
	private static void singleton(final boolean present){
		final KmerTable table=table();
		check(table.set(0, 5)==1 && table.size()==1 && table.contains(0) && !table.contains(3),
			"Singleton fixture must contain only key zero");
		final int[] buffer={99};
		final int[] result=table.getValues(present ? 0 : 3, buffer);
		check(result==buffer, "KmerTable must reuse the supplied scalar buffer");
		final int expected=present ? 5 : -1;
		check(buffer[0]==expected, "KMR005/006: singleton result expected "+expected+", got "+buffer[0]);
	}

	/** Dumps a two-node positive chain or an empty table, closing the writer before checking bytes. */
	private static void dump(final Path dir, final int head, final int tail, final int min, final int max,
			final String expected, final boolean empty) throws Exception{
		final KmerTable table=table();
		if(!empty){
			check(head>0 && tail>0, "Positive fixture nodes isolate KMR007 from nonpositive-link pruning");
			check(table.set(0, head)==1 && table.set(3, tail)==1, "Insert head AAA before colliding tail AAT");
			check(table.size()==2 && table.contains(0) && table.contains(3), "Dump fixture must retain both keys");
		}else{check(table.size()==0, "Empty control must have no entries");}
		AbstractKmerTable.FASTA_DUMP=true;
		AbstractKmerTable.NUMERIC_DUMP=false;
		Shared.AMINO_IN=false;
		Files.createDirectories(dir);
		final Path output=dir.resolve("dump.fa");
		check(!Files.exists(output), "Dump output must be new; stale bytes cannot satisfy the oracle");
		final ByteStreamWriter writer=new ByteStreamWriter(output.toString(), false, false, false);
		final ByteBuilder buffer=new ByteBuilder();
		writer.start();
		final boolean success;
		final boolean writerError;
		try{
			success=table.dumpKmersAsBytes_MT(writer, buffer, 3, min, max, null);
			writer.print(buffer);
		}finally{
			writerError=writer.poisonAndWait();
		}
		check(success && !writerError, "Table traversal and writer completion must both succeed");
		final byte[] actual=Files.readAllBytes(output);
		check(Arrays.equals(actual, expected.getBytes(StandardCharsets.US_ASCII)),
			"KMR007/count bounds: expected "+expected.replace("\n", "\\n")+", got "+
			new String(actual, StandardCharsets.US_ASCII).replace("\n", "\\n"));
	}

	/** Fixes the modulus at three so keys zero and three share a bucket without resizing. */
	private static KmerTable table(){
		final KmerTable table=new KmerTable(3, false);
		check(table.arrayLength()==3, "Three buckets are required for the collision fixture");
		return table;
	}

	/** Keeps test oracles active under -da and supplies a concrete failed contract. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
