package test.kmer;

import java.nio.charset.StandardCharsets;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.Shared;
import structures.LongList;

/**
 * Checks replacement of a reused forward-kmer list for inputs around k=3.
 * Literal two-bit windows provide the oracle; no hash tables or reads are loaded.
 * Run one case per fresh JVM with either assertion mode.
 * @author Sucrose
 */
public final class KmerTableSetFillKmersTest{

	/** Accepts case=short-input|empty-input|valid-control. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=name argument.");
		}
		final String mode=args[0].substring(5);
		final String input;
		switch(mode){
			case "short-input": input="AC"; break;
			case "empty-input": input=""; break;
			case "valid-control": input="TTA"; break;
			default: throw new IllegalArgumentException("Unknown fillKmers case: "+mode);
		}

		Shared.AMINO_IN=false;
		AbstractKmerTableSet.DISPLAY_PROGRESS=false;
		final KmerTableSet tables=new KmerTableSet(new String[] {
			"k=3", "ways=1", "initialsize=17", "threads=1", "prefilter=0",
			"prealloc=f", "rcomp=f", "minprob=0"}, 12);
		check(tables.k==3 && tables.ways==1 && tables.initialSize==17 && !tables.amino
			&& !tables.prefilter && tables.prefilterArray==null && tables.tables()==null,
			mode, "Fixture requires three-base DNA without allocated tables or a prefilter.");
		final LongList kmers=new LongList(2);
		final int seeded=tables.fillKmers(bytes("ACGT"), kmers);
		//A=0, C=1, G=2, T=3: ACG=6 and CGT=27, independent of table lookup.
		check(seeded==2 && kmers.size==2 && kmers.get(0)==6 && kmers.get(1)==27,
			mode, "Fixture must start with the literal windows [6,27].");

		final int valid=tables.fillKmers(bytes(input), kmers);
		if(mode.equals("valid-control")){
			check(valid==1 && kmers.size==1 && kmers.get(0)==60,
				mode, "Exactly-k input TTA must replace the old list with [60].");
		}else{
			check(valid==0, mode, "Short input must report zero valid windows: "+valid);
			check(kmers.size==0, mode, "Expected empty destination, got size="+kmers.size);
		}
		check(tables.tables()==null && tables.prefilterArray==null,
			mode, "Window extraction must not allocate hash tables or a prefilter.");
		System.out.println("PASS KMR010 case="+mode);
	}

	/** Encodes only the literal ASCII inputs; no production kmer encoder supplies expected values. */
	private static byte[] bytes(final String bases){return bases.getBytes(StandardCharsets.US_ASCII);}

	/** Keeps fixture and oracle checks active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR010 case="+mode+": "+message);}
	}
}
