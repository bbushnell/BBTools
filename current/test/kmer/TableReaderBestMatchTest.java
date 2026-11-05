package test.kmer;

import java.nio.charset.StandardCharsets;

import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.TableReader;
import stream.Read;

/**
 * Finite reused-reader checks for TableReader.findBestMatch sparse counters.
 * Each case uses two literal three-base keys and an independent count array.
 * Run one case per fresh JVM under both assertion modes.
 * @author Jean
 */
public final class TableReaderBestMatchTest{

	/** Selects a direct cleanup reproduction, stale-state result check or unaffected control. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=partial-cleanup|reused-reader|fresh-control");
		}
		final String mode=args[0].substring(5);
		if(mode.equals("partial-cleanup")){partialCleanup();}
		else if(mode.equals("reused-reader")){reusedReader();}
		else if(mode.equals("fresh-control")){freshControl();}
		else{throw new IllegalArgumentException("Unknown case: "+mode);}
		System.out.println("PASS case="+mode);
	}

	/** A failed below-threshold lookup must not leave its touched ID counter dirty. */
	private static void partialCleanup(){
		final Fixture fixture=new Fixture();
		final int result=fixture.reader.findBestMatch(read("AAA"), fixture.tables);
		check(result==-1, "A single hit must remain below minHits=2: "+result);
		check(fixture.reader.countArray[1]==0,
			"KMR112: the partial ID 1 hit was not cleared: "+fixture.reader.countArray[1]);
		check(fixture.reader.countArray[2]==0, "Unmatched ID 2 must stay zero");
	}

	/** A failed lookup must not change the later winner relative to a fresh reader. */
	private static void reusedReader(){
		final Fixture reused=new Fixture();
		check(reused.reader.findBestMatch(read("AAA"), reused.tables)==-1,
			"The first lookup must inject only a below-threshold partial hit");
		final int actual=reused.reader.findBestMatch(read("AAAACG"), reused.tables);
		final Fixture fresh=new Fixture();
		final int expected=fresh.reader.findBestMatch(read("AAAACG"), fresh.tables);
		check(expected==1, "The fresh known-answer read must prefer ID 1: "+expected);
		check(actual==expected, "KMR112: stale state changed best match from "+expected+" to "+actual);
	}

	/** Control: the later read has two AAA hits and one ACG hit in a fresh reader. */
	private static void freshControl(){
		final Fixture fixture=new Fixture();
		final int result=fixture.reader.findBestMatch(read("AAAACG"), fixture.tables);
		check(result==1, "The two-hit ID 1 must beat the one-hit ID 2: "+result);
		check(fixture.reader.countArray[1]==0 && fixture.reader.countArray[2]==0,
			"A successful lookup must clear all packed counters");
	}

	/** Constructs one reader and seven fixed shards containing AAA->1 and ACG->2. */
	private static final class Fixture{
		Fixture(){
			reader=new TableReader(3, 0, 0, 0, 0, false, false);
			reader.minHits=2;
			reader.countArray=new int[3];
			tables=new AbstractKmerTable[7];
			for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
			check(tables[1].set(64, 1)==1, "AAA must be installed as ID 1");
			check(tables[0].set(70, 2)==1, "ACG must be installed as ID 2");
			for(int i=0; i<tables.length; i++){
				final HashArray1D table=(HashArray1D)tables[i];
				check(table.arrayLength()==77 && table.victims().size()==0,
					"Fixture shards must stay fixed and victim-free");
			}
			check(tables[1].getValue(64)==1 && tables[0].getValue(70)==2,
				"Literal encoded keys must retrieve their intended IDs");
		}
		final TableReader reader;
		final AbstractKmerTable[] tables;
	}

	/** Creates an unpaired read over the supplied literal bases. */
	private static Read read(final String bases){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, bases, 0);
	}

	/** Keeps the independent oracle enabled when production assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
