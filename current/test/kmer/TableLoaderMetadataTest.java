package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.HashArray1D;
import kmer.TableLoaderLockFree;
import stream.FastaReadInputStream;

/**
 * Checks exact reference keys, scaffold IDs and retained names with tiny real inputs.
 * Each case uses a fresh loader/JVM and seven independent, nonresizing table shards.
 * Literal keys and metadata are independent of production encoding/remapping helpers.
 * Worker-failure injection and repeated loader calls are outside this fixture.
 * @author Sucrose
 */
public final class TableLoaderMetadataTest{

	/** Accepts one of six case names and a directory for retained, never-overwritten input files. */
	public static void main(String[] args) throws Exception{
		String name=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){name=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		check(name!=null && dir!=null, "Expected case= and dir= arguments");
		final boolean files, keepNames, useRefNames, mixed;
		switch(name){
			case "literal-control": files=false; keepNames=false; useRefNames=false; mixed=false; break;
			case "literal-names": files=false; keepNames=true; useRefNames=false; mixed=false; break;
			case "file-control": files=true; keepNames=false; useRefNames=false; mixed=false; break;
			case "file-names": files=true; keepNames=true; useRefNames=false; mixed=false; break;
			case "file-refnames": files=true; keepNames=true; useRefNames=true; mixed=false; break;
			case "mixed-refnames": files=true; keepNames=true; useRefNames=true; mixed=true; break;
			default: throw new IllegalArgumentException("Unknown case: "+name);
		}
		Files.createDirectories(dir);
		final String[] refs;
		if(files){
			final Path alpha=dir.resolve("alpha.fa"), beta=dir.resolve("beta.fa");
			Files.write(alpha, ">first\nAAA\n>second\nACG\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			Files.write(beta, ">third\nGAA\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			refs=new String[] {alpha.toString(), beta.toString()};
		}else{refs=null;}
		final String[] literals=mixed ? new String[] {"AAC"} : files ? null : new String[] {"AAA", "ACG"};
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		FastaReadInputStream.MIN_READ_LEN=1;
		FastaReadInputStream.SPLIT_READS=false;
		TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
		TableLoaderLockFree.DISPLAY_PROGRESS=false;
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 3, 0, 0, 0, 0, true, false);
		loader.setRefSkip(0);
		loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
		final long added=loader.processData(refs, literals, keepNames, useRefNames, false);
		check(!loader.errorState, "Valid tiny inputs must not report an input-close error");
		System.out.println("case="+name+", added="+added+", scaffoldNames="+loader.scaffoldNames+
			", refScafCounts="+Arrays.toString(loader.refScafCounts));

		//Diagnose each baseline defect before another metadata mismatch can hide it.
		if(keepNames && !files){
			check(tables[1].contains(127) && tables[0].contains(91),
				"KMR114: retaining literal names must not drop AAA or ACG");
		}else if(mixed){
			check(tables[6].contains(111), "KMR114: the mixed file/literal load must retain AAC");
		}else if(useRefNames){
			check(loader.scaffoldNames.contains("alpha") && loader.scaffoldNames.contains("beta"),
				"KMR115: reference-name remapping must run after scaffold counts are populated");
		}else if(keepNames){
			check(tables[1].getValue(127)==1, "KMR115: the first real scaffold must have ID 1, not reserved ID 0");
		}

		final int expected=files ? (mixed ? 4 : 3) : 2;
		verifyKeys(tables, keepNames, files, mixed, expected);
		check(added==expected, "Each distinct reference key must be counted once: "+added+" versus "+expected);
		verifyMetadata(loader, refs, keepNames, useRefNames, files, mixed, expected);
		System.out.println("PASS case="+name+" keys="+expected);
	}

	/** Checks every possible three-base encoding/tag combination, every shard and physical occupancy. */
	private static void verifyKeys(final HashArray1D[] tables, final boolean keepNames,
			final boolean files, final boolean mixed, final int expectedSize){
		//Canonical two-bit values plus length bit64: AAA127, ACG91, GAA125, AAC111.
		for(int key=0; key<128; key++){
			final int id=key==127 ? 1 : key==91 ? 2 : key==125 && files ? 3 : key==111 && mixed ? 4 : -1;
			for(int shard=0; shard<tables.length; shard++){
				final boolean present=id>0 && shard==key%7;
				final int value=present ? (keepNames ? id : 1) : -1;
				check(tables[shard].contains(key)==present, "Unexpected presence for key "+key+" in shard "+shard);
				check(tables[shard].getValue(key)==value, "Wrong scaffold ID for key "+key+" in shard "+shard+
					": "+tables[shard].getValue(key)+" versus "+value);
			}
		}
		long total=0;
		for(int shard=0; shard<tables.length; shard++){
			final HashArray1D table=tables[shard];
			check(table.arrayLength()==77 && table.victims().size()==0, "Tiny fixed tables must not grow or overflow");
			int occupied=0;
			for(long key : table.array()){
				if(key>=0){
					check(key<128 && key%7==shard, "Unexpected physical key or shard: "+key);
					occupied++;
				}
			}
			check(occupied==table.size(), "Physical and reported primary occupancy must agree");
			total+=occupied;
		}
		check(total==expectedSize, "Wrong total distinct key count: "+total+" versus "+expectedSize);
	}

	/** Checks literal metadata order, the reserved slot, lengths and reference grouping. */
	private static void verifyMetadata(final TableLoaderLockFree loader, final String[] refs,
			final boolean keepNames, final boolean useRefNames, final boolean files,
			final boolean mixed, final int expectedSize){
		if(!keepNames){
			check(loader.scaffoldNames==null && loader.scaffoldLengths==null &&
				loader.refNames==null && loader.refScafCounts==null, "Names-off loading must leave metadata absent");
			return;
		}
		final String[] names=!files ? new String[] {"", "1", "2"} :
			!useRefNames ? new String[] {"", "first", "second", "third"} :
			mixed ? new String[] {"", "alpha", "alpha", "beta", "literal"} :
			new String[] {"", "alpha", "alpha", "beta"};
		check(loader.scaffoldNames.equals(Arrays.asList(names)), "Wrong ordered scaffold names: "+loader.scaffoldNames);
		check(loader.scaffoldLengths!=null && loader.scaffoldLengths.size==expectedSize+1,
			"Lengths must include exactly one reserved slot and every real scaffold");
		for(int id=0; id<=expectedSize; id++){
			check(loader.scaffoldLengths.get(id)==(id==0 ? 0 : 3), "Wrong scaffold length at ID "+id);
		}
		final ArrayList<String> expectedRefs=new ArrayList<String>();
		if(files){expectedRefs.addAll(Arrays.asList(refs));}
		if(!files || mixed){expectedRefs.add("literal");}
		check(loader.refNames.equals(expectedRefs), "Reference paths/literal group must retain input order");
		final int[] counts=!files ? new int[] {2} : mixed ? new int[] {2, 1, 1} : new int[] {2, 1};
		check(Arrays.equals(loader.refScafCounts, counts), "Wrong scaffold counts per reference");
	}

	/** Keeps all literal oracles active under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
