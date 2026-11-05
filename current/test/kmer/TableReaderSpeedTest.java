package test.kmer;

import java.lang.reflect.Field;

import kmer.HashArray1D;
import kmer.TableLoaderLockFree;
import kmer.TableReader;

/** Finite checks for the reader's agreed modulo17 sampling contract.
 * Literal three-base encodings and expected IDs are independent of production
 * filtering and encoding helpers. Prepopulated keys isolate reader filtering;
 * a separate real-loader control preserves unfiltered mutation storage.
 * Run one case per fresh JVM under either assertion mode.
 * @author Sucrose
 */
public final class TableReaderSpeedTest{

	/** Runs three mismatched-selector witnesses and six unchanged controls. */
	public static void main(String[] args) throws ReflectiveOperationException{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=<reader-speed case> argument.");
		}
		final String mode=args[0].substring(5);
		final int speed, skip;
		switch(mode){
			case "tac-speed-one": speed=1; skip=0; break;
			case "tcg-speed-sixteen": speed=16; skip=0; break;
			case "gcg-rejected": speed=1; skip=0; break;
			case "zero-control": speed=0; skip=0; break;
			case "negative-control": speed=-1; skip=0; break;
			case "seventeen-control": speed=17; skip=0; break;
			case "qskip-reject-control": speed=0; skip=2; break;
			case "qskip-pass-control": speed=0; skip=2; break;
			case "mutation-storage-control":
				mutationControl(mode);
				System.out.println("PASS reader-speed case="+mode+" keys=10");
				return;
			default: throw new IllegalArgumentException("Unknown reader-speed case: "+mode);
		}
		final HashArray1D[] tables=emptyTables(mode);
		for(long key : QUERY_KEYS){
			fixture(tables[(int)(key%7)].set(key, 1)==1, mode, "Each literal query key must be newly stored");
		}
		verifyStorage(tables, QUERY_KEYS, mode);
		final TableReader reader=new TableReader(3, 0, speed, skip, 0, true, false);
		switch(mode){
			case "tac-speed-one":
				witness(reader.getValue(49, 44, 3, 0, 0, tables), 1, mode, "TAC key113 must pass speed1"); break;
			case "tcg-speed-sixteen":
				witness(reader.getValue(54, 24, 3, 0, 0, tables), 1, mode, "TCG key118 must pass speed16"); break;
			case "gcg-rejected":
				witness(reader.getValue(38, 25, 3, 0, 0, tables), -1, mode, "GCG key102 must fail speed1 even when stored"); break;
			case "qskip-reject-control":
				fixture(reader.getValue(38, 25, 3, 0, 3, tables)==-1, mode, "qSkip2 must skip query position3 with speed0"); break;
			case "qskip-pass-control":
				fixture(reader.getValue(38, 25, 3, 0, 4, tables)==1, mode, "qSkip2 must retain query position4 with speed0"); break;
			default:
				final int expected=mode.equals("seventeen-control") ? -1 : 1;
				fixture(reader.getValue(49, 44, 3, 0, 0, tables)==expected
					&& reader.getValue(54, 24, 3, 0, 0, tables)==expected
					&& reader.getValue(38, 25, 3, 0, 0, tables)==expected, mode,
					"Nonpositive speed must retain all three literal keys; speed17 must reject all three");
		}
		verifyStorage(tables, QUERY_KEYS, mode);
		System.out.println("PASS reader-speed case="+mode);
	}

	/** Retains the complete literal GCG substitution neighborhood even at nonzero speed. */
	private static void mutationControl(final String mode) throws ReflectiveOperationException{
		final HashArray1D[] tables=emptyTables(mode);
		final boolean oldProgress=TableLoaderLockFree.DISPLAY_PROGRESS;
		final boolean oldAmbiguous=TableLoaderLockFree.REPLICATE_AMBIGUOUS;
		try{
			TableLoaderLockFree.DISPLAY_PROGRESS=false;
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
			final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 3, 0, 1, 1, 0, true, false);
			loader.setRefSkip(0);
			loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
			final long added=loader.processData(null, new String[] {"GCG"}, false, false, false);
			fixture(!loader.errorState && added==10 && statistic(loader, "refReads")==1
				&& statistic(loader, "refBases")==3 && statistic(loader, "refKmers")==1, mode,
				"One GCG read/window must create exactly ten canonical substitution keys");
			verifyStorage(tables, MUTATION_KEYS, mode);
			fixture(tables[4].getValue(102)==1, mode, "Mutation storage must retain key102 despite speed1");
			final TableReader unsampled=new TableReader(3, 0, 0, 0, 0, true, false);
			fixture(unsampled.getValue(38, 25, 3, 0, 0, tables)==1, mode,
				"An unsampled reader must still retrieve the retained original GCG key");
		}finally{
			TableLoaderLockFree.DISPLAY_PROGRESS=oldProgress;
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=oldAmbiguous;
		}
	}

	/** Allocates seven fixed tiny shards and checks their initial physical state. */
	private static HashArray1D[] emptyTables(final String mode){
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		verifyStorage(tables, EMPTY, mode);
		return tables;
	}

	/** Verifies the exact literal set, IDs, physical size and absence of resize/overflow. */
	private static void verifyStorage(final HashArray1D[] tables, final long[] keys, final String mode){
		long previous=-1;
		for(long key : keys){
			fixture(key>previous && tables[(int)(key%7)].getValue(key)==1, mode,
				"Independent literal keys must be sorted, unique and stored with ID1: "+key);
			previous=key;
		}
		long total=0;
		for(int shard=0; shard<7; shard++){
			final HashArray1D table=tables[shard];
			fixture(table.arrayLength()==77 && table.victims().size()==0, mode, "Tiny shards must not grow or overflow");
			long occupied=0;
			for(long key : table.array()){
				if(key>=0){
					fixture(key%7==shard, mode, "Every stored key must occupy its modulo-seven shard");
					occupied++;
				}
			}
			fixture(occupied==table.size(), mode, "Physical occupancy must match reported primary size");
			total+=occupied;
		}
		fixture(total==keys.length, mode, "No keys outside the independent literal set may be stored");
	}

	/** Reads existing aggregate counters without changing production visibility. */
	private static long statistic(final TableLoaderLockFree loader, final String name) throws ReflectiveOperationException{
		final Field field=TableLoaderLockFree.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(loader);
	}

	/** Distinguishes the intended selector mismatch from a failed fixture precondition. */
	private static void witness(final int actual, final int expected, final String mode, final String message){
		if(actual!=expected){throw new AssertionError("KMR125 case="+mode+": "+message+"; expected="+expected+", actual="+actual);}
	}

	/** Keeps all control and geometry checks enabled under -da. */
	private static void fixture(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("Fixture qualification case="+mode+": "+message);}
	}

	private static final long[] EMPTY={};
	/** GCG:38/25->102, TAC:49/44->113, TCG:54/24->118; canonical payload plus length marker64. */
	private static final long[] QUERY_KEYS={102,113,118};
	/** CCG,ACG,GAG,GCG,GCT,GCC,GGG,GTG,TCG,GCA after canonicalization and length marker64. */
	private static final long[] MUTATION_KEYS={90,91,98,102,103,105,106,110,118,121};
}
