package test.kmer;

import java.lang.reflect.Field;
import java.util.Arrays;

import kmer.AbstractKmerTable;
import kmer.HashArray;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.TableLoaderLockFree;
import structures.IntList;

/**
 * Checks repeat-load statistics and retained-metadata reuse boundaries.
 * Names-off calls may append to existing caller-owned values. Retained metadata
 * calls require all borrowed destination tables to be empty, including HashArray
 * victim forests, and must reject before replacing old metadata or adding keys.
 * @author Jean
 */
public final class TableLoaderReuseTest{

	/** Runs one KMR121 case in a fresh JVM. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=names-off-counters|retained-existing|"+
				"prepopulated-primary|prepopulated-victim|cleared-control.");
		}
		final boolean oldReplicate=TableLoaderLockFree.REPLICATE_AMBIGUOUS;
		final boolean oldProgress=TableLoaderLockFree.DISPLAY_PROGRESS;
		final String mode=args[0].substring(5);
		try{
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
			TableLoaderLockFree.DISPLAY_PROGRESS=false;
			if("names-off-counters".equals(mode)){checkNamesOffCounters(mode);}
			else if("retained-existing".equals(mode)){checkRetainedExistingRejection(mode);}
			else if("prepopulated-primary".equals(mode)){checkPrepopulatedPrimaryRejection(mode);}
			else if("prepopulated-victim".equals(mode)){checkPrepopulatedVictimRejection(mode);}
			else if("cleared-control".equals(mode)){checkRetainedNamesAfterExplicitClear(mode);}
			else{throw new IllegalArgumentException("Unknown case: "+mode);}
			System.out.println("PASS table-loader-reuse case="+mode);
		}finally{
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=oldReplicate;
			TableLoaderLockFree.DISPLAY_PROGRESS=oldProgress;
		}
	}

	/** Eight AAA roots followed by one CCC root must report only the second invocation. */
	private static void checkNamesOffCounters(final String mode) throws ReflectiveOperationException{
		final AbstractKmerTable[] tables=hashArrays();
		final TableLoaderLockFree loader=loader(tables);
		final long first=loader.processData(null, copies("AAA", 8), false, false, false);
		fixture(first==1, mode, "Eight identical AAA roots should create one key, not "+first);
		checkCounters(loader, mode, 8, 24, 8, "first call");

		final long second=loader.processData(null, new String[] {"CCC"}, false, false, false);
		kmr121(second==1, mode, "One new CCC root should create one key, not "+second);
		checkCounters(loader, mode, 1, 3, 1, "second call");
		kmr121(countKeys(tables)==2, mode, "Names-off append should retain AAA and add CCC");
		kmr121(tables[(int)(KEY_AAA%7)].getValue(KEY_AAA)==1, mode, "AAA key missing after append");
		kmr121(tables[(int)(KEY_CCC%7)].getValue(KEY_CCC)==1, mode, "CCC key missing after append");
		kmr121(!loader.errorState, mode, "Names-off append should not set errorState");
		checkNoMetadata(loader, mode);
	}

	/** A second retained-metadata call must fail before metadata replacement or table mutation. */
	private static void checkRetainedExistingRejection(final String mode){
		final AbstractKmerTable[] tables=hashArrays();
		final TableLoaderLockFree loader=loader(tables);
		final long first=loader.processData(null, new String[] {"AAA"}, true, false, false);
		fixture(first==1 && countKeys(tables)==1, mode, "Initial retained-name load failed");
		loader.errorState=true;
		final Object[] state=metadataState(loader);
		fixture(loader.scaffoldNames.equals(Arrays.asList("", "1")), mode,
			"Wrong initial names: "+loader.scaffoldNames);
		fixture(loader.scaffoldLengths.get(1)==3 && loader.refScafCounts[0]==1, mode,
			"Wrong initial retained metadata");

		expectRetainedNamesRejected(loader, new String[] {"CCC"}, mode);
		checkMetadataState(loader, state, mode);
		kmr121(loader.errorState, mode, "Rejected retained-name call must preserve errorState");
		kmr121(countKeys(tables)==1 && tables[(int)(KEY_AAA%7)].getValue(KEY_AAA)==1 &&
			!tables[(int)(KEY_CCC%7)].contains(KEY_CCC),
			mode, "Rejected retained-name call must preserve old keys and avoid CCC");
	}

	/** Retained-name loading must reject a fresh loader with primary prepopulation. */
	private static void checkPrepopulatedPrimaryRejection(final String mode){
		final AbstractKmerTable[] tables=hashArrays();
		tables[(int)(KEY_AAA%7)].set(KEY_AAA, 7);
		final TableLoaderLockFree loader=loader(tables);
		loader.errorState=true;
		final Object[] state=metadataState(loader);
		expectRetainedNamesRejected(loader, new String[] {"CCC"}, mode);
		checkMetadataState(loader, state, mode);
		kmr121(loader.errorState, mode, "Rejected fresh-loader call must preserve errorState");
		kmr121(tables[(int)(KEY_AAA%7)].getValue(KEY_AAA)==7 &&
			!tables[(int)(KEY_CCC%7)].contains(KEY_CCC),
			mode, "Rejected fresh-loader call must preserve preloaded primary keys");
	}

	/** Retained-name loading must reject a fresh loader with victim-only prepopulation. */
	private static void checkPrepopulatedVictimRejection(final String mode){
		final AbstractKmerTable[] tables=hashArrays();
		final HashArray victimShard=(HashArray)tables[0];
		victimShard.victims().set(KEY_VICTIM, 9);
		fixture(victimShard.size()==0 && victimShard.victims().size()==1,
			mode, "Expected no primary keys and one victim key");
		final TableLoaderLockFree loader=loader(tables);
		final Object[] state=metadataState(loader);

		expectRetainedNamesRejected(loader, new String[] {"CCC"}, mode);
		checkMetadataState(loader, state, mode);
		kmr121(victimShard.size()==0 && victimShard.victims().size()==1 &&
			victimShard.victims().getValue(KEY_VICTIM)==9,
			mode, "Rejected fresh-loader call must preserve victim-only occupancy");
	}

	/** When the same tables are explicitly emptied, retained IDs can safely restart at one. */
	private static void checkRetainedNamesAfterExplicitClear(final String mode) throws ReflectiveOperationException{
		final HashForest[] tables=forests();
		final TableLoaderLockFree loader=loader(tables);
		loader.processData(null, new String[] {"AAA"}, true, false, false);
		fixture(countKeys(tables)==1, mode, "Initial retained-name load should store one key");
		for(HashForest table : tables){table.clear();}
		fixture(countKeys(tables)==0, mode, "Explicit forest clear must empty every shard");

		final long added=loader.processData(null, new String[] {"CCC"}, true, false, false);
		kmr121(added==1, mode, "Retained-name load after explicit clear should add CCC, not "+added);
		kmr121(countKeys(tables)==1 && !tables[(int)(KEY_AAA%7)].contains(KEY_AAA) &&
			tables[(int)(KEY_CCC%7)].getValue(KEY_CCC)==1,
			mode, "Cleared tables should contain only fresh CCC with scaffold ID 1");
		kmr121(loader.scaffoldNames.equals(Arrays.asList("", "1")) &&
			loader.scaffoldLengths.get(1)==3 && loader.refScafCounts[0]==1,
			mode, "Fresh retained metadata should restart after explicit clear");
		checkCounters(loader, mode, 1, 3, 1, "after explicit clear");
		kmr121(!loader.errorState, mode, "Explicit-clear control should not set errorState");
	}

	/** Constructs a loader with exact canonical k=3 keys and caller-owned destination tables. */
	private static TableLoaderLockFree loader(final AbstractKmerTable[] tables){
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 3);
		loader.setRefSkip(0);
		loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
		return loader;
	}

	/** Allocates seven fixed HashArray shards. */
	private static HashArray1D[] hashArrays(){
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		return tables;
	}

	/** Allocates seven clearable fixed forest shards. */
	private static HashForest[] forests(){
		final HashForest[] tables=new HashForest[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashForest(17, false);}
		return tables;
	}

	/** Counts primary plus victim keys because HashArray.size() excludes overflow. */
	private static long countKeys(final AbstractKmerTable[] tables){
		long sum=0;
		for(AbstractKmerTable table : tables){
			sum+=table.size();
			if(table instanceof HashArray){sum+=((HashArray)table).victims().size();}
		}
		return sum;
	}

	/** Expects retained metadata to reject nonempty destination tables. */
	private static void expectRetainedNamesRejected(final TableLoaderLockFree loader,
			final String[] literal, final String mode){
		try{
			loader.processData(null, literal, true, false, false);
			throw new AssertionError("KMR121 case="+mode+
				": Expected retained-name reuse on nonempty tables to be rejected");
		}catch(IllegalStateException expected){
			// Expected.
		}
	}

	/** Captures identity and content for metadata fields before an expected rejection. */
	private static Object[] metadataState(final TableLoaderLockFree loader){
		final int[] counts=loader.refScafCounts;
		final int[] countsCopy=counts==null ? null : counts.clone();
		final String[] names=loader.scaffoldNames==null ? null : loader.scaffoldNames.toArray(new String[0]);
		final String[] refs=loader.refNames==null ? null : loader.refNames.toArray(new String[0]);
		final int[] lengths=new int[loader.scaffoldLengths.size];
		for(int i=0; i<lengths.length; i++){lengths[i]=loader.scaffoldLengths.get(i);}
		return new Object[] {loader.scaffoldNames, names, loader.refNames, refs,
			counts, countsCopy, loader.scaffoldLengths, lengths, loader.errorState};
	}

	/** Checks that a rejected call changed neither metadata identities nor contents. */
	private static void checkMetadataState(final TableLoaderLockFree loader, final Object[] state,
			final String mode){
		final int[] expectedCounts=(int[])state[5];
		kmr121(loader.scaffoldNames==state[0] && loader.refNames==state[2] &&
			loader.refScafCounts==state[4] && loader.scaffoldLengths==state[6],
			mode, "Rejected call must preserve metadata field identities");
		kmr121(Arrays.equals(loader.scaffoldNames==null ? null : loader.scaffoldNames.toArray(new String[0]),
			(String[])state[1]), mode, "Rejected call changed scaffoldNames");
		kmr121(Arrays.equals(loader.refNames==null ? null : loader.refNames.toArray(new String[0]),
			(String[])state[3]), mode, "Rejected call changed refNames");
		kmr121(Arrays.equals(loader.refScafCounts, expectedCounts),
			mode, "Rejected call changed refScafCounts");
		final IntList lengths=loader.scaffoldLengths;
		final int[] oldLengths=(int[])state[7];
		kmr121(lengths.size==oldLengths.length, mode, "Rejected call changed scaffoldLengths size");
		for(int i=0; i<oldLengths.length; i++){
			kmr121(lengths.get(i)==oldLengths[i], mode, "Rejected call changed scaffold length "+i);
		}
		kmr121(loader.errorState==((Boolean)state[8]).booleanValue(),
			mode, "Rejected call changed errorState");
	}

	/** Reads existing package-private aggregate counters without production visibility changes. */
	private static long statistic(final TableLoaderLockFree loader, final String name)
			throws ReflectiveOperationException{
		final Field field=TableLoaderLockFree.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(loader);
	}

	/** Checks the three per-invocation aggregate counters. */
	private static void checkCounters(final TableLoaderLockFree loader, final String mode,
			final long reads, final long bases, final long kmers, final String phase)
			throws ReflectiveOperationException{
		kmr121(statistic(loader, "refReads")==reads &&
			statistic(loader, "refBases")==bases &&
			statistic(loader, "refKmers")==kmers,
			mode, "Wrong "+phase+" counters: "+statistic(loader, "refReads")+"/"+
			statistic(loader, "refBases")+"/"+statistic(loader, "refKmers"));
	}

	/** Returns a literal array filled with one sequence. */
	private static String[] copies(final String s, final int count){
		final String[] array=new String[count];
		Arrays.fill(array, s);
		return array;
	}

	/** Checks that a names-off post-load state has no retained metadata. */
	private static void checkNoMetadata(final TableLoaderLockFree loader, final String mode){
		kmr121(loader.scaffoldNames==null && loader.refNames==null && loader.refScafCounts==null &&
			loader.scaffoldLengths==null,
			mode, "Names-off loading must leave retained metadata absent");
	}

	/** Expected failures of the proposed KMR121 production repair. */
	private static void kmr121(final boolean condition, final String mode, final String message){
		if(!condition){throw new AssertionError("KMR121 case="+mode+": "+message);}
	}

	/** Fixture setup checks that are not themselves KMR121 production witnesses. */
	private static void fixture(final boolean condition, final String mode, final String message){
		if(!condition){throw new AssertionError("Fixture qualification case="+mode+": "+message);}
	}

	private static final long KEY_AAA=127;
	private static final long KEY_CCC=106;
	private static final long KEY_VICTIM=0;
}
