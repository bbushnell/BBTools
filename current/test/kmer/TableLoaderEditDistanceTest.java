package test.kmer;

import java.lang.reflect.Field;

import kmer.HashArray1D;
import kmer.TableLoaderLockFree;

/** Separates full edit dispatch from short-key indel eligibility with literal key sets.
 * Seven tiny shards and one literal read exercise the actual loader. Expected keys
 * use forward two-bit bases plus the length marker; no production mutation helper
 * generates an oracle. Run one case per fresh JVM under either assertion mode.
 * @author Sucrose
 */
public final class TableLoaderEditDistanceTest{

	/** Runs three independent defect cases and three unchanged nearby controls. */
	public static void main(String[] args) throws ReflectiveOperationException{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=<edit-distance case> argument.");
		}
		final String mode=args[0].substring(5);
		final String literal;
		final int k, mink, hdist, edist, shortHdist, shortEdist, expectedCount;
		final long[] fullKeys, shortKeys;
		switch(mode){
			case "full-edits-without-hamming":
				literal="AAA"; k=3; mink=0; hdist=0; edist=1; shortHdist=shortEdist=0;
				fullKeys=AAA_EDITS; shortKeys=EMPTY; expectedCount=10; break;
			case "short-edits-without-full-edits":
				literal="ACGT"; k=4; mink=3; hdist=edist=0; shortHdist=shortEdist=1;
				fullKeys=ACGT_EXACT; shortKeys=SHORT_EDITS; expectedCount=26; break;
			case "short-hamming-with-full-edits":
				literal="ACGT"; k=4; mink=3; hdist=edist=1; shortHdist=1; shortEdist=0;
				fullKeys=ACGT_EDITS; shortKeys=SHORT_HAMMING; expectedCount=39; break;
			case "exact-control":
				literal="AAA"; k=3; mink=0; hdist=edist=shortHdist=shortEdist=0;
				fullKeys=AAA_EXACT; shortKeys=EMPTY; expectedCount=1; break;
			case "short-hamming-control":
				literal="ACGT"; k=4; mink=3; hdist=edist=shortEdist=0; shortHdist=1;
				fullKeys=ACGT_EXACT; shortKeys=SHORT_HAMMING; expectedCount=21; break;
			case "both-edits-control":
				literal="ACGT"; k=4; mink=3; hdist=edist=shortHdist=shortEdist=1;
				fullKeys=ACGT_EDITS; shortKeys=SHORT_EDITS; expectedCount=44; break;
			default: throw new IllegalArgumentException("Unknown edit-distance case: "+mode);
		}
		check(fullKeys.length+shortKeys.length==expectedCount, mode,
			"Literal expected arrays must agree with the independently specified total");
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		check(storedKeys(tables, mode)==0, mode, "Fixture requires seven empty tiny shards");
		TableLoaderLockFree.DISPLAY_PROGRESS=false;
		TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, k, mink, 0, hdist, edist, false, false);
		loader.hammingDistance2=shortHdist;
		loader.editDistance2=shortEdist;
		loader.setRefSkip(0);
		loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
		final long added=loader.processData(null, new String[] {literal}, false, false, false);
		final long stored=storedKeys(tables, mode);
		check(!loader.errorState && statistic(loader, "refReads")==1
			&& statistic(loader, "refBases")==literal.length() && statistic(loader, "refKmers")==1,
			mode, "The real loader must process exactly one literal and one full window without failure");
		check(value(tables, k==3 ? 64 : 283)==1, mode, "Original full key must retain ID1");
		if(mink>0){check(value(tables, 70)==1 && value(tables, 91)==1, mode, "Original ACG/CGT terminal keys must retain ID1");}

		//Witnesses precede size checks so baseline failures name the dispatch/indel defect.
		if(mode.equals("full-edits-without-hamming") || mode.equals("short-edits-without-full-edits")){
			check(value(tables, 65)==1, mode, "Edit request must store literal AAC key65");
		}else if(mode.equals("short-hamming-with-full-edits")){
			check(value(tables, 65)==-1, mode, "Short hamming-only request must exclude insertion-only AAC key65");
		}
		checkKeys(tables, fullKeys, mode);
		checkKeys(tables, shortKeys, mode);
		check(stored==expectedCount && added==expectedCount, mode,
			"Exact literal key set and creation return must agree: expected="+expectedCount+", stored="+stored+", added="+added);
		System.out.println("PASS loader-edits case="+mode+" keys="+expectedCount);
	}

	/** Counts primary entries only after proving no resize or victim overflow occurred. */
	private static long storedKeys(final HashArray1D[] tables, final String mode){
		long stored=0;
		for(HashArray1D table : tables){
			check(table.arrayLength()==77 && table.victims().size()==0, mode,
				"Fixture must keep 77 cells per shard and no victim entries: cells="+table.arrayLength()+", victims="+table.victims().size());
			stored+=table.size();
		}
		return stored;
	}

	/** Verifies literal membership and stored ID without deriving keys from loader code. */
	private static void checkKeys(final HashArray1D[] tables, final long[] keys, final String mode){
		long previous=-1;
		for(long key : keys){
			check(key>previous, mode, "Literal key arrays must be sorted and duplicate-free");
			check(value(tables, key)==1, mode, "Expected ID1 for literal key "+key);
			previous=key;
		}
	}

	/** Queries the production loader's fixed key-modulo-seven shard assignment. */
	private static int value(final HashArray1D[] tables, final long key){return tables[(int)(key%7)].getValue(key);}

	/** Reads existing package-private aggregate counters without changing production visibility. */
	private static long statistic(final TableLoaderLockFree loader, final String name) throws ReflectiveOperationException{
		final Field field=TableLoaderLockFree.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(loader);
	}

	/** Keeps all fixture and expected-key checks active under -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR122 case="+mode+": "+message);}
	}

	private static final long[] EMPTY={};
	private static final long[] AAA_EXACT={64};
	/** AAA and its nine substitutions; one-edit insert/truncate results add no new keys. */
	private static final long[] AAA_EDITS={64,65,66,67,68,72,76,80,96,112};
	private static final long[] ACGT_EXACT={283};
	/** ACGT plus twelve substitutions and six additional internal insertion/truncation results. */
	private static final long[] ACGT_EDITS={262,267,274,275,278,279,280,281,282,283,286,287,294,299,310,315,347,411,475};
	/** Union of the one-substitution neighborhoods of terminal ACG and CGT. */
	private static final long[] SHORT_HAMMING={66,68,69,70,71,74,75,78,83,86,87,88,89,90,91,95,102,107,118,123};
	/** Adds AAC, AGC, ATC, CAG and CTG to the short substitution union. */
	private static final long[] SHORT_EDITS={65,66,68,69,70,71,73,74,75,77,78,82,83,86,87,88,89,90,91,94,95,102,107,118,123};
}
