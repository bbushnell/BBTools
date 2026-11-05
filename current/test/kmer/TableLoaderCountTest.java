package test.kmer;

import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.TableLoaderLockFree;

/** Checks that increment-mode loading reports new keys rather than accumulated counts.
 * Literal AAA windows and their nine single-base substitutions provide independent
 * key/count oracles. Each case uses seven tiny shards in a fresh JVM.
 * @author Shinobu
 */
public final class TableLoaderCountTest{

	/** Runs a repeated/previously stored key case or a nearby successful control. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=exact-repeat|exact-existing|mutation-existing-zero|mutation-control|single-control.");
		}
		final String mode=args[0].substring(5);
		final boolean existing=mode.equals("exact-existing"), zero=mode.equals("mutation-existing-zero");
		final boolean mutation=zero || mode.equals("mutation-control");
		final boolean single=mode.equals("single-control");
		if(!existing && !mutation && !single && !mode.equals("exact-repeat")){
			throw new IllegalArgumentException("Unknown loader-count case: "+mode);
		}
		final AbstractKmerTable[] tables=new AbstractKmerTable[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		final long[] keys=mutation ? new long[] {64, 65, 66, 67, 68, 72, 76, 80, 96, 112} : new long[] {64};
		if(existing){tables[1].set(64, 5);}
		if(zero){
			for(long key : keys){
				tables[(int)(key%7)].set(key, 0);
				check(tables[(int)(key%7)].getValue(key)==0, mode, "Zero-count key was not present before loading");
			}
		}
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 3, 0, 0, mutation ? 1 : 0, 0, false, false);
		loader.storeMode(TableLoaderLockFree.INCREMENT);
		TableLoaderLockFree.DISPLAY_PROGRESS=false;
		TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
		final int copies=existing ? 2 : single ? 1 : 4;
		final String[] literal=new String[copies];
		for(int i=0; i<literal.length; i++){literal[i]="AAA";}
		final long added=loader.processData(null, literal, false, false, false);
		final int expectedCount=existing ? 7 : copies;
		long stored=0;
		for(AbstractKmerTable table : tables){stored+=table.size();}
		check(!loader.errorState && stored==keys.length, mode, "Unexpected failure or number of stored keys");
		for(long key : keys){
			check(tables[(int)(key%7)].getValue(key)==expectedCount, mode, "Unexpected count for literal key "+key);
		}
		final long expectedAdded=(existing || zero) ? 0 : keys.length;
		check(added==expectedAdded, mode, "Expected newly created keys="+expectedAdded+", got="+added);
		System.out.println("PASS loader-count case="+mode);
	}

	/** Keeps the independent oracle active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR124 case="+mode+": "+message);}
	}
}
