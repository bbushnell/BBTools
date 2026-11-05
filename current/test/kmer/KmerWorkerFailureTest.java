package test.kmer;

import java.util.Arrays;

import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.HistogramMaker;
import shared.Shared;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HistogramMakerU;
import ukmer.Kmer;

/** Checks that table-worker failures reach the calling thread, with real-table controls.
 * A null table deterministically causes the same delegated failure in sequential
 * and worker paths; the test never depends on timing or worker assignment order.
 * @author Shinobu */
public final class KmerWorkerFailureTest {

	/** Runs one explicit control or rejection group in a fresh JVM. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=controls|short-hist|long-hist|short-init|short-clear|long-init|long-clear.");
		}
		final String name=args[0].substring(5);
		final int oldThreads=Shared.threads();
		try{
			Shared.setThreads(4);
			if(name.equals("controls")){controls();}
			else{
				boolean rejected=false;
				try{fail(name);}
				catch(RuntimeException e){
					for(Throwable cause=e; cause!=null; cause=cause.getCause()){
						if(cause instanceof NullPointerException){rejected=true; break;}
					}
					if(!rejected){throw e;}
				}
				check(rejected, "Worker failure was hidden from the caller: "+name);
			}
		}finally{Shared.setThreads(oldThreads);}
		System.out.println("PASS "+name);
	}

	/** Positive controls exercise both histogram strategies and all ownership operations. */
	private static void controls(){
		final HashArray1D[] shortTables=shortTables();
		final HashArrayU1D[] longTables=longTables();
		final long[] expected={0, 0, 4, 4, 0, 0};
		for(int threads : new int[] {1, 4}){
			Shared.setThreads(threads);
			check(Arrays.equals(expected, HistogramMaker.fillHistogram(shortTables, 5)), "Short histogram control differs.");
			check(Arrays.equals(expected, HistogramMakerU.fillHistogram(longTables, 5)), "Long histogram control differs.");
			kmer.OwnershipThread.initialize(shortTables);
			ukmer.OwnershipThread.initialize(longTables);
			for(HashArray1D table : shortTables){
				check(table.getOwner(1)==-1 && table.setOwner(1, 2)==2, "Short ownership control differs.");
			}
			for(HashArrayU1D table : longTables){
				final Kmer key=key(0);
				check(table.getOwner(key)==-1 && table.setOwner(key, 2)==2, "Long ownership control differs.");
			}
			kmer.OwnershipThread.clear(shortTables);
			ukmer.OwnershipThread.clear(longTables);
			kmer.OwnershipThread.initialize(shortTables);
			ukmer.OwnershipThread.initialize(longTables);
			for(HashArray1D table : shortTables){check(table.getOwner(1)==-1 && table.getValue(1)==2, "Short reset changed data.");}
			for(HashArrayU1D table : longTables){check(table.getOwner(key(0))==-1 && table.getValue(key(0))==2, "Long reset changed data.");}
			kmer.OwnershipThread.clear(shortTables);
			ukmer.OwnershipThread.clear(longTables);
		}
		// The single-table fast paths already reject a bad delegate synchronously.
		int rejected=0;
		try{kmer.OwnershipThread.initialize(new AbstractKmerTable[] {null});}catch(NullPointerException e){rejected++;}
		try{ukmer.OwnershipThread.initialize(new AbstractKmerTableU[] {null});}catch(NullPointerException e){rejected++;}
		Shared.setThreads(1);
		try{HistogramMaker.fillHistogram(new AbstractKmerTable[] {null}, 5);}catch(NullPointerException e){rejected++;}
		try{HistogramMakerU.fillHistogram(new AbstractKmerTableU[] {null}, 5);}catch(NullPointerException e){rejected++;}
		check(rejected==4, "Sequential controls must expose the same null-table error.");
	}

	/** A fixed four-table shape forces worker dispatch; exactly one delegate is invalid. */
	private static void fail(final String name){
		final HashArray1D[] shorts=shortTables();
		final HashArrayU1D[] longs=longTables();
		shorts[2]=null;
		longs[2]=null;
		if(name.equals("short-hist")){HistogramMaker.fillHistogram(shorts, 5);}
		else if(name.equals("long-hist")){HistogramMakerU.fillHistogram(longs, 5);}
		else if(name.equals("short-init")){kmer.OwnershipThread.initialize(shorts);}
		else if(name.equals("short-clear")){kmer.OwnershipThread.clear(shorts);}
		else if(name.equals("long-init")){ukmer.OwnershipThread.initialize(longs);}
		else if(name.equals("long-clear")){ukmer.OwnershipThread.clear(longs);}
		else{throw new IllegalArgumentException("Unknown case: "+name);}
	}

	/** Four independent small tables avoid introducing shared-table mutation in the fixture. */
	private static HashArray1D[] shortTables(){
		final HashArray1D[] tables=new HashArray1D[4];
		for(int i=0; i<tables.length; i++){
			tables[i]=new HashArray1D(new int[] {101}, -1L);
			tables[i].set(1, 2);
			tables[i].set(2, 3);
		}
		return tables;
	}

	/** Long-key equivalent with literal count bins2 and3 in each table. */
	private static HashArrayU1D[] longTables(){
		final HashArrayU1D[] tables=new HashArrayU1D[4];
		final Kmer a=key(0), c=key(1);
		for(int i=0; i<tables.length; i++){
			tables[i]=new HashArrayU1D(new int[] {101}, a.k, a.kbig);
			tables[i].set(a, 2);
			tables[i].set(c, 3);
		}
		return tables;
	}

	/** Independent homopolymer keys are sufficient for this orchestration-only fixture. */
	private static Kmer key(final int base){
		final Kmer key=new Kmer(63);
		for(int i=0; i<key.kbig; i++){key.addRightNumeric(base);}
		return key;
	}

	/** Remains effective under both assertion modes. */
	private static void check(final boolean ok, final String message){
		if(!ok){throw new AssertionError(message);}
	}
}
