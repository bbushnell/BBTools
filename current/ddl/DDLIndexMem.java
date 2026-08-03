package ddl;

import java.util.ArrayList;

import cardinality.DynamicDemiLog;
import parse.Parse;
import shared.Timer;

/**
 * Directly measures the heap footprint of a DDL sketch index. Loads a sketch file, snapshots
 * used heap, builds the index (matrix or CSR via {@code csr=t/f}), snapshots again, and reports
 * the delta -- so matrix vs CSR memory can be compared apples-to-apples on the same data. Uses
 * System.gc() around the snapshots (a hint, but adequate for a coarse comparative number).
 *
 * Usage: DDLIndexMem <sketch.tsv.gz> [csr=f] [k=25] [t=32]
 *
 * @author Noire
 * @date August 3, 2026
 */
public class DDLIndexMem {

	public static void main(String[] args){
		String inFile=null;
		int k=25, threads=32;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("csr") || a.equals("ddlcsr") || a.equals("packedindex")){DDLIndexBase.USE_CSR=Parse.parseBoolean(b);}
			else if(a.equals("csr2") || a.equals("index2") || a.equals("packed2")){DDLIndexBase.USE_CSR2=Parse.parseBoolean(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(inFile==null){inFile=arg;}
		}
		if(inFile==null){System.err.println("Usage: DDLIndexMem <sketch.tsv.gz> [csr=f] [k=25] [t=32]"); System.exit(1);}

		final Runtime rt=Runtime.getRuntime();
		System.err.println("Loading "+inFile+"...");
		final ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(inFile, k, threads);
		final int buckets=records.get(0).ddl.buckets;
		System.err.println("Loaded "+records.size()+" records ("+buckets+" buckets).");

		gc();
		final long before=rt.totalMemory()-rt.freeMemory();

		final Timer bt=new Timer();
		final DDLIndexBase index=DDLIndexBase.create(buckets);
		index.addAll(records, threads);
		bt.stop();

		gc();
		final long after=rt.totalMemory()-rt.freeMemory();

		System.err.println("=== index="+index.getClass().getSimpleName()+"  csr="+DDLIndexBase.USE_CSR
			+"  csr2="+DDLIndexBase.USE_CSR2+"  buckets="+buckets+" ===");
		System.err.println("populatedCells:     "+index.populatedCells());
		System.err.println("Index build time:   "+bt);
		System.err.println("Heap w/ records:    "+mb(before)+" MB");
		System.err.println("Heap w/ +index:     "+mb(after)+" MB");
		System.err.println("INDEX HEAP DELTA:   "+mb(after-before)+" MB");
		//Reference the index so JIT can't elide the build before measurement.
		if(index.numClades()<0){System.err.println("unreachable "+records.size());}
	}

	private static void gc(){
		for(int i=0; i<3; i++){System.gc(); try{Thread.sleep(200);}catch(InterruptedException e){}}
	}

	private static long mb(long bytes){return bytes/(1L<<20);}
}
