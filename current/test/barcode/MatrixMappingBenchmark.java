package test.barcode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import barcode.Barcode;
import barcode.DemuxData;
import barcode.PCRMatrixHDist;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Measures repeated Hamming assignment of a fixed 513-name, eight-base input.
 * Includes map allocation, candidate searches, count updates and worker startup/join.
 * This is a helper benchmark, not demultiplexer or file-throughput measurement.
 * Baseline workers accumulate state; fixed mapping replaces it. Both must return
 * the same map and aggregate totals. Use a fresh JVM for each driver and repeat.
 * @author Brian Bushnell, Shinobu
 */
public final class MatrixMappingBenchmark{

	/** Parses mode=baseline|fixed, driver=serial|local|shared and loops; warms up before timing. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MatrixMappingBenchmark.class, false);
		try{
			String mode=null, driver=null;
			int loops=0, warmup=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("driver")){driver=value;}
				else if(key.equals("loops")){loops=Integer.parseInt(value);}
				else if(key.equals("warmup")){warmup=Integer.parseInt(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or fixed.");
			check("serial".equals(driver) || "local".equals(driver) || "shared".equals(driver),
					"Set driver=serial, local or shared.");
			check(loops>0 && loops<=1000000, "Use 1..1000000 mapping calls per measurement.");
			if(warmup<0){warmup=Math.max(64, loops/8);}
			check(warmup>0 && warmup<=1000000, "Use 1..1000000 warmup calls before measurement.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Mapping benchmark requires -ea for production invariants.");
			final boolean serial="serial".equals(driver), baseline="baseline".equals(mode);
			Shared.setThreads(serial ? 1 : 2);
			PCRMatrixHDist.parseStatic("maxhdist=2", "maxhdist", "2");
			PCRMatrixHDist.parseStatic("clearzone=1", "clearzone", "1");
			final Matrix matrix=new Matrix(!"shared".equals(driver));
			final long warmupChecksum=matrix.run(warmup);
			check(warmupChecksum==512L*warmup, "Warmup lost map entries.");
			matrix.verify(baseline ? (serial ? 0 : warmup) : 1);
			final long start=System.nanoTime();
			final long checksum=matrix.run(loops);
			final long elapsed=System.nanoTime()-start;
			matrix.verify(baseline ? (serial ? 0 : warmup+(long)loops) : 1);
			check(checksum==512L*loops && elapsed>0, "Timed calls lost map entries or clock did not advance.");
			blackhole=checksum+matrix.lastMap.hashCode();
			System.out.println("MATRIX_MAPPING_BENCHMARK mode="+mode+" driver="+driver+
					" loops="+loops+" warmup="+warmup+" queries=513 references=256"+
					" checksum="+checksum+" expected="+(512L*loops)+" elapsed_nanos="+elapsed);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private MatrixMappingBenchmark(){}

	/** Owns references independently of observation objects and inspects protected counts. */
	private static final class Matrix extends PCRMatrixHDist{
		Matrix(final boolean local){
			super(8, 0, 0, false);
			localCounts=local;
			matrixThreads=2;
			final char[] bases={'A', 'C', 'G', 'T'};
			for(int code=0; code<256; code++){
				final char[] sequence={'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
				for(int pos=4; pos<8; pos++){sequence[pos]=bases[(code>>((7-pos)*2))&3];}
				final String ref=new String(sequence);
				names.add(ref);
				queries.add(new Barcode(ref, 1));
				queries.add(new Barcode("T"+ref.substring(1), 1));
			}
			queries.add(new Barcode(POLY_C, 5));
			populateExpected(names);
			final DemuxData data=new DemuxData(8, 0, 0);
			data.addPolyA=data.addPolyG=data.addPolyT=false;
			data.addPolyC=true;
			populateSplitCodes(data);
			check(queries.size()==513 && expectedList.size()==256 &&
					queries.size()*(long)expectedList.size()>=100000,
					"Input must cross the worker-dispatch threshold; threads selects serial versus workers.");
		}

		/** Includes one whole mapping operation per pass, consuming every returned map's size. */
		long run(final int loops){
			assert(loops>0) : "The final map must come from at least one measured call.";
			long checksum=0;
			for(int i=0; i<loops; i++){
				lastMap=makeAssignmentMap(queries, 1);
				checksum+=lastMap.size();
			}
			return checksum;
		}

		/** Checks the entire final map, independent positional oracle and borrowed weights outside timing. */
		void verify(final long copies){
			check(lastMap.size()==512 && !lastMap.containsKey(POLY_C), "Unexpected map membership.");
			for(final String ref : names){
				check(ref.equals(lastMap.get(ref)) && ref.equals(lastMap.get("T"+ref.substring(1))),
						"Exact/mutant pair must map to its independently known reference: "+ref);
			}
			check(totalCounted==517 && totalAssigned==517 && totalAssignedToExpected==512 && !errorState,
					"Aggregate weights must describe the last call, not accumulated calls.");
			for(int pos=0; pos<8; pos++){
				for(int call=0; call<5; call++){
					for(int ref=0; ref<5; ref++){
						final long unit=pos==0 ? (ref==0 && (call==0 || call==3) ? 256 : 0) :
								pos<4 ? (call==0 && ref==0 ? 512 : 0) : (call==ref && call<4 ? 128 : 0);
						check(counts[pos][call][ref]==(unit+(call==1 && ref==1 ? 5 : 0))*copies,
								"Unexpected positional state at "+pos+","+call+","+ref+" after "+copies+" sets.");
					}
				}
			}
			for(final Barcode ref : expectedList){check(ref.count()==2*copies, "Expected reference weight changed.");}
			check(allCodesMap.get(POLY_C).count()==5*copies, "Synthetic reference weight changed.");
			for(final Barcode query : queries){
				check(query.count()==(POLY_C.equals(query.name) ? 5 : 1), "Borrowed observation weight changed.");
			}
		}

		private final ArrayList<String> names=new ArrayList<String>(256);
		private final ArrayList<Barcode> queries=new ArrayList<Barcode>(513);
		private HashMap<String, String> lastMap;
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private static final String POLY_C="CCCCCCCC";
	private static volatile long blackhole;
}
