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
 * Checks assignment maps, aggregate weights, positional counts and reference counts
 * across serial/local/shared worker paths and repeated calls. Fixed mode expresses
 * the proposed replacement-count contract; production has not yet been changed.
 * @author Brian Bushnell, Shinobu
 */
public final class MatrixMappingStateTest{

	/** Requires mode=baseline|fixed, driver=serial|local|shared and seeded=t|f. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MatrixMappingStateTest.class, false);
		try{
			String mode=null, driver=null;
			boolean seeded=false;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("driver")){driver=value;}
				else if(key.equals("seeded")){seeded=Parse.parseBoolean(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("serial".equals(driver) || "local".equals(driver) || "shared".equals(driver),
					"Set driver=serial, local or shared.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Mapping-state checks require -ea.");
			final ArrayList<String> expected=expectedNames();
			final ArrayList<Barcode> queries=queries(expected);
			PCRMatrixHDist.parseStatic("maxhdist=2", "maxhdist", "2");
			PCRMatrixHDist.parseStatic("clearzone=1", "clearzone", "1");
			Shared.setThreads(1);
			final Matrix matrix=new Matrix(expected, !"shared".equals(driver));
			matrix.initializeData();
			matrix.refine(queries, 1);
			matrix.checkCounts(0);
			if(seeded){matrix.populateCounts(queries, 1); matrix.checkCounts(1);}
			final boolean serial="serial".equals(driver), baseline="baseline".equals(mode);
			Shared.setThreads(serial ? 1 : 2);
			for(int pass=1; pass<=2; pass++){
				final HashMap<String, String> map=matrix.makeAssignmentMap(queries, 1);
				check(map.size()==512, "The 512 expected-matching names must be mapped; size="+map.size());
				for(final String ref : expected){
					check(ref.equals(map.get(ref)) && ref.equals(map.get(mutant(ref))),
							"Exact and one-mismatch observations must map to "+ref);
				}
				matrix.checkTotals();
				final int copies=baseline ? (seeded ? 1 : 0)+(serial ? 0 : pass) : 1;
				matrix.checkCounts(copies);
				for(final Barcode query : queries){
					check(query.count()==(POLY_C.equals(query.name) ? 5 : 1), "Mapping mutated a borrowed observation count.");
				}
				System.out.println("MATRIX_MAPPING mode="+mode+" driver="+driver+" seeded="+seeded+
						" pass="+pass+" map="+map.size()+" total=517 positional_copies="+copies+
						" expected_reference_count_each="+(2*copies)+" unexpected_reference_count="+(5*copies));
			}
		}finally{Shared.closeStream(pp.outstream);}
	}

	private MatrixMappingStateTest(){}

	/** All 256 four-base suffixes, each prefixed with AAAA. */
	private static ArrayList<String> expectedNames(){
		final ArrayList<String> names=new ArrayList<String>(256);
		final char[] bases={'A', 'C', 'G', 'T'};
		for(int code=0; code<256; code++){
			final char[] sequence={'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
			int value=code;
			for(int pos=7; pos>=4; pos--){sequence[pos]=bases[value&3]; value>>>=2;}
			names.add(new String(sequence));
		}
		return names;
	}

	/** Exact/one-mismatch observations plus one unexpected homopolymer of weight five. */
	private static ArrayList<Barcode> queries(final ArrayList<String> expected){
		final ArrayList<Barcode> queries=new ArrayList<Barcode>(513);
		for(final String ref : expected){
			queries.add(new Barcode(ref, 1));
			queries.add(new Barcode(mutant(ref), 1));
		}
		queries.add(new Barcode(POLY_C, 5));
		return queries;
	}

	private static String mutant(final String ref){return "T"+ref.substring(1);}

	/** Exposes protected state only for literal independent count checks. */
	private static final class Matrix extends PCRMatrixHDist{
		Matrix(final ArrayList<String> names, final boolean local){
			super(8, 0, 0, false);
			localCounts=local;
			matrixThreads=2;
			populateExpected(names);
			final DemuxData data=new DemuxData(8, 0, 0);
			data.addPolyA=data.addPolyG=data.addPolyT=false;
			data.addPolyC=true;
			populateSplitCodes(data);
			check(513L*expectedList.size()>=100000, "Fixture must cross the production worker-dispatch threshold.");
		}

		void checkTotals(){
			check(totalCounted==517 && totalAssigned==517 && totalAssignedToExpected==512 && !errorState,
					"Mapping totals must describe one completed observation set: "+totalCounted+", "+
					totalAssigned+", "+totalAssignedToExpected+", error="+errorState);
		}

		/** Expected pairs supply 512 observations; the unexpected C homopolymer adds five per position. */
		void checkCounts(final int copies){
			check(counts.length==8 && expectedList.size()==256, "Fixture matrix dimensions changed.");
			for(int pos=0; pos<8; pos++){
				for(int call=0; call<5; call++){
					for(int ref=0; ref<5; ref++){
						final long unit=pos==0 ? (ref==0 && (call==0 || call==3) ? 256 : 0) :
								pos<4 ? (call==0 && ref==0 ? 512 : 0) : (call==ref && call<4 ? 128 : 0);
						final long expected=(unit+(call==1 && ref==1 ? 5 : 0))*copies;
						check(counts[pos][call][ref]==expected,
								"Positional count mismatch at "+pos+","+call+","+ref+": "+counts[pos][call][ref]+
								" expected="+expected);
					}
				}
			}
			for(final Barcode ref : expectedList){
				check(ref.count()==2L*copies, "Reference count mismatch for "+ref.name+": "+ref.count());
			}
			final Barcode unexpected=allCodesMap.get(POLY_C);
			check(unexpected!=null && unexpected.expected==0 && unexpected.count()==5L*copies,
					"Unexpected candidate count must also follow the selected state contract: "+unexpected);
		}
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private static final String POLY_C="CCCCCCCC";
}
