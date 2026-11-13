package test.barcode;

import java.util.Locale;

import barcode.BarcodeStats;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite full-distance histogram checks with literal weights and unique counts.
 * Uses calcEdist=false; the edit histogram mirrors Hamming distance here.
 * Storage and report selectors use separate cases; one mode/distance per fresh JVM.
 * No files, workers or production changes are needed.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeDistanceHistogramTest{

	/** Requires mode=baseline|fixed; case=storage (default) uses8/29/30/31, report uses3/1000. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BarcodeDistanceHistogramTest.class, false);
		try{
			String mode=null;
			String testCase="storage";
			int distance=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("distance")){distance=Integer.parseInt(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("storage".equals(testCase) || "report".equals(testCase), "Set case=storage or report.");
			final boolean report=testCase.equals("report");
			check(report ? distance==3 || distance==1000 : distance==8 || distance==29 || distance==30 || distance==31,
					"Storage distances are8/29/30/31; report distances are3/1000.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Run the production and fixture checks with -ea.");
			Shared.setThreads(1);
			BarcodeStats.calcEdist=false;
			checkDistance(mode, distance, report);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Disjoint alternating alphabets give exactly the selected distance without N or homopolymer shortcuts. */
	private static void checkDistance(final String mode, final int distance, final boolean report){
		final char[] reference=new char[distance], observed=new char[distance];
		for(int i=0; i<distance; i++){
			reference[i]=((i&1)==0 ? 'G' : 'T');
			observed[i]=((i&1)==0 ? 'A' : 'C');
		}
		final String expectedName=new String(reference), unexpectedName=new String(observed);
		//The public constructor/increment API does not require optional loader length metadata.
		final BarcodeStats stats=new BarcodeStats(0, 1, null);
		check(stats.addExpectedCode(expectedName)==1, "The fresh expected set must retain its one name.");
		stats.increment(expectedName, 5);
		stats.increment(unexpectedName, 7);
		final boolean expectFailure=!report && mode.equals("baseline") && distance>=30;
		ArrayIndexOutOfBoundsException observedFailure=null;
		try{stats.calcStats();}
		catch(final ArrayIndexOutOfBoundsException failure){observedFailure=failure;}
		check((observedFailure!=null)==expectFailure, "Unexpected distance-bound outcome for "+mode+" distance="+distance+": "+observedFailure);
		if(observedFailure!=null){
			final StackTraceElement[] frames=observedFailure.getStackTrace();
			check(observedFailure.getClass()==ArrayIndexOutOfBoundsException.class && frames.length>0 &&
					frames[0].getClassName().equals("barcode.BarcodeStats") && frames[0].getMethodName().equals("calcStats"),
					"The baseline must fail at the histogram write in BarcodeStats.calcStats: "+observedFailure);
			//calcStats is incremental and HashMap iteration order is unspecified; discard this partial instance.
			System.out.println("DISTANCE_HISTOGRAM mode="+mode+" distance="+distance+" bounds_failure=true calcStats_frame=true");
			return;
		}
		check(stats.totalCodes==12 && stats.totalCodesU==2 && stats.expectedCodes==5 && stats.expectedCodesU==1,
				"Literal totals must be weight12/unique2 and expected weight5/unique1: "+stats.totalCodes+","+
				stats.totalCodesU+","+stats.expectedCodes+","+stats.expectedCodesU);
		check(stats.nCodes==0 && stats.nCodesU==0, "Fully defined names must not enter the undefined-base bucket.");
		for(final long value : stats.polymerArray){check(value==0, "Alternating names must not enter a homopolymer bucket.");}
		for(final long value : stats.polymerArrayU){check(value==0, "Alternating names must not enter a unique homopolymer bucket.");}
		check(stats.codeMap.size()==2 && stats.codeMap.get(expectedName).count()==5 && stats.codeMap.get(unexpectedName).count()==7,
				"Observed names and their original weights must remain intact.");
		check(stats.expectedCodeMap.size()==1 && stats.expectedCodeMap.get(expectedName).count()==5,
				"The one expected name must receive exactly its five observations.");
		checkHistogram(stats.hdistArray, distance, 7, "Hamming weight");
		checkHistogram(stats.edistArray, distance, 7, "Mirrored edit weight");
		checkHistogram(stats.hdistArrayU, distance, 1, "Hamming unique");
		checkHistogram(stats.edistArrayU, distance, 1, "Mirrored edit unique");
		check(stats.hdistSum(0, Integer.MAX_VALUE)==7 && stats.edistSum(0, Integer.MAX_VALUE)==7,
				"Distance-range reporting must preserve all seven unexpected observations.");
		if(report){checkReport(stats, mode, distance); return;}
		System.out.println("DISTANCE_HISTOGRAM mode="+mode+" distance="+distance+
				" bounds_failure=false total_weight=12 total_unique=2 expected_weight=5 expected_unique=1 bin_weight=7 bin_unique=1");
	}

	/** Checks the actual toStats text after storing the full distance; only the new1000 case differs. */
	private static void checkReport(final BarcodeStats stats, final String mode, final int distance){
		final long displayed=mode.equals("baseline") && distance==1000 ? 0 : 7;
		final String expected="-hDist 1,2,3+\t0\t0\t"+displayed;
		int found=0;
		for(final String line : stats.toStats("Barcodes:").split("\n")){
			if(line.startsWith("-hDist 1,2,3+\t")){
				check(line.equals(expected), "Exact three-plus report mismatch: "+line+" expected="+expected);
				found++;
			}
		}
		check(found==1, "The single-index report must contain exactly one Hamming-distance row; found="+found);
		System.out.println("DISTANCE_REPORT mode="+mode+" distance="+distance+
				" stored_weight=7 stored_unique=1 displayed_3plus="+displayed);
	}

	/** Requires the actual distance bin and zero in every other bin; clamping or dropping cannot pass. */
	private static void checkHistogram(final long[] histogram, final int distance, final long weight, final String name){
		check(histogram.length>distance, name+" must represent the full distance "+distance+"; length="+histogram.length);
		for(int bin=0; bin<histogram.length; bin++){
			final long expected=bin==distance ? weight : 0;
			check(histogram[bin]==expected, name+" bin="+bin+" observed="+histogram[bin]+" expected="+expected);
		}
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private BarcodeDistanceHistogramTest(){}
}
