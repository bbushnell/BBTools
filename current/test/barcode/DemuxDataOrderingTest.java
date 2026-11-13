package test.barcode;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Locale;

import barcode.Barcode;
import barcode.DemuxData;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite DemuxData ordering regression with a sorted control and unsorted counts.
 * Baseline mode requires the current silent delta-count corruption; fixed mode
 * requires a sorted round trip of the same unsorted input.
 * @author Brian Bushnell, Jean
 */
public final class DemuxDataOrderingTest{

	/** Requires mode=baseline|fixed; assertions must be enabled. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DemuxDataOrderingTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("DemuxData ordering regression requires -ea.");}

			final ArrayList<Barcode> control=roundTrip(sortedCounts());
			checkRows("Sorted control", control, "AC", 2, "GT", 1);
			final ArrayList<Barcode> decoded=roundTrip(unsortedCounts());
			if("baseline".equals(mode)){
				checkRows("Baseline unsorted", decoded, "GT", 1, "AC", 1);
			}else{
				checkRows("Fixed unsorted", decoded, "AC", 2, "GT", 1);
			}
			System.out.println("DEMUXDATA_ORDERING mode="+mode+" control="+control.size()+
					" decoded="+decoded.size()+" secondCount="+decoded.get(1).count());
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each process owns one finite ordering check. */
	private DemuxDataOrderingTest(){}

	/** Already sorted by Barcode.compareTo: descending count, then barcode name. */
	private static ArrayList<Barcode> sortedCounts(){
		final ArrayList<Barcode> list=new ArrayList<Barcode>(2);
		list.add(new Barcode("AC", 2));
		list.add(new Barcode("GT", 1));
		return list;
	}

	/** Intentionally opposite of Barcode.compareTo to exercise encode's ENSURE_SORTED path. */
	private static ArrayList<Barcode> unsortedCounts(){
		final ArrayList<Barcode> list=new ArrayList<Barcode>(2);
		list.add(new Barcode("GT", 1));
		list.add(new Barcode("AC", 2));
		return list;
	}

	/** Encodes a single-index A48 count set and decodes it in the same fresh JVM. */
	private static ArrayList<Barcode> roundTrip(final ArrayList<Barcode> counts){
		final DemuxData data=new DemuxData(2, 0, 0);
		data.expectedList=new LinkedHashSet<String>();
		data.expectedList.add("AC");
		data.expectedList.add("GT");
		data.codeCounts=counts;
		return new ArrayList<Barcode>(new DemuxData(data.encode()).codeCounts);
	}

	/** Confirms literal expected rows and counts without hiding ordering changes. */
	private static void checkRows(final String label, final ArrayList<Barcode> observed,
			final String name0, final long count0, final String name1, final long count1){
		check(observed.size()==2, label+" decoded "+observed.size()+" rows.");
		check(name0.equals(observed.get(0).name) && observed.get(0).count()==count0,
				label+" row 0: expected "+name0+"/"+count0+", observed "+observed.get(0));
		check(name1.equals(observed.get(1).name) && observed.get(1).count()==count1,
				label+" row 1: expected "+name1+"/"+count1+", observed "+observed.get(1));
	}

	/** Always enforce oracle checks, independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
