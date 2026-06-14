package sort;

import java.io.File;
import java.util.HashMap;

import fileIO.TextFile;
import shared.Shared;
import shared.Tools;
import stream.Read;


/**
 * Comparator for sorting reads based on a predefined order list.
 * Orders reads according to their identifiers' positions in a reference list,
 * which can be loaded from a file or provided as a comma-separated string.
 * Reads not found in the list are sorted to the end.
 *
 * Unlike the other comparators this is not a shared singleton (it carries a per-instance
 * order map), so direction is supplied via the constructor / {@link #getComparator}.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public final class ReadComparatorList extends ReadComparator {

	/**
	 * Constructs an ascending comparator from a file or comma-separated list of read IDs.
	 * If the path exists, reads one identifier per line; otherwise splits the string on commas.
	 * @param fname File path or comma-separated identifier list
	 */
	public ReadComparatorList(String fname){
		String[] array;
		if(new File(fname).exists()){
			array=TextFile.toStringLines(fname);
		}else{
			array=fname.split(",");
		}
		int mapSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, (array.length*3L)/2);
		map=new HashMap<String, Integer>(mapSize);
		for(int i=0; i<array.length; i++){
			map.put(array[i], i);
		}
		mult=1;
	}

	/** Private constructor sharing an existing order map with a chosen direction. */
	private ReadComparatorList(HashMap<String, Integer> map_, int mult_){
		map=map_;
		mult=mult_;
	}

	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return mult*x;
	}

	/** Compares two reads by their position in the predefined id list; ids absent from the list
	 * sort to the end, with pair number as a tiebreaker. */
	public int compareInner(Read r1, Read r2) {

		Integer a=(r1.id==null ? null : map.get(r1.id));
		Integer b=(r2.id==null ? null : map.get(r2.id));

		if(a==null && b==null){return r1.pairnum()-r2.pairnum();}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int dif=a-b;
		if(dif==0){return r1.pairnum()-r2.pairnum();}
		return dif;
	}

	@Override
	public ReadComparator getComparator(boolean asc){
		final int m=(asc ? 1 : -1);
		return m==mult ? this : new ReadComparatorList(map, m);
	}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "List";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Maps each read id to its index (sort rank) in the predefined order list. */
	private final HashMap<String, Integer> map;

}
