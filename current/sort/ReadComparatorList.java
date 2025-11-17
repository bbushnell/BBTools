package sort;

import java.io.File;
import java.util.HashMap;

import fileIO.TextFile;
import shared.Shared;
import shared.Tools;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorList extends ReadComparator {
	
	/**
	 * Constructs a ReadComparatorList from a file or comma-separated string.
	 * If fname corresponds to an existing file, reads one identifier per line.
	 * Otherwise, treats fname as a comma-separated list of identifiers.
	 * @param fname File path or comma-separated string of read identifiers
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
	}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	/**
	 * Internal comparison logic for ordering reads by reference list position.
	 * Reads with null IDs or IDs not in the map are sorted after mapped reads.
	 * For reads with equal positions, sorts by pair number.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative, zero, or positive value for less than, equal, or greater than
	 */
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
	
	/** Multiplier for sort direction: 1 for ascending, -1 for descending */
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
	/** Maps read identifiers to their sort order positions */
	private HashMap<String, Integer> map;
	
}
