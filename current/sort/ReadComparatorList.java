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
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public final class ReadComparatorList extends ReadComparator {
	
	/**
	 * Constructs a comparator from a file or comma-separated list of read IDs.
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
	}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
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
	
	private int ascending=1;
	
	/** Sets the sort direction for read comparison.
	 * @param asc true for ascending order, false for descending order */
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
	@Override
	public final int ascendingMult() {return ascending;}
	@Override
	public final String name() {return "List";}
	
	private HashMap<String, Integer> map;
	
}
