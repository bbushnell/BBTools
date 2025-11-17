package barcode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteStreamWriter;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * Tracks barcode mapping statistics for sequence analysis.
 * Maintains counts of barcode occurrences and cross-references between
 * read barcodes and their mapped reference locations to detect contamination.
 * @author Brian Bushnell
 */
public class BarcodeMappingStats {
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a new barcode mapping statistics tracker with empty collections */
	public BarcodeMappingStats() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges statistics from another BarcodeMappingStats instance into this one.
	 * Combines barcode counts and source mapping data by adding counts together.
	 * @param bs The other statistics instance to merge into this one
	 */
	public void merge(BarcodeMappingStats bs) {
		for(Entry<String, Barcode> e : bs.codeMap.entrySet()) {
			Barcode b=e.getValue();
			incrementCodeMap(b.name, b.count());
		}
		for(Entry<String, HashMap<String, Barcode>> e : bs.sourceMap.entrySet()) {
			final String readKey=e.getKey();
			final HashMap<String, Barcode> map=e.getValue();
			for(Entry<String, Barcode> ee : map.entrySet()) {
				final Barcode b=ee.getValue();
				final String refKey=ee.getKey();
				incrementSourceMap(readKey, refKey, b.count());
			}
		}
	}
	
	/**
	 * Increments barcode statistics for a read and its mapped reference.
	 * Updates both the overall barcode count and the barcode-to-reference mapping.
	 * @param r The read containing barcode information
	 * @param refKey Reference sequence identifier where the read mapped, or null for unknown
	 */
	public void increment(Read r, String refKey){
		String barcode=r.barcode(true);
		incrementCodeMap(barcode, r.pairCount());
		incrementSourceMap(barcode, refKey==null ? "UNKNOWN" : refKey, r.pairCount());
	}
	
	/**
	 * Increments the count for a specific barcode in the overall barcode map.
	 * Creates a new barcode entry if one doesn't exist for the given key.
	 * @param key The barcode string to increment
	 * @param amt The amount to add to the barcode's count
	 */
	public void incrementCodeMap(String key, long amt) {
		Barcode b=codeMap.get(key);
		if(b==null){
			b=new Barcode(key);
			codeMap.put(key, b);
		}
		b.increment(amt);
	}
	
	/**
	 * Increments the count for a barcode-reference mapping in the source map.
	 * Creates nested map structures as needed for new barcode-reference pairs.
	 *
	 * @param readKey The barcode from the read
	 * @param refKey The reference sequence identifier
	 * @param amt The amount to add to this mapping's count
	 */
	public void incrementSourceMap(String readKey, String refKey, long amt) {
		HashMap<String, Barcode> map=sourceMap.get(readKey);
		if(map==null){
			map=new HashMap<String, Barcode>();
			sourceMap.put(readKey, map);
		}
		Barcode b=map.get(refKey);
		if(b==null){
			b=new Barcode(refKey);
			map.put(refKey, b);
		}
		b.increment(amt);
	}

	/**
	 * Writes barcode mapping statistics to a tab-delimited output file.
	 * Output includes barcode, source reference, count, and fraction columns
	 * sorted by barcode count in descending order.
	 *
	 * @param outbarcodes Output file path for the statistics report
	 * @param overwrite Whether to overwrite existing output files
	 */
	public void writeStats(String outbarcodes, boolean overwrite) {
		ByteBuilder bb=new ByteBuilder();
		ArrayList<Barcode> codeList=toSortedList(codeMap);
		final long sum=sum(codeList);
		final double invSum=1.0/(Tools.max(1, sum));

		ByteStreamWriter bsw=new ByteStreamWriter(outbarcodes, overwrite, false, true);
		bsw.start();

		bsw.println("#Reads\t"+sum);
		bsw.println("#Barcode\tSource\tCount\tFraction");
		for(Barcode bc : codeList) {
			HashMap<String, Barcode> map=sourceMap.get(bc.name);
			ArrayList<Barcode> sourceList=toSortedList(map);
			final long sum2=sum(sourceList);
			final double invSum2=1.0/(Tools.max(1, sum2));
			for(Barcode source : sourceList) {
				bb.append(bc.name).tab().append(source.name).tab().append(source.count()).tab().append(source.count()*invSum2, 6).nl();
				bsw.print(bb);
				bb.clear();
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/**
	 * Calculates the total count across all barcodes in the list.
	 * @param list List of barcodes to sum
	 * @return Total count of all barcodes in the list
	 */
	private static long sum(ArrayList<Barcode> list) {
		long sum=0;
		for(Barcode bc : list) {
			sum+=bc.count();
		}
		return sum;
	}
	
	/**
	 * Converts a barcode map to a sorted list ordered by count (descending).
	 * Returns null if the input map is null or empty.
	 * @param map Map of barcode names to Barcode objects
	 * @return Sorted list of barcodes by count, or null if map is empty
	 */
	private static ArrayList<Barcode> toSortedList(HashMap<String, Barcode> map){
		if(map==null || map.isEmpty()){return null;}
		ArrayList<Barcode> list=new ArrayList<Barcode>(map.size());
		for(Entry<String, Barcode> e : map.entrySet()) {
			list.add(e.getValue());
		}
		Collections.sort(list);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Raw counts of barcodes */
	public HashMap<String, Barcode> codeMap=new HashMap<String, Barcode>();
	
	/** 
	 * A table of tables.  Key1 is the barcode of a read; key2 is where the read mapped to.
	 * The barcode for key2 tracks the number of times reads with  key1 barcode mapped to key2's reference. 
	 * E.G. if key1 is ABC-DEF then the top barcode in its table would be expected to
	 * be ABC-DEF, and other entries would indicate contamination of that library.
	 * 
	 */
	public HashMap<String, HashMap<String, Barcode>> sourceMap=new HashMap<String, HashMap<String, Barcode>>();
	
	/** Tracks whether any errors occurred during statistics processing or output */
	public boolean errorState=false;
	
}
