package barcode;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import hiseq.IlluminaHeaderParser2;
import shared.LineParser2;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.SuperLongList;

/**
 * Static utility class for barcode counting and related file operations.
 * Provides methods for reading sequence files, extracting barcodes from read headers,
 * maintaining barcode count statistics, and writing results and mappings to tab-separated output.
 * @author Brian Bushnell
 * @date April 9, 2024
 */
public class BarcodeCounter {
	//TODO: Make non-static and put barcode-length/number/delimiter detection here.

	/*--------------------------------------------------------------*/
	/*----------------           File I/O           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Counts barcodes from a reads file with automatic format detection.
	 * Detects barcode delimiter and barcodes-per-read from FileFormat, then delegates to the main counting method.
	 * @param readsFile Path to the reads file (FASTQ or other supported format)
	 * @param maxReads Maximum number of reads to process (0 = unlimited)
	 * @param addTile Whether to include tile information in barcode keys
	 * @return HashMap mapping barcode strings to Barcode objects with counts
	 */
	public static final HashMap<String, Barcode> countBarcodes(String readsFile, long maxReads, boolean addTile) {
		FileFormat ff=FileFormat.testInput(readsFile, FileFormat.FASTQ, null, true, true);
		final byte delimiter=(byte)ff.barcodeDelimiter();
		final int barcodesPerRead=ff.barcodesPerRead();
		return countBarcodes(ff, maxReads, addTile);
	}
	
	/**
	 * Main barcode counting method that processes reads from a FileFormat.
	 * Creates a concurrent input stream, processes reads in batches, and uses IlluminaHeaderParser2 to extract barcodes.
	 * @param ff FileFormat object containing input file configuration
	 * @param maxReads Maximum number of reads to process (0 = unlimited)
	 * @param addTile Whether to include tile information in barcode keys for uniqueness
	 * @return HashMap mapping barcode strings (optionally with tile) to Barcode objects with counts
	 */
	public static final HashMap<String, Barcode> countBarcodes(FileFormat ff, long maxReads, boolean addTile) {
		final ConcurrentReadInputStream cris=makeCris(ff, maxReads);
		HashMap<String, Barcode> map=new HashMap<String, Barcode>();
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//Check to ensure pairing is as expected
		if(ln!=null && !ln.isEmpty()){
			Read r=ln.get(0);
			assert(r.samline!=null || (r.mate!=null)==cris.paired());
		}
		
		IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
			processList(ln, cris, map, ihp, addTile);

			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}

		errorState|=ReadWrite.closeStream(cris);
		
		return map;
	}
	
	private static final ConcurrentReadInputStream makeCris(FileFormat ff1, long maxReads){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(
				maxReads, true, ff1, null);
		cris.start(); //Start the stream
		return cris;
	}
	
	private static final void processList(ListNum<Read> ln, final ConcurrentReadInputStream cris, 
			HashMap<String, Barcode> codeMap, IlluminaHeaderParser2 ihp, boolean addTile){

		//Grab the actual read list from the ListNum
		final ArrayList<Read> reads=ln.list;
		
		//Loop through each read in the list
		for(int idx=0; idx<reads.size(); idx++){
			final Read r1=reads.get(idx);
			if(!r1.validated()){r1.validate(true);}
			processRead(r1, codeMap, ihp, addTile);
		}

		//Notify the input stream that the list was used
		cris.returnList(ln);
	}
	
	private static final void processRead(final Read r1, HashMap<String, Barcode> codeMap, IlluminaHeaderParser2 ihp, boolean addTile){
		ihp.parse(r1.id);
		final String key=ihp.barcode();
		if(key==null){return;}
		String key2=key;
		int tile=0;
		if(addTile) {
			tile=ihp.tile();
			assert(tile>10 && tile<10000);
			key2+=tile;
		}
		Barcode b=codeMap.get(key2);
		if(b==null){
			b=new Barcode(key);
			b.tile=tile;
			codeMap.put(key2, b);
		}
		b.increment(1);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads expected barcode sequences from a file or comma-separated string.
	 * If the input is a valid file path, reads barcodes from file (one per line); otherwise treats it as a comma-separated list.
	 * Filters out comment lines starting with '#' and normalizes barcodes via BarcodeStats.fixBarcode.
	 * @param fname File path or comma-separated barcode string
	 * @param forceDelimiter Delimiter character to force in barcode formatting
	 * @return ArrayList of formatted barcode strings, or null if input is null
	 */
	public static final ArrayList<String> loadBarcodes(String fname, int forceDelimiter){
		if(fname==null){return null;}
		String[] codes;
		if(new File(fname).exists()) {
			codes=TextFile.toStringLines(fname);
		}else {
			codes=fname.split(",");
		}
		ArrayList<String> expected=new ArrayList<String>();
		for(int i=0; i<codes.length; i++){
			String s=codes[i];
			if(!Tools.startsWith(s, '#')){
//				s=removeTab(s);
//				assert(s.indexOf('\t')<0) : "Barcodes should not contain a tab: '"+s+"'";//Although it's fine if they do; just use -da
				s=BarcodeStats.fixBarcode(s, forceDelimiter, false, false);
				expected.add(s);
			}
		}
		return expected;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads barcode count data from a tab-separated file.
	 * Each non-comment line must contain barcode name and count separated by tab; filters by a minimum count threshold.
	 * @param fname Path to tab-separated barcode count file
	 * @param minCount Minimum count threshold for inclusion
	 * @return ArrayList of Barcode objects meeting the count threshold
	 */
	public static final ArrayList<Barcode> loadCounts(String fname, long minCount){
		LineParser2 lp=new LineParser2('\t');
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		ArrayList<Barcode> list=new ArrayList<Barcode>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line[0]=='#') {continue;}
			lp.set(line);
			String a=lp.parseString();
			assert(lp.hasMore()) : new String(line);
			long b=lp.parseLong();
			if(b>=minCount) {
				Barcode c=new Barcode(a, b);
				list.add(c);
			}
		}
		errorState=bf.close()|errorState;
		return list;
	}
	
	/**
	 * Writes barcode counts to file using FileFormat settings.
	 * Delegates to the String-filename writeCounts method, using overwrite and append flags from the FileFormat.
	 * @param counts Collection of Barcode objects to write
	 * @param minCount Minimum count threshold for inclusion
	 * @param ff FileFormat specifying output settings
	 * @param sort Whether to sort barcodes before writing
	 * @param overwrite Whether to overwrite existing files
	 * @param append Whether to append to existing files
	 * @return true if writing succeeded, false otherwise
	 */
	public static final boolean writeCounts(Collection<Barcode> counts, long minCount, FileFormat ff, 
			boolean sort, boolean overwrite, boolean append){
		return writeCounts(counts, minCount, ff, sort, ff.overwrite(), ff.append());
	}
	
	/**
	 * Writes barcode counts to a tab-separated file with header information.
	 * Optionally sorts barcodes, writes summary header lines, and filters output by a minimum count threshold.
	 * @param counts Collection of Barcode objects to write
	 * @param minCount Minimum count threshold for inclusion
	 * @param fname Output file path
	 * @param sort Whether to sort barcodes by count before writing
	 * @param overwrite Whether to overwrite existing files
	 * @param append Whether to append to existing files
	 * @return true if writing succeeded, false otherwise
	 */
	public static final boolean writeCounts(Collection<Barcode> counts, long minCount, String fname, 
			boolean sort, boolean overwrite, boolean append){
		if(sort) {
			ArrayList<Barcode> list=new ArrayList<Barcode>(counts);
			Collections.sort(list);
			counts=list;
		}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		
		long sum=0;
		for(Barcode bc : counts) {sum+=bc.count();}
		bsw.println("#Barcodes\t"+sum);
		bsw.println("#Unique\t"+counts.size());
		ByteBuilder bb=new ByteBuilder(128);
		for(Barcode b : counts) {
			if(b.count()>=minCount) {
				b.appendTo(bb).nl();
				bsw.print(bb);
				bb.clear();
			}
//			if(b.count()>=minCount) {
//				bsw.print(b.name).tab().print(b.count()).nl();
//			}
		}
		boolean b=bsw.poisonAndWait();
		errorState|=b;
		return b;
	}
	
	/**
	 * Writes barcode counts from a HashMap to a tab-separated file.
	 * Outputs barcode name and count per line without additional filtering or sorting.
	 * @param counts HashMap of barcode names to Barcode objects
	 * @param fname Output file path
	 * @param overwrite Whether to overwrite existing files
	 * @param append Whether to append to existing files
	 * @return true if writing succeeded, false otherwise
	 */
	public static final boolean writeCounts(HashMap<String, Barcode> counts, String fname, 
			boolean overwrite, boolean append){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		for(Entry<String, Barcode> e : counts.entrySet()) {
			Barcode b=e.getValue();
			bsw.print(b.name).tab().print(b.count()).nl();
		}
		boolean b=bsw.poisonAndWait();
		errorState|=b;
		return b;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads a barcode assignment map from a tab-separated file.
	 * Each non-comment line must contain two tab-separated values representing key-value pairs.
	 * @param mapIn Path to tab-separated assignment map file
	 * @return HashMap mapping first column values to second column values
	 */
	public static HashMap<String, String> loadAssignmentMap(String mapIn) {
		LineParser2 lp=new LineParser2('\t');
		ByteFile bf=ByteFile.makeByteFile(mapIn, true);
		HashMap<String, String> map=new HashMap<String, String>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line[0]=='#') {continue;}
			lp.set(line);
			String a=lp.parseString();
			assert(lp.hasMore()) : new String(line);
			String b=lp.parseString();
			map.put(a, b);
		}
		errorState=bf.close()|errorState;
		return map;
	}
	
	/**
	 * Writes a barcode assignment map to a tab-separated file.
	 * Emits one line per entry with key and value separated by a tab.
	 * @param assignmentMap HashMap containing key-value assignments to write
	 * @param mapOut Output file path
	 * @param overwrite Whether to overwrite existing files
	 * @param append Whether to append to existing files
	 * @return true if writing succeeded, false otherwise
	 */
	public static boolean writeAssignmentMap(HashMap<String, String> assignmentMap,
			String mapOut, boolean overwrite, boolean append) {

		ByteStreamWriter bsw=new ByteStreamWriter(mapOut, overwrite, append, true);
		bsw.start();
		for(Entry<String, String> e : assignmentMap.entrySet()) {
			String a=e.getKey(), b=e.getValue();
			bsw.print(a).tab().print(b);
			bsw.nl();
		}
		boolean b=bsw.poisonAndWait();
		errorState|=b;
		return b;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a sorted SuperLongList of barcode counts for statistical analysis.
	 * Increments the list by each barcode count and sorts the list in ascending order.
	 * @param barcodes Collection of Barcode objects
	 * @return Sorted SuperLongList containing all barcode counts
	 */
	public static SuperLongList makeCountList(Collection<Barcode> barcodes) {
		SuperLongList sll=new SuperLongList(100000);
		for(Barcode b : barcodes) {
			sll.increment(b.count());
		}
		sll.sort();
		return sll;
	}
	
	/**
	 * Returns barcode count at a specified percentile threshold using cumulative count.
	 * Computes the value from the sorted SuperLongList and currently asserts for testing rather than production use.
	 * @param barcodes Collection of Barcode objects for analysis
	 * @param percentile Percentile threshold (0.0 to 1.0)
	 * @return Count value at the specified percentile (asserts in current implementation)
	 */
	public static long barcodeCountPercentile(Collection<Barcode> barcodes, float percentile) {
		//Just an example; don't call it directly.
		SuperLongList sll=makeCountList(barcodes);
		long x=sll.percentileValueBySum(percentile);//TODO: Note - this is untested and gives assertion error.
		assert(false) : x; //Make sure this returns the count described in the annotation, not +1 or -1.
		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean errorState=false;
	
}
