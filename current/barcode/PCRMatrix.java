package barcode;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import shared.Parse;
import shared.Tools;
import sketch.Sketch;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Tracks data about bar code mismatches by position.
 * Used for demultiplexing.
 * 
 * @author Brian Bushnell
 * @contributor Chloe
 * @date March 7, 2024
 *
 */
public abstract class PCRMatrix {

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a PCRMatrix with specified barcode structure parameters.
	 * Initializes counting matrix and validates licensing.
	 *
	 * @param length1_ Length of first barcode segment
	 * @param length2_ Length of second barcode segment (0 for single barcodes)
	 * @param delimiter_ Delimiter character between segments (0 for none)
	 * @param hdistSum_ Whether to sum hamming distances or take maximum
	 */
	protected PCRMatrix(int length1_, int length2_, int delimiter_, boolean hdistSum_) {
		writelock();
		if(!valid()) {
			KillSwitch.killTraceless("This software does not have a valid license.");
		}
		length1=length1_;
		length2=length2_;
		delimiter=delimiter_;
		hdistSum=hdistSum_ || length2<1;
		
		length=length1+length2+(delimiter>0 ? 1 : 0);
		delimiterPos=(delimiter>0 ? length1 : -1);
		start2=(length2>0 ? length1+(delimiter>0 ? 1 : 0) : -1);
		letters=length1+length2;
		
		assert(length>0) : length+", "+length1+", "+length2+", "+delimiter;
		assert(length1>0) : length+", "+length1+", "+length2+", "+delimiter;
		assert(delimiter==0 || length2>0) : length+", "+length1+", "+length2+", "+delimiter;
		
		counts=new long[length][5][5];
		writeunlock();
	}
	
	/**
	 * Factory method to create PCRMatrix with default settings.
	 *
	 * @param length1_ Length of first barcode segment
	 * @param length2_ Length of second barcode segment
	 * @param delimiter_ Delimiter character between segments
	 * @return New PCRMatrix instance of the configured type
	 */
	public static PCRMatrix create(int length1_, int length2_, int delimiter_) {
		return create(matrixType0, length1_, length2_, delimiter_, hdistSum0);
	}
	
	/**
	 * Factory method to create PCRMatrix of specified type.
	 *
	 * @param type Matrix implementation type (HDIST_TYPE, PROB_TYPE, or TILE_TYPE)
	 * @param length1_ Length of first barcode segment
	 * @param length2_ Length of second barcode segment
	 * @param delimiter_ Delimiter character between segments
	 * @param hdistSum_ Whether to sum hamming distances or take maximum
	 * @return New PCRMatrix instance of the specified type
	 * @throws RuntimeException if type is unknown
	 */
	public static PCRMatrix create(int type, int length1_, int length2_, int delimiter_, boolean hdistSum_) {
		if(type==HDIST_TYPE) {
			return new PCRMatrixHDist(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==PROB_TYPE) {
			return createProbMatrix(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==TILE_TYPE) {
			return createTileMatrix(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==HDIST_OLD_TYPE) {
//			return new PCRMatrixHDist_old(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==PROB_OLD_TYPE) {
//			return new PCRMatrixProb_old(length1_, length2_, delimiter_, hdistSum_);
		}
		throw new RuntimeException("Unknown PCRMatrix type "+type);
	}

	private static synchronized PCRMatrix createProbMatrix(int length1, int length2, int delimiter, boolean hdistSum){
		try{return PMPConstructor.newInstance(length1, length2, delimiter, hdistSum);}
		catch(Exception e){throw new RuntimeException("Failed to instantiate PCRMatrixProb", e);}
	}

	private static synchronized PCRMatrix createTileMatrix(int length1, int length2, int delimiter, boolean hdistSum){
		try{return PMTConstructor.newInstance(length1, length2, delimiter, hdistSum);}
		catch(Exception e){throw new RuntimeException("Failed to instantiate PCRMatrixTile", e);}
	}

	static synchronized boolean probLoaded() {return probLoaded;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses static configuration parameters for PCRMatrix classes.
	 * Handles parameters like localcounts, coding, matrixtype, and various flags.
	 *
	 * @param arg Full argument string
	 * @param a Parameter name
	 * @param b Parameter value
	 * @return true if parameter was recognized and parsed
	 */
	public static boolean parseStatic(String arg, String a, String b){
		
		if(a.equals("localcounts")){
			localCounts=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("coding")){
			if(Tools.startsWithDigit(b)) {
				DemuxData.DEFAULT_CODING=Integer.parseInt(b);
			}else if("raw".equalsIgnoreCase(b)) {
				DemuxData.DEFAULT_CODING=Sketch.RAW;
			}else if("a48".equalsIgnoreCase(b)) {
				DemuxData.DEFAULT_CODING=Sketch.A48;
			}else {
				assert(false) : "Unknown coding format: "+b;
			}
		}else if(a.equalsIgnoreCase("raw")){
			DemuxData.DEFAULT_CODING=Sketch.RAW;
		}else if(a.equalsIgnoreCase("a48")){
			DemuxData.DEFAULT_CODING=Sketch.A48;
		}else if(a.equals("writesent")){
			if(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")) {
				DemuxClient.writeSent="sent.txt";
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")) {
				DemuxClient.writeSent=null;
			}else {
				DemuxClient.writeSent=b;
			}
		}else if(a.equals("deltacounts")){
			DemuxData.deltaCountsDefault=Parse.parseBoolean(b);
		}else if(a.equals("deltacodes") || a.equals("deltabarcodes")){
			DemuxData.deltaBarcodesDefault=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("delta")){
			DemuxData.deltaBarcodesDefault=DemuxData.deltaCountsDefault=Parse.parseBoolean(b);
		}
		
		
		else if(a.equalsIgnoreCase("polya") || a.equalsIgnoreCase("addpolya")){
			addPolyA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyc") || a.equalsIgnoreCase("addpolyc")){
			addPolyC=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyg") || a.equalsIgnoreCase("addpolyg")){
			addPolyG=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyt") || a.equalsIgnoreCase("addpolyt")){
			addPolyT=Parse.parseBoolean(b);
		}else if(a.equals("hdistsum") || a.equals("pairhdist") || a.equals("hdistpair") || a.equals("sumhdist")){
			hdistSum0=Parse.parseBoolean(b);
		}else if(a.equals("matrixtype") || a.equals("type") || a.equals("mode") || a.equals("pcrmatrixtype")){
			if("hdist".equals(b)) {matrixType0=HDIST_TYPE;}
			else if("prob".equals(b) || "probability".equals(b)) {matrixType0=PROB_TYPE;}
			else if("tile".equals(b) || "bytile".equals(b)) {matrixType0=TILE_TYPE;}
			else {matrixType0=Integer.parseInt(b);}
		}else if(a.equals("probability") || a.equals("prob")){
			matrixType0=(Parse.parseBoolean(b) ? PROB_TYPE : (matrixType0==PROB_TYPE ? HDIST_TYPE : matrixType0));
		}else if(a.equals("bytile") || a.equals("tile")){
			matrixType0=(Parse.parseBoolean(b) ? TILE_TYPE : (matrixType0==TILE_TYPE ? PROB_TYPE : matrixType0));
		}else if(a.equalsIgnoreCase("matrixThreads")){
			matrixThreads=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("devmode")){
			devMode=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("ensuresorted")){
			DemuxData.ENSURE_SORTED=Parse.parseBoolean(b);
		}else if(PCRMatrixHDist.parseStatic(arg, a, b)) {
			callParseStatic(arg, a, b);//In case of shared flags with different defaults
		}else if(callParseStatic(arg, a, b)) {
			PCRMatrixHDist.parseStatic(arg, a, b);
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses instance-specific configuration parameters.
	 * Implementation varies by concrete PCRMatrix subclass.
	 *
	 * @param arg Full argument string
	 * @param a Parameter name
	 * @param b Parameter value
	 * @return true if parameter was recognized and parsed
	 */
	public abstract boolean parse(String arg, String a, String b);
	
	/** Finalizes static configuration after all parameters have been parsed.
	 * Calls postParse methods on subclasses and sets tile mode flag. */
	public static void postParseStatic(){
		PCRMatrixHDist.postParseStatic();
		callPostParseStatic();
		if(matrixType0==TILE_TYPE) {byTile=true;}
		else {
			assert(byTile==false);
			byTile=false;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           File I/O           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Represents a mapping between observed and assigned barcodes with count and distance data.
	 * Used for sorting and outputting barcode assignment results. */
	public static class MapLine implements Comparable<MapLine> {
		
		/**
		 * Constructs a MapLine with observed/assigned barcode mapping data.
		 *
		 * @param observed_ The observed barcode sequence
		 * @param assigned_ The assigned reference barcode
		 * @param count_ Number of observations
		 * @param prob_ Assignment probability (-1 if not calculated)
		 */
		public MapLine(String observed_, String assigned_, long count_, float prob_){
			observed=observed_;
			assigned=assigned_;
			count=count_;
			hdist=Barcode.hdist(observed, assigned);
			prob=prob_;
		}
		
		@Override
		public int compareTo(MapLine b) {
			int x=assigned.compareTo(b.assigned);
			if(x!=0) {return x;}
			if(count!=b.count) {return count>b.count ? -1 : 1;}
			if(hdist!=b.hdist) {return hdist-b.hdist;}
			return observed.compareTo(b.observed);
		}
		
		/** Converts MapLine to tab-separated byte representation.
		 * @return ByteBuilder containing formatted output */
		public ByteBuilder toBytes() {
			ByteBuilder bb=new ByteBuilder();
			bb.append(observed).tab().append(assigned);
			bb.tab().append(count).tab().append(hdist);
//			if(prob>=0) {bb.tab().append(prob, 5);}
			return bb;
		}
		
		@Override
		public String toString() {
			return toBytes().toString();
		}
		
		/** The observed barcode sequence */
		final String observed;
		/** The assigned reference barcode sequence */
		final String assigned;
		/** Number of times this mapping was observed */
		final long count;
		/** Hamming distance between observed and assigned sequences */
		final int hdist;
		/** Assignment probability (-1 if not calculated) */
		float prob=-1;
	}
	
	/**
	 * Prints barcode assignment map to file using barcode collection for counts.
	 *
	 * @param assignmentMap Map from observed to assigned barcodes
	 * @param mapOut Output file path
	 * @param counts Barcode collection with count data
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 */
	public final void printAssignmentMap(HashMap<String, String> assignmentMap,
			String mapOut, Collection<Barcode> counts, boolean overwrite, boolean append) {
		readlock();
		HashMap<String, Barcode> countMap=Barcode.barcodesToMap(counts);
		printAssignmentMap(assignmentMap, mapOut, countMap, overwrite, append);
		readunlock();
	}
	
	/**
	 * Prints barcode assignment map to file using barcode map for counts.
	 *
	 * @param assignmentMap Map from observed to assigned barcodes
	 * @param mapOut Output file path
	 * @param counts Barcode map with count data
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 */
	public void printAssignmentMap(HashMap<String, String> assignmentMap,
			String mapOut, HashMap<String, Barcode> counts, boolean overwrite, boolean append) {
		readlock();
		ArrayList<MapLine> lines=new ArrayList<MapLine>();
		for(Entry<String, String> e : assignmentMap.entrySet()) {
			String a=e.getKey(), b=e.getValue();
			Barcode v=(counts==null ? null : counts.get(a));
			lines.add(new MapLine(a, b, v==null ? 0 : v.count(), -1));
		}
		Collections.sort(lines);
		printAssignmentMap(lines, mapOut, overwrite, append);
		readunlock();
	}
	
	/**
	 * Static version of printAssignmentMap using barcode collection.
	 *
	 * @param assignmentMap Map from observed to assigned barcodes
	 * @param mapOut Output file path
	 * @param counts Barcode collection with count data
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 */
	public static final void printAssignmentMapStatic(HashMap<String, String> assignmentMap,
			String mapOut, Collection<Barcode> counts, boolean overwrite, boolean append) {
		HashMap<String, Barcode> countMap=Barcode.barcodesToMap(counts);
		printAssignmentMapStatic(assignmentMap, mapOut, countMap, overwrite, append);
	}
	
	/**
	 * Static version of printAssignmentMap using barcode map.
	 *
	 * @param assignmentMap Map from observed to assigned barcodes
	 * @param mapOut Output file path
	 * @param counts Barcode map with count data
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 */
	public static void printAssignmentMapStatic(HashMap<String, String> assignmentMap,
			String mapOut, HashMap<String, Barcode> counts, boolean overwrite, boolean append) {
		ArrayList<MapLine> lines=new ArrayList<MapLine>();
		for(Entry<String, String> e : assignmentMap.entrySet()) {
			String a=e.getKey(), b=e.getValue();
			Barcode v=(counts==null ? null : counts.get(a));
			lines.add(new MapLine(a, b, v==null ? 0 : v.count(), -1));
		}
		Collections.sort(lines);
		printAssignmentMap(lines, mapOut, overwrite, append);
	}
	
	/**
	 * Prints pre-formatted assignment map lines to file.
	 *
	 * @param lines List of MapLine objects to output
	 * @param mapOut Output file path
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 */
	public static final void printAssignmentMap(ArrayList<MapLine> lines, String mapOut, 
			boolean overwrite, boolean append) {
		ByteStreamWriter bsw=new ByteStreamWriter(mapOut, overwrite, append, true);
		bsw.start();
		for(MapLine line : lines) {bsw.println(line.toBytes());}
		bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/
	/*----------------            HDist             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Fills IntList with hamming distances between query and reference sequences.
	 * @param q Query sequence
	 * @param list Array of reference sequences
	 * @param hdList List to populate with distances
	 */
	protected static final void fillHDistList(byte[] q, byte[][] list, IntList hdList){
		hdList.clear();
		for(byte[] r : list) {
			hdList.add(hdist(q, r));
		}
	}

	/**
	 * Calculates hamming distance between two byte sequences.
	 * @param q First sequence
	 * @param r Second sequence
	 * @return Hamming distance between sequences
	 */
	public static final int hdist(byte[] q, byte[] r) {return Barcode.hdist(q, r);}
	
	/**
	 * Combines two distance values based on hdistSum setting.
	 * @param d1 First distance value
	 * @param d2 Second distance value
	 * @return Combined distance (sum or maximum)
	 */
	public final int hdist(int d1, int d2) {return hdistSum ? d1+d2 : Tools.max(d1, d2);}
	
	/**
	 * Calculates distance between barcode and string.
	 * @param q Query barcode
	 * @param r Reference string
	 * @return Distance between sequences
	 */
	public final int hdist(Barcode q, String r) {
		return hdist(q.name, r);
	}
	
	/**
	 * Calculates distance between two strings based on hdistSum setting.
	 * For dual barcodes, computes left and right distances separately.
	 *
	 * @param q Query string
	 * @param r Reference string
	 * @return Combined distance
	 */
	public final int hdist(String q, String r) {
		return hdistSum ? Barcode.hdist(q, r) : Tools.max(Barcode.hdistL(q, r), Barcode.hdistR(q, r));
	}
	
	/**
	 * Encodes index and two distance values into a single long.
	 * Uses bit packing for efficient storage and retrieval.
	 *
	 * @param idx Index value
	 * @param hdist First distance
	 * @param hdist2 Second distance
	 * @return Packed long containing all three values
	 */
	public static final long encodeHDist(long idx, long hdist, long hdist2) {
		return idx|(hdist<<32)|(hdist2<<48);
	}
	
	/**
	 * Filters barcodes to only include those meeting minimum count threshold.
	 * @param codeCounts Input barcode collection
	 * @param minCount Minimum count threshold
	 * @return Filtered list of barcodes above threshold
	 */
	protected final ArrayList<Barcode> highpass(Collection<Barcode> codeCounts, long minCount) {
		if(minCount<=1 || codeCounts.size()<minSizeToFilter) {
			return codeCounts instanceof ArrayList ? (ArrayList<Barcode>)codeCounts 
					: new ArrayList<Barcode>(codeCounts);
		}
		ArrayList<Barcode> list=new ArrayList<Barcode>(256);
		for(Barcode b : codeCounts) {
			if(b.count()>=minCount) {list.add(b);}
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------            HDist             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Finds closest matching barcode using hamming distance.
	 * Dispatches to single or dual barcode implementation.
	 *
	 * @param s Query sequence
	 * @param maxHDist Maximum allowed distance
	 * @param clearzone Minimum distance margin to second-best match
	 * @return Closest matching barcode or null if none found
	 */
	protected final Barcode findClosestHDist(String s, int maxHDist, int clearzone) {
		return length2<1 ? findClosestSingleHDist(s, maxHDist, clearzone) : 
			findClosestDualHDist(s, maxHDist, clearzone);
	}
	
	/**
	 * Finds closest matching single barcode using hamming distance.
	 *
	 * @param q Query sequence
	 * @param maxHDist Maximum allowed distance
	 * @param clearzone Minimum distance margin to second-best match
	 * @return Closest matching barcode or null if none found
	 */
	protected final Barcode findClosestSingleHDist(String q, int maxHDist, int clearzone) {
		assert(q.length()==length1) : q+", "+length1+", "+length2;
		final byte[] left=q.getBytes();
		long packed=findClosestHDist(left, leftBytes, maxHDist, clearzone);
		int idx=(int)(packed&0xFFFFFFFFL);
		int hdist=(int)((packed>>32)&0xFFFFL);
		int hdist2=(int)((packed>>48)&0xFFFFL);
//		System.err.println("q="+q+/*", packed="+packed+*/", idx="+idx+", hd="+hdist+", hd2="+hdist2);
		assert(idx<leftCodes.length);
		assert(idx<0 || hdist+clearzone<=hdist2);
		assert(idx<0 || hdist<=maxHDist) : idx+", "+hdist+", "+maxHDist;
		return idx<0 ? null : leftCodes[idx];
	}
	
	/**
	 * Finds closest matching dual barcode using hamming distance.
	 * Splits query into left and right segments and finds best matches for each.
	 *
	 * @param q Query sequence
	 * @param maxHDist Maximum allowed distance
	 * @param clearzone Minimum distance margin to second-best match
	 * @return Closest matching barcode or null if none found
	 */
	protected final Barcode findClosestDualHDist(String q, int maxHDist, int clearzone) {
		//if(verbose) {System.err.println("Looking for "+q);}
		byte[] left=new byte[length1];
		byte[] right=new byte[length2];
		for(int i=0; i<length1; i++) {left[i]=(byte) q.charAt(i);}
		for(int i=length2-1, j=q.length()-1; i>=0; i--, j--) {
			right[i]=(byte) q.charAt(j);
		}
		final long lpacked=findClosestHDist(left, leftBytes, maxHDist, clearzone);
		final long rpacked=findClosestHDist(right, rightBytes, maxHDist, clearzone);
		
		final int lidx=(int)(lpacked&0xFFFFFFFFL);
		final int lhdist=(int)((lpacked>>32)&0xFFFFL);
		final int lhdist2=(int)((lpacked>>48)&0xFFFFL);
		final int lmargin=lhdist2-lhdist;
		assert(lidx<leftCodes.length);
		assert(lidx<0 || lmargin>=clearzone);
		assert(lidx<0 || lhdist<=maxHDist) : lidx+", "+lhdist+", "+lhdist2+", "+lmargin+", "+maxHDist;
		
		final int ridx=(int)(rpacked&0xFFFFFFFFL);
		final int rhdist=(int)((rpacked>>32)&0xFFFFL);
		final int rhdist2=(int)((rpacked>>48)&0xFFFFL);
		final int rmargin=rhdist2-rhdist;
		assert(ridx<rightCodes.length);
		assert(ridx<0 || rmargin>=clearzone) : ridx+", "+rhdist+", "+rhdist2+", "+rmargin;
		assert(ridx<0 || rhdist<=maxHDist) : ridx+", "+rhdist+", "+rhdist2+", "+rmargin;
		
		if(lidx<0 || ridx<0) {return null;}
		
		if(hdistSum) {
//			int idx=lidx*rightBytes.length+ridx;
//			Barcode bc=allCodes[idx];
			int hdist=lhdist+rhdist;
			int hdist2=lhdist2+rhdist2;
			int margin=hdist2-hdist;
//			System.err.println("\n"+q+"\n"+bc+"\ntrue="+bc.hdist(q)+", max="+maxHDist+", margin="+margin+", h="+hdist+", h2="+hdist2+"\n");
//			System.err.println((hdist>maxHDist)+" "+(margin>=clearzone)+" "+clearzone);
			if(hdist>maxHDist || margin<clearzone) {return null;}
//			System.err.println("PASS");
		}
		
		int idx=lidx*rightBytes.length+ridx;
		Barcode bc=allCodes[idx];
		//if(verbose) {System.err.println(q+" -> "+bc.name+" ("+bc.expected+")");}
		if(bc.expected==1) {assert(expectedMap.containsKey(bc.name)) : "\n"+expectedMap;}
		else {assert(!expectedMap.containsKey(bc.name));}
//		assert(!hdistSum || bc.hdist(q)<=maxHDist) : "\n"+q+"\n"+bc+"\n"+bc.hdist(q)+", "+maxHDist+"\n";
		assert(hdist(bc.name, q)<=maxHDist) : "\n"+q+"\n"+bc+"\n"+hdist(bc.name, q)+", "+maxHDist+"\n";
		return bc;
	}
	
	/**
	 * Finds closest matching sequence from a list using hamming distance.
	 * Returns packed long containing index and distance information.
	 *
	 * @param q Query sequence
	 * @param list Array of reference sequences
	 * @param maxHDist Maximum allowed distance
	 * @param clearzone Minimum distance margin to second-best match
	 * @return Packed long with index and distances, or -1 if none found
	 */
	protected final long findClosestHDist(byte[] q, byte[][] list, int maxHDist, int clearzone) {
		assert(!expectedList.isEmpty());
		int best=-1, best2=-1;
		
		int hdist=Tools.min(q.length, letters);
		int hdist2=hdist;
		
		for(int i=0; i<list.length; i++) {
			byte[] ref=list[i];
			final int d=hdist(q, ref);
			if(d<hdist2) {
				best2=i;
				hdist2=d;
				if(d<hdist) {
					best2=best;
					hdist2=hdist;
					best=i;
					hdist=d;
				}
				if(hdist2<clearzone) {return -1;}
			}
		}
		//if(verbose) {System.err.println("best="+best+", hdist="+hdist+", hdist2="+hdist2);}
		if(best<0 || hdist>maxHDist) {return -1;}
		if(hdist+clearzone>hdist2) {return -1;}
		//if(verbose) {System.err.println("q="+new String(q)+" -> best="+new String(list[best]));}
		return encodeHDist(best, hdist, hdist2);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Various            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Finds closest matching barcode for the given sequence.
	 * Implementation varies by concrete subclass.
	 * @param s Query sequence
	 * @return Closest matching barcode or null
	 */
	public abstract Barcode findClosest(String s);
	
	/**
	 * Refines barcode assignments based on observed counts.
	 * Implementation varies by concrete subclass.
	 * @param codeCounts Observed barcode counts
	 * @param minCount Minimum count threshold
	 */
	public abstract void refine(Collection<Barcode> codeCounts, long minCount);
	
	/**
	 * Creates mapping from observed to assigned barcodes.
	 * Implementation varies by concrete subclass.
	 *
	 * @param codeCounts Observed barcode counts
	 * @param minCount Minimum count threshold
	 * @return Map from observed to assigned barcode sequences
	 */
	public abstract HashMap<String, String> makeAssignmentMap(Collection<Barcode> codeCounts, long minCount);
	
	/**
	 * Populates count data for barcode list.
	 * Implementation varies by concrete subclass.
	 * @param list List of barcodes to populate
	 * @param minCount Minimum count threshold
	 */
	public abstract void populateCounts(ArrayList<Barcode> list, long minCount);
	
	/** Calculates probability matrices from count data.
	 * Implementation varies by concrete subclass. */
	public abstract void makeProbs();

	/** Initializes data structures for the matrix.
	 * Implementation varies by concrete subclass. */
	public abstract void initializeData();
	
	/**
	 * Validates licensing or other requirements.
	 * Implementation varies by concrete subclass.
	 * @return true if matrix is valid for use
	 */
	protected abstract boolean valid();// {return true;}

	/*--------------------------------------------------------------*/
	/*----------------          Populating          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Populates the matrix with expected barcode sequences.
	 * Creates internal data structures for tracking expected barcodes.
	 * @param expected Collection of expected barcode sequences
	 */
	@SuppressWarnings("unchecked")
	public final void populateExpected(Collection<String> expected) {
		writelock();
		assert(expectedList==null);
		
		expectedList=new ArrayList<Barcode>(expected.size());
		expectedMap=new HashMap<String, Barcode>(expectedList.size()*2);
		for(String s : expected) {
			assert(s.length()==counts.length) : "Expected barcode lengths do not match actual barcode lengths."
				+"\n"+s+", "+s.length()+", "+counts.length+", "+expected.size();
//			assert(!expectedMap.containsKey(s)) : "Duplicate key: "+s;
			if(!expectedMap.containsKey(s)) {
				Barcode bc=new Barcode(s);
				expectedList.add(bc);
				expectedMap.put(s, bc);
			}
		}
//		System.err.println("size="+expectedList.size());
		writeunlock();
	}
	/**
	 * Populates the matrix with expected barcode objects.
	 * Creates internal data structures from existing Barcode objects.
	 * @param expected Collection of expected Barcode objects
	 */
	public final void populateExpectedFromBarcodes(Collection<Barcode> expected) {
		writelock();
		assert(expectedList==null);
		assert(expected!=null && expected.size()>0) : expected;
		expectedList=new ArrayList<Barcode>(expected.size());
		expectedMap=new HashMap<String, Barcode>(expectedList.size()*2);
		for(Barcode bc0 : expected) {
			String s=bc0.name;
			Barcode bc=new Barcode(s);
			assert(s.length()==counts.length);
//			assert(!expectedMap.containsKey(s)) : "Duplicate key: "+s;
			if(!expectedMap.containsKey(s)) {
				expectedList.add(bc);
				expectedMap.put(s, bc);
			}
		}
//		System.err.println("size="+expectedList.size());
		writeunlock();
	}
	
	/** Populates data for unexpected (novel) barcodes.
	 * Implementation varies by concrete subclass. */
	public abstract void populateUnexpected();
	
	/** Populates split barcode data using default DemuxData settings.
	 * Creates separate arrays for left and right barcode segments. */
	public final void populateSplitCodes() {
		DemuxData dd=new DemuxData(length1, length2, delimiter);
		populateSplitCodes(dd);
	}
	
	/**
	 * Populates split barcode data using provided DemuxData configuration.
	 * Creates arrays for left/right segments and all possible combinations.
	 * @param dd DemuxData configuration object
	 */
	public final void populateSplitCodes(DemuxData dd) {
		writelock();
		assert(leftBytes==null) : "Already populated.";
		LinkedHashSet<String> set1=new LinkedHashSet<String>();
		LinkedHashSet<String> set2=new LinkedHashSet<String>();
		@SuppressWarnings("unchecked")
		LinkedHashSet<String>[] sets=new LinkedHashSet[] {set1, set2};
		for(Barcode b : expectedList) {
			String code1=b.name.substring(0, length1);
			set1.add(code1);
			if(length2>0) {
				String code2=b.name.substring(start2);
				assert(code2.length()==length2);
				set2.add(code2);
			}
		}

		if(dd.addPolyA && length1>0) {set1.add(poly('A', length1));}
		if(dd.addPolyC && length1>0) {set1.add(poly('C', length1));}
		if(dd.addPolyG && length1>0) {set1.add(poly('G', length1));}
		if(dd.addPolyT && length1>0) {set1.add(poly('T', length1));}
		
		if(dd.addPolyA && length2>0) {set2.add(poly('A', length2));}
		if(dd.addPolyC && length2>0) {set2.add(poly('C', length2));}
		if(dd.addPolyG && length2>0) {set2.add(poly('G', length2));}
		if(dd.addPolyT && length2>0) {set2.add(poly('T', length2));}
		
		final int size1=set1.size();
		final int size2=set2.size();
		leftBytes=new byte[size1][length1];
		rightBytes=new byte[size2][length2];
		splitBytes=new byte[][][] {leftBytes, rightBytes};
		
		leftCodes=new Barcode[size1];
		rightCodes=new Barcode[size2];
		allCodes=new Barcode[size1*size2];
		allCodesMap=new HashMap<String, Barcode>((size1*size2*3)/2);
		
		for(int i=0; i<sets.length; i++){
			LinkedHashSet<String> set=sets[i];
			int j=0;
			for(String s : set) {
				Tools.copy(s, splitBytes[i][j]);
				j++;
			}
		}
		
		ByteBuilder bb=new ByteBuilder();
		if(length2>0) {//dual
			for(int i=0; i<leftBytes.length; i++) {
				bb.clear().append(leftBytes[i]);
				if(delimiter>0) {bb.append((byte)delimiter);}
				int len=bb.length();
				leftCodes[i]=new Barcode(new String(leftBytes[i]));
				for(int j=0; j<rightBytes.length; j++) {
					bb.setLength(len);
					bb.append(rightBytes[j]);
					String s=bb.toString();
					Barcode bc=expectedMap.get(s);
					if(bc==null) {
						bc=new Barcode(s, 0, 0);
					}else {
						assert(bc.expected==1);
					}
					int idx=i*size2+j;
					assert(allCodes[idx]==null);
					allCodes[idx]=bc;
					allCodesMap.put(s, bc);
					rightCodes[j]=new Barcode(new String(rightBytes[j]));
				}
			}
		}else {//single
//			assert(leftCodes.length==expectedList.size()) : leftCodes.length+", "+expectedList.size()+
//			"\n"+expectedList+"\n"+Arrays.toString(leftCodes);
			assert(leftCodes.length==leftBytes.length) : leftCodes.length+", "+expectedList.size()+
				"\n"+expectedList+"\n"+Arrays.toString(leftCodes);
			
			for(int i=0; i<leftBytes.length; i++) {
				String s=new String(leftBytes[i]);
				Barcode bc=expectedMap.get(s);
				if(bc==null) {
					bc=new Barcode(s, 0, 0);
				}else {
					assert(bc.expected==1);
				}
				leftCodes[i]=bc;
				allCodesMap.put(bc.name, bc);
			}
			allCodes=leftCodes;
		}
		writeunlock();
	}
	
	/** Resets all count data to zero.
	 * Clears count matrix and barcode counts. */
	public void clearCounts() {
		writelock();
		Tools.fill(counts, 0);
		totalCounted=totalAssigned=totalAssignedToExpected=0;
		for(Barcode bc : expectedList) {bc.setCount(0);}
		writeunlock();
	}
	
	/** Calculates fraction of reads that were assigned to any barcode.
	 * @return Fraction of total reads that were assigned */
	public final float assignedFraction() {
		return (totalAssigned/(1.0f*totalCounted));
	}
	
	/** Calculates fraction of reads assigned to expected barcodes.
	 * @return Fraction of total reads assigned to expected barcodes */
	public final float expectedFraction() {
		return (totalAssignedToExpected/(1.0f*totalCounted));
	}
	
	/** Calculates fraction of reads assigned to unexpected (chimeric) barcodes.
	 * @return Fraction of reads assigned to chimeric barcodes */
	public final float chimericFraction() {
		return ((totalAssigned-totalAssignedToExpected)/(1.0f*totalCounted));
	}

	/**
	 * Adds observed barcode data with reference barcode match.
	 * @param query Observed barcode
	 * @param ref Reference barcode match
	 */
	public final void add(Barcode query, Barcode ref) {add(query.name, ref, query.count());}
	/**
	 * Adds observed barcode data with specified count.
	 * @param query Observed barcode sequence
	 * @param ref Reference barcode match
	 * @param count Number of observations
	 */
	public final void add(String query, Barcode ref, long count) {add(query, ref, count, 0);}
	/**
	 * Adds observed vs reference base counts to the matrix at each position.
	 * Updates total counts and barcode-specific counts.
	 *
	 * @param query Observed barcode sequence
	 * @param ref Reference barcode (null for unassigned)
	 * @param count Number of observations
	 * @param pos Starting position offset
	 */
	public final void add(String query, Barcode ref, long count, int pos) {
		assert(ref==null || ref.length()==counts.length || matrixType0>=5);
		for(int i=0; i<query.length(); i++, pos++) {
			final int q=query.charAt(i), r=(ref==null ? 'N' : ref.charAt(i));
			final byte xq=baseToNumber[q], xr=baseToNumber[r];
			counts[pos][xq][xr]+=count;
		}
		totalCounted+=count;
		if(ref!=null) {
			ref.incrementSync(count);
			totalAssigned+=count;
			totalAssignedToExpected+=ref.expected*count;
		}
	}
	
	/**
	 * Converts count matrix to tab-separated byte format.
	 * @param bb ByteBuilder to append to (created if null)
	 * @return ByteBuilder containing formatted matrix data
	 */
	public ByteBuilder toBytes(ByteBuilder bb) {
		return toBytes(counts, bb);
	}
	
	/**
	 * Converts count matrix to accuracy statistics format.
	 * @param bb ByteBuilder to append to (created if null)
	 * @param excludeUnknown Whether to exclude unknown base calls
	 * @return ByteBuilder containing accuracy data
	 */
	public ByteBuilder toAccuracy(ByteBuilder bb, boolean excludeUnknown) {
		return toAccuracy(counts, bb, excludeUnknown);
	}
	
	/** Returns whether this matrix processes data by sequencing tile.
	 * @return false for base implementation */
	public boolean byTile() {return false;}
	/**
	 * Converts probability data to byte format.
	 * Implementation varies by concrete subclass.
	 * @param bb ByteBuilder to append to
	 * @return ByteBuilder containing probability data
	 */
	public abstract ByteBuilder toBytesProb(ByteBuilder bb);
	
	/**
	 * Adds another PCRMatrix's data to this matrix.
	 * Combines count arrays and totals from both matrices.
	 * @param p PCRMatrix to add data from
	 */
	final public void add(PCRMatrix p) {//Unused?
		p.readlock();
		writelock();
		Tools.add(counts, p.counts);
		totalCounted+=p.totalCounted;
		totalAssigned+=p.totalAssigned;
		totalAssignedToExpected+=p.totalAssignedToExpected;
		writeunlock();
		p.readunlock();
	}
	
	/**
	 * Retrieves expected barcode by sequence string.
	 * @param s Barcode sequence
	 * @return Barcode object or null if not found
	 */
	protected final Barcode getBarcode(String s) {return expectedMap.get(s);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Static method to convert count array to tab-separated format.
	 * Creates detailed position-by-base count table.
	 *
	 * @param counts 3D count array [position][call][reference]
	 * @param bb ByteBuilder to append to (created if null)
	 * @return ByteBuilder containing formatted count data
	 */
	public static final ByteBuilder toBytes(long[][][] counts, ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		bb.append("#pos\tcall\tA\tC\tG\tT\tN\tSum\n");
		for(int pos=0; pos<counts.length; pos++) {
			for(int xq=0; xq<5; xq++) {
				final byte q=numberToBase[xq];
				bb.append(pos).tab().append(q);
				long sum=0;
				for(int xr=0; xr<5; xr++) {
//					final byte r=numberToBase[xr];
					final long count=counts[pos][xq][xr];
					sum+=count;
					bb.tab().append(count);
				}
				bb.tab().append(sum).nl();
			}
		}
		return bb;
	}
	
	/**
	 * Static method to convert count array to accuracy statistics.
	 * Calculates call and reference accuracy by position.
	 *
	 * @param counts 3D count array [position][call][reference]
	 * @param bb ByteBuilder to append to (created if null)
	 * @param excludeUnknown Whether to exclude unknown bases from accuracy calculation
	 * @return ByteBuilder containing accuracy statistics
	 */
	public static final ByteBuilder toAccuracy(long[][][] counts, ByteBuilder bb, boolean excludeUnknown) {
		if(bb==null) {bb=new ByteBuilder();}

		bb.append("#pos\tA_call\tC_call\tG_call\tT_call\tA_ref\tC_ref\tG_ref\tT_ref\n");
		for(int pos=0; pos<counts.length; pos++) {
			long[] refCount=new long[5];
			long[] callCount=new long[5];
			long[] unknown=new long[5];
			long[] callCorrect=new long[5];
			long[] refCorrect=new long[5];
			long sum=0;
			for(int xc=0; xc<5; xc++) {
				for(int xr=0; xr<5; xr++) {
					final long count=counts[pos][xc][xr];
					refCount[xr]+=count;
					callCount[xc]+=count;
					sum+=count;
					if(xc==xr) {
						refCorrect[xr]+=count;
						callCorrect[xc]+=count;
					}
					if(xr==4) {unknown[xc]=count;}
				}
			}
			if(sum>refCount[4]) {
				bb.append(pos);
				for(int i=0; i<4; i++) {
					long count=callCount[i], correct=callCorrect[i];
					if(excludeUnknown) {count-=unknown[i];}
					double accuracy=(count<1 ? 1 : correct/(double)count);
					bb.tab().append(accuracy, 4);
//					assert(false) : "\ncount="+count+", correct="+correct
//						+", cc4="+callCount[4]+", acc="+accuracy
//						+"\ncounts="+Arrays.toString(callCount)
//						+"\ncorrects="+Arrays.toString(callCorrect);
				}
				for(int i=0; i<4; i++) {
					long count=refCount[i], correct=refCorrect[i];
					double accuracy=(count<1 ? 1 : correct/(double)count);
					bb.tab().append(accuracy, 4);
				}
				bb.nl();
			}
		}
		return bb;
	}
	
	/**
	 * Creates a string of repeated characters.
	 * @param c Character to repeat
	 * @param len Number of repetitions
	 * @return String containing repeated character
	 */
	private static String poly(char c, int len) {
		ByteBuilder bb=new ByteBuilder(len);
		for(int i=0; i<len; i++) {bb.append(c);}
		return bb.toString();
	}

	/** Acquires write lock for thread-safe matrix modification. */
	protected final void writelock() {rwlock.writeLock().lock();}
	/** Releases write lock after matrix modification. */
	protected final void writeunlock() {rwlock.writeLock().unlock();}
	/** Acquires read lock for thread-safe matrix access. */
	protected final void readlock() {rwlock.readLock().lock();}
	/** Releases read lock after matrix access. */
	protected final void readunlock() {rwlock.readLock().unlock();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the ReadWriteLock used for thread synchronization.
	 * @return ReadWriteLock instance for this matrix */
	public final ReadWriteLock rwlock() {return rwlock;}
	/** ReadWriteLock for thread-safe access to matrix data */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/** Total length of barcode including delimiter */
	protected final int length;
	/** Length of first barcode segment */
	protected final int length1;
	/** Length of second barcode segment (0 for single barcodes) */
	protected final int length2;
	/** Delimiter character between barcode segments (0 for none) */
	protected final int delimiter;
	/** Position of delimiter in barcode (-1 if none) */
	protected final int delimiterPos;
	/** Starting position of second barcode segment (-1 if single) */
	protected final int start2;
	/** Total number of letter positions (excluding delimiter) */
	protected final int letters;
	/** Whether to sum hamming distances rather than take maximum */
	protected final boolean hdistSum;
	
	/** [position][Call: A,C,G,T,N][Ref: A,C,G,T,N,Sum] */
	protected final long[][][] counts;
	/** Total number of barcodes processed */
	protected long totalCounted;
	/** Total number of barcodes assigned to any reference */
	protected long totalAssigned;
	/** Total number of barcodes assigned to expected references */
	protected long totalAssignedToExpected;
	
	/** List of expected barcode sequences */
	protected ArrayList<Barcode> expectedList;
	/** Map from barcode sequence to expected Barcode object */
	protected HashMap<String, Barcode> expectedMap;

	/** Byte arrays for left barcode segments */
	protected byte[][] leftBytes;
	/** Byte arrays for right barcode segments */
	protected byte[][] rightBytes;
	/** Combined array containing leftBytes and rightBytes */
	protected byte[][][] splitBytes;

	/** Barcode objects for left segments */
	protected Barcode[] leftCodes;
	/** Barcode objects for right segments */
	protected Barcode[] rightCodes;
	/** Array of all possible barcode combinations */
	protected Barcode[] allCodes;
	/** Map from barcode sequence to Barcode object for all combinations */
	protected HashMap<String, Barcode> allCodesMap;

	/** Whether to output verbose debugging information */
	public boolean verbose=true;
	/** Global verbose flag for additional debugging output */
	public static boolean verbose2=false;
	/** Whether this matrix instance is in an error state */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global error state flag */
	public static boolean errorStateS=false;
	
	public static int minSizeToFilter=100000;

	/** Whether to add poly-A sequences to barcode sets */
	protected static boolean addPolyA=false;
	/** Whether to add poly-C sequences to barcode sets */
	protected static boolean addPolyC=false;
	/** Whether to add poly-G sequences to barcode sets */
	protected static boolean addPolyG=false;
	/** Whether to add poly-T sequences to barcode sets */
	protected static boolean addPolyT=false;
	
	//This is just a maximum; will not go over Shared.threads()
	/** Maximum threads to use for matrix operations (capped by Shared.threads()) */
	protected static int matrixThreads=64;
	
	/** Whether to use local count storage */
	protected static boolean localCounts=true;
	/** Whether to process data by sequencing tile */
	public static boolean byTile=false;
	/** Whether development mode is enabled */
	protected static boolean devMode=false;
	
	/** Default setting for hamming distance calculation mode */
	protected static boolean hdistSum0=false;
	/** Constant for legacy hamming distance matrix type */
	/** Constant for tile-based matrix type */
	/** Constant for probability matrix type */
	/** Constant for hamming distance matrix type */
	public static final int HDIST_TYPE=0, PROB_TYPE=1, TILE_TYPE=2, 
			HDIST_OLD_TYPE=7, PROB_OLD_TYPE=8;
	/** Default matrix type for new instances */
	protected static int matrixType0=HDIST_TYPE;
	
	/** Lookup array for converting bases to numeric codes */
	protected static final byte[] baseToNumber=AminoAcid.baseToNumber4;
	/** Lookup array for converting numeric codes to bases */
	protected static final byte[] numberToBase=AminoAcid.numberToBase;

	/*--------------------------------------------------------------*/
	/*----------------          Reflection          ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean probLoaded=false;

	private static Class<? extends PCRMatrix> pcrMatrixProbAbstractClass=getPMPAClass();
	private static Class<? extends PCRMatrix> pcrMatrixProbClass=getPMPClass();
	private static Class<? extends PCRMatrix> pcrMatrixTileClass=getPMTClass();
	private static Method PMPA_parseStatic=getParseStatic(pcrMatrixProbAbstractClass);
	private static Method PMPA_postParseStatic=getPostParseStatic(pcrMatrixProbAbstractClass);
	private static Constructor<? extends PCRMatrix> PMPConstructor=getConstructor(pcrMatrixProbClass);
	private static Constructor<? extends PCRMatrix> PMTConstructor=getConstructor(pcrMatrixTileClass);

	private static synchronized Class<? extends PCRMatrix> getPMPAClass(){
		try{
			Class<? extends PCRMatrix> c=(Class<? extends PCRMatrix>)Class.forName("barcode.prob.PCRMatrixProbAbstract");
			probLoaded=true;
			return c;
		}catch(ClassNotFoundException e) {return barcode.stub.PCRMatrixProbAbstract.class;}
	}

	private static Class<? extends PCRMatrix> getPMPClass(){
		try{return (Class<? extends PCRMatrix>)Class.forName("barcode.prob.PCRMatrixProb");}
		catch(ClassNotFoundException e) {return barcode.stub.PCRMatrixProb.class;}
	}

	private static Class<? extends PCRMatrix> getPMTClass(){
		try{return (Class<? extends PCRMatrix>)Class.forName("barcode.prob.PCRMatrixTile");}
		catch(ClassNotFoundException e) {return barcode.stub.PCRMatrixTile.class;}
	}

	private static boolean callParseStatic(String arg, String a, String b) {
		try{return (Boolean)PMPA_parseStatic.invoke(null, arg, a, b);}
		catch(Exception e){throw new RuntimeException("Error calling parseStatic via reflection", e);}
	}

	private static void callPostParseStatic() {
		try{PMPA_postParseStatic.invoke(null);}
		catch(Exception e){throw new RuntimeException("Error calling postParseStatic via reflection", e);}
	}

	private static Method getParseStatic(Class<?> c) {
		try{return c.getMethod("parseStatic", String.class, String.class, String.class);}
		catch(Exception e){throw new RuntimeException("Error calling parseStatic via reflection", e);}
	}

	private static Method getPostParseStatic(Class<?> c) {
		try{return c.getMethod("postParseStatic");}
		catch(Exception e){throw new RuntimeException("Error calling postParseStatic via reflection", e);}
	}

	private static Constructor<? extends PCRMatrix> getConstructor(Class<?> c) {
		try{return (Constructor<? extends PCRMatrix>)c.getConstructor(int.class, int.class, int.class, boolean.class);}
		catch(Exception e){throw new RuntimeException("Error calling getConstructor via reflection", e);}
	}

}
