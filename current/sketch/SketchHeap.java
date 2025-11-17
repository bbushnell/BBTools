package sketch;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import structures.LongHashMap;
import structures.LongHeap;
import structures.LongHeapMap;
import structures.LongHeapSet;
import structures.LongHeapSetInterface;

/**
 * A bounded heap-based data structure for maintaining k-mer sketches.
 * Supports both set-based and map-based storage for k-mer frequency tracking.
 * Used for MinHash sketching in sequence comparison and taxonomic classification.
 * @author Brian Bushnell
 */
public class SketchHeap {
	
	/**
	 * Constructs a SketchHeap with specified capacity and counting behavior.
	 * Automatically selects between set-based or map-based storage depending
	 * on whether frequency counting is required.
	 *
	 * @param limit Maximum number of k-mers to store in the heap
	 * @param minKeyOccuranceCount_ Minimum occurrence count for k-mer retention
	 * @param trackCounts Whether to track individual k-mer frequencies
	 */
	SketchHeap(int limit, int minKeyOccuranceCount_, boolean trackCounts){
		minKeyOccuranceCount=minKeyOccuranceCount_;
		setMode=minKeyOccuranceCount<2 && !trackCounts;
		if(setMode){
			setOrMap=set=new LongHeapSet(limit);
			map=null;
			heap=set.heap;
		}else{
			if(minKeyOccuranceCount>1){limit=(int)Tools.min(10000000, limit*SketchObject.sketchHeapFactor);}
			setOrMap=map=new LongHeapMap(limit);
			set=null;
			heap=map.heap;
		}
	}
	
	/**
	 * Resets all metadata fields and clears the k-mer data structure.
	 * Preserves filename if requested for reuse with different sequence data.
	 * @param clearFname Whether to also clear the filename field
	 */
	public void clear(boolean clearFname){
		taxID=-1;
		imgID=-1;
		genomeSizeBases=0;
		genomeSizeKmers=0;
		genomeSequences=0;
		if(baseCounts!=null){Arrays.fill(baseCounts, 0);}
		r16S=null;
		r18S=null;
		probSum=0;
		taxName=null;
		name0=null;
		if(clearFname){fname=null;}
		setOrMap.clear();
	}
	
	/**
	 * Merges another SketchHeap into this one.
	 * Combines k-mer data and merges metadata fields, preferring non-null values.
	 * @param b The SketchHeap to merge into this one
	 */
	public void add(SketchHeap b){
		if(taxID<0){taxID=b.taxID;}
		if(imgID<0){imgID=b.imgID;}
		if(taxName==null){taxName=b.taxName;}
		if(name0==null){name0=b.name0;}
		if(fname==null){fname=b.fname;}
		genomeSizeBases+=b.genomeSizeBases;
		genomeSizeKmers+=b.genomeSizeKmers;
		genomeSequences+=b.genomeSequences;
		if(baseCounts!=null){Tools.add(baseCounts, b.baseCounts);}
		set16S(b.r16S);
		set18S(b.r18S);
		probSum+=b.probSum;
		if(setMode){
			set.add(b.set);
		}else{
			map.add(b.map);
		}
	}
	
	/**
	 * Adds k-mers from a Sketch object into this SketchHeap.
	 * Converts sketch keys and merges metadata information.
	 * @param b The Sketch to add k-mers from
	 */
	public void add(Sketch b){
		if(taxID<0){taxID=b.taxID;}
		if(imgID<0){imgID=b.imgID;}
		if(taxName==null){taxName=b.taxName();}
		if(name0==null){name0=b.name0();}
		if(fname==null){fname=b.fname();}
		genomeSizeBases+=b.genomeSizeBases;
		genomeSizeKmers+=b.genomeSizeKmers;
		genomeSequences+=b.genomeSequences;
		if(baseCounts!=null && b.baseCounts!=null){Tools.add(baseCounts, b.baseCounts);}
		set16S(b.r16S());
		set18S(b.r18S());
		
		long[] keys=b.keys;
		int[] counts=b.keyCounts;
		assert(keys.length==b.length()) : keys.length+", "+b.length(); //Otherwise, change to loop through the size
		for(int i=0; i<keys.length; i++){
			long key=Long.MAX_VALUE-keys[i];
			int count=(counts==null ? 1 : counts[i]);
			assert((key>=SketchObject.minHashValue)==(count>0));
			increment(key, count);
		}
	}
	
	/**
	 * Generates a header string containing all metadata for this sketch.
	 * Includes size, coding parameters, k-mer settings, genome statistics,
	 * taxonomic information, and ribosomal RNA sequences.
	 * @return StringBuilder containing formatted header information
	 */
	public StringBuilder toHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#SZ:"+setOrMap.size());
		
		sb.append("\tCD:");
		sb.append(SketchObject.codingArray[SketchObject.CODING]);
		if(SketchObject.deltaOut){sb.append('D');}
		if(SketchObject.aminoOrTranslate()){sb.append('M');}
		if(SketchObject.amino8){sb.append('8');}
		
		sb.append("\tK:").append(SketchObject.k);
		if(SketchObject.k2>0){sb.append(",").append(SketchObject.k2);}
		if(SketchObject.HASH_VERSION>1){sb.append("\tH:").append(SketchObject.HASH_VERSION);}

		if(genomeSizeBases>0){sb.append("\tGS:"+genomeSizeBases);}
		if(genomeSizeKmers>0){sb.append("\tGK:"+genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){sb.append("\tGE:").append(ge);}
		if(genomeSequences>0){sb.append("\tGQ:"+genomeSequences);}
		if(baseCounts!=null && !SketchObject.aminoOrTranslate()){
			sb.append("\tBC:").append(baseCounts[0]).append(',').append(baseCounts[1]).append(',');
			sb.append(baseCounts[2]).append(',').append(baseCounts[3]);
		}
		if(probSum>0){sb.append("\tPC:"+Tools.format("%.4f",probSum/genomeSizeKmers));}
		if(taxID>=0){sb.append("\tID:"+taxID);}
		if(imgID>=0){sb.append("\tIMG:"+imgID);}
		if(taxName!=null){sb.append("\tNM:"+taxName);}
		if(name0!=null){sb.append("\tNM0:"+name0);}
		if(fname!=null){sb.append("\tFN:"+fname);}

		if(r16S!=null){sb.append("\t16S:"+r16S.length);}
		if(r18S!=null){sb.append("\t18S:"+r18S.length);}
		if(r16S!=null){
			sb.append('\n').append("#16S:");
			for(byte b : r16S){sb.append((char)b);}
		}
		if(r18S!=null){
			sb.append('\n').append("#18S:");
			for(byte b : r18S){sb.append((char)b);}
		}
		return sb;
	}
	
	/**
	 * Adds a k-mer value after checking blacklist/whitelist constraints.
	 * Applies filtering rules before attempting to add to the data structure.
	 * @param value The k-mer hash value to potentially add
	 * @return true if the value was added, false if rejected or filtered
	 */
	public boolean checkAndAdd(long value){
		assert(value>=SketchObject.minHashValue);
		
//		if(!heap.hasRoom() && value<=heap.peek()){return false;}
//		if(Blacklist.contains(value)){return false;}
//		if(!Whitelist.contains(value)){return false;}
		
		if(Blacklist.exists() || Whitelist.exists()){
			if(!heap.hasRoom() && value<=heap.peek()){return false;}
			if(Blacklist.contains(value)){return false;}
			if(!Whitelist.containsRaw(value)){return false;}
		}
		
		return add(value);
	}
	
	/**
	 * Calculates the maximum sketch length based on genome size estimates.
	 * Used to determine optimal sketch size for this genome.
	 * @return Maximum recommended sketch length
	 */
	public final int maxLen(){
		return SketchObject.toSketchSize(genomeSizeBases, genomeSizeKmers, genomeSizeEstimate(), SketchObject.targetSketchSize);
	}
	
	/**
	 * Converts the heap contents to a sorted array of k-mer values.
	 * Uses default maximum length and minimum count parameters.
	 * @return Sorted array of k-mer hash values
	 */
	public final long[] toSketchArray(){
		int maxLen=maxLen();
		return toSketchArray(maxLen, minKeyOccuranceCount);
	}
	
	/**
	 * Converts heap to sketch array using specified minimum count threshold.
	 * @param minKeyOccuranceCount_ Minimum occurrence count for k-mer inclusion
	 * @return Sorted array of k-mer hash values meeting the count threshold
	 */
	public final long[] toSketchArray_minCount(int minKeyOccuranceCount_){
		int maxLen=maxLen();
		return toSketchArray(maxLen, minKeyOccuranceCount_);
	}
	
	/**
	 * Converts heap to sketch array with specified maximum length.
	 * @param maxLen Maximum number of k-mers to include in output array
	 * @return Sorted array of k-mer hash values, truncated to maxLen
	 */
	final long[] toSketchArray_maxLen(int maxLen){
		return toSketchArray(maxLen, minKeyOccuranceCount);
	}
	
	/**
	 * Legacy destructive method to convert heap to sketch array.
	 * Empties the heap while creating the output array.
	 * @param maxLen Maximum number of k-mers to extract
	 * @return Sorted array of k-mer hash values
	 */
	private final long[] toSketchArrayOld(int maxLen){//Destructive
		final int initial=heap.size();
		final int len=Tools.min(maxLen, initial);
		final long[] array=KillSwitch.allocLong1D(len);
		
		int toSkip=heap.size()-len;
		for(int i=0; i<toSkip; i++){heap.poll();}
		for(int i=0; i<len; i++){
			array[i]=Long.MAX_VALUE-heap.poll();
		}
		Tools.reverseInPlace(array);
		assert(heap.size()==0) : heap.size()+", "+len+", "+maxLen+", "+initial;
		return array;
	}
	
	/**
	 * Non-destructive conversion of heap/map to sketch array.
	 * Preserves original data structure while creating output array.
	 *
	 * @param maxLen Maximum number of k-mers to include
	 * @param minKeyOccuranceCount_ Minimum occurrence count threshold
	 * @return Sorted array of k-mer hash values meeting criteria
	 */
	private final long[] toSketchArray(int maxLen, int minKeyOccuranceCount_){//Non-destructive
		if(minKeyOccuranceCount_<0){minKeyOccuranceCount_=minKeyOccuranceCount;}
		if(setMode){return toSketchArrayOld(maxLen);}
		long[] keys=map().toArray(minKeyOccuranceCount_);
		for(int i=0; i<keys.length; i++){
//			assert(keys[i]>0) : Arrays.toString(keys);
			keys[i]=Long.MAX_VALUE-keys[i];
//			assert(keys[i]>0) : Arrays.toString(keys);
		}
		Arrays.sort(keys);
		if(keys.length>maxLen){
			keys=Arrays.copyOf(keys, maxLen);
		}
		
//		final LongHeap heap=heap;
//		heap.clear();
//		assert(heap.size()==0) : heap.size()+", "+maxLen;
		return keys;
	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
		return code;
	}
	
	/**
	 * Estimates genome size based on minimum k-mer hash value and sketch size.
	 * Uses statistical properties of MinHash sketching for size estimation.
	 * @return Estimated genome size in k-mers
	 */
	public long genomeSizeEstimate() {
		int size=size();
		if(size==0){return 0;}
		long min=peek();
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min, size));
//		assert(est<30000000) : min+", "+(Long.MAX_VALUE-min)+", "+size+", "+genomeSizeKmers+", "+Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min, size));
		return est;
	}
	
	/**
	 * Estimates genome size using k-mers with specified minimum count.
	 * Filters out low-frequency k-mers for more robust size estimation.
	 * @param minCount Minimum k-mer occurrence count to consider
	 * @return Estimated genome size based on filtered k-mers
	 */
	public long genomeSizeEstimate(int minCount) {
		if(minCount<2){return genomeSizeEstimate();}
		if(size()==0){return 0;}
		long[] min=map.map.getMin(minCount);
		if(min[1]==0){return 0;}
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min[0], (int)min[1]));
		return est;
	}
	
	/** Estimates optimal sketch size for this genome.
	 * @return Recommended sketch size based on genome size estimates */
	public long sketchSizeEstimate(){
		return SketchObject.toSketchSize(genomeSizeBases, genomeSizeKmers, genomeSizeEstimate(), SketchObject.targetSketchSize);
	}

	/**
	 * Checks if a k-mer hash value exists in this sketch.
	 * @param key The k-mer hash value to search for
	 * @return true if the key exists in the sketch
	 */
	public boolean contains(long key) {
		return setOrMap.contains(key);
	}
	
	@Override
	public String toString(){return toHeader().toString();}

	/** Gets the primary name, preferring taxonomic name over original name.
	 * @return The taxonomic name if available, otherwise the original name */
	public String name(){return taxName==null ? name0 : taxName;}
	/** Gets the taxonomic name */
	public String taxName(){return taxName;}
	/** Gets the original sequence name */
	public String name0(){return name0;}
	/** Gets the source filename */
	public String fname(){return fname;}
	/**
	 * Gets base composition counts for A, C, G, T bases.
	 * @param original Whether to return original array or a clone
	 * @return Array of base counts [A, C, G, T] or null if not tracked
	 */
	public long[] baseCounts(boolean original){return baseCounts==null ? null : original ? baseCounts : baseCounts.clone();}
	/** Sets the taxonomic name.
	 * @param s The taxonomic name to set */
	public void setTaxName(String s){taxName=s;}
	/** Sets the original sequence name.
	 * @param s The original name to set */
	public void setName0(String s){name0=s;}
	/** Sets the source filename.
	 * @param s The filename to set */
	public void setFname(String s){fname=s;}
	
	/** Gets the 16S ribosomal RNA sequence */
	public byte[] r16S(){return r16S;}
	/** Gets the length of the 16S rRNA sequence */
	public int r16SLen(){return r16S==null ? 0 : r16S.length;}
	/**
	 * Sets the 16S ribosomal RNA sequence, keeping the highest-scoring one.
	 * Evaluates sequence quality based on length and base composition.
	 * @param b The 16S sequence to potentially store
	 */
	public void set16S(byte[] b){
		if(b==null || b.length<SketchObject.min_SSU_len){return;}
		
		if(r16S==null || score16S(b)>score16S(r16S)){
			r16S=b;
		}
	}

	/**
	 * Scores a 16S sequence based on length similarity to ideal (1533 bp).
	 * @param seq The sequence to score
	 * @return Quality score for the 16S sequence
	 */
	private float score16S(byte[] seq){return scoreSSU(seq, 1533);}
	/**
	 * Scores an 18S sequence based on length similarity to ideal (1858 bp).
	 * @param seq The sequence to score
	 * @return Quality score for the 18S sequence
	 */
	private float score18S(byte[] seq){return scoreSSU(seq, 1858);}
	
	/**
	 * Scores small subunit ribosomal RNA based on length and base quality.
	 * Combines length similarity score with fraction of defined bases.
	 *
	 * @param seq The ribosomal RNA sequence to score
	 * @param idealLen The ideal length for this rRNA type
	 * @return Combined quality score (0.0 to 1.0)
	 */
	private float scoreSSU(byte[] seq, int idealLen){
		float lengthScore=lengthScore(seq.length, idealLen);
		float definedScore=(seq.length-AminoAcid.countUndefined(seq))/(float)seq.length;
		return lengthScore*definedScore;
	}
	/**
	 * Calculates length similarity score between actual and ideal lengths.
	 * @param len Actual sequence length
	 * @param ideal Ideal sequence length
	 * @return Length similarity score (0.0 to 1.0)
	 */
	private float lengthScore(int len, int ideal){return Tools.min(len, ideal)/(float)Tools.max(len, ideal);}
	
	/** Gets the 18S ribosomal RNA sequence */
	public byte[] r18S(){return r18S;}
	/** Gets the length of the 18S rRNA sequence */
	public int r18SLen(){return r18S==null ? 0 : r18S.length;}
	/**
	 * Sets the 18S ribosomal RNA sequence, keeping the highest-scoring one.
	 * Evaluates sequence quality based on length and base composition.
	 * @param b The 18S sequence to potentially store
	 */
	public void set18S(byte[] b){
		if(b==null || b.length<SketchObject.min_SSU_len){return;}
		if(r18S==null || score18S(b)>score16S(r18S)){
			r18S=b;
		}
	}
	
	/**
	 * Determines if this organism is classified as a eukaryote.
	 * Uses taxonomic tree lookup based on the taxonomic ID.
	 * @return true if organism is eukaryotic based on taxonomy
	 */
	boolean isEukaryote(){
		if(taxID<1 || taxID>=SketchObject.minFakeID){return false;}
		if(SketchObject.taxtree==null){return false;}
		return SketchObject.taxtree.isEukaryote((int)taxID);
	}
	
	/** Taxonomic name of the organism */
	private String taxName;
	/** Original sequence identifier */
	private String name0;
	/** Source filename */
	private String fname;
	/** NCBI taxonomic identifier */
	public long taxID=-1;
	/** IMG genome identifier */
	public long imgID=-1;
	/** Genome size in base pairs */
	public long genomeSizeBases=0;
	/** Genome size in k-mers */
	public long genomeSizeKmers=0;
	/** Number of sequences in the genome */
	public long genomeSequences=0;
	/** Base composition counts [A, C, G, T] */
	public final long[] baseCounts=(SketchObject.aminoOrTranslate() ? null : new long[4]);
	/** 16S ribosomal RNA sequence */
	private byte[] r16S;
	/** 18S ribosomal RNA sequence */
	private byte[] r18S;
	/** Sum of k-mer probability scores */
	double probSum=0;

	/** Calculates average probability correctness across all k-mers.
	 * @return Average probability of k-mer correctness (0.0 to 1.0) */
	public float probCorrect(){return probSum<=0 ? 0f : (float)(probSum/Tools.max(genomeSizeKmers, 1f));}
	/** Gets the maximum capacity of the underlying heap */
	public int capacity(){return heap.capacity();}
	/** Checks if the heap has room for more elements */
	public boolean hasRoom(){return heap.hasRoom();}
	/** Gets the minimum k-mer value without removing it */
	public long peek(){return heap.peek();}
	/** Gets the current number of k-mers stored */
	public int size(){return heap.size();}
	/** Gets the underlying hash map for k-mer frequency access */
	public LongHashMap map(){return map.map;}

	/** Clears all k-mer data while preserving metadata */
	public void clear(){
		setOrMap.clear();
	}
	/** Clears the underlying set or map data structure */
	public void clearSet(){
		if(set==null){map.map.clear();}
		else{set.set.clear();}
	}
	/**
	 * Adds a k-mer to the data structure.
	 * @param key The k-mer hash value to add
	 * @return true if the key was added successfully
	 */
	public boolean add(long key){return setOrMap.add(key);}
	/**
	 * Increments the count for a k-mer by the specified amount.
	 * @param key The k-mer hash value to increment
	 * @param incr The amount to increment by
	 * @return The new count for this k-mer
	 */
	public int increment(long key, int incr){return setOrMap.increment(key, incr);}

	/** Set-based storage for k-mers without frequency tracking */
	private final LongHeapSet set;
	/** Map-based storage for k-mers with frequency tracking */
	private final LongHeapMap map;
	/** Unified interface to either set or map storage */
	private final LongHeapSetInterface setOrMap;
	/** Underlying heap structure for maintaining k-mer order */
	public final LongHeap heap;
	/** Minimum occurrence count required for k-mer retention */
	public final int minKeyOccuranceCount;
	/** Determines whether to use LongHeapSet or LongHeapMap */
	public final boolean setMode;
}
