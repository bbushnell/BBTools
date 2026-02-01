package sketch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ReadWrite;
import json.JsonObject;
import map.AbstractBitSet;
import map.LongHashMap;
import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import structures.LongPair;
import tax.ImgRecord;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch extends SketchObject implements Comparable<Sketch>, Cloneable {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Sketch clone() {
		Sketch copy=null;
		try {
			copy = (Sketch)super.clone();
		} catch (CloneNotSupportedException e) {
//			e.printStackTrace();
			throw new RuntimeException(e);
		}
		copy.compareBitSet=null;
		copy.indexBitSet=null;
		return copy;
	}
	
	//Array should already be hashed, sorted, unique, subtracted from Long.MAX_VALUE, then reversed.
	/**
	 * Constructs a sketch with k-mer keys and metadata.
	 * Uses default values for taxonomy and genome size information.
	 *
	 * @param keys_ Sorted k-mer hash keys
	 * @param keyCounts_ Count of each k-mer (may be null)
	 * @param baseCounts_ Base composition counts [A,C,G,T] (may be null)
	 * @param r16S_ 16S rRNA sequence (may be null)
	 * @param r18S_ 18S rRNA sequence (may be null)
	 * @param meta_ Additional metadata strings (may be null)
	 */
	public Sketch(long[] keys_, int[] keyCounts_, long[] baseCounts_, byte[] r16S_, byte[] r18S_, ArrayList<String> meta_){
		this(keys_, keyCounts_, baseCounts_, r16S_, r18S_, -1, -1, -1, -1, -1, -1, null, null, null, meta_);
	}
	
	/**
	 * Constructs a sketch from a SketchHeap with minimum count filtering.
	 *
	 * @param heap The SketchHeap containing k-mer data
	 * @param clearFname Whether to clear the filename from the heap
	 * @param keepCounts Whether to preserve k-mer count information
	 * @param meta_ Additional metadata strings
	 */
	public Sketch(SketchHeap heap, boolean clearFname, boolean keepCounts, ArrayList<String> meta_){
		this(heap, clearFname, keepCounts, meta_, -1);
	}
	
	/**
	 * Constructs a sketch from a SketchHeap with specified minimum occurrence count.
	 * Filters k-mers to include only those appearing at least minKeyOccuranceCount times.
	 * Preserves k-mer counts if requested and not in set mode.
	 *
	 * @param heap The SketchHeap containing k-mer data
	 * @param clearFname Whether to clear the filename from the heap
	 * @param keepCounts Whether to preserve k-mer count information
	 * @param meta_ Additional metadata strings
	 * @param minKeyOccuranceCount Minimum occurrence count for k-mers to include
	 */
	public Sketch(SketchHeap heap, boolean clearFname, boolean keepCounts, ArrayList<String> meta_, int minKeyOccuranceCount){
		this(heap.toSketchArray_minCount(minKeyOccuranceCount), null, heap.baseCounts(false), heap.r16S(), heap.r18S(), (int)heap.taxID, heap.imgID,
				heap.genomeSizeBases, heap.genomeSizeKmers, heap.genomeSequences,
				heap.probCorrect(), heap.taxName(), heap.name0(), heap.fname(), meta_);
		assert(keyCounts==null);
		if(!heap.setMode && keepCounts){
			LongHashMap map=heap.map();
			keyCounts=new int[keys.length];
			for(int i=0; i<keys.length; i++){
				int count=map.get(Long.MAX_VALUE-keys[i]);
				assert(count>0) : keys[i]+" -> "+count+"\n"+Arrays.toString(map.values())+"\n"+Arrays.toString(map.keys());
				keyCounts[i]=count;
			}
		}
		if(heap.setMode){heap.clearSet();}
		heap.clear(clearFname);
//		System.err.println("size="+size+", genome="+this.genomeSize+", m"); : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (counts==null)+", "+heap.setMode+", "+keepCounts;
	}

	/**
	 * Full constructor for creating a sketch with all possible metadata.
	 * Initializes all fields including taxonomy information, genome size estimates,
	 * and various identifiers. Handles IMG record lookup if available.
	 *
	 * @param keys_ Sorted k-mer hash keys
	 * @param keyCounts_ Count of each k-mer (may be null)
	 * @param baseCounts_ Base composition counts [A,C,G,T] (may be null)
	 * @param r16S_ 16S rRNA sequence (may be null)
	 * @param r18S_ 18S rRNA sequence (may be null)
	 * @param taxID_ NCBI taxonomy ID
	 * @param imgID_ IMG database ID
	 * @param gSizeBases_ Genome size in bases
	 * @param gSizeKmers_ Genome size in k-mers
	 * @param gSequences_ Number of sequences in genome
	 * @param probCorrect_ Probability of correct base calls
	 * @param taxName_ Taxonomic name
	 * @param name0_ Primary sequence name
	 * @param fname_ Source filename
	 * @param meta_ Additional metadata strings
	 */
	public Sketch(long[] keys_, int[] keyCounts_, long[] baseCounts_, byte[] r16S_, byte[] r18S_, int taxID_, long imgID_, long gSizeBases_, 
			long gSizeKmers_, long gSequences_, double probCorrect_, String taxName_, String name0_, String fname_, ArrayList<String> meta_){
//		Exception e=new Exception(baseCounts_+", "+Arrays.toString(baseCounts_));
//		e.printStackTrace();
		keys=keys_;
		keyCounts=keyCounts_;
//		fixKeyCounts();
		baseCounts=baseCounts_;
		r16S=r16S_;
		r18S=r18S_;
		assert(keyCounts==null || keys==null || keyCounts.length==keys.length) : (keys==null ? "null" : keys.length)+", "+(keyCounts==null ? "null" : keyCounts.length);
		taxID=taxID_;
		imgID=imgID_;
		sketchID=nextSketch.getAndIncrement();
		genomeSizeBases=gSizeBases_;
		genomeSizeKmers=gSizeKmers_;
		genomeSequences=gSequences_;
		probCorrect=probCorrect_<=0 ? 0f : (float)probCorrect_;
		meta=meta_;
		
		taxName=fix(taxName_);
		name0=fix(name0_);
		fname=fix(fname_);
		if(fname!=null && (fname.startsWith("stdin.") || fname.endsWith(".sketch") || fname.endsWith(".sketch.gz"))){fname=null;}
		
//		if(k2>0){
//			if(useToValue2) {
//				int count=0;
//				for(long x : array){
//					count+=(int)(x&1);
//				}
//				k1Count=count;
//			}else {
//				k1Count=array.length/2;
//			}
//		}else{
//			k1Count=array.length;
//		}
		
		if(ImgRecord.imgMap!=null && imgID>=0 && taxID<0){
			ImgRecord record=ImgRecord.imgMap.get(imgID);
			if(record!=null){
				if(record.name!=null && taxName==null){taxName=record.name;}
				taxID=record.taxID;
			}
		}
	}
	
	/** Loads SSU rRNA sequences from the SSU map if available.
	 * Only loads for valid taxonomy IDs below the fake ID threshold. */
	void loadSSU(){
		if(taxID>0 && taxID<minFakeID){
			if(r16S==null && SSUMap.r16SMap!=null){
				r16S=SSUMap.r16SMap.get(taxID);
			}
			if(r18S==null && SSUMap.r18SMap!=null){
				r18S=SSUMap.r18SMap.get(taxID);
			}
		}
//		new Exception().printStackTrace();
//		System.err.println(taxID+", "+ssu+", "+SSUMap.ssuMap.size());
	}
	
	/**
	 * Adds a metadata string to this sketch after fixing formatting.
	 * Creates the metadata list if it doesn't exist.
	 * @param s Metadata string to add
	 */
	void addMeta(String s){
		s=fixMeta(s);
		if(s==null){return;}
		if(meta==null){meta=new ArrayList<String>(1);}
		meta.add(s);
	}
	
	/**
	 * Sets the metadata list for this sketch.
	 * Should only be called when no metadata exists.
	 * @param list The metadata list to set
	 */
	void setMeta(ArrayList<String> list){
		assert(meta==null);
		meta=list;
	}
	
	/**
	 * Fixes string formatting by replacing tabs with spaces.
	 * @param s String to fix
	 * @return Fixed string or null if input was null
	 */
	public static String fix(String s){
		if(s==null){return null;}
		return s.replace('\t', ' ');
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges k-mer keys from another sketch into this one.
	 * Maintains sorted order and limits result to maxlen keys.
	 * Both sketches must not have key counts.
	 *
	 * @param other The sketch to merge with
	 * @param maxlen Maximum number of keys to retain
	 */
	public void add(Sketch other, int maxlen){
		assert(keyCounts==null && other.keyCounts==null);
		final long[] a=keys;
		final long[] b=other.keys;
		if(maxlen<1){
			assert(false);
			maxlen=1000000;
		}
		LongList list=new LongList(Tools.min(maxlen, a.length+b.length));
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//match
				list.add(ka);
				i++;
				j++;
			}else if(ka<kb){
				list.add(ka);
				i++;
			}else{
				list.add(kb);
				j++;
			}
			if(list.size()>=maxlen){break;}
		}
		
		if(keys.length==list.size()){
			for(int i=0; i<list.size; i++){
				keys[i]=list.array[i];
			}
		}else{
			keys=list.toArray();
		}
	}

	/**
	 * Resizes the key arrays to the specified length.
	 * Only reduces size, does not expand.
	 * @param newSize New size for key arrays
	 */
	public void resize(int newSize) {
		if(newSize>=length()){return;}
		keys=Arrays.copyOf(keys, newSize);
//		fixKeyCounts(newSize);
		if(keyCounts!=null){
			for(int i=0; i<newSize; i++)
			keyCounts=Arrays.copyOf(keyCounts, newSize);
		}
	}
	
	/**
	 * Fixes key counts by removing count array if maximum count is less than 2.
	 * Updates maxCount field with the highest count found.
	 * @return The maximum count value
	 */
	int fixKeyCounts(){
		return fixKeyCounts(keyCounts==null ? maxCount : keyCounts.length);
	}
	
	/**
	 * Fixes key counts checking up to maxLen entries.
	 * Removes count array if maximum count is less than 2.
	 * @param maxLen Maximum number of entries to check
	 * @return The maximum count value
	 */
	int fixKeyCounts(int maxLen){
		int max=0;
		if(keyCounts!=null){
			for(int i=0; i<maxLen && max<2; i++){
				max=Tools.max(keyCounts[i], max);
			}
		}
		if(max<2){keyCounts=null;}
		maxCount=Tools.max((maxCount>0 ? 1 : 0), max);
		return maxCount;
	}

	/**
	 * Removes blacklisted k-mers from this sketch.
	 * Filters both keys and counts arrays to exclude blacklisted entries.
	 * @return Number of k-mers removed
	 */
	public int applyBlacklist() {
		assert(blacklist!=null);
		LongList keylist=new LongList(keys.length);
		IntList countlist=(keyCounts==null ? null : new IntList(keys.length));
		int removed=0;
		for(int i=0; i<keys.length; i++){
			long key=keys[i];
			if(!Blacklist.contains(Long.MAX_VALUE-key)){
				keylist.add(key);
				if(keyCounts!=null){countlist.add(keyCounts[i]);}
			}else{
				removed++;
			}
		}
		if(keylist.size()!=keys.length){
			keys=keylist.toArray();
			if(keyCounts!=null){keyCounts=countlist.toArray();}
		}
//		if(removed>0){
//			fixKeyCounts();
//		}
		return removed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Set Operations        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a new sketch containing the intersection of two sketches.
	 * Result inherits metadata from the second sketch.
	 *
	 * @param sa First sketch
	 * @param sb Second sketch
	 * @return New sketch with shared k-mers, or null if no intersection
	 */
	public static final Sketch intersection(Sketch sa, Sketch sb){
		Sketch shared=intersection(sa.keys, sb.keys, sa.keyCounts);
		if(shared!=null){
			shared.taxID=sb.taxID;
			shared.taxName=sb.taxName;
			shared.name0=sb.name0;
			shared.fname=sb.fname;
			shared.meta=sb.meta;
			shared.imgID=sb.imgID;
			shared.spid=sb.spid;
		}
		return shared;
	}
	
	/**
	 * Creates a sketch from the intersection of two k-mer key arrays.
	 * Preserves counts from the first array if available.
	 *
	 * @param a First k-mer array
	 * @param b Second k-mer array
	 * @param aCounts Count array for first k-mer array (may be null)
	 * @return New sketch with shared k-mers, or null if no matches
	 */
	public static final Sketch intersection(long[] a, long[] b, int[] aCounts){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=new IntList();
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(aCounts!=null){
					il.add(aCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				i++;
			}else{
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null, null, null, null, null);
	}
	
	/**
	 * Creates a new sketch containing the union of two sketches.
	 * Result inherits metadata from the first sketch.
	 *
	 * @param sa First sketch
	 * @param sb Second sketch
	 * @return New sketch with combined k-mers from both inputs
	 */
	public static final Sketch union(Sketch sa, Sketch sb){
		Sketch shared=union(sa.keys, sb.keys, sa.keyCounts, sb.keyCounts, sa.baseCounts, sb.baseCounts, sa.r16S, sb.r16S, sa.r18S, sb.r18S);
		if(shared!=null){
			shared.taxID=sa.taxID;
			shared.taxName=sa.taxName;
			shared.name0=sa.name0;
			shared.fname=sa.fname;
			shared.meta=sa.meta;
			shared.imgID=sa.imgID;
			shared.spid=sa.spid;
		}
		return shared;
	}
	
	/**
	 * Creates a sketch from the union of two k-mer key arrays with full metadata.
	 * Combines counts and base counts, selects longer SSU sequences.
	 *
	 * @param a First k-mer array
	 * @param b Second k-mer array
	 * @param aCounts Count array for first k-mer array
	 * @param bCounts Count array for second k-mer array
	 * @param aBaseCounts Base counts for first sketch
	 * @param bBaseCounts Base counts for second sketch
	 * @param a16S 16S sequence from first sketch
	 * @param b16S 16S sequence from second sketch
	 * @param a18S 18S sequence from first sketch
	 * @param b18S 18S sequence from second sketch
	 * @return New sketch with combined data from both inputs
	 */
	public static final Sketch union(long[] a, long[] b, int[] aCounts, int[] bCounts, long[] aBaseCounts, long[] bBaseCounts, 
			byte[] a16S, byte[] b16S, byte[] a18S, byte[] b18S){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=(aCounts==null || bCounts==null ? null : new IntList());
		byte[] r16S=(b16S==null ? a16S : a16S==null ? b16S : a16S.length>b16S.length ? a16S : b16S);
		byte[] r18S=(b18S==null ? a18S : a18S==null ? b18S : a18S.length>b18S.length ? a18S : b18S);
		long[] baseCounts=(aBaseCounts==null || bBaseCounts==null ? null : new long[aBaseCounts.length]);
		if(baseCounts!=null){
			for(int k=0; k<baseCounts.length; k++) {baseCounts[k]+=(aBaseCounts[k]+bBaseCounts[k]);}
		}
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(il!=null){
					il.add(aCounts[i]+bCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				ll.add(ka);
				if(il!=null){
					il.add(aCounts[i]);
				}
				i++;
			}else{
				ll.add(kb);
				if(il!=null){
					il.add(bCounts[i]);
				}
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null, baseCounts, r16S, r18S, null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Filtering          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks if this sketch passes metadata filtering criteria.
	 * @param params Display parameters containing filter criteria
	 * @return true if sketch passes all metadata filters
	 */
	boolean passesMeta(DisplayParams params) {
		return passesMeta(params.requiredMeta, params.bannedMeta, /*params.requiredTaxid, params.bannedTaxid,*/ params.requiredMetaAnd);
	}

	/**
	 * Checks if this sketch passes specified metadata filtering criteria.
	 * Evaluates required and banned metadata tags with AND/OR logic.
	 *
	 * @param requiredMeta List of required metadata tags (may be null)
	 * @param bannedMeta List of banned metadata tags (may be null)
	 * @param requiredMetaAnd Whether all required tags must be present (AND) vs any (OR)
	 * @return true if sketch passes all specified filters
	 */
	boolean passesMeta(ArrayList<String> requiredMeta, ArrayList<String> bannedMeta, /*IntList requiredTaxid, IntList bannedTaxid,*/ boolean requiredMetaAnd) {
		assert(requiredMeta!=null || bannedMeta!=null /*|| requiredTaxid!=null || bannedTaxid!=null*/);
		assert(requiredMeta==null || requiredMeta.size()>0);
		assert(bannedMeta==null || bannedMeta.size()>0);

//		if(requiredTaxid!=null && !requiredTaxid.contains(taxID)){return false;}
//		if(bannedTaxid!=null && bannedTaxid.contains(taxID)){return false;}
		
		if(requiredMeta!=null){
			if(meta==null){return false;}
			for(String tag : requiredMeta){
				if(meta.contains(tag)){
					if(!requiredMetaAnd){break;}
				}else if(requiredMetaAnd){
					return false;
				}
			}
		}
		if(bannedMeta!=null && meta!=null){
			for(String tag : bannedMeta){
				if(!meta.contains(tag)){
					return false;
				}
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Counts matching k-mers between this sketch and another sketch.
	 * Wrapper method that delegates to the static version.
	 *
	 * @param other Sketch to compare against
	 * @param buffer Buffer to store comparison results
	 * @param present BitSet tracking which k-mers have been seen
	 * @param fillPresent Whether to populate the present BitSet
	 * @param taxHits Array mapping k-mer indices to taxonomy hits
	 * @param contamLevel Taxonomy level for contamination analysis
	 * @return Number of matching k-mers
	 */
	public int countMatches(Sketch other, CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){
		return countMatches(keys, other.keys, keyCounts, other.keyCounts, refHitCounts, other.taxID, buffer, present, fillPresent, taxHits, contamLevel);
	}
	
	/**
	 * Counts matching k-mers between two sorted k-mer arrays with detailed statistics.
	 * Performs comprehensive analysis including contamination detection, taxonomy hits,
	 * and depth calculations. Uses two-pointer technique for efficient comparison.
	 *
	 * @param a First k-mer array (query)
	 * @param b Second k-mer array (reference)
	 * @param aCounts Count array for first k-mer array
	 * @param bCounts Count array for second k-mer array
	 * @param refHitCounts Array tracking hits per k-mer across references
	 * @param bid Reference taxonomy ID for contamination analysis
	 * @param buffer Buffer to store detailed comparison statistics
	 * @param present BitSet tracking which k-mers have been seen
	 * @param fillPresent Whether to populate or read from the present BitSet
	 * @param taxHits Array mapping k-mer indices to taxonomy hits
	 * @param contamLevel Taxonomy level for contamination grouping
	 * @return Number of matching k-mers between the arrays
	 */
	public static final int countMatches(long[] a, long[] b, int[] aCounts, int[] bCounts, int[] refHitCounts, int bid,
			CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){

		
//		if(verbose2){System.err.println("fillPresent: "+fillPresent+", "+present);}
//		assert(fillPresent) : bid+", "+minFakeID+", "+(taxHits!=null);
		
		if(bid>0 && bid<minFakeID && taxHits!=null){
			bid=taxtree.getIdAtLevelExtended(bid, contamLevel);
		}else{
			bid=-1;
		}
		
//		assert(false) : (buffer==null)+", "+fillPresent+", "+present.cardinality();
		assert(a.length>0 && b.length>0);
		
		//Kmers hitting this reference
		int matches=0;
		
		//Kmers hitting this reference and others
		int multiMatches=0;
		
		//Kmers hitting nothing
		int noHits=0;
		
		//Kmers hitting some organism but not this reference
		int contamHits=0;
		
		//Kmers hitting something in this taxa, but not this reference
		int sameTax=0;
		
		//Kmers hitting this organism and no other taxa
		int unique2=0;
		
		//Kmers hitting only this taxa but not this organism (this count may not include everything due to early exit)
		int unique3_temp=0;
		
		//Kmers hitting multiple organisms but not this reference
		int multiContamHits=0;
		
		//Sum of query counts for shared kmers
		long depthSum=0;
		
		//Sum of query counts for shared kmers divided by ref counts for those kmers
		double depthSum2=0;//Slow, but necessary.
		
		//Number of times matching query keys occurred in all references
		long refHitSum=0;

		//Matches from k1
		int k1hits=0;
		
		//Query kmers from k1
		int k1seenQ=0;
		
		//Reference kmers from k1
		int k1seenR=0;
		int i=0, j=0;
		assert(present==null || present.capacity()==a.length);
//		assert(false) : buffer.rbs.capacity()+", "+buffer.rbs+", "+present;
		if(present!=null){
			if(fillPresent){
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					final int bit=(int)(ka&1);
					if(ka==kb){
						present.increment(i);
						matches++;
						k1hits+=bit;
						k1seenQ+=bit;
						k1seenR+=bit;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
						i++;
						j++;
					}else if(ka<kb){
						i++;
						k1seenQ+=bit;
					}else{
						j++;
						k1seenR+=bit;
					}
				}
			}else{
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					final int bit=(int)(ka&1);
					if(ka==kb){
						final int count=present.getCount(i);
						if(count>1){
							multiMatches++;
						}
						
						matches++;
						k1hits+=bit;
						k1seenQ+=bit;
						k1seenR+=bit;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique2++;}
						}
						
						i++;
						j++;
					}else if(ka<kb){
						k1seenQ+=bit;
						final int count=present.getCount(i);
						if(count>0){
							contamHits++;
							if(count>1){
								multiContamHits++;
							}
						}else{
							noHits++;
						}
						
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null){
								if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
								for(int tid : taxHitsRow){
									if(tid==bid){
										sameTax++;
										break;
									}
								}
							}
						}
						
						i++;
					}else{
						k1seenR+=bit;
						j++;
					}
				}
				
				//For the remaining query kmers, we don't know whether the reference sketch would have shared them had it been longer.
				//This section can be disabled to prevent them from being displayed.
				if(bid>0 && i<a.length-1){
					for(; i<a.length; i++){
						int[] taxHitsRow=taxHits[i];
						if(taxHitsRow!=null){
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
						}
					}
				}
			}
		}else{
			for(; i<a.length && j<b.length; ){
				final long ka=a[i], kb=b[j];
				final int bit=(int)(ka&1);
//				final int count=present.getCount(i);
				if(ka==kb){
					matches++;
					k1hits+=bit;
					k1seenQ+=bit;
					k1seenR+=bit;
					if(aCounts!=null){depthSum+=aCounts[i];}
					if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
					i++;
					j++;
				}else if(ka<kb){
					i++;
					k1seenQ+=bit;
				}else{
					j++;
					k1seenR+=bit;
				}
			}
		}
		
		if(k2<1){
			k1hits=matches;
			k1seenQ=i;
			k1seenR=j;
		}
		
//		if(taxHits!=null){
//			System.err.println("matches="+matches+", noHits="+noHits+", contamHits="+contamHits+", sameTax="+sameTax+", multiContamHits="+multiContamHits);
//		}

//		assert(bid<1 || unique2>=(matches-multiMatches)) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches;
//		assert(matches<1000 || multiMatches==0) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches+", "+fillPresent;
		
		if(buffer!=null){
//			System.err.println("*A) "+matches+", "+multiMatches+", "+unique2+", "+unique3_temp);
//			new Exception().printStackTrace();
//			assert(k1hits<=(matches-k1hits)) : k1hits+", "+matches;
			buffer.set(matches, multiMatches, unique2, unique2+unique3_temp, noHits, 
					contamHits, contamHits-sameTax, multiContamHits, i, j, 
					a.length, b.length, depthSum, depthSum2, refHitSum, k1hits, k1seenQ, k1seenR);
		}
		return matches;
	}
	
//	public float identity(Sketch b, float[] ret){
//		if(ret!=null){Arrays.fill(ret, 0);}
//		return identityWeighted(array, b.array, ret);
//	}
//
//	public static float identity(long[] a, long[] b){
//		int matches=countMatches(a, b);
//		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
//	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
//		System.err.println(code+", "+gSize+", "+taxID+", "+imgID+", "+name0);
		return code;
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		int x=taxName.compareTo(b.taxName);
		if(x!=0){return x;}
		if(name0!=null && b.name0!=null){return name0.compareTo(b.name0);}
		return name0!=null ? 1 : b.name0!=null ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	/**
	 * Tests equality with another sketch using compareTo.
	 * @param b Sketch to compare against
	 * @return true if sketches are equal
	 */
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a formatted header string containing sketch metadata.
	 * @return ByteBuilder containing the header */
	public ByteBuilder toHeader(){
		ByteBuilder sb=new ByteBuilder();
		return toHeader(sb);
	}
	
	/**
	 * Appends a formatted header to the provided ByteBuilder.
	 * Header includes size, coding, k-mer parameters, genome information,
	 * taxonomy data, and SSU sequences.
	 *
	 * @param bb ByteBuilder to append header to
	 * @return The same ByteBuilder with header appended
	 */
	public ByteBuilder toHeader(ByteBuilder bb){
		bb.append("#SZ:").append(keys.length);
		bb.append("\tCD:");
		bb.append(codingArray[CODING]);
		if(deltaOut){bb.append('D');}
		if(keyCounts!=null){bb.append('C');}
		if(aminoOrTranslate()){bb.append('M');}
		if(amino8){bb.append('8');}

		bb.append("\tK:").append(k);
		if(k2!=0){bb.append(",").append(k2);}
		if(HASH_VERSION>1){bb.append("\tH:").append(HASH_VERSION);}
//		if(maxCount>0){bb.append("\tMC:").append(maxCount);}
		
		if(genomeSizeBases>0){bb.append("\tGS:").append(genomeSizeBases);}
		if(genomeSizeKmers>0){bb.append("\tGK:").append(genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){bb.append("\tGE:").append(ge);}
		if(genomeSequences>0){bb.append("\tGQ:"+genomeSequences);}
		if(baseCounts!=null && !aminoOrTranslate()){
			bb.append("\tBC:").append(baseCounts[0]).append(',').append(baseCounts[1]).append(',');
			bb.append(baseCounts[2]).append(',').append(baseCounts[3]);
		}
		if(probCorrect>0){bb.append("\tPC:"+Tools.format("%.4f", probCorrect));}
		if(taxID>=0){bb.append("\tID:").append(taxID);}
		if(imgID>=0){bb.append("\tIMG:").append(imgID);}
		if(spid>0){bb.append("\tSPID:").append(spid);}
		if(fname!=null){bb.append("\tFN:").append(fname);}
		if(taxName!=null){bb.append("\tNM:").append(taxName);}
		if(name0!=null){bb.append("\tNM0:").append(name0);}
		if(meta!=null){
			for(String s : meta){
				bb.append("\tMT_").append(s);
			}
		}
		
		if(r16S!=null){
			bb.append("\t16S:").append(r16S.length);
		}
		if(r18S!=null){
			bb.append("\t18S:").append(r18S.length);
		}
		if(r16S!=null){
			bb.append('\n');
			bb.append("#16S:").append(r16S);
		}
		if(r18S!=null){
			bb.append('\n');
			bb.append("#18S:").append(r18S);
		}
		return bb;
	}
	
	/** Converts the sketch to its serialized byte representation.
	 * @return ByteBuilder containing the complete sketch data */
	public ByteBuilder toBytes(){
		return toBytes(new ByteBuilder());
	}
	
	/**
	 * Appends the sketch data to the provided ByteBuilder.
	 * Uses delta encoding if enabled for compact representation.
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder with sketch data appended
	 */
	public ByteBuilder toBytes(ByteBuilder bb){
		if(CODING==A48 && deltaOut){return toBytesA48D(bb);}
		long prev=0;
		toHeader(bb);
		bb.append("\n");
		byte[] temp=null;
		if(CODING==A48){temp=KillSwitch.allocByte1D(12);}
		for(int i=0; i<keys.length; i++){
			long key=keys[i];
			int count=(keyCounts==null ? 1 : keyCounts[i]);
			long x=key-prev;
			if(CODING==A48){
				appendA48(x, bb, temp);
				if(count>1){
					bb.append('\t');
					appendA48(count-1, bb, temp);
				}
				bb.append('\n');
			}else if(CODING==HEX){
				bb.append(Long.toHexString(x)).append('\n');
			}else if(CODING==RAW){
				bb.append(x).append('\n');
			}else{
				assert(false);
			}
			if(deltaOut){prev=key;}
		}
		return bb;
	}
	
	//This is to make the common case fast
	/**
	 * Fast path for serializing sketch data using A48 coding with delta encoding.
	 * Optimized version of toBytes for the most common encoding format.
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder with sketch data appended
	 */
	private ByteBuilder toBytesA48D(ByteBuilder bb){
		assert(CODING==A48 && deltaOut);
		long prev=0;
		toHeader(bb);
		bb.append("\n");
		final byte[] temp=KillSwitch.allocByte1D(12);

		if(keyCounts==null){
			for(int i=0; i<keys.length; i++){
				long key=keys[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, bb, temp);
					bb.append('\n');
				}
				prev=key;
			}
		}else{
			for(int i=0; i<keys.length; i++){
				long key=keys[i];
				int count=keyCounts[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, bb, temp);
					if(count>1){
						bb.append('\t');
						appendA48(count-1, bb, temp);
					}
					bb.append('\n');
				}
				prev=key;
			}
		}
		return bb;
	}
	
	/**
	 * Encodes a long value using base-48 encoding and appends to ByteBuilder.
	 * Uses a temporary byte array to avoid allocations.
	 *
	 * @param value Value to encode
	 * @param bb ByteBuilder to append to
	 * @param temp Temporary byte array for encoding
	 */
	public static final void appendA48(long value, ByteBuilder bb, byte[] temp){
		int i=0;
//		long value=value0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			bb.append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				bb.append((char)(temp[i]+48));
			}
		}
	}
	
	/**
	 * Converts a long value to its base-48 string representation.
	 * @param value Value to convert
	 * @return Base-48 encoded string
	 */
	public static final String toA48(long value){
		int i=0;
//		long value=value0;
		StringBuilder sb=new StringBuilder(12);
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			sb.append((char)(b+48));
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			sb.reverse();
		}
		return sb.toString();
	}
	
	@Override
	public String toString(){
		return toBytes().toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       String Parsing         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a base-48 encoded string to a long value.
	 * @param line Base-48 encoded string
	 * @return Decoded long value
	 */
	public static long parseA48(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=6;
			long c=line.charAt(i);
			x|=(c-48);
		}
		return x;
	}
	
	/** Parses coverage too */
	public static long parseA48C(String line, IntList covList){
		if(line.length()==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length();
		for(; i<len; i++){
			long c=line.charAt(i);
			if(c<48){break;}
			key<<=6;
			key|=(c-48);
		}
		for(i++; i<len; i++){
			long c=line.charAt(i);
			cov<<=6;
			cov|=(c-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	/**
	 * Parses a hexadecimal string to a long value.
	 * Handles negative values indicated by leading minus sign.
	 * @param line Hexadecimal string
	 * @return Decoded long value
	 */
	public static long parseHex(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=4;
			x|=hexTable[line.charAt(i)];
		}
		if(line.charAt(0)=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Array Parsing         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a base-48 encoded byte array to a long value.
	 * @param line Base-48 encoded byte array
	 * @return Decoded long value
	 */
	public static long parseA48(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	/**
	 * Parses a nucleotide sequence string to find the maximum k-mer hash.
	 * @param line Nucleotide sequence string
	 * @return Maximum k-mer hash value from the sequence
	 */
	public static long parseNuc(String line){
		return parseNuc(line.getBytes());
	}
	
	/** Returns the maximal key in the sequence */
	public static long parseNuc(byte[] bases){
		if(bases.length<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		long kmer=0, rkmer=0;
		int len=0;
		
		long key=-1;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
//				long z=Tools.max(kmer, rkmer);
				final long hashcode=hash(kmer, rkmer);
				key=Tools.max(key, hashcode);
			}
		}
		return key<minHashValue ? -1 : Long.MAX_VALUE-key;
	}
	
	/** Parses coverage too */
	public static long parseA48C(byte[] line, IntList covList){
		if(line.length==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length;
		for(; i<len; i++){
			long b=line[i];
			if(b<48){break;}
			key<<=6;
			key|=(b-48);
		}
		for(i++; i<len; i++){
			long b=line[i];
			cov<<=6;
			cov|=(b-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	/**
	 * Parses a hexadecimal byte array to a long value.
	 * @param line Hexadecimal byte array
	 * @return Decoded long value
	 */
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a base-48 encoded ByteBuilder to a long value.
	 * @param bb ByteBuilder containing base-48 encoded data
	 * @return Decoded long value
	 */
	public static long parseA48(ByteBuilder bb){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(len==0){return 0;}
		long x=0;
		for(int i=0; i<len; i++){
			x<<=6;
			x|=(((long)line[i])-48);
		}
		return x;
	}
	
	/**
	 * Parses a nucleotide sequence from ByteBuilder to find maximum k-mer hash.
	 * @param bb ByteBuilder containing nucleotide sequence
	 * @return Maximum k-mer hash value from the sequence
	 */
	public static long parseNuc(ByteBuilder bb){
		return parseNuc(bb.toBytes());
	}
	
	/** Parses coverage too */
	public static long parseA48C(ByteBuilder bb, IntList covList){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(len==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0;
		for(; i<len; i++){
			long b=line[i];
			if(b<48){break;}
			key<<=6;
			key|=(b-48);
		}
		for(i++; i<len; i++){
			long b=line[i];
			cov<<=6;
			cov|=(b-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	/**
	 * Parses a hexadecimal ByteBuilder to a long value.
	 * @param bb ByteBuilder containing hexadecimal data
	 * @return Decoded long value
	 */
	public static long parseHex(ByteBuilder bb){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(line.length==0){return 0;}
		long x=0;
		for(int i=0; i<len; i++){
			x<<=4;
			x|=hexTable[line[i]];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Estimates genome size based on the largest k-mer hash and sketch size.
	 * Returns the minimum of stored genome size and calculated estimate.
	 * @return Estimated genome size in k-mers
	 */
	public long genomeSizeEstimate() {
		return keys.length==0 ? 0 : Tools.min(genomeSizeKmers, genomeSizeEstimate(keys[keys.length-1], keys.length));
	}
	
	/**
	 * Estimates genome size using only k-mers with at least the specified count.
	 * Filters to high-confidence k-mers for more accurate size estimation.
	 * @param minCount Minimum k-mer count to include in estimation
	 * @return Estimated genome size in k-mers based on filtered data
	 */
	public long genomeSizeEstimate(int minCount) {
		if(minCount<2){return genomeSizeEstimate();}
		if(length()==0){return 0;}
		long max=0;
		int num=0;
		for(int i=0; i<keyCounts.length; i++){
			if(keyCounts[i]>=minCount){
				max=keys[i];
				num++;
			}
		}
		if(max==0){return 0;}
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(max, num));
		return est;
	}
	
	/** Calculates GC content from base composition counts.
	 * @return GC fraction (0.0 to 1.0), or 0 if no base counts available */
	public float gc(){
		if(baseCounts==null){return 0;}
		long at=baseCounts[0]+baseCounts[3];
		long gc=baseCounts[1]+baseCounts[2];
		return gc/(Tools.max(1.0f, at+gc));
	}
	
	/** Extracts the core filename without path and extension.
	 * @return Core filename, or null if no filename set */
	public String filePrefix(){return ReadWrite.stripToCore(fname);}
	/**
	 * Gets the best available name for this sketch.
	 * Prefers taxonomy name, then sequence name, then filename.
	 * @return The most appropriate name for display
	 */
	public String name(){return taxName!=null ? taxName : name0!=null ? name0 : fname;}
	/** Gets the taxonomic name */
	public String taxName(){return taxName;}
	/** Gets the primary sequence name */
	public String name0(){return name0;}
	/** Gets the source filename */
	public String fname(){return fname;}
	/** Gets the number of k-mer keys in this sketch */
	public int length(){return keys.length;}
	/** Sets the taxonomic name.
	 * @param s Taxonomic name to set */
	public void setTaxName(String s){taxName=s;}
	/** Sets the primary sequence name.
	 * @param s Sequence name to set */
	public void setName0(String s){name0=s;}
	/** Sets the source filename, stripping path information.
	 * @param s Filename to set */
	public void setFname(String s){
//		assert(!s.endsWith("sketch")) : s; //123
		fname=(s==null ? s : ReadWrite.stripPath(s));
	}
	
//	public float k1Fraction(){return k1Count/(float)array.length;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Assorted           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a k-mer count histogram from this sketch.
	 * Groups k-mers by count value and returns sorted pairs of (count, frequency).
	 * @return List of (count, frequency) pairs sorted by count value
	 */
	public ArrayList<LongPair> toKhist() {
		HashMap<Long, LongPair> map=new HashMap<Long, LongPair>();
		for(int i=0; i<keyCounts.length; i++){
			int a=keyCounts[i];
			Long key=(long)a;
			LongPair value=map.get(key);
			if(value==null){
				value=new LongPair();
				value.a=a;
				map.put(key, value);
			}
			value.b++;
		}
		ArrayList<LongPair> list=new ArrayList<LongPair>(map.size());
		list.addAll(map.values());
		Collections.sort(list);
		return list;
	}
	
//	static boolean warned=false;
	/**
	 * Adds SSU rRNA sequences from the SSU map based on taxonomy.
	 * Handles prokaryote/eukaryote preferences and restrictions for 16S/18S sequences.
	 * Requires taxonomy tree for organism classification.
	 */
	public void addSSU(){
		if(useSSUMapOnly){r16S=r18S=null;}
		if(taxID<1 || taxID>=minFakeID){return;}
		if(!SSUMap.hasMap()){return;}
		if(taxtree==null/* && !warned*/){
//			warned=true;
//			System.err.println("*** Warning - no taxtree loaded. A taxtree is recommended when adding SSUs. ***");
			assert(false) : "Please set a path to the taxtree when adding SSUs.";
		}
		final boolean euk=taxtree!=null && (preferSSUMapForEuks || useSSUMapOnlyForEuks || ban16SForEuks) ? taxtree.isEukaryote(taxID) : false;
		final boolean prok=taxtree!=null && ban18SForProks ? taxtree.isProkaryote(taxID) : false;
		
		final boolean mapOnly=(useSSUMapOnly || (useSSUMapOnlyForEuks && euk));
		final boolean preferMap=(mapOnly || preferSSUMap || (preferSSUMapForEuks && euk));
		if(mapOnly){r16S=r18S=null;}
		
		if(SSUMap.r16SMap!=null && (r16S==null || preferMap)){
			byte[] x=SSUMap.r16SMap.get(taxID);
			if(x!=null){r16S=x;}
		}
		if(SSUMap.r18SMap!=null && (r18S==null || preferMap)){
			byte[] x=SSUMap.r18SMap.get(taxID);
			if(x!=null){r18S=x;}
		}
		if(ban18SForProks && prok){r18S=null;}
		if(ban16SForEuks && euk){r16S=null;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            BitSet            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates bit sets for contamination tracking if enabled.
	 * @param printContam Whether contamination printing is enabled
	 * @param index Whether to create an index bit set
	 */
	public void makeBitSets(boolean printContam, boolean index){
		assert(compareBitSet==null && indexBitSet==null);
		if(!printContam){return;}
		compareBitSet=AbstractBitSet.make(length(), bitSetBits);
		if(index){indexBitSet=AbstractBitSet.make(length(), bitSetBits);}
	}
	
	/** Adds another bit set to the comparison bit set.
	 * @param rbs Bit set to add */
	public void addToBitSet(AbstractBitSet rbs){
		compareBitSet.add(rbs);
	}
	
	/** Gets the comparison bit set */
	public AbstractBitSet compareBitSet(){return compareBitSet;}
	
	/** Gets the index bit set */
	public AbstractBitSet indexBitSet(){return indexBitSet;}
	
	/** Merges the index bit set into the comparison bit set.
	 * Clears the index bit set after merging. */
	public void mergeBitSets(){
		assert(!mergedBitSets);
		if(compareBitSet!=null && indexBitSet!=null){
			compareBitSet.setToMax(indexBitSet);
		}
		indexBitSet=null;
		mergedBitSets=true;
	}
	
	/** Returns whether bit sets have been merged */
	public boolean merged(){return mergedBitSets;}

	/** Gets the length of the 16S rRNA sequence */
	public int r16SLen(){return r16S==null ? 0 : r16S.length;}
	/** Gets the length of the 18S rRNA sequence */
	public int r18SLen(){return r18S==null ? 0 : r18S.length;}
	/** Gets the 16S rRNA sequence */
	public byte[] r16S(){return r16S;}
	/** Gets the 18S rRNA sequence */
	public byte[] r18S(){return r18S;}
	/** Returns whether this sketch has any SSU rRNA sequences */
	public boolean hasSSU(){return r16S!=null || r18S!=null;}
	/**
	 * Checks if this sketch shares SSU sequences with another sketch.
	 * @param b Sketch to compare SSU sequences with
	 * @return true if both sketches have compatible SSU sequences
	 */
	public boolean sharesSSU(Sketch b){
		return (r16S!=null && b.r16S!=null) || (r18S!=null && b.r18S!=null);
	}
	
	public void set16S(byte[] seq) {
		assert(r16S==null || r16S.length<seq.length) : this+"\n"+new String(seq);
		r16S=seq;
	}
	
	public void set18S(byte[] seq) {
		assert(r18S==null || r18S.length<seq.length) : this+"\n"+new String(seq);
		r18S=seq;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Stores sorted hashcodes, ascending, as Long.MAX_VALUE-(raw hashcode) */
	public long[] keys;
	/** Stores kmer (hashcode) observation counts */
	public int[] keyCounts;
	/** Stores base (ACGTN) counts */
	final long[] baseCounts;
	/** 16S rRNA sequence data */
	private byte[] r16S;
	/** 18S rRNA sequence data */
	private byte[] r18S;
	/** NCBI taxonomy identifier */
	public int taxID;
	/** Maximum k-mer observation count in this sketch */
	public int maxCount;
	/** Unique sketch identifier based on loading order */
	int sketchID;//Based on loading order
	/** Number of sequences in the genome */
	public final long genomeSequences;
	/** Genome size in bases */
	public final long genomeSizeBases;
	/** Genome size in k-mers */
	public final long genomeSizeKmers;
	/** Probability of correct base calls */
	public final float probCorrect;
//	public final int k1Count; //Number of keys made from k1 rather than k2
	public String taxName;
	public String name0;
	/** Source filename */
	private String fname;
	/** Additional metadata strings */
	ArrayList<String> meta;
	
	//TODO: These should move to SketchResults.
	/** BitSet used for comparison operations */
	private AbstractBitSet compareBitSet; //Used for comparison
	/** BitSet used for indexing operations */
	private AbstractBitSet indexBitSet;
	
	
	//Extended information
	/** IMG database identifier */
	public long imgID=-1;
	/** Species identifier */
	public long spid=-1;
//	public String seqUnitName=null;
	
	/** Flag indicating whether bit sets have been merged */
	private boolean mergedBitSets=false; //TODO: Temporary for debugging
	/** Tracks the number of reference sketches sharing each kmer.
	 * Should be set to null when no longer needed. */
	private int[] refHitCounts;

	/** Gets the reference hit counts array */
	public int[] refHitCounts(){return refHitCounts;}
	/** Clears the reference hit counts array to free memory */
	public void clearRefHitCounts(){
		refHitCounts=null;
	}
	/** Sets the reference hit counts array.
	 * @param x Array of hit counts per k-mer */
	public void setRefHitCounts(int[] x) {
		refHitCounts=x;
		assert(x!=null);
	}
	
	/** Atomic counter for generating unique sketch IDs */
	private static AtomicInteger nextSketch=new AtomicInteger(1);
	

	/*--------------------------------------------------------------*/
	
	/** TODO: Move up */
	public final void setFrom(JsonObject all) {
		JsonObject top=null;
		if(all!=null && all.jmapSize()>0) {
			for(String key : all.jmap.keySet()){
				JsonObject hit=all.jmap.get(key);
				if(top==null) {top=hit;}
			}
		}
		if(top!=null) {
			taxID=top.getLong("TaxID").intValue();
			taxName=top.getString("taxName");
		}
	}
	
}
