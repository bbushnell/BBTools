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

public class SketchHeap {
	
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
	
	public final int maxLen(){
		return SketchObject.toSketchSize(genomeSizeBases, genomeSizeKmers, genomeSizeEstimate(), SketchObject.targetSketchSize);
	}
	
	public final long[] toSketchArray(){
		int maxLen=maxLen();
		return toSketchArray(maxLen, minKeyOccuranceCount);
	}
	
	public final long[] toSketchArray_minCount(int minKeyOccuranceCount_){
		int maxLen=maxLen();
		return toSketchArray(maxLen, minKeyOccuranceCount_);
	}
	
	final long[] toSketchArray_maxLen(int maxLen){
		return toSketchArray(maxLen, minKeyOccuranceCount);
	}
	
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
	
	public long genomeSizeEstimate() {
		int size=size();
		if(size==0){return 0;}
		long min=peek();
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min, size));
//		assert(est<30000000) : min+", "+(Long.MAX_VALUE-min)+", "+size+", "+genomeSizeKmers+", "+Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min, size));
		return est;
	}
	
	public long genomeSizeEstimate(int minCount) {
		if(minCount<2){return genomeSizeEstimate();}
		if(size()==0){return 0;}
		long[] min=map.map.getMin(minCount);
		if(min[1]==0){return 0;}
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min[0], (int)min[1]));
		return est;
	}
	
	public long sketchSizeEstimate(){
		return SketchObject.toSketchSize(genomeSizeBases, genomeSizeKmers, genomeSizeEstimate(), SketchObject.targetSketchSize);
	}

	public boolean contains(long key) {
		return setOrMap.contains(key);
	}
	
	@Override
	public String toString(){return toHeader().toString();}

	public String name(){return taxName==null ? name0 : taxName;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public String fname(){return fname;}
	public long[] baseCounts(boolean original){return baseCounts==null ? null : original ? baseCounts : baseCounts.clone();}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	public void setFname(String s){fname=s;}
	
	public byte[] r16S(){return r16S;}
	public int r16SLen(){return r16S==null ? 0 : r16S.length;}
	public void set16S(byte[] b){
		if(b==null || b.length<SketchObject.min_SSU_len){return;}
		
		if(r16S==null || score16S(b)>score16S(r16S)){
			r16S=b;
		}
	}

	private float score16S(byte[] seq){return scoreSSU(seq, 1533);}
	private float score18S(byte[] seq){return scoreSSU(seq, 1858);}
	
	private float scoreSSU(byte[] seq, int idealLen){
		float lengthScore=lengthScore(seq.length, idealLen);
		float definedScore=(seq.length-AminoAcid.countUndefined(seq))/(float)seq.length;
		return lengthScore*definedScore;
	}
	private float lengthScore(int len, int ideal){return Tools.min(len, ideal)/(float)Tools.max(len, ideal);}
	
	public byte[] r18S(){return r18S;}
	public int r18SLen(){return r18S==null ? 0 : r18S.length;}
	public void set18S(byte[] b){
		if(b==null || b.length<SketchObject.min_SSU_len){return;}
		if(r18S==null || score18S(b)>score16S(r18S)){
			r18S=b;
		}
	}
	
	boolean isEukaryote(){
		if(taxID<1 || taxID>=SketchObject.minFakeID){return false;}
		if(SketchObject.taxtree==null){return false;}
		return SketchObject.taxtree.isEukaryote((int)taxID);
	}
	
	private String taxName;
	private String name0;
	private String fname;
	public long taxID=-1;
	public long imgID=-1;
	public long genomeSizeBases=0;
	public long genomeSizeKmers=0;
	public long genomeSequences=0;
	public final long[] baseCounts=(SketchObject.aminoOrTranslate() ? null : new long[4]);
	private byte[] r16S;
	private byte[] r18S;
	double probSum=0;

	public float probCorrect(){return probSum<=0 ? 0f : (float)(probSum/Tools.max(genomeSizeKmers, 1f));}
	public int capacity(){return heap.capacity();}
	public boolean hasRoom(){return heap.hasRoom();}
	public long peek(){return heap.peek();}
	public int size(){return heap.size();}
	public LongHashMap map(){return map.map;}

	public void clear(){
		setOrMap.clear();
	}
	public void clearSet(){
		if(set==null){map.map.clear();}
		else{set.set.clear();}
	}
	public boolean add(long key){return setOrMap.add(key);}
	public int increment(long key, int incr){return setOrMap.increment(key, incr);}

	private final LongHeapSet set;
	private final LongHeapMap map;
	private final LongHeapSetInterface setOrMap;
	public final LongHeap heap;
	public final int minKeyOccuranceCount;
	/** Determines whether to use LongHeapSet or LongHeapMap */
	public final boolean setMode;
}
