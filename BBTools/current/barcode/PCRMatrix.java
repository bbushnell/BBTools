package barcode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import barcode.stub.PCRMatrixProb;
import barcode.stub.PCRMatrixProbAbstract;
import barcode.stub.PCRMatrixTile;
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
 * @date March 7, 2024
 *
 */
public abstract class PCRMatrix {

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public static PCRMatrix create(int length1_, int length2_, int delimiter_) {
		return create(matrixType0, length1_, length2_, delimiter_, hdistSum0);
	}
	
	public static PCRMatrix create(int type, int length1_, int length2_, int delimiter_, boolean hdistSum_) {
		if(type==HDIST_TYPE) {
			return new PCRMatrixHDist(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==PROB_TYPE) {
			return new PCRMatrixProb(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==TILE_TYPE) {
			return new PCRMatrixTile(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==HDIST_OLD_TYPE) {
//			return new PCRMatrixHDist_old(length1_, length2_, delimiter_, hdistSum_);
		}else if(type==PROB_OLD_TYPE) {
//			return new PCRMatrixProb_old(length1_, length2_, delimiter_, hdistSum_);
		}
		throw new RuntimeException("Unknown PCRMatrix type "+type);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
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
		}else if(a.equals("probability")){
			matrixType0=(Parse.parseBoolean(b) ? PROB_TYPE : (matrixType0==PROB_TYPE ? HDIST_TYPE : matrixType0));
		}else if(a.equals("bytile") || a.equals("tile")){
			matrixType0=(Parse.parseBoolean(b) ? TILE_TYPE : (matrixType0==TILE_TYPE ? PROB_TYPE : matrixType0));
		}else if(a.equalsIgnoreCase("matrixThreads")){
			matrixThreads=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("devmode")){
			devMode=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("ensuresorted")){
			DemuxData.ENSURE_SORTED=Parse.parseBoolean(b);
		}else if(PCRMatrixHDist.parseStatic(arg, a, b) || 
				PCRMatrixProbAbstract.parseStatic(arg, a, b)){
			//In case of shared flags with different defaults
			PCRMatrixHDist.parseStatic(arg, a, b);
			PCRMatrixProbAbstract.parseStatic(arg, a, b);
		}else{
			return false;
		}
		return true;
	}
	
	public abstract boolean parse(String arg, String a, String b);
	
	public static void postParseStatic(){
		PCRMatrixHDist.postParseStatic();
		PCRMatrixProbAbstract.postParseStatic();
		if(matrixType0==TILE_TYPE) {byTile=true;}
		else {
			assert(byTile==false);
			byTile=false;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           File I/O           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static class MapLine implements Comparable<MapLine> {
		
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
		
		final String observed;
		final String assigned;
		final long count;
		final int hdist;
		float prob=-1;
	}
	
	public final void printAssignmentMap(HashMap<String, String> assignmentMap,
			String mapOut, Collection<Barcode> counts, boolean overwrite, boolean append) {
		readlock();
		HashMap<String, Barcode> countMap=Barcode.barcodesToMap(counts);
		printAssignmentMap(assignmentMap, mapOut, countMap, overwrite, append);
		readunlock();
	}
	
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
	
	public static final void printAssignmentMapStatic(HashMap<String, String> assignmentMap,
			String mapOut, Collection<Barcode> counts, boolean overwrite, boolean append) {
		HashMap<String, Barcode> countMap=Barcode.barcodesToMap(counts);
		printAssignmentMapStatic(assignmentMap, mapOut, countMap, overwrite, append);
	}
	
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
	
	protected static final void fillHDistList(byte[] q, byte[][] list, IntList hdList){
		hdList.clear();
		for(byte[] r : list) {
			hdList.add(hdist(q, r));
		}
	}

	public static final int hdist(byte[] q, byte[] r) {return Barcode.hdist(q, r);}
	
	public final int hdist(int d1, int d2) {return hdistSum ? d1+d2 : Tools.max(d1, d2);}
	
	public final int hdist(Barcode q, String r) {
		return hdist(q.name, r);
	}
	
	public final int hdist(String q, String r) {
		return hdistSum ? Barcode.hdist(q, r) : Tools.max(Barcode.hdistL(q, r), Barcode.hdistR(q, r));
	}
	
	public static final long encodeHDist(long idx, long hdist, long hdist2) {
		return idx|(hdist<<32)|(hdist2<<48);
	}
	
	protected final ArrayList<Barcode> highpass(Collection<Barcode> codeCounts, long minCount) {
		if(minCount<=1) {
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
	
	protected final Barcode findClosestHDist(String s, int maxHDist, int clearzone) {
		return length2<1 ? findClosestSingleHDist(s, maxHDist, clearzone) : 
			findClosestDualHDist(s, maxHDist, clearzone);
	}
	
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
	
	public abstract Barcode findClosest(String s);
	
	public abstract void refine(Collection<Barcode> codeCounts, long minCount);
	
	public abstract HashMap<String, String> makeAssignmentMap(Collection<Barcode> codeCounts, long minCount);
	
	public abstract void populateCounts(ArrayList<Barcode> list, long minCount);
	
	public abstract void makeProbs();

	public abstract void initializeData();
	
	protected abstract boolean valid();// {return true;}

	/*--------------------------------------------------------------*/
	/*----------------          Populating          ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public abstract void populateUnexpected();
	
	public final void populateSplitCodes() {
		DemuxData dd=new DemuxData(length1, length2, delimiter);
		populateSplitCodes(dd);
	}
	
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
	
	public void clearCounts() {
		writelock();
		Tools.fill(counts, 0);
		totalCounted=totalAssigned=totalAssignedToExpected=0;
		for(Barcode bc : expectedList) {bc.setCount(0);}
		writeunlock();
	}
	
	public final float assignedFraction() {
		return (totalAssigned/(1.0f*totalCounted));
	}
	
	public final float expectedFraction() {
		return (totalAssignedToExpected/(1.0f*totalCounted));
	}
	
	public final float chimericFraction() {
		return ((totalAssigned-totalAssignedToExpected)/(1.0f*totalCounted));
	}

	public final void add(Barcode query, Barcode ref) {add(query.name, ref, query.count());}
	public final void add(String query, Barcode ref, long count) {add(query, ref, count, 0);}
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
	
	public ByteBuilder toBytes(ByteBuilder bb) {
		return toBytes(counts, bb);
	}
	
	public ByteBuilder toAccuracy(ByteBuilder bb, boolean excludeUnknown) {
		return toAccuracy(counts, bb, excludeUnknown);
	}
	
	public boolean byTile() {return false;}
	public abstract ByteBuilder toBytesProb(ByteBuilder bb);
	
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
	
	protected final Barcode getBarcode(String s) {return expectedMap.get(s);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	private static String poly(char c, int len) {
		ByteBuilder bb=new ByteBuilder(len);
		for(int i=0; i<len; i++) {bb.append(c);}
		return bb.toString();
	}

	protected final void writelock() {rwlock.writeLock().lock();}
	protected final void writeunlock() {rwlock.writeLock().unlock();}
	protected final void readlock() {rwlock.readLock().lock();}
	protected final void readunlock() {rwlock.readLock().unlock();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	protected final int length;
	protected final int length1;
	protected final int length2;
	protected final int delimiter;
	protected final int delimiterPos;
	protected final int start2;
	protected final int letters;
	protected final boolean hdistSum;
	
	/** [position][Call: A,C,G,T,N][Ref: A,C,G,T,N,Sum] */
	protected final long[][][] counts;
	protected long totalCounted;
	protected long totalAssigned;
	protected long totalAssignedToExpected;
	
	protected ArrayList<Barcode> expectedList;
	protected HashMap<String, Barcode> expectedMap;

	protected byte[][] leftBytes;
	protected byte[][] rightBytes;
	protected byte[][][] splitBytes;

	protected Barcode[] leftCodes;
	protected Barcode[] rightCodes;
	protected Barcode[] allCodes;
	protected HashMap<String, Barcode> allCodesMap;

	public boolean verbose=true;
	public static boolean verbose2=false;
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean errorStateS=false;

	protected static boolean addPolyA=false;
	protected static boolean addPolyC=false;
	protected static boolean addPolyG=false;
	protected static boolean addPolyT=false;
	
	//This is just a maximum; will not go over Shared.threads()
	protected static int matrixThreads=64;
	
	protected static boolean localCounts=true;
	public static boolean byTile=false;
	protected static boolean devMode=false;
	
	protected static boolean hdistSum0=false;
	public static final int HDIST_TYPE=0, PROB_TYPE=1, TILE_TYPE=2, 
			HDIST_OLD_TYPE=7, PROB_OLD_TYPE=8;
	protected static int matrixType0=HDIST_TYPE;
	
	protected static final byte[] baseToNumber=AminoAcid.baseToNumber4;
	protected static final byte[] numberToBase=AminoAcid.numberToBase;
	
}
