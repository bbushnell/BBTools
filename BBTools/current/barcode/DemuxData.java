package barcode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;

import shared.LineParser1;
import shared.LineParser2;
import shared.Parse;
import shared.Tools;
import sketch.Sketch;
import structures.ByteBuilder;

public class DemuxData {
	
	public static void main(String[] args) {
		DemuxData dd=new DemuxData(10, 10, '+');
		dd.deltaBarcodes=true;
		dd.deltaCounts=true;
//		dd.coding=Sketch.A48;
		
		dd.expectedList=new LinkedHashSet<String>();
		dd.codeCounts=new ArrayList<Barcode>();

		dd.expectedList.add("ACGTACGTAC+CGAACGAACG");
		dd.expectedList.add("GGGGGGGGGG+CGAACGAACG");

		dd.codeCounts.add(new Barcode("ACGTACGTAC+CGAACGAACG", 17));
		dd.codeCounts.add(new Barcode("ACGTACGTAC+NNNNNNNNNN", 3));
		dd.codeCounts.add(new Barcode("GGGGGGGGGG+GNNNNNNNNN", 3));
		dd.codeCounts.add(new Barcode("GGGGGGGGGG+NNNNNNNNNN", 3));
		dd.codeCounts.add(new Barcode("GGGGGGGGGG+GNNNNNNNNN", 2));
		dd.codeCounts.add(new Barcode("GGGGGGGGGT+GNNNNNNNNN", 2));
		dd.codeCounts.add(new Barcode("GGGGGGGGGT+TNNNNNNNNN", 2));
		dd.codeCounts.add(new Barcode("GGGGGGGGGG+TCGTACGTTT", 1));
		dd.codeCounts.add(new Barcode("NNNNNNNNNN+NNNNNNNNNN", 1));
		System.err.println("Expected:\n"+dd.expectedList);
		
		
		ArrayList<Barcode> list=(ArrayList<Barcode>)(dd.codeCounts);

		System.err.println("Encoding.");
		ArrayList<byte[]> chunks=dd.encode();
		System.err.println(new String(chunks.get(0)));

		System.err.println("Decoding.");
		DemuxData dd2=new DemuxData(chunks);
		System.err.println("Counts:\n"+dd.codeCounts);
		System.err.println("Result:\n"+dd2.codeCounts);
	}
	
	public DemuxData(int length1_, int length2_, int delimiter_){
		length1=length1_;
		length2=length2_;
		barcodeDelimiter=delimiter_;
	}
	
	public DemuxData(ArrayList<byte[]> chunks) {
		decode(chunks);
	}
	
	public ArrayList<byte[]> encode() {
		ByteBuilder bb=new ByteBuilder();

		ArrayList<byte[]> chunks=new ArrayList<byte[]>(2);
		
		bb.append("#Flags").nl();
		bb.append("length1=").append(length1).nl();
		bb.append("length2=").append(length2).nl();
		bb.append("barcodedelimiter=").append(barcodeDelimiter).nl();
		bb.append("coding=").append(coding).nl();
		bb.append("deltacounts=").append(deltaCounts).nl();
		bb.append("deltabarcodes=").append(deltaBarcodes).nl();
		
		bb.append("maxhdist0=").append(maxHDist0).nl();
		bb.append("minratio0=").append(minRatio0, 6).nl();
		bb.append("minprob0=").append(minProb0, 6).nl();
		bb.append("maxhdist1=").append(maxHDist1).nl();
		bb.append("minratio1=").append(minRatio1, 6).nl();
		bb.append("minprob1=").append(minProb1, 6).nl();
		bb.append("addpolya=").append(addPolyA).nl();
		bb.append("addpolyc=").append(addPolyC).nl();
		bb.append("addpolyg=").append(addPolyG).nl();
		bb.append("addpolyt=").append(addPolyT).nl();
		
		
		bb.nl().append("#Expected\t").append(expectedList.size()).nl();
		for(String s : expectedList) {
			bb.append(s).nl();
		}
		
		if(ENSURE_SORTED) {
			Barcode prev=null;
			boolean sorted=true;
			for(Barcode bc : codeCounts) {
				if(prev!=null) {
					if(bc.compareTo(prev)<0) {
						sorted=false;
						break;
					}
				}
				prev=bc;
				if(!sorted) {
					ArrayList<Barcode> list=new ArrayList<Barcode>(codeCounts);
					Collections.sort(list);
					codeCounts=list;
				}
			}
		}
		
		long total=0;
		bb.nl().append("#Counts\t").append(codeCounts.size()).nl();
		long prevCount=-1;
		String prevName=null;
		final byte[] buffer=new byte[64];
		prevCode1=prevCode2=0;
		for(Barcode bc : codeCounts) {
			appendBarcode(bc, bb, prevCount, prevName, buffer);
			prevCount=bc.count();
			prevName=bc.name;
			if(bb.length>=SEND_BUFFER_MAX_BYTES) {
				total+=bb.length;
				byte[] chunk=bb.toBytes();
				chunks.add(chunk);
				bb.clear();
//				System.err.println("Created chunk "+toStringChunk(chunk, 50));
			}
		}
		bb.append("#End").nl();
		if(bb.length>0) {
			total+=bb.length;
			byte[] chunk=bb.toBytes();
			chunks.add(chunk);
			bb.clear();
//			System.err.println("Created chunk "+toStringChunk(chunk, 50));
		}
		
		System.err.println("Encoded "+total+" bytes as "+chunks.size()+" chunks.");
//		for(byte[] b : chunks) {
//			System.err.println(new String(b));
//		}
//		assert(false);
		
		return chunks;
	}
	
	private ByteBuilder appendBarcode(Barcode bc, ByteBuilder bb, long prevCount, 
			String prevName, final byte[] buffer) {
		final long count=bc.count();
		if(coding==Sketch.RAW) {
			bb.append(bc.name);
		}else if(coding==Sketch.HEX) {
			assert(false) : coding;
		}else if(coding==Sketch.A48) {
			final boolean sameCount=(count==prevCount);
			final long x1=encodeACGTN(bc.name, 0, length1);
			assert(!deltaBarcodes || !sameCount || x1>=prevCode1) : 
				prevCode1+", "+x1+", "+prevName+", "+bc+", "+deltaBarcodes+", "+count+", "+prevCount;
			final boolean match1=(x1==prevCode1);
			if(!SKIP_DUPLICATE || !match1) {
				final long x=x1;
				final long code=(deltaBarcodes && sameCount ? x-prevCode1 : x);
				appendA48(code, bb, buffer);
				prevCode1=x;
			}
			if(length2>0) {
				final long x2=encodeACGTN(bc.name, length1+1, length1+1+length2);
				assert(!deltaBarcodes || !sameCount || !match1 || x2>=prevCode2) : 
					"\npc2="+prevCode2+", x2="+x2+", pn="+prevName+", n="+bc.name+", "
							+ "pcount="+prevCount+", count="+bc.count();
				final boolean match2=(x2==prevCode2);
				final long code=(deltaBarcodes && sameCount && match1 ? x2-prevCode2 : x2);
				bb.tab();
				if(!SKIP_DUPLICATE || !match2) {appendA48(code, bb, buffer);}
				prevCode2=x2;
			}
		}else {
			assert(false) : coding;
		}
		if(count!=prevCount) {
			bb.tab().append(prevCount<1 || !deltaCounts ? count : prevCount-count);
		}
		return bb.nl();
	}
	
	public static final long encodeACGTN(String s, int a, int b) {
		final byte[] table=baseToNumberACGNT;
		//Starting with 1 would simplify things by encoding the length, at the cost of ~1 bit.
		long code=0;
		for(int i=a; i<b; i++) {
			final char c=s.charAt(i);
			final int x=table[c];
			assert(x>=0) : x+", "+c+", "+s+"["+i+"] is not ACGTN";
			code=code*5+x;
			assert(code>=0) : s+", "+a+", "+b;
		}
		return code;
	}
	
	public static final void decodeACGTN(final long code0, final ByteBuilder bb, final int len) {
		final byte[] table=numberToBaseACGNT;
		long code=code0;
		final int from=bb.length();
		for(int i=0; i<len; i++) {
			final int x=(int)(code%5);
			code/=5;
			final byte c=table[x];
			bb.append(c);
		}
		final int to=bb.length();
		bb.reverseInPlace(from, to);
//		assert(false);
		assert(code<=1) : "Code too big for expected length: len="+
			len+", code0="+code0+", code="+code+", term="+bb;
	}
	
	public static final void appendA48(long value, ByteBuilder bb, byte[] temp){
		int i=0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
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
	
	public static long parseA48(byte[] line, int from){
		if(line.length<=from){return 0;}
		long x=0;
		for(int i=from; i<line.length; i++){
			final byte b=line[i];
			if(b<48) {break;}
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	private static String toStringChunk(byte[] chunk, int x) {
		return "'"+new String(chunk, 0, x)+"','"+new String(chunk, chunk.length-x, x)+"'";
	}
	
	public void decode(Collection<byte[]> chunks) {
		System.err.println("Decoding "+chunks.size()+" chunks.");
		int mode=0;
		LineParser2 lp=new LineParser2('\n');
//		expectedList=new ArrayList<Barcode>();
//		codeCounts=new ArrayList<Barcode>();
		int numExpected=0, expectedSeen=0;
		int numCounts=0, countsSeen=0;
		int chunksSeen=0;
		
		ByteBuilder bbChunk=new ByteBuilder(), bbBuffer=new ByteBuilder();
		long prevCount=-1;
		prevCode1=prevCode2=0;
		String prevName=null;
		for(byte[] chunk0 : chunks) {
			chunksSeen++;

			byte[] chunk=chunk0;
			if(!bbChunk.endsWith('\n')) {//This block fixes each chunk so they end with a linebreak
				bbChunk.append(chunk0);
				int offset=0;
				while(!bbChunk.endsWith('\n')) {
					offset++;
					bbChunk.length--;
				}
				chunk=bbChunk.toBytes();
				bbChunk.clear();
				for(; offset>0; offset--) {
					bbChunk.append(chunk0[chunk0.length-offset]);
				}
			}
			
//			System.err.println("Processing chunk "+chunksSeen+": "+new String(chunk, 0, 60)+"..."+
//					new String(chunk, chunk.length-60, 60));
			lp.set(chunk);
			
			//Inverted the above loop for faster processing
			while(mode<1 && lp.hasMore()) {
				lp.advance();
				if(lp.currentFieldLength()>0) {
					final byte[] line=lp.parseByteArrayFromCurrentField();
					if(line[0]=='#') {
						int x=parseMode(line);
						assert(x>mode);
						mode=x;
						assert(mode==FLAGS) : mode+", "+new String(line);
					}else {assert(false) : new String(line);}
				}
			}
			while(mode==FLAGS && lp.hasMore()) {
				lp.advance();
				if(lp.currentFieldLength()>0) {
					final byte[] line=lp.parseByteArrayFromCurrentField();
					tabParser.set(line);
					if(line[0]=='#') {
						int x=parseMode(line);
						assert(x>mode);
						mode=x;
						assert(mode==EXPECTED) : mode+", "+new String(line);
						numExpected=tabParser.parseInt(1);
						expectedList=new LinkedHashSet<String>(numExpected);
					}else{
						boolean x=parseFlag(line);
					}
				}
			}
			while(mode==EXPECTED && lp.hasMore()) {
				lp.advance();
				if(lp.currentFieldLength()>0) {
					final byte[] line=lp.parseByteArrayFromCurrentField();
					tabParser.set(line);
					if(line[0]=='#') {
						int x=parseMode(line);
						assert(x>mode);
						mode=x;
						assert(mode==COUNTS) : mode+", "+new String(line);
						numCounts=tabParser.parseInt(1);
						codeCounts=new ArrayList<Barcode>(numCounts);
					}else{
						expectedSeen++;
						Barcode bc=new Barcode(tabParser.parseString(0));
						expectedList.add(bc.name);
						if(length1<0) {
							length1=bc.length1();
							length2=bc.length2();
						}
					}
				}
			}
			while(mode==COUNTS && lp.hasMore()) {
				lp.advance();
				if(lp.currentFieldLength()>0) {
					final byte[] line=lp.parseByteArrayFromCurrentField();
					if(lp.parseByteFromCurrentField()=='#') {
						int x=parseMode(line);
						assert(x>mode);
						mode=x;
						assert(mode==END) : mode+", "+new String(line);
					}else{
						countsSeen++;
//						tabParser.set(line);
//						Barcode bc=new Barcode(tabParser.parseString(0), tabParser.parseLong(1));
						Barcode bc=parseCountLine(line, prevCount, prevName, bbBuffer);
						prevCount=bc.count();
						prevName=bc.name;
						codeCounts.add(bc);
						if(length1<0) {
							length1=bc.length1();
							length2=bc.length2();
						}
					}
				}
			}
			if(mode==END){break;}
		}
		assert(countsSeen==numCounts) : countsSeen+", "+numCounts+", "+mode;
		assert(mode==END) : mode;
	}
	
	private Barcode parseCountLine(byte[] line, long prevCount, String prevName, ByteBuilder bb) {
		bb.clear();
		tabParser.set(line);
		int term=0;
		tabParser.setBounds(term);
		if(coding==Sketch.A48) {
			final long x1, x2, x3;
			x1=(tabParser.currentFieldLength()<1 ? -1 : parseA48(line, tabParser.a()));
			if(length2>0) {
				term++;
				tabParser.setBounds(term);
				x2=(tabParser.currentFieldLength()<1 ? -1 : parseA48(line, tabParser.a()));
			}else {x2=-1;}
			if(tabParser.hasMore()) {
				term++;
				x3=tabParser.parseLong(term);
			}else {x3=-1;}
			
			final long count=(x3<0 ? prevCount : deltaCounts && prevCount>0 ? prevCount-x3 : x3);
			assert(count>0 && (count<=prevCount || prevCount<0)) : x3+", "+count+", "+prevCount+"\n"+new String(line);
			final long code1=(x1<0 ? prevCode1 : 
				deltaBarcodes && count==prevCount ? prevCode1+x1 : x1);
			final long code2=(x2<0 ? prevCode2 : 
				deltaBarcodes && count==prevCount && code1==prevCode1 ? prevCode2+x2 : x2);
			
			decodeACGTN(code1, bb, length1);
			if(length2>0) {
				if(barcodeDelimiter>0) {bb.append((byte)barcodeDelimiter);}
				decodeACGTN(code2, bb, length2);
			}
			prevCode1=code1;
			prevCode2=code2;
			return new Barcode(bb.toString(), count);
		}else {
			bb.append(line, tabParser.a(), tabParser.b());
		}
		long count=prevCount;
		if(tabParser.hasMore()) {
			term++;
			long x=tabParser.parseLong(term);
			count=prevCount>0 && deltaCounts ? prevCount-x : x;
			assert(count>0) : prevCount+", "+x+", "+count;
		}
		return new Barcode(bb.toString(), count);
	}
	
	public int parseMode(byte[] line) {
		if(Tools.startsWith(line, "#Flags")) {return FLAGS;}
		if(Tools.startsWith(line, "#Expected")) {return EXPECTED;}
		if(Tools.startsWith(line, "#Counts")) {return COUNTS;}
		if(Tools.startsWith(line, "#End")) {return END;}
		return -1;
	}
	
	public boolean parseFlag(byte[] line) {
		equalsParser.set(line);
		String a=equalsParser.parseString(0);
		String b=equalsParser.terms()>1 ? equalsParser.parseString(1) : null;
		boolean c=parseFlag(a, b);
		assert(c) : a+", "+b;
		return c;
	}
	

	
	
	
	public boolean parseFlag(String a, String b){
		
		if(a.equalsIgnoreCase("polya") || a.equalsIgnoreCase("addpolya")){
			addPolyA=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyc") || a.equalsIgnoreCase("addpolyc")){
			addPolyC=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyg") || a.equalsIgnoreCase("addpolyg")){
			addPolyG=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("polyt") || a.equalsIgnoreCase("addpolyt")){
			addPolyT=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("coding") || a.equalsIgnoreCase("encoding")){
			coding=Tools.startsWithDigit(b) ? Integer.parseInt(b) : Tools.indexOf(Sketch.codingArray, b.charAt(0));
		}else if(a.equalsIgnoreCase("length1") || a.equalsIgnoreCase("len1")){
			length1=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("length2") || a.equalsIgnoreCase("len2")){
			length2=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("delimiter") || a.equalsIgnoreCase("barcodedelimiter")){
			barcodeDelimiter=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("deltacounts")){
			deltaCounts=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("deltacodes") || a.equalsIgnoreCase("deltabarcodes")){
			deltaBarcodes=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("delta")){
			deltaBarcodes=deltaCounts=Parse.parseBoolean(b);
		}
		
		else if(a.equals("maxhdist0") || a.equals("hdist0")){
			maxHDist0=Integer.parseInt(b);
		}else if(a.equals("maxhdist") || a.equals("hdist") || a.equals("maxhdist1") || a.equals("hdist1")){
			maxHDist1=Integer.parseInt(b);
		}else if(a.equals("minratio0") || a.equals("ratio0")){
			minRatio0=(float)Parse.parseDoubleKMG(b);
			assert(minRatio0>=0) : minRatio0;
		}else if(a.equals("minratio1") || a.equals("ratio1") || a.equals("ratio") || a.equals("minratio")){
			minRatio1=(float)Parse.parseDoubleKMG(b);
			assert(minRatio1>=0) : minRatio1;
		}else if(a.equals("minprob0")){
			minProb0=(float)Parse.parseDoubleKMG(b);
			minProb0=(minProb0>=0 ? minProb0 : (float)Math.pow(10, minProb0));
			assert(minProb0<=1 && minProb0>=0) : minProb0;
		}else if(a.equals("minprob1") || a.equals("minprob")){
			minProb1=(float)Parse.parseDoubleKMG(b);
			minProb1=(minProb1>=0 ? minProb1 : (float)Math.pow(10, minProb1));
			assert(minProb1<=1 && minProb1>=0) : minProb1;

//		}else if(a.equals("hybrid") || a.equalsIgnoreCase("hybridhdist")){
//			hybridHDist=Tools.startsWithDigit(b) ? Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : -1;
//		}else if(a.equalsIgnoreCase("hybridClearzone") || a.equals("hybridcz")){
//			hybridClearzone=Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}
	
	public Collection<Barcode> codeCounts;
	public LinkedHashSet<String> expectedList;

	private long prevCode1=0;
	private long prevCode2=0;
	
	public int length1=-1;
	public int length2=-1;
	public int barcodeDelimiter='+';
	
	public int maxHDist0=6;
	public float minRatio0=20f;
	public float minProb0=-12;

	public int maxHDist1=6;
	public float minRatio1=1_000_000f; //Possibly 5k for single-ended (10bp), or based on length
	public float minProb1=-5.6f;
	
	public boolean addPolyA=PCRMatrix.addPolyA;
	public boolean addPolyC=PCRMatrix.addPolyC;
	public boolean addPolyG=PCRMatrix.addPolyG;
	public boolean addPolyT=PCRMatrix.addPolyT;
	
	public int SEND_BUFFER_MAX_BYTES=8000000;
	public boolean deltaCounts=deltaCountsDefault;
	public boolean deltaBarcodes=deltaBarcodesDefault;
	
	public static boolean deltaCountsDefault=true;
	public static boolean deltaBarcodesDefault=false;
	
	private LineParser1 equalsParser=new LineParser1('=');
	private LineParser1 tabParser=new LineParser1('\t');
	
	public int type=PCRMatrix.PROB_TYPE;
	public boolean hdistSum=false;
	public int coding=DEFAULT_CODING;
	
	public static int DEFAULT_CODING=Sketch.A48;
	public static boolean ENSURE_SORTED=true;
	public static boolean SKIP_DUPLICATE=true;
	private static final int FLAGS=1, EXPECTED=2, COUNTS=3, END=4;
	
	public static final byte[] numberToBaseACGNT=new byte[] {'A', 'C', 'G', 'N', 'T'};
	public static final byte[] baseToNumberACGNT=makeBaseToNumberACGNT();
	
	private static final byte[] makeBaseToNumberACGNT() {
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)(-1));
		array['A']=array['a']=0;
		array['C']=array['c']=1;
		array['G']=array['g']=2;
		array['N']=array['n']=3;
		array['T']=array['t']=4;
		array['U']=array['u']=4;
		return array;
	}
	
}
