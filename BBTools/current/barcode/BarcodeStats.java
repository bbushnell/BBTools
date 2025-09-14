package barcode;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;

import align2.BandedAlignerConcrete;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.TextFile;
import shared.LineParser2;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

public class BarcodeStats {
	
	public static void main(String[] args) {
		
		String fname=args[0];
		int codesPerRead=(args.length>1 ? Integer.parseInt(args[1]) : 1);
		char delimiter=(args.length>2 ? args[2].charAt(0) : 0);
		String expected=(args.length>3 ? args[3] : null);
		
		BarcodeStats bs=new BarcodeStats(delimiter, codesPerRead, null);
		if(expected!=null) {bs.loadBarcodeList(expected, delimiter, false, false);}
		bs.loadBarcodeCounts(fname);
		if(codesPerRead>1) {
			bs.leftStats=bs.makeLeft();
			bs.rightStats=bs.makeRight();
			bs.leftStats.calcStats();
			bs.rightStats.calcStats();
		}
		bs.calcStats();
		String s=bs.toStats("Barcodes:");
		System.err.println(s);
		if(bs.leftStats!=null){
			System.err.println(bs.leftStats.toStats("Left:     "));
		}
		if(bs.rightStats!=null){
			System.err.println(bs.rightStats.toStats("Right:     "));
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public BarcodeStats(int delimiter_, int barcodesPerRead_, String label_) {
		delimiter=(char)delimiter_;
		assert(delimiter==delimiter_) : "Invalid delimiter; character value "+delimiter_;
		barcodesPerRead=barcodesPerRead_;
		label=label_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Methods    ----------------*/
	/*--------------------------------------------------------------*/
	

	
	/**
	 * Loads a set of barcodes from a set of barcodes or filenames.
	 * They can have a tab after them with more stuff, which will be trimmed.
	 * Modifies the input set.
	 */
	public static LinkedHashSet<String> loadBarcodeSet(LinkedHashSet<String> barcodeSet,
			byte barcodeDelimiter, boolean rcIndex1, boolean rcIndex2){
		if(barcodeSet==null || barcodeSet.isEmpty()) {return null;}
		int barcodeLength=0;
		
		//Handle names
		@SuppressWarnings("unchecked")
		LinkedHashSet<String> tempSet=new LinkedHashSet<String>(barcodeSet);
		barcodeSet.clear();
		for(String s : tempSet) {
			File f=new File(s);
			if(f.exists() && f.isFile()){//Read the file
				TextFile tf=new TextFile(s);
				String[] lines=tf.toStringLines();
				for(String line : lines){
//					assert(false) : "'"+Character.toString(barcodeDelimiter)+"'="+((int)(barcodeDelimiter));
					if(barcodeDelimiter>0 || rcIndex1 || rcIndex2) {
						line=BarcodeStats.fixBarcode(line, barcodeDelimiter, rcIndex1, rcIndex2);
						assert(barcodeDelimiter<1 || line.indexOf(barcodeDelimiter)>=0) : 
							line+", "+barcodeDelimiter;
					}
					assert(barcodeLength==0 || barcodeLength==line.length()) : barcodeLength+", "+line.length()+
						", "+"'"+Character.toString(barcodeDelimiter<0 ? '0' : barcodeDelimiter)+"'="+((int)(barcodeDelimiter))+"\n"+line;
					barcodeLength=line.length();
					barcodeSet.add(line);
				}
			}else{
				if(barcodeDelimiter>0 || rcIndex1 || rcIndex2) {
					s=BarcodeStats.fixBarcode(s, barcodeDelimiter, rcIndex1, rcIndex2);
					assert(barcodeDelimiter<1 || s.indexOf(barcodeDelimiter)>=0) : 
						s+", "+barcodeDelimiter;
				}
				assert(barcodeLength==0 || barcodeLength==s.length());
				barcodeLength=s.length();
				barcodeSet.add(s);
			}//Re-add each key as a literal name
		}
		return barcodeSet;
	}
	
	public int loadBarcodeList(String fname, int forceDelimiter, boolean rcIndex1, boolean rcIndex2){
		if(fname==null){return 0;}
		String[] codes;
		if(new File(fname).exists()) {
			codes=TextFile.toStringLines(fname);
		}else {
			codes=fname.split(",");
		}
		for(int i=0; i<codes.length; i++){
			String s=codes[i];
			linesProcessed++;
			bytesProcessed+=s.length();
			if(!Tools.startsWith(s, '#')){
//				s=removeTab(s);
				//assert(s.indexOf('\t')<0) : "Barcodes should not contain a tab: '"+s+"'";//Although it's fine if they do; just use -da
				s=fixBarcode(s, forceDelimiter, rcIndex1, rcIndex2);
				Barcode b=new Barcode(s);
				expectedCodeList.add(b);
				expectedCodeMap.put(b.name, b);
			}
		}
		return expectedCodeList.size();
	}
	
//	//For the special case of barcode <tab> count, like CountBarcodes2's output
//	//Now handled by fixBarcode
//	private static String removeTab(String s) {
//		int tab=s.indexOf('\t');
//		if(tab<0) {return s;}
//		for(int i=tab+1; i<s.length(); i++) {
//			char c=s.charAt(i);
//			if(!Tools.isDigit(c) && c!='\t') {return s;}
//		}
//		return s.substring(0, tab);
//	}
	
	public static String fixBarcode(String s, int delimiter, boolean rcIndex1, boolean rcIndex2) {
		if(s==null) {return null;}
		final int tab=s.indexOf('\t');
		if(delimiter!='\t' && tab>0) {s=s.substring(0, tab);}
		if(delimiter>0) {s=fixDelimiter(s, delimiter);}
		s=rcompBarcode(s, delimiter, rcIndex1, rcIndex2);
		return s;
	}
	
	public static String rcompBarcode(String s, int delimiter, boolean rcIndex1, boolean rcIndex2) {
		if(s==null || (!rcIndex1 && !rcIndex2)) {return s;}
		int idx=s.indexOf((char)delimiter);
		if(idx<0) {
			return rcIndex1 ? AminoAcid.reverseComplementBases(s) : s;
		}else {
			String a=s.substring(0, idx);
			String b=s.substring(idx+1);
			a=(rcIndex1 ? AminoAcid.reverseComplementBases(a) : a);
			b=(rcIndex2 ? AminoAcid.reverseComplementBases(b) : b);
			return a+((char)delimiter)+b;
		}
	}
	
	public static String fixDelimiter(String s, int forceDelimiter) {
		if(s==null || forceDelimiter<1 || s.indexOf(forceDelimiter)>=0) {return s;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(!Character.isLetter(c)) {
				return s.replace(c, (char)forceDelimiter);
			}
		}
		assert(false) : "Can't find delimiter in '"+s+"'";
		return s;
	}
	
	public long loadBarcodeCounts(String fname){
		LineParser2 lp=new LineParser2('\t');
		if(fname==null){
			assert(false) : "Null filename when loading barcode counts.";
			return 0;
		}
		long added=0;
		long sum=0;
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			linesProcessed++;
			bytesProcessed+=line.length;
			if(!Tools.startsWith(line, '#')){
				lp.set(line);
				String name=lp.parseString();
				long count=lp.parseLong();
				increment(name, count);
				added++;
				sum+=count;
			}
		}
		bf.close();
		return added;
	}
	
	//TODO: Add load from fastq
	
	public BarcodeStats split(int pairnum) {
		assert(barcodesPerRead==2) : barcodesPerRead;
		BarcodeStats split=new BarcodeStats((char)0, 1, pairnum==0 ? "Left" : "Right");
		split.length1=(pairnum==0 ? length1 : length2);
		for(Barcode b : expectedCodeList) {
			split.addExpectedCode(b.getStringForPairnum(pairnum, delimiter));
		}
		for(Entry<String, Barcode> e : codeMap.entrySet()){
			Barcode b=e.getValue();
			String code=b.getStringForPairnum(pairnum, delimiter);
			split.increment(code, b.count());
//			assert(false) : b.count()+", "+codeMap.size()+", "+split.codeMap.get(code);
		}
		return split;
	}
	
	public static int findDelimiter(String fname) {
		TextFile tf=new TextFile(fname);
		int delimiter=-1;
		for(String line=tf.nextLine(); line!=null && delimiter==-1; line=tf.nextLine()){
			if(!Tools.startsWith(line, '#')){
				String b=Tools.tabPattern.split(line)[0];
				delimiter=FileFormat.barcodeDelimiter(b);
				break;
			}
		}
		tf.close();
		return delimiter;
	}
	
	public static BarcodeStats loadStatic(String fname, String expectedBarcodeFile, int expectedBarcodeCount) {
		BarcodeStats bs=loadStatic(fname);
		if(expectedBarcodeFile!=null){
			bs.loadBarcodeList(expectedBarcodeFile, bs.delimiter, false, false);
		}else if(expectedBarcodeCount>0) {
			ArrayList<Barcode> list=bs.toList();
			for(int i=0, max=Tools.min(expectedBarcodeCount, list.size()); i<max; i++) {
				Barcode b=list.get(i);
				bs.addExpectedCode(b.name);
			}
		}
		return bs;
	}
	
	public static BarcodeStats loadStatic(String fname) {
		int delimiter=findDelimiter(fname);
		BarcodeStats bs=new BarcodeStats(delimiter>0 ? delimiter : 0, delimiter>0 ? 2 : 1, null);
		bs.loadBarcodeCounts(fname);
		return bs;
	}
	
	public BarcodeStats makeLeft() {return split(0);}
	public BarcodeStats makeRight() {return split(1);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addExpectedCodes(BarcodeStats bs) {
		assert(expectedCodeList.isEmpty());
		expectedCodeList.addAll(bs.expectedCodeList);
		expectedCodeMap.putAll(bs.expectedCodeMap);
	}
	
	public int addExpectedCode(String s) {
		if(expectedCodeMap.containsKey(s)){return 0;}
		return addExpectedCode(new Barcode(s));
	}
	
	public int addExpectedCode(Barcode b) {
		if(expectedCodeMap.containsKey(b.name)){return 0;}
		expectedCodeMap.put(b.name, b);
		expectedCodeList.add(b);
		return 1;
	}
	
	public void increment(String key, long amt) {
		Barcode b=codeMap.get(key);
		if(b==null){
			b=new Barcode(key);
			codeMap.put(key, b);
		}
		b.increment(amt);
	}
	
	public void increment(String key, long amt, int tile) {
		final String key2=(tile==0 ? key : key+tile);
		Barcode b=codeMap.get(key2);
		if(b==null){
			b=new Barcode(key);
			b.tile=tile;
			codeMap.put(key2, b);
		}else {
			assert(b.tile==tile);
		}
		b.increment(amt);
	}
	
	public void incrementBad(String key, long amt) {
		Barcode b=badPairMap.get(key);
		if(b==null){
			b=new Barcode(key);
			badPairMap.put(key, b);
		}
		b.increment(amt);
	}
	
	public void incrementGood(String key, long amt) {
		Barcode b=goodPairMap.get(key);
		if(b==null){
			b=new Barcode(key);
			goodPairMap.put(key, b);
		}
		b.increment(amt);
	}
	
	public void incrementMismatches(final String s, final long count, LongList transitions) {
		final Barcode b=findClosest(s, 6, 1);
		if(b==null) {return;}
		final String s2=b.name;
		for(int i=0; i<s.length(); i++) {
			final char q=s.charAt(i), r=s2.charAt(i);
			if(Tools.isLetter(q) && r!='N') {
				if(q!=r) {
					incrementMismatch(i, r, q, count, transitions);
				}
			}
		}
//		assert(false) : transitions;
	}
	
	public void incrementMismatch(int pos, int ref, int query, long incr, LongList transitions) {
		int idx=Transition.encode(pos, ref, query);
		transitions.increment(idx, incr);
	}
	
	public ByteBuilder printTransitions() {return printTransitions(transitions, null);}
	
	public static ByteBuilder printTransitions(LongList transitions, ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		ArrayList<Transition> list=new ArrayList<Transition>();
		for(int i=0; i<transitions.size(); i++) {
			long count=transitions.get(i);
			if(count>0) {
				Transition t=Transition.decode(i);
				t.count=count;
				list.add(t);
			}
		}
		Collections.sort(list);
		bb.append("#pos\tref\tquery\tcount\n");
		for(Transition t : list) {
			t.appendTo(bb).nl();
		}
//		assert(false) : transitions+"\n"+bb;
		return bb;
	}
	
	public void calcStats() {//Should only be called after merging, only on master copy
//		System.err.println(codeMap);
//		assert(false) : codeMap.size();
		for(Entry<String, Barcode> e : codeMap.entrySet()) {
			Barcode b=e.getValue();
			totalCodes+=b.count();
			totalCodesU++;
			
			Barcode eb=expectedCodeMap.get(b.name);
			byte mer=-1;
			if(eb!=null){
				eb.increment(b.count());
				expectedCodes+=b.count();
				expectedCodesU++;
				//validArray[2]+=b.count();
			}else{
				if(!expectedCodeList.isEmpty()) {
					incrementMismatches(b.name, b.count(), transitions);
				}
				if(b.countUndefined()>0){
					nCodes+=b.count();
					nCodesU++;
				}else if((mer=b.checkHomopolymer())>=0){
					polymerArray[mer]+=b.count();
					polymerArrayU[mer]++;
				}else if(!expectedCodeList.isEmpty() && barcodesPerRead<2){
					int hdist=calcHdist(b);
					int edist=(hdist>1 && calcEdist ? calcEdist(b) : hdist);
					hdistArray[hdist]+=b.count();
					edistArray[edist]+=b.count();
					hdistArrayU[hdist]++;
					edistArrayU[edist]++;
				}
			}
			
			if(leftStats!=null && !leftStats.expectedCodeList.isEmpty()) {
				assert(rightStats!=null);
				String left=b.leftString(delimiter);
				String right=b.rightString(delimiter);
//				leftStats.increment(left, b.count());
//				rightStats.increment(right, b.count());
				if(eb==null) { 
					final boolean leftMatched=leftStats.expectedCodeMap.containsKey(left), rightMatched=rightStats.expectedCodeMap.containsKey(right);
					final int matched=(leftMatched ? 1 : 0)+(rightMatched ? 1 : 0);
					validArray[matched]+=b.count();
					validArrayU[matched]++;
					if(matched==2){
						leftStats.incrementBad(left, b.count());
						rightStats.incrementBad(right, b.count());
					}
				}else{
					leftStats.incrementGood(left, b.count());
					rightStats.incrementGood(right, b.count());
				}
			}
		}
		
		if(leftStats!=null && !leftStats.expectedCodeList.isEmpty()) {
			badPairFraction=(float)(validArray[2]/(1.0*totalCodes));
			goodPairFraction=(float)(expectedCodes/(1.0*totalCodes));
		}else{
			long bad=0, good=0;//This should just give the same result as above, but for left and right...
			for(Barcode e : expectedCodeList) {
				{
					Barcode bc=badPairMap.get(e.name);
					bad+=(bc==null ? 0 : bc.count());
				}
				{
					Barcode bc=goodPairMap.get(e.name);
					good+=(bc==null ? 0 : bc.count());
				}
			}
			badPairFraction=(float)(bad/(1.0*totalCodes));
			goodPairFraction=(float)(good/(1.0*totalCodes));
		}
	}
	
	public Barcode findClosest(String s, int maxHDist) {
		assert(!expectedCodeList.isEmpty());
		Barcode best=expectedCodeMap.get(s);
		if(best!=null || maxHDist<1) {return best;}
		int min=s.length();
		for(int i=0; i<expectedCodeList.size(); i++) {
			Barcode b=expectedCodeList.get(i);
			int dist=b.hdist(s);
			if(dist<min) {
				best=b;
				min=dist;
			}
		}
		assert(min>0);
		return min<=maxHDist ? best : null;
	}
	
	public Barcode findClosest(String s, int maxHDist, int clearzone) {
		assert(!expectedCodeList.isEmpty());
		Barcode best=expectedCodeMap.get(s);
//		Barcode best2=null;
		
		//Note:  This ignores cases where you have a perfect match to the wrong barcode,
		//or two barcodes are within clearzone of each other.
		if(best!=null || maxHDist<1) {return best;}
		int hdist=s.length();
		int hdist2=hdist;
		assert(best==null); //Otherwise hdist should be 0.
		
		for(int i=0; i<expectedCodeList.size(); i++) {
			Barcode b=expectedCodeList.get(i);
			int d=b.hdist(s);
			if(d<hdist2) {
//				best2=b;
				hdist2=d;
				if(d<hdist) {
//					best2=best;
					hdist2=hdist;
					best=b;
					hdist=d;
				}
			}
		}
		if(hdist<=maxHDist && hdist+clearzone<=hdist2) {return best;}
		return null;
	}
	
	//Slightly different than the below method; it uses a String and has no preconditions
	public int calcHdist(String s) {
		assert(!expectedCodeList.isEmpty());
		if(expectedCodeMap.containsKey(s)) {return 0;}
		int min=s.length();
		for(int i=0; i<expectedCodeList.size(); i++) {
			Barcode b=expectedCodeList.get(i);
			min=Tools.min(min, b.hdist(s));
		}
		return min;
	}
	
	public int calcHdist(Barcode b) {
		assert(!expectedCodeList.isEmpty());
		assert(!expectedCodeMap.containsKey(b.name)) : "Check this first.";
		int min=b.length();
		for(Barcode expected : expectedCodeList) {
			min=Tools.min(min, b.hdist(expected));
			if(min<=1) {
				assert(min==1);
				return min;
			}
		}
		return min;
	}
	
	public int calcEdist(Barcode b) {
		assert(!expectedCodeList.isEmpty());
		assert(!expectedCodeMap.containsKey(b.name)) : "Check this first.";
		int min=b.length();
		for(Barcode expected : expectedCodeList) {
			min=Tools.min(min, b.edist(expected, bandy));
			if(min<=1) {
				assert(min==1);
				return min;
			}
		}
		return min;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Distribution         ----------------*/
	/*--------------------------------------------------------------*/
	
	public BarcodeStats copySpecial() {
		BarcodeStats bs=new BarcodeStats(delimiter, barcodesPerRead, label);
		bs.addExpectedCodes(this);
		return bs;
	}
	
	public void merge(BarcodeStats bs) {
		for(Entry<String, Barcode> e : bs.codeMap.entrySet()) {
			Barcode b=e.getValue();
			increment(b.name, b.count());
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reporting          ----------------*/
	/*--------------------------------------------------------------*/
	
//	private float badPairFraction() {
//		if(leftStats!=null) {
//			return (float)(validArray[2]/(1.0*totalCodes));
//		}
//	}
	
	public boolean isBadPair(String s){
		if(expectedCodeMap.containsKey(s)) {return false;}
		String[] split=split(s, delimiter);
		return leftStats.expectedCodeMap.containsKey(split[0]) && rightStats.expectedCodeMap.containsKey(split[1]);
	}
	
	public boolean hasErrors(String s){
		if(expectedCodeMap.containsKey(s)) {return false;}
		String[] split=split(s, delimiter);
		return !leftStats.expectedCodeMap.containsKey(split[0]) || !rightStats.expectedCodeMap.containsKey(split[1]);
	}
	
	private static String[] split(String s, int delimiter){
		String a=null, b=null;
		if(delimiter<=0) {
			a=s.substring(0, s.length()/2);
			b=s.substring(s.length()/2, s.length());
		}else{
			for(int i=0; a==null && i<s.length(); i++) {
				if(s.charAt(i)==delimiter){
					a=s.substring(0, i);
					b=s.substring(i+1, s.length());
				}
			}
		}
		assert(a!=null && a.length()>0);
		assert(b!=null && b.length()>0);
		return new String[] {a, b};
	}
	
	public int size() {return codeMap.size();}
	
	public ArrayList<Barcode> toList(){
		ArrayList<Barcode> list=new ArrayList<Barcode>(codeMap.size());
		list.addAll(codeMap.values());
		Collections.sort(list);
		return list;
	}
	
	public long hdistSum(int min, int max){
		long sum=0;
		for(int i=min; i<=max && i<hdistArray.length; i++) {
			sum+=hdistArray[i];
		}
		return sum;
	}
	
	public long edistSum(int min, int max){
		long sum=0;
		for(int i=min; i<=max && i<edistArray.length; i++) {
			sum+=edistArray[i];
		}
		return sum;
	}
	
	String pad(long x) {
		String s=""+x;
		while(s.length()<10) {s=s+" ";}
		return s;
	}
	
	public String toStats(String hdr) {
		StringBuilder sb=new StringBuilder();
		sb.append(hdr+"\t"+pad(totalCodes)+"\t("+totalCodesU+" unique)");
		if(barcodesPerRead>1) {
			sb.append("\n-codesPerRead\t"+pad(barcodesPerRead));
			sb.append("\n-delimiter\t"+(delimiter>0 ? Character.toString(delimiter) : ""));
			if(length1>0) {sb.append("\n-bcLength1\t"+pad(length1));}
			if(length2>0) {sb.append("\n-bcLength2\t"+pad(length2));}
		}else {
			if(length1>0) {sb.append("\n-bcLength\t"+pad(length1));}
		}
		if(!expectedCodeList.isEmpty()) {
			sb.append("\n-expected\t"+pad(expectedCodes)+"\t("+expectedCodesU+" unique)");
			sb.append("\n-unexpected\t"+pad(totalCodes-expectedCodes)+"\t("+(totalCodesU-expectedCodesU)+" unique)");
		}
		sb.append("\n-codesWithNs\t"+pad(nCodes)+"\t("+nCodesU+" unique)");
		sb.append("\n-homopolymers\t"+pad(shared.Vector.sum(polymerArray))+"\t("+shared.Vector.sum(polymerArrayU)+" unique)");
		if(polymerArray[0]>0) {sb.append("\n--polyA\t"+polymerArray[0]);}
		if(polymerArray[1]>0) {sb.append("\n--polyC\t"+polymerArray[1]);}
		if(polymerArray[2]>0) {sb.append("\n--polyG\t"+polymerArray[2]);}
		if(polymerArray[3]>0) {sb.append("\n--polyT\t"+polymerArray[3]);}
		if(barcodesPerRead>1 && !expectedCodeList.isEmpty()){
			sb.append("\n-badPair\t"+pad(validArray[2])+"\t("+validArrayU[2]+" unique)"+"\t"+String.format("(%.4f)",badPairFraction));
			sb.append("\n-singleMatch\t"+pad(validArray[1])+"\t("+validArrayU[1]+" unique)"+"\t"+String.format("(%.4f)",(float)(validArray[1]/(1.0*totalCodes))));
			sb.append("\n-neitherMatch\t"+pad(validArray[0])+"\t("+validArrayU[0]+" unique)"+"\t"+String.format("(%.4f)",(float)(validArray[0]/(1.0*totalCodes))));
		}
		if(!expectedCodeList.isEmpty() && barcodesPerRead<2) {
//			sb.append("\n-hDist 1,2,3+\t"+hdistSum(1,1)+", "+hdistSum(2,2)+", "+hdistSum(3,999));
//			if(calcEdist) {
//				sb.append("\n-eDist 1,2,3+\t"+edistSum(1,1)+", "+edistSum(2,2)+", "+edistSum(3,999));
//			}
			sb.append("\n-hDist 1,2,3+\t"+hdistSum(1,1)+"\t"+hdistSum(2,2)+"\t"+hdistSum(3,999));
			if(calcEdist) {
				sb.append("\n-eDist 1,2,3+\t"+edistSum(1,1)+"\t"+edistSum(2,2)+"\t"+edistSum(3,999));
			}
		}
		return sb.toString();
	}
	
	public void printToFile(String fname, boolean overwrite){
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TXT, null, true, overwrite, false, false);
		printToFile(ff);
	}
	
	public void printToFile(FileFormat ff){printToFile(ff, 0);}
	
	public void printToFile(FileFormat ff, long minCount){
		if(ff==null) {return;}
		ArrayList<Barcode> list=toList();
		long sum=0;
		for(Barcode bc : list) {sum+=bc.count();}
		
		try {
			ByteStreamWriter bsw=new ByteStreamWriter(ff);
			bsw.start();
			bsw.println("#Barcodes\t"+sum);
			bsw.println("#Unique\t"+list.size());
			ByteBuilder bb=new ByteBuilder(128);
			for(Barcode bc : list){
				if(bc.count()>=minCount) {
					bc.appendTo(bb).nl();
					bsw.print(bb);
					bb.clear();
				}
			}
			errorState|=bsw.poisonAndWait();
		} catch (Throwable e) {
			System.err.println("ERROR - Could not write barcode file "+ff.name()+": "+e.toString());
			errorState=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final HashMap<String, Barcode> codeMap=new HashMap<String, Barcode>();
	public HashMap<String, Barcode> badPairMap=new HashMap<String, Barcode>();//Only for left or right codes
	public HashMap<String, Barcode> goodPairMap=new HashMap<String, Barcode>();//Only for left or right codes
	public final ArrayList<Barcode> expectedCodeList=new ArrayList<Barcode>();
	public final HashMap<String,Barcode> expectedCodeMap=new HashMap<String,Barcode>();

	public final char delimiter;
	public final int barcodesPerRead;
	public int length1;
	public int length2;
	
	public long linesProcessed=0, bytesProcessed=0;
	
	public long nCodes=0;
	public long expectedCodes=0;
	public long totalCodes=0;
	public float badPairFraction=-1;
	public float goodPairFraction=-1;
	/** 0 is nonmatches, 1 is single matches, 2 is both barcodes match but pair doesn't */
	public long[] validArray=new long[3];
	public long[] polymerArray=new long[4];
	public long[] hdistArray=new long[30];
	public long[] edistArray=new long[30];
	
	public long nCodesU=0;
	public long expectedCodesU=0;
	public long totalCodesU=0;
	public long[] validArrayU=new long[3];
	public long[] polymerArrayU=new long[4];
	public long[] hdistArrayU=new long[30];
	public long[] edistArrayU=new long[30];
	
	public LongList transitions=new LongList();

	public BarcodeStats leftStats=null;
	public BarcodeStats rightStats=null;
	private final BandedAlignerConcrete bandy=new BandedAlignerConcrete(31);
	private String label;
	
	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean calcEdist=false;
	public static boolean errorState=false;
	
}
