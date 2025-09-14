package barcode;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import align2.BandedAligner;
import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

public class Barcode implements Comparable<Barcode> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public Barcode(String s){this(s, 0, 1, 0);}
	public Barcode(String s, long c){this(s, c, 1, 0);}
	public Barcode(String s, long c, int e){this(s, c, e, 0);}
	public Barcode(String s, long c, int e, int t){
		name=s;
		count=c;
		expected=e;
		tile=t;
		assert(expected==0 || expected==1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets one of two barcodes for dual indexes */
	public Barcode getBarcodeForPairnum(int pairnum, char delimiter) {
		return new Barcode(getStringForPairnum(pairnum, delimiter), count);
	}
	public String getStringForPairnum(int pairnum, char delimiter) {
		assert(pairnum==0 || pairnum==1);
		final int pos, pos2;
		if(delimiter>0){//If there is a delimiter
			pos=name.indexOf(delimiter);
			pos2=pos+1;
			assert(pos>=0) : pos+", "+Character.toString(delimiter)+", "+name;
		}else{
			boolean even=((name.length()&1)==0);
			pos=(name.length()/2)+1;
			pos2=(even ? pos : pos+1);
		}
		String s=(pairnum==0 ? name.substring(0, pos) : name.substring(pos2, name.length()));
		return s;
	}

	public Barcode left(char delimiter){return getBarcodeForPairnum(0, delimiter);}
	public Barcode right(char delimiter){return getBarcodeForPairnum(1, delimiter);}
	public String leftString(char delimiter){return getStringForPairnum(0, delimiter);}
	public String rightString(char delimiter){return getStringForPairnum(1, delimiter);}
	
	/*--------------------------------------------------------------*/
	
	public int countUndefined() {
		int sum=0;
		for(int i=0; i<length(); i++){
			char c=charAt(i);
			sum+=(Tools.isLetter(c) && !AminoAcid.isFullyDefined(c)) ? 1 : 0;
		}
		return sum;
	}
	
	public boolean isHomopolymer() {return isHomopolymer(name, name.charAt(0));}
	public boolean isHomopolymer(char h) {return isHomopolymer(name);}
	public static boolean isHomopolymer(String name) {return isHomopolymer(name, name.charAt(0));}
	public static boolean isHomopolymer(String name, char h) {
		for(int i=0; i<name.length(); i++) {
			char c=name.charAt(i);
			if(c!=h && Tools.isLetter(c)) {return false;}
		}
		return true;
	}
	
	/** Returns -1 if false, and the homopolymer character 2-bit encoding if true */
	public byte checkHomopolymer() {
		int match=0;
		int nonletter=0;
		final byte mer=AminoAcid.baseToNumber[charAt(0)];
		if(mer<0) {return mer;}
		for(int i=1; i<length(); i++){
			char c=charAt(i);
			byte n=AminoAcid.baseToNumber[c];
			if(n==mer) {
				match++;
			}else if(Tools.isLetter(c)) {
				return -1;
			}else {
				nonletter++;
			}
		}
		return (match>0 && nonletter<=1 ? mer : 0);
	}

	public static final boolean nonLetter(int a, int b) {
		return !Tools.isLetter(a) && !Tools.isLetter(b);
	}
	
	public int hdist(final Barcode b){return hdist(name, b.name);}
	public int hdist(final String b){return hdist(name, b);}
	public int hdist(final byte[] b){return hdist(name, b);}
	public int hdistL(final Barcode b){return hdistL(name, b.name);}
	public int hdistL(final String b){return hdistL(name, b);}
	public int hdistR(final Barcode b){return hdistR(name, b.name);}
	public int hdistR(final String b){return hdistR(name, b);}
	public static int hdist(final byte[] a, final byte[] b){
		assert(a.length==b.length) : "'"+new String(a)+"', '"+new String(b)+"'";
		final int min=Tools.min(a.length, b.length);
		int subs=0;
		for(int i=0; i<min; i++){
			subs+=(a[i]==b[i] ? 0 : 1);
//			subs+=(a[i]==b[i] || nonLetter(a[i], b[i]) ? 0 : 1);
		}
		return subs;
	}
	public static int hdist(final String a, final String b){
		//The second clause is to allow tile numbers to be present
		assert(a.length()==b.length() || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
				&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		final int min=Tools.min(a.length(), b.length());
		int subs=0;
		for(int i=0; i<min; i++){
			final char ca=a.charAt(i), cb=b.charAt(i);
			subs+=(ca==cb ? 0 : 1);
//			subs+=(ca==cb || nonLetter(ca, cb) ? 0 : 1);
		}
		return subs;
	}
	public static int hdist(final String a, final byte[] b){
		assert(a.length()==b.length || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
				&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		final int min=Tools.min(a.length(), b.length);
		int subs=0;
		for(int i=0; i<min; i++){
			final int ca=a.charAt(i), cb=b[i];
			subs+=(ca==cb ? 0 : 1);
//			subs+=(ca==cb || nonLetter(ca, cb) ? 0 : 1);
		}
		return subs;//+Tools.absdif(length(), b.length());
	}
	public static int hdistL(final String a, final String b, int len1){
		assert(a.length()==b.length() || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
				&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		int subs=0;
		for(int i=0; i<len1; i++){
			final char ca=a.charAt(i), cb=b.charAt(i);
			subs+=(ca==cb ? 0 : 1);
		}
		return subs;
	}
	public static int hdistR(final String a, final String b, int len2){
		assert(a.length()==b.length() || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
				&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		int subs=0;
		final int minlen=Tools.min(a.length(), b.length());
		for(int i=minlen-1, min=minlen-len2; i>=min; i--){
			final char ca=a.charAt(i), cb=b.charAt(i);
			subs+=(ca==cb ? 0 : 1);
		}
		return subs;
	}
	public static int hdistL(final String a, final String b){
		assert(a.length()==b.length() || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
				&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		int subs=0;
		for(int i=0, max=a.length(); i<max; i++){
			final char ca=a.charAt(i), cb=b.charAt(i);
			if(!Tools.isLetter(ca)) {break;}
			subs+=(ca==cb ? 0 : 1);
		}
		return subs;
	}
	public static int hdistR(final String a, final String b){
		assert(a.length()==b.length() || (Tools.endsWithLetter(a)!=Tools.endsWithLetter(b) 
			&& Tools.startsWithLetter(a) && Tools.startsWithLetter(b)))
			: "'"+a+"', '"+b+"'";
		final int minlen=Tools.min(a.length(), b.length());
		int subs=0;
		for(int i=minlen-1; i>=0; i--){
			final char ca=a.charAt(i), cb=b.charAt(i);
			if(!Tools.isLetter(ca)) {
				//At least one of them should not be a tile number
				assert(!Tools.isDigit(ca));
				break;
			}
			subs+=(ca==cb ? 0 : 1);
		}
		return subs;
	}

	public int edist(final Barcode b, BandedAligner bandy){return edist(b.name, bandy);}
	public int edist(final String b, BandedAligner bandy){
		int dist=hdist(name, b);
		if(dist>1){dist=bandy.alignForward(getBytes(), b.getBytes(), 0, 0, length(), true);}
		return dist;
	}
	//Bandy is, e.g., new BandedAlignerConcrete(21);
	
	/*--------------------------------------------------------------*/
	
	public static HashMap<String, Barcode> barcodesToMap(Collection<Barcode> codes) {
		HashMap<String, Barcode> countMap=new HashMap<String, Barcode>();
		for(Barcode b : codes) {countMap.put(b.name, b);}
		return countMap;
	}
	
	public static ArrayList<Barcode> summateAssignments(HashMap<String, String> assignmentMap,
			ArrayList<Barcode> expectedCodeList, HashMap<String, Barcode> countMap) {
		ArrayList<String> list=new ArrayList<String>(expectedCodeList.size());
		for(Barcode b : expectedCodeList) {
			if(b.expected==1) {list.add(b.name);}
		}
		return summateAssignments(assignmentMap, list, countMap);
	}
	
	public static ArrayList<Barcode> summateAssignments(HashMap<String, String> assignmentMap,
			Collection<String> expectedCodeList, Collection<Barcode> counts) {
		HashMap<String, Barcode> countMap=barcodesToMap(counts);
		return summateAssignments(assignmentMap, expectedCodeList, countMap);
	}
	
	public static ArrayList<Barcode> summateAssignments(HashMap<String, String> assignmentMap,
			Collection<String> expectedCodeList, HashMap<String, Barcode> countMap) {

		HashMap<String, Barcode> sumMap=new HashMap<String, Barcode>();
		for(String s : expectedCodeList) {
			sumMap.put(s, new Barcode(s));
		}
		for(Entry<String, String> e : assignmentMap.entrySet()) {
			String observed=e.getKey(), assigned=e.getValue();
			Barcode sum=sumMap.get(assigned);
			if(sum==null) {
				sumMap.put(assigned, sum=new Barcode(assigned));//Optional.  This shouldn't really happen.
			}
			Barcode count=countMap.get(observed);
			if(count!=null) {
				sum.increment(count);
			}
		}
		ArrayList<Barcode> list=new ArrayList<Barcode>(sumMap.size());
		list.addAll(sumMap.values());
		Collections.sort(list);
		return list;
	}

	//52866.4.475040.GAGGCCGCCA-TTATCTAGCT.filter-DNA.fastq.gz
	public static String parseBarcodeFromFname(String fname) {
		String[] split=Tools.dotPattern.split(fname);
		for(String s : split) {
			if(isBarcode(s)) {return s;}
			else if("UNKNOWN".equalsIgnoreCase(s)) {return s;}
		}
//		assert(false) : "Can't find barcode in filename "+fname+"\n"+Arrays.toString(split);
		return null;
	}

	public static final boolean isBarcodeSymbol(char x) {return AminoAcid.isACGTN(x) || x=='-' || x=='+';}
	public static final boolean isBarcodeSymbol(byte x) {return AminoAcid.isACGTN(x) || x=='-' || x=='+';}
	
	public static boolean isBarcode(String s) {
		int bases=0;
		int delimiters=0;
		if(s.length()<6) {return false;}
		for(int i=0; i<s.length(); i++) {
			char c=s.charAt(i);
			if(AminoAcid.isACGTN(c)) {bases++;}
			else if(c=='-' || c=='+') {delimiters++;}
			else {return false;}
		}
		return bases>=6 && delimiters<=1;
	}
	
	public static final boolean containsOnlyBarcodeSymbols(String s) {
		for(int i=0; i<s.length(); i++) {
			if(!isBarcodeSymbol(s.charAt(i))) {return false;}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/

	public long count(){return count;}
	public void setCount(long x){count=x;}
	public int length() {return name.length();}
	public char charAt(int i) {return name.charAt(i);}
	public byte[] getBytes() {return name.getBytes();}
	
	public int length1() {
		for(int i=0; i<name.length(); i++) {
			if(!Tools.isLetter(name.charAt(i))) {return i;}
		}
		return name.length();
	}
	public int length2() {
		int i=name.length()-1;
		while(i>=0 && Tools.isDigit(name.charAt(i))){
			i--;
		}
		int len=0;
		while(i>=0 && Tools.isLetter(name.charAt(i))){
			i--;
			len++;
		}
		return i<0 ? 0 : len;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void increment() {count++;}
	public void increment(Barcode b) {count+=b.count;}
	public void increment(long x) {count+=x;}
	public synchronized void incrementSync(long x) {count+=x;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Overrides          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int hashCode(){
		return name.hashCode()^(tile*65);
	}
	
	@Override
	public int compareTo(Barcode b) {
		if(count!=b.count){return count>b.count ? -1 : 1;}
		if(tile!=b.tile){return tile>b.tile ? 1 : -1;}
		return name.compareTo(b.name);
	}
	
	@Override
	public boolean equals(Object b) {
//		return(b!=null && getClass()==b.getClass() && equals((Barcode)b));//Proper way
		return(equals((Barcode)b));//Faster way
	}
	public boolean equals(Barcode b) {return(b!=null && name.equals(b.name) && tile==b.tile);}//Ignores count
	
	@Override
	public String toString(){
		return name+"\t"+count+(expected!=1 ? "\te"+expected : "")+
				(frequency!=1 ? "\tf"+frequency : "")+(tile>0 ? "\tt"+tile : "");
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(name);
		if(tile>0) {bb.append(tile);}
		return bb.tab().append(count);
	}
	
	public ByteBuilder appendIndex(ByteBuilder bb, byte delimiter, int indexNum) {
		return appendIndex(bb, delimiter, indexNum, name);
	}
	
	public static ByteBuilder appendIndex(ByteBuilder bb, byte delimiter, int indexNum, String name) {
		for(int i=0, currentIndex=1; i<name.length() && currentIndex<=indexNum; i++) {
			final char c=name.charAt(i);
			if(c==delimiter) {currentIndex++;}
			else if(currentIndex==indexNum){bb.append(c);}
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static byte delimiter(String barcode){
		if(barcode==null || barcode.length()<3) {return 0;}
		int letters=0, nonletters=0;
		byte delimiter=0;
		for(int i=0; i<barcode.length(); i++){
			char c=barcode.charAt(i);
			if(Tools.isLetter(c)){
				letters++;
			}else{
				if(nonletters==0) {delimiter=(byte)c;}
				nonletters++;
			}
		}
		if(nonletters==1 && letters>1){
			return delimiter;
		}
		return 0;//No delimiter or multiple delimiters
	}
	
	public static byte delimiter(byte[] barcode){
		if(barcode==null || barcode.length<3) {return 0;}
		int letters=0, numbers=0, other=0;
		byte delimiter=0;
		for(int i=0; i<barcode.length; i++){
			byte c=barcode[i];
			if(Tools.isLetter(c)){
				letters++;
			}else if(Tools.isDigit(c)){
				numbers++;
			}else{
				if(other==0) {delimiter=(byte)c;}
				other++;
			}
		}
		if(other==1 && letters>1){
			return delimiter;
		}
		return 0;//No delimiter or multiple delimiters
	}
	
	public static int countPolymers(String code) {
		int delimiter=code.length();
		for(int i=0; i<code.length(); i++) {
			if(!Tools.isLetter(code.charAt(i))){delimiter=i; break;}
		}
		
		int polymers=0;
		{
			int same=0, len=0;
			char last=0;
			for(int i=0; i<delimiter; i++) {
				char c=code.charAt(i);
				same+=(c==last ? 1 : 0);
				len++;
				last=c;
			}
			polymers+=(len>1 && same>=len-1) ? 1 : 0;
		}
		{
			int same=0, len=0;
			char last=0;
			for(int i=delimiter+1; i<code.length(); i++) {
				char c=code.charAt(i);
				same+=(c==last ? 1 : 0);
				len++;
				last=c;
			}
			polymers+=(len>1 && same>=len-1) ? 1 : 0;
		}
		return polymers; //0, 1, or 2.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String name;
	private long count=0;
	public final int expected;
	public float frequency=1f;
	public int tile=0;
	
}
