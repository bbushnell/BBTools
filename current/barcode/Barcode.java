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

/**
 * Represents a DNA sequencing barcode with comprehensive analysis and comparison capabilities.
 * Supports dual-indexed barcodes, distance calculations, validation, and count tracking.
 * Used for barcode demultiplexing, quality control, and sequencing data analysis.
 *
 * @author Brian Bushnell
 * @date December 30, 2013
 */
public class Barcode implements Comparable<Barcode> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a barcode with default count and expected values.
	 * @param s The barcode sequence string */
	public Barcode(String s){this(s, 0, 1, 0);}
	/**
	 * Creates a barcode with specified count and default expected value.
	 * @param s The barcode sequence string
	 * @param c The barcode count
	 */
	public Barcode(String s, long c){this(s, c, 1, 0);}
	/**
	 * Creates a barcode with specified count and expected status.
	 * @param s The barcode sequence string
	 * @param c The barcode count
	 * @param e Expected status (0 or 1)
	 */
	public Barcode(String s, long c, int e){this(s, c, e, 0);}
	/**
	 * Creates a barcode with all parameters specified.
	 *
	 * @param s The barcode sequence string
	 * @param c The barcode count
	 * @param e Expected status (0 or 1)
	 * @param t Tile number for spatial tracking
	 */
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
	/**
	 * Extracts one component string from a dual-indexed barcode.
	 * Uses delimiter if present, otherwise splits at midpoint.
	 *
	 * @param pairnum Which component to extract (0=left/first, 1=right/second)
	 * @param delimiter Character separating barcode components
	 * @return String containing the specified barcode component
	 */
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

	/**
	 * Gets the left/first component of a dual-indexed barcode.
	 * @param delimiter Character separating barcode components
	 * @return Barcode object containing left component
	 */
	public Barcode left(char delimiter){return getBarcodeForPairnum(0, delimiter);}
	/**
	 * Gets the right/second component of a dual-indexed barcode.
	 * @param delimiter Character separating barcode components
	 * @return Barcode object containing right component
	 */
	public Barcode right(char delimiter){return getBarcodeForPairnum(1, delimiter);}
	/**
	 * Gets the left/first component string of a dual-indexed barcode.
	 * @param delimiter Character separating barcode components
	 * @return String containing left component
	 */
	public String leftString(char delimiter){return getStringForPairnum(0, delimiter);}
	/**
	 * Gets the right/second component string of a dual-indexed barcode.
	 * @param delimiter Character separating barcode components
	 * @return String containing right component
	 */
	public String rightString(char delimiter){return getStringForPairnum(1, delimiter);}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Counts ambiguous or undefined bases in the barcode sequence.
	 * Identifies bases that are letters but not fully defined (e.g., N, Y, R).
	 * @return Number of undefined bases in the sequence
	 */
	public int countUndefined() {
		int sum=0;
		for(int i=0; i<length(); i++){
			char c=charAt(i);
			sum+=(Tools.isLetter(c) && !AminoAcid.isFullyDefined(c)) ? 1 : 0;
		}
		return sum;
	}
	
	/**
	 * Checks if the barcode is a homopolymer sequence.
	 * Uses the first character as the expected homopolymer base.
	 * @return true if all bases match the first character
	 */
	public boolean isHomopolymer() {return isHomopolymer(name, name.charAt(0));}
	/**
	 * Checks if the barcode is a homopolymer of the specified base.
	 * @param h Expected homopolymer base character
	 * @return true if all bases match the specified character
	 */
	public boolean isHomopolymer(char h) {return isHomopolymer(name);}
	/**
	 * Checks if a sequence is a homopolymer.
	 * Uses the first character as the expected homopolymer base.
	 * @param name Sequence to check
	 * @return true if all bases match the first character
	 */
	public static boolean isHomopolymer(String name) {return isHomopolymer(name, name.charAt(0));}
	/**
	 * Checks if a sequence is a homopolymer of the specified base.
	 * Non-letter characters are ignored in the comparison.
	 *
	 * @param name Sequence to check
	 * @param h Expected homopolymer base character
	 * @return true if all letter bases match the specified character
	 */
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

	/**
	 * Checks if both characters are non-letters.
	 * @param a First character
	 * @param b Second character
	 * @return true if both are non-letters
	 */
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
	
	/**
	 * Converts a collection of barcodes to a name-indexed HashMap.
	 * @param codes Collection of barcodes to map
	 * @return HashMap with barcode names as keys
	 */
	public static HashMap<String, Barcode> barcodesToMap(Collection<Barcode> codes) {
		HashMap<String, Barcode> countMap=new HashMap<String, Barcode>();
		for(Barcode b : codes) {countMap.put(b.name, b);}
		return countMap;
	}
	
	/**
	 * Aggregates barcode counts based on assignment mapping.
	 * Sums counts from observed barcodes to their assigned expected barcodes.
	 *
	 * @param assignmentMap Mapping from observed to assigned barcode names
	 * @param expectedCodeList List of expected barcodes to aggregate to
	 * @param countMap Counts for observed barcodes
	 * @return List of expected barcodes with aggregated counts
	 */
	public static ArrayList<Barcode> summateAssignments(HashMap<String, String> assignmentMap,
			ArrayList<Barcode> expectedCodeList, HashMap<String, Barcode> countMap) {
		ArrayList<String> list=new ArrayList<String>(expectedCodeList.size());
		for(Barcode b : expectedCodeList) {
			if(b.expected==1) {list.add(b.name);}
		}
		return summateAssignments(assignmentMap, list, countMap);
	}
	
	/**
	 * Aggregates barcode counts based on assignment mapping.
	 * Converts collection inputs to required format and delegates.
	 *
	 * @param assignmentMap Mapping from observed to assigned barcode names
	 * @param expectedCodeList Expected barcode names to aggregate to
	 * @param counts Observed barcode counts
	 * @return List of expected barcodes with aggregated counts
	 */
	public static ArrayList<Barcode> summateAssignments(HashMap<String, String> assignmentMap,
			Collection<String> expectedCodeList, Collection<Barcode> counts) {
		HashMap<String, Barcode> countMap=barcodesToMap(counts);
		return summateAssignments(assignmentMap, expectedCodeList, countMap);
	}
	
	/**
	 * Aggregates barcode counts based on assignment mapping.
	 * Creates new barcodes for expected codes and sums assigned counts.
	 *
	 * @param assignmentMap Mapping from observed to assigned barcode names
	 * @param expectedCodeList Expected barcode names to aggregate to
	 * @param countMap Counts for observed barcodes
	 * @return Sorted list of expected barcodes with aggregated counts
	 */
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
	/**
	 * Extracts barcode sequence from filename.
	 * Splits filename by dots and identifies barcode-like segments.
	 * @param fname Filename to parse
	 * @return Barcode sequence string, or null if not found
	 */
	public static String parseBarcodeFromFname(String fname) {
		String[] split=Tools.dotPattern.split(fname);
		for(String s : split) {
			if(isBarcode(s)) {return s;}
			else if("UNKNOWN".equalsIgnoreCase(s)) {return s;}
		}
//		assert(false) : "Can't find barcode in filename "+fname+"\n"+Arrays.toString(split);
		return null;
	}

	/**
	 * Checks if character is valid barcode symbol.
	 * Valid symbols are ACGTN plus delimiter characters - and +.
	 * @param x Character to check
	 * @return true if valid barcode symbol
	 */
	public static final boolean isBarcodeSymbol(char x) {return AminoAcid.isACGTN(x) || x=='-' || x=='+';}
	/**
	 * Checks if byte is valid barcode symbol.
	 * Valid symbols are ACGTN plus delimiter characters - and +.
	 * @param x Byte to check
	 * @return true if valid barcode symbol
	 */
	public static final boolean isBarcodeSymbol(byte x) {return AminoAcid.isACGTN(x) || x=='-' || x=='+';}
	
	/**
	 * Determines if string represents a valid barcode sequence.
	 * Requires at least 6 bases and at most 1 delimiter.
	 * @param s String to validate
	 * @return true if valid barcode format
	 */
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
	
	/**
	 * Checks if string contains only valid barcode symbols.
	 * @param s String to validate
	 * @return true if all characters are valid barcode symbols
	 */
	public static final boolean containsOnlyBarcodeSymbols(String s) {
		for(int i=0; i<s.length(); i++) {
			if(!isBarcodeSymbol(s.charAt(i))) {return false;}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/

	/** Gets the count value for this barcode */
	public long count(){return count;}
	/** Sets the count value for this barcode.
	 * @param x New count value */
	public void setCount(long x){count=x;}
	/** Gets the length of the barcode sequence */
	public int length() {return name.length();}
	/**
	 * Gets character at specified position in barcode sequence.
	 * @param i Position index
	 * @return Character at position i
	 */
	public char charAt(int i) {return name.charAt(i);}
	/** Gets byte array representation of barcode sequence */
	public byte[] getBytes() {return name.getBytes();}
	
	/**
	 * Gets length of first component (letters from start).
	 * Counts letters until first non-letter character.
	 * @return Length of first letter sequence
	 */
	public int length1() {
		for(int i=0; i<name.length(); i++) {
			if(!Tools.isLetter(name.charAt(i))) {return i;}
		}
		return name.length();
	}
	/**
	 * Gets length of second component (letters from end).
	 * Skips trailing digits, then counts letters backwards.
	 * @return Length of second letter sequence
	 */
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
	
	/** Increments the count by 1 */
	public void increment() {count++;}
	/** Increments the count by another barcode's count.
	 * @param b Barcode whose count to add */
	public void increment(Barcode b) {count+=b.count;}
	/** Increments the count by specified amount.
	 * @param x Amount to add to count */
	public void increment(long x) {count+=x;}
	/** Thread-safe increment of count by specified amount.
	 * @param x Amount to add to count */
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
	/**
	 * Tests equality with another barcode.
	 * Based on name and tile, ignores count.
	 * @param b Barcode to compare
	 * @return true if name and tile match
	 */
	public boolean equals(Barcode b) {return(b!=null && name.equals(b.name) && tile==b.tile);}//Ignores count
	
	@Override
	public String toString(){
		return name+"\t"+count+(expected!=1 ? "\te"+expected : "")+
				(frequency!=1 ? "\tf"+frequency : "")+(tile>0 ? "\tt"+tile : "");
	}
	
	/**
	 * Appends barcode representation to ByteBuilder.
	 * Includes name with optional tile number, then tab and count.
	 * @param bb ByteBuilder to append to
	 * @return The same ByteBuilder for chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(name);
		if(tile>0) {bb.append(tile);}
		return bb.tab().append(count);
	}
	
	/**
	 * Appends specific index component to ByteBuilder.
	 *
	 * @param bb ByteBuilder to append to
	 * @param delimiter Character separating index components
	 * @param indexNum Which index component to extract (1-based)
	 * @return The same ByteBuilder for chaining
	 */
	public ByteBuilder appendIndex(ByteBuilder bb, byte delimiter, int indexNum) {
		return appendIndex(bb, delimiter, indexNum, name);
	}
	
	/**
	 * Appends specific index component from name to ByteBuilder.
	 * Extracts the indexNum-th component separated by delimiter.
	 *
	 * @param bb ByteBuilder to append to
	 * @param delimiter Character separating index components
	 * @param indexNum Which index component to extract (1-based)
	 * @param name Barcode name to parse
	 * @return The same ByteBuilder for chaining
	 */
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
	
	/**
	 * Identifies delimiter character in barcode string.
	 * Returns delimiter only if exactly one non-letter present with multiple letters.
	 * @param barcode Barcode string to analyze
	 * @return Delimiter byte, or 0 if none or multiple delimiters
	 */
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
	
	/**
	 * Identifies delimiter character in barcode byte array.
	 * Returns delimiter only if exactly one non-letter/non-digit present.
	 * @param barcode Barcode byte array to analyze
	 * @return Delimiter byte, or 0 if none or multiple delimiters
	 */
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
	
	/**
	 * Counts homopolymer segments in dual-indexed barcode.
	 * Checks left and right components separately for homopolymer patterns.
	 * @param code Barcode string to analyze
	 * @return Number of homopolymer segments (0, 1, or 2)
	 */
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
	
	/** The barcode sequence string */
	public final String name;
	/** Count of observations for this barcode */
	private long count=0;
	/** Expected status flag (0 or 1) */
	public final int expected;
	/** Frequency value for this barcode */
	public float frequency=1f;
	/** Tile number for spatial tracking */
	public int tile=0;
	
}
