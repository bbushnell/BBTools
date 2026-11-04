package dna;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.Range;

/** Mutable chromosome bases with zero-based inclusive logical bounds and spare capacity.
 * Empty instances use minIndex=Integer.MAX_VALUE and maxIndex=-1. Setters copy and
 * normalize input, grow storage when needed, and widen the logical interval; gaps
 * between writes are not initialized to N. Padded accessors return N outside the
 * logical interval but preserve bytes, including zero, inside it.
 *
 * The public backing array is owned by this object but directly accessible to callers;
 * direct mutations must maintain bounds and valid base-table indices. Numeric helpers
 * read storage directly rather than applying logical-bound padding. Instances and
 * shared normalization switches are not safe for concurrent mutation.
 */
public class ChromosomeArray implements Serializable{

	/** Serial version UID for Serializable interface */
	private static final long serialVersionUID=3199182397853127842L;

	/** Legacy serialization conversion driver: chromosome number followed by input filename. */
	public static void main(String[] args){
		translateFile(args[1], Byte.parseByte(args[0]));
	}

	/** Reads a serialized array, assigns its chromosome number, writes it and reloads it.
	 * The legacy destination name is derived from the last .fa in fname; this is not a FASTA parser. */
	private static void translateFile(String fname, int chrom){

		long time1=System.nanoTime();

		ChromosomeArray cha=read(fname, chrom);
		cha.chromosome=chrom;
		long time2=System.nanoTime();

		int dot=fname.lastIndexOf(".fa");
		String outfile=fname.substring(0, dot).replace("hs_ref_", "")+".chrom";

		System.out.println("Writing to "+outfile);

		System.out.println("minIndex="+cha.minIndex+", maxIndex="+cha.maxIndex+", length="+cha.array.length+
				"; time="+Tools.format("%.3f seconds", (time2-time1)/1000000000d));

		long time3=System.nanoTime();
		ReadWrite.write(cha, outfile, false);
		cha=null;
		System.gc();
		cha=read(outfile);
		long time4=System.nanoTime();

		System.out.println("minIndex="+cha.minIndex+", maxIndex="+cha.maxIndex+", length="+cha.array.length+
				"; time="+Tools.format("%.3f seconds", (time4-time3)/1000000000d));
	}

	/** Deserializes an unassigned chromosome and replaces its number.
	 * @param fname Serialized ChromosomeArray input
	 * @param chrom New chromosome number
	 * @return Deserialized array, normalized if the shared read flag requests it */
	public static ChromosomeArray read(String fname, int chrom){
		ChromosomeArray cha=read(fname);
		assert(cha.chromosome<1);
		cha.chromosome=chrom;
		return cha;
	}

	/** Deserializes an array and optionally applies changeUndefinedToN to its complete storage.
	 * @param fname Serialized input path, with compression handled by ReadWrite
	 * @return Loaded mutable array */
	public static ChromosomeArray read(String fname){

//		if(fname.endsWith(".chrom") || fname.endsWith(".chrom.gz")){}
		ChromosomeArray ca=ReadWrite.read(ChromosomeArray.class, fname, false);
		if(CHANGE_UNDEFINED_TO_N_ON_READ){
			ca.changeUndefinedToN();
		}
		return ca;
	}

	/** Replaces bytes rejected by AminoAcid.isACGTN with N across the complete backing array.
	 * Accepted lowercase bases and U are preserved; logical bounds do not change. */
	public void changeUndefinedToN(){
		for(int i=0; i<array.length; i++){
//			array[i]=AminoAcid.numberToBase[AminoAcid.baseToNumberACGTother[array[i]]];
			if(!AminoAcid.isACGTN(array[i])){array[i]='N';}
		}
	}

	/** Creates an empty plus-strand array with unassigned chromosome number -1. */
	public ChromosomeArray(){
		this((byte)-1, Shared.PLUS);
	}

	/**
	 * Creates the reverse complement of this ChromosomeArray.
	 * Reverses raw storage coordinates zero through maxIndex, including positions below
	 * minIndex, using the extended complement table. The source is unchanged.
	 * @return Independent array with opposite strand, minIndex zero and the same maximum
	 */
	public ChromosomeArray complement(){
		byte otherStrand=(strand==Shared.MINUS ? Shared.PLUS : Shared.MINUS);
		ChromosomeArray ca=new ChromosomeArray(chromosome, otherStrand, 0, maxIndex);
		for(int i=0; i<=maxIndex; i++){
			int pos=maxIndex-i;
			byte b=AminoAcid.baseToComplementExtended[array[i]];
			ca.array[pos]=b;
		}
		return ca;
	}

	/** Creates empty storage with the requested chromosome and strand; initial capacity is at least 1000. */
	public ChromosomeArray(int chrom, byte strnd){
		this(chrom, strnd, Integer.MAX_VALUE, -1);
	}

	/** Allocates zero-filled storage and installs the supplied logical bounds without filling bases.
	 * @param chrom Chromosome number
	 * @param strnd Strand marker
	 * @param min Inclusive logical minimum, or Integer.MAX_VALUE for the usual empty sentinel
	 * @param max Inclusive logical maximum, or -1 for an empty array */
	public ChromosomeArray(int chrom, byte strnd, int min, int max){
		chromosome=chrom;
		strand=strnd;
		array=KillSwitch.allocByte1D(Tools.max(1000, max+1));
		minIndex=min;
		maxIndex=max;
	}

	/** Normalizes one table-supported base and widens logical bounds to include its location.
	 * @param loc Nonnegative destination coordinate
	 * @param val ASCII base code supported by the active normalization table */
	public void set(int loc, int val){

		if(loc>=array.length){//Increase size
			int newlen=(int)(1+(3L*max(array.length, loc))/2);
			assert(newlen>loc) : newlen+", "+loc+", "+array.length;
			resize(newlen);
			assert(array.length==newlen);
//			System.err.println("Resized array to "+newlen);
		}
		if(CHANGE_U_TO_T && CHANGE_DEGENERATE_TO_N){
			val=AminoAcid.baseToACGTN[val];
		}else{
			val=Tools.toUpperCase((char)val);
			if(AminoAcid.baseToNumberExtended[val]<0){val='N';}
		}
		array[loc]=(val>Byte.MAX_VALUE ? Byte.MAX_VALUE : (byte)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}

	/** Copies and normalizes the sequence starting at loc; an empty sequence changes nothing.
	 * @param loc Nonnegative first destination coordinate; the full write must fit int indexing
	 * @param s Nonnull sequence of table-supported characters, not modified by this method */
	public void set(int loc, CharSequence s){
		if(s.length()==0){return;}
		int loc2=loc+s.length();
		if(loc2>array.length){//Increase size
			int newlen=(int)(1+(3L*max(array.length, loc2))/2);
			assert(newlen>loc2) : newlen+", "+loc2+", "+array.length;
			resize(newlen);
			assert(array.length==newlen);
//			System.err.println("Resized array to "+newlen);
		}

		if(CHANGE_U_TO_T && CHANGE_DEGENERATE_TO_N){
			for(int i=0; i<s.length(); i++, loc++){
				array[loc]=AminoAcid.baseToACGTN[s.charAt(i)];
			}
		}else{
			for(int i=0; i<s.length(); i++, loc++){
				char c=Tools.toUpperCase(s.charAt(i));
				if(AminoAcid.baseToNumberExtended[c]<0){c='N';}
				assert(Tools.isLetter(c));
				assert(c<=Byte.MAX_VALUE);
				array[loc]=(byte)c;
			}
		}

		loc--;
		assert(loc==loc2-1) : "loc="+loc+", loc2="+loc2+", s.len="+s.length();
		minIndex=min(loc2-s.length(), minIndex); //first written index, not last
		maxIndex=max(loc, maxIndex);
	}

	/** Copies and normalizes the complete nonnull byte array starting at loc. */
	public void set(int loc, byte[] s){
		set(loc, s, s.length);
	}

	/** Copies and normalizes the builder's used prefix without retaining or modifying the builder. */
	public void set(int loc, ByteBuilder bb){
		set(loc, bb.array, bb.length());
	}

	/** Copies the first slen input bytes into storage; zero length leaves bounds and capacity alone.
	 * @param loc Nonnegative first destination coordinate; the write must fit int indexing
	 * @param s Nonnull source storage
	 * @param slen Prefix length from zero through s.length */
	public void set(int loc, byte[] s, final int slen){
		assert(slen<=s.length && slen>=0);
		if(slen==0){return;}
		int loc2=loc+slen;
		if(loc2>array.length){//Increase size
			int newlen=(int)(1+(3L*max(array.length, loc2))/2);
			assert(newlen>loc2) : newlen+", "+loc2+", "+array.length;
			resize(newlen);
			assert(array.length==newlen);
//			System.err.println("Resized array to "+newlen);
		}

		if(CHANGE_U_TO_T && CHANGE_DEGENERATE_TO_N){
			for(int i=0; i<slen; i++, loc++){
				byte b=(byte)Tools.max(0, s[i]);
				array[loc]=AminoAcid.baseToACGTN[b];
			}
		}else{
			for(int i=0; i<slen; i++, loc++){
				char c=Tools.max((char)0, Tools.toUpperCase((char)s[i]));
				if(AminoAcid.baseToNumberExtended[c]<0){c='N';}
				assert(Tools.isLetter(c));
				assert(c<=Byte.MAX_VALUE);
				array[loc]=(byte)c;
			}
		}
		loc--;
		assert(loc==loc2-1) : "loc="+loc+", loc2="+loc2+", s.len="+slen;
		minIndex=min(loc2-slen, minIndex); //first written index, not last
		maxIndex=max(loc, maxIndex);
	}

	/**
	 * Calculates GC content for a genomic region.
	 * @param loc Starting position
	 * @param length Region length
	 * @param counts Optional array to store base counts (reused if provided)
	 * @return GC fraction among A/C/G/T/U bases, or zero when that denominator is empty
	 */
	public float calcGC(int loc, int length, int[] counts){
		counts=countACGTINOC(loc, length, counts);
		long at=counts[0]+counts[3];
		long gc=counts[1]+counts[2];
		return gc/(float)Tools.max(at+gc, 1);
	}

	/**
	 * Counts bases in a genomic region by category.
	 * Uses padded get access within the nonnegative interval ending no later than maxIndex.
	 * @param loc Nonnegative starting position
	 * @param length Nonnegative region length; loc+length must not exceed maxIndex+1
	 * @param counts Optional eight-element array, cleared before reuse
	 * @return Counts in bins {A, C, G, T/U, Iupac/X, N, Other, Control}
	 */
	public int[] countACGTINOC(final int loc, final int length, int[] counts){
		final int lim=loc+length;
		assert(loc>=0 && lim<=maxIndex+1 && loc<=lim);
		if(counts==null){counts=new int[8];}
		else{Arrays.fill(counts, 0);}
		assert(counts.length==8);
		for(int i=loc; i<lim; i++){
			byte b=get(i);
			int num=charToNum[b];
			counts[num]++;
		}
		return counts;
	}

	/**
	 * Gets the base at the specified genomic location.
	 * Returns 'N' if loc is before minIndex or after maxIndex.
	 * @param loc Genomic coordinate
	 * @return Base as byte value, or 'N' if outside valid range
	 */
	public byte get(int loc){
		return loc<minIndex || loc>maxIndex ? (byte)'N' : array[loc]; //maxIndex is inclusive
	}

	/** Copies the inclusive interval using get, with N outside logical bounds.
	 * The requested length must be representable as an int; stored zero bytes are preserved. */
	public String getString(int a, int b){
		StringBuilder sb=new StringBuilder(b-a+1);
		for(int i=a; i<=b; i++){
			sb.append((char)get(i));
		}
		return sb.toString();
	}

	/**
	 * Extracts a genomic sequence as a byte array.
	 * Copies in-range requests in bulk and fills positions outside logical bounds with N.
	 * Returns independent storage and preserves every byte inside the logical interval.
	 * @param a Starting coordinate (inclusive)
	 * @param b Ending coordinate (inclusive)
	 * @return Byte array containing sequence from coordinates a to b
	 * @throws IllegalArgumentException If the requested length is negative or exceeds an int
	 */
	public byte[] getBytes(int a, int b){
		if(a>=minIndex && b<=maxIndex){return KillSwitch.copyOfRange(array, a, b+1);}
		final long length=(long)b-a+1;
		if(length<0 || length>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Invalid genomic interval ["+a+", "+b+"]: length="+length);
		}
		final byte[] out=KillSwitch.allocByte1D((int)length);
		Arrays.fill(out, (byte)'N');
		final int from=max(a, minIndex), to=min(b, maxIndex);
		if(from<=to){
			System.arraycopy(array, from, out, (int)((long)from-a), to-from+1);
		}
		return out;
	}

	/** Returns the ACGTN lookup code for raw storage at loc, without logical-bound padding. */
	public byte getNumberACGTN(int loc){
		return AminoAcid.baseToNumberACGTN[array[loc]];
	}

	/** Returns the A/C/G/T code (U maps to T), or -1 for an undefined stored base. */
	public byte getNumber(int loc){
		return AminoAcid.baseToNumber[array[loc]];
	}

	/** Tests raw storage in the inclusive interval for nonnegative baseToNumber codes.
	 * An empty interval returns true; no logical-bound padding is applied. */
	public boolean isFullyDefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x<0){return false;}
		}
		return true;
	}

	/** Tests raw storage in the inclusive interval for negative baseToNumber codes.
	 * An empty interval returns true; no logical-bound padding is applied. */
	public boolean isFullyUndefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){return false;}
		}
		return true;
	}

	/** Counts table-defined bases in the logical interval, or zero for the empty sentinel. */
	public int countDefinedBases(){
		return countDefinedBases(minIndex, maxIndex);
	}

	/** Counts nonnegative baseToNumber codes in the inclusive raw-storage interval. */
	public int countDefinedBases(int a, int b){
		int sum=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){sum++;}
		}
		return sum;
	}

	/** Packs an inclusive stored interval through the byte-array toNumber helper. */
	public int getNumber(int a, int b){
		return toNumber(a, b, array);
	}

	/** Packs bases left-to-right using two bits per base; an undefined base returns -1.
	 * At most 16 bases fit in the int representation; 16-base results may be negative
	 * and the all-T encoding is indistinguishable from the -1 undefined sentinel.
	 * @param a First source offset, inclusive
	 * @param b Last source offset, inclusive
	 * @param bases Source bytes indexed directly by the base-code table
	 * @return Packed bits, or -1 if any base is undefined */
	public static int toNumber(int a, int b, byte[] bases){
		assert(b>=a);
		//TODO: Probable bug (DNA-021) - this assertion permits 17 bases, whose 34 bits
		//truncate in an int. KeyRing.makeKeys requires chunksize<16; a live 17-base caller is unconfirmed.
		assert(b-a<17);
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}

	/** String variant of two-bit packing with the same sixteen-base representation limit.
	 * @return Packed bits, or -1 for an undefined base; bounds/length are not validated here */
	public static int toNumber(int a, int b, String bases){
		//TODO: Probable bug (DNA-021) - oversized intervals silently discard high packed bits.
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases.charAt(i)];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}

	/** Changes capacity without discarding the base at the inclusive maximum index.
	 * @throws IllegalArgumentException If the capacity is negative or would discard stored bases */
	public void resize(int newlen){
		if(newlen<0 || newlen<=maxIndex){
			throw new IllegalArgumentException("Capacity "+newlen+" cannot preserve inclusive maximum index "+maxIndex);
		}
		byte[] temp=KillSwitch.allocByte1D(newlen);
		int lim=min(array.length, newlen);
		assert(lim>maxIndex) : "Copy must include the stored base at maxIndex: copied="+lim+", maxIndex="+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}

	/** Decodes the entire backing capacity using the platform charset, including unused storage. */
	public String toBaseString(){
		String s=new String(array);
		return s;
	}

	/** Returns saturated distances to the nearest uppercase A/C/G/T across all allocated storage.
	 * Requires nonempty storage. Character.MAX_VALUE denotes no nearer defined base;
	 * logical bounds and lowercase/other table-defined symbols are not consulted. */
	public char[] nearestDefinedBase(){
		char[] r=new char[array.length];
		final char max=Character.MAX_VALUE;

		char dist=max;
		for(int i=0; i<r.length; i++){
			byte b=array[i];
			if(b=='A' || b=='C' || b=='G' || b=='T'){
				dist=0;
			}else{
				dist=(dist==max ? max : (char)(dist+1));
			}
			r[i]=dist;
		}

		dist=r[r.length-1];
		for(int i=r.length-1; i>=0; i--){
			byte b=array[i];
			if(b=='A' || b=='C' || b=='G' || b=='T'){
				dist=0;
			}else{
				dist=(dist==max ? max : (char)(dist+1));
			}
			r[i]=Tools.min(dist, r[i]);
		}
		return r;
	}

	/** Finds inclusive contig spans in the logical interval, dropping terminal N/X runs.
	 * Uppercase X always splits; uppercase N splits at nBlockSize, while shorter gaps remain.
	 * Other stored symbols are treated as contig bases by this method.
	 * @param nBlockSize Positive N-run splitting threshold
	 * @return Newly allocated ordered ranges, empty when no contig was found */
	public ArrayList<Range> toContigRanges(final int nBlockSize){
		assert(nBlockSize>0);
		ArrayList<Range> list=new ArrayList<Range>();

		//State machine: 'contig' is open while inside defined bases; stop tracks the last
		//defined base, so emitted Range(start,stop) is inclusive and excludes trailing Ns.
		//An 'X' breaks immediately; an 'N'-run breaks only once it reaches nBlockSize, so
		//gaps shorter than nBlockSize are bridged within one contig. ns seeds to nBlockSize+1
		//(as if preceded by a full gap) so a leading defined base opens cleanly.
		int start=-1;
		int stop=-1;
		int ns=nBlockSize+1;

		boolean contig=false;

		for(int i=minIndex; i<=maxIndex; i++){
			byte b=array[i];
			if(b=='N' || b=='X'){
				ns++;
				if(contig && (b=='X' || ns>=nBlockSize)){
					Range r=new Range(start, stop);
					list.add(r);
					contig=false;
				}
			}else{
				ns=0;
				if(!contig){start=i;}
				contig=true;
				stop=i;
			}
		}
		if(contig){
			Range r=new Range(start, stop);
			list.add(r);
		}
		return list;
	}

	/** Compares chromosome number, logical bounds, capacity and case-folded logical bases.
	 * Ignores strand and bytes outside logical bounds. A mismatch prints a legacy
	 * diagnostic letter to stderr; this does not override Object.equals.
	 * @param other Nonnull array to compare
	 * @return Whether the inspected fields and bases match */
	public boolean equalsIgnoreCase(ChromosomeArray other){
		if(minIndex!=other.minIndex){System.err.println("a");return false;}
		if(maxIndex!=other.maxIndex){System.err.println("b");return false;}
		if(chromosome!=other.chromosome){System.err.println("c");return false;}
		if(array.length!=other.array.length){System.err.println("d");return false;}
		for(int i=minIndex; i<=maxIndex; i++){
			if(Tools.toLowerCase(array[i])!=Tools.toLowerCase(other.array[i])){
				System.err.println("e");
				return false;
			}
		}
		return true;
	}

	/** Returns the smaller long value. */
	private static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger long value. */
	private static final long max(long x, long y){return x>y ? x : y;}
	/** Returns the smaller int value. */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger int value. */
	private static final int max(int x, int y){return x>y ? x : y;}

	/** Strand marker fixed at construction; complement creates an opposite-strand object. */
	public final byte strand;
	/** Mutable chromosome identifier; -1 is the usual unassigned value. */
	public int chromosome;
	/** Owned storage, directly exposed for callers that maintain bounds and base-table invariants. */
	public byte[] array;
	/** Inclusive logical maximum coordinate; -1 for the usual empty sentinel. */
	public int maxIndex=-1;
	/** Inclusive logical minimum coordinate; Integer.MAX_VALUE for the usual empty sentinel. */
	public int minIndex=Integer.MAX_VALUE;

	/** Enables complete-storage normalization after deserialization. */
	public static boolean CHANGE_UNDEFINED_TO_N_ON_READ=false;
	/** Together with CHANGE_DEGENERATE_TO_N selects ACGTN normalization, including U-to-T. */
	public static boolean CHANGE_U_TO_T=true;
	/** Together with CHANGE_U_TO_T selects ACGTN normalization; otherwise the extended path is used. */
	public static boolean CHANGE_DEGENERATE_TO_N=true;

	/** Translation array for converting characters to numeric base categories */
	private static final byte[] charToNum=AssemblyStats2.makeCharToNum();

}
