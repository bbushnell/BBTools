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


public class ChromosomeArray implements Serializable {
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 3199182397853127842L;

	public static void main(String[] args){
		translateFile(args[1], Byte.parseByte(args[0]));
	}
	
	
	private static void translateFile(String fname, int chrom){
		
		long time1=System.nanoTime();
		
		ChromosomeArray cha=read(fname, chrom);
		cha.chromosome=chrom;
		long time2=System.nanoTime();
		
		int dot=fname.lastIndexOf(".fa");
		String outfile=fname.substring(0,dot).replace("hs_ref_", "")+".chrom";
		
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
	
	public static ChromosomeArray read(String fname, int chrom){
		ChromosomeArray cha=read(fname);
		assert(cha.chromosome<1);
		cha.chromosome=chrom;
		return cha;
	}
	
	public static ChromosomeArray read(String fname){
		
//		if(fname.endsWith(".chrom") || fname.endsWith(".chrom.gz")){}
		ChromosomeArray ca=ReadWrite.read(ChromosomeArray.class, fname, false);
		if(CHANGE_UNDEFINED_TO_N_ON_READ){
			ca.changeUndefinedToN();
		}
		return ca;
	}
	
	public void changeUndefinedToN(){
		for(int i=0; i<array.length; i++){
//			array[i]=AminoAcid.numberToBase[AminoAcid.baseToNumberACGTother[array[i]]];
			if(!AminoAcid.isACGTN(array[i])){array[i]='N';}
		}
	}
	
	public ChromosomeArray(){
		this((byte)-1, Shared.PLUS);
	}
	
	/** Actually does reverse complement */
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
	
	public ChromosomeArray(int chrom, byte strnd){
		this(chrom, strnd, Integer.MAX_VALUE, -1);
	}
	
	public ChromosomeArray(int chrom, byte strnd, int min, int max){
		chromosome=chrom;
		strand=strnd;
		array=KillSwitch.allocByte1D(Tools.max(1000, max+1));
		minIndex=min;
		maxIndex=max;
	}
	
	
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
	
	
	public void set(int loc, CharSequence s){
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
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	public void set(int loc, byte[] s){
		set(loc, s, s.length);
	}
	
	public void set(int loc, ByteBuilder bb){
		set(loc, bb.array, bb.length());
	}
	
	public void set(int loc, byte[] s, final int slen){
		assert(slen<=s.length && slen>=0);
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
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}

	/**
	 * @param loc
	 * @param length
	 * @param counts
	 * @return gc fraction
	 */
	public float calcGC(int loc, int length, int[] counts) {
		counts=countACGTINOC(loc, length, counts);
		long at=counts[0]+counts[3];
		long gc=counts[1]+counts[2];
		return gc/(float)Tools.max(at+gc, 1);
	}

	/**
	 * @param loc
	 * @param length
	 * @return counts: {A, C, G, T, Iupac, N, Other, Control}
	 */
	public int[] countACGTINOC(final int loc, final int length, int[] counts) {
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
	
	
	/** Returns the letter (IUPAC) representation of the base, as a byte */
	public byte get(int loc){
		return loc<minIndex || loc>=maxIndex ? (byte)'N' : array[loc];
	}
	
	public String getString(int a, int b){
		StringBuilder sb=new StringBuilder(b-a+1);
		for(int i=a; i<=b; i++){
			sb.append((char)get(i));
		}
		return sb.toString();
	}
	
	/** Returns FASTA format bytes.  Same as getString, but faster. */
	public byte[] getBytes(int a, int b){
		byte[] out=KillSwitch.copyOfRange(array, a, b+1);
//		assert(out[0]>0 && out[out.length-1]>0) : a+", "+b+", "+minIndex+", "+maxIndex+", "+array.length;
		if(a<minIndex || b>maxIndex){
			for(int i=0; i<out.length; i++){
				if(out[i]==0){out[i]='N';}
			}
		}
		return out;
	}
	
	public byte getNumberACGTN(int loc){
		return AminoAcid.baseToNumberACGTN[array[loc]];
	}
	
	public byte getNumber(int loc){
		return AminoAcid.baseToNumber[array[loc]];
	}
	
	public boolean isFullyDefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x<0){return false;}
		}
		return true;
	}
	
	public boolean isFullyUndefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){return false;}
		}
		return true;
	}
	
	public int countDefinedBases(){
		return countDefinedBases(minIndex, maxIndex);
	}
	
	public int countDefinedBases(int a, int b){
		int sum=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){sum++;}
		}
		return sum;
	}
	
	public int getNumber(int a, int b){
		return toNumber(a, b, array);
	}
	
	public static int toNumber(int a, int b, byte[] bases){
		assert(b>=a);
		assert(b-a<17); //<17 for unsigned, <16 for signed
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}
	
	public static int toNumber(int a, int b, String bases){
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases.charAt(i)];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}
	
	public void resize(int newlen){
		byte[] temp=KillSwitch.allocByte1D(newlen);
		int lim=min(array.length, newlen);
		assert(lim>=maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}
	
	public String toBaseString(){
		String s=new String(array);
		return s;
	}
	
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
	
	public ArrayList<Range> toContigRanges(final int nBlockSize){
		assert(nBlockSize>0);
		ArrayList<Range> list=new ArrayList<Range>();
		
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
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public final byte strand;
	public int chromosome;
	public byte[] array;
	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	
	public static boolean CHANGE_UNDEFINED_TO_N_ON_READ=false;
	public static boolean CHANGE_U_TO_T=true;
	public static boolean CHANGE_DEGENERATE_TO_N=true;
	
	/** Translation array for tracking base counts */
	private static final byte[] charToNum=AssemblyStats2.makeCharToNum();
	
	
}
