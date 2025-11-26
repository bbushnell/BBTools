package structures;

import java.io.Serializable;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;

import bin.Cluster;
import bin.Contig;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.LineParser;
import shared.LineParser1;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.FASTQ;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Oct 8, 2013
 *
 */
public final class ByteBuilder implements Serializable, CharSequence {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -4786450129730831665L;
	
	/** Constructs ByteBuilder with initial capacity of 32 bytes */
	public ByteBuilder(){
		array=KillSwitch.allocByte1D(32);
	}
	
	/** Constructs ByteBuilder with specified initial capacity.
	 * @param initial Initial capacity in bytes (must be >= 1) */
	public ByteBuilder(int initial){
		assert(initial>=1) : initial;
		array=KillSwitch.allocByte1D(initial);
	}
	
	/** Constructs ByteBuilder initialized with string representation of object.
	 * @param o Object to convert to string and append */
	public ByteBuilder(Object o){
		String s=o.toString();
		array=KillSwitch.allocByte1D(s.length()+1);
		append(s);
	}
	
	/**
	 * Constructs ByteBuilder wrapping existing byte array.
	 * Length is set to full array length.
	 * @param array_ Byte array to wrap (not copied)
	 */
	public ByteBuilder(byte[] array_){
		assert(array_!=null);
		array=array_;
		length=array.length;
	}
	
	/** Copy constructor that creates independent copy of another ByteBuilder.
	 * @param bb ByteBuilder to copy */
	public ByteBuilder(ByteBuilder bb){
		array=bb.toBytes();
		length=bb.length();
	}

	@Override
	public CharSequence subSequence(int start, int end) throws IndexOutOfBoundsException{
		if(start<0 || end>=length()){throw new IndexOutOfBoundsException();}
		return new ByteBuilder(KillSwitch.copyOfRange(array, start, end));
	}
	
	/**
	 * Finds first occurrence of character.
	 * @param c Character to search for
	 * @return Index of first occurrence, or -1 if not found
	 */
	public int indexOf(char c) {return indexOf((byte)c);}
	/**
	 * Finds last occurrence of character.
	 * @param c Character to search for
	 * @return Index of last occurrence, or -1 if not found
	 */
	public int lastIndexOf(char c) {return lastIndexOf((byte)c);}
	/**
	 * Counts occurrences of character.
	 * @param c Character to count
	 * @return Number of occurrences
	 */
	public int count(char c) {return count((byte)c);}
	
	/**
	 * Finds first occurrence of byte value.
	 * @param c Byte value to search for
	 * @return Index of first occurrence, or -1 if not found
	 */
	public int indexOf(byte c) {
		for(int i=0; i<length; i++) {
			if(array[i]==c) {return i;}
		}
		return -1;
	}
	
	/**
	 * Finds last occurrence of byte value.
	 * @param c Byte value to search for
	 * @return Index of last occurrence, or -1 if not found
	 */
	public int lastIndexOf(byte c) {
		for(int i=length-1; i>=0; i--) {
			if(array[i]==c) {return i;}
		}
		return -1;
	}

	/**
	 * Counts occurrences of byte value.
	 * @param c Byte value to count
	 * @return Number of occurrences
	 */
	public int count(byte c) {
		int count=0;
		for(int i=0; i<length; i++) {
			count+=(array[i]==c) ? 1 : 0;
		}
		return count;
	}

//	public ByteBuilder append(float x, int places){return append(Tools.format(decimalFormat[places], x));}
//	public ByteBuilder append(double x, int places){return append(Tools.format(decimalFormat[places], x));}

	/**
	 * Appends float using Float.toString() (slower but more precise).
	 * @param x Float value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendSlow(float x){return append(Float.toString(x));}
	/**
	 * Appends double using Double.toString() (slower but more precise).
	 * @param x Double value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendSlow(double x){return append(Double.toString(x));}
	/**
	 * Appends float with specified decimal places using String.format().
	 * @param x Float value to append
	 * @param decimals Number of decimal places
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendSlow(float x, int decimals){return append(Tools.format("%."+decimals+"f", x));}
	/**
	 * Appends double with specified decimal places using String.format().
	 * @param x Double value to append
	 * @param decimals Number of decimal places
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendSlow(double x, int decimals){return append(Tools.format("%."+decimals+"f", x));}
	/**
	 * Appends boolean value as "true" or "false".
	 * @param x Boolean value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(boolean x){return append(x ? tbool : fbool);}

	/**
	 * Appends two characters.
	 * @param a First character
	 * @param b Second character
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char a, char b) {
		append(a);append(b);
		return this;
	}
	/**
	 * Appends two bytes.
	 * @param a First byte
	 * @param b Second byte
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte a, byte b) {
		append(a);append(b);
		return this;
	}
	/**
	 * Appends three characters.
	 * @param a First character
	 * @param b Second character
	 * @param c Third character
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char a, char b, char c) {
		append(a);append(b);append(c);
		return this;
	}
	/**
	 * Appends three bytes.
	 * @param a First byte
	 * @param b Second byte
	 * @param c Third byte
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte a, byte b, byte c) {
		append(a);append(b);append(c);
		return this;
	}
	/**
	 * Appends four characters.
	 * @param a First character
	 * @param b Second character
	 * @param c Third character
	 * @param d Fourth character
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char a, char b, char c, char d) {
		append(a);append(b);append(c);append(d);
		return this;
	}
	/**
	 * Appends four bytes.
	 * @param a First byte
	 * @param b Second byte
	 * @param c Third byte
	 * @param d Fourth byte
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte a, byte b, byte c, byte d) {
		append(a);append(b);append(c);append(d);
		return this;
	}
	/**
	 * Appends character followed by integer.
	 * @param x Character to append
	 * @param y Integer to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char x, int y){
		append(x); append(y); return this;
	}
	/**
	 * Appends character followed by long integer.
	 * @param x Character to append
	 * @param y Long integer to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char x, long y){
		append(x); append(y); return this;
	}
	
	/**
	 * Appends single character, expanding array if needed.
	 * @param x Character to append (cast to byte)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char x){
		if(length>=array.length){expand();}
		array[length++]=(byte)x;
		return this;
	}
	/**
	 * Appends single byte, expanding array if needed.
	 * @param x Byte value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte x){
		if(length>=array.length){expand();}
		array[length++]=x;
		return this;
	}
	/** For sam/bam */
	public ByteBuilder appendColon(char x){
		if(length+1>=array.length){expand();}
		array[length++]=(byte)x;
		array[length++]=(byte)':';
		return this;
	}
	/** For sam/bam */
	public ByteBuilder appendColon(byte x){
		if(length+1>=array.length){expand();}
		array[length++]=x;
		array[length++]=(byte)':';
		return this;
	}
	
	/**
	 * Appends k-mer as nucleotide sequence.
	 * @param kmer K-mer object to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendKmer(Kmer kmer) {
		return appendKmer(kmer.array1(), kmer.k);
	}
	
	/**
	 * Appends array of k-mers as nucleotide sequences.
	 * @param kmer Array of k-mer values
	 * @param k K-mer length
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendKmer(long[] kmer, int k) {
		for(long subkmer : kmer){
			appendKmer(subkmer, k);
		}
		return this;
	}
	
	/**
	 * @param kmer
	 * @param k
	 */
	public ByteBuilder appendKmer(long kmer, int k) {
		kmer=AminoAcid.reverseComplementBinaryFast(~kmer, k);
		for(int i=0; i<k; i++){
			int x=(int)(kmer&3);
			append((char)AminoAcid.numberToBase[x]);
			kmer>>=2;
		}
		return this;
	}
	
	/**
	 * Appends specified term from LineParser.
	 * @param lp LineParser containing terms
	 * @param term Term index to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendTerm(LineParser lp, int term) {
		return lp.appendTerm(this, term);
	}
	
	/**
	 * Appends integer followed by newline.
	 * @param x Integer value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendln(int x){
		expand(12);
		append(x);
		return append(newline);
	}
	
	public ByteBuilder appendln(char x){
		expand(2);
		array[length++]=(byte)x;
		array[length++]=newline;
		return this;
	}
	
	public ByteBuilder appendln(byte x){
		expand(2);
		array[length++]=(byte)x;
		array[length++]=newline;
		return this;
	}
	
	public ByteBuilder appendln(long x){
		expand(12);
		append(x);
		return append(newline);
	}
	
	/**
	 * Appends integer using optimized formatting with lookup tables.
	 * Handles negatives and uses two-digit chunking for efficiency.
	 * @param x Integer value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(int x){
		expand(11);
		if(x<0){
			if(x<=Integer.MIN_VALUE){
				return append((long)x);
			}else{
				array[length]='-';
				length++;
				x=-x;
			}
		}else if(x==0){
			array[length]='0';
			length++;
			return this;
		}
		
		int pos=0;
		while(x>9){
			int y=x%100;
			x=x/100;
			numbuffer[pos]=ones100[y];
			pos++;
			numbuffer[pos]=tens100[y];
			pos++;
		}
		while(x>0){
			int y=x%10;
			x=x/10;
			numbuffer[pos]=ones100[y];
			pos++;
		}
		
		assert(pos>0) : pos+", "+x;
		while(pos>0){
			pos--;
			array[length]=numbuffer[pos];
			length++;
		}
		
		return this;
	}
	
	public ByteBuilder appendOpUnsafe(int x, byte b){
		if(x<0){
			if(x<=Integer.MIN_VALUE){
				return append((long)x).append(b);
			}else{
				array[length++]='-';
				x=-x;
			}
		}else if(x==0){
			array[length++]='0';
			array[length++]=b;
			return this;
		}
		
		int pos=0;
		while(x>9){
			int y=x%100;
			x=x/100;
			numbuffer[pos++]=ones100[y];
			numbuffer[pos++]=tens100[y];
		}
		while(x>0){
			int y=x%10;
			x=x/10;
			numbuffer[pos++]=ones100[y];
		}
		
		assert(pos>0) : pos+", "+x;
		while(pos>0){
			pos--;
			array[length++]=numbuffer[pos];
		}
		array[length++]=b;
		return this;
	}
	
	/**
	 * Appends long integer using optimized formatting with lookup tables.
	 * Uses int path for values that fit, otherwise uses two-digit chunking.
	 * @param x Long value to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(long x){
		if(x>Integer.MIN_VALUE && x<=Integer.MAX_VALUE){return append((int)x);}
		expand(20);
		if(x<0){
			if(x==Long.MIN_VALUE){
				return append(Long.toString(Long.MIN_VALUE));
			}else{
				array[length]='-';
				length++;
				x=-x;
			}
		}else if(x==0){
			array[length]='0';
			length++;
			return this;
		}
		
		int pos=0;
		while(x>9){
			int y=(int)(x%100);
			x=x/100;
			numbuffer[pos]=ones100[y];
			pos++;
			numbuffer[pos]=tens100[y];
			pos++;
		}
		while(x>0){
			int y=(int)(x%10);
			x=x/10;
			numbuffer[pos]=ones100[y];
			pos++;
		}
		
		assert(pos>0) : x;
		while(pos>0){
			pos--;
			array[length]=numbuffer[pos];
			length++;
		}
		
		return this;
	}
	
	/**
	 * Appends double with specified decimal places using fast formatting.
	 * @param x0 Double value to append
	 * @param decimals0 Number of decimal places
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(final double x0, final int decimals0){
		return append(x0, decimals0, false);
	}
	
	public ByteBuilder append(final double x0, final int decimals0, boolean concise){
//		if(true){return append(x0, decimals0);}
		if(x0==(long)x0) {return append((long)x0);}
		if(decimals0<1){return append((long)(x0+0.5));}
		expand(21+decimals0);
		double x=x0;
		int decimals=decimals0;
		if(x<0){
			array[length]='-';
			length++;
			x=-x;
		}
		x=x+(0.5*decimalInvMult[decimals]);
		long upper=(long)x;
		long lower=(long)((x-upper)*decimalMult[decimals]);
		x*=decimalMult[decimals];
		x=x+0.5;
//		long longRep=(long)x;
//		assert(longRep==(long)(decimalMult[decimals]*Double.parseDouble((Tools.format("%."+decimals0+"f", x0))))) : 
//			"\n"+longRep+"\n"+decimalMult[decimals]*Double.parseDouble((Tools.format("%."+decimals0+"f", x0)))+"\n"+x0;
		
//		long upper=longRep/longMult[decimals];
//		long lower=longRep%longMult[decimals];
//		
//		append(upper);
		
		int pos=0;
		
		{
			//Lower digits
			for(; decimals>1; decimals-=2){
				int y=(int)(lower%100);
				lower=lower/100;
				numbuffer[pos]=ones100[y];
				pos++;
				numbuffer[pos]=tens100[y];
				pos++;
			}
			for(; decimals>0; decimals--){
				int y=(int)(lower%10);
				lower=lower/10;
				numbuffer[pos]=ones100[y];
				pos++;
			}
			numbuffer[pos]='.';
			pos++;

			//Upper digits
			if(upper==0){
				numbuffer[pos]='0';
				pos++;
			}else{
				while(upper>9){
					int y=(int)(upper%100);
					upper=upper/100;
					numbuffer[pos]=ones100[y];
					pos++;
					numbuffer[pos]=tens100[y];
					pos++;
				}
				while(upper>0){
					int y=(int)(upper%10);
					upper=upper/10;
					numbuffer[pos]=ones100[y];
					pos++;
				}
			}
		}
		
//		assert(Tools.format("%."+decimals0+"f", x0).equals(new String(Tools.reverseAndCopy(KillSwitch.copyOf(numbuffer, pos))))) : 
//			Tools.format("%."+decimals0+"f", x0)+"\n"+new String(Tools.reverseAndCopy(KillSwitch.copyOf(numbuffer, pos)));

//		longRep=(long)x;
//		assert(longRep==(long)(decimalMult[decimals0]*Double.parseDouble((Tools.format("%."+decimals0+"f", x0))))) : 
//			"\n"+longRep+"\n"+decimalMult[decimals0]*Double.parseDouble((Tools.format("%."+decimals0+"f", x0)))+
//			"\n"+x0+"\n"+Tools.format("%."+decimals0+"f", x0)+"\n"+new String(Tools.reverseAndCopy(KillSwitch.copyOf(numbuffer, pos)));
		
		int lim=0;
		if(concise) {//eliminates trailing zeros.  This would also be a good place to limit sig figs.
			while(concise && numbuffer[lim]=='0' && lim<decimals0-1) {lim++;}
		}
		
		while(pos>lim){
			pos--;
			array[length]=numbuffer[pos];
			length++;
		}
		
		return this;
	}
	
	public ByteBuilder append(String s){//Baseline
		if(s==null){return append(nullBytes);}
		final int len=s.length();
		expand(len);
		for(int i=0; i<len; i++, length++){
			array[length]=(byte)s.charAt(i);
		}
		return this;
	}
	
	//1.5x faster on cluster microbenchmark
	//4.5x faster on laptop micorbenchmark
	public ByteBuilder appendByCopy(String s){
		if(s==null){return append(nullBytes);}
		byte[] x=s.getBytes(StandardCharsets.ISO_8859_1);		
		expand(x.length);
		System.arraycopy(x, 0, array, length, x.length);
		length+=x.length;
		return this;
	}
	
	//11x faster on cluster and laptop microbenchmark
	//Not faster in practice though, with fastq or sam
	public ByteBuilder appendByVH(String s) {
		if(s==null){return append(nullBytes);}
		if(Vector.varHandles) {return Vector.append(this, s);}
//		else{return appendByCopy(s);}
		else{return append(s);}
	}
	
	/**
	 * Appends substring by converting characters to bytes.
	 * @param x String to append from
	 * @param from Starting index (inclusive)
	 * @param toExclusive Ending index (exclusive)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(String x, int from, int toExclusive){
		int len=toExclusive-from;
		if(x==null || len<1){return append(nullBytes);}
		expand(len);
		for(int i=from; i<toExclusive; i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	/**
	 * Appends StringBuilder by converting characters to bytes.
	 * @param x StringBuilder to append (null becomes "null")
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(StringBuilder x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	/**
	 * Appends CharSequence by converting characters to bytes.
	 * @param x CharSequence to append (null becomes "null")
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(CharSequence x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	/**
	 * Appends CharSequence followed by newline.
	 * @param x CharSequence to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendln(CharSequence x){
		if(x==null){return appendln(nullBytes);}
		expand(x.length()+1);
		append(x);
		array[length]=newline;
		length++;
		return this;
	}
	
//	public ByteBuilder append(Object x){
//		if(x==null){return append(nullBytes);}
//		return append(x.toString());
//	}
	
	/**
	 * Appends entire byte array.
	 * @param x Byte array to append (null becomes "null")
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte[] x){
		if(x==null){x=nullBytes;}
		expand(x.length);
		System.arraycopy(x, 0, array, length, x.length);
		length+=x.length;
//		for(int i=0; i<x.length; i++){
//			array[length]=x[i];
//			length++;
//		}
		return this;
	}

	/**
	 * Appends string followed by tab character.
	 * @param x String to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(String x){return append(x).tab();}
	/**
	 * Appends integer followed by tab character.
	 * @param x Integer to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(int x){return append(x).tab();}
	/**
	 * Appends long followed by tab character.
	 * @param x Long to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(long x){return append(x).tab();}
	/**
	 * Appends byte array followed by tab character.
	 * @param x Byte array to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(byte[] x){return append(x).tab();}
	/**
	 * Appends float with decimals followed by tab character.
	 * @param x Float to append
	 * @param decimals Number of decimal places
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(float x, int decimals){return append(x, decimals).tab();}
	/**
	 * Appends double with decimals followed by tab character.
	 * @param x Double to append
	 * @param decimals Number of decimal places
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendt(double x, int decimals){return append(x, decimals).tab();}
	
	/**
	 * Appends byte array followed by newline.
	 * @param x Byte array to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendln(byte[] x){
		expand(x.length+1);
		append(x);
		array[length]=newline;
		length++;
		return this;
	}
	
	public ByteBuilder appendQuality(byte[] x){
		return appendQuality(x, FASTQ.ASCII_OFFSET_OUT);
	}
	
	/**
	 * Appends quality scores as ASCII by adding 33 to each byte.
	 * Converts Phred quality scores to FASTQ ASCII representation.
	 * @param x Quality score array (null is ignored)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendQuality(byte[] x, int offset){
		if(x==null){return this;}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=(byte)(x[i]+offset);
			length++;
		}
		return this;
	}
	
	/**
	 * Appends differential quality encoding starting from ASCII 77.
	 * Encodes difference from previous quality score plus 77.
	 * @param x Quality score array (null is ignored)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendQualityDif(byte[] x){
		if(x==null){return this;}
		expand(x.length);
		int last=0;
		for(int i=0; i<x.length; i++){
			final int q=x[i];
			assert(q<=44);
			array[length]=(byte)((q-last)+77);
			length++;
			last=q;
		}
		return this;
	}
	
	/**
	 * Appends long array with delimiter between values.
	 * @param array Array of long values
	 * @param delimiter Character to place between values
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(long[] array, char delimiter){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i]);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	public ByteBuilder appendA48(long[] array, char delimiter){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			appendA48(array[i]);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	private ByteBuilder appendA48_old(long value, byte[] temp){
		int i=0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				append((char)(temp[i]+48));
			}
		}
		return this;
	}
	
	public ByteBuilder appendA48(long value){
		if(value==0){
			return append((byte)'0');
		}
		
		int highBit=63-Long.numberOfLeadingZeros(value);
		int symbols=(highBit/6)+1;
		
		expand(symbols);
		for(int shift=(symbols-1)*6; shift>=0; shift-=6){
			byte b=(byte)((value>>shift)&0x3F);
			array[length++]=(byte)(b+48);
		}
		return this;
	}
	
	public ByteBuilder appendA48(int value){
		if(value==0){
			return append((byte)'0');
		}
		
		int highBit=31-Integer.numberOfLeadingZeros(value);
		int symbols=(highBit/6)+1;
		
		expand(symbols);
		for(int shift=(symbols-1)*6; shift>=0; shift-=6){
			byte b=(byte)((value>>shift)&0x3F);
			array[length++]=(byte)(b+48);
		}
		return this;
	}
	
	/**
	 * Appends int array with delimiter between values.
	 * @param array Array of int values
	 * @param delimiter Character to place between values
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(int[] array, char delimiter){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i]);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	/**
	 * Appends float array with specified decimals and delimiter between values.
	 * @param array Array of float values
	 * @param delimiter Character to place between values
	 * @param decimals Number of decimal places for each float
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(float[] array, char delimiter, int decimals){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i], decimals);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	/**
	 * Appends contents of another ByteBuilder.
	 * @param bb ByteBuilder to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(ByteBuilder bb){
		return append(bb.array, 0, bb.length);
	}
	
	/**
	 * Appends ByteBuilder followed by newline.
	 * @param x ByteBuilder to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendln(ByteBuilder x){
		expand(x.length+1);
		append(x);
		array[length++]=newline;
		return this;
	}
	
	/**
	 * Appends byte array up to specified length.
	 * @param x Byte array to append
	 * @param len Maximum number of bytes to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte[] x, int len){
		return append(x, 0, len);
	}
	
	/**
	 * Appends portion of byte array from start position for specified length.
	 * @param x Byte array to append from
	 * @param start Starting position in array
	 * @param len Number of bytes to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(byte[] x, int start, int len){
//		if(x==null){x=nullBytes;}
		len=Tools.min(len, x.length-start);
		expand(len);
//		final int lim=(int)Tools.min(start+(long)len, x.length);
//		for(int i=start; i<lim; i++){
//			array[length]=x[i];
//			length++;
//		}
		System.arraycopy(x, start, array, length, len);
		length+=len;
		return this;
	}
	
	/**
	 * Appends char array by converting characters to bytes.
	 * @param x Char array to append (null becomes "null")
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder append(char[] x){
		if(x==null){return append(nullBytes);}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=(byte)x[i];
			length++;
		}
		return this;
	}
	
	/**
	 * Appends char array followed by newline.
	 * @param x Char array to append
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder appendln(char[] x){
		expand(x.length+1);
		append(x);
		array[length]=newline;
		length++;
		return this;
	}
	
	public void replaceLast(char c){array[length-1]=(byte)c;}
	
	public void replaceLast(byte b) {array[length-1]=b;}
	
	/**
	 * Append a newline.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder nl(){return append(newline);}
	
	/**
	 * Append a tab.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder tab(){return append(tab);}
	
	/**
	 * Append a space.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder space(){return append(' ');}
	
	/**
	 * Append an underscore.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder under(){return append('_');}
	
	/**
	 * Append a dash (hyphen).
	 * @return This ByteBuilder.
	 */
	public ByteBuilder dash(){return append('-');}
	
	/**
	 * Append a period.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder dot(){return append('.');}
	
	/**
	 * Append a comma.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder comma(){return append(',');}
	
	/**
	 * Append a colon.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder colon(){return append(':');}
	
	/**
	 * Append a colon.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder equals(){return append('=');}
	
	/**
	 * Append a semicolon.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder semi(){return append(';');}
	
	/**
	 * Append a plus.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder plus(){return append('+');}
	
	/**
	 * Append a percent.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder percent(){return append('%');}
	
	/**
	 * Gets byte at specified index.
	 * @param i Index position
	 * @return Byte value at index
	 */
	public byte get(int i){
		assert(i<length);
		return array[i];
	}
	
	/**
	 * Sets byte at specified index.
	 * @param i Index position
	 * @param b Byte value to set
	 */
	public void set(int i, byte b){
		assert(i<length);
		array[i]=b;
	}
	
	/**
	 * Sets character at specified index (converted to byte).
	 * @param i Index position
	 * @param b Character value to set
	 */
	public void set(int i, char b){
		assert(i<length) : i+", "+b+", "+length;
		array[i]=(byte)b;
	}
	
	@Override
	public char charAt(int i){
		assert(i<length);
		return (char)array[i];
	}

	/**
	 * Tests if sequence ends with specified character.
	 * @param c Character to test
	 * @return true if sequence ends with character
	 */
	public boolean endsWith(char c) {
		return length<1 ? false : array[length-1]==c;
	}

	/**
	 * Tests if sequence starts with specified string.
	 * @param s String to test
	 * @return true if sequence starts with string
	 */
	public boolean startsWith(String s) {
		return length>=s.length() && Tools.startsWith(array, s);
	}

	/**
	 * @param left Amount to trim from the left
	 * @param right Amount to trim from the right
	 */
	public void trimByAmount(int left, int right) {
		assert(left>=0 && right>=0);
		int newlen=length-left-right;
		if(newlen==length){return;}
		length=Tools.max(newlen, 0);
		if(length<1){return;}
		for(int i=0, j=left; i<newlen; i++, j++){
			array[i]=array[j];
		}
	}
	
	@Override
	public final String toString(){
		return new String(array, 0, length, StandardCharsets.US_ASCII);
	}
	
	/** Returns copy of internal array trimmed to current length.
	 * @return New byte array containing current content */
	public final byte[] toBytes(){
		return KillSwitch.copyOf(array, length);
	}
	
	/**
	 * Returns copy of internal array slice from specified range.
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 * @return New byte array containing specified range
	 */
	public final byte[] toBytes(int from, int to){
		return KillSwitch.copyOfRange(array, from, to);
	}
	
	/**
	 * Extracts prefix and shifts remaining content with optional overlap.
	 * Returns extracted bytes while compacting remaining content.
	 * @param len Number of bytes to extract from beginning
	 * @param overlap Number of bytes to preserve from extracted portion
	 * @return Extracted byte array of specified length
	 */
	public final byte[] expelAndShift(int len, int overlap){
		assert(overlap<=len) : overlap+", "+len;
		assert(len<=length);
		byte[] expel=KillSwitch.copyOf(array, len);
		if(len>0){
			for(int i=0, j=len-overlap; j<length; i++, j++){
				array[i]=array[j];
			}
		}
		length=length-len+overlap;
		return expel;
	}

	/**
	 * Shrinks internal array to specified maximum length if larger.
	 * Preserves content up to current length or maxLen, whichever is smaller.
	 * @param maxLen Maximum desired array length
	 */
	public void shrinkTo(int maxLen) {
		assert(maxLen>=0);
		if(array.length<=maxLen){return;}
//		assert(length<=maxLen) : length+", "+array.length+", "+maxLen;
		if(length<1){
			array=KillSwitch.allocByte1D(maxLen);
		}else{
			array=KillSwitch.copyOf(array, maxLen);
			length=Tools.min(length, maxLen);
		}
	}
	
	/**
	 * Tests if internal array has room for specified additional bytes.
	 * @param x Number of additional bytes needed
	 * @return true if sufficient space available
	 */
	private final boolean isRoom(int x){
		return array.length-length>=x;
	}
	
	/** Doubles internal array size up to maximum allowed length */
	private final void expand(){
		long x=Tools.min(MAXLEN, array.length*2L);
		if(x<=array.length){
			throw new RuntimeException("Array overflow: "+x+"<="+array.length);
		}
		assert(((int)x)>array.length) : "Array overflow: "+x+"<="+array.length;
		array=KillSwitch.copyOf(array, (int)x);
	}

	/**
	 * Inserts byte at specified position, shifting existing content right.
	 * @param pos Position to insert at
	 * @param c Byte value to insert
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder insert(int pos, byte c) {
		assert(pos<=length);
		if(pos==length){return append(c);}
		
		expand(1);
		for(int i=length-1; i>=pos; i--){
			array[i+1]=array[i];
		}
		array[pos]=c;
		length++;
		return this;
	}
	
	/**
	 * Expands internal array to accommodate additional bytes if needed.
	 * Doubles size repeatedly until sufficient space available.
	 * @param extra Number of additional bytes required
	 */
	public final void expand(int extra){//Doubles length until it can accommodate this much more.
		long x=array.length;
		if(x>=length+extra){return;}
//		System.err.println("x="+array.length+", extra="+extra+", length="+length);
		while(x<=length+extra){
//			System.err.println("*\t"+x+"-"+length+"<"+extra);
			x<<=1;
		}
		x=Tools.min(MAXLEN, x);
		assert(x>0 && ((int)x)>=array.length) : "Array overflow: "+x+"<array.length";
		assert(x>array.length) : "Resizing to an non-longer array ("+array.length+"); probable array size overflow.";
		array=KillSwitch.copyOf(array, (int)x);
	}
	
	/** Reverses the order of bytes in current content.
	 * @return This ByteBuilder for chaining */
	public ByteBuilder reverse() {
		Tools.reverseInPlace(array, 0, length);
		return this;
	}
	
	/** Reverses the order of bytes in current content in-place.
	 * @return This ByteBuilder for chaining */
	public ByteBuilder reverseInPlace() {
		Tools.reverseInPlace(array, 0, length);
		return this;
	}
	
	/**
	 * Reverses bytes within specified range in-place.
	 * @param from Starting position (inclusive)
	 * @param toExclusive Ending position (exclusive)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder reverseInPlace(int from, int toExclusive) {
		Tools.reverseInPlace(array, from, toExclusive);
		return this;
	}
	
	/** Converts nucleotide bases to their complements in-place.
	 * @return This ByteBuilder for chaining */
	public ByteBuilder complementInPlace() {
		AminoAcid.complementBasesInPlace(array, length);
		return this;
	}
	
	/** Converts to reverse complement of nucleotide sequence in-place.
	 * @return This ByteBuilder for chaining */
	public ByteBuilder reverseComplementInPlace() {
		Vector.reverseComplementInPlace(array, length);
		return this;
	}
	
	/** Ensures internal array has space for specified additional bytes.
	 * @param extra Number of additional bytes that will be needed */
	public final void ensureExtra(int extra){
		if(array.length-length<extra){expand(extra);}
	}
	
	/** Tests if ByteBuilder contains no content */
	public boolean isEmpty(){return length==0;}
	@Override
	public int length(){return length;}
	/** Resets length to zero without deallocating array.
	 * @return This ByteBuilder for chaining */
	public ByteBuilder clear(){
		length=0;
		return this;
	}
	/**
	 * Removes specified number of characters from end.
	 * @param x Number of characters to remove
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder trimLast(int x){
		assert(x>=0);
		setLength(length-x);
		return this;
	}
	/**
	 * Sets content length directly (may truncate or expand).
	 * @param x New length (must be within array bounds)
	 * @return This ByteBuilder for chaining
	 */
	public ByteBuilder setLength(int x){
		assert(x>=0 && x<=array.length);
		length=x;
		return this;
	}

	public byte[][] split(char c){
		LineParser1 lp=new LineParser1(c);
		lp.set(toBytes());
		final int terms=lp.terms();
		byte[][] out=new byte[terms][];
		for(int i=0; i<terms; i++) {out[i]=lp.parseByteArray(i);}
		return out;
	}
	

	/*--------------------------------------------------------------*/
	
	/** Appends unsigned 8-bit value */
	public ByteBuilder appendU8(int value){
		array[length++]=(byte)(value&0xFF);
		return this;
	}
	
	/** Appends unsigned 16-bit little-endian value */
	public ByteBuilder appendU16LE(int value){
		array[length++]=(byte)(value&0xFF);
		array[length++]=(byte)((value>>>8)&0xFF);
		return this;
	}
	
	/** Appends signed 32-bit little-endian value */
	public ByteBuilder setI32LE(int value, int offset){
		array[offset++]=(byte)(value&0xFF);
		array[offset++]=(byte)((value>>>8)&0xFF);
		array[offset++]=(byte)((value>>>16)&0xFF);
		array[offset++]=(byte)((value>>>24)&0xFF);
		return this;
	}
	
	/** Appends signed 32-bit little-endian value */
	public ByteBuilder appendI32LE(int value){
		array[length++]=(byte)(value&0xFF);
		array[length++]=(byte)((value>>>8)&0xFF);
		array[length++]=(byte)((value>>>16)&0xFF);
		array[length++]=(byte)((value>>>24)&0xFF);
		return this;
	}
	
	/** Appends unsigned 32-bit little-endian value */
	public ByteBuilder appendU32LE(long value){
		return appendI32LE((int)(value&0xFFFFFFFFL));
	}
	
	/** Appends signed 64-bit little-endian value */
	public ByteBuilder appendI64LE(long value){
		array[length++]=(byte)(value&0xFF);
		array[length++]=(byte)((value>>>8)&0xFF);
		array[length++]=(byte)((value>>>16)&0xFF);
		array[length++]=(byte)((value>>>24)&0xFF);
		array[length++]=(byte)((value>>>32)&0xFF);
		array[length++]=(byte)((value>>>40)&0xFF);
		array[length++]=(byte)((value>>>48)&0xFF);
		array[length++]=(byte)((value>>>56)&0xFF);
		return this;
	}
	
	/** Appends 32-bit little-endian float */
	public ByteBuilder appendFloatLE(float value){
		return appendI32LE(Float.floatToRawIntBits(value));
	}
	
	/** Appends 64-bit little-endian double */
	public ByteBuilder appendDoubleLE(double value){
		return appendI64LE(Double.doubleToRawLongBits(value));
	}
	
	/** Appends null-terminated string */
	public ByteBuilder appendNullTerminated(String s){
		byte[] bytes=s.getBytes();
		append(bytes);
		append((byte)0);
		return this;
	}
	
	/*--------------------------------------------------------------*/
	
	/** something */
	public byte[] array;
	/** something else */
	public int length=0;
	/** Temporary buffer for numeric formatting operations */
	private final byte[] numbuffer=KillSwitch.allocByte1D(40);

	private static final byte tab='\t';
	private static final byte newline='\n';
	/** Byte array containing ASCII digits 0-9 */
	public static final byte[] numbers=new byte[] {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	/** Byte representation of "null" string */
	public static final byte[] nullBytes="null".getBytes();
	/** Byte representation of "false" string */
	public static final byte[] fbool="false".getBytes();
	/** Byte representation of "true" string */
	public static final byte[] tbool="true".getBytes();

	public static final byte[] ones100, tens100;

	public static final double[] decimalMult, decimalInvMult;
	/** Powers of 10 as long values for numeric formatting */
	public static final long[] longMult;
	/** Format strings for different decimal place counts */
	public static final String[] decimalFormat;
	/** Maximum allowed array length to prevent overflow */
	public static final int MAXLEN=Integer.MAX_VALUE-20;
	
	static{
		ones100=new byte[100];
		tens100=new byte[100];
		for(int i=0; i<100; i++){
			ones100[i]=(byte)('0'+i%10);
			tens100[i]=(byte)('0'+i/10);
		}
		longMult=new long[50];
		decimalMult=new double[50];
		decimalInvMult=new double[50];
		decimalFormat=new String[50];
		decimalMult[0]=decimalInvMult[0]=longMult[0]=1;
		decimalFormat[0]="%.0f";
		for(int i=1; i<50; i++){
			longMult[i]=longMult[i-1]*10;
			decimalMult[i]=decimalMult[i-1]*10;
			decimalInvMult[i]=1/decimalMult[i];
			decimalFormat[i]="%."+i+"f";
		}
	}
	
	public ByteBuilder appendUint32(long val) {
	    append((byte)(val & 0xFF));
	    append((byte)((val >> 8) & 0xFF));
	    append((byte)((val >> 16) & 0xFF));
	    append((byte)((val >> 24) & 0xFF));
	    return this;
	}
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String in=args[0];
		String out=(args.length<2 || "null".equalsIgnoreCase(args[1])) ? null : args[1];
		int type=(args.length>2 ? Integer.parseInt(args[2]) : 4);
		FileFormat ffin=FileFormat.testInput(in, null, false);
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.FASTQ, null, false, true, false, false);
		ByteFile bf=ByteFile.makeByteFile(ffin, type);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		ByteBuilder bb=new ByteBuilder(32768);
//		assert(false) : out+", "+ffout+", "+bf+", "+bsw;
		long lines=0, bytes=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			bb.append(line);
			lines++;
			bytes+=line.length;
//			assert(false) : "BB "+bb.length;
			if(bb.length>=16384) {
				if(bsw!=null) {
					bsw.print(bb);
//					assert(false) : "Wrote "+bb.length;
				}
				bb.clear();
			}
		}
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
		bf.close();
		if(bsw!=null) {
			boolean error=bsw.poisonAndWait();
			assert(!error);
		}
		t.stop();
		System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
	}
	
}
