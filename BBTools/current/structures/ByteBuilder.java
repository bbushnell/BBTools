package structures;

import java.io.Serializable;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.LineParser;
import shared.Tools;
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
	
	public ByteBuilder(){
		array=KillSwitch.allocByte1D(32);
	}
	
	public ByteBuilder(int initial){
		assert(initial>=1) : initial;
		array=KillSwitch.allocByte1D(initial);
	}
	
	public ByteBuilder(Object o){
		String s=o.toString();
		array=KillSwitch.allocByte1D(s.length()+1);
		append(s);
	}
	
	public ByteBuilder(byte[] array_){
		assert(array_!=null);
		array=array_;
		length=array.length;
	}
	
	public ByteBuilder(ByteBuilder bb){
		array=bb.toBytes();
		length=bb.length();
	}

	@Override
	public CharSequence subSequence(int start, int end) throws IndexOutOfBoundsException{
		if(start<0 || end>=length()){throw new IndexOutOfBoundsException();}
		return new ByteBuilder(KillSwitch.copyOfRange(array, start, end));
	}
	
	public int indexOf(char c) {return indexOf((byte)c);}
	public int lastIndexOf(char c) {return lastIndexOf((byte)c);}
	public int count(char c) {return count((byte)c);}
	
	public int indexOf(byte c) {
		for(int i=0; i<length; i++) {
			if(array[i]==c) {return i;}
		}
		return -1;
	}
	
	public int lastIndexOf(byte c) {
		for(int i=length-1; i>=0; i--) {
			if(array[i]==c) {return i;}
		}
		return -1;
	}

	public int count(byte c) {
		int count=0;
		for(int i=0; i<length; i++) {
			count+=(array[i]==c) ? 1 : 0;
		}
		return count;
	}

//	public ByteBuilder append(float x, int places){return append(Tools.format(decimalFormat[places], x));}
//	public ByteBuilder append(double x, int places){return append(Tools.format(decimalFormat[places], x));}

	public ByteBuilder appendSlow(float x){return append(Float.toString(x));}
	public ByteBuilder appendSlow(double x){return append(Double.toString(x));}
	public ByteBuilder appendSlow(float x, int decimals){return append(Tools.format("%."+decimals+"f", x));}
	public ByteBuilder appendSlow(double x, int decimals){return append(Tools.format("%."+decimals+"f", x));}
	public ByteBuilder append(boolean x){return append(x ? tbool : fbool);}

	public ByteBuilder append(char a, char b) {
		append(a);append(b);
		return this;
	}
	public ByteBuilder append(byte a, byte b) {
		append(a);append(b);
		return this;
	}
	public ByteBuilder append(char a, char b, char c) {
		append(a);append(b);append(c);
		return this;
	}
	public ByteBuilder append(byte a, byte b, byte c) {
		append(a);append(b);append(c);
		return this;
	}
	public ByteBuilder append(char a, char b, char c, char d) {
		append(a);append(b);append(c);append(d);
		return this;
	}
	public ByteBuilder append(byte a, byte b, byte c, byte d) {
		append(a);append(b);append(c);append(d);
		return this;
	}
	public ByteBuilder append(char x, int y){
		append(x); append(y); return this;
	}
	public ByteBuilder append(char x, long y){
		append(x); append(y); return this;
	}
	
	public ByteBuilder append(char x){
		if(length>=array.length){expand();}
		array[length]=(byte)x;
		length++;
		return this;
	}
	public ByteBuilder append(byte x){
		if(length>=array.length){expand();}
		array[length]=x;
		length++;
		return this;
	}
	
	public ByteBuilder appendKmer(Kmer kmer) {
		return appendKmer(kmer.array1(), kmer.k);
	}
	
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
	
	public ByteBuilder appendTerm(LineParser lp, int term) {
		return lp.appendTerm(this, term);
	}
	
	public ByteBuilder appendln(int x){
		expand(12);
		append(x);
		return append('\n');
	}
	
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
	
	public ByteBuilder append(String x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
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
	
	public ByteBuilder append(StringBuilder x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(CharSequence x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	public ByteBuilder appendln(CharSequence x){
		expand(x.length()+1);
		append(x);
		array[length]='\n';
		length++;
		return this;
	}
	
//	public ByteBuilder append(Object x){
//		if(x==null){return append(nullBytes);}
//		return append(x.toString());
//	}
	
	public ByteBuilder append(byte[] x){
		if(x==null){x=nullBytes;}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=x[i];
			length++;
		}
		return this;
	}

	public ByteBuilder appendt(String x){return append(x).tab();}
	public ByteBuilder appendt(int x){return append(x).tab();}
	public ByteBuilder appendt(long x){return append(x).tab();}
	public ByteBuilder appendt(byte[] x){return append(x).tab();}
	public ByteBuilder appendt(float x, int decimals){return append(x, decimals).tab();}
	public ByteBuilder appendt(double x, int decimals){return append(x, decimals).tab();}
	
	public ByteBuilder appendln(byte[] x){
		expand(x.length+1);
		append(x);
		array[length]='\n';
		length++;
		return this;
	}
	
	public ByteBuilder appendQuality(byte[] x){
		if(x==null){return this;}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=(byte)(x[i]+33);
			length++;
		}
		return this;
	}
	
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
	
	public ByteBuilder append(long[] array, char delimiter){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i]);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	public ByteBuilder appendA48(long[] array, char delimiter, byte[] temp){
		if(array==null || array.length<1){return this;}
		if(temp==null) {temp=new byte[12];}
		for(int i=0; i<array.length; i++){
			appendA48(array[i], temp);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	public ByteBuilder appendA48(long value, byte[] temp){
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
	
	public ByteBuilder append(int[] array, char delimiter){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i]);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	public ByteBuilder append(float[] array, char delimiter, int decimals){
		if(array==null || array.length<1){return this;}
		for(int i=0; i<array.length; i++){
			append(array[i], decimals);
			append(delimiter);
		}
		length--;
		return this;
	}
	
	public ByteBuilder append(ByteBuilder bb){
		return append(bb.array, 0, bb.length);
	}
	
	public ByteBuilder appendln(ByteBuilder x){
		expand(x.length+1);
		append(x);
		array[length]='\n';
		length++;
		return this;
	}
	
	public ByteBuilder append(byte[] x, int len){
		return append(x, 0, len);
	}
	
	public ByteBuilder append(byte[] x, int start, int len){
//		if(x==null){x=nullBytes;}
		final int lim=(int)Tools.min(start+(long)len, x.length);
		len=Tools.min(len, x.length-start);
		expand(len);
		for(int i=start; i<lim; i++){
			array[length]=x[i];
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(char[] x){
		if(x==null){return append(nullBytes);}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=(byte)x[i];
			length++;
		}
		return this;
	}
	
	public ByteBuilder appendln(char[] x){
		expand(x.length+1);
		append(x);
		array[length]='\n';
		length++;
		return this;
	}
	
	/**
	 * Append a newline.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder nl(){return append('\n');}
	
	/**
	 * Append a tab.
	 * @return This ByteBuilder.
	 */
	public ByteBuilder tab(){return append('\t');}
	
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
	
	public byte get(int i){
		assert(i<length);
		return array[i];
	}
	
	public void set(int i, byte b){
		assert(i<length);
		array[i]=b;
	}
	
	public void set(int i, char b){
		assert(i<length) : i+", "+b+", "+length;
		array[i]=(byte)b;
	}
	
	@Override
	public char charAt(int i){
		assert(i<length);
		return (char)array[i];
	}

	public boolean endsWith(char c) {
		return length<1 ? false : array[length-1]==c;
	}

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
		return new String(array, 0, length);
	}
	
	public final byte[] toBytes(){
		return KillSwitch.copyOf(array, length);
	}
	
	public final byte[] toBytes(int from, int to){
		return KillSwitch.copyOfRange(array, from, to);
	}
	
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
	
	private final boolean isRoom(int x){
		return array.length-length>=x;
	}
	
	private final void expand(){
		long x=Tools.min(MAXLEN, array.length*2L);
		if(x<=array.length){
			throw new RuntimeException("Array overflow: "+x+"<="+array.length);
		}
		assert(((int)x)>array.length) : "Array overflow: "+x+"<="+array.length;
		array=KillSwitch.copyOf(array, (int)x);
	}

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
	
	private final void expand(int extra){
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
	
	public ByteBuilder reverse() {
		Tools.reverseInPlace(array, 0, length);
		return this;
	}
	
	public ByteBuilder reverseInPlace() {
		Tools.reverseInPlace(array, 0, length);
		return this;
	}
	
	public ByteBuilder reverseInPlace(int from, int toExclusive) {
		Tools.reverseInPlace(array, from, toExclusive);
		return this;
	}
	
	public ByteBuilder complementInPlace() {
		AminoAcid.complementBasesInPlace(array, length);
		return this;
	}
	
	public ByteBuilder reverseComplementInPlace() {
		AminoAcid.reverseComplementBasesInPlace(array, length);
		return this;
	}
	
	public final void ensureExtra(int extra){
		if(array.length-length<extra){expand(extra);}
	}
	
	public boolean isEmpty(){return length==0;}
	@Override
	public int length(){return length;}
	public ByteBuilder clear(){
		length=0;
		return this;
	}
	public ByteBuilder trimLast(int x){
		assert(x>=0);
		setLength(length-x);
		return this;
	}
	public ByteBuilder setLength(int x){
		assert(x>=0 && x<=array.length);
		length=x;
		return this;
	}
	
	/** something */
	public byte[] array;
	/** something else */
	public int length=0;
	private final byte[] numbuffer=KillSwitch.allocByte1D(40);

	public static final byte[] numbers=new byte[] {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	public static final byte[] nullBytes="null".getBytes();
	public static final byte[] fbool="false".getBytes();
	public static final byte[] tbool="true".getBytes();

	public static final byte[] ones100, tens100;

	public static final double[] decimalMult, decimalInvMult;
	public static final long[] longMult;
	public static final String[] decimalFormat;
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
	
}
