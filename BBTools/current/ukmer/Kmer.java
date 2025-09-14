package ukmer;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Jul 9, 2015
 *
 */
public class Kmer implements Cloneable {
	
	public Kmer(Kmer o){
		this(o.k, o.mult);
		setFrom(o);
	}
	
	public Kmer(int kbig_){
		this(getK(kbig_), getMult(kbig_));
	}
	
	@Override
	public Kmer clone(){
		return new Kmer(this);
	}
	
	public Kmer(int k_, int mult_){
		k=k_;
		mult=mult_;
		maxindex=mult-1;
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		coreMask=toCoreMask(k);
		
		kbig=k*mult;
		array1=new long[mult];
		array2=new long[mult];
		key=null;
	}
	
	public static final long toCoreMask(int k){
//		System.err.println(k+", "+MASK_CORE);
		return MASK_CORE ? ((~((-1L)<<(2*k)))>>2)&~(3L) : -1L;
	}
	
	public static int getMult(int kbig){
		final int mult=getMult0(kbig);
//		assert(mult==getMult0(kbig*(mult/mult))) : mult+", "+getMult0(mult*(kbig/mult));
		return mult;
	}
	
	public static int getKbig(int kbig){
		int x=getMult(kbig)*getK(kbig);
		assert(x<=kbig) : x+", "+kbig;
		assert(kbig>31 || x==kbig);
		return x;
	}
	
	private static int getMult0(int kbig){
//		if(true){return 2;}//TODO: 123 //Enable to allow multi-word arrays for k<32
		final int word=31;
		
		final int mult1=(kbig+word-1)/word;
		final int mult2=Tools.max(1, kbig/word);
		if(mult1==mult2){return mult1;}

		final int k1=Tools.min(word, kbig/mult1);
		final int k2=Tools.min(word, kbig/mult2);
		
		final int kbig1=k1*mult1;
		final int kbig2=k2*mult2;
		
		assert(kbig1<=kbig);
		assert(kbig2<=kbig);
		assert(mult2<=mult1);

//		assert(false) : mult1+", "+mult2+", "+k1+", "+k2;
		
		final int mult=kbig2>=kbig1 ? mult2 : mult1;
		
		return mult;
	}
	
	public static int getK(int kbig){
		int mult=getMult(kbig);
		int x=kbig/mult;
		assert(x*mult<=kbig) : x+", "+kbig;
		assert(x<=31) : kbig+", "+mult+", "+x;
		return x;
	}
	
	public Kmer setFrom(Kmer o){
		for(int i=0; i<mult; i++){
			array1[i]=o.array1[i];
			array2[i]=o.array2[i];
			len=o.len;
		}
		incarnation++;
		return this;
	}
	
	public Kmer setFrom(long[] array){
		for(int i=0; i<mult; i++){
			array1[i]=array[i];
		}
		fillArray2();
		incarnation++;
		return this;
	}
	
	public void clear() {
		len=0;
		for(int i=0; i<mult; i++){
			array1[i]=0;
			array2[i]=0;
		}
		lastIncarnation=-1;
		incarnation=0;
		//incarnation++;
	}
	
	public void clearFast() {
		len=0;
		lastIncarnation=-1;
		incarnation=0;
		//incarnation++;
	}
	
	public boolean verify(boolean update){
//		boolean b=verify();
//		if(b){
//			if(update){update();}
//			b=verify();
//			assert(len<kbig || incarnation==lastIncarnation);
//		}
		if(update){
			update();
			assert(len<kbig || incarnation==lastIncarnation) : "incarnation="+incarnation+", last="+lastIncarnation+", len="+len+", kbig="+kbig;
		}
		boolean b=verify();
		return b;
	}
	
	private boolean verify(){
		if(len<kbig){return true;}
		for(int i=maxindex, j=0; j<mult; j++, i--){
			long kmer=array1[i];
			long rkmer=array2[j];
			if(kmer!=AminoAcid.reverseComplementBinaryFast(rkmer, k)){
//				assert(false) : Arrays.toString(array1);
				return false;
			}
		}
		assert(incarnation==lastIncarnation) : "incarnation="+incarnation+", last="+lastIncarnation+", len="+len+", kbig="+kbig;
		return true;
	}
	
	public byte addRight(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
	public byte addRight(final char b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addRightNumeric(x)];
	}
	
	public byte addLeft(final byte b){
		long x=AminoAcid.baseToNumber[b];
		return AminoAcid.numberToBase[(int)addLeftNumeric(x)];
	}
	
	public long addRightNumeric(long x){
		long x2;
		
		if(x<0){
			x=0;
			x2=3;
			len=0;
		}else{
			x2=AminoAcid.numberToComplement[(int)x];
			len++;
		}
		
		for(int i=maxindex, j=0; j<mult; j++, i--){
			
			long y=(array1[i]>>>shift2)&3L;
			long y2=array2[j]&3L;
			
			//Update kmers
			array1[i]=((array1[i]<<2)|x)&mask;
			array2[j]=((array2[j]>>>2)|(x2<<shift2))&mask;
			
			x=y;
			x2=y2;
		}
		incarnation++;
		return x;
	}
	
	public long addLeftNumeric(long x){
		assert(x>=0 && x<4) : x;
		long x2=AminoAcid.numberToComplement[(int)x];
		
		assert(x>=0);
		assert(len>=kbig);
		
		for(int i=0, j=maxindex; i<mult; i++, j--){

			long y=array1[i]&3L;
			long y2=(array2[j]>>>shift2)&3L;
			
			//Update kmers
			array1[i]=(array1[i]>>>2)|(x<<shift2);
			array2[j]=((array2[j]<<2)|x2)&mask;
			
			x=y;
			x2=y2;
		}
		incarnation++;
		return x;
	}
	
	public void fillArray2() {
		for(int i=maxindex, j=0; j<mult; j++, i--){
			array2[j]=AminoAcid.reverseComplementBinaryFast(array1[i], k);
		}
		len=kbig;
		incarnation++;
	}
	
	@Override
	public String toString(){
//		update();
		assert(verify(true));
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<mult; i++){
			bb.appendKmer(array1[i], k);
//			bb.append(" ");
		}
////		bb.append("~");
//		for(int i=0; i<mult; i++){
//			bb.appendKmer(array2[i], k);
////			bb.append(" ");
//		}
		return bb.toString();
	}
	
	public boolean equals(Kmer x){
		if(xor()!=x.xor()){return false;}
		return AbstractKmerTableU.equals(key(), x.key());
	}
	
	public boolean sameOrientation(Kmer x){
		if(xor()!=x.xor()){return false;}
		return Tools.equals(array1, array2);
	}
	
	public int compareTo(Kmer x){
		return compare(key(), x.key());
	}
	
	public int compareTo(long[] key2){
		assert(false);
		return compare(key(), key2);
	}
	
	public static int compare(long[] key1, long[] key2){
//		assert(false); //Why was this here?
		return AbstractKmerTableU.compare(key1, key2);
	}
	
	public static boolean equals(long[] key1, long[] key2){
		assert(false);
		return AbstractKmerTableU.equals(key1, key2);
	}
	
	public long[] array1(){return array1;}
		
	public long[] array2(){return array2;}
	
	/** WARNING!
	 * Do not confuse this with xor()! */
	public long[] key(){
		update();
//		assert(verify(false));
		return key;
	}
	
	public boolean corePalindrome(){//TODO: This can be set as a flag from setKey0
		update();
		return corePalindrome;
	}
	
	private void setKey0(){
		corePalindrome=false;
		key=array1;
		if(!rcomp){return;}
		for(int i=0; i<mult; i++){
			final long a=array1[i]&coreMask, b=array2[i]&coreMask;
			if(a>b){return;}
			else if(a<b){
				key=array2;
				return;
			}
		}
		corePalindrome=true;
		setKey0safe();
	}
	
	private void setKey0safe(){
		key=array1;
		for(int i=0; i<mult; i++){
			final long a=array1[i], b=array2[i];
			if(a>b){break;}
			else if(a<b){
				key=array2;
				break;
			}
		}
	}
	
	public static long xor(long[] key, long coreMask){
		long xor=key[0]&coreMask;
		for(int i=1; i<key.length; i++){
			xor=(Long.rotateLeft(xor, 25))^(key[i]&coreMask);
		}
		return xor&mask63;
	}
	
	/** WARNING!
	 * Do not confuse this with key()! */
	public long xor(){
		update();
		return lastXor;
	}

	/**
	 * @param divisor
	 * @return This kmer's xor modulo the divisor
	 */
	public int mod(int divisor) {
		int x=(int)(xor()%divisor);
//		System.err.println(xor()+"%"+value+"="+x);
		return x;
	}
	
	public void rcomp() {
		long[] temp=array1;
		array1=array2;
		array2=temp;
	}
	
	private void update(){
		if(verbose){System.err.println("update() - len="+len);}
		assert(TESTMODE || len>=kbig) : len+", "+kbig+", "+array1[0];
		if(incarnation==lastIncarnation){return;}
		setKey0();
		lastXor=xor0();
		lastIncarnation=incarnation;
		if(verbose){System.err.println("After update - kmer "+this+"; key="+Arrays.toString(key)+"; a1="+Arrays.toString(array1())+"; a2="+Arrays.toString(array2()));}
	}
	
	private long xor0(){
		return xor(key, coreMask);
	}
	
	public String arraysToString() {
		return "key="+Arrays.toString(key)+", a1="+Arrays.toString(array1)+", a2="+Arrays.toString(array2);
	}
	
	public final int gc(){
		int gc=0;
		for(long kmer : array1){
			while(kmer>0){
				long x=kmer&3;
				kmer>>>=2;
				if(x==1 || x==2){gc++;}
			}
		}
		return gc;
	}
	
	static boolean rcomp=true;
	
	private long lastXor=-1;
	private long incarnation=0;
	private long lastIncarnation=-1;
	private boolean corePalindrome=false;
	private long[] key=null;
	
	private long[] array1;
	private long[] array2;
	public final int kbig;
	public final int k;
	final int mult, maxindex;
	
	private final int shift;
	private final int shift2;
	private final long mask;
	private final long coreMask;
	
	public int len=0; //TODO: Make private; use getter.
	public final int len(){return len;}
	
	public static boolean MASK_CORE=false;
	private static final long mask63=Long.MAX_VALUE;
	private static final boolean TESTMODE=false; //123
	private static final boolean verbose=false;
}
