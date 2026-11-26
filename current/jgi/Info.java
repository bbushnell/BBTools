package jgi;

import java.util.Arrays;

import shared.Tools;

/**
 * Calculates information content of DNA sequences using combinatorics.
 * Computes sequence complexity based on base composition and frequency distribution.
 * Used for analyzing sequence entropy and finding prefixes with specific information content.
 *
 * @author Brian Bushnell
 * @date Sep 19, 2012
 */
public final class Info {
	
	/** Prints built-in demonstrations or command-line sequence summaries. */
	public static void main(String[] args){
		
		if(args.length>0){
			if(args.length==2 && Tools.isDigit(args[1].charAt(0))){
				byte[] s=args[0].getBytes();
				int b=Integer.parseInt(args[1]);
				int len=prefixForInfoBits(s, b);
				if(len<0){
					System.out.println("Input string only contains "+Tools.format("%.2f",infoInBitsDouble(s, 0, s.length))+" bits.");
				}else{
					System.out.println("Prefix needed for "+b+" bits is length "+len+": "+args[0].substring(0, len));
//					assert(false) : "TODO: This is clearly broken.";
				}
			}else{
				for(String s : args){
					printInfo(s);
					System.out.println();
				}
			}
			System.exit(0);
		}
		
		System.out.println();
		printInfo("");
		System.out.println();
		printInfo("A");
		System.out.println();
		printInfo("AG");
		System.out.println();
		printInfo("AGT");
		System.out.println();
		printInfo("AANAA");
		System.out.println();
		printInfo("GGGGGGGCGGG");
		System.out.println();
		printInfo("CGGGGGGGGGG");
		System.out.println();
		printInfo("AGTCAGTCCTAGNGTACGT");
		System.out.println();
		printInfo("AGTCAGTCAGTCAGTC");
		System.out.println();
		printInfo("GCGCGCGCGCGCGCGC");
		System.out.println();
		
		String[] s=new String[] {"A", "G", "C", "T", ""};
		for(int i=0; i<40; i++){
			System.out.println();
			s[4]=s[4]+s[i%4];
			printInfo(s[4]);
		}
		
		System.out.println("PrefixForBits for AAAATATATGAAATGCATGCAATATGTTATGAAA");
		for(int i=0; i<60; i+=2){
			System.out.println(i+"\t"+prefixForInfoBits("AAAATATATGAAATGCATGCAATATGTTATGAAA".getBytes(), i));
		}

		
		System.out.println("PrefixForBits for GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");
		for(int i=0; i<60; i+=2){
			System.out.println(i+"\t"+prefixForInfoBits("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC".getBytes(), i));
		}

		
		System.out.println("PrefixForBits for ACGTACGTACGTACGTACGTACGTACGTACGTAC");
		for(int i=0; i<63; i+=2){
			System.out.println(i+"\t"+prefixForInfoBits("ACGTACGTACGTACGTACGTACGTACGTACGTAC".getBytes(), i));
		}
	}
	
	/** Prints one sequence, its length in characters, its approximate information in bits,
	 * and its capped multinomial count.
	 * @param s ASCII DNA sequence to print; the count is delegated to info(String)
	 */
	public static void printInfo(String s){
		long r=info(s);
		double bits=Math.log(r)/Math.log(2);
		System.out.println(s+"\nlen="+s.length()+" \tinfo = "+Tools.format("%.2f", bits)+" bits. \t("+r+")");
	}
	
	/** Counts arrangements of recognized bases in an ASCII string.
	 * @param s ASCII DNA sequence; converted with the default charset, then the
	 * first s.length() bytes are scored
	 * @return Multinomial count of recognized bases, capped at Long.MAX_VALUE
	 */
	public static long info(String s){
		return info(s.getBytes(), 0, s.length());
	}
	
	/** Returns the integer floor of a slice's capped count in log2 units.
	 * @param array ASCII base array
	 * @param from First absolute index to score
	 * @param len Maximum number of bytes to inspect
	 * @return floor(log2(info(array, from, len)))
	 */
	public static int infoInBits(final byte[] array, final int from, final int len){return 63-Long.numberOfLeadingZeros(info(array, from, len));}
	/** Returns the approximate log2 of a slice's capped count.
	 * @param array ASCII base array
	 * @param from First absolute index to score
	 * @param len Maximum number of bytes to inspect
	 * @return log2(info(array, from, len)) as a double
	 */
	public static double infoInBitsDouble(final byte[] array, final int from, final int len){return Math.log(info(array, from, len))*invlog2;}
	/** Counts arrangements of recognized bases in a complete byte array.
	 * @param array ASCII base array; A/C/G/T/U are recognized and other bytes ignored
	 * @return Multinomial count of recognized bases, capped at Long.MAX_VALUE
	 */
	public static long info(final byte[] array){return info(array, 0, array.length);}
	/** Counts arrangements of recognized bases in a bounded byte slice.
	 * Supported calls use nonnegative from and len values whose sum does not overflow.
	 * The high endpoint is clipped to array.length.
	 * @param array ASCII base array; not modified
	 * @param from First absolute index to score
	 * @param len Maximum number of bytes to inspect
	 * @return Multinomial count of recognized bases in the slice, capped at Long.MAX_VALUE
	 */
	public static long info(final byte[] array, final int from, final int len){
		int[] counts=new int[4];
		long r=1;
		int used=0;
		for(int i=from, lim=min(from+len, array.length); i<lim; i++){
//			System.out.print(((char)array[i])+" -> ");
			byte num=baseToNumber[array[i]];
//			System.out.println(num);
			if(num>=0){
				counts[num]++;
				used++;
				
				// JGI-216: Preserve a representable quotient when the product would overflow.
				if(used>32 && used>MAX/r){//overflow
					r=nextInfo(r, used, counts[num]);
					if(r==MAX){return MAX;}
					continue;
				}
				r=r*used;
				
				/* alternate method */
//				long temp=r*used;
//			    if(used>32 && temp/used!=r){//overflow
//			    	return MAX;
//			    }
//			    r=temp;
			    
			    r=r/counts[num];
			}
		}
		return r;
	}

	/** Finds the first prefix whose count reaches a power-of-two threshold.
	 * @param array ASCII base array to scan from index 0
	 * @param bits Power-of-two threshold in bits; assertions require bits in [0, 63)
	 * @return Absolute exclusive endpoint reaching 2^bits, or -1 if unreachable;
	 * a threshold of 1 can return the starting index before reading any byte
	 */
	public static int prefixForInfoBits(final byte[] array, final int bits){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, 0);}
	/** Finds the first prefix whose count reaches a power-of-two threshold.
	 * @param array ASCII base array
	 * @param bits Power-of-two threshold in bits; assertions require bits in [0, 63)
	 * @param from First absolute index to scan
	 * @return Absolute exclusive endpoint reaching 2^bits, or -1 if unreachable;
	 * a threshold of 1 can return from before reading any byte
	 */
	public static int prefixForInfoBits(final byte[] array, final int bits, final int from){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, from);}
	/** Finds the first prefix whose count reaches a raw count threshold.
	 * @param array ASCII base array to scan from index 0
	 * @param info Count threshold; assertions require a nonnegative value
	 * @return Absolute exclusive endpoint reaching info, or -1 if unreachable;
	 * thresholds 0 and 1 return the starting index before reading any byte
	 */
	public static int prefixForInfo(final byte[] array, final long info){return prefixForInfo(array, info, 0);}
	
	/** Finds the first prefix whose count reaches a raw count threshold.
	 * Counts restart at from. Ignored ASCII bytes advance the returned endpoint but do
	 * not contribute to base counts.
	 * @param array ASCII base array; not modified
	 * @param info Count threshold; assertions require a nonnegative value
	 * @param from First absolute index to scan
	 * @return Absolute exclusive endpoint reaching info, or -1 if unreachable;
	 * thresholds 0 and 1 return from before reading any byte
	 */
	public static int prefixForInfo(final byte[] array, final long info, final int from){
		assert(info>=0);
		int[] counts=new int[4];
		long r=1;
		int used=0;
		int i=from;
		for(; i<array.length && r<info; i++){
//			System.out.print(((char)array[i])+" -> ");
			byte num=baseToNumber[array[i]];
//			System.out.println(num);
			if(num>=0){
				counts[num]++;
				used++;
				
				if(used>32 && used>MAX/r){//overflow
					r=nextInfo(r, used, counts[num]);
					continue;// JGI-217: Let the loop update produce the exclusive endpoint.
				}
				r=r*used;
				
				/* alternate method */
//				long temp=r*used;
//			    if(used>32 && temp/used!=r){//overflow
//			    	return MAX;
//			    }
//			    r=temp;
			    
			    r=r/counts[num];
//
//			    {
//			    	String s=new String(array).substring(0, i+1);
//			    	System.out.println("\n"+s);
//			    	System.out.println("For len "+i+": r="+r+", bits="+(63-Long.numberOfLeadingZeros(r))+"\t->\t"+(Math.log(r)*invlog2));
//			    	System.out.println(infoInBitsDouble(s.getBytes(), 0, i+1));
//			    	System.out.println(info(s.getBytes(), 0, i+1));
//			    }
			}
		}
		return r<info ? -1 : i;
	}
	
	/** Computes a capped multinomial count without overflowing the intermediate product.
	 * Cancelling used/count makes the remaining denominator divide the exact prior count.
	 * @param prior Exact previous count, before saturation
	 * @param used Number of recognized bases including the new base
	 * @param count Number of occurrences of the new base including this occurrence
	 * @return Next count, capped at Long.MAX_VALUE
	 */
	private static long nextInfo(final long prior, final int used, final int count){
		assert(prior>0 && count>0 && count<=used) :
				"The multinomial recurrence needs a positive exact count and one updated bin: "+prior+", "+used+", "+count;
		int a=used, b=count;
		while(b!=0){
			final int remainder=a%b;
			a=b;
			b=remainder;
		}
		final int numerator=used/a, denominator=count/a;
		assert(prior%denominator==0) :
				"Coprime cancellation makes the denominator divide the prior multinomial count: "+prior+", "+used+", "+count;
		final long reduced=prior/denominator;
		return reduced>MAX/numerator ? MAX : reduced*numerator;
	}

	private static final byte[] numberToBase={
		'A','C','G','T','N'
	};
	
	/** Initial mutable mapping from 7-bit ASCII bytes to base numbers.
	 * A/C/G/T map to 0/1/2/3 case-insensitively, U/u map to T, and all other
	 * 7-bit ASCII bytes including N are ignored with -1.
	 */
	public static final byte[] baseToNumber=new byte[128];
	
	static{
		Arrays.fill(baseToNumber, (byte)-1);
		for(int i=0; i<numberToBase.length; i++){
			char x=(char)numberToBase[i];
			if(x=='A' || x=='C' || x=='G' || x=='T'){
				baseToNumber[x]=(byte)i;
				baseToNumber[Tools.toLowerCase(x)]=(byte)i;
			}
		}
		baseToNumber['U']=3;
		baseToNumber['u']=3;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	private static final long MAX=Long.MAX_VALUE;
	private static final double invlog2=1.0/Math.log(2);
}
