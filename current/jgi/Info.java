package jgi;

import java.util.Arrays;

import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Sep 19, 2012
 *
 */
public final class Info {
	
	/**
	 * Program entry point with dual modes for information content calculation.
	 * Mode 1: Two arguments (sequence, bits) finds prefix length for target bits.
	 * Mode 2: Variable arguments prints information content for each sequence.
	 * @param args Command-line arguments - either sequence strings or sequence+bits pair
	 */
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
	
	/** Prints sequence information including length, information content in bits, and raw value.
	 * @param s The DNA sequence to analyze */
	public static void printInfo(String s){
		long r=info(s);
		double bits=Math.log(r)/Math.log(2);
		System.out.println(s+"\nlen="+s.length()+" \tinfo = "+Tools.format("%.2f", bits)+" bits. \t("+r+")");
	}
	
	/**
	 * Calculates information content for a string sequence.
	 * @param s The DNA sequence string
	 * @return Information content as combinatorial value
	 */
	public static long info(String s){
		return info(s.getBytes(), 0, s.length());
	}
	
	/**
	 * Calculates information content in bits using bit position of highest set bit.
	 *
	 * @param array The sequence byte array
	 * @param from Starting position in array
	 * @param len Length of sequence to analyze
	 * @return Information content in bits (integer approximation)
	 */
	public static int infoInBits(final byte[] array, final int from, final int len){return 63-Long.numberOfLeadingZeros(info(array, from, len));}
	/**
	 * Calculates information content in bits using logarithm for precise floating-point result.
	 *
	 * @param array The sequence byte array
	 * @param from Starting position in array
	 * @param len Length of sequence to analyze
	 * @return Information content in bits (floating-point)
	 */
	public static double infoInBitsDouble(final byte[] array, final int from, final int len){return Math.log(info(array, from, len))*invlog2;}
	/**
	 * Calculates information content for entire byte array.
	 * @param array The sequence byte array
	 * @return Information content as combinatorial value
	 */
	public static long info(final byte[] array){return info(array, 0, array.length);}
	/**
	 * Calculates information content using combinatorial formula based on base frequencies.
	 * Uses multinomial coefficient: n! / (c1! * c2! * c3! * c4!) where n is sequence length
	 * and c1-c4 are counts of each base type. Includes overflow protection.
	 *
	 * @param array The sequence byte array
	 * @param from Starting position in array
	 * @param len Length of sequence to analyze
	 * @return Information content as combinatorial value, MAX if overflow occurs
	 */
	public static long info(final byte[] array, final int from, final int len){
		short[] counts=new short[4];
		long r=1;
		int used=0;
		for(int i=from, lim=min(from+len, array.length); i<lim; i++){
//			System.out.print(((char)array[i])+" -> ");
			byte num=baseToNumber[array[i]];
//			System.out.println(num);
			if(num>=0){
				counts[num]++;
				used++;
				
				if(used>32 && used>MAX/r){//overflow
//					System.out.println("***");
					return MAX;
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

	/**
	 * Finds shortest prefix length that contains at least the specified bits of information.
	 * @param array The sequence byte array
	 * @param bits Target information content in bits
	 * @return Prefix length needed, or -1 if sequence has insufficient information
	 */
	public static int prefixForInfoBits(final byte[] array, final int bits){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, 0);}
	/**
	 * Finds shortest prefix length starting from specified position that contains
	 * at least the specified bits of information.
	 *
	 * @param array The sequence byte array
	 * @param bits Target information content in bits
	 * @param from Starting position in array
	 * @return Prefix length needed, or -1 if sequence has insufficient information
	 */
	public static int prefixForInfoBits(final byte[] array, final int bits, final int from){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, from);}
	/**
	 * Finds shortest prefix length that contains at least the specified raw information value.
	 * @param array The sequence byte array
	 * @param info Target information content as combinatorial value
	 * @return Prefix length needed, or -1 if sequence has insufficient information
	 */
	public static int prefixForInfo(final byte[] array, final long info){return prefixForInfo(array, info, 0);}
	
	/**
	 * Finds shortest prefix length starting from specified position that contains
	 * at least the specified raw information value. Incrementally calculates information
	 * content using the same combinatorial approach as info() method.
	 *
	 * @param array The sequence byte array
	 * @param info Target information content as combinatorial value
	 * @param from Starting position in array
	 * @return Prefix length needed, or -1 if sequence has insufficient information
	 */
	public static int prefixForInfo(final byte[] array, final long info, final int from){
		assert(info>=0);
		short[] counts=new short[4];
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
//					System.out.println("***");
					return i;
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
	
	/** Array mapping base numbers (0-4) to ASCII characters (A, C, G, T, N) */
	private static final byte[] numberToBase={
		'A','C','G','T','N'
	};
	
	/** Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise */
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
	
	/**
	 * Returns the smaller of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return The smaller value
	 */
	private static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns the larger of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return The larger value
	 */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/** Maximum value used for overflow detection in information calculations */
	private static final long MAX=Long.MAX_VALUE;
	/** Inverse of natural logarithm of 2, used for converting to bits: 1/ln(2) */
	private static final double invlog2=1.0/Math.log(2);
}
