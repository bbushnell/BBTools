package jgi;

import java.util.Arrays;

import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Sep 19, 2012
 *
 */
public final class Info {
	
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
	
	public static void printInfo(String s){
		long r=info(s);
		double bits=Math.log(r)/Math.log(2);
		System.out.println(s+"\nlen="+s.length()+" \tinfo = "+Tools.format("%.2f", bits)+" bits. \t("+r+")");
	}
	
	public static long info(String s){
		return info(s.getBytes(), 0, s.length());
	}
	
	public static int infoInBits(final byte[] array, final int from, final int len){return 63-Long.numberOfLeadingZeros(info(array, from, len));}
	public static double infoInBitsDouble(final byte[] array, final int from, final int len){return Math.log(info(array, from, len))*invlog2;}
	public static long info(final byte[] array){return info(array, 0, array.length);}
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

	public static int prefixForInfoBits(final byte[] array, final int bits){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, 0);}
	public static int prefixForInfoBits(final byte[] array, final int bits, final int from){assert(bits>=0 && bits<63);return prefixForInfo(array, 1L<<bits, from);}
	public static int prefixForInfo(final byte[] array, final long info){return prefixForInfo(array, info, 0);}
	
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
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	private static final long MAX=Long.MAX_VALUE;
	private static final double invlog2=1.0/Math.log(2);
}
