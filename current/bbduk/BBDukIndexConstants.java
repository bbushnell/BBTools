package bbduk;

import dna.AminoAcid;
import shared.Tools;

/**
 * Stores all k-mer hashing and manipulation constants derived from
 * command-line arguments. This class is static.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public final class BBDukIndexConstants{

	private BBDukIndexConstants(){}
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	public static boolean amino;
	public static int bitsPerBase;
	public static int maxSymbol;
	public static int symbols;
	public static int symbolArrayLen;
	public static int symbolSpace;
	public static long symbolMask;
	public static int k;
	public static int k2;
	public static int shift;
	public static int shift2;
	public static long mask;
	public static long kmask;
	public static int mink;
	public static long middleMask;

	public static long[] clearMasks;
	public static long[][] setMasks;
	public static long[] leftMasks;
	public static long[] rightMasks;
	public static long[] lengthMasks;
	
	public static byte[] symbolToNumber;
	public static byte[] symbolToNumber0;
	public static byte[] symbolToComplementNumber0;
	
	public static long refReads=0;
	public static long refBases=0;
	public static long refKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Utility        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Computes reverse complement of kmer. */
	public static long rcomp(long kmer, int len){return amino?kmer:AminoAcid.reverseComplementBinaryFast(kmer,len);}
	
	/** Transforms a kmer into a canonical value stored in the table. */
	public static long toValue(long kmer,long rkmer,long lengthMask,boolean rcompFlag,long middleMask){
		final long value=(rcompFlag?Tools.max(kmer,rkmer):kmer);
		return (value&middleMask)|lengthMask;
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	public static final boolean isFullyDefined(byte symbol){return symbol>=0&&symbolToNumber[symbol]>=0;}
	
	/** Initializes all k-mer hashing constants based on k and alphabet type. */
	public static void initialize(int k_,int mink_,boolean amino_,int midMaskLen_){
		k=k_;
		mink=mink_;
		amino=amino_;
		k2=k-1;
		
		bitsPerBase=(amino?5:2);
		maxSymbol=(amino?20:3);
		symbols=maxSymbol+1;
		symbolArrayLen=(64+bitsPerBase-1)/bitsPerBase;
		symbolSpace=(1<<bitsPerBase);
		symbolMask=symbolSpace-1;
		
		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63?-1L:~((-1L)<<shift));
		
		symbolToNumber=AminoAcid.symbolToNumber(amino);
		symbolToNumber0=AminoAcid.symbolToNumber0(amino);
		symbolToComplementNumber0=AminoAcid.symbolToComplementNumber0(amino);
		
		clearMasks=new long[symbolArrayLen];
		leftMasks=new long[symbolArrayLen];
		rightMasks=new long[symbolArrayLen];
		lengthMasks=new long[symbolArrayLen];
		setMasks=new long[symbols][symbolArrayLen];
		
		for(int i=0;i<symbolArrayLen;i++){
			clearMasks[i]=~(symbolMask<<(bitsPerBase*i));
			leftMasks[i]=((-1L)<<(bitsPerBase*i));
			rightMasks[i]=~((-1L)<<(bitsPerBase*i));
			lengthMasks[i]=((1L)<<(bitsPerBase*i));
			for(long j=0;j<symbols;j++){
				setMasks[(int)j][i]=(j<<(bitsPerBase*i));
			}
		}
		
		kmask=lengthMasks[k]; 
		
		if(midMaskLen_>0){
			assert(k>midMaskLen_+1);
			int bits=midMaskLen_*bitsPerBase;
			int shift=((k-midMaskLen_)/2)*bitsPerBase; 
			middleMask=~((~((-1L)<<bits))<<shift);
		}else{
			middleMask=-1L;
		}
	}
}