package bbduk;

import shared.Parser;
import shared.Shared;

/**
 * Abstract base class for BBDuk components (Processor, Index, etc.).
 * Injects configuration constants into local fields.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 19, 2025
 */
public abstract class BBDukObject{
	
	protected static final int BBDUK_SPEED_DIVISOR=17;
	protected static int speed=0;
	protected static int qSkip=1;
	
	protected BBDukObject(BBDukParser bbdp, Parser parser){
		initLocalFields(bbdp, parser);
		initConstants();
	}
	
	/** Copy constructor */
	protected BBDukObject(BBDukObject s){
		k=s.k;
		mink=s.mink;
		amino=s.amino;
		k2=s.k2; 
		rcomp=s.rcomp;
		forbidNs=s.forbidNs;
		minlen=s.minlen;
		midMaskLen=s.midMaskLen;
		minlen2=s.minlen2;
		
		bitsPerBase=s.bitsPerBase;
		maxSymbol=s.maxSymbol;
		symbols=s.symbols;
		shift=s.shift;
		shift2=s.shift2;
		mask=s.mask;
		symbolMask=s.symbolMask;
		symbolToNumber=s.symbolToNumber;
		symbolToNumber0=s.symbolToNumber0;
		symbolToComplementNumber0=s.symbolToComplementNumber0;
		clearMasks=s.clearMasks;
		setMasks=s.setMasks;
		leftMasks=s.leftMasks;
		rightMasks=s.rightMasks;
		lengthMasks=s.lengthMasks;
		middleMask=s.middleMask;
	}
	
	protected final void initLocalFields(BBDukParser bbdp, Parser parser){
		k=bbdp.k;
		mink=bbdp.mink;
		amino=Shared.AMINO_IN;
		k2=bbdp.k2; 
		rcomp=bbdp.rcomp;
		forbidNs=bbdp.forbidNs;
		minlen=k-1;
		midMaskLen=bbdp.midMaskLen; 
	}
	
	protected final void initConstants(){
		bitsPerBase=BBDukIndexConstants.bitsPerBase;
		maxSymbol=BBDukIndexConstants.maxSymbol;
		symbols=BBDukIndexConstants.symbols;
		shift=BBDukIndexConstants.shift;
		shift2=BBDukIndexConstants.shift2;
		mask=BBDukIndexConstants.mask;
		symbolMask=BBDukIndexConstants.symbolMask;
		symbolToNumber=BBDukIndexConstants.symbolToNumber;
		symbolToNumber0=BBDukIndexConstants.symbolToNumber0;
		symbolToComplementNumber0=BBDukIndexConstants.symbolToComplementNumber0;
		clearMasks=BBDukIndexConstants.clearMasks;
		setMasks=BBDukIndexConstants.setMasks;
		leftMasks=BBDukIndexConstants.leftMasks;
		rightMasks=BBDukIndexConstants.rightMasks;
		lengthMasks=BBDukIndexConstants.lengthMasks;
		middleMask=BBDukIndexConstants.middleMask;
		minlen2=(middleMask==-1L ? k : (k-midMaskLen)/2);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public final boolean passesSpeed(long key){return speed<1 || ((key&Long.MAX_VALUE)%BBDUK_SPEED_DIVISOR)>=speed;}

	public final boolean failsSpeed(long key){return speed>0 && ((key&Long.MAX_VALUE)%BBDUK_SPEED_DIVISOR)<speed;}
	
	public final boolean isFullyDefined(byte symbol){return symbol>=0 && symbolToNumber[symbol]>=0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Local Fields          ----------------*/
	/*--------------------------------------------------------------*/

	// Kmer/Alphabet/Flags
	protected int k;
	protected int mink;
	protected boolean amino;
	protected int k2;
	protected boolean rcomp;
	protected boolean forbidNs;
	protected int midMaskLen;
	
	// Derived Lengths/Mins
	protected int minlen;
	protected int minlen2;
	
	// Hashing Constants
	protected int bitsPerBase;
	protected int maxSymbol;
	protected int symbols;
	protected int shift;
	protected int shift2;
	protected long mask;
	protected long symbolMask;
	
	// Lookup/Manipulation Arrays
	protected byte[] symbolToNumber;
	protected byte[] symbolToNumber0;
	protected byte[] symbolToComplementNumber0;
	protected long[] clearMasks;
	protected long[][] setMasks;
	protected long[] leftMasks;
	protected long[] rightMasks;
	protected long[] lengthMasks;
	protected long middleMask;
}