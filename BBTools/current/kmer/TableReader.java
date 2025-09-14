package kmer;

import java.io.PrintStream;
import java.util.BitSet;

import dna.AminoAcid;
import jgi.Dedupe;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date Mar 5, 2015
 *
 */
public class TableReader {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		AbstractKmerTable[] tables=TableLoaderLockFree.makeTables(AbstractKmerTable.ARRAY1D, 12, -1L, false, 1.0);
		
		int k=31;
		int mink=0;
		int speed=0;
		int hdist=0;
		int edist=0;
		boolean rcomp=true;
		boolean maskMiddle=false;
		
		//Create a new Loader instance
		TableLoaderLockFree loader=new TableLoaderLockFree(tables, k, mink, speed, hdist, edist, rcomp, maskMiddle);
		loader.setRefSkip(0);
		loader.hammingDistance2=0;
		loader.editDistance2=0;
		loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
		
		///And run it
		String[] refs=args;
		String[] literals=null;
		boolean keepNames=false;
		boolean useRefNames=false;
		long kmers=loader.processData(refs, literals, keepNames, useRefNames, false);
		t.stop();

		outstream.println("Load Time:\t"+t);
		outstream.println("Return:   \t"+kmers);
		outstream.println("refKmers: \t"+loader.refKmers);
		outstream.println("refBases: \t"+loader.refBases);
		outstream.println("refReads: \t"+loader.refReads);
		
		int qskip=0;
		int qhdist=0;
		TableReader tr=new TableReader(k, mink, speed, qskip, qhdist, rcomp, maskMiddle);
		
		//TODO: Stuff...
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public TableReader(int k_){
		this(k_, 0, 0, 0, 0, true, false);
	}
	
	public TableReader(int k_, int mink_, int speed_, int qskip_, int qhdist_, boolean rcomp_, boolean maskMiddle_){
		k=k_;
		k2=k-1;
		mink=mink_;
		rcomp=rcomp_;
		useShortKmers=(mink>0 && mink<k);
		speed=speed_;
		qSkip=qskip_;
		qHammingDistance=qhdist_;
		middleMask=maskMiddle ? ~(3L<<(2*(k/2))) : -1L;
		
		noAccel=(speed<1 && qSkip<2);
		accel=!noAccel;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Mask a read to cover matching kmers.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases masked
	 */
	public final int kMask(final Read r, final AbstractKmerTable[] sets){
		if(r==null){return 0;}
		if(verbose){outstream.println("KMasking read "+r.id);}
		
		BitSet bs=markBits(r, sets);
		if(verbose){outstream.println("Null bitset.");}
		if(bs==null){return 0;}

		final byte[] bases=r.bases, quals=r.quality;
		final int cardinality=bs.cardinality();
		assert(cardinality>0);
		
		//Replace kmer hit zone with the trim symbol
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				if(kmaskLowercase){
					bases[i]=(byte)Tools.toLowerCase(bases[i]);
				}else{
					bases[i]=trimSymbol;
					if(quals!=null && trimSymbol=='N'){quals[i]=0;}
				}
			}
		}
		return cardinality;
	}
	
	
	/**
	 * Counts the number of kmer hits for a read.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of hits
	 */
	public final int countKmerHits(final Read r, final AbstractKmerTable[] sets){
		if(r==null || r.length()<k){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		final byte[] bases=r.bases;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		
		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber0[b];
			long x2=AminoAcid.baseToComplementNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(b=='N' && forbidNs){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=getValue(kmer, rkmer, k, qHammingDistance, i, sets);
				if(verbose){outstream.println("Testing kmer "+kmer+"; id="+id);}
				if(id>0){
					if(verbose){outstream.println("Found = "+(found+1)+"/"+minHits);}
					if(found>=minHits){
						return (found=found+1); //Early exit
					}
					found++;
				}
			}
		}
		
		return found;
	}
	
	/**
	 * Returns the id of the sequence with the most kmer matches to this read, or -1 if none are at least minHits.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return id of best match
	 */
	public final int findBestMatch(final Read r, final AbstractKmerTable[] sets){
		idList.size=0;
		if(r==null || r.length()<k){return -1;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return -1;}
		final byte[] bases=r.bases;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		int len=0;
		int found=0;
		
		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber0[b];
			long x2=AminoAcid.baseToComplementNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(b=='N' && forbidNs){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=getValue(kmer, rkmer, k, qHammingDistance, i, sets);
				if(id>0){
					countArray[id]++;
					if(countArray[id]==1){idList.add(id);}
					found++;
					if(verbose){outstream.println("Found = "+found+"/"+minHits);}
				}
			}
		}
		
		final int id, max;
		if(found>=minHits){
			max=condenseLoose(countArray, idList, countList);
			int id0=-1;
			for(int i=0; i<countList.size; i++){
				if(countList.get(i)==max){
					id0=idList.get(i); break;
				}
			}
			id=id0;
		}else{
			max=0;
			id=-1;
		}
		
		return id;
	}
	
	
	/**
	 * Mask a read to cover matching kmers.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases masked
	 */
	public final BitSet markBits(final Read r, final AbstractKmerTable[] sets){
		if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k))){
			if(verbose){outstream.println("Read too short.");}
			return null;
		}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){
			if(verbose){outstream.println("Skipping read.");}
			return null;
		}
		if(verbose){outstream.println("Marking bitset for read "+r.id);}
		final byte[] bases=r.bases;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		int id0=-1; //ID of first kmer found.
		
		BitSet bs=new BitSet(bases.length+trimPad+1);
		
		final int minus=k-1-trimPad;
		final int plus=trimPad+1;
		
		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));
		
		//Scan for normal kmers
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber0[b];
			long x2=AminoAcid.baseToComplementNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(b=='N' && forbidNs){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=getValue(kmer, rkmer, k, qHammingDistance, i, sets);
				if(id>0){
					if(id0<0){id0=id;}
					if(verbose){
						outstream.println("a: Found "+kmer);
						outstream.println("Setting "+Tools.max(0, i-minus)+", "+(i+plus));
						outstream.println("i="+i+", minus="+minus+", plus="+plus+", trimpad="+trimPad+", k="+k);
					}
					bs.set(Tools.max(0, i-minus), i+plus);
					found++;
				}
			}
		}
		
		//If nothing was found, scan for short kmers.
		if(useShortKmers){
			assert(!maskMiddle && middleMask==-1) : maskMiddle+", "+middleMask+", k="+", mink="+mink;
			
			//Look for short kmers on left side
			{
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.min(k, stop);
				for(int i=start; i<lim; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=rkmer|(x2<<(2*len));
					len++;
					if(verbose){outstream.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						
						if(verbose){
							outstream.println("Looking for left kmer  "+AminoAcid.kmerToString(kmer, len));
							outstream.println("Looking for left rkmer "+AminoAcid.kmerToString(rkmer, len));
						}
						final int id=getValue(kmer, rkmer, len, qHammingDistance2, i, sets);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("b: Found "+kmer);
								outstream.println("Setting "+0+", "+(i+plus));
							}
							bs.set(0, i+plus);
							found++;
						}
					}
				}
			}

			//Look for short kmers on right side
			{
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.max(-1, stop-k);
				for(int i=stop-1; i>lim; i--){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=kmer|(x<<(2*len));
					rkmer=((rkmer<<2)|x2)&mask;
					len++;
					if(verbose){outstream.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						if(verbose){
							outstream.println("Looking for right kmer "+
									AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+"; value="+toValue(kmer, rkmer, lengthMasks[len])+"; kmask="+lengthMasks[len]);
						}
						final int id=getValue(kmer, rkmer, len, qHammingDistance2, i, sets);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("c: Found "+kmer);
								outstream.println("Setting "+Tools.max(0, i-trimPad)+", "+bases.length);
							}
							bs.set(Tools.max(0, i-trimPad), bases.length);
							found++;
						}
					}
				}
			}
		}
		
		
		if(verbose){outstream.println("found="+found+", bitset="+bs);}
		
		if(found==0){return null;}
		assert(found>0) : "Overflow in 'found' variable.";
		
		int cardinality=bs.cardinality();
		assert(cardinality>0);
		
		return bs;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	/**
	 * Transforms a kmer into all canonical values for a given Hamming distance.
	 * Returns the related id stored in the tables.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len kmer length
	 * @param qHDist Hamming distance
	 * @param qPos Position of kmer in query
	 * @param sets Kmer hash tables
	 * @return Value stored in table, or -1
	 */
	public final int getValue(final long kmer, final long rkmer, final int len, final int qHDist, final int qPos, final AbstractKmerTable[] sets){
		if(qSkip>1 && (qPos%qSkip!=0)){return -1;}
		return qHDist<1 ? getValue(kmer, rkmer, len, sets) : getValue(kmer, rkmer, len, qHDist, sets);
	}
	
	/**
	 * Transforms a kmer into all canonical values for a given Hamming distance.
	 * Returns the related id stored in the tables.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len kmer length
	 * @param qHDist Hamming distance
	 * @param sets Kmer hash tables
	 * @return Value stored in table, or -1
	 */
	public final int getValue(final long kmer, final long rkmer, final int len, final int qHDist, final AbstractKmerTable[] sets){
		int id=getValue(kmer, rkmer, len, sets);
		if(id<1 && qHDist>0){
			final int qHDistMinusOne=qHDist-1;
			
			//Sub
			for(int j=0; j<4 && id<1; j++){
				for(int i=0; i<len && id<1; i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
						id=getValue(temp, rtemp, len, qHDistMinusOne, sets);
					}
				}
			}
		}
		return id;
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table and search.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param len kmer length
	 * @param sets Kmer hash tables
	 * @return Value stored in table
	 */
	public final int getValue(final long kmer, final long rkmer, final int len, final AbstractKmerTable[] sets){
		return getValueWithMask(kmer, rkmer, lengthMasks[len], sets);
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table and search.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @param sets Kmer hash tables
	 * @return Value stored in table
	 */
	public final int getValueWithMask(final long kmer, final long rkmer, final long lengthMask, final AbstractKmerTable[] sets){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		
//		final long max=(rcomp ? Tools.max(kmer, rkmer) : kmer);
//		final long key=(max&middleMask)|lengthMask;
		
		final long key=toValue(kmer, rkmer, lengthMask);
		
		if(noAccel || ((key/WAYS)&15)>=speed){
			if(verbose){outstream.println("Testing key "+key);}
			AbstractKmerTable set=sets[(int)(key%WAYS)];
			final int id=set.getValue(key);
			return id;
		}
		return -1;
	}
	
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	
	/**
	 * Pack a list of counts from an array to an IntList.
	 * @param loose Counter array
	 * @param packed Unique values
	 * @param counts Counts of values
	 * @return Highest observed count
	 */
	public static int condenseLoose(int[] loose, IntList packed, IntList counts){
		counts.size=0;
		if(packed.size<1){return 0;}

		int max=0;
		for(int i=0; i<packed.size; i++){
			final int p=packed.get(i);
			final int c=loose[p];
			counts.add(c);
			loose[p]=0;
			max=Tools.max(max, c);
		}
		return max;
	}
	
	public final int kmerToWay(final long kmer){
//		final int way=(int)((kmer&coreMask)%WAYS);
//		return way;
		return (int)(kmer%WAYS);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	public final boolean maskMiddle=false;
	
	/** Search for query kmers with up to this many substitutions */
	private final int qHammingDistance;
	/** Search for short query kmers with up to this many substitutions */
	public int qHammingDistance2=-1;
	
	/** Trim this much extra around matched kmers */
	public int trimPad=0;
	
	/** If positive, only look for kmer matches in the leftmost X bases */
	public int restrictLeft=0;
	/** If positive, only look for kmer matches the rightmost X bases */
	public int restrictRight=0;
	
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	public boolean forbidNs=false;
	
	/** Replace bases covered by matched kmers with this symbol */
	public byte trimSymbol='N';
	
	/** Convert masked bases to lowercase */
	public boolean kmaskLowercase=false;
	
	/** Don't look for kmers in read 1 */
	public boolean skipR1=false;
	/** Don't look for kmers in read 2 */
	public boolean skipR2=false;

	/** A read must contain at least this many kmer hits before being considered a match.  Default: 1 */
	public int minHits=1;
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
//	public long storedKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------      Per-Thread Fields       ----------------*/
	/*--------------------------------------------------------------*/
	
	public int[] countArray;
	
	private final IntList idList=new IntList();
	private final IntList countList=new IntList();
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Shortest kmer to use for trimming */
	private final int mink;
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private final boolean useShortKmers;
	
	/** Fraction of kmers to skip, 0 to 15 out of 16 */
	private final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int qSkip;
	
	/** noAccel is true if speed and qSkip are disabled, accel is the opposite. */
	private final boolean noAccel, accel;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	
	/** x&clearMasks[i] will clear base i */
	private static final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	private static final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private static final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] lengthMasks;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
	static{
		clearMasks=new long[32];
		leftMasks=new long[32];
		rightMasks=new long[32];
		lengthMasks=new long[32];
		setMasks=new long[4][32];
		for(int i=0; i<32; i++){
			clearMasks[i]=~(3L<<(2*i));
		}
		for(int i=0; i<32; i++){
			leftMasks[i]=((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			rightMasks[i]=~((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			lengthMasks[i]=((1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			for(long j=0; j<4; j++){
				setMasks[(int)j][i]=(j<<(2*i));
			}
		}
	}
	
}
