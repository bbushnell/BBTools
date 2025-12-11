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
 * Reads and processes k-mer tables for sequence matching and masking.
 * Performs k-mer based read processing, including matching, counting, and masking
 * operations with support for forward and reverse complement k-mer matching.
 *
 * @author Brian Bushnell
 * @date Mar 5, 2015
 */
public class TableReader {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line for testing TableReader functionality.
	 * Initializes k-mer tables, loads reference data, and demonstrates basic usage.
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
	 * Masks bases in a read that correspond to matching k-mers in the reference tables.
	 * Identified bases are replaced with the trim symbol or converted to lowercase.
	 *
	 * @param r Read to process and mask
	 * @param sets K-mer hash tables containing reference k-mers
	 * @return Number of bases masked in the read
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
	 * Counts the number of k-mer matches between a read and reference k-mer tables.
	 * Scans through read k-mers and counts hits against the provided tables.
	 *
	 * @param r Read to analyze for k-mer matches
	 * @param sets K-mer hash tables to search against
	 * @return Number of k-mer hits found
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
	 * Finds the reference sequence ID with the most k-mer matches to the given read.
	 * Returns the ID of the best matching sequence, or -1 if no matches meet the minimum threshold.
	 *
	 * @param r Read to find matches for
	 * @param sets K-mer hash tables containing reference sequences
	 * @return ID of best matching reference sequence, or -1 if insufficient matches
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
	 * Creates a BitSet marking all bases in the read that should be masked based on k-mer matches.
	 * Handles both normal k-mers and short k-mers at read ends when configured.
	 *
	 * @param r Read to analyze and mark
	 * @param sets K-mer hash tables to search against
	 * @return BitSet with positions to mask, or null if no matches found
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
	 * Retrieves the stored value for a k-mer from hash tables, considering query position skipping.
	 * Delegates to Hamming distance-aware search if configured.
	 *
	 * @param kmer Forward k-mer sequence
	 * @param rkmer Reverse complement k-mer sequence
	 * @param len Length of the k-mer
	 * @param qHDist Maximum Hamming distance for fuzzy matching
	 * @param qPos Position of k-mer in query sequence
	 * @param sets K-mer hash tables to search
	 * @return Stored value for the k-mer, or -1 if not found
	 */
	public final int getValue(final long kmer, final long rkmer, final int len, final int qHDist, final int qPos, final AbstractKmerTable[] sets){
		if(qSkip>1 && (qPos%qSkip!=0)){return -1;}
		return qHDist<1 ? getValue(kmer, rkmer, len, sets) : getValue(kmer, rkmer, len, qHDist, sets);
	}
	
	/**
	 * Searches for k-mer value with Hamming distance tolerance.
	 * First searches for exact match, then generates variants with single substitutions if needed.
	 *
	 * @param kmer Forward k-mer sequence
	 * @param rkmer Reverse complement k-mer sequence
	 * @param len Length of the k-mer
	 * @param qHDist Maximum Hamming distance for fuzzy matching
	 * @param sets K-mer hash tables to search
	 * @return Stored value for the k-mer or variant, or -1 if not found
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
	 * Retrieves stored value for a k-mer from hash tables using exact matching.
	 *
	 * @param kmer Forward k-mer sequence
	 * @param rkmer Reverse complement k-mer sequence
	 * @param len Length of the k-mer
	 * @param sets K-mer hash tables to search
	 * @return Stored value for the k-mer, or -1 if not found
	 */
	public final int getValue(final long kmer, final long rkmer, final int len, final AbstractKmerTable[] sets){
		return getValueWithMask(kmer, rkmer, lengthMasks[len], sets);
	}
	
	/**
	 * Searches hash tables for a canonical k-mer value using the provided length mask.
	 * Applies speed filtering and middle masking before table lookup.
	 *
	 * @param kmer Forward k-mer sequence
	 * @param rkmer Reverse complement k-mer sequence
	 * @param lengthMask Bitmask indicating k-mer length
	 * @param sets K-mer hash tables to search
	 * @return Stored value for the k-mer, or -1 if not found
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
	 * Converts k-mer pair to canonical representation for hash table storage.
	 * Selects larger of forward/reverse k-mer, applies middle masking, and adds length marker.
	 *
	 * @param kmer Forward k-mer sequence
	 * @param rkmer Reverse complement k-mer sequence
	 * @param lengthMask Bitmask marking k-mer length
	 * @return Canonical k-mer value for hash table lookup
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	
	/**
	 * Transfers counts from a sparse array to compact IntLists and resets the array.
	 * Used to consolidate k-mer hit counts from array-based tracking.
	 *
	 * @param loose Sparse counter array to read from
	 * @param packed List of unique sequence IDs
	 * @param counts List to store corresponding counts
	 * @return Maximum count observed
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
	
	/** Tracks whether this TableReader has encountered processing errors */
	public boolean errorState=false;
	
	/** Whether to mask middle bases in k-mers for fuzzy matching (always false) */
	public final boolean maskMiddle=false;
	
	/** Maximum Hamming distance allowed for query k-mer matching */
	private final int qHammingDistance;
	/** Maximum Hamming distance for short k-mer matching at read ends */
	public int qHammingDistance2=-1;
	
	/** Additional bases to mask around k-mer matches */
	public int trimPad=0;
	
	/** If positive, only search for k-mers in leftmost X bases of read */
	public int restrictLeft=0;
	/** If positive, only search for k-mers in rightmost X bases of read */
	public int restrictRight=0;
	
	/** Whether to prevent ambiguous bases ('N') from matching reference bases */
	public boolean forbidNs=false;
	
	/** Character used to replace bases covered by matched k-mers */
	public byte trimSymbol='N';
	
	/**
	 * Whether to convert masked bases to lowercase instead of replacing with trimSymbol
	 */
	public boolean kmaskLowercase=false;
	
	/** Whether to skip k-mer processing for read 1 in paired reads */
	public boolean skipR1=false;
	/** Whether to skip k-mer processing for read 2 in paired reads */
	public boolean skipR2=false;

	/** Minimum number of k-mer hits required to consider a read a match */
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
	
	/**
	 * Whether to search for reverse complement k-mers in addition to forward k-mers
	 */
	private final boolean rcomp;
	/** Bitmask with zeros at middle base positions for fuzzy k-mer matching */
	private final long middleMask;
	
	/** Primary k-mer length for sequence processing */
	private final int k;
	/** K-mer length minus 1, used in various calculations */
	private final int k2;
	/** Minimum k-mer length for short k-mer matching at read ends */
	private final int mink;
	/** Whether to attempt matching k-mers shorter than k at read ends */
	private final boolean useShortKmers;
	
	/** Speed setting (0-15) controlling fraction of k-mers to process */
	private final int speed;
	
	/**
	 * Number of k-mers to skip when scanning reads (1=every k-mer, 2=every other, etc.)
	 */
	private final int qSkip;
	
	/** True when speed and qSkip acceleration are disabled */
	private final boolean noAccel, accel;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of hash table partitions for parallel processing */
	private static final int WAYS=7; //123
	/** Whether to print verbose debugging messages */
	public static final boolean verbose=false; //123
	
	/** Output stream for printing messages and results */
	private static PrintStream outstream=System.err;
	
	/** Bitmasks for clearing individual bases in k-mers (position i) */
	private static final long[] clearMasks;
	/** Bitmasks for setting base i to nucleotide j in k-mer manipulation */
	private static final long[][] setMasks;
	/** Bitmasks for clearing all bases to the right of position i */
	private static final long[] leftMasks;
	/** Bitmasks for clearing all bases to the left of position i */
	private static final long[] rightMasks;
	/** Bitmasks marking k-mer lengths with single bit set to left of k-mer */
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
