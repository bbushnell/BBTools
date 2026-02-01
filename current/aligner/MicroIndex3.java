package aligner;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import map.LongHashMap;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.Read;
import stream.SamHeader;

/**
 * K-mer index for MicroAligner providing fast sequence mapping.
 * Builds a hash map from reference sequences using canonical k-mers with middle-masking
 * for approximate matching. Optimized for single reference alignment scenarios.
 *
 * @author Brian Bushnell
 * @date November 15, 2024
 */
public class MicroIndex3 {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public MicroIndex3(int k_, int midMaskLen_, String path, boolean setSamStatics) {
		this(k_, midMaskLen_, loadRef(path, setSamStatics));
	}

	/**
	 * Constructs a MicroIndex3 from a Read object.
	 * Initializes k-mer indexing parameters and builds the hash map index.
	 *
	 * @param k_ K-mer length for indexing
	 * @param midMaskLen_ Number of middle bases to mask in k-mers for approximate matching
	 * @param r Reference sequence as Read object
	 */
	public MicroIndex3(int k_, int midMaskLen_, Read r) {
//		this(k_, minIdentity_, ref_.length);
		ref=r.bases;
		refname=r.id;
		k=k_;
		k2=k-1;
		midMaskLen=midMaskLen_;
		middleMask=makeMidMask(k, midMaskLen);
		map=new LongHashMap(ref.length*2);
		
		index();
		assert(map.size()>0) : ref.length+", "+k;
	}
	
	/**
	 * Creates a bit mask that zeros out middle bases of k-mers for fuzzy matching.
	 * Allows k-mers with different middle bases but matching flanking regions to match.
	 *
	 * @param k K-mer length
	 * @param midMaskLen Number of middle bases to mask out
	 * @return Bit mask with middle region zeroed
	 */
	public static long makeMidMask(int k, int midMaskLen) {
		assert(k>midMaskLen+1);
		int bitsPerBase=2;
		int bits=midMaskLen*bitsPerBase;
		int shift=((k-midMaskLen)/2)*bitsPerBase;
		long middleMask=~((~((-1L)<<bits))<<shift);
		return middleMask;
	}
	
	/**
	 * Loads reference sequence from file and optionally initializes SAM header data.
	 * Sets up genome data structures including scaffold names, locations, and lengths
	 * required for SAM output formatting.
	 *
	 * @param path File path to reference sequence
	 * @param setSamStatics Whether to initialize Data class static fields for SAM output
	 * @return Reference sequence as Read object
	 */
	public static Read loadRef(String path, boolean setSamStatics) {
		ArrayList<Read> list=ConcurrentGenericReadInputStream.getReads(1, false, 
				FileFormat.testInput(path, null, false), null, null, null);
//		assert(false && !list.isEmpty() && list.get(0)!=null) : list;
		Read r=list.get(0);
		
		int numChroms=Data.numChroms=1;
		SamHeader.PN="MicroAligner3";
		Data.scaffoldNames=new byte[numChroms+1][][];
		Data.scaffoldLocs=new int[numChroms+1][];
		Data.scaffoldLengths=new int[numChroms+1][];
		Data.chromScaffolds=new int[] {0, 1};
		
		Data.scaffoldNames[1]=new byte[][] {Tools.trimToWhitespace(r.id.getBytes())};
		Data.scaffoldLocs[1]=new int[] {0};
		Data.scaffoldLengths[1]=new int[] {r.length()};
		
		return r;
	}
	
	public void index() {indexRef(k, middleMask, ref, map);}
	
	/**
	 * Builds a k-mer hash map index from reference sequence.
	 * Uses canonical k-mers (max of forward/reverse) with middle masking applied.
	 * Stores position with strand information encoded using MINUS_CODE.
	 *
	 * @param k K-mer length
	 * @param midMask Bit mask for middle base masking
	 * @param ref Reference sequence bases
	 * @param map Hash map to populate with k-mer positions
	 * @return The populated hash map
	 */
	private static LongHashMap indexRef(int k, long midMask, byte[] ref, LongHashMap map) {
		map.clear();
		byte[] bases=ref;
		
		assert(k<=32);
		
		if(bases==null || bases.length<k){return map;}
		final int bitsPerBase=2;
		final int shift=bitsPerBase*k;
		final int shift2=shift-bitsPerBase;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
				if(len>=k){
					long key=Tools.max(kmer, rkmer)&midMask;
					int value=(kmer>=rkmer ? i : i+MINUS_CODE);
//					if(!map.containsKey(key)) {//Not needed, since old values are not replaced
					map.put(key, value);
//					}
				}
			}
		}
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mapping            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Maps a read to the reference using k-mer lookup.
	 * Returns first matching position found by scanning read k-mers against the index.
	 * Sets read alignment fields (chrom, strand, start, stop) if mapping succeeds.
	 *
	 * @param r Read to map against the reference
	 * @return Encoded position and strand ((offset<<1)|strand) or NO_HIT if no match
	 */
	public long map(Read r) {
		assert(!r.mapped() && r.match==null && r.samline==null);
		r.chrom=-1;
		if(r==null || r.length()<k || r.match!=null || r.samline!=null) {return 0;}
		byte[] bases=r.bases;
		if(bases==null || bases.length<k){return 0;}
		final int bitsPerBase=2;
		final int shift=bitsPerBase*k;
		final int shift2=shift-bitsPerBase;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		int offset=-1;
		int orientation=-1;
		int ivalue=-1;
		int vvalue=-1;
		
		//TODO: This loop could be replaced by the kmer list which already exists.
		for(int i=0; i<bases.length && orientation<0; i++){
			final byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
				if(len>=k /*&& ((i&1)==0)*/){
					long key=Tools.max(kmer, rkmer)&middleMask;
					int value=map.get(key);
					if(value>=0) {
						ivalue=i;
						vvalue=value;
						
						if(value>=MINUS_CODE) {
							value-=MINUS_CODE;
							if(kmer>=rkmer) {
								orientation=1;
								offset=value-k2-(bases.length-i-1);
							}else {
//								assert(key==rkmer);
								orientation=2;
								offset=value-i;
							}
						}else {
							if(kmer>=rkmer) {
								orientation=0;
								offset=value-i;
							}else {
//								assert(key==rkmer);
								orientation=3;
								offset=value-k2-(bases.length-i-1);
							}
						}
					}
				}
			}
		}
//		assert(map.size()>0) : orientation+", "+offset+", "+map.size();
		if(orientation<0) {return NO_HIT;}
		r.chrom=1;
		int strand=orientation&1;
		r.setStrand(strand);
		r.start=offset;
		r.stop=offset+r.length()-1;
//		r.setMapped(true);
		return (offset<<1)|strand;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public LongHashMap getMap() {return map;}
	
	final int k;
	final int k2;
	final int midMaskLen;
	final long middleMask;
	final private String refname;
	final byte[] ref;
	final LongHashMap map;
	
	//Indicates the position is on the minus strand
	static final int MINUS_CODE=1000000000;
	static final int NO_HIT=Integer.MIN_VALUE;
}
