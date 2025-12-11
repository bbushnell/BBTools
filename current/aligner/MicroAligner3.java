package aligner;

import prok.GeneCaller;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.LongHashMap;

/**
 * Aligns reads to a small, single sequence reference like PhiX.
 * Uses k-mer lookup for initial positioning, then applies alignment at that position.
 * Designed for references without duplicate k-mers, performing alignment only once
 * at the first matching k-mer location.
 *
 * @author Brian Bushnell
 * @date November 15, 2024
 */
public class MicroAligner3 implements MicroAligner {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a MicroAligner3 with specified index and alignment parameters.
	 * Initializes internal state including k-mer parameters and thread-local resources.
	 *
	 * @param index_ The k-mer index for the reference sequence
	 * @param minIdentity_ Minimum identity threshold for accepting alignments
	 * @param shared_ Whether to use shared thread-local resources or private instances
	 */
	public MicroAligner3(MicroIndex3 index_, float minIdentity_, boolean shared_) {
		index=index_;
		minIdentity=minIdentity_;
		maxSubFraction=1-minIdentity;
		k=index.k;
		k2=k-1;
		middleMask=index.middleMask;
		shared=shared_;
		myBuffer=(shared ? null : new ByteBuilder());
		mySSA=(shared ? null : new SingleStateAlignerFlat2());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Alignment           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Maps a read to the reference using the default minimum identity threshold.
	 * @param r The read to align
	 * @return Identity score of the alignment, or 0 if no mapping found
	 */
	public float map(Read r) {
		return map(r, minIdentity);
	}
	
	/**
	 * Maps a read to the reference with specified minimum identity threshold.
	 * Uses k-mer lookup to find initial position, then performs alignment.
	 * Handles both forward and reverse complement orientations.
	 *
	 * @param r The read to align
	 * @param minid Minimum identity threshold for accepting the alignment
	 * @return Identity score of the alignment, or 0 if no mapping found above threshold
	 */
	public float map(Read r, float minid) {
		if(r==null || r.length()<k || r.match!=null || r.samline!=null) {return 0;}
		
		final long ret=index.map(r);
		if(ret==index.NO_HIT) {
			assert(!r.mapped());
			return 0;
		}
		int strand=(int)(ret&1);
		int offset=(int)(ret>>1);
		assert(offset==r.start);
		
		final float id;
		int pad=5;
		if(strand==1) {
			r.reverseComplement();
			id=align(r, index.ref, offset, offset+r.length(), pad, minid);
			r.reverseComplement();
			if(r.mapped()) {r.setStrand(1);}
		}else {
			id=align(r, index.ref, offset, offset+r.length(), pad, minid);
		}
		assert(id>=minid || !r.mapped()) : "\nid="+id+"<"+minid+"\n"+r+"\n"+r.mate.toFastq()+"\n";
		return id;
	}
	
	/**
	 * Performs alignment of read to reference sequence within specified bounds.
	 * First attempts quick alignment, then falls back to full dynamic programming
	 * alignment if quick alignment fails to meet identity threshold.
	 *
	 * @param r The read to align
	 * @param ref Reference sequence bytes
	 * @param a Start position in reference
	 * @param b End position in reference
	 * @param pad Padding bases to extend alignment window
	 * @param minid Minimum identity threshold
	 * @return Identity score of the alignment
	 */
	public float align(Read r, byte[] ref, int a, int b, int pad, float minid){
		assert(!r.mapped());
		{
			final float id=quickAlign(r, ref, a, minid);
			if(id>minid) {
				assert(r.mapped());
				return id;
			}
		}
		assert(!r.mapped());
		
		SingleStateAlignerFlat2 ssa=getSSA();
		a=Tools.max(0, a-pad);
		b=Tools.min(ref.length-1, b+pad);
		int[] max=ssa.fillUnlimited(r.bases, ref, a, b, 0);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		int[] score=ssa.score(r.bases, ref, a, b, rows, maxCol, maxState);
		int rstart=score[1];
		int rstop=score[2];
		r.start=rstart;
		r.stop=rstop;
		r.chrom=1;
		
		final byte[] match=ssa.traceback(r.bases, ref, a, b, rows, maxCol, maxState);
		final float id=Read.identity(match);
		
		if(id<minid) {return id;}//Probably not phix
		r.setMapped(true);
		r.match=match;

		assert(id>=minid || !r.mapped()) : "\nid="+id+"<"+minid+"\n"+r+"\n"+r.mate.toFastq()+"\n";
		return id;
	}
	
	/**
	 * Performs fast approximate alignment by direct base comparison.
	 * Builds match string and calculates identity without dynamic programming.
	 * Used as first-pass alignment attempt before full alignment.
	 *
	 * @param read The read to align
	 * @param ref Reference sequence bytes
	 * @param a Start position in reference for alignment
	 * @param minid Minimum identity threshold
	 * @return Approximate identity score, accounting for N bases and substitutions
	 */
	public float quickAlign(Read read, byte[] ref, int a, float minid) {
		byte[] bases=read.bases;
		ByteBuilder buffer=getBuffer();
		buffer.clear();
		int subs=0, ns=0;
		for(int i=0, j=a; i<bases.length; i++, j++) {
			if(j<0 || j>=ref.length) {
				buffer.append('C');
			}else {
				final byte q=bases[i], r=ref[j];
				if(q=='N') {
					buffer.append('N');
					ns++;
				}else if(r=='N' || r==q) {
					buffer.append('m');
				}else {
					buffer.append('S');
					subs++;
				}
			}
		}
		int matches=bases.length-subs-ns;
		if(subs>3 || matches*4<bases.length) {return 0;}
		float id=(subs+ns*0.0625f)/Tools.max(1f, matches+0.25f*ns+subs);
		if(id>=minid) {
			read.match=buffer.toBytes();
			read.start=a;
			read.stop=a+read.length()-1;
			read.chrom=1;
			read.setMapped(true);
		}
		return id;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Gets SingleStateAligner instance for dynamic programming alignment.
	 * Returns private instance if not shared, otherwise uses shared thread-local instance.
	 * @return SingleStateAlignerFlat2 instance for alignment operations
	 */
	private final SingleStateAlignerFlat2 getSSA() {
		if(!shared) {return mySSA;}
		return GeneCaller.getSSA();
	}
	
	/**
	 * Gets ByteBuilder instance for building match strings.
	 * Returns private instance if not shared, otherwise uses shared thread-local instance.
	 * @return ByteBuilder instance for constructing alignment strings
	 */
	private final ByteBuilder getBuffer() {
		if(!shared) {return myBuffer;}
		return buffer();
	}
	
	/**
	 * Gets thread-local ByteBuilder instance for shared usage.
	 * Creates new instance if none exists for current thread.
	 * @return Thread-local ByteBuilder instance
	 */
	private static final ByteBuilder buffer() {
		ByteBuilder buffer=bufferHolder.get();
		if(buffer==null) {
			buffer=new ByteBuilder();
			bufferHolder.set(buffer);
		}
		return buffer;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public LongHashMap getMap() {return index.map;}
	
	final float minIdentity;
	final float maxSubFraction;
	final int k;
	final int k2;
	final long middleMask;
	final MicroIndex3 index;
	final SingleStateAlignerFlat2 mySSA;
	final ByteBuilder myBuffer;
	final boolean shared;
	
	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/

	private static final float nMult=1024;
	private static final float nMultInv=1.0f/nMult;
	private static final ThreadLocal<ByteBuilder> bufferHolder=new ThreadLocal<ByteBuilder>();
//	final ByteBuilder buffer=new ByteBuilder();
	
}
