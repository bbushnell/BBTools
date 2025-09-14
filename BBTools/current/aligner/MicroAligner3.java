package aligner;

import prok.GeneCaller;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.LongHashMap;

/**
 * Aligns reads to a small, single sequence reference like PhiX.
 * The reference should not have any duplicate kmers.
 * Alignment is only attempted once, at the first matching kmer.
 * This will generate a match string and return the identity.
 * Mapped reads can be printed as sam output.
 * 
 * @author Brian Bushnell
 * @date November 15, 2024
 *
 */
public class MicroAligner3 implements MicroAligner {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	/** Returns identity */
	public float map(Read r) {
		return map(r, minIdentity);
	}
	
	/** Returns identity */
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
	
	/** Returns identity (approx, in the case of Ns) */
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
	
	private final SingleStateAlignerFlat2 getSSA() {
		if(!shared) {return mySSA;}
		return GeneCaller.getSSA();
	}
	
	private final ByteBuilder getBuffer() {
		if(!shared) {return myBuffer;}
		return buffer();
	}
	
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
