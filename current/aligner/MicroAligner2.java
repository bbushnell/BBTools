package aligner;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import prok.GeneCaller;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.Read;
import structures.LongHashMap;

/**
 * Similar to MicroAligner, but designed to only test whether something matches,
 * not generate match strings.
 * 
 * @author Brian Bushnell
 * @date May 24, 2024
 *
 */
public class MicroAligner2 implements MicroAligner {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public MicroAligner2(int k_, float minIdentity_, String path) {
		this(k_, minIdentity_, loadRef(path));
	}

	public MicroAligner2(int k_, float minIdentity_, byte[] ref_) {
		this(k_, minIdentity_, ref_, indexRef(k_, ref_));
	}

	public MicroAligner2(int k_, float minIdentity_, byte[] ref_, LongHashMap map_) {
		k=k_;
		k2=k-1;
		ref=ref_;
		map=map_;
		minIdentity=minIdentity_;
		maxSubFraction=1-minIdentity;
	}
	
	private static byte[] loadRef(String path) {
		if(path==null) {return null;}
		if(!Tools.isReadableFile(path)) {
			if(path.equalsIgnoreCase("phix")) {
				path=Data.findPath("?phix2.fa.gz");
			}
		}
		ArrayList<Read> list=ConcurrentGenericReadInputStream.getReads(1, false, 
				FileFormat.testInput(path, null, false), null, null, null);
		return list.get(0).bases;
	}
	
	//TODO: Add MaskMiddle and/or hdist
	private static LongHashMap indexRef(int k, byte[] ref) {
		LongHashMap map=new LongHashMap(ref.length*2);
		byte[] bases=ref;
		
		assert(k<=32);
		
		if(bases==null || bases.length<k){return null;}
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
					long key=Tools.max(kmer, rkmer);
					int value=(key==kmer ? i : i+MINUS_CODE);
					map.put(key, value);
				}
			}
		}
		return map;
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
		if(r==null || r.length()<k) {return 0;}
		mapCount++;
		byte[] bases=r.bases;
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
				if(len>=k && ((i&skipmask)==0)){
					long key=Tools.max(kmer, rkmer);
					int value=map.get(key);
					if(value>=0) {
						ivalue=i;
						vvalue=value;
						
						if(value>=MINUS_CODE) {
							value-=MINUS_CODE;
							if(key==kmer) {
								orientation=1;
								offset=value-k2-(bases.length-i-1);
							}else {
								assert(key==rkmer);
								orientation=2;
								offset=value-i;
							}
						}else {
							if(key==kmer) {
								orientation=0;
								offset=value-i;
							}else {
								assert(key==rkmer);
								orientation=3;
								offset=value-k2-(bases.length-i-1);
							}
						}
					}
				}
			}
		}
		
		if(orientation<0) {return 0;}
		final float id;
		int pad=5;
		if(orientation==1 || orientation==3) {
			r.reverseComplement();
			id=align(r, ref, offset, offset+r.length(), pad, minIdentity, null);
			r.reverseComplement();
		}else {
			id=align(r, ref, offset, offset+r.length(), pad, minIdentity, null);
		}
		metCutoff+=(id>=minIdentity ? 1 : 0);
		idSum+=(id>=minIdentity ? id : 0);
		return id;
	}
	
	public float align(Read r, byte[] ref, int a, int b, int pad, float minIdentity, int[] extra){
		float id=quickAlign(r, ref, a);
		if(id>=minIdentity) {return id;}
		
		slowAligns++;
		SingleStateAlignerFlat3 ssa=GeneCaller.getSSA3();
		a=Tools.max(0, a-pad);
		b=Tools.min(ref.length-1, b+pad);
		int[] max=ssa.fillUnlimited(r.bases, ref, a, b, 0);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		
		id=ssa.tracebackIdentity(r.bases, ref, a, b, rows, maxCol, maxState, extra);
		return id;
	}
	
	public float quickAlign(Read read, byte[] ref, int a) {
		quickAligns++;
		byte[] bases=read.bases;
		int subs=0, ns=0, matches=0;
		int maxSubs=(int)(maxSubFraction*read.length());
		//This starts it at the lowest index preventing either from going off the left end
		for(int i=Tools.max(0, -a), j=Tools.max(0, a); 
				i<bases.length && j<ref.length && subs<maxSubs; i++, j++) {
				final byte q=bases[i], r=ref[j];
				final boolean qr=(q==r), qn=(q=='N'), rn=(r=='N');
				matches+=(qr && !qn) ? 1 :0;
				subs+=(!qr && !qn && !rn) ? 1 :0;
				ns+=(qn || rn) ? 1 :0;
		}
		return matches/(float)read.length();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final float minIdentity;
	final float maxSubFraction;
	final int k;
	final int k2;
	final byte[] ref;
	final LongHashMap map;
	public int skipmask=0;

	public long mapCount=0;
	public long quickAligns=0;
	public long slowAligns=0;
	public long metCutoff=0;
	public double idSum=0;
	
	//Indicates the position is on the minus strand
	private static final int MINUS_CODE=1000000000;
}
