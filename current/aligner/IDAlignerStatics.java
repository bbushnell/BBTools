package aligner;

import shared.Shared;

public class IDAlignerStatics{

	/**
	 * Counts substitutions up to maxBandwidth.
	 * Updates pos[0] to optimal offset if SIMD is used.
	 */
	public static int countSubs(byte[] query, byte[] ref, int[] pos, int maxBandwidth){
		assert(pos!=null && pos.length>=2);
		if(Shared.SIMD){
			return simd.SIMDAlignByte.countSubs(query, ref, pos, maxBandwidth);
		}
		int offset=pos[0];
		if(offset==Integer.MIN_VALUE){offset=0; pos[0]=0;}
		int subs=0, qLen=query.length, rLen=ref.length;
		
		for(int i=0; i<qLen && subs<maxBandwidth; i++){
			int rIdx=i+offset;
			if(rIdx>=0 && rIdx<rLen){
				if(query[i]!=ref[rIdx]){subs++;}
			}else{
				subs++;
			}
		}
		return subs;
	}
	
	/**
	 * Counts substitutions up to maxBandwidth.
	 * Returns min(subs+1, maxBandwidth);
	 * Updates pos[0] to optimal offset if SIMD is used.
	 */
	public static int decideBandwidth(final byte[] query, final byte[] ref, 
			final int[] pos, final int maxBandwidth){
		assert(pos!=null && pos.length>=2);
		if(Shared.SIMD){
			final int subs=simd.SIMDAlignByte.countSubs(query, ref, pos, maxBandwidth);
			return Math.min(subs+1, maxBandwidth);
		}
		
		final int qLen=query.length, rLen=ref.length;
		int offset=pos[0];
		if(offset==Integer.MIN_VALUE){offset=0; pos[0]=0;}
		
		// 1. Calculate Intersection Bounds (Relative to Query)
		// Left Overhang: Query indices 0 to (-offset - 1) are blind
		final int qStart=Math.max(0, -offset);
		// Right Overhang: Query indices beyond (rLen - offset) are blind
		final int qEnd=Math.min(qLen, rLen-offset);

		// 2. Pre-calculate Overhang Penalties
		// If qStart >= qEnd, there is no intersection (pure overhang)
		if(qStart>=qEnd){return Math.min(qLen+1, maxBandwidth);}
		
		int subs=qStart + Math.max(0, qLen-qEnd);
		if(subs>=maxBandwidth){return maxBandwidth;} 

		// 3. The "Hot" Loop - No Branching, No Bounds Checks!
		// We only iterate where BOTH query and ref exist.
		for(int q=qStart; q<qEnd && subs<maxBandwidth; q++){
			// Using your wicked XOR logic, but masked to safe positive int
			// If q!=r, xor is >0. min(1, >0) is 1. 
			// If q==r, xor is 0. min(1, 0) is 0.
			subs+=Math.min(1, (query[q]^ref[q+offset])&0xFF);
		}
		return Math.min(subs+1, maxBandwidth);
	}

//	/**
//	 * Counts substitutions up to maxBandwidth.
//	 * Returns min(subs+1, maxBandwidth);
//	 * Updates pos[0] to optimal offset if SIMD is used.
//	 */
//	public static int decideBandwidth(final byte[] query, final byte[] ref, 
//			final int[] pos, final int maxBandwidth){
//		assert(pos!=null && pos.length>=2);
//		if(Shared.SIMD){
//			return simd.SIMDAlignByte.countSubs(query, ref, pos, maxBandwidth);
//		}
//		final int qLen=query.length, rLen=ref.length;
//		int offset=pos[0], subs=0;
//		if(offset==Integer.MIN_VALUE){offset=0; pos[0]=0;}
//		
//		for(int qpos=0, rpos=offset; qpos<qLen && subs<maxBandwidth; qpos++, rpos++){
//			final byte q=query[qpos], r=(rpos>=0 && rpos<rLen ? ref[rpos] : (byte)'$');
////			subs+=(q==r ? 0 : 1);
//			subs+=Math.min(1, (q^r)&0xFF);
//		}
//		return Math.min(subs+1, maxBandwidth);
//	}
}
