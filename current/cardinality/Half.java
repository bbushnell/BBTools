package cardinality;

/**
 * IEEE 754 half-precision (binary16) conversion, Java-8-safe (no
 * Float.floatToFloat16, which is Java 20+).  Round-to-nearest-even.
 * Full correctness including subnormals and overflow-to-infinity;
 * NaN input is rejected by assertion (crash loud, never wrong).
 *
 * Purpose: compact per-sketch chronicle state (Fll53TrajC) — checkpoint
 * scalars stored as 16-bit floats so the memory saved can buy organisms.
 *
 * @author Amber (priority: Brian)
 * @date July 2026
 */
public final class Half {

	/** float -> binary16 bits, round-to-nearest-even. */
	public static short f2h(float f){
		assert(!Float.isNaN(f)) : "NaN reaching fp16 storage";
		final int b=Float.floatToRawIntBits(f);
		final int s=(b>>>16)&0x8000;
		final int e=(b>>>23)&0xFF;
		final int m=b&0x7FFFFF;
		if(e>=143){return (short)(s|0x7C00);}   // overflow (incl. Inf) -> Inf
		if(e>=113){   // normal half range
			int h=s|((e-112)<<10)|(m>>>13);
			final int rem=m&0x1FFF;
			if(rem>0x1000 || (rem==0x1000 && (h&1)==1)){h++;}   // RNE; carry ok
			return (short)h;
		}
		if(e>=101){   // subnormal half range
			final int shift=126-e;   // 14..25
			final int mant=m|0x800000;
			int h=mant>>>shift;
			final int rem=mant&((1<<shift)-1);
			final int mid=1<<(shift-1);
			if(rem>mid || (rem==mid && (h&1)==1)){h++;}
			return (short)(s|h);
		}
		return (short)s;   // underflow -> signed zero
	}

	/** binary16 bits -> float (exact). */
	public static float h2f(short h){
		final int s=(h&0x8000)<<16;
		final int em=h&0x7FFF;
		if(em==0){return Float.intBitsToFloat(s);}
		final int e=em>>>10, m=em&0x3FF;
		if(e==0x1F){return Float.intBitsToFloat(s|0x7F800000|(m<<13));}
		if(e==0){   // half subnormal -> float normal
			int k=0, mm=m;
			while((mm&0x400)==0){mm<<=1; k++;}
			return Float.intBitsToFloat(s|((113-k)<<23)|((mm&0x3FF)<<13));
		}
		return Float.intBitsToFloat(s|((e+112)<<23)|(m<<13));
	}

	/** Round-trips f through fp16: what a stored value becomes. */
	public static float round(float f){return h2f(f2h(f));}

	/** Exhaustive self-test: h2f/f2h idempotent on all non-NaN bit patterns,
	 * and round(x) is within half an fp16 ulp on random normal-range floats. */
	public static void main(String[] args){
		int tested=0;
		for(int i=0; i<65536; i++){
			final short h=(short)i;
			if(((h&0x7FFF)>>>10)==0x1F && (h&0x3FF)!=0){continue;}   // NaN patterns
			final short back=f2h(h2f(h));
			if(back!=h && !(h==(short)0x8000 && back==0)){   // -0 folds to +0? keep signed
				if(back!=h){throw new AssertionError(
					String.format("idempotence: %04x -> %g -> %04x", i, h2f(h), back&0xFFFF));}
			}
			tested++;
		}
		final java.util.Random r=new java.util.Random(42);
		for(int i=0; i<1000000; i++){
			final float f=(float)((r.nextDouble()*2-1)*Math.pow(2, r.nextInt(30)-15));
			final float g=round(f);
			final float ulp=Math.max(Math.ulp(g)*8192, 6.0e-8f);   // fp16 ulp ~ 2^13 float ulps
			if(Math.abs(f-g)>ulp*0.5f*1.0001f && Math.abs(f)<65504){
				throw new AssertionError(f+" -> "+g+" err "+(f-g)+" > half-ulp "+ulp*0.5);
			}
			tested++;
		}
		System.out.println("Half: "+tested+" cases OK");
	}
}
