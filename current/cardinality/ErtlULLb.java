package cardinality;

import shared.Tools;

/**
 * Copy of ErtlULL using BBTools bit conventions:
 * - Bucket index from LOW bits of hash (key & bucketMask)
 * - NLZ from leading zeros of full 64-bit hash
 *
 * Everything else is identical to ErtlULL. Used to validate that
 * the bit source doesn't affect accuracy before building UDLL6.
 */
public final class ErtlULLb extends CardinalityTracker {

    final byte[] registers;

    ErtlULLb(){this(2048, 31, -1, 0);}
    ErtlULLb(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        registers=new byte[buckets];
    }
    @Override public ErtlULLb copy(){return new ErtlULLb(buckets, k, -1, minProb);}

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);

        // BBTools convention: low bits = bucket, NLZ of full hash
        final int idx=(int)(key&bucketMask);
        final int nlz=Long.numberOfLeadingZeros(key);
        final int p=bucketBits; // log2(buckets)

        final byte oldState=registers[idx];
        long hashPrefix=ErtlULL.unpack(oldState);
        hashPrefix|=1L<<(nlz+p-1); // same bit position formula as Ertl
        final byte newState=ErtlULL.pack(hashPrefix);
        if((newState&0xFF)>(oldState&0xFF)){
            registers[idx]=newState;
            lastCardinality=-1;
        }
    }

    /** Raw FGRA estimate using the static method from ErtlULL. */
    public double fgraEstimate(){
        return ErtlULL.fgraEstimateStatic(registers, bucketBits);
    }

    @Override
    public final long cardinality(){
        if(lastCardinality>=0) return lastCardinality;
        long card=Math.max(0, Math.round(fgraEstimate()));
        card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
        lastCardinality=card;
        return card;
    }

    @Override public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
    @Override public final float[] compensationFactorLogBucketsArray(){return null;}

    public int filledBuckets(){
        int c=0; for(byte b:registers){if(b!=0)c++;} return c;
    }
    public double occupancy(){return (double)filledBuckets()/buckets;}

    @Override
    public double[] rawEstimates(){
        final int total=11+6+CardinalityStats.NUM_DLC_TIERS;
        final double[] r=new double[total];
        final double fgra=fgraEstimate();
        r[0]=fgra; r[1]=fgra; r[4]=fgra; r[6]=fgra; r[8]=fgra;
        return r;
    }
}
