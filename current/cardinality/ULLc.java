package cardinality;

import shared.Tools;

/**
 * ULLc: ErtlULLb + minZeros tracking and early exit.
 * Registers remain absolute 8-bit (identical to ErtlULLb).
 * Early exit skips elements whose NLZ is too low to affect any register.
 * <p>
 * With HISTORY_MARGIN=2, elements with NLZ < minZeros-2 cannot affect any
 * register (the 2-bit history only tracks 2 levels below the max). So early
 * exit rejects those elements via eeMask, producing IDENTICAL register states
 * to ErtlULLb — just faster.
 */
public final class ULLc extends CardinalityTracker {

    final byte[] registers;
    private int minZeros=0;
    /** Number of buckets with NLZ <= minZeros. When 0, advance floor. */
    private int floorCount;
    /** Early exit mask. Rejects elements with NLZ < minZeros - HISTORY_MARGIN. */
    private long eeMask=-1L;
    private int filledBuckets=0;

    ULLc(){this(2048, 31, -1, 0);}
    ULLc(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        registers=new byte[buckets];
        floorCount=buckets; // all empty = all at or below floor
    }
    @Override public ULLc copy(){return new ULLc(buckets, k, -1, minProb);}

    /** Extract NLZ from an absolute 8-bit Ertl register. Returns -1 for empty. */
    private int nlzFromRegister(int reg){
        if(reg==0){return -1;}
        return ((reg&0xFF)>>>2)-(bucketBits-1);
    }

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);

        // Early exit: reject elements whose NLZ is too low to affect any register.
        if(Long.compareUnsigned(key, eeMask)>0){return;}
        branch1++;

        final int idx=(int)(key&bucketMask);
        final int nlz=Long.numberOfLeadingZeros(key);
        final int p=bucketBits;

        final byte oldState=registers[idx];
        long hashPrefix=ErtlULL.unpack(oldState);
        hashPrefix|=1L<<(nlz+p-1);
        final byte newState=ErtlULL.pack(hashPrefix);
        if((newState&0xFF)<=(oldState&0xFF)){return;} // no change
        branch2++;

        registers[idx]=newState;
        lastCardinality=-1;

        // Track floor transitions
        final int oldNlz=nlzFromRegister(oldState);
        final int newNlz=nlzFromRegister(newState);

        if(oldState==0){filledBuckets++;}

        // If old bucket was at or below floor and new is above floor: decrement floorCount
        if(oldNlz<=minZeros && newNlz>minZeros){
            floorCount--;
            if(floorCount<=0){
                advanceFloor();
            }
        }
    }

    /**
     * Scans all registers to find the new minimum NLZ, advances minZeros,
     * updates eeMask, and recounts buckets at the new floor.
     */
    private void advanceFloor(){
        // Find minimum NLZ across all registers
        int minNlz=Integer.MAX_VALUE;
        for(int i=0; i<buckets; i++){
            int r=registers[i]&0xFF;
            if(r==0){
                minNlz=-1;
                break; // can't advance past empty buckets
            }
            int nlz=((r&0xFF)>>>2)-(bucketBits-1);
            if(nlz<minNlz){minNlz=nlz;}
        }
        if(minNlz<=minZeros){
            // Can't advance — some bucket is still at or below current floor
            floorCount=0;
            for(int i=0; i<buckets; i++){
                if(nlzFromRegister(registers[i]&0xFF)<=minZeros){floorCount++;}
            }
            return;
        }
        minZeros=minNlz;
        // Update eeMask: reject NLZ < minZeros - HISTORY_MARGIN
        int exitThreshold=Math.max(0, minZeros-HISTORY_MARGIN);
        eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;

        // Count buckets at the new floor
        floorCount=0;
        for(int i=0; i<buckets; i++){
            if(nlzFromRegister(registers[i]&0xFF)<=minZeros){floorCount++;}
        }
    }

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

    public int filledBuckets(){return filledBuckets;}
    public double occupancy(){return (double)filledBuckets/buckets;}
    public int getMinZeros(){return minZeros;}

    @Override
    public double[] rawEstimates(){
        final int total=11+4+CardinalityStats.NUM_DLC_TIERS;
        final double[] r=new double[total];
        final double fgra=fgraEstimate();
        r[0]=fgra; r[1]=fgra; r[4]=fgra; r[6]=fgra; r[8]=fgra;
        return r;
    }

    private static final int HISTORY_MARGIN=2;

    public long branch1=0, branch2=0;
    public double branch1Rate(){return branch1/(double)Math.max(1, added);}
    public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}
}
