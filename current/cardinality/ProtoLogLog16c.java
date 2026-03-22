package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * PLL16c: 16-bit register prototype using ULL8's proven correction pipeline.
 * Multiplicative per-tier corrections, single summarize loop, hybridDLL estimator.
 * Supports all 5 sub-NLZ modes: history, luck, mantissa, andtissa, nlz2.
 */
public final class ProtoLogLog16c extends CardinalityTracker {

    public static final int MODE_NONE=0, MODE_MANTISSA=1, MODE_NLZ2=2, MODE_HISTORY=4, MODE_LUCK=8, MODE_ANDTISSA=16;
    public static int MODE=MODE_MANTISSA;
    public static int MANTISSA_BITS=2, NLZ2_BITS=2, HISTORY_BITS=2, LUCK_BITS=2, ANDTISSA_BITS=2, NLZ_BITS=6;

    // Per-state CFs from MantissaCompare2 tier 8 (131072 trials)
    // These are additive NLZ corrections: multiplier = 2^(-CF)
    public static final double[] CF_MANTISSA_2={-0.4928, -0.2633, -0.0342, +0.2640};
    public static final double[] CF_ANDTISSA_2={-0.3295, +0.0017, +0.2241, +0.3588};
    public static final double[] CF_NLZ2_2={-0.3562, -0.0342, +0.1648, +0.3382};
    public static final double[] CF_LUCK_1={+0.1686, -1.3682};
    public static final double[] CF_LUCK_2={+0.1686, -1.1662, -2.3051, -3.6058};
    public static final double[] CF_LUCK_3={+0.1686, -1.1662, -2.3051, -3.3914, -4.4281, -5.5035, -6.5513, -7.9929};
    public static final double[] CF_HISTORY_1={-1.3682, +0.1686};
    public static final double[] CF_HISTORY_2={-2.5208, -1.1662, -1.6570, +0.2078};
    public static final double[] CF_HISTORY_3={-3.6058, -2.3051, -2.7423, -1.1164, -2.8410, -1.5820, -2.1327, +0.2137};

    // Per-tier CFs for history 2-bit (tiers 0-2)
    public static final double[][] CF_HISTORY_2_TIERS={
        {+0.0000,  0.0000,  0.0000,  0.0000},
        {-1.8472,  0.0000, +0.2051,  0.0000},
        {-3.1800, -1.3051, -1.9303, +0.2266},
    };
    // Combined history+mantissa CFs from simulator (tier 8, 124M samples)
    // State index = mantissa_val | (history_val << mbits)
    public static final double[] CF_HISTMANT_2H2M={
        -2.5995, -2.5518, -2.4953, -2.4397,  // h=0: m=0,1,2,3
        -1.3324, -1.2396, -1.1559, -1.0402,  // h=1
        -1.7472, -1.6938, -1.6204, -1.5905,  // h=2
        -0.2122, -0.0389, +0.1694, +0.4223}; // h=3
    public static final double[] CF_HISTMANT_1H3M={
        -1.5812, -1.5239, -1.4818, -1.4209, -1.3730, -1.3330, -1.2667, -1.1835,  // h=0: m=0..7
        -0.3085, -0.2349, -0.1372, -0.0452, +0.0691, +0.1802, +0.3047, +0.4605}; // h=1: m=0..7
    public static final double[] CF_HISTMANT_3H1M={
        -3.6133, -3.5716,  // h=0: m=0,1
        -2.3542, -2.2627,  // h=1
        -2.7558, -2.7068,  // h=2
        -1.2263, -1.0449,  // h=3
        -2.8923, -2.8413,  // h=4
        -1.6402, -1.5328,  // h=5
        -2.1591, -2.1219,  // h=6
        -0.0999, +0.3347}; // h=7

    // Per-tier CFs for history 3-bit (tiers 0-3, from simulator 500k trials)
    public static final double[][] CF_HISTORY_3_TIERS={
        {0, 0, 0, 0, 0, 0, 0, 0},
        {-1.8456, 0, 0, 0, +0.2051, 0, 0, 0},
        {-3.1813, 0, -1.2962, 0, -1.9285, 0, +0.2280, 0},
        {-4.3479, -2.5276, -3.0717, -1.2338, -3.2677, -1.6811, -2.3671, +0.2221},
    };
    // Per-tier CFs for combined 1h+3m (tiers 0-1, from simulator 500k trials)
    // State = mantissa_val(3bit) | (history_val(1bit) << 3)
    public static final double[][] CF_HISTMANT_1H3M_TIERS={
        {-0.8227, -0.6814, -0.5253, -0.3557, -0.1542, +0.0391, +0.2449, +0.4656,
               0,       0,       0,       0,       0,       0,       0,       0},
        {-2.1487, -2.0697, -2.0043, -1.9300, -1.8541, -1.7696, -1.6953, -1.5896,
         -0.3752, -0.2654, -0.1551, -0.0362, +0.0784, +0.2230, +0.3644, +0.5282},
    };

    // Per-tier CFs for luck 2-bit (tiers 0-2)
    public static final double[][] CF_LUCK_2_TIERS={
        { 0.0000,  0.0000,  0.0000, +0.0000},
        {+0.2051,  0.0000,  0.0000, -1.8472},
        {+0.1855, -1.3051,  0.0000, -3.1800},
    };

    // Active correction state — set by setMode()
    public static boolean USE_EEMASK=true;
    public static boolean SKIP_UNCHANGED_LUCK=true;
    public static double[] CORRECTION_TABLE=null;
    public static double[][] TIER_TABLES=null;
    /** Additive offset applied to all per-state CFs: tierMult = 2^(-(cf+offset)). */
    public static double CF_OFFSET=0;

    final short[] maxArray;
    private final byte[] luckSecond;
    private int filledBuckets=0;
    private long eeMask=-1L;
    private int minNlzStored=0; // minimum nlzStored (absNlz+1) across filled buckets, 0 before all filled
    private int minCount=0; // buckets at minNlzStored; 0 means not yet tracking
    int[] lastRawNlz;

    ProtoLogLog16c(){this(2048, 31, -1, 0);}
    ProtoLogLog16c(Parser p){ super(p); maxArray=new short[buckets]; luckSecond=usesLuck()?new byte[buckets]:null; }
    ProtoLogLog16c(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        maxArray=new short[buckets];
        luckSecond=usesLuck()?new byte[buckets]:null;
    }
    @Override public ProtoLogLog16c copy(){return new ProtoLogLog16c(buckets, k, -1, minProb);}

    private static boolean usesLuck(){return (MODE & MODE_LUCK)!=0;}
    private static boolean usesHistory(){return (MODE & MODE_HISTORY)!=0;}
    private static int nlzShift(){return 16-NLZ_BITS;}
    private static int maxNlzStored(){return (1<<NLZ_BITS)-1;}
    private static int getAbsNlz(int stored){return (stored>>>nlzShift())-1;}
    private static int getExtra(int stored){return stored & ((1<<nlzShift())-1);}

    public static void setMode(int mode){
        MODE=mode;
        if(mode==MODE_HISTORY){
            CORRECTION_TABLE=(HISTORY_BITS==3?CF_HISTORY_3:HISTORY_BITS==1?CF_HISTORY_1:CF_HISTORY_2);
            TIER_TABLES=(HISTORY_BITS==3?CF_HISTORY_3_TIERS:HISTORY_BITS==2?CF_HISTORY_2_TIERS:null);
        }else if(mode==MODE_LUCK){
            CORRECTION_TABLE=(LUCK_BITS==3?CF_LUCK_3:LUCK_BITS==1?CF_LUCK_1:CF_LUCK_2);
            TIER_TABLES=(LUCK_BITS==2?CF_LUCK_2_TIERS:null);
        }else if(mode==MODE_MANTISSA){
            CORRECTION_TABLE=CF_MANTISSA_2; TIER_TABLES=null;
        }else if(mode==MODE_ANDTISSA){
            CORRECTION_TABLE=CF_ANDTISSA_2; TIER_TABLES=null;
        }else if(mode==MODE_NLZ2){
            CORRECTION_TABLE=CF_NLZ2_2; TIER_TABLES=null;
        }else if(mode==(MODE_HISTORY|MODE_MANTISSA)){
            // Combined: use simulator-derived 16-state CFs when available
            if(HISTORY_BITS==2 && MANTISSA_BITS==2){CORRECTION_TABLE=CF_HISTMANT_2H2M;}
            else if(HISTORY_BITS==1 && MANTISSA_BITS==3){CORRECTION_TABLE=CF_HISTMANT_1H3M;}
            else if(HISTORY_BITS==3 && MANTISSA_BITS==1){CORRECTION_TABLE=CF_HISTMANT_3H1M;}
            else{CORRECTION_TABLE=null;}
            // Per-tier tables: use simulator-derived when available, else outer product
            if(HISTORY_BITS==1 && MANTISSA_BITS==3){TIER_TABLES=CF_HISTMANT_1H3M_TIERS;}
            else{
                // Outer product of history tier CFs × mantissa CFs
                final double[][] htiers=(HISTORY_BITS==2?CF_HISTORY_2_TIERS:HISTORY_BITS==3?CF_HISTORY_3_TIERS:null);
                if(htiers!=null){
                    final int mstates=1<<MANTISSA_BITS, hstates=1<<HISTORY_BITS;
                    TIER_TABLES=new double[htiers.length][];
                    for(int t=0; t<htiers.length; t++){
                        TIER_TABLES[t]=new double[mstates*hstates];
                        for(int h=0; h<hstates; h++){
                            for(int m=0; m<mstates; m++){
                                TIER_TABLES[t][m+(h<<MANTISSA_BITS)]=CF_MANTISSA_2[Math.min(m,CF_MANTISSA_2.length-1)]+htiers[t][h];
                            }
                        }
                    }
                }else{TIER_TABLES=null;}
            }
        }else{CORRECTION_TABLE=null; TIER_TABLES=null;}
    }

    // ---- computeExtra and hashAndStore: identical to PLL16b ----

    private int computeExtra(long key, int nlz, int bucket, int oldStored){
        int extra=0, shift=0;
        final int oldNlzS=(oldStored>0) ? (oldStored>>>nlzShift()) : 0;
        final int newNlzS=Math.min(nlz+1, maxNlzStored());

        if((MODE & MODE_MANTISSA)!=0){
            final int mbits=MANTISSA_BITS, mshift=63-nlz-mbits;
            int val=(mshift<0)?0:(int)((~(key>>>mshift))&((1<<mbits)-1));
            extra|=(val<<shift); shift+=mbits;
        }
        if((MODE & MODE_ANDTISSA)!=0){
            final int abits=ANDTISSA_BITS, ashift=63-nlz-abits;
            long anded=key&(key<<2);
            int val=(ashift<0)?0:(int)((anded>>>ashift)&((1<<abits)-1));
            extra|=(val<<shift); shift+=abits;
        }
        if((MODE & MODE_NLZ2)!=0){
            final int nbits=NLZ2_BITS, cap2=(1<<nbits)-1, consumed=nlz+1;
            int val=(consumed>=64)?cap2:Math.min(cap2, Long.numberOfLeadingZeros(key<<consumed));
            extra|=(val<<shift); shift+=nbits;
        }
        if((MODE & MODE_HISTORY)!=0){
            final int hbits=HISTORY_BITS;
            int hist=0;
            if(oldStored>0){
                if(newNlzS>oldNlzS){
                    int oldH=(oldStored>>shift)&((1<<hbits)-1);
                    hist=((oldH|(1<<hbits))>>(newNlzS-oldNlzS))&((1<<hbits)-1);
                }else{
                    hist=(oldStored>>shift)&((1<<hbits)-1);
                }
            }
            extra|=(hist<<shift); shift+=hbits;
        }
        if((MODE & MODE_LUCK)!=0){
            final int lbits=LUCK_BITS, cap=(1<<lbits)-1;
            final int oldM=(oldStored>0)?getAbsNlz(oldStored):-1;
            final int effM=Math.max(nlz, oldM);
            final int sec=(luckSecond[bucket]&0xFF)-1;
            int val=(effM>=0 && sec>=0) ? Math.min(cap, Math.max(0, effM-sec-1)) : cap;
            extra|=(val<<shift); shift+=lbits;
        }
        return extra;
    }

    /** Margin below global min NLZ for early exit.
     *  History needs max-1/max-2; luck needs to track second-best. */
    private static int eeMargin(){
        if(usesHistory()) return HISTORY_BITS;
        if(usesLuck()) return LUCK_BITS;
        return 0;
    }

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);

        // Early exit: key has too few leading zeros to affect any bucket
        if(USE_EEMASK && Long.compareUnsigned(key, eeMask)>0){return;}

        final int nlz=Long.numberOfLeadingZeros(key);
        final int bucket=(int)(key&bucketMask);

        final long micro=(key>>bucketBits)&0x3FL;
        microIndex|=(1L<<micro);
        if(USE_MICRO && Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}

        final int newNlzS=Math.min(nlz+1, maxNlzStored());
        final int oldStored=maxArray[bucket]&0xFFFF;
        final int oldNlzS=(oldStored>0) ? (oldStored>>>nlzShift()) : 0;

        // Update luck leaderboard; track whether it changed
        boolean luckChanged=false;
        if(usesLuck()){
            final int curM=(oldStored>0)?getAbsNlz(oldStored):-1;
            final int sec=(luckSecond[bucket]&0xFF)-1;
            if(nlz>curM){
                luckSecond[bucket]=(byte)(curM+1);
                luckChanged=true;
            }else if(nlz>sec && nlz<curM){
                luckSecond[bucket]=(byte)(nlz+1);
                luckChanged=true;
            }
        }

        if(newNlzS>oldNlzS){
            lastCardinality=-1;
            final boolean wasEmpty=(oldStored==0);
            if(wasEmpty){filledBuckets++;}
            maxArray[bucket]=(short)((newNlzS<<nlzShift())|computeExtra(key, nlz, bucket, oldStored));
            if(wasEmpty && filledBuckets>=buckets){
                // All buckets now filled — initialize min tracking
                minNlzStored=0; minCount=0;
                advanceEeMask();
            }else if(minCount>0 && oldNlzS==minNlzStored){
                if(--minCount==0){advanceEeMask();}
            }
        }else if(newNlzS==oldNlzS && newNlzS>0){
            final int nEx=computeExtra(key, nlz, bucket, oldStored);
            final int oEx=getExtra(oldStored);
            final boolean improved=usesLuck() ? (nEx<oEx) : (nEx>oEx);
            if(improved){
                maxArray[bucket]=(short)((newNlzS<<nlzShift())|nEx);
                lastCardinality=-1;
            }
        }else if(oldNlzS>newNlzS){
            if(usesHistory()){
                final int diff=oldNlzS-newNlzS;
                final int hshift=getHistoryShift();
                if(diff>=1 && diff<=HISTORY_BITS){
                    final int bit=1<<(HISTORY_BITS-diff);
                    final int newStored=oldStored|(bit<<hshift);
                    if(newStored!=oldStored){ maxArray[bucket]=(short)newStored; lastCardinality=-1; }
                }
            }
            if(usesLuck() && (SKIP_UNCHANGED_LUCK ? luckChanged : true)){
                // Only recompute gap when leaderboard actually changed
                final int nEx=computeExtra(key, nlz, bucket, oldStored);
                final int oEx=getExtra(oldStored);
                if(nEx<oEx){
                    maxArray[bucket]=(short)((oldNlzS<<nlzShift())|nEx);
                    lastCardinality=-1;
                }
            }
        }
    }

    /** Advance eeMask when all buckets have moved past the current minimum. */
    private void advanceEeMask(){
        // Count buckets at new candidate minimum
        while(minCount==0 && minNlzStored<63){
            minNlzStored++;
            // Count how many buckets are at the new minNlzStored
            int c=0;
            for(int i=0; i<buckets; i++){
                final int nlzS=(maxArray[i]&0xFFFF)>>>nlzShift();
                if(nlzS==minNlzStored){c++;}
            }
            minCount=c;
        }
        // Update eeMask: skip elements with NLZ < (minNlzStored-1) - margin
        // minNlzStored is stored as absNlz+1, so absNlz = minNlzStored-1
        final int absMin=minNlzStored-1;
        final int effective=Math.max(0, absMin-eeMargin());
        eeMask=(effective>=63) ? 0L : (-1L>>>(effective));
    }

    // ---- Summarize: ULL8's approach — multiplicative corrections, single loop ----

    private CardinalityStats summarize(){
        final int states=1<<getActiveExtraBits();
        final int mask=(1<<nlzShift())-1;
        // minTier: below this, use per-tier tables (if available) or no correction.
        // With tier tables: tiers 0..len-1 use tier CFs, tier len+ use steady-state.
        // Without tier tables: history needs HISTORY_BITS+1, luck needs 2^LUCK_BITS.
        // For combined modes, minTier is driven by the history/luck component alone
        // (mantissa reaches steady state immediately — it's just hash bits).
        final int minTier;
        if(TIER_TABLES!=null){minTier=TIER_TABLES.length;}
        else if(usesLuck()){minTier=(1<<LUCK_BITS);}
        else if(usesHistory()){minTier=HISTORY_BITS+1;}
        else{minTier=0;}

        // Build per-tier multiplier table (like ULL8)
        final double[][] tierMult=new double[64][states];
        for(int t=0; t<64; t++){
            double[] cf=null;
            if(TIER_TABLES!=null && t<TIER_TABLES.length){
                cf=TIER_TABLES[t]; // per-tier correction
            }else if(t>=minTier){
                cf=CORRECTION_TABLE; // steady-state correction
            }
            // else: cf stays null → multiplier = 1.0 (no correction for low tiers without tier tables)
            for(int s=0; s<states; s++){
                tierMult[t][s]=(cf!=null && s<cf.length) ? Math.pow(2.0, -(cf[s]+CF_OFFSET)) : 1.0;
            }
        }

        double difSum=0, hllSumFilled=0, gSum=0;
        int count=0;
        final int[] nlzCounts=new int[64];

        for(int i=0; i<buckets; i++){
            final int stored=maxArray[i]&0xFFFF;
            if(stored>0){
                final int absNlz=getAbsNlz(stored);
                final int extra=stored&mask;
                if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
                final double m=tierMult[Math.min(absNlz, 63)][Math.min(extra, states-1)];
                final double base=Math.pow(2.0, -absNlz)*m;
                final double dif=(absNlz==0 ? (double)Long.MAX_VALUE :
                    (absNlz<64 ? (double)(1L<<(63-absNlz)) : 1.0))*m;
                difSum+=dif;
                hllSumFilled+=base;
                gSum+=Math.log(Tools.max(1, dif));
                count++;
            }
        }
        lastRawNlz=nlzCounts;
        return new CardinalityStats(difSum, hllSumFilled, hllSumFilled,
            gSum, count, buckets, null, CF_MATRIX, CF_BUCKETS,
            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex,
            nlzCounts, 0);
    }

    @Override
    public final long cardinality(){
        if(lastCardinality>=0) return lastCardinality;
        final CardinalityStats s=summarize();
        long card=(long)s.hybridDLL();
        card=Math.max(card, s.microCardinality());
        card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
        lastCardinality=card;
        return card;
    }

    @Override public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
    public int filledBuckets(){return filledBuckets;}
    public double occupancy(){return (double)filledBuckets/buckets;}
    @Override public final float[] compensationFactorLogBucketsArray(){return null;}

    @Override
    public double[] rawEstimates(){
        final CardinalityStats s=summarize();
        return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
    }

    private static int getActiveExtraBits(){
        int b=0;
        if((MODE&MODE_MANTISSA)!=0) b+=MANTISSA_BITS;
        if((MODE&MODE_ANDTISSA)!=0) b+=ANDTISSA_BITS;
        if((MODE&MODE_NLZ2)!=0) b+=NLZ2_BITS;
        if((MODE&MODE_HISTORY)!=0) b+=HISTORY_BITS;
        if((MODE&MODE_LUCK)!=0) b+=LUCK_BITS;
        return b;
    }

    private static int getHistoryShift(){
        int s=0;
        if((MODE&MODE_MANTISSA)!=0) s+=MANTISSA_BITS;
        if((MODE&MODE_ANDTISSA)!=0) s+=ANDTISSA_BITS;
        if((MODE&MODE_NLZ2)!=0) s+=NLZ2_BITS;
        return s;
    }

    private static int CF_BUCKETS=2048;
    private static float[][] CF_MATRIX=null;
    public static void setCFMatrix(float[][] matrix, int buckets){
        CF_MATRIX=matrix; CF_BUCKETS=buckets;
    }
}
