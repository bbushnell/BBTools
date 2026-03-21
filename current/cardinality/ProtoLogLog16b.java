package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Gemini's rewrite of ProtoLogLog16. Clean implementation.
 */
public final class ProtoLogLog16b extends CardinalityTracker {

    public static final int MODE_NONE=0, MODE_MANTISSA=1, MODE_NLZ2=2, MODE_HISTORY=4, MODE_LUCK=8, MODE_ANDTISSA=16;
    public static int MODE=MODE_MANTISSA;
    public static int MANTISSA_BITS=2, NLZ2_BITS=2, HISTORY_BITS=2, LUCK_BITS=2, ANDTISSA_BITS=2, NLZ_BITS=6;

    public static double[] CORRECTION_TABLE=null;
    public static double[][] TIER_TABLES=null; // per-tier CFs for tiers 0..TIER_TABLES.length-1
    public static double CF_OFFSET=0.0;

    // Per-state CFs from MantissaCompare2 (131072 trials)
    // Mantissa/Andtissa/NLZ2: same CFs for all tiers, no minimum tier needed
    public static final double[] CF_MANTISSA_2={-0.4928, -0.2633, -0.0342, +0.2640};
    public static final double[] CF_ANDTISSA_2={-0.3295, +0.0017, +0.2241, +0.3588};
    public static final double[] CF_NLZ2_2={-0.3562, -0.0342, +0.1648, +0.3382};

    // History/Luck: steady-state CFs (tier 8) plus per-tier tables for tiers 0-2
    public static final double[] CF_LUCK_1={+0.1686, -1.3682};
    public static final double[] CF_LUCK_2={+0.1686, -1.1662, -2.3051, -3.6058};
    public static final double[] CF_LUCK_3={+0.1686, -1.1662, -2.3051, -3.3914, -4.4281, -5.5035, -6.5513, -7.9929};
    public static final double[] CF_HISTORY_1={-1.3682, +0.1686};
    public static final double[] CF_HISTORY_2={-2.5208, -1.1662, -1.6570, +0.2078};
    public static final double[] CF_HISTORY_3={-3.6058, -2.3051, -2.7423, -1.1164, -2.8410, -1.5820, -2.1327, +0.2137};

    // Per-tier CFs for history 2-bit (tiers 0-2, N/A states use 0.0)
    public static final double[][] CF_HISTORY_2_TIERS={
        {+0.0000,  0.0000,  0.0000,  0.0000},  // tier 0: all state 0
        {-1.8472,  0.0000, +0.2051,  0.0000},  // tier 1: states 0,2 only
        {-3.1800, -1.3051, -1.9303, +0.2266},  // tier 2: all states
    };
    // Per-tier CFs for luck 2-bit (tiers 0-2, N/A states use 0.0)
    public static final double[][] CF_LUCK_2_TIERS={
        { 0.0000,  0.0000,  0.0000, +0.0000},  // tier 0: all state cap
        {+0.2051,  0.0000,  0.0000, -1.8472},  // tier 1: states 0,3 only
        {+0.1855, -1.3051,  0.0000, -3.1800},  // tier 2: states 0,1,3
    };

    final short[] maxArray;
    private final byte[] luckSecond;
    private int filledBuckets=0;
    int[] lastRawNlz;

    ProtoLogLog16b(){this(2048, 31, -1, 0);}
    ProtoLogLog16b(Parser p){ super(p); maxArray=new short[buckets]; luckSecond=usesLuck()?new byte[buckets]:null; }
    ProtoLogLog16b(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        maxArray=new short[buckets];
        luckSecond=usesLuck()?new byte[buckets]:null;
    }
    @Override public ProtoLogLog16b copy(){return new ProtoLogLog16b(buckets, k, -1, minProb);}

    private static boolean usesLuck(){return (MODE & MODE_LUCK)!=0;}
    private static boolean usesHistory(){return (MODE & MODE_HISTORY)!=0;}
    private static int nlzShift(){return 16-NLZ_BITS;}
    private static int maxNlzStored(){return (1<<NLZ_BITS)-1;}
    private static int getAbsNlz(int stored){return (stored>>>nlzShift())-1;}
    private static int getExtra(int stored){return stored & ((1<<nlzShift())-1);}

    public static void setMode(int mode){
        MODE=mode;
        CF_OFFSET=0.0;
        TIER_TABLES=null; // per-tier tables need more calibration work; disabled for now
        if(mode==MODE_HISTORY){
            CORRECTION_TABLE=(HISTORY_BITS==3?CF_HISTORY_3:HISTORY_BITS==1?CF_HISTORY_1:CF_HISTORY_2);
        }else if(mode==MODE_LUCK){
            CORRECTION_TABLE=(LUCK_BITS==3?CF_LUCK_3:LUCK_BITS==1?CF_LUCK_1:CF_LUCK_2);
        }else if(mode==MODE_MANTISSA){
            CORRECTION_TABLE=CF_MANTISSA_2;
        }else if(mode==MODE_ANDTISSA){
            CORRECTION_TABLE=CF_ANDTISSA_2;
        }else if(mode==MODE_NLZ2){
            CORRECTION_TABLE=CF_NLZ2_2;
        }else{CORRECTION_TABLE=null;}
    }

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
                    hist=(oldStored>>shift)&((1<<hbits)-1); // preserve on tie
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

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);
        final int nlz=Long.numberOfLeadingZeros(key);
        final int bucket=(int)(key&bucketMask);

        final long micro=(key>>bucketBits)&0x3FL;
        microIndex|=(1L<<micro);
        if(USE_MICRO && Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}

        final int newNlzS=Math.min(nlz+1, maxNlzStored());
        final int oldStored=maxArray[bucket]&0xFFFF;
        final int oldNlzS=(oldStored>0) ? (oldStored>>>nlzShift()) : 0;

        // Update luck leaderboard before anything else
        if(usesLuck()){
            final int curM=(oldStored>0)?getAbsNlz(oldStored):-1;
            final int sec=(luckSecond[bucket]&0xFF)-1;
            if(nlz>curM){
                luckSecond[bucket]=(byte)(curM+1); // old max becomes second
            }else if(nlz>sec && nlz<curM){
                luckSecond[bucket]=(byte)(nlz+1); // silver medalist (strictly between)
            }
        }

        if(newNlzS>oldNlzS){
            // Promotion
            lastCardinality=-1;
            if(oldStored==0){filledBuckets++;}
            maxArray[bucket]=(short)((newNlzS<<nlzShift())|computeExtra(key, nlz, bucket, oldStored));
        }else if(newNlzS==oldNlzS && newNlzS>0){
            // Same NLZ: update extra if improved
            final int nEx=computeExtra(key, nlz, bucket, oldStored);
            final int oEx=getExtra(oldStored);
            // Luck: smaller gap is better (min-update). Others: bigger is better (max-update).
            final boolean improved=usesLuck() ? (nEx<oEx) : (nEx>oEx);
            if(improved){
                maxArray[bucket]=(short)((newNlzS<<nlzShift())|nEx);
                lastCardinality=-1;
            }
        }else if(oldNlzS>newNlzS){
            // Sub-max element: update extra bits if improved
            if(usesHistory()){
                final int diff=oldNlzS-newNlzS;
                final int hshift=getHistoryShift();
                if(diff>=1 && diff<=HISTORY_BITS){
                    final int bit=1<<(HISTORY_BITS-diff);
                    final int newStored=oldStored|(bit<<hshift);
                    if(newStored!=oldStored){ maxArray[bucket]=(short)newStored; lastCardinality=-1; }
                }
            }
            if(usesLuck()){
                // Recompute gap from updated leaderboard
                final int nEx=computeExtra(key, nlz, bucket, oldStored);
                final int oEx=getExtra(oldStored);
                if(nEx<oEx){
                    maxArray[bucket]=(short)((oldNlzS<<nlzShift())|nEx);
                    lastCardinality=-1;
                }
            }
        }
    }

    /**
     * Restore the approximate original value from leading zeros + fractional offset.
     * restore(nlz) = 2^(63-nlz). With float nlz, gives fractional precision.
     */
    private static double restore(double nlz){
        return Math.pow(2.0, 63.0-nlz);
    }

    /**
     * Summarize using Brian's Mean formula with per-bucket CF correction.
     * CF is added to the NLZ before restore: restore(NLZ + CF[state]).
     * Negative CF (outlier) → lower effective NLZ → larger val → smaller estimate.
     */
    private CardinalityStats summarize(){
        final int states=1<<getActiveExtraBits();
        final int mask=(1<<nlzShift())-1;
        final int extraBits=getActiveExtraBits();
        // minTier: history needs bits+1, luck needs 2^bits, others can correct all tiers
        final int minTier;
        if(TIER_TABLES!=null){minTier=TIER_TABLES.length;}
        else if(usesLuck()){minTier=(1<<extraBits);}
        else if(usesHistory()){minTier=extraBits+1;}
        else{minTier=0;}

        double sum=0;
        double cardSum=0, logCardSum=0;
        int count=0;
        final int[] nlzCounts=new int[64];

        for(int i=0; i<buckets; i++){
            final int stored=maxArray[i]&0xFFFF;
            if(stored>0){
                final int absNlz=getAbsNlz(stored);
                final int extra=stored&mask;
                if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}

                // Compute effective NLZ = absNlz + CF[state]
                double effNlz=absNlz;
                if(extra<states){
                    if(absNlz>=minTier && CORRECTION_TABLE!=null && CORRECTION_TABLE.length==states){
                        effNlz+=CORRECTION_TABLE[extra]+CF_OFFSET;
                    }else if(TIER_TABLES!=null && absNlz<TIER_TABLES.length
                        && TIER_TABLES[absNlz]!=null && TIER_TABLES[absNlz].length==states){
                        effNlz+=TIER_TABLES[absNlz][extra]+CF_OFFSET;
                    }
                }

                double val=restore(effNlz);
                double dif=Long.MAX_VALUE-val;
                if(val>0){
                    sum+=dif;
                    // Per-bucket cardinality estimate
                    double perCard=2.0*(Long.MAX_VALUE/val);
                    cardSum+=perCard;
                    logCardSum+=Math.log(perCard);
                    count++;
                }
            }
        }

        final int subsets=count;
        final double mantissaFactor=0.7213428177;
        final double emptyBucketModifier=((count+buckets)/(float)(buckets+buckets));

        // Original: average restored values, then convert
        final double meanDif=sum/Tools.max(subsets, 1);
        final double meanVal=(double)Long.MAX_VALUE-meanDif;
        final double estimatePerSet=2.0*(Long.MAX_VALUE/meanVal);
        final double total=estimatePerSet*subsets*mantissaFactor*emptyBucketModifier;
        lastCorrectedMean=total;

        // Arithmetic mean of per-bucket cardinality estimates
        final double arithCard=cardSum/Tools.max(subsets, 1);
        lastArithMean=arithCard*subsets*mantissaFactor*emptyBucketModifier;

        // Geometric mean of per-bucket cardinality estimates
        final double geoCard=Math.exp(logCardSum/Tools.max(subsets, 1));
        lastGeoMean=geoCard*subsets*mantissaFactor*emptyBucketModifier;

        // Build corrected sums for CardinalityStats.
        // Use the corrected 'sum' (from effNlz) as difSum so DThHyb/Mean benefit
        // from per-state corrections. Also build corrected hllSum and gSum.
        double hllSum=0, gSum=0;
        for(int i=0; i<buckets; i++){
            final int stored=maxArray[i]&0xFFFF;
            if(stored>0){
                final int absNlz=getAbsNlz(stored);
                final int extra=stored&mask;
                double effNlz=absNlz;
                if(extra<states){
                    if(absNlz>=minTier && CORRECTION_TABLE!=null && CORRECTION_TABLE.length==states){
                        effNlz+=CORRECTION_TABLE[extra]+CF_OFFSET;
                    }else if(TIER_TABLES!=null && absNlz<TIER_TABLES.length
                        && TIER_TABLES[absNlz]!=null && TIER_TABLES[absNlz].length==states){
                        effNlz+=TIER_TABLES[absNlz][extra]+CF_OFFSET;
                    }
                }
                hllSum+=Math.pow(2.0, -effNlz);
                double dif=restore(effNlz);
                gSum+=Math.log(Tools.max(1, dif));
            }
        }
        lastRawNlz=nlzCounts;
        // Pass corrected sums so Mean/DThHyb/HMean all use sub-NLZ corrections
        return new CardinalityStats(sum, hllSum, hllSum,
            gSum, count, buckets, null, CF_MATRIX, CF_BUCKETS,
            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex,
            nlzCounts, 0);
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

    /** Cached estimates from last summarize(). */
    private double lastCorrectedMean=0;
    private double lastArithMean=0;
    private double lastGeoMean=0;
    public static int ESTIMATE_MODE=0; // 0=original, 1=arithMean, 2=geoMean

    @Override public final long cardinality(){
        if(lastCardinality>=0) return lastCardinality;
        summarize();
        long card=(long)(ESTIMATE_MODE==2 ? lastGeoMean : ESTIMATE_MODE==1 ? lastArithMean : lastCorrectedMean);
        card=Math.max(card, 0);
        card=Math.min(clampToAdded?added:Long.MAX_VALUE, card);
        lastCardinality=card;
        return card;
    }

    @Override public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
    public int filledBuckets(){return filledBuckets;}
    public double occupancy(){return (double)filledBuckets/buckets;}
    @Override public final float[] compensationFactorLogBucketsArray(){return null;}
    @Override public double[] rawEstimates(){
        final CardinalityStats s=summarize();
        final double est=(ESTIMATE_MODE==2 ? lastGeoMean : ESTIMATE_MODE==1 ? lastArithMean : lastCorrectedMean);
        return s.toArray(est);
    }

    private static int CF_BUCKETS=2048;
    private static float[][] CF_MATRIX=null; // must be set via setCFMatrix or cffile
    public static void setCFMatrix(float[][] matrix, int buckets){
        CF_MATRIX=matrix; CF_BUCKETS=buckets;
    }
}
