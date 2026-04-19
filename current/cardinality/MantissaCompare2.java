package cardinality;

import rand.FastRandomXoshiro;
import parse.Parse;
import parse.PreParser;
import shared.Tools;

/**
 * Gemini's rewrite of MantissaCompare. Single-bucket 16-bit simulation.
 * State 0 = Established/Small Gap. State Cap = Outlier/Big Gap.
 */
public class MantissaCompare2 {

    static final int MODE_MANTISSA=0, MODE_ANDTISSA=1, MODE_NLZ2=2, MODE_HISTORY=3, MODE_LUCK=4, MODE_HISTMANT=5, MODE_TWINTAIL=6, MODE_MASTERSLAVE=7, MODE_SUPER1BIT=8, MODE_AVGNLZ=9, MODE_MEDIAN3=10, MODE_MINPLUS1=11, MODE_CTLL=12;
    static final String[] MODE_NAMES={"Mantissa", "Andtissa", "NLZ2", "History", "Luck", "HistMant", "TwinTail", "MasterSlave", "Super1Bit", "AvgNlz", "Median3", "MinPlus1", "CTLL"};
    static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

    static final int SAMPLE_ALL=0, SAMPLE_ENTRY=1, SAMPLE_BOTH=2;
    static final String[] SAMPLE_NAMES={"all", "entry", "both"};

    public static void main(String[] args) throws InterruptedException {
        {PreParser pp=new PreParser(args, null, false); args=pp.args;}
        int inner=32768, outer=131072, maxTier=11, mode=MODE_MANTISSA, bits=2, threads=1;
        int sampleMode=SAMPLE_ALL;
        int hbits=-1, mbits=-1; // for combined mode; -1 = use 'bits'

        for(String arg : args){
            final String[] ab=arg.split("=");
            final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
            if(a.equals("inner")){inner=Parse.parseIntKMG(b);}
            else if(a.equals("outer")){outer=Parse.parseIntKMG(b);}
            else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
            else if(a.equals("bits")){bits=Integer.parseInt(b);}
            else if(a.equals("hbits")){hbits=Integer.parseInt(b);}
            else if(a.equals("mbits")){mbits=Integer.parseInt(b);}
            else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
            else if(a.equals("sample") || a.equals("s")){
                if(b.equals("all")){sampleMode=SAMPLE_ALL;}
                else if(b.equals("entry")){sampleMode=SAMPLE_ENTRY;}
                else if(b.equals("both")){sampleMode=SAMPLE_BOTH;}
                else{throw new RuntimeException("Unknown sample mode '"+b+"'");}
            }
            else if(a.equals("mode")){
                if(b.equals("mantissa")){mode=MODE_MANTISSA;}
                else if(b.equals("andtissa")){mode=MODE_ANDTISSA;}
                else if(b.equals("nlz2")){mode=MODE_NLZ2;}
                else if(b.equals("history")){mode=MODE_HISTORY;}
                else if(b.equals("luck")){mode=MODE_LUCK;}
                else if(b.equals("histmant") || b.equals("historymantissa")){mode=MODE_HISTMANT;}
                else if(b.equals("twintail") || b.equals("ttll")){mode=MODE_TWINTAIL;}
                else if(b.equals("masterslave") || b.equals("ms")){mode=MODE_MASTERSLAVE;}
                else if(b.equals("super1bit") || b.equals("s1b") || b.equals("dualhist")){mode=MODE_SUPER1BIT;}
                else if(b.equals("avgnlz") || b.equals("avg")){mode=MODE_AVGNLZ;}
                else if(b.equals("median3") || b.equals("med3")){mode=MODE_MEDIAN3;}
                else if(b.equals("minplus1") || b.equals("mp1")){mode=MODE_MINPLUS1;}
                else if(b.equals("ctll") || b.equals("compressedtier") || b.equals("manttier")){mode=MODE_CTLL;}
                else{throw new RuntimeException("Unknown mode '"+b+"'");}
            }
            else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
        }
        if(mode==MODE_HISTMANT){
            if(hbits<0){hbits=2;}
            if(mbits<0){mbits=2;}
        }

        // TwinTail and MasterSlave always have 16 states: (h1:2)(h0:2) = 4 bits
        final int numStates=(mode==MODE_TWINTAIL || mode==MODE_MASTERSLAVE) ? 16 :
            (mode==MODE_HISTMANT) ? (1<<(hbits+mbits)) : (1<<bits);
        final int cap=numStates-1;
        final String modeDesc=(mode==MODE_TWINTAIL || mode==MODE_MASTERSLAVE) ?
            MODE_NAMES[mode]+" (h1:2)(h0:2)=16 states" :
            (mode==MODE_HISTMANT) ?
            MODE_NAMES[mode]+" hbits="+hbits+" mbits="+mbits :
            MODE_NAMES[mode]+" bits="+bits;
        System.err.println("MantissaCompare2: "+modeDesc
            +" states="+numStates+" inner="+inner+" outer="+outer
            +" maxtier="+maxTier+" threads="+threads+" sample="+SAMPLE_NAMES[sampleMode]);
        final long t0=System.nanoTime();

        final long[][] count=new long[64][numStates];
        final double[][] sum=new double[64][numStates];
        final double[][] geoSum=new double[64][numStates]; // sum of log(card)
        final double[][] harmSum=new double[64][numStates]; // sum of 1/card

        // Capture effectively final copies for lambda
        final int fInner=inner, fOuter=outer, fMaxTier=maxTier, fMode=mode;
        final int fBits=bits, fHbits=hbits, fMbits=mbits, fNumStates=numStates, fCap=cap;
        final int fSampleMode=sampleMode;

        final Thread[] workers=new Thread[threads];
        final long[][][] tCount=new long[threads][64][numStates];
        final double[][][] tSum=new double[threads][64][numStates];
        final double[][][] tGeo=new double[threads][64][numStates];
        final double[][][] tHarm=new double[threads][64][numStates];

        for(int ti=0; ti<threads; ti++){
            final int threadIdx=ti;
            final int trialStart=ti*(outer/threads);
            final int trialEnd=(ti==threads-1) ? outer : (ti+1)*(outer/threads);
            workers[ti]=new Thread(()->{
                final long[][] lc=tCount[threadIdx];
                final double[][] ls=tSum[threadIdx];
                final double[][] lg=tGeo[threadIdx];
                final double[][] lh=tHarm[threadIdx];
                for(int trial=trialStart; trial<trialEnd; trial++){
                    final FastRandomXoshiro rng=new FastRandomXoshiro(trial);
                    char stored=0;
                    int maxNlz=-1, extra=0, luckSecond=-1;
                    int prevTier=-1, prevState=-1;
                    int entryCard=0; // card when entering current state

                    for(int card=1; card<=fInner; card++){
                final long number=rng.nextLong();
                final long key=Tools.hash64shift(number);
                int nlz=Long.numberOfLeadingZeros(key);

                // Multi-hash NLZ override for dual/triple hash modes
                if(fMode==MODE_SUPER1BIT){
                    // Unsigned max of 2 hashes -> min NLZ. P(NLZ>=k) = 4^(-k).
                    final long key2=Tools.hash64shift(number^123456789L);
                    final long dualKey=(Long.compareUnsigned(key, key2)>=0 ? key : key2);
                    nlz=Long.numberOfLeadingZeros(dualKey);
                }else if(fMode==MODE_AVGNLZ){
                    // Average of 2 NLZs. Softer compression than min.
                    final long key2=Tools.hash64shift(number^123456789L);
                    final int nlz2=Long.numberOfLeadingZeros(key2);
                    nlz=(nlz+nlz2)/2;
                }else if(fMode==MODE_MEDIAN3){
                    // Median of 3 hash codes (unsigned). Suppresses both extremes.
                    final long key2=Tools.hash64shift(number^123456789L);
                    final long key3=Tools.hash64shift(number^987654321L);
                    final long median;
                    if(Long.compareUnsigned(key, key2)<=0){
                        if(Long.compareUnsigned(key2, key3)<=0){median=key2;}
                        else if(Long.compareUnsigned(key, key3)<=0){median=key3;}
                        else{median=key;}
                    }else{
                        if(Long.compareUnsigned(key, key3)<=0){median=key;}
                        else if(Long.compareUnsigned(key2, key3)<=0){median=key3;}
                        else{median=key2;}
                    }
                    nlz=Long.numberOfLeadingZeros(median);
                }else if(fMode==MODE_MINPLUS1){
                    // min(NLZ1, NLZ2+1): asymmetric dual hash.
                    final long key2=Tools.hash64shift(number^123456789L);
                    final int nlz2=Long.numberOfLeadingZeros(key2)+1;
                    nlz=Math.min(nlz, nlz2);
                }else if(fMode==MODE_CTLL){
                    // CompressedTierLogLog: mantissa threshold maps rawNlz to tier.
                    // halfNlz = 2*rawNlz + mantissa, tier = halfNlz/3.
                    // Tier ratio = (sqrt(2))^3 = 2*sqrt(2) ≈ 2.828.
                    final int rawNlz=nlz;
                    final int mBits;
                    if(rawNlz>=43){
                        mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
                    }else{
                        mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
                    }
                    final int mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
                    nlz=(2*rawNlz+mant)/3; // tier number
                }

                if(fMode==MODE_HISTMANT){
                    // Combined history+mantissa: mantissa in low bits, history in high bits
                    final int mcap=(1<<fMbits)-1, hcap=(1<<fHbits)-1;
                    final int mshift=63-nlz-fMbits;
                    final int mval=(mshift<0) ? 0 : (int)((~(key>>>mshift))&mcap);
                    final int newNlzStored=Math.min(nlz+1, 63);
                    final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;
                    final int oldMant=(stored>0) ? (stored&mcap) : 0;
                    final int oldHist=(stored>0) ? ((stored>>>fMbits)&hcap) : 0;
                    if(newNlzStored>oldNlzStored){
                        // New max: update mantissa, shift history register
                        final int k=newNlzStored-oldNlzStored;
                        final int carry=(stored>0) ? (1<<fHbits) : 0;
                        final int newHist=((oldHist|carry)>>k)&hcap;
                        extra=mval|(newHist<<fMbits);
                        maxNlz=newNlzStored-1;
                    }else if(newNlzStored==oldNlzStored && newNlzStored>0){
                        // Same max: update mantissa if better, keep history
                        if(mval>oldMant){ extra=mval|(oldHist<<fMbits); }
                    }else if(newNlzStored<oldNlzStored){
                        // Sub-max: update history bits, keep mantissa
                        final int diff=oldNlzStored-newNlzStored;
                        if(diff>=1 && diff<=fHbits){
                            final int newHist=oldHist|(1<<(fHbits-diff));
                            extra=oldMant|(newHist<<fMbits);
                        }
                    }
                }else if(fMode==MODE_MANTISSA || fMode==MODE_ANDTISSA || fMode==MODE_NLZ2){
                    int val=0;
                    if(fMode==MODE_MANTISSA){
                        final int shift=63-nlz-fBits;
                        val=(shift<0) ? 0 : (int)((~(key>>>shift))&fCap);
                    }else if(fMode==MODE_ANDTISSA){
                        final long anded=key&(key<<2);
                        final int shift=63-nlz-fBits;
                        val=(shift<0) ? 0 : (int)((anded>>>shift)&fCap);
                    }else{
                        final int consumed=nlz+1;
                        val=(consumed>=64) ? fCap : Math.min(fCap, Long.numberOfLeadingZeros(key<<consumed));
                    }
                    if(nlz>maxNlz){ maxNlz=nlz; extra=val; }
                    else if(nlz==maxNlz && val>extra){ extra=val; }

                }else if(fMode==MODE_HISTORY || fMode==MODE_SUPER1BIT || fMode==MODE_AVGNLZ || fMode==MODE_MEDIAN3 || fMode==MODE_MINPLUS1 || fMode==MODE_CTLL){
                    // All history-based modes use the same update logic;
                    // the NLZ was already overridden above for multi-hash modes.
                    final int newNlzStored=Math.min(nlz+1, 63);
                    final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;
                    if(newNlzStored>oldNlzStored){
                        final int k=newNlzStored-oldNlzStored;
                        final int oldHist=(stored>0) ? (stored&fCap) : 0;
                        final int carry=(stored>0) ? (1<<fBits) : 0;
                        extra=((oldHist|carry)>>k)&fCap;
                        maxNlz=newNlzStored-1;
                    }else if(newNlzStored<oldNlzStored){
                        final int diff=oldNlzStored-newNlzStored;
                        if(diff>=1 && diff<=fBits){ extra|=(1<<(fBits-diff)); }
                    }
                }else if(fMode==MODE_LUCK){
                    if(nlz>maxNlz){ luckSecond=maxNlz; maxNlz=nlz; }
                    else if(nlz>luckSecond && nlz<maxNlz){ luckSecond=nlz; }
                    extra=(maxNlz>=0 && luckSecond>=0) ? Math.min(fCap, Math.max(0, maxNlz-luckSecond-1)) : fCap;
                }else if(fMode==MODE_TWINTAIL){
                    // TwinTail: two independent 2-bit shift registers (h0, h1).
                    // histBit selects which tail to update. In real TTLL this is
                    // (key>>>bucketBits)&1; here we use bit 10 of the key.
                    final int histBit=(int)((key>>>10)&1);
                    final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;
                    final int newNlzStored=Math.min(nlz+1, 63);
                    final int delta=newNlzStored-oldNlzStored;

                    if(stored==0 && nlz>=0){
                        // First element: set maxNlz, set MSB of triggered tail, other bits 0
                        maxNlz=nlz;
                        extra=(histBit==0) ? 0x2 : 0x8; // MSB of h0=bit1, MSB of h1=bit3
                    }else if(delta>0){
                        // New max NLZ: shift both tails right, set MSB of triggered tail
                        maxNlz=newNlzStored-1;
                        final int oldH=(stored&0xF);
                        final int shiftAmt=Math.min(delta, 2);
                        final int h0=(oldH&0x3)>>>shiftAmt;
                        final int h1=((oldH>>>2)&0x3)>>>shiftAmt;
                        extra=(h1<<2)|h0;
                        // Set MSB of triggered tail
                        final int bitPos=(histBit==0) ? 1 : 3;
                        extra|=(1<<bitPos);
                    }else if(delta==0 && stored>0){
                        // Equal to current max: set MSB of triggered tail
                        final int bitPos=(histBit==0) ? 1 : 3;
                        extra=(stored&0xF)|(1<<bitPos);
                    }else if(delta==-1){
                        // One below current max: set LSB of triggered tail
                        final int bitPos=(histBit==0) ? 0 : 2;
                        extra=(stored&0xF)|(1<<bitPos);
                    }else{
                        // delta < -1: ignore (too far below ceiling)
                        extra=(stored&0xF);
                    }
                    // Keep maxNlz as the max seen so far
                    if(stored>0 && newNlzStored-1<=maxNlz){ /* maxNlz unchanged */ }
                }else if(fMode==MODE_MASTERSLAVE){
                    // Master/Slave: h0=master (pure NLZ history like UDLL6),
                    // h1=slave (bit-filtered, only responds to histBit==1).
                    // h0 MSB = "saw NLZ-1", h0 LSB = "saw NLZ-2" (unconditional).
                    // h1 MSB = "saw current NLZ with bit=1", h1 LSB = "saw NLZ-1 with bit=1".
                    final int histBit=(int)((key>>>10)&1);
                    final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;
                    final int newNlzStored=Math.min(nlz+1, 63);
                    final int delta=newNlzStored-oldNlzStored;
                    final int oldH0=(stored>0) ? (stored&0x3) : 0;
                    final int oldH1=(stored>0) ? ((stored>>>2)&0x3) : 0;
                    int h0, h1;

                    if(stored==0 && nlz>=0){
                        // First element: h0 gets carry (via UDLL6 formula), h1 gets MSB if histBit==1
                        maxNlz=nlz;
                        h0=0; // no old history, carry shifts out for delta=nlz+1
                        // Actually: carry = (1<<2) = 4, but stored==0 means no carry
                        // First element is always at NLZ=nlz, so h0=00 (no sub-tiers seen)
                        h1=(histBit==1) ? 0x2 : 0; // MSB of h1 if bit=1
                    }else if(delta>0){
                        // Advance: h0 uses UDLL6 carry mechanism, h1 just shifts
                        maxNlz=newNlzStored-1;
                        final int carry=(stored>0) ? (1<<2) : 0; // bit above 2-bit field
                        h0=((oldH0|carry)>>>delta)&0x3;
                        h1=(oldH1>>>Math.min(delta, 2))&0x3;
                        // h1: set MSB if advancing hash has histBit==1
                        if(histBit==1){ h1|=0x2; }
                    }else if(delta==0 && stored>0){
                        // Same NLZ: h0 unchanged, h1 sets MSB if histBit==1
                        h0=oldH0;
                        h1=oldH1;
                        if(histBit==1){ h1|=0x2; }
                    }else if(delta==-1){
                        // One below: h0 sets MSB, h1 sets LSB if histBit==1
                        h0=oldH0|0x2; // set MSB = "saw NLZ-1"
                        h1=oldH1;
                        if(histBit==1){ h1|=0x1; } // set LSB = "saw NLZ-1 with bit=1"
                    }else if(delta==-2){
                        // Two below: h0 sets LSB, h1 nothing
                        h0=oldH0|0x1; // set LSB = "saw NLZ-2"
                        h1=oldH1;
                    }else{
                        // delta < -2: too far below, no change
                        h0=oldH0;
                        h1=oldH1;
                    }
                    extra=(h1<<2)|h0;
                    if(stored>0 && newNlzStored-1<=maxNlz){ /* maxNlz unchanged */ }
                }

                if(maxNlz>=0){
                    stored=(char)(((maxNlz+1)<<10)|extra);
                    if(maxNlz<=fMaxTier && extra<fNumStates){
                        final boolean stateChanged=(maxNlz!=prevTier || extra!=prevState);
                        if(fSampleMode==SAMPLE_ALL){
                            lc[maxNlz][extra]++;
                            ls[maxNlz][extra]+=card;
                            lg[maxNlz][extra]+=Math.log(card);
                            lh[maxNlz][extra]+=1.0/card;
                        }else if(fSampleMode==SAMPLE_ENTRY){
                            if(stateChanged){
                                lc[maxNlz][extra]++;
                                ls[maxNlz][extra]+=card;
                                lg[maxNlz][extra]+=Math.log(card);
                                lh[maxNlz][extra]+=1.0/card;
                            }
                        }else{ // SAMPLE_BOTH
                            if(stateChanged){
                                // Record exit from previous state
                                if(prevTier>=0 && prevTier<=fMaxTier && prevState<fNumStates){
                                    lc[prevTier][prevState]++;
                                    ls[prevTier][prevState]+=card;
                                    lg[prevTier][prevState]+=Math.log(card);
                                    lh[prevTier][prevState]+=1.0/card;
                                }
                                // Record entry to new state
                                lc[maxNlz][extra]++;
                                ls[maxNlz][extra]+=card;
                                lg[maxNlz][extra]+=Math.log(card);
                                lh[maxNlz][extra]+=1.0/card;
                            }
                        }
                        if(stateChanged){
                            prevTier=maxNlz; prevState=extra; entryCard=card;
                        }
                    }
                }
            }
        }
            });
            workers[ti].start();
        }
        for(Thread w : workers){w.join();}
        // Merge per-thread arrays
        for(int ti=0; ti<threads; ti++){
            for(int t=0; t<64; t++){
                for(int s=0; s<numStates; s++){
                    count[t][s]+=tCount[ti][t][s];
                    sum[t][s]+=tSum[ti][t][s];
                    geoSum[t][s]+=tGeo[ti][t][s];
                    harmSum[t][s]+=tHarm[ti][t][s];
                }
            }
        }

        final double elapsed=(System.nanoTime()-t0)*1e-9;
        System.err.println(String.format("Elapsed: %.1fs", elapsed));

        // Report per-tier
        StringBuilder header=new StringBuilder("Tier\tTotal\tAvgCard\tRatio");
        for(int s=0; s<numStates; s++){header.append("\tP("+s+")\tLinCF("+s+")\tGeoCF("+s+")\tHarmCF("+s+")");}
        System.out.println(header);

        double prevAvg=0;
        for(int t=0; t<=maxTier; t++){
            long tierTotal=0; double tierSum=0, tierGeo=0, tierHarm=0;
            for(int s=0; s<numStates; s++){
                tierTotal+=count[t][s]; tierSum+=sum[t][s];
                tierGeo+=geoSum[t][s]; tierHarm+=harmSum[t][s];
            }
            if(tierTotal<100) continue;
            double tierLinAvg=tierSum/tierTotal;
            double tierGeoAvg=Math.exp(tierGeo/tierTotal);
            double tierHarmAvg=tierTotal/tierHarm;
            StringBuilder sb=new StringBuilder();
            sb.append(t).append('\t').append(tierTotal);
            sb.append('\t').append(String.format("%.2f", tierLinAvg));
            sb.append('\t').append(prevAvg>0 ? String.format("%.4f", tierLinAvg/prevAvg) : "-");
            for(int s=0; s<numStates; s++){
                sb.append('\t').append(String.format("%.6f", count[t][s]/(double)tierTotal));
                if(count[t][s]>10){
                    double stateLinAvg=sum[t][s]/count[t][s];
                    double stateGeoAvg=Math.exp(geoSum[t][s]/count[t][s]);
                    double stateHarmAvg=count[t][s]/harmSum[t][s];
                    sb.append('\t').append(String.format("%+.8f", Math.log(stateLinAvg/tierLinAvg)/Math.log(2)));
                    sb.append('\t').append(String.format("%+.8f", Math.log(stateGeoAvg/tierGeoAvg)/Math.log(2)));
                    sb.append('\t').append(String.format("%+.8f", Math.log(stateHarmAvg/tierHarmAvg)/Math.log(2)));
                }else{sb.append("\tN/A\tN/A\tN/A");}
            }
            System.out.println(sb);
            prevAvg=tierLinAvg;
        }

        // Weighted summary
        long gTotal=0; double[] gCF=new double[numStates];
        for(int t=3; t<=maxTier; t++){
            long tt=0; double ts=0;
            for(int s=0; s<numStates; s++){ tt+=count[t][s]; ts+=sum[t][s]; }
            if(tt<100) continue;
            double ta=ts/tt;
            for(int s=0; s<numStates; s++){
                if(count[t][s]>10){ gCF[s]+=Math.log(sum[t][s]/count[t][s]/ta)/Math.log(2)*count[t][s]; }
                gTotal+=count[t][s];
            }
        }
        final String summaryLabel=(mode==MODE_TWINTAIL || mode==MODE_MASTERSLAVE) ? MODE_NAMES[mode]+" 4-bit" :
            (mode==MODE_HISTMANT) ? MODE_NAMES[mode]+" "+hbits+"h+"+mbits+"m" :
            MODE_NAMES[mode]+" "+bits+"-bit";
        System.out.println("\n# "+summaryLabel+" weighted CFs (tiers 3-"+maxTier+"):");
        StringBuilder cfLine=new StringBuilder("# CF=[");
        for(int s=0; s<numStates; s++){
            if(s>0) cfLine.append(", ");
            cfLine.append(String.format("%+.8f", gCF[s]/(gTotal>0?gTotal:1)));
        }
        cfLine.append("]");
        System.out.println(cfLine);
    }
}
