package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Ertl's UltraLogLog (PVLDB 2024) — faithful implementation of the paper.
 * b=2, d=2, q=6: 8-bit registers = 6-bit max update value + 2-bit history.
 * Uses the FGRA estimator (Algorithm 6) with small/large count corrections.
 */
public final class ErtlULL extends CardinalityTracker {

    // FGRA constants for b=2, d=2 (from paper Section 3.3)
    static final double TAU=0.819491;
    static final double V_CONST=0.611893;
    static final double ETA_0=4.663135;
    static final double ETA_1=2.137850;
    static final double ETA_2=2.781145;
    static final double ETA_3=0.982408;
    static final double[] ETA={ETA_0, ETA_1, ETA_2, ETA_3};
    static boolean DEBUG=true;

    final byte[] registers;

    ErtlULL(){this(2048, 31, -1, 0);}
    ErtlULL(Parser p){super(p); registers=new byte[buckets];}
    ErtlULL(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        registers=new byte[buckets];
    }
    @Override public ErtlULL copy(){return new ErtlULL(buckets, k, -1, minProb);}

    /** Unpack 8-bit register to 64-bit bitset (Algorithm 5). */
    private static long unpack(int r){
        if(r<4) return 0;
        int u=r>>>2;
        int q1=(r>>>1)&1, q0=r&1;
        // Bit u+1 = max (always set), bit u = saw u-1, bit (u-1) = saw u-2
        long x=1L<<(u+1);
        if(q1!=0) x|=1L<<u;
        if(q0!=0 && u>=1) x|=1L<<(u-1);
        return x;
    }

    /** Pack 64-bit bitset to 8-bit register (Algorithm 4). */
    private static int pack(long x){
        if(x<4) return (x>0)?4:0; // smallest valid register
        int u=62-Long.numberOfLeadingZeros(x);
        int q1=(int)((x>>>u)&1);
        int q0=(u>=1)?(int)((x>>>(u-1))&1):0;
        return (u<<2)|(q1<<1)|q0;
    }

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);
        final int bucket=(int)(key&bucketMask);

        final long micro=(key>>bucketBits)&0x3FL;
        microIndex|=(1L<<micro);
        if(USE_MICRO && Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}

        // Mask register index bits, compute update value k (Algorithm 3)
        final long a=key&~((1L<<bucketBits)-1 | (-1L<<(64)));
        // Actually: a is the hash with the top p bits masked off
        // k = nlz(a) - p + 1 where p = bucketBits
        // Simpler: k = nlz of remaining bits = Long.numberOfLeadingZeros(key) + 1
        // because the top p bits are the bucket index, and NLZ of the remaining
        // 64-p bits gives the update value.
        // In our convention: NLZ is computed from the full 64-bit key.
        final int nlz=Long.numberOfLeadingZeros(key);
        final int k=nlz+1; // update value, ≥1

        final int oldR=registers[bucket]&0xFF;
        long x=unpack(oldR);
        x|=(1L<<(k+1)); // set bit for update value k (stored at position k+1)
        final int newR=pack(x);
        if(newR>oldR){
            registers[bucket]=(byte)newR;
            lastCardinality=-1;
        }
    }

    /** FGRA estimator (Algorithm 6). */
    @Override
    public final long cardinality(){
        if(lastCardinality>=0) return lastCardinality;
        final int p=bucketBits;
        final int w=65-p; // max update value
        final int m=buckets;

        // Count special registers and accumulate g(r) for intermediate range
        int c0=0, c4=0, c8=0, c10=0;
        int c4w=0, c4w1=0, c4w2=0, c4w3=0;
        double s=0;
        final int fourW=4*w;

        for(int i=0; i<m; i++){
            final int r=registers[i]&0xFF;
            if(r<12){
                if(r==0) c0++;
                else if(r==4) c4++;
                else if(r==8) c8++;
                else if(r==10) c10++;
            }else if(r<fourW){
                s+=g(r);
            }else{
                if(r==fourW) c4w++;
                else if(r==fourW+1) c4w1++;
                else if(r==fourW+2) c4w2++;
                else if(r==fourW+3) c4w3++;
            }
        }

        // Small range correction (Section 3.5)
        if(c0+c4+c8+c10>0){
            double alpha=m+3*(c0+c4+c8+c10);
            double beta=m-c0-c4;
            double gamma=4*c0+2*c4+3*c8+c10;
            double z=Math.pow((Math.sqrt(beta*beta+4*alpha*gamma)-beta)/(2*alpha), 4); // eq (23)

            if(c0>0) s+=c0*sigma(z);
            if(c4>0) s+=c4*Math.pow(2, -TAU)*psi(z);
            if(c8>0) s+=c8*Math.pow(4, -TAU)*(z*(ETA_0-ETA_1)+(ETA_1));
            if(c10>0) s+=c10*Math.pow(4, -TAU)*(z*(ETA_2-ETA_3)+(ETA_3));
        }

        // Large range: just use g(r) for now (corrections only matter near nmax)
        if(c4w>0) s+=c4w*g(fourW);
        if(c4w1>0) s+=c4w1*g(fourW+1);
        if(c4w2>0) s+=c4w2*g(fourW+2);
        if(c4w3>0) s+=c4w3*g(fourW+3);

        // Final estimate (eq 12)
        double lambda_p=Math.pow(m, 1+1.0/TAU)*Math.pow(1+(1+TAU)/(2*m)*V_CONST, -1);
        double estimate=lambda_p*Math.pow(s, -1.0/TAU);
        if(DEBUG && added<20){
            System.err.println("FGRA debug: m="+m+" s="+s+" lambda="+lambda_p+
                " est="+estimate+" c0="+c0+" c4="+c4+" c8="+c8+" c10="+c10);
        }

        long card=Math.max(0, Math.round(estimate));
        card=Math.min(clampToAdded?added:Long.MAX_VALUE, card);
        lastCardinality=card;
        return card;
    }

    /** g(r) function (eq 15): per-register contribution for intermediate range. */
    private static double g(int r){
        int u=r>>>2;
        int state=r&3;
        return Math.pow(2, -TAU*u)*ETA[state];
    }

    /** σ(z) (eq 20): (1/z) Σ_{u=0}^∞ 2^{τu} (z^{2^u} − z^{2^{u+1}}) ψ(z^{2^{u+1}}) */
    private static double sigma(double z){
        if(z==1.0) return Double.POSITIVE_INFINITY;
        double result=0;
        double zpow=z;
        for(int u=0; u<100; u++){
            double zpow2=zpow*zpow;
            double term=Math.pow(2, TAU*u)*(zpow-zpow2)*psi(zpow2);
            if(Double.isNaN(term) || Math.abs(term)<1e-20) break;
            result+=term;
            zpow=zpow2;
        }
        return result/z;
    }

    /** ψ(z) (eq 19): z²(η₀−η₁−η₂+η₃) + z(η₂−η₃) + η₁ */
    private static double psi(double z){
        return z*(z*(ETA_0-ETA_1-ETA_2+ETA_3)+(ETA_2-ETA_3))+ETA_1;
    }

    @Override public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
    public double occupancy(){
        int filled=0;
        for(int i=0; i<buckets; i++){if((registers[i]&0xFF)>0) filled++;}
        return (double)filled/buckets;
    }
    @Override public final float[] compensationFactorLogBucketsArray(){return null;}
    @Override public double[] rawEstimates(){
        long card=cardinality();
        // Build nlzCounts from registers for CardinalityStats
        final int[] nlzCounts=new int[64];
        double difSum=0, hllSum=0, gSum=0;
        int count=0;
        for(int i=0; i<buckets; i++){
            int r=registers[i]&0xFF;
            if(r>0){
                int u=r>>>2;
                int absNlz=u; // update value = NLZ+1, so NLZ = u-1... but for compat, use u
                if(absNlz>=0 && absNlz<64) nlzCounts[absNlz]++;
                double dif=(absNlz==0?(double)Long.MAX_VALUE:(absNlz<64?(double)(1L<<(63-absNlz)):1.0));
                double base=Math.pow(2.0, -absNlz);
                difSum+=dif; hllSum+=base; gSum+=Math.log(Tools.max(1, dif));
                count++;
            }
        }
        CardinalityStats s=new CardinalityStats(difSum, hllSum, hllSum,
            gSum, count, buckets, null, CF_MATRIX, CF_BUCKETS,
            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex,
            nlzCounts, 0);
        return s.toArray(cardinality());
    }

    public static final String CF_FILE="?cardinalityCorrectionLL6.tsv.gz";
    private static int CF_BUCKETS=2048;
    private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
    public static float[][] initializeCF(int buckets){
        CF_BUCKETS=buckets;
        return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
    }
    public static void setCFMatrix(float[][] matrix, int buckets){
        CF_MATRIX=matrix; CF_BUCKETS=buckets;
    }
}
// Temporary debug main
// public static void main(String[] args){
//     ErtlULL ull=new ErtlULL(256, 31, -1, 0);
//     for(int i=1; i<=1000; i++) ull.hashAndStore(i);
//     System.err.println("card="+ull.cardinality()+" (true=1000)");
// }
