package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Ertl's UltraLogLog (PVLDB 2024) — ported from hash4j production Java code.
 * 8-bit registers: 6-bit NLZ + 2-bit sub-NLZ history.
 * Uses the optimal FGRA estimator (closed-form, no CF table needed).
 *
 * Register format uses the mapping r -> r + 4*p - 8 from the paper,
 * so the maximum register value is always 255 regardless of p.
 *
 * Ported from: https://github.com/dynatrace-oss/hash4j
 * Original author: Otmar Ertl, Dynatrace (Apache 2.0 license)
 */
public final class ErtlULL extends CardinalityTracker {

    final byte[] registers;

    ErtlULL(){this(2048, 31, -1, 0);}
    ErtlULL(Parser p){super(p); registers=new byte[buckets];}
    ErtlULL(int buckets_, int k_, long seed, float minProb_){
        super(buckets_, k_, seed, minProb_);
        registers=new byte[buckets];
    }
    @Override public ErtlULL copy(){return new ErtlULL(buckets, k, -1, minProb);}

    /*--------------------------------------------------------------*/
    /*----------------     Pack / Unpack (Ertl)     ----------------*/
    /*--------------------------------------------------------------*/

    /** Unpack 8-bit register to 64-bit hash prefix bitset. */
    static long unpack(byte register){
        return (4L | (register & 3)) << ((register >>> 2) - 2);
    }

    /** Pack 64-bit hash prefix bitset to 8-bit register. */
    static byte pack(long hashPrefix){
        int nlz=Long.numberOfLeadingZeros(hashPrefix)+1;
        return (byte)((-nlz << 2) | ((hashPrefix << nlz) >>> 62));
    }

    /*--------------------------------------------------------------*/
    /*----------------        Hash & Store          ----------------*/
    /*--------------------------------------------------------------*/

    @Override
    public final void hashAndStore(final long number){
        final long key=Tools.hash64shift(number^hashXor);
        // Ertl convention: top p bits = bucket index, NLZ of remaining (64-p) bits
        final int q=Long.numberOfLeadingZeros(buckets-1L); // q = 64 - p
        final int idx=(int)(key >>> q);
        final int nlz=Long.numberOfLeadingZeros(~(~key << -q)); // nlz in {0, ..., 64-p}
        final byte oldState=registers[idx];
        long hashPrefix=unpack(oldState);
        hashPrefix|=1L << (nlz + ~q); // bit at position (nlz + p - 1)
        final byte newState=pack(hashPrefix);
        if((newState&0xFF)>(oldState&0xFF)){
            registers[idx]=newState;
            lastCardinality=-1;
        }
    }

    /*--------------------------------------------------------------*/
    /*----------------      FGRA Estimator          ----------------*/
    /*--------------------------------------------------------------*/

    // Constants from Ertl's optimal FGRA estimator (b=2, d=2)
    static final double ETA_0=4.663135422063788;
    static final double ETA_1=2.1378502137958524;
    static final double ETA_2=2.781144650979996;
    static final double ETA_3=0.9824082545153715;
    static final double TAU=0.8194911375910897;
    static final double V_CONST=0.6118931496978437;

    private static final double POW_2_TAU=Math.pow(2., TAU);
    private static final double POW_2_MINUS_TAU=Math.pow(2., -TAU);
    private static final double POW_4_MINUS_TAU=Math.pow(4., -TAU);
    private static final double MINUS_INV_TAU=-1.0/TAU;
    private static final double ETA_X=ETA_0-ETA_1-ETA_2+ETA_3;
    private static final double ETA23X=(ETA_2-ETA_3)/ETA_X;
    private static final double ETA13X=(ETA_1-ETA_3)/ETA_X;
    private static final double ETA3012XX=(ETA_3*ETA_0-ETA_1*ETA_2)/(ETA_X*ETA_X);
    private static final double POW_4_MINUS_TAU_ETA_01=POW_4_MINUS_TAU*(ETA_0-ETA_1);
    private static final double POW_4_MINUS_TAU_ETA_23=POW_4_MINUS_TAU*(ETA_2-ETA_3);
    private static final double POW_4_MINUS_TAU_ETA_1=POW_4_MINUS_TAU*ETA_1;
    private static final double POW_4_MINUS_TAU_ETA_3=POW_4_MINUS_TAU*ETA_3;
    private static final double POW_2_MINUS_TAU_ETA_X=POW_2_MINUS_TAU*ETA_X;
    private static final double POW_2_MINUS_TAU_ETA_02=POW_2_MINUS_TAU*(ETA_0-ETA_2);
    private static final double POW_2_MINUS_TAU_ETA_13=POW_2_MINUS_TAU*(ETA_1-ETA_3);
    private static final double POW_2_MINUS_TAU_ETA_2=POW_2_MINUS_TAU*ETA_2;
    private static final double POW_2_MINUS_TAU_ETA_3=POW_2_MINUS_TAU*ETA_3;
    private static final double PHI_1=ETA_0/(POW_2_TAU*(2.*POW_2_TAU-1));
    private static final double P_INITIAL=ETA_X*(POW_4_MINUS_TAU/(2-POW_2_MINUS_TAU));

    /** Pre-computed estimation factors per precision p (index = p - 3). */
    static final double[] ESTIMATION_FACTORS={
        94.59941722950778, 455.6358404615186, 2159.476860400962,
        10149.51036338182, 47499.52712820488, 221818.76564766388,
        1034754.6840013304, 4824374.384717942, 2.2486750611989766E7,
        1.0479810199493326E8, 4.8837185623048025E8, 2.275794725435168E9,
        1.0604938814719946E10, 4.9417362104242645E10, 2.30276227770117E11,
        1.0730444972228585E12, 5.0001829613164E12, 2.329988778511272E13,
        1.0857295240912981E14, 5.059288069986326E14, 2.3575295235667005E15,
        1.0985627213141412E16, 5.119087674515589E16, 2.3853948339571715E17
    };

    /** Pre-computed per-register contributions for intermediate range (236 entries). */
    static final double[] REGISTER_CONTRIBUTIONS={
        0.8484061093359406, 0.38895829052007685, 0.5059986252327467, 0.17873835725405993,
        0.48074234060273024, 0.22040001471443574, 0.2867199572932749, 0.10128061935935387,
        0.2724086914332655, 0.12488785473931466, 0.16246750447680292, 0.057389829555353204,
        0.15435814343988866, 0.0707666752272979, 0.09206087452057209, 0.03251947467566813,
        0.08746577181824695, 0.0400993542020493, 0.05216553700867983, 0.018426892732996067,
        0.04956175987398336, 0.022721969094305374, 0.029559172293066274, 0.01044144713836362,
        0.02808376340530896, 0.012875216815740723, 0.01674946174724118, 0.005916560101748389,
        0.015913433441643893, 0.0072956356627506685, 0.009490944673308844, 0.0033525700962450116,
        0.009017216113341773, 0.004134011914931561, 0.0053779657012946284, 0.0018997062578498703,
        0.005109531310944485, 0.002342503834183061, 0.00304738001114257, 0.001076452918957914,
        0.0028952738727082267, 0.0013273605219527246, 0.0017267728074345586, 6.09963188753462E-4,
        0.0016405831157217021, 7.521379173550258E-4, 9.78461602292084E-4, 3.4563062172237723E-4,
        9.2962292270938E-4, 4.2619276177576713E-4, 5.544372155028133E-4, 1.958487477192352E-4,
        5.267631795945699E-4, 2.4149862146135835E-4, 3.141672858847145E-4, 1.1097608132071735E-4,
        2.9848602115777116E-4, 1.3684320663902123E-4, 1.7802030736817869E-4, 6.288368329501905E-5,
        1.6913464774658265E-4, 7.754107700464113E-5, 1.0087374230011362E-4, 3.563252169014952E-5,
        9.583875639268212E-5, 4.393801322487549E-5, 5.715927601779108E-5, 2.0190875207520577E-5,
        5.430624268457414E-5, 2.4897113642537945E-5, 3.2388833410757184E-5, 1.144099329232623E-5,
        3.0772185549154786E-5, 1.4107744575453657E-5, 1.8352865935237916E-5, 6.482944704957522E-6,
        1.7436805727319977E-5, 7.99403737572986E-6, 1.0399500462555932E-5, 3.67350727106242E-6,
        9.880422483694849E-6, 4.529755498675165E-6, 5.892791363067244E-6, 2.081562667074589E-6,
        5.5986600976661345E-6, 2.5667486794686803E-6, 3.339101736056405E-6, 1.1795003568090263E-6,
        3.1724346748254955E-6, 1.4544270182973653E-6, 1.8920745223756656E-6, 6.683541714686068E-7,
        1.7976340035771381E-6, 8.241391019206623E-7, 1.072128458850476E-6, 3.7871739159788393E-7,
        1.0186145159929963E-6, 4.6699164053601817E-7, 6.075127690181302E-7, 2.1459709360913574E-7,
        5.77189533646426E-7, 2.6461697039041317E-7, 3.442421115430427E-7, 1.2159967724530947E-7,
        3.27059699739513E-7, 1.4994302882644454E-7, 1.9506195985170504E-7, 6.890345650764188E-8,
        1.853256875916027E-7, 8.49639834530526E-8, 1.1053025444979778E-7, 3.904357664636507E-8,
        1.0501327589016596E-7, 4.814414208323267E-8, 6.263105916717392E-8, 2.2123721430020238E-8,
        5.9504908663745294E-8, 2.7280481949286693E-8, 3.548937430686624E-8, 1.2536224699555158E-8,
        3.371796684815404E-8, 1.545826061452554E-8, 2.0109761920695445E-8, 7.103548569567803E-9,
        1.910600846054063E-8, 8.759296176321385E-9, 1.139503111580109E-8, 4.0251673442004705E-9,
        1.082626247715867E-8, 4.963383100969499E-9, 6.456900615837058E-9, 2.28082795382416E-9,
        6.134612546958812E-9, 2.812460192131048E-9, 3.65874960227048E-9, 1.292412391857717E-9,
        3.476127720042246E-9, 1.5936574250689536E-9, 2.0732003554895977E-9, 7.323348470132607E-10,
        1.9697191686598677E-9, 9.030328662369446E-10, 1.2872203525845489E-9, 4.54696198313039E-10,
        1.2229703685228866E-9, 5.606801491206791E-10, 7.293928206826874E-10, 2.5764985922987735E-10,
        6.92986095905959E-10, 3.1770479284824887E-10, 4.1330443990824427E-10, 1.4599517261737423E-10,
        3.926748688923721E-10, 1.8002480658009348E-10, 2.3419555992885186E-10, 8.272696321778206E-11,
        2.225059832666067E-10, 1.0200957528418621E-10, 1.327049869160979E-10, 4.687655297461429E-11,
        1.2608118449008524E-10, 5.780288643182276E-11, 7.519618885068399E-11, 2.656221301145837E-11,
        7.144286571105751E-11, 3.2753529955811655E-11, 4.2609301647742677E-11, 1.5051259431302017E-11,
        4.0482511975524363E-11, 1.8559518231526075E-11, 2.4144210160882415E-11, 8.528672304925501E-12,
        2.293908229376684E-11, 1.0516598285774437E-11, 1.3681118012966618E-11, 4.832701981970378E-12,
        1.2998242223663023E-11, 5.959143881034847E-12, 7.752292944665042E-12, 2.7384108113817744E-12,
        7.365346997814574E-12, 3.376699844369893E-12, 4.392773006047039E-12, 1.5516979527951759E-12,
        4.173513269314059E-12, 1.9133791810691354E-12, 2.4891286772044455E-12, 8.792568765435867E-13,
        2.1582778227746425E-12, 9.89479027998621E-13, 1.2872203525845489E-12, 4.54696198313039E-13,
        1.2229703685228866E-12, 5.606801491206791E-13, 7.293928206826874E-13, 2.5764985922987735E-13,
        6.92986095905959E-13, 3.1770479284824887E-13, 4.1330443990824427E-13, 1.4599517261737423E-13,
        3.926748688923721E-13, 1.8002480658009348E-13, 2.3419555992885186E-13, 8.272696321778206E-14,
        2.225059832666067E-13, 1.0200957528418621E-13, 1.327049869160979E-13, 4.687655297461429E-14,
        1.2608118449008524E-13, 5.780288643182276E-14, 7.519618885068399E-14, 2.656221301145837E-14,
        7.144286571105751E-14, 3.2753529955811655E-14, 4.2609301647742677E-14, 1.5051259431302017E-14,
        4.0482511975524363E-14, 1.8559518231526075E-14, 2.4144210160882415E-14, 8.528672304925501E-15,
        2.293908229376684E-14, 1.0516598285774437E-14, 1.3681118012966618E-14, 4.832701981970378E-15,
        1.2998242223663023E-14, 5.959143881034847E-15, 7.752292944665042E-15, 2.7384108113817744E-15,
        7.365346997814574E-15, 3.376699844369893E-15, 4.392773006047039E-15, 1.5516979527951759E-15,
        4.173513269314059E-15, 1.9133791810691354E-15, 2.4891286772044455E-15, 8.792568765435867E-16
    };

    // ψ'(z) = psiPrime (Ertl's notation, divided by ETA_X)
    private static double psiPrime(double z, double zSquare){
        return (z+ETA23X)*(zSquare+ETA13X)+ETA3012XX;
    }

    // σ(z): small-range series (eq 20)
    private static double sigma(double z){
        if(z<=0.) return ETA_3;
        if(z>=1.) return Double.POSITIVE_INFINITY;
        double powZ=z, nextPowZ=powZ*powZ;
        double s=0, powTau=ETA_X;
        while(true){
            double oldS=s;
            double nextNextPowZ=nextPowZ*nextPowZ;
            s+=powTau*(powZ-nextPowZ)*psiPrime(nextPowZ, nextNextPowZ);
            if(!(s>oldS)) return s/z;
            powZ=nextPowZ; nextPowZ=nextNextPowZ; powTau*=POW_2_TAU;
        }
    }

    // φ(z): large-range series
    private static double phi(double z, double zSquare){
        if(z<=0.) return 0.;
        if(z>=1.) return PHI_1;
        double previousPowZ=zSquare, powZ=z;
        double nextPowZ=Math.sqrt(powZ);
        double p=P_INITIAL/(1.+nextPowZ);
        double ps=psiPrime(powZ, previousPowZ);
        double s=nextPowZ*(ps+ps)*p;
        while(true){
            previousPowZ=powZ; powZ=nextPowZ;
            double oldS=s;
            nextPowZ=Math.sqrt(powZ);
            double nextPs=psiPrime(powZ, previousPowZ);
            p*=POW_2_MINUS_TAU/(1.+nextPowZ);
            s+=nextPowZ*((nextPs+nextPs)-(powZ+nextPowZ)*ps)*p;
            if(!(s>oldS)) return s;
            ps=nextPs;
        }
    }

    private static double smallRangeEstimate(long c0, long c4, long c8, long c10, long m){
        long alpha=m+3*(c0+c4+c8+c10);
        long beta=m-c0-c4;
        long gamma=4*c0+2*c4+3*c8+c10;
        double quadRootZ=(Math.sqrt((double)(beta*beta+4*alpha*gamma))-beta)/(2*alpha);
        double rootZ=quadRootZ*quadRootZ;
        return rootZ*rootZ;
    }

    private static double largeRangeEstimate(long c4w0, long c4w1, long c4w2, long c4w3, long m){
        long alpha=m+3*(c4w0+c4w1+c4w2+c4w3);
        long beta=c4w0+c4w1+2*(c4w2+c4w3);
        long gamma=m+2*c4w0+c4w2-c4w3;
        return Math.sqrt((Math.sqrt((double)(beta*beta+4*alpha*gamma))-beta)/(2*alpha));
    }

    /*--------------------------------------------------------------*/
    /*----------------       Cardinality            ----------------*/
    /*--------------------------------------------------------------*/

    /** Raw FGRA estimate (unclamped, for fair comparison in rawEstimates). */
    public double fgraEstimatePublic(){return fgraEstimate();}
    private double fgraEstimate(){
        final int m=buckets;
        final int p=bucketBits;

        int c0=0, c4=0, c8=0, c10=0;
        int c4w0=0, c4w1=0, c4w2=0, c4w3=0;
        double sum=0;
        final int off=(p<<2)+4;

        for(int i=0; i<m; i++){
            int r=registers[i]&0xFF;
            int r2=r-off;
            if(r2<0){
                if(r2<-8) c0++;
                if(r2==-8) c4++;
                if(r2==-4) c8++;
                if(r2==-2) c10++;
            }else if(r<252){
                sum+=REGISTER_CONTRIBUTIONS[r2];
            }else{
                if(r==252) c4w0++;
                if(r==253) c4w1++;
                if(r==254) c4w2++;
                if(r==255) c4w3++;
            }
        }

        // Small range correction
        if(c0>0 || c4>0 || c8>0 || c10>0){
            double z=smallRangeEstimate(c0, c4, c8, c10, m);
            if(c0>0) sum+=c0*sigma(z);
            if(c4>0) sum+=c4*POW_2_MINUS_TAU_ETA_X*psiPrime(z, z*z);
            if(c8>0) sum+=c8*(z*POW_4_MINUS_TAU_ETA_01+POW_4_MINUS_TAU_ETA_1);
            if(c10>0) sum+=c10*(z*POW_4_MINUS_TAU_ETA_23+POW_4_MINUS_TAU_ETA_3);
        }

        // Large range correction
        if(c4w0>0 || c4w1>0 || c4w2>0 || c4w3>0){
            double z=largeRangeEstimate(c4w0, c4w1, c4w2, c4w3, m);
            double rootZ=Math.sqrt(z);
            double s2=phi(rootZ, z)*(c4w0+c4w1+c4w2+c4w3);
            s2+=z*(1+rootZ)*(c4w0*ETA_0+c4w1*ETA_1+c4w2*ETA_2+c4w3*ETA_3);
            s2+=rootZ*((c4w0+c4w1)*(z*POW_2_MINUS_TAU_ETA_02+POW_2_MINUS_TAU_ETA_2)
                      +(c4w2+c4w3)*(z*POW_2_MINUS_TAU_ETA_13+POW_2_MINUS_TAU_ETA_3));
            sum+=s2*Math.pow(POW_2_MINUS_TAU, 65-p)/((1+rootZ)*(1+z));
        }

        return ESTIMATION_FACTORS[p-3]*Math.pow(sum, MINUS_INV_TAU);
    }

    /**
     * Static FGRA estimator for use by other classes (e.g., UDLL6).
     * Takes an array of 8-bit Ertl-format registers and precision parameter p.
     */
    public static double fgraEstimateStatic(byte[] regs, int p){
        final int m=regs.length;
        int c0=0, c4=0, c8=0, c10=0;
        int c4w0=0, c4w1=0, c4w2=0, c4w3=0;
        double sum=0;
        final int off=(p<<2)+4;
        for(int i=0; i<m; i++){
            int r=regs[i]&0xFF;
            int r2=r-off;
            if(r2<0){
                if(r2<-8) c0++;
                if(r2==-8) c4++;
                if(r2==-4) c8++;
                if(r2==-2) c10++;
            }else if(r<252){
                sum+=REGISTER_CONTRIBUTIONS[r2];
            }else{
                if(r==252) c4w0++;
                if(r==253) c4w1++;
                if(r==254) c4w2++;
                if(r==255) c4w3++;
            }
        }
        if(c0>0 || c4>0 || c8>0 || c10>0){
            double z=smallRangeEstimate(c0, c4, c8, c10, m);
            if(c0>0) sum+=c0*sigma(z);
            if(c4>0) sum+=c4*POW_2_MINUS_TAU_ETA_X*psiPrime(z, z*z);
            if(c8>0) sum+=c8*(z*POW_4_MINUS_TAU_ETA_01+POW_4_MINUS_TAU_ETA_1);
            if(c10>0) sum+=c10*(z*POW_4_MINUS_TAU_ETA_23+POW_4_MINUS_TAU_ETA_3);
        }
        if(c4w0>0 || c4w1>0 || c4w2>0 || c4w3>0){
            double z=largeRangeEstimate(c4w0, c4w1, c4w2, c4w3, m);
            double rootZ=Math.sqrt(z);
            double s2=phi(rootZ, z)*(c4w0+c4w1+c4w2+c4w3);
            s2+=z*(1+rootZ)*(c4w0*ETA_0+c4w1*ETA_1+c4w2*ETA_2+c4w3*ETA_3);
            s2+=rootZ*((c4w0+c4w1)*(z*POW_2_MINUS_TAU_ETA_02+POW_2_MINUS_TAU_ETA_2)
                      +(c4w2+c4w3)*(z*POW_2_MINUS_TAU_ETA_13+POW_2_MINUS_TAU_ETA_3));
            sum+=s2*Math.pow(POW_2_MINUS_TAU, 65-p)/((1+rootZ)*(1+z));
        }
        if(p-3<0 || p-3>=ESTIMATION_FACTORS.length){return 0;}
        return ESTIMATION_FACTORS[p-3]*Math.pow(sum, MINUS_INV_TAU);
    }

    @Override
    public final long cardinality(){
        if(lastCardinality>=0) return lastCardinality;
        long card=Math.max(0, Math.round(fgraEstimate()));
        card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
        lastCardinality=card;
        return card;
    }

    /*--------------------------------------------------------------*/
    /*----------------       Framework Glue         ----------------*/
    /*--------------------------------------------------------------*/

    @Override
    public final void add(CardinalityTracker log){
        assert(log.getClass()==this.getClass());
        final ErtlULL other=(ErtlULL)log;
        added+=other.added;
        lastCardinality=-1;
        // ULL merge: unpack → OR hashPrefixes → repack.
        // Per-register max is WRONG because sub-bits encode history
        // that may differ between instances for the same bucket.
        for(int i=0; i<registers.length; i++){
            long hp=unpack(registers[i]) | unpack(other.registers[i]);
            registers[i]=pack(hp);
        }
    }

    public byte[] getRegisters(){return registers;}
    public double occupancy(){
        int filled=0;
        for(int i=0; i<buckets; i++){if(registers[i]!=0) filled++;}
        return (double)filled/buckets;
    }

    @Override public final float[] compensationFactorLogBucketsArray(){return null;}

    @Override
    public double[] rawEstimates(){
        // ErtlULL uses FGRA estimator directly — CardinalityStats estimators don't apply
        // because the register encoding is incompatible with standard NLZ conventions.
        // Fill key slots with unclamped FGRA estimate for fair comparison.
        final int total=11+4+CardinalityStats.NUM_DLC_TIERS;
        final double[] r=new double[total];
        final double fgra=fgraEstimate();
        r[0]=fgra;  // Mean
        r[1]=fgra;  // HMean
        r[4]=fgra;  // HLL
        r[6]=fgra;  // Hybrid
        r[8]=fgra;  // DThHyb
        return r;
    }

    // No static CF loading — ErtlULL uses its own FGRA estimator, not CF tables
    public static void setCFMatrix(float[][] matrix, int buckets){}
}
