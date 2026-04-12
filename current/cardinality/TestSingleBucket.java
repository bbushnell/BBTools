package cardinality;
import shared.Tools;

/** Minimal single-bucket simulation of the 3-rule history mechanism. */
public class TestSingleBucket {
    static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*65536);

    public static void main(String[] args){
        // Single bucket, feed random hashes, track history
        int storedTier=0; // absolute tier
        int hist=0;
        boolean occupied=false;
        long h0count=0, h1count=0, totalObs=0;
        long advSet=0, advClr=0, dNeg1Set=0, dNeg1Skip=0;
        java.util.Random rng=new java.util.Random(42);

        for(int i=0; i<50_000_000; i++){
            long key=Tools.hash64shift(rng.nextLong());
            int rawNlz=Long.numberOfLeadingZeros(key);
            int mant=0;
            if(rawNlz<47){
                int mBits=(int)((key>>>(46-rawNlz))&0xFFFF);
                mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
            }
            int elemTier=(2*rawNlz+mant)/3;

            if(!occupied){
                storedTier=elemTier; hist=0; occupied=true;
            }else{
                int delta=elemTier-storedTier;
                if(delta>0){
                    // Rule 1 & 2: advance
                    hist=(delta==1) ? 1 : 0;
                    if(delta==1) advSet++; else advClr++;
                    storedTier=elemTier;
                }else if(delta==-1){
                    // Rule 3: set history
                    if(hist==0){dNeg1Set++; hist=1;}
                    else{dNeg1Skip++;}
                }
                // delta==0 or delta<=-2: no action
            }

            if(occupied){
                totalObs++;
                if(hist==0) h0count++; else h1count++;
            }
        }
        System.err.println("=== Single-Bucket Simulation ===");
        System.err.println("Final tier: "+storedTier+" hist: "+hist);
        System.err.println("P(0) = "+String.format("%.4f", (double)h0count/totalObs));
        System.err.println("advance set (delta=1): "+advSet);
        System.err.println("advance clear (delta>=2): "+advClr);
        System.err.println("delta=-1 set: "+dNeg1Set);
        System.err.println("delta=-1 skip: "+dNeg1Skip);
        long totalAdv=advSet+advClr;
        System.err.println("advance delta=1 rate: "+String.format("%.1f%%", 100.0*advSet/totalAdv));
        System.err.println("delta=-1 per advance: "+String.format("%.2f", (double)(dNeg1Set+dNeg1Skip)/totalAdv));
    }
}
