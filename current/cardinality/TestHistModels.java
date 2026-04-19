package cardinality;
import shared.Tools;

/** Compare 3-rule carry model vs bitmap model for 1-bit history. */
public class TestHistModels {
    static final int MT=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

    static int toTier(long key){
        int rawNlz=Long.numberOfLeadingZeros(key);
        final int mBits;
        if(rawNlz>=43){
            mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
        }else{
            mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
        }
        int mant=(mBits>=MT) ? 1 : 0;
        return (2*rawNlz+mant)/3;
    }

    public static void main(String[] args){
        java.util.Random rng=new java.util.Random(42);
        // Model A: 3-rule carry
        int tierA=0; int histA=0; boolean occA=false;
        long a0=0, a1=0;
        // Model B: bitmap - track which half-NLZ positions observed
        // For 1-bit history: history = have we seen halfNlz at (3*storedTier - 1)?
        // i.e., top position of previous tier
        int tierB=0; boolean occB=false;
        long bitmapB=0; // bit positions = halfNlz values observed
        long b0=0, b1=0;

        int N=50_000_000;
        for(int i=0; i<N; i++){
            long key=Tools.hash64shift(rng.nextLong());
            int elemTier=toTier(key);
            int rawNlz=Long.numberOfLeadingZeros(key);
            final int mBits2;
            if(rawNlz>=43){
                mBits2=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
            }else{
                mBits2=(int)((key>>>(42-rawNlz))&0xFFFFF);
            }
            int mant=(mBits2>=MT) ? 1 : 0;
            int halfNlz=2*rawNlz+mant;

            // Model A: 3 rules
            if(!occA){
                tierA=elemTier; histA=0; occA=true;
            }else{
                int delta=elemTier-tierA;
                if(delta>0){
                    histA=(delta==1) ? 1 : 0;
                    tierA=elemTier;
                }else if(delta==-1){
                    histA=1;
                }
            }
            if(occA){if(histA==0) a0++; else a1++;}

            // Model B: bitmap-like
            // Set bit for this halfNlz. Extract history from sub-tier boundary.
            if(!occB){
                tierB=elemTier; bitmapB=(1L<<halfNlz); occB=true;
            }else{
                bitmapB|=(1L<<Math.min(halfNlz,63));
                if(elemTier>tierB){
                    tierB=elemTier;
                }
            }
            if(occB){
                // History: have we seen any halfNlz in the previous tier?
                // Previous tier = halfNlz range [3*(tierB-1), 3*tierB - 1]
                int prevLo=3*(tierB-1);
                int prevHi=3*tierB-1;
                boolean histBitmapSet=false;
                for(int h=Math.max(0,prevLo); h<=Math.min(63,prevHi); h++){
                    if((bitmapB&(1L<<h))!=0){histBitmapSet=true; break;}
                }
                if(histBitmapSet) b1++; else b0++;
            }
        }
        System.err.println("Model A (3-rule carry):");
        System.err.println("  P(0) = "+String.format("%.4f", (double)a0/(a0+a1)));
        System.err.println("Model B (bitmap):");
        System.err.println("  P(0) = "+String.format("%.4f", (double)b0/(b0+b1)));
    }
}
