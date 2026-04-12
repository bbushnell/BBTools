package cardinality;
import shared.Tools;
public class TestCDLL4Hist {
    public static void main(String[] args) {
        CompressedDynamicLogLog4.DISABLE_EEMASK = false;
        if(args.length > 0 && args[0].equals("noeemask")) CompressedDynamicLogLog4.DISABLE_EEMASK = true;
        CompressedDynamicLogLog4.DEBUG_HIST = false;
        CompressedDynamicLogLog4 c = new CompressedDynamicLogLog4(512, 31, 12345L, 0);
        java.util.Random r = new java.util.Random(42);
        // Feed lots of elements
        for(int i = 0; i < 5000000; i++) {
            c.add(r.nextLong());
        }
        // Force summarize to get hist distribution
        double[] est = c.rawEstimates();
        int h0=0, h1=0, phantom=0;
        for(int i=0; i<512; i++){
            // Can't read nibbles from outside... use summarize debug
        }
        System.err.println("minZeros=" + c.getMinZeros() + " card=" + c.cardinality());
        CompressedDynamicLogLog4.printDebugCounters();
        // Print ratios
        long advSet = CompressedDynamicLogLog4.cntAdvanceSet;
        long advClr = CompressedDynamicLogLog4.cntAdvanceClear;
        long d1Set = CompressedDynamicLogLog4.cntDeltaNeg1Set;
        long d1Skip = CompressedDynamicLogLog4.cntDeltaNeg1Skip;
        long totalAdv = advSet + advClr;
        System.err.println("--- Ratios ---");
        System.err.println("  advance delta=1 rate: " + String.format("%.1f%%", 100.0*advSet/totalAdv));
        System.err.println("  advance delta>=2 rate: " + String.format("%.1f%%", 100.0*advClr/totalAdv));
        System.err.println("  delta=-1 hit rate: " + String.format("%.1f%%", 100.0*d1Set/(d1Set+d1Skip)));
    }
}
