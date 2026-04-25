package cardinality;
public class TestAVDLL32Hist {
    public static void main(String[] args) {
        ArithmeticVariableDynamicLogLog32 e = new ArithmeticVariableDynamicLogLog32(1536, 31, 42L, 0f);
        for (long i = 0; i < 5000; i++) {
            e.hashAndStore(i);
        }
        e.rawEstimates();
        CardStats cs = e.consumeLastSummarized();
        System.out.println("--- AVDLL32 ---");
        System.out.println("globalNLZ+1: " + e.getMinZeros());
        System.out.println("DLC: " + cs.dlcRaw() + " HC: " + cs.hc() + " LDLC: " + cs.ldlc());
        System.out.println("SBS: " + cs.dlcSbs() + " LC: " + cs.lcRaw());
        System.out.println("Mean+H: " + cs.meanHistCF() + " HybPlus2: " + cs.hybridPlus2());
        
        ArithmeticUltraDynamicLogLog32 e2 = new ArithmeticUltraDynamicLogLog32(1536, 31, 42L, 0f);
        for (long i = 0; i < 5000; i++) {
            e2.hashAndStore(i);
        }
        e2.rawEstimates();
        CardStats cs2 = e2.consumeLastSummarized();
        System.out.println("\n--- AUDLL32 ---");
        System.out.println("globalNLZ+1: " + e2.getMinZeros());
        System.out.println("DLC: " + cs2.dlcRaw() + " HC: " + cs2.hc() + " LDLC: " + cs2.ldlc());
        System.out.println("SBS: " + cs2.dlcSbs() + " LC: " + cs2.lcRaw());
        System.out.println("Mean+H: " + cs2.meanHistCF() + " HybPlus2: " + cs2.hybridPlus2());
    }
}
