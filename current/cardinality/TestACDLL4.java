package cardinality;
public class TestACDLL4 {
    public static void main(String[] args) {
        ArithmeticCompressedDynamicLogLog4 t = new ArithmeticCompressedDynamicLogLog4(1856, 31, 42L, 0f);
        java.util.Random rng = new java.util.Random(123);
        for (int i = 0; i < 50000; i++) {
            t.add(rng.nextLong());
        }
        // Call rawEstimates to trigger summarize + lastSummarized
        double[] est = t.rawEstimates();
        CardStats cs = t.consumeLastSummarized();
        if (cs == null) {
            System.out.println("CardStats is NULL!");
            return;
        }
        System.out.println("hc=" + cs.hc());
        System.out.println("ldlc=" + cs.ldlc());
        System.out.println("meanHistCF=" + cs.meanHistCF());
        System.out.println("hybridPlus2=" + cs.hybridPlus2());
        System.out.println("dlcSbs=" + cs.dlcSbs());
        System.out.println("dlcRaw=" + cs.dlcRaw());
        System.out.println("historyCoverage=" + cs.historyCoverage());
    }
}
