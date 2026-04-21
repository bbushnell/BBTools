package cardinality;

/** Quick test: verify ACDLL5 word packing round-trips and history preservation. */
public class TestACDLL5 {
    public static void main(String[] args) {
        ArithmeticCompressedDynamicLogLog5 a = new ArithmeticCompressedDynamicLogLog5(64, 31, -1, 0);
        
        // Test 1: write and read back all 6 positions in a word
        System.out.println("=== Test 1: Read/Write round-trip ===");
        for (int pos = 0; pos < 6; pos++) {
            int bucket = pos; // all in word 0
            // Write exponent and history
            a.writeBucket(bucket, pos + 1, pos % 4);
        }
        for (int pos = 0; pos < 6; pos++) {
            int bucket = pos;
            int exp = a.readExponent(bucket);
            int hist = a.readHistory(bucket);
            System.out.printf("  bucket %d: exp=%d (expect %d), hist=%d (expect %d) %s%n",
                bucket, exp, pos + 1, hist, pos % 4,
                (exp == pos + 1 && hist == pos % 4) ? "OK" : "FAIL");
        }
        
        // Test 2: write exponent without clobbering history
        System.out.println("=== Test 2: writeExponent preserves history ===");
        ArithmeticCompressedDynamicLogLog5 b = new ArithmeticCompressedDynamicLogLog5(64, 31, -1, 0);
        b.writeBucket(0, 3, 2); // exp=3, hist=2
        b.writeExponent(0, 5);  // change exp to 5, should keep hist=2
        int hist = b.readHistory(0);
        System.out.printf("  After writeExponent: exp=%d hist=%d (expect 5, 2) %s%n",
            b.readExponent(0), hist, (b.readExponent(0)==5 && hist==2) ? "OK" : "FAIL");
        
        // Test 3: writeHistory preserves exponent
        b.writeHistory(0, 3); // change hist to 3, should keep exp=5
        System.out.printf("  After writeHistory: exp=%d hist=%d (expect 5, 3) %s%n",
            b.readExponent(0), b.readHistory(0),
            (b.readExponent(0)==5 && b.readHistory(0)==3) ? "OK" : "FAIL");
        
        // Test 4: countAndDecrement preserves history
        System.out.println("=== Test 3: countAndDecrement history preservation ===");
        ArithmeticCompressedDynamicLogLog5 c = new ArithmeticCompressedDynamicLogLog5(64, 31, -1, 0);
        // Set up: bucket 0: exp=5 hist=3, bucket 1: exp=2 hist=1, bucket 2: exp=0 hist=2
        c.writeBucket(0, 5, 3);
        c.writeBucket(1, 2, 1);
        c.writeBucket(2, 0, 2); // floor-level bucket with history
        System.out.printf("  Before: b0=(exp=%d,hist=%d) b1=(exp=%d,hist=%d) b2=(exp=%d,hist=%d)%n",
            c.readExponent(0), c.readHistory(0),
            c.readExponent(1), c.readHistory(1),
            c.readExponent(2), c.readHistory(2));
        
        // Make methods accessible for test
        // Actually countAndDecrement is private... let me just run hashAndStore and check
        
        // Test 5: Hash many values, compare CDLL5 vs ACDLL5 bucket states
        System.out.println("=== Test 4: CDLL5 vs ACDLL5 bucket comparison ===");
        CompressedDynamicLogLog5 cdll5 = new CompressedDynamicLogLog5(64, 31, 42, 0);
        ArithmeticCompressedDynamicLogLog5 acdll5 = new ArithmeticCompressedDynamicLogLog5(64, 31, 42, 0);
        
        for (long i = 0; i < 10000; i++) {
            cdll5.hashAndStore(i);
            acdll5.hashAndStore(i);
        }
        
        int mismatches = 0;
        for (int i = 0; i < 64; i++) {
            int cReg = cdll5.readBucket(i);
            int cExp = (cReg >>> 2) & 0x7;
            int cHist = cReg & 0x3;
            int aExp = acdll5.readExponent(i);
            int aHist = acdll5.readHistory(i);
            if (cExp != aExp || cHist != aHist) {
                mismatches++;
                if (mismatches <= 20) {
                    System.out.printf("  MISMATCH bucket %d: CDLL5(exp=%d,hist=%d) ACDLL5(exp=%d,hist=%d)%n",
                        i, cExp, cHist, aExp, aHist);
                }
            }
        }
        System.out.printf("  Total mismatches: %d / 64 buckets (minZeros: CDLL5=%d ACDLL5=%d)%n",
            mismatches, cdll5.getMinZeros(), acdll5.getMinZeros());
    }
}
