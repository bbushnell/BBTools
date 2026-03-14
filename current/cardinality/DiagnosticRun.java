package cardinality;

import rand.FastRandomXoshiro;

/**
 * Single-DDL diagnostic: runs DLL4 at 512 buckets using the 2048-bucket v5 CF table,
 * and for each threshold where trueCard > 1000, asserts corrected Mean is within 15%.
 * On failure, prints the raw Mean, DLC estimate, CF factor looked up, and corrected Mean.
 */
public class DiagnosticRun {

    public static void main(String[] args) throws Exception {
        final int buckets = 512;
        final long maxTrue = buckets * 200L;
        final double reportFrac = 0.01;
        final String cfFile = "/mnt/c/playground/Chloe/dll4_cf_v5b.tsv";

        // Force class init of DLL4 FIRST (triggers its default v4 CF load, resets v1Buckets=0),
        // then load the v5 CF table so it correctly overwrites the default.
        final CardinalityTracker ddl = DDLCalibrationDriver.makeInstance("dll4", buckets, 31, 42L, 0);
        CorrectionFactor.initialize(cfFile, buckets);
        CorrectionFactor.USE_CORRECTION = true;
        System.err.println("v1Buckets=" + CorrectionFactor.v1Buckets
                + "  tableVersion=" + CorrectionFactor.tableVersion
                + "  keyScale=" + ((double)CorrectionFactor.v1Buckets / buckets));
        final FastRandomXoshiro rng = new FastRandomXoshiro(99999L);

        // Compute thresholds
        final long[] thresholds = DDLCalibrationDriver.computeThresholds(maxTrue, reportFrac);
        int ti = 0;

        int failures = 0;
        int checked = 0;

        for (long trueCard = 1; trueCard <= maxTrue; trueCard++) {
            ddl.hashAndStore(rng.nextLong());

            if (ti < thresholds.length && trueCard >= thresholds[ti]) {

                // Get raw estimates (CF off)
                CorrectionFactor.USE_CORRECTION = false;
                final double[] raw = ddl.rawEstimates();

                // Get corrected estimates (CF on)
                CorrectionFactor.USE_CORRECTION = true;
                final double[] corr = ddl.rawEstimates();

                final double rawMean   = raw[0];   // r[0] = Mean (no CF)
                final double corrMean  = corr[0];  // r[0] = Mean * cf
                final double rawDLC    = raw[11];  // r[11] = dlcLogSpace025() (no CF)
                final double corrDLC   = corr[11]; // r[11] = dlcLogSpace025() * cf
                final double cfFactor  = (rawMean > 0) ? corrMean / rawMean : 0;
                final double keyUsed   = rawDLC * ((double)CorrectionFactor.v1Buckets / buckets);
                final double relErr    = (corrMean - trueCard) / (double)trueCard;

                if (trueCard > 1000) {
                    checked++;
                    if (Math.abs(relErr) > 0.15) {
                        failures++;
                        System.out.printf(
                            "FAIL trueCard=%6d  rawMean=%8.1f  corrMean=%8.1f  relErr=%+7.3f" +
                            "  rawDLC=%8.1f  keyUsed=%8.1f  cfFactor=%.5f  corrDLC=%8.1f%n",
                            trueCard, rawMean, corrMean, relErr,
                            rawDLC, keyUsed, cfFactor, corrDLC);
                    }
                }

                ti++;
                if (ti >= thresholds.length) break;
            }
        }

        System.out.printf("%nChecked %d thresholds above trueCard=1000.  Failures (>15%%): %d%n",
                checked, failures);
    }
}
