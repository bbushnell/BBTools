package clade;

/**
 * Represents PreClade v2.0 format data and converts it to proper Clade objects.
 * PreClade format contains raw k-mer counts from k=1 through k=5 that need
 * canonical mapping and GC compensation to match the database format.
 *
 * PreClade v2.0 Format:
 * //PreClade Format 2.0
 * #
 * sequence_name
 * 1-mers (5 values: A,C,G,T,N)
 * 2-mers (16 values)
 * 3-mers (64 values)
 * 4-mers (256 values)
 * 5-mers (1024 values)
 *
 * @author Chloe
 * @date September 28, 2025
 */
public class PreClade {

    /** Sequence name from PreClade entry */
    public final String name;

    /** Raw k-mer counts indexed by k-mer length [1-5] */
    public final long[][] rawCounts;

    /** Total bases calculated from 1-mer counts */
    public final long bases;

    /**
     * Constructs a PreClade from pre-parsed k-mer count arrays.
     * This constructor avoids String operations for maximum efficiency.
     *
     * @param name_ Sequence name
     * @param countsArray Pre-parsed k-mer counts [1-5][counts]
     *                    Index 0 unused, 1=1-mers(5), 2=2-mers(16), 3=3-mers(64), 4=4-mers(256), 5=5-mers(1024)
     */
    public PreClade(String name_, long[][] countsArray) {
        this.name = name_;
        this.rawCounts = countsArray;
        // Calculate total bases from 1-mers (A,C,G,T only, exclude N)
        this.bases = countsArray[1][0] + countsArray[1][1] + countsArray[1][2] + countsArray[1][3];
    }

    /**
     * Constructs a PreClade from parsed data lines.
     *
     * @param name_ Sequence name
     * @param monomersLine 1-mers line (5 comma-separated values: A,C,G,T,N)
     * @param dimersLine 2-mers line (16 comma-separated values)
     * @param trimersLine 3-mers line (64 comma-separated values)
     * @param tetramersLine 4-mers line (256 comma-separated values)
     * @param pentamersLine 5-mers line (1024 comma-separated values)
     * @throws IllegalArgumentException if format is invalid
     */
    public PreClade(String name_, String monomersLine, String dimersLine,
                   String trimersLine, String tetramersLine, String pentamersLine) {
        this.name = name_;
        this.rawCounts = new long[6][]; // Index 0 unused, 1-5 for k=1 through k=5

        // Parse 1-mers (A,C,G,T,N counts)
        this.rawCounts[1] = parseCountsLine(monomersLine, 5, "1-mer");

        // Parse 2-mers through 5-mers
        this.rawCounts[2] = parseCountsLine(dimersLine, 16, "2-mer");
        this.rawCounts[3] = parseCountsLine(trimersLine, 64, "3-mer");
        this.rawCounts[4] = parseCountsLine(tetramersLine, 256, "4-mer");
        this.rawCounts[5] = parseCountsLine(pentamersLine, 1024, "5-mer");

        // Calculate total bases from 1-mers (A,C,G,T only, exclude N)
        this.bases = rawCounts[1][0] + rawCounts[1][1] + rawCounts[1][2] + rawCounts[1][3];
    }

    /**
     * Parses a comma-delimited line of k-mer counts.
     *
     * @param line Comma-delimited counts
     * @param expectedCount Expected number of counts
     * @param kmerType Description for error messages
     * @return Array of parsed counts
     * @throws IllegalArgumentException if wrong number of counts
     */
    private long[] parseCountsLine(String line, int expectedCount, String kmerType) {
        String[] countStrs = line.split(",");
        if(countStrs.length != expectedCount) {
            throw new IllegalArgumentException("Expected " + expectedCount + " " + kmerType +
                " counts, got " + countStrs.length);
        }

        long[] counts = new long[expectedCount];
        for(int i = 0; i < expectedCount; i++) {
            counts[i] = Long.parseLong(countStrs[i].trim());
        }
        return counts;
    }

    /**
     * Converts this PreClade to a proper Clade object with canonical k-mer mapping
     * and GC compensation applied. This ensures PreClade signatures match the
     * canonical database format used by standard Clade processing.
     *
     * @return Clade object with properly processed k-mer signatures
     */
    public Clade toClade() {
        // Create Clade with dummy taxonomic info (will be determined by classification)
        Clade clade = new Clade(-1, -1, name);

        // Set basic statistics
        clade.bases = this.bases;
        clade.contigs = 1; // Assume single sequence

        // Copy 1-mers directly (no canonical mapping needed for ACGTN)
        System.arraycopy(rawCounts[1], 0, clade.counts[1], 0, 5);

        // Copy 2-mers directly (special case - no canonical mapping in CladeObject)
        System.arraycopy(rawCounts[2], 0, clade.counts[2], 0, 16);

        // Apply canonical k-mer mapping for k=3 through k=5
        for(int k = 3; k <= 5; k++) {
            applyCanonicalMapping(rawCounts[k], clade.counts[k], k);
        }

        // Calculate derived statistics and complete the Clade
        clade.finish();

        return clade;
    }

    /**
     * Applies canonical k-mer mapping using CladeObject.remapMatrix.
     * This maps raw k-mer indices to canonical indices, combining
     * forward and reverse complement k-mers.
     *
     * @param rawKmerCounts Raw k-mer counts indexed by lexicographic order
     * @param canonicalCounts Output array for canonical k-mer counts
     * @param k K-mer length (3, 4, or 5)
     */
    private void applyCanonicalMapping(long[] rawKmerCounts, long[] canonicalCounts, int k) {
        int[] remap = CladeObject.remapMatrix[k];
        if(remap == null) {
            throw new IllegalStateException("No remap matrix for k=" + k);
        }

        // Initialize canonical counts to zero
        java.util.Arrays.fill(canonicalCounts, 0);

        // Map raw k-mers to canonical indices and accumulate counts
        for(int i = 0; i < rawKmerCounts.length && i < remap.length; i++) {
            int canonicalIndex = remap[i];
            if(canonicalIndex >= 0 && canonicalIndex < canonicalCounts.length) {
                canonicalCounts[canonicalIndex] += rawKmerCounts[i];
            }
        }
    }

    @Override
    public String toString() {
        return "PreClade{name='" + name + "', bases=" + bases + "}";
    }
}