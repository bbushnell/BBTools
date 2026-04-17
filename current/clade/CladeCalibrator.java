package clade;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import bin.BinObject;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.PreParser;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tax.TaxTree;

/**
 * Builds QuickClade confidence-calibration tables from hit TSVs.
 *
 * Input: one or more QuickClade machine-format TSVs produced via
 *   quickclade.sh ... percontig hits=N format=machine out=hits.tsv
 * Each line is one query-hit pair with columns (tab-separated):
 *   QueryName, Q_GC, Q_Bases, Q_Contigs, RefName, R_TaxID, R_GC, R_Bases,
 *   R_Contigs, R_Level, GCdif, STRdif, k3dif, k4dif, k5dif, [lineage]
 * Query TaxID is parsed from QueryName via bin.BinObject.parseTaxID
 * (patterns tid_NNN or tid|NNN| both accepted).
 *
 * Output: one TSV per taxonomic level (species, genus, family, order, class,
 * phylum, kingdom, domain). Columns: len, k5difBin, k5difLow, correct, total, prob.
 *
 * @author Ady
 */
public class CladeCalibrator {

    public static void main(String[] args) {
        Timer t = new Timer();
        CladeCalibrator cc = new CladeCalibrator(args);
        cc.process(t);
    }

    public CladeCalibrator(String[] args) {
        {
            PreParser pp = new PreParser(args, CladeCalibrator.class, false);
            args = pp.args;
            outstream = pp.outstream;
        }

        ArrayList<String> in = new ArrayList<>();
        String outPrefix = "probCorrect";
        int k5difBins = 500;
        float k5difMax = 1.0f;
        String treeFile = null;
        rawMode = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String[] split = arg.split("=", 2);
            String a = split[0].toLowerCase();
            String b = split.length > 1 ? split[1] : null;

            if (a.equals("in") || a.equals("in1")) {
                for (String s : b.split(",")) { in.add(s); }
            } else if (a.equals("out")) {
                outPrefix = b;
            } else if (a.equals("bins") || a.equals("k5difbins")) {
                k5difBins = Integer.parseInt(b);
            } else if (a.equals("k5difmax") || a.equals("maxk5dif")) {
                k5difMax = Float.parseFloat(b);
            } else if (a.equals("tree") || a.equals("taxtree")) {
                treeFile = b;
            } else if (a.equals("raw")) {
                rawMode = parse.Parse.parseBoolean(b);
            } else if (b == null && !arg.startsWith("-")) {
                in.add(arg);
            } else {
                throw new RuntimeException("Unknown parameter '" + arg + "'");
            }
        }

        if (in.isEmpty()) {
            throw new RuntimeException("No input files specified (use in=file.tsv)");
        }

        this.inputs = in.toArray(new String[0]);
        this.outPrefix = outPrefix;
        this.k5difBins = k5difBins;
        this.k5difMax = k5difMax;
        this.binWidth = k5difMax / k5difBins;
        this.treeFile = treeFile;
        if (rawMode) { outstream.println("Raw mode: outputting per-hit features."); }
    }

    public void process(Timer t) {
        outstream.println("Loading TaxTree...");
        if (treeFile != null) {
            tree = TaxTree.loadTaxTree(treeFile, null, true, false);
        } else {
            tree = TaxTree.getTree();
        }
        if (tree == null) {
            throw new RuntimeException("Failed to load TaxTree. Use tree=path.taxtree.gz");
        }
        outstream.println("TaxTree loaded: " + tree.nodes.length + " nodes");

        // Length bucket schedule (matches QuickCladeDesign.md §4).
        lengthSchedule = new int[]{2500, 5000, 10000, 20000, 40000, 80000,
                160000, 320000, 640000, 1280000, 2560000, 5120000};
        int nLens = lengthSchedule.length;

        // tallies[metric][level][lenIdx][bin] = {correct, total}
        // metric: 0=k5, 1=k4. k5dif=1.0 is a "not computed" flag; use k4 then.
        long[][][][][] tallies = rawMode ? null
                : new long[2][LEVELS.length][nLens][k5difBins + 1][2];

        // Raw mode: write per-hit features directly.
        ByteStreamWriter rawOut = null;
        if (rawMode) {
            String rawFile = outPrefix + "_raw.tsv";
            rawOut = new ByteStreamWriter(rawFile, true, false, false);
            rawOut.start();
            rawOut.print(("#queryName\thitRank\tlength\tk3dif\tk4dif\tk5dif\tcaLevel\n").getBytes());
        }
        String prevQuery = null;
        int hitRank = 0;

        long rowsParsed = 0, rowsSkipped = 0, rowsUsed = 0;
        long[] skipReason = new long[5]; // [0]=badTid, [1]=numFmt, [2]=badLen, [3]=caFail, [4]=shortCols

        for (String f : inputs) {
            outstream.println("Reading " + f + "...");
            ByteFile bf = ByteFile.makeByteFile(f, true);
            byte[] line;
            while ((line = bf.nextLine()) != null) {
                if (line.length == 0 || line[0] == '#') { continue; }
                rowsParsed++;
                // Columns: 0=QueryName, 1=Q_GC, 2=Q_Bases, 3=Q_Contigs, 4=RefName,
                // 5=R_TaxID, 6=R_GC, 7=R_Bases, 8=R_Contigs, 9=R_Level, 10=GCdif,
                // 11=STRdif, 12=k3dif, 13=k4dif, 14=k5dif, 15=lineage (optional)
                String s = new String(line);
                String[] cols = s.split("\t", -1);
                if (cols.length < 15) { rowsSkipped++; skipReason[4]++; continue; }

                int qTid = BinObject.parseTaxID(cols[0]);
                int qBases;
                int rTid;
                float k3dif;
                float k4dif;
                float k5dif;
                try {
                    qBases = Integer.parseInt(cols[2]);
                    rTid = Integer.parseInt(cols[5]);
                    k3dif = Float.parseFloat(cols[12]);
                    k4dif = Float.parseFloat(cols[13]);
                    k5dif = Float.parseFloat(cols[14]);
                } catch (NumberFormatException e) { rowsSkipped++; skipReason[1]++; continue; }

                if (qTid < 1 || rTid < 1) { rowsSkipped++; skipReason[0]++; continue; }

                int lenIdx = lengthBucket(qBases);
                if (lenIdx < 0) { rowsSkipped++; skipReason[2]++; continue; }

                int caLevel = safeCommonAncestorLevel(qTid, rTid);
                if (caLevel < 0) { rowsSkipped++; skipReason[3]++; continue; }

                rowsUsed++;

                if (rawMode) {
                    String queryName = cols[0];
                    if (!queryName.equals(prevQuery)) {
                        prevQuery = queryName;
                        hitRank = 0;
                    }
                    hitRank++;
                    ByteBuilder bb = new ByteBuilder();
                    bb.append(queryName).tab();
                    bb.append(hitRank).tab();
                    bb.append(qBases).tab();
                    bb.append(k3dif, 5).tab();
                    bb.append(k4dif, 5).tab();
                    bb.append(k5dif, 5).tab();
                    bb.append(caLevel).append('\n');
                    rawOut.print(bb.toBytes());
                } else {
                    int k5bin = binIdx(k5dif);
                    int k4bin = binIdx(k4dif);
                    for (int li = 0; li < LEVELS.length; li++) {
                        int L = LEVELS[li];
                        boolean correct = (caLevel > 0 && caLevel <= L);
                        tallies[1][li][lenIdx][k4bin][1]++;
                        if (correct) { tallies[1][li][lenIdx][k4bin][0]++; }
                        if (k5dif < 1.0f) {
                            tallies[0][li][lenIdx][k5bin][1]++;
                            if (correct) { tallies[0][li][lenIdx][k5bin][0]++; }
                        }
                    }
                }
            }
            bf.close();
        }

        outstream.println("Rows parsed: " + rowsParsed + ", used: " + rowsUsed
                + ", skipped: " + rowsSkipped);
        outstream.println("Skip reasons: badTid=" + skipReason[0]
                + " numFmt=" + skipReason[1] + " badLen=" + skipReason[2]
                + " caFail=" + skipReason[3] + " shortCols=" + skipReason[4]);

        if (rawMode) {
            rawOut.poisonAndWait();
            outstream.println("Raw output: " + outPrefix + "_raw.tsv (" + rowsUsed + " rows)");
            t.stop();
            outstream.println("Time: " + t);
            return;
        }

        // Write TSVs: one per (metric, level) pair. 2 * 8 = 16 files.
        String[] metricNames = {"k5", "k4"};
        for (int mi = 0; mi < 2; mi++) {
            for (int li = 0; li < LEVELS.length; li++) {
                String outFile = outPrefix + "_" + metricNames[mi] + "_"
                        + LEVEL_NAMES[li] + ".tsv";
                outstream.println("Writing " + outFile);
                ByteStreamWriter bsw = new ByteStreamWriter(outFile, true, false, false);
                bsw.start();
                ByteBuilder bb = new ByteBuilder();
                bb.append("#len\tbin\tdifLow\tcorrect\ttotal\tprob\n");
                bsw.print(bb.toBytes());
                for (int lenIdx = 0; lenIdx < nLens; lenIdx++) {
                    for (int bin = 0; bin <= k5difBins; bin++) {
                        long correct = tallies[mi][li][lenIdx][bin][0];
                        long total = tallies[mi][li][lenIdx][bin][1];
                        if (total == 0) { continue; }
                        float low = bin * binWidth;
                        float prob = correct / (float) total;
                        bb.clear();
                        bb.append(lengthSchedule[lenIdx]).tab();
                        bb.append(bin).tab();
                        bb.append(low, 5).tab();
                        bb.append(correct).tab();
                        bb.append(total).tab();
                        bb.append(prob, 5).append('\n');
                        bsw.print(bb.toBytes());
                    }
                }
                bsw.poisonAndWait();
            }
        }

        t.stop();
        outstream.println("Time: " + t);
    }

    private int binIdx(float x) {
        int b = (int) Math.floor(x / binWidth);
        if (b < 0) { b = 0; }
        if (b > k5difBins) { b = k5difBins; }
        return b;
    }

    private int lengthBucket(int qBases) {
        if (rawMode) { return 0; } // raw mode doesn't need bucketing
        int best = -1;
        int bestDist = Integer.MAX_VALUE;
        for (int i = 0; i < lengthSchedule.length; i++) {
            int dist = Math.abs(qBases - lengthSchedule[i]);
            if (dist < bestDist) { bestDist = dist; best = i; }
        }
        // Accept if within 50% of the target length
        if (best >= 0 && bestDist <= lengthSchedule[best] * 0.5) { return best; }
        // Also accept anything larger than the last bucket
        if (qBases > lengthSchedule[lengthSchedule.length - 1]) {
            return lengthSchedule.length - 1;
        }
        return -1;
    }

    private int safeCommonAncestorLevel(int a, int b) {
        try {
            return tree.commonAncestorLevel(a, b);
        } catch (Throwable e) {
            return -1;
        }
    }

    // Taxonomic levels we care about. Skip SUBSPECIES(1) and SUPERKINGDOM(9) and LIFE(11).
    // TaxTree constants: SPECIES=2, GENUS=3, FAMILY=4, ORDER=5, CLASS=6,
    // PHYLUM=7, KINGDOM=8, DOMAIN=10.
    private static final int[] LEVELS = {2, 3, 4, 5, 6, 7, 8, 10};
    private static final String[] LEVEL_NAMES = {"species", "genus", "family",
            "order", "class", "phylum", "kingdom", "domain"};

    private final String[] inputs;
    private final String outPrefix;
    private final int k5difBins;
    private final float k5difMax;
    private final float binWidth;
    private int[] lengthSchedule;
    private TaxTree tree;
    private final String treeFile;
    private boolean rawMode;
    private PrintStream outstream = System.err;
}
