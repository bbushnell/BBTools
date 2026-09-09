package prot;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Synthetic fixture for {@link MarkerSalvageDiagnostic} and the {@link MarkerSelector} refactor
 * it depends on (Yoimiya's review, 2026-09-08). No real corpus, no cluster execution, no marker
 * set (re)generated for production use. Covers, per her explicit list: hand-derived P/N and S/P
 * arithmetic, the zero-copy denominator, small-population threshold boundaries, monotonic yield
 * under relaxation, absence of any diagnostic rank-file write, and before/after byte parity of
 * {@link MarkerSelector}'s own existing single-threshold output across the accumulation-extraction
 * refactor.
 *
 * <p>The synthetic universe (fixed, hand-derivable): a 20-organism "TestPhylum" with 6 family
 * ranks at controlled P/S values (rank1 P=20/S=20; rank2 P=19/S=19; rank3 P=15/S=15; rank4
 * P=20/S=14; rank5 P=8/S=6 with 12 organisms never mentioning it at all -- the zero-copy-
 * denominator case; rank6 P=17/S=17); a 5-organism "TinyPhylum" with rank8 at P=4/S=4 (the
 * small-n boundary case: {@code ceil(0.95*5)==ceil(0.90*5)==5}, so P=4 fails identically at both
 * thresholds); and a 2-organism "ExcludedPhylum" (below {@code minOrgs=3}) to prove exclusion
 * happens at output, not accumulation.
 *
 * @author Eru
 */
public final class MarkerSalvageDiagnosticTest {

	/** Golden hashes of {@link MarkerSelector}'s real output on this fixture's synthetic input at
	 *  minprev=0.95/minsc=0.90/minorgs=3, captured from the git-HEAD PRE-REFACTOR source (commit
	 *  ba0dccd..be0ecb8-era `MarkerSelector.java`, before {@code accumulateGroups} existed) via a
	 *  real build+run, then independently reproduced byte-for-byte by the CURRENT (refactored)
	 *  source before these constants were written -- see
	 *  `results/marker_salvage_diagnostic_fixture_20260908.md` for the exact commands. */
	static final String EXPECTED_MARKERSETS_SHA256="e1f4e47d5fdceccd94f22d733011b9cf99ca21e4bc34d830b664507d75216061";
	static final String EXPECTED_TESTPHYLUM_RANKFILE_SHA256="a6e2b7a040683432de03a18fd8a1939a2fdf82585b364bfc874bdd4095c4cae1";
	static final String EXPECTED_TINYPHYLUM_RANKFILE_SHA256="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";

	public static void main(String[] args) throws IOException {
		int pass=0, total=0;

		final File tmpDir=Files.createTempDirectory("msdtest_").toFile();
		tmpDir.deleteOnExit();
		final String dir=tmpDir.getAbsolutePath();

		final String perorg=buildPerorgFixture(dir);
		final String taxpgm=buildTaxpgmFixture(dir);

		// ============================================================
		// Before/after MarkerSelector byte parity across the refactor
		// ============================================================
		final String msOutDir=dir+"/ms_out";
		new File(msOutDir+"/subsets").mkdirs();
		MarkerSelector.main(new String[]{
			"perorg="+perorg, "taxpgm="+taxpgm,
			"out="+msOutDir+"/markersets.tsv", "outdir="+msOutDir+"/subsets",
			"tag=marker", "minprev=0.95", "minsc=0.90", "minorgs=3"
		});
		pass+=check(true, EXPECTED_MARKERSETS_SHA256.equals(sha256HexOfRawFile(new File(msOutDir+"/markersets.tsv"))),
			"MarkerSelector's markersets.tsv on this fixture matches the pre-refactor golden hash (byte parity)"); total++;
		pass+=check(true, EXPECTED_TESTPHYLUM_RANKFILE_SHA256.equals(sha256HexOfRawFile(new File(msOutDir+"/subsets/marker_TestPhylum.txt"))),
			"MarkerSelector's marker_TestPhylum.txt matches the pre-refactor golden hash"); total++;
		pass+=check(true, EXPECTED_TINYPHYLUM_RANKFILE_SHA256.equals(sha256HexOfRawFile(new File(msOutDir+"/subsets/marker_TinyPhylum.txt"))),
			"MarkerSelector's marker_TinyPhylum.txt matches the pre-refactor golden hash"); total++;
		// Real qualification result, independently hand-derived (see class javadoc P/S values):
		// TestPhylum: rank1 (P/N=20/20=1.00>=0.95, S/P=20/20=1.00>=0.90) and rank2
		// (P/N=19/20=0.95>=0.95 exactly, S/P=19/19=1.00>=0.90) qualify; ranks 3-6 do not. 2 markers.
		final ArrayList<String> testPhylumRanks=readLines(new File(msOutDir+"/subsets/marker_TestPhylum.txt"));
		pass+=check(true, testPhylumRanks.size()==2 && testPhylumRanks.contains("1") && testPhylumRanks.contains("2"),
			"MarkerSelector's real qualification at 0.95/0.90 matches hand-derivation: TestPhylum={1,2}"); total++;

		// ============================================================
		// Run the diagnostic on the SAME fixture
		// ============================================================
		final String diagPrefix=dir+"/diag";
		MarkerSalvageDiagnostic.main(new String[]{
			"perorg="+perorg, "taxpgm="+taxpgm, "out="+diagPrefix, "minorgs=3"
		});

		final File rawFile=new File(diagPrefix+"_raw.tsv");
		final File gridFile=new File(diagPrefix+"_grid.tsv");
		pass+=check(true, rawFile.exists(), "diagnostic writes the raw P/S/N table"); total++;
		pass+=check(true, gridFile.exists(), "diagnostic writes the grid-count table"); total++;

		final ArrayList<String[]> raw=readTsvRows(rawFile, 6);
		final ArrayList<String[]> grid=readTsvRows(gridFile, 6);

		// ============================================================
		// Hand-derived P/N and S/P arithmetic (raw table), independent of MarkerSelector's own
		// aggregate n_markers -- reads the raw per-rank numbers directly.
		// ============================================================
		final HashMap<String, int[]> testPhylumRaw=rawByRank(raw, "TestPhylum");
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("1"), new int[]{20, 20}), "raw: TestPhylum rank1 P=20,S=20"); total++;
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("2"), new int[]{19, 19}), "raw: TestPhylum rank2 P=19,S=19"); total++;
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("3"), new int[]{15, 15}), "raw: TestPhylum rank3 P=15,S=15"); total++;
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("4"), new int[]{20, 14}), "raw: TestPhylum rank4 P=20,S=14"); total++;
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("6"), new int[]{17, 17}), "raw: TestPhylum rank6 P=17,S=17"); total++;

		// ============================================================
		// Zero-copy denominator: rank5 present in only 8/20 organisms (12 never mention it at
		// all) -- P must be 8 and the GROUP's n_orgs must still be the full 20, not 8. A bug that
		// used "organisms mentioning any family" as the denominator would not show up here since
		// n_orgs lives on the GROUP not the rank -- checked directly: the raw row's n_orgs column
		// (index 2) is 20, matching the group total, independent of rank5's own P=8.
		// ============================================================
		pass+=check(true, java.util.Arrays.equals(testPhylumRaw.get("5"), new int[]{8, 6}), "raw: TestPhylum rank5 P=8,S=6 (zero-copy denominator case)"); total++;
		final int testPhylumNOrgs=findNOrgs(raw, "TestPhylum");
		pass+=check(true, testPhylumNOrgs==20, "raw: TestPhylum group n_orgs is the full 20 (population denominator unaffected by rank5's absence in 12 organisms)"); total++;

		// ============================================================
		// Small-n boundary: TinyPhylum rank8, P=4 of N=5. ceil(0.95*5)=5 and ceil(0.90*5)=5 are
		// IDENTICAL thresholds at this n -- P=4 fails both, a real discreteness artifact (matches
		// the real Asgardarchaeota finding in the source audit), not a bug.
		// ============================================================
		final HashMap<String, int[]> tinyPhylumRaw=rawByRank(raw, "TinyPhylum");
		pass+=check(true, java.util.Arrays.equals(tinyPhylumRaw.get("8"), new int[]{4, 4}), "raw: TinyPhylum rank8 P=4,S=4"); total++;
		final int tinyAt095=gridCount(grid, "TinyPhylum", 0.95, 0.90);
		final int tinyAt090=gridCount(grid, "TinyPhylum", 0.90, 0.90);
		final int tinyAt085=gridCount(grid, "TinyPhylum", 0.85, 0.90);
		final int tinyAt080=gridCount(grid, "TinyPhylum", 0.80, 0.90);
		pass+=check(true, tinyAt095==0 && tinyAt090==0 && tinyAt085==0,
			"small-n boundary: TinyPhylum n_markers IDENTICAL (0) at minPrev=0.95, 0.90, AND 0.85 (P=4 fails all three -- ceil(x*5) doesn't distinguish them)"); total++;
		pass+=check(true, tinyAt080==1,
			"small-n boundary: TinyPhylum n_markers becomes 1 exactly at minPrev=0.80 (P=4>=ceil(0.80*5)=4)"); total++;

		// ============================================================
		// Cross-check: diagnostic's grid at the two REAL data points must equal MarkerSelector's
		// own actual n_markers count -- proves the diagnostic's counting logic matches the real
		// tool, not just an independent guess.
		// ============================================================
		pass+=check(true, gridCount(grid, "TestPhylum", 0.95, 0.90)==2,
			"diagnostic grid at (0.95,0.90) matches MarkerSelector's real n_markers=2 for TestPhylum"); total++;

		// ============================================================
		// group kind classification: ALL_domain vs phylum vs excluded_minorgs, and confirmation
		// that ExcludedPhylum (n_orgs=2<minOrgs=3) IS present in both tables (accumulated, then
		// only filtered at MarkerSelector's OWN output stage -- not "never formed").
		// ============================================================
		pass+=check(true, findKind(raw, "ALL_bacteria").equals("all_domain"), "kind: ALL_bacteria classified all_domain"); total++;
		pass+=check(true, findKind(raw, "TestPhylum").equals("phylum"), "kind: TestPhylum classified phylum"); total++;
		pass+=check(true, findKind(raw, "ExcludedPhylum").equals("excluded_minorgs"), "kind: ExcludedPhylum classified excluded_minorgs"); total++;
		pass+=check(true, findNOrgs(raw, "ExcludedPhylum")==2, "ExcludedPhylum IS present with n_orgs=2 in the raw table (accumulated, not skipped)"); total++;
		pass+=check(true, findKind(grid, "ExcludedPhylum").equals("excluded_minorgs"), "ExcludedPhylum IS present in the grid table too, classified excluded_minorgs"); total++;

		// ============================================================
		// Monotonic yields: relaxing (lowering) minPrev at fixed minSc, or lowering minSc at
		// fixed minPrev, must never DECREASE n_markers for any group -- checked across the whole
		// grid, all 4 groups.
		// ============================================================
		pass+=check(true, monotonicAcrossGrid(grid), "n_markers is monotonically non-decreasing as either threshold relaxes, for every group in the grid"); total++;

		// ============================================================
		// NaN/nonfinite rejection (Yoimiya/Root's review, 2026-09-08): Double.parseDouble("NaN")
		// succeeds, and every comparison with NaN (including <0 and >1) is false, so a bare range
		// check silently admits it. parseGrid must reject NaN, +Infinity, -Infinity explicitly.
		// ============================================================
		pass+=check(true, parseGridThrows("0.95,NaN,0.70"), "parseGrid throws on a NaN grid value (silently passed the old <0||>1 check)"); total++;
		pass+=check(true, parseGridThrows("0.95,Infinity,0.70"), "parseGrid throws on +Infinity"); total++;
		pass+=check(true, parseGridThrows("0.95,-Infinity,0.70"), "parseGrid throws on -Infinity"); total++;
		pass+=check(false, parseGridThrows("0.95,0.90,0.70"), "parseGrid still accepts an ordinary finite grid"); total++;

		// ============================================================
		// Custom-precision threshold: label text must match EXACTLY what countQualifying used,
		// not a rounded display (Yoimiya/Root's review, 2026-09-08). Direct unit-level proof first
		// (bypasses file I/O): a hand-built group with N=1000, P=949 (P/N=0.949 exactly) qualifies
		// at the TRUE threshold 0.949 (949>=949) but would NOT qualify if the code silently used a
		// rounded "0.95" (949<950). This is a real boundary, not a contrived one.
		// ============================================================
		final MarkerSelector.Group precisionGroup=new MarkerSelector.Group();
		precisionGroup.nOrgs=1000;
		precisionGroup.present.put(20, 949);
		precisionGroup.single.put(20, 949);
		pass+=check(1, MarkerSalvageDiagnostic.countQualifying(precisionGroup, 0.949, 0.0),
			"countQualifying at the EXACT unrounded threshold 0.949: P=949>=949 qualifies (count=1)"); total++;
		pass+=check(0, MarkerSalvageDiagnostic.countQualifying(precisionGroup, 0.95, 0.0),
			"countQualifying at the (different, rounded-looking) threshold 0.95: P=949<950 does NOT qualify (count=0) -- proves 0.949 and 0.95 are genuinely different thresholds, not the same value rounded for display"); total++;

		// End-to-end: run the diagnostic with a custom prevgrid=0.949 and confirm the WRITTEN
		// table's label text round-trips to exactly 0.949 (not "0.95"), and its n_markers count
		// for PrecisionPhylum matches the true-threshold unit result above (1, not 0).
		final String precisionPerorg=writeLines(dir, "perorg_precision.tsv",
			buildPrecisionOrgLines());
		final String precisionTaxpgm=writeLines(dir, "taxpgm_precision.tsv",
			buildPrecisionTaxpgmLines());
		final String precisionPrefix=dir+"/precision_diag";
		MarkerSalvageDiagnostic.main(new String[]{
			"perorg="+precisionPerorg, "taxpgm="+precisionTaxpgm, "out="+precisionPrefix,
			"minorgs=3", "prevgrid=0.949", "scgrid=0.00"
		});
		final ArrayList<String[]> precisionGrid=readTsvRows(new File(precisionPrefix+"_grid.tsv"), 6);
		boolean foundExactLabel=false, labelWasWrong=false;
		int precisionNMarkers=-1;
		for(String[] row : precisionGrid){
			if(row[0].equals("PrecisionPhylum")){
				if(row[3].equals("0.949")){foundExactLabel=true; precisionNMarkers=Integer.parseInt(row[5]);}
				if(row[3].equals("0.95")){labelWasWrong=true;}
			}
		}
		pass+=check(true, foundExactLabel && !labelWasWrong,
			"end-to-end grid table's threshold LABEL for a custom prevgrid=0.949 reads back as exactly \"0.949\", never rounded to \"0.95\""); total++;
		pass+=check(true, precisionNMarkers==1,
			"end-to-end grid table's n_markers at the 0.949 row is 1 (matches the true-threshold unit result, not the rounded-threshold result of 0)"); total++;

		// ============================================================
		// No rank-file writes anywhere under the diagnostic's output prefix.
		// ============================================================
		final File[] diagFiles=tmpDir.listFiles();
		boolean noRankFile=true;
		for(File f : diagFiles==null ? new File[0] : diagFiles){
			if(f.getName().startsWith("diag") && !f.getName().endsWith("_raw.tsv") && !f.getName().endsWith("_grid.tsv")){
				noRankFile=false;
			}
		}
		final File diagSubsetsDir=new File(dir, "subsets");
		pass+=check(true, noRankFile && !diagSubsetsDir.exists(),
			"no rank-file (or subsets/ dir) was written anywhere by the diagnostic -- only _raw.tsv and _grid.tsv exist"); total++;

		System.err.println("MarkerSalvageDiagnosticTest: "+pass+"/"+total+" passed.");
		if(pass!=total){throw new RuntimeException((total-pass)+" regression case(s) failed.");}
	}

	// ---------------- fixture construction ----------------

	/** 27 organisms, 16-field sparse rows (13 placeholder fields at 2..14, famcounts at 15),
	 *  matching this file's javadoc P/S design exactly. */
	static String buildPerorgFixture(String dir) throws IOException {
		final ArrayList<String> lines=new ArrayList<String>();
		for(int i=1; i<=20; i++){
			final int tid=1000+i;
			final StringBuilder fam=new StringBuilder("1:1");
			if(i<=19){fam.append(";2:1");}
			if(i<=15){fam.append(";3:1");}
			fam.append(i<=14 ? ";4:1" : ";4:2");
			if(i<=8){fam.append(i<=6 ? ";5:1" : ";5:2");}
			if(i<=17){fam.append(";6:1");}
			lines.add(perorgLine(tid, "bacteria", fam.toString()));
		}
		for(int i=1; i<=5; i++){
			final int tid=2000+i;
			lines.add(perorgLine(tid, "bacteria", i<=4 ? "8:1" : ""));
		}
		lines.add(perorgLine(3001, "bacteria", "9:1"));
		lines.add(perorgLine(3002, "bacteria", ""));
		return writeLines(dir, "perorg_sparse_test.tsv", lines.toArray(new String[0]));
	}

	static String perorgLine(int tid, String domain, String famcounts){
		final StringBuilder sb=new StringBuilder();
		sb.append(tid).append('\t').append(domain);
		for(int i=0; i<13; i++){sb.append("\t0");}
		sb.append('\t').append(famcounts);
		return sb.toString();
	}

	/** 1000 organisms, "PrecisionPhylum": rank20 present (single-copy) in exactly 949 of them
	 *  (P=949, N=1000 -> P/N=0.949 exactly), absent in the remaining 51. */
	static String[] buildPrecisionOrgLines(){
		final String[] lines=new String[1000];
		for(int i=1; i<=1000; i++){
			lines[i-1]=perorgLine(4000+i, "bacteria", i<=949 ? "20:1" : "");
		}
		return lines;
	}

	static String[] buildPrecisionTaxpgmLines(){
		final String[] lines=new String[1000];
		for(int i=1; i<=1000; i++){lines[i-1]=(4000+i)+"\tPrecisionPhylum";}
		return lines;
	}

	static boolean parseGridThrows(String s){
		try{MarkerSalvageDiagnostic.parseGrid(s); return false;}
		catch(RuntimeException e){return true;}
	}

	static String buildTaxpgmFixture(String dir) throws IOException {
		final ArrayList<String> lines=new ArrayList<String>();
		for(int i=1; i<=20; i++){lines.add((1000+i)+"\tTestPhylum");}
		for(int i=1; i<=5; i++){lines.add((2000+i)+"\tTinyPhylum");}
		lines.add("3001\tExcludedPhylum");
		lines.add("3002\tExcludedPhylum");
		return writeLines(dir, "taxpgm_test.tsv", lines.toArray(new String[0]));
	}

	static String writeLines(String dir, String name, String... lines) throws IOException {
		final File f=new File(dir, name);
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			for(String line : lines){w.write(line); w.write('\n');}
		}
		return f.getAbsolutePath();
	}

	// ---------------- TSV / assertion helpers ----------------

	static ArrayList<String> readLines(File f) throws IOException {
		final ArrayList<String> out=new ArrayList<String>();
		for(String line : Files.readAllLines(f.toPath())){
			if(!line.isEmpty()){out.add(line);}
		}
		return out;
	}

	/** Reads all non-comment rows with exactly {@code expectedFields} tab-separated fields. */
	static ArrayList<String[]> readTsvRows(File f, int expectedFields) throws IOException {
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		for(String line : Files.readAllLines(f.toPath())){
			if(line.isEmpty() || line.startsWith("#")){continue;}
			final String[] parts=line.split("\t", -1);
			if(parts.length!=expectedFields){throw new RuntimeException("Expected "+expectedFields+" fields, got "+parts.length+": "+line);}
			rows.add(parts);
		}
		return rows;
	}

	/** Raw table columns: group, kind, n_orgs, rank, P, S. Returns rank->{P,S} for one group. */
	static HashMap<String, int[]> rawByRank(ArrayList<String[]> raw, String group){
		final HashMap<String, int[]> out=new HashMap<String, int[]>();
		for(String[] row : raw){
			if(row[0].equals(group)){out.put(row[3], new int[]{Integer.parseInt(row[4]), Integer.parseInt(row[5])});}
		}
		return out;
	}

	static int findNOrgs(ArrayList<String[]> rows, String group){
		for(String[] row : rows){if(row[0].equals(group)){return Integer.parseInt(row[2]);}}
		throw new RuntimeException("Group not found: "+group);
	}

	static String findKind(ArrayList<String[]> rows, String group){
		for(String[] row : rows){if(row[0].equals(group)){return row[1];}}
		throw new RuntimeException("Group not found: "+group);
	}

	/** Grid table columns: group, kind, n_orgs, minPrev, minSc, n_markers. */
	static int gridCount(ArrayList<String[]> grid, String group, double minPrev, double minSc){
		for(String[] row : grid){
			if(row[0].equals(group) && closeTo(Double.parseDouble(row[3]), minPrev) && closeTo(Double.parseDouble(row[4]), minSc)){
				return Integer.parseInt(row[5]);
			}
		}
		throw new RuntimeException("Grid cell not found: "+group+" "+minPrev+" "+minSc);
	}

	static boolean closeTo(double a, double b){return Math.abs(a-b)<1e-9;}

	/** For every group, for every fixed minSc, n_markers must be non-decreasing as minPrev
	 *  decreases (sorted descending->ascending order of grid values); symmetric check for minSc. */
	static boolean monotonicAcrossGrid(ArrayList<String[]> grid){
		final HashMap<String, ArrayList<double[]>> byGroup=new HashMap<String, ArrayList<double[]>>();
		for(String[] row : grid){
			final String g=row[0];
			ArrayList<double[]> list=byGroup.get(g);
			if(list==null){list=new ArrayList<double[]>(); byGroup.put(g, list);}
			list.add(new double[]{Double.parseDouble(row[3]), Double.parseDouble(row[4]), Integer.parseInt(row[5])});
		}
		for(ArrayList<double[]> cells : byGroup.values()){
			//Fix minSc, sweep minPrev descending (0.95->0.70): as minPrev decreases, n_markers must not decrease.
			final ArrayList<Double> scValues=new ArrayList<Double>();
			for(double[] c : cells){if(!scValues.contains(c[1])){scValues.add(c[1]);}}
			for(double sc : scValues){
				final ArrayList<double[]> row=new ArrayList<double[]>();
				for(double[] c : cells){if(closeTo(c[1], sc)){row.add(c);}}
				row.sort((x, y) -> Double.compare(y[0], x[0]));//descending minPrev
				for(int i=1; i<row.size(); i++){
					if(row.get(i)[2]<row.get(i-1)[2]){return false;}
				}
			}
			//Fix minPrev, sweep minSc descending: same requirement.
			final ArrayList<Double> prevValues=new ArrayList<Double>();
			for(double[] c : cells){if(!prevValues.contains(c[0])){prevValues.add(c[0]);}}
			for(double pv : prevValues){
				final ArrayList<double[]> row=new ArrayList<double[]>();
				for(double[] c : cells){if(closeTo(c[0], pv)){row.add(c);}}
				row.sort((x, y) -> Double.compare(y[1], x[1]));//descending minSc
				for(int i=1; i<row.size(); i++){
					if(row.get(i)[2]<row.get(i-1)[2]){return false;}
				}
			}
		}
		return true;
	}

	static String sha256HexOfRawFile(File f) throws IOException {
		final MessageDigest md;
		try{md=MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
		try(InputStream in=new FileInputStream(f)){
			final byte[] buf=new byte[65536];
			int n;
			while((n=in.read(buf))>=0){md.update(buf, 0, n);}
		}
		final byte[] digest=md.digest();
		final StringBuilder sb=new StringBuilder(digest.length*2);
		for(byte b : digest){sb.append(String.format("%02x", b));}
		return sb.toString();
	}

	static int check(boolean expected, boolean actual, String label){
		final boolean ok=(expected==actual);
		System.err.println((ok?"PASS":"FAIL")+": "+label+" (expected="+expected+" actual="+actual+")");
		return ok ? 1 : 0;
	}

	static int check(int expected, int actual, String label){
		final boolean ok=(expected==actual);
		System.err.println((ok?"PASS":"FAIL")+": "+label+" (expected="+expected+" actual="+actual+")");
		return ok ? 1 : 0;
	}
}
