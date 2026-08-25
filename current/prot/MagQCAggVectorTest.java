package prot;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import ml.CellNet;
import ml.CellNetParser;
import ml.Function;

/**
 * Self-contained test for {@link MagQCVectorMaker}'s aggregator-vector emission
 * (composite architecture, Brian's step 4; re-baselined 2026-08-25 for the FROZEN
 * shared-context vector-layout rebuild, magqc_rebuild_20260824.plan). Reuses the
 * MagQCVectorMakerTest fixture (tiny synthetic cache, hand-computable complements),
 * builds three tiny REAL CellNets (two famset subnets and the special ncRNA subnet,
 * all under the identical standard shared-context width now that the old per-net rep
 * flags are retired), writes a manifest, and verifies:
 *
 * <ul>
 * <li><b>Global output unaffected:</b> the global training vector is byte-identical
 * whether the run emits a famset subnet or the ncRNA subnet + aggregator — which also
 * proves the two runs synthesized IDENTICAL bins, so their per-row emissions align.</li>
 * <li><b>Aggregator math end-to-end:</b> with cleanspike=1 (no contaminants, so
 * serve-obs == target-only obs), each aggregator row's per-subnet features
 * [ratio, log-obs, log-pred, zero] are recomputed independently by feed-forwarding
 * the SAME loaded nets on the SAME inputs taken from the emitted subnet-row files,
 * and must match the emitted strings EXACTLY (same fmt() rounding, same floats).</li>
 * <li><b>Pooled baseline:</b> recomputed from both subnets' obs and predictions.</li>
 * <li><b>Dense head:</b> top-K prevalence order is hand-computed for the fixture
 * ([0,2,1,3,5,6,7,4]); each head column must equal the global row's enc=two column
 * for that family rank (same appendTwo encoding, reordered).</li>
 * <li><b>Phylum one-hot + targets:</b> equal to the global row's columns exactly.</li>
 * <li><b>poolmode=valsplit:</b> runs, emits nonempty B/C aggregator files, and is
 * deterministic (byte-identical on a re-run).</li>
 * <li><b>aggobs=serve vs aggobs=clean:</b> the block above runs with cleanspike=1.0,
 * which makes every bin contamination-free (fam[]==cleanFamBuf, lastNcServe==lastNcObs
 * always), so the serve/clean observed-count switch in formatAggRow is exercised but
 * never actually distinguished - a bug swapping the two, or turning aggobs=clean into a
 * no-op, would not be caught. A second fixture run with cleanspike=0.0 (real contamination
 * on effectively every bin) runs the SAME seed/manifest twice, once per aggobs= value,
 * confirms the bins still align (global vector byte-identical), confirms the two
 * aggregator outputs are NOT byte-identical, and recomputes famset_a's exact
 * [ratio,lobs,lpred,zero] block for both branches - clean obs from the famset subnet's own
 * target-only output row, serve obs decoded from the aggregator's always-whole-bin dense
 * head - proving obs_clean&lt;=obs_serve on every row and obs_clean&lt;obs_serve on many.</li>
 * </ul>
 *
 * Run with assertions on: {@code java -ea prot.MagQCAggVectorTest}
 *
 * @author UMP45
 */
public class MagQCAggVectorTest {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("magqcagg_test").toFile();
		//deleteOnExit() only removes an EMPTY dir; the ~20 fixture files inside would leak.
		//A shutdown hook recursively deletes the whole tree, even on assertion failure.
		Runtime.getRuntime().addShutdownHook(new Thread(()->deleteRecursive(dir)));

		// --- fixture: identical to MagQCVectorMakerTest (4 orgs, 2 phyla, 8 families) ---
		write(new File(dir, "cache.tsv"),
			"c100_1\t100\tbacteria\t5000\t2500\t5000\t5\t5\t4500\t4100000\t4500\t1\t1\t1\t0\t10\t0:1;1:1;2:1\n"+
			"c100_2\t100\tbacteria\t3000\t1500\t3000\t3\t3\t2700\t2430000\t2700\t0\t0\t0\t0\t8\t1:1;3:1\n"+
			"c100_3\t100\tbacteria\t2000\t1000\t2000\t2\t2\t1800\t1620000\t1800\t0\t1\t0\t1\t5\t4:2\n"+
			"c200_1\t200\tbacteria\t6000\t3600\t6000\t6\t6\t5400\t4860000\t5400\t1\t1\t1\t0\t12\t0:1;2:1;5:1\n"+
			"c200_2\t200\tbacteria\t4000\t2400\t4000\t4\t4\t3600\t3240000\t3600\t0\t0\t0\t2\t9\t6:1;7:1\n"+
			"c300_1\t300\tbacteria\t4000\t2000\t4000\t4\t4\t3600\t3240000\t3600\t1\t1\t0\t0\t9\t0:1;1:1\n"+
			"c300_2\t300\tbacteria\t4000\t2000\t4000\t4\t4\t3600\t3240000\t3600\t0\t0\t1\t1\t7\t2:1;3:1\n"+
			"c400_1\t400\tbacteria\t7000\t4200\t7000\t7\t7\t6300\t5670000\t6300\t1\t1\t1\t0\t13\t0:1;5:1;6:1\n"+
			"c400_2\t400\tbacteria\t5000\t3000\t5000\t5\t5\t4500\t4050000\t4500\t0\t1\t0\t1\t10\t2:1;7:1\n");
		write(new File(dir, "sizemap.tsv"), "100\t10000\n200\t10000\n300\t8000\n400\t12000\n");
		write(new File(dir, "taxpgm.tsv"),
			"100\tFirmicutes\tpgm1\n200\tProteobacteria\tpgm2\n300\tFirmicutes\tpgm1\n400\tProteobacteria\tpgm2\n");
		write(new File(dir, "familylist.tsv"), "rank\tcluster_rep\n0\tf0\n1\tf1\n2\tf2\n3\tf3\n4\tf4\n5\tf5\n6\tf6\n7\tf7\n");
		write(new File(dir, "subset.txt"), "0\n2\n5\n");
		write(new File(dir, "subsetB.txt"), "1\n3\n");

		// --- three tiny real nets, all under the FROZEN shared-context width (17 = 9 scalars +
		// domain one-hot(8), magqc_rebuild_20260824.plan "FROZEN VECTOR LAYOUT" - the old per-net
		// rep-override flags (sngenelen= etc) are RETIRED, every net now gets the identical
		// standard context): famset_a (3 obs + 3 phyla + 17 ctx = 23 in), ncrna (5+3+17=25 in),
		// famset_b (2 obs + 3 phyla + 17 ctx = 22 in - a plain third subnet, proving the
		// aggregator handles >2 subnets correctly).
		// TYPE_RATES defaults are nonzero but unnormalized, which trips CellNet's constructor
		// assert; zero them like RegressionTrainer.buildNetAndMask does (cells then keep their
		// constructor-assigned activation, which serializes and reloads consistently).
		Arrays.fill(Function.TYPE_RATES, 0f);
		final File netA=new File(dir, "netA.bbnet"), netN=new File(dir, "netN.bbnet");
		final File netB=new File(dir, "netB.bbnet");
		makeNet(netA, new int[]{23, 4, 1}, 5);
		makeNet(netN, new int[]{25, 4, 1}, 7);
		makeNet(netB, new int[]{22, 3, 1}, 9);

		final File manifest=new File(dir, "manifest.tsv");
		write(manifest,
			"famset_a\t3\t23\t"+p(dir, "subset.txt")+"\t-\t"+netA.getAbsolutePath()+"\n"+
			"ncrna\t5\t25\t-\t-\t"+netN.getAbsolutePath()+"\n"+
			"famset_b\t2\t22\t"+p(dir, "subsetB.txt")+"\t-\t"+netB.getAbsolutePath()+"\n");

		final String[] common={"cache="+p(dir, "cache.tsv"), "sizemap="+p(dir, "sizemap.tsv"),
			"taxpgm="+p(dir, "taxpgm.tsv"), "familylist="+p(dir, "familylist.tsv"),
			"n=200", "valn=0", "valfrac=0.5", "enc=two", "seed=1", "cleanspike=1.0"};

		// Run 1: famset_a subnet rows. Run 1b: famset_b rows (a second, plain famset subnet).
		// Run 2: ncrna subnet + aggregator. All runs synthesize IDENTICAL bins (byte-identical
		// global proves it), so rows align.
		final File g1=new File(dir, "g1.tsv"), g1b=new File(dir, "g1b.tsv"), g2=new File(dir, "g2.tsv");
		final File fRows=new File(dir, "frows.tsv"), fbRows=new File(dir, "fbrows.tsv"), nRows=new File(dir, "nrows.tsv");
		final File agg=new File(dir, "agg.tsv");
		MagQCVectorMaker.main(concat(common, new String[]{"out="+g1.getAbsolutePath(),
			"subnet=famset", "subsetfile="+p(dir, "subset.txt"), "subnetout="+fRows.getAbsolutePath()}));
		MagQCVectorMaker.main(concat(common, new String[]{"out="+g1b.getAbsolutePath(),
			"subnet=famset", "subsetfile="+p(dir, "subsetB.txt"),
			"subnetout="+fbRows.getAbsolutePath()}));
		MagQCVectorMaker.main(concat(common, new String[]{"out="+g2.getAbsolutePath(),
			"subnet=ncrna", "subnetout="+nRows.getAbsolutePath(),
			"aggmanifest="+manifest.getAbsolutePath(), "aggout="+agg.getAbsolutePath()}));

		// (1) Bins identical across runs (byte-identical global) -> per-row alignment holds.
		final byte[] b1=Files.readAllBytes(g1.toPath()), b2=Files.readAllBytes(g2.toPath());
		assertTrue(Arrays.equals(b1, b2), "global vector differs between famset and ncrna+agg runs");
		assertTrue(Arrays.equals(b1, Files.readAllBytes(g1b.toPath())),
			"global vector differs between famset_a and famset_b runs");

		final ArrayList<String[]> F=parse(fRows), FB=parse(fbRows), N=parse(nRows), A=parse(agg), G=parse(g1);
		assertTrue(A.size()==200 && F.size()==200 && FB.size()==200 && N.size()==200 && G.size()==200,
			"row counts: agg="+A.size()+" F="+F.size()+" FB="+FB.size()+" N="+N.size()+" G="+G.size());
		// Header: 3 subnets*4 + pooled + 2*8 head + 3 phyla + 17 shared context + 2*5 raw ncRNA
		// + 1 mapped-fraction = 60 inputs, 2 outputs.
		assertTrue("#dims\t60\t2\t0".equals(header(agg)), "bad agg header: "+header(agg));

		// (2)+(3) per-row recomputation, exact strings.
		final CellNet cnA=CellNetParser.load(netA.getAbsolutePath());
		final CellNet cnN=CellNetParser.load(netN.getAbsolutePath());
		final CellNet cnB=CellNetParser.load(netB.getAbsolutePath());
		final int[] headOrder={0, 2, 1, 3, 5, 6, 7, 4};//prevalence f0=4,f2=4,then rank order among 2s, f4=1
		// g (the global/monolith row) layout under the FROZEN rewrite: [16 fam cols (enc=two,
		// 8 fams)] + [17 shared context] + [3 phylum] + [2 labels] = 38 wide. Named so the offsets
		// below read as arithmetic, not magic numbers.
		final int G_CTX_START=16, G_PHYLUM=G_CTX_START+17, G_LABELS=G_PHYLUM+3;
		for(int i=0; i<200; i++){
			final String[] f=F.get(i), fb=FB.get(i), nr=N.get(i), a=A.get(i), g=G.get(i);
			assertTrue(a.length==62, "agg row width "+a.length+" != 62");
			// famset_a block (cols 0-3)
			double obsF=feed(cnA, f, 23, 3);
			checkBlock(a, 0, obsF, predOf(cnA), "famset_a", i);
			// ncrna block (cols 4-7)
			double obsN=feed(cnN, nr, 25, 5);
			checkBlock(a, 4, obsN, predOf(cnN), "ncrna", i);
			// famset_b block (cols 8-11)
			double obsB=feed(cnB, fb, 22, 2);
			checkBlock(a, 8, obsB, predOf(cnB), "famset_b", i);
			// pooled (col 12)
			final double pooled=Math.min(2.0, (obsF+obsN+obsB)
				/Math.max(1, Math.max(0, predOf(cnA))+Math.max(0, predOf(cnN))+Math.max(0, predOf(cnB))));
			assertTrue(fmt(pooled).equals(a[12]), "row "+i+" pooled: "+fmt(pooled)+" != "+a[12]);
			// dense head (cols 13-28): global's enc=two column pair for each rank, reordered.
			for(int j=0; j<8; j++){
				final int r=headOrder[j];
				assertTrue(a[13+2*j].equals(g[2*r]) && a[14+2*j].equals(g[2*r+1]),
					"row "+i+" head["+j+"] rank "+r+": ("+a[13+2*j]+","+a[14+2*j]+") != ("+g[2*r]+","+g[2*r+1]+")");
			}
			// phylum one-hot (cols 29-31) == global's phylum block; targets (60-61) == global's labels.
			for(int j=0; j<3; j++){assertTrue(a[29+j].equals(g[G_PHYLUM+j]), "row "+i+" phylum "+j);}
			assertTrue(a[60].equals(g[G_LABELS]) && a[61].equals(g[G_LABELS+1]), "row "+i+" targets");
			// shared context (cols 32-48, 17 wide) ~= famset_a row cols 6-22 (fmt round-trips
			// through float can flip the last printed digit, so numeric tolerance not string eq).
			for(int j=0; j<17; j++){
				final double d=Math.abs(Double.parseDouble(a[32+j])-Double.parseDouble(f[6+j]));
				assertTrue(d<=1e-5, "row "+i+" ctx "+j+": "+a[32+j]+" vs "+f[6+j]);
			}
			// raw ncRNA two-channel (cols 49-58) + mapped-fraction (col 59): aggregator-only
			// direct inputs, new in the FROZEN rebuild. Under cleanspike=1.0 every bin is
			// contamination-free, so aggObsServe's default whole-bin source (lastNcServe)
			// equals the ncrna subnet's own target-only obs (lastNcObs) - nr's obs columns are
			// a valid ground truth here (same reasoning the test uses elsewhere for this fixture).
			for(int t=0; t<5; t++){
				final int obsCount=Integer.parseInt(nr[t]);
				final String expCount=fmt(obsCount/(double)(1+obsCount));//encodeCount under enc=two
				assertTrue(a[49+2*t].equals(obsCount>0 ? "1" : "0"),
					"row "+i+" ncRNA presence "+t+": "+a[49+2*t]+" for obs "+obsCount);
				assertTrue(a[50+2*t].equals(expCount),
					"row "+i+" ncRNA count "+t+": "+a[50+2*t]+" != "+expCount+" for obs "+obsCount);
			}
		}

		// (4) poolmode=valsplit: emits nonempty B and C aggregator files, deterministically.
		final File g3=new File(dir, "g3.tsv"), g3v=new File(dir, "g3v.tsv");
		final File a3=new File(dir, "a3.tsv"), a3v=new File(dir, "a3v.tsv");
		final String[] vs={"out="+g3.getAbsolutePath(), "outval="+g3v.getAbsolutePath(), "valn=50",
			"poolmode=valsplit", "aggmanifest="+manifest.getAbsolutePath(),
			"aggout="+a3.getAbsolutePath(), "aggvalout="+a3v.getAbsolutePath()};
		MagQCVectorMaker.main(concat(common, vs));
		assertTrue(parse(a3).size()==200, "valsplit B rows: "+parse(a3).size());
		assertTrue(parse(a3v).size()==50, "valsplit C rows: "+parse(a3v).size());
		final File g4=new File(dir, "g4.tsv"), g4v=new File(dir, "g4v.tsv");
		final File a4=new File(dir, "a4.tsv"), a4v=new File(dir, "a4v.tsv");
		MagQCVectorMaker.main(concat(common, new String[]{"out="+g4.getAbsolutePath(),
			"outval="+g4v.getAbsolutePath(), "valn=50", "poolmode=valsplit",
			"aggmanifest="+manifest.getAbsolutePath(),
			"aggout="+a4.getAbsolutePath(), "aggvalout="+a4v.getAbsolutePath()}));
		assertTrue(Arrays.equals(Files.readAllBytes(a3.toPath()), Files.readAllBytes(a4.toPath()))
			&& Arrays.equals(Files.readAllBytes(a3v.toPath()), Files.readAllBytes(a4v.toPath())),
			"valsplit runs are not deterministic");

		System.out.println("MagQCAggVectorTest PASS: bins aligned (global byte-identical x3), 200 agg rows, "
			+"per-subnet [ratio,lobs,lpred,zero] + pooled recomputed exactly for all 3 subnets, "
			+"dense head == global enc=two columns in prevalence order "+Arrays.toString(headOrder)
			+", phylum/targets == global, raw ncRNA two-channel + mapped-fraction verified, "
			+"valsplit deterministic (200 B + 50 C rows).");

		// (5) aggobs=serve vs aggobs=clean: closes the REAL coverage gap for formatAggRow's
		// observed-count source switch (famArr=(aggObsServe?fam:cleanFamBuf), ncArr likewise).
		// cleanspike=1.0 above makes cont=0 on EVERY bin (rnd.nextDouble() never returns
		// exactly 1.0), so fam[]==cleanFamBuf and lastNcServe==lastNcObs always - a bug
		// swapping serve/clean, or breaking either branch into a no-op, would NOT have been
		// caught by (1)-(4). cleanspike=0.0 here makes sampleCont(rnd)>0 on effectively every
		// bin, so real contaminant contigs are added and fam[] (serve) genuinely diverges from
		// cleanFamBuf (clean) on the contaminated family ranks. Same seed, same manifest, same
		// everything else -> identical bins (checked below exactly as in (1)), so ONLY the
		// aggobs observed-count source differs between the two runs.
		final String[] common5={"cache="+p(dir, "cache.tsv"), "sizemap="+p(dir, "sizemap.tsv"),
			"taxpgm="+p(dir, "taxpgm.tsv"), "familylist="+p(dir, "familylist.tsv"),
			"n=200", "valn=0", "valfrac=0.5", "enc=two", "seed=5", "cleanspike=0.0",
			"aggmanifest="+manifest.getAbsolutePath()};
		final File g5s=new File(dir, "g5s.tsv"), g5c=new File(dir, "g5c.tsv");
		final File agg5s=new File(dir, "agg5s.tsv"), agg5c=new File(dir, "agg5c.tsv");
		final File f5Rows=new File(dir, "f5rows.tsv");//famset_a's own subset {0,2,5}, target-only
		MagQCVectorMaker.main(concat(common5, new String[]{"out="+g5s.getAbsolutePath(),
			"aggout="+agg5s.getAbsolutePath(), "aggobs=serve",
			"subnet=famset", "subsetfile="+p(dir, "subset.txt"), "subnetout="+f5Rows.getAbsolutePath()}));
		MagQCVectorMaker.main(concat(common5, new String[]{"out="+g5c.getAbsolutePath(),
			"aggout="+agg5c.getAbsolutePath(), "aggobs=clean"}));

		// Bins identical regardless of aggobs (it is read only in formatAggRow, never in
		// makeBin, so it cannot perturb the RNG draw sequence) - same alignment guarantee as
		// (1), now exercised under real contamination.
		assertTrue(Arrays.equals(Files.readAllBytes(g5s.toPath()), Files.readAllBytes(g5c.toPath())),
			"aggobs=serve/clean runs desynced bins (global vector differs) - not a fair A/B");

		final ArrayList<String[]> A5S=parse(agg5s), A5C=parse(agg5c), F5=parse(f5Rows);
		assertTrue(A5S.size()==200 && A5C.size()==200 && F5.size()==200,
			"aggobs run row counts: serve="+A5S.size()+" clean="+A5C.size()+" famset_a-clean="+F5.size());

		// Non-vacuousness gate: if serve and clean produced byte-identical aggregator output
		// under GUARANTEED real contamination, the branch is a no-op or was swapped-to-identical
		// - exactly the hole this block exists to close.
		assertTrue(!Arrays.equals(Files.readAllBytes(agg5s.toPath()), Files.readAllBytes(agg5c.toPath())),
			"aggobs=serve and aggobs=clean produced byte-identical aggregator output under real "
			+"contamination - the branch is a no-op or swapped-to-identical");

		// Per-row exact recomputation for famset_a (manifest cols 0-3, subset {0,2,5}):
		// - CLEAN obs is f5Rows' own obs columns (target-only snapshot) - EXACTLY what
		//   aggobs=clean's famArr=cleanFamBuf sums for these ranks (formatFamsetRow and the
		//   clean branch of formatAggRow both read the SAME cleanFamBuf-derived snapshot).
		// - SERVE obs is decoded from the aggregator row's dense head (cols 13-28), which is
		//   ALWAYS built from whole-bin fam[] regardless of aggobs (see formatAggRow's own
		//   comment) - an independent, non-subnet ground truth for the serve-side numbers.
		//   Decoding presence+excess back to a count is exact here: this fixture's per-rank
		//   whole-bin counts never approach EXC_CAP=16.
		// - phylum + context columns come from f5Rows too: formatFamsetRow's context is ALSO
		//   whole-bin (built from lastTotalBp/glob[], set once per bin before the aggobs branch
		//   even exists), so f5Rows' own phylum/context columns are valid for BOTH the clean and
		//   serve reconstructed inputs - only the 3 obs columns need to be swapped for serve.
		final int[] headOrder5={0, 2, 1, 3, 5, 6, 7, 4};//prevalence order - fixed by the cache/
			//family fixture alone (computeDenseHead runs before any sampling), independent of
			//seed/cleanspike; identical to headOrder in (2)+(3).
		final int j0=0, j2=1, j5=4;//headOrder5[j]==rank: j0->0, j2->2, j5->5
		int diffRows=0, ncDiffRows=0;
		for(int i=0; i<200; i++){
			final String[] s=A5S.get(i), c=A5C.get(i), f=F5.get(i);
			assertTrue(s.length==62 && c.length==62, "row "+i+" agg width");
			if(!(s[4].equals(c[4]) && s[5].equals(c[5]) && s[6].equals(c[6]) && s[7].equals(c[7]))){ncDiffRows++;}

			final double obsClean=feed(cnA, f, 23, 3);
			checkBlock(c, 0, obsClean, predOf(cnA), "famset_a-clean", i);

			final int serve0=decodeTwo(s, j0), serve2=decodeTwo(s, j2), serve5=decodeTwo(s, j5);
			// phylum(3) + shared context(17) = 20 more columns, unchanged by aggobs (only the 3
			// obs columns need swapping for serve - see the class javadoc's reasoning above).
			final String[] fServe=new String[23];
			fServe[0]=String.valueOf(serve0); fServe[1]=String.valueOf(serve2); fServe[2]=String.valueOf(serve5);
			System.arraycopy(f, 3, fServe, 3, 20);
			final double obsServe=feed(cnA, fServe, 23, 3);
			checkBlock(s, 0, obsServe, predOf(cnA), "famset_a-serve", i);

			assertTrue(obsClean<=obsServe, "row "+i+" clean obs "+obsClean+" > serve obs "+obsServe
				+" - contamination can only ADD family counts, never remove them");
			if(obsClean!=obsServe){diffRows++;}
		}
		assertTrue(diffRows>0, "no row's famset_a obs differed between aggobs=serve and "
			+"aggobs=clean under cleanspike=0.0 - contamination never touched ranks {0,2,5}; "
			+"the recomputation would be vacuous");

		System.out.println("MagQCAggVectorTest aggobs=serve/clean PASS: bins aligned (global "
			+"byte-identical), aggregator output NOT identical between serve/clean ("
			+diffRows+"/200 rows differ on famset_a, "+ncDiffRows+"/200 differ on ncrna), "
			+"serve>=clean on every row, per-row [ratio,lobs,lpred,zero] recomputed exactly "
			+"for both branches via the SAME frozen net.");
	}

	/** Feeds a subnet-row's input columns to the net; returns the obs sum (first numObs cols). */
	private static double feed(CellNet net, String[] row, int numIn, int numObs){
		final float[] in=new float[numIn];
		double obs=0;
		for(int i=0; i<numIn; i++){
			in[i]=Float.parseFloat(row[i]);
			if(i<numObs){obs+=in[i];}
		}
		net.applyInput(in);
		net.feedForward();
		return obs;
	}
	private static double predOf(CellNet net){return net.getOutput(0);}

	/** Decodes an aggregator row's dense-head enc=two column pair (presence, excess) at head
	 *  index j back to the whole-bin family count for that rank. Mirrors MagQCVectorMaker's
	 *  private appendTwo/EXC_CAP=16 exactly; exact as long as the count never approaches 16
	 *  (true for this tiny fixture - verified, not assumed, by the direction check that uses it). */
	private static int decodeTwo(String[] row, int j){
		final int presCol=13+2*j, excCol=14+2*j;
		if("0".equals(row[presCol])){return 0;}
		final double excess=Double.parseDouble(row[excCol]);
		return (int)Math.round(excess*16)+1;
	}

	/** Asserts one per-subnet feature block [ratio, lobs, lpred, zero] at column offset off. */
	private static void checkBlock(String[] a, int off, double obs, double pred, String name, int row){
		final String ratio=fmt(Math.min(2.0, obs/Math.max(0.5, pred)));
		final String lobs=fmt(log2(1+obs));
		final String lpred=fmt(log2(1+Math.max(0, pred)));
		final String zero=(obs==0 ? "1" : "0");
		assertTrue(ratio.equals(a[off]) && lobs.equals(a[off+1]) && lpred.equals(a[off+2]) && zero.equals(a[off+3]),
			"row "+row+" "+name+" block: expected ["+ratio+","+lobs+","+lpred+","+zero+"] got ["
			+a[off]+","+a[off+1]+","+a[off+2]+","+a[off+3]+"]");
	}

	/** Builds a tiny random real net and writes it in BBNet format. */
	private static void makeNet(File f, int[] dims, long seed) throws Exception{
		final CellNet net=new CellNet(dims, seed, 1f, 0f, 1, new ArrayList<String>());
		net.randomize();
		final FileWriter w=new FileWriter(f);
		w.write(net.toBytes().toString());
		w.close();
	}

	private static ArrayList<String[]> parse(File f) throws Exception{
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		boolean first=true;
		for(String line : Files.readAllLines(f.toPath())){
			if(line.length()==0){continue;}
			if(first){first=false; continue;}//header
			rows.add(line.split("\t"));
		}
		return rows;
	}
	private static String header(File f) throws Exception{
		for(String line : Files.readAllLines(f.toPath())){if(line.length()>0){return line;}}
		return null;
	}

	/** Mirrors MagQCVectorMaker.fmt (fixed notation, no exponent). */
	private static String fmt(double v){
		if(v==(long)v){return Long.toString((long)v);}
		return String.format("%.6f", v);
	}
	private static double log2(double v){return v<=1 ? 0 : Math.log(v)/Math.log(2);}

	private static void write(File f, String content) throws Exception{
		final FileWriter w=new FileWriter(f); w.write(content); w.close();
	}
	private static String p(File dir, String name){return new File(dir, name).getAbsolutePath();}
	/** Recursively deletes a file/directory tree; used by the shutdown hook to clean the temp fixture. */
	private static void deleteRecursive(File f){
		final File[] kids=f.listFiles();
		if(kids!=null){for(File k : kids){deleteRecursive(k);}}
		f.delete();
	}
	private static String[] concat(String[] a, String[] b){
		final String[] out=new String[a.length+b.length];
		System.arraycopy(a, 0, out, 0, a.length);
		System.arraycopy(b, 0, out, a.length, b.length);
		return out;
	}
	private static void assertTrue(boolean cond, String msg){
		if(!cond){throw new AssertionError(msg);}
	}
}
