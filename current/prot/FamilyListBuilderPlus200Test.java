package prot;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;

/**
 * Self-contained test for {@link FamilyListBuilder}'s global single-copy backbone
 * extension (Rebuild Point 2, magqc_rebuild_20260824.plan / records/FAMILYLIST_PLUS200_SPEC.md).
 * Builds a tiny synthetic corpus (2 phyla x 6 genomes) with hand-computable occ_total/occ_single
 * for every family, and verifies:
 *
 * <ul>
 * <li><b>Backward-compat:</b> globalsingle=0 reproduces the EXACT pre-extension per-phylum set
 * (gate criterion 5) -- {@code F_p1}/{@code F_p2}, each single-copy in all 6 genomes of its
 * phylum, occ_total=6.</li>
 * <li><b>Correct pick with an adequate candidate pool:</b> {@code F_sneaky} (single-copy in 2
 * genomes of EACH phylum, occ_total=4, occ_single=4) never clears either phylum's 50% gate
 * (2/6=33%) so the per-phylum step never sees it -- exactly the case the extension exists for.
 * With globalsinglecandidates=4 (wide enough to include it), globalsingle=1 must select it.</li>
 * <li><b>Disjointness</b> (gate criterion 1): the addition is never one of the per-phylum
 * selected reps.</li>
 * <li><b>occSingleOut merge correctness:</b> {@code F_sneaky} lands ABOVE nothing but alongside
 * the originals in occ_total-desc order without the sanity-print's occSingle lookup NPEing --
 * a real gap this session's first version had (occSingle was only ever populated for the
 * original per-phylum candidate universe).</li>
 * <li><b>The safety assert actually fires, not just decoratively</b>: three junk families
 * ({@code J1,J2,J3}, present 2-copies-per-genome in ALL 12 genomes -- occ_total=12, occ_single=0
 * since never single-copy) rank ABOVE {@code F_sneaky} by raw occ_total. With
 * globalsinglecandidates=3 (just enough to admit only the junk), the candidate pool never sees
 * F_sneaky at all, so s1=0 &lt; T=12 -- KillSwitch.assertDie MUST fire (verified via a subprocess:
 * assertDie halts the JVM, so this cannot be caught in-process).</li>
 * </ul>
 *
 * Run with assertions on: {@code java -ea prot.FamilyListBuilderPlus200Test}
 *
 * @author Eru
 */
public class FamilyListBuilderPlus200Test {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("flb_plus200_test").toFile();
		Runtime.getRuntime().addShutdownHook(new Thread(()->deleteRecursive(dir)));

		// taxpgm: tid 1-6 = P1, tid 7-12 = P2 (minPhylumGenomes default 3, both qualify at 6).
		final StringBuilder taxpgm=new StringBuilder();
		for(int t=1; t<=6; t++){taxpgm.append(t).append("\tP1\n");}
		for(int t=7; t<=12; t++){taxpgm.append(t).append("\tP2\n");}
		write(new File(dir, "taxpgm.tsv"), taxpgm.toString());
		write(new File(dir, "excluded.tsv"), "");

		final StringBuilder cluster=new StringBuilder();
		// F_p1: single-copy in ALL 6 P1 genomes -> occ_total=6, occ_in_phylum(P1)=6/6=100%.
		for(int t=1; t<=6; t++){cluster.append("F_p1\ttid_").append(t).append("_g1\n");}
		// F_p2: single-copy in ALL 6 P2 genomes -> occ_total=6, occ_in_phylum(P2)=6/6=100%.
		for(int t=7; t<=12; t++){cluster.append("F_p2\ttid_").append(t).append("_g1\n");}
		// J1,J2,J3: 2 copies in ALL 12 genomes -> occ_total=12, occ_single=0 (never single-copy) --
		// high-occ_total "junk" that fills a small candidate pool without being a real backbone family.
		for(String j : new String[]{"J1", "J2", "J3"}){
			for(int t=1; t<=12; t++){
				cluster.append(j).append("\ttid_").append(t).append("_g1\n");
				cluster.append(j).append("\ttid_").append(t).append("_g2\n");
			}
		}
		// F_sneaky: single-copy in tid1,2 (P1, 2/6=33%<50% gate) and tid7,8 (P2, 2/6=33%<50% gate)
		// -> occ_total=4, occ_single=4. Clears NEITHER phylum's gate, so the per-phylum step never
		// sees it -- the exact scenario the global extension exists to catch.
		for(int t : new int[]{1, 2, 7, 8}){cluster.append("F_sneaky\ttid_").append(t).append("_g1\n");}
		write(new File(dir, "cluster.tsv"), cluster.toString());

		final StringBuilder reps=new StringBuilder();
		for(String r : new String[]{"F_p1", "F_p2", "J1", "J2", "J3", "F_sneaky"}){
			reps.append(">").append(r).append("\nMKV\n");
		}
		write(new File(dir, "reps.fasta"), reps.toString());

		final String cl=p(dir, "cluster.tsv"), rp=p(dir, "reps.fasta"), tp=p(dir, "taxpgm.tsv"), ex=p(dir, "excluded.tsv");
		final String[] common={"cluster="+cl, "reps="+rp, "taxpgm="+tp, "excluded="+ex,
			"ntop=1", "nbot=1", "minoccfrac=0.5", "minphylumgenomes=3"};

		// (1) Backward-compat: globalsingle=0 (or omitted) reproduces the exact per-phylum set.
		final File baseOut=new File(dir, "base.tsv");
		FamilyListBuilder.main(concat(common, new String[]{"out="+baseOut.getAbsolutePath(), "ow=t"}));
		final ArrayList<String[]> baseRows=parse(baseOut);
		assertTrue(baseRows.size()==2, "base set size "+baseRows.size()+" != 2");
		assertTrue(hasRow(baseRows, "F_p1", 6) && hasRow(baseRows, "F_p2", 6),
			"base set is not exactly {F_p1:6, F_p2:6}: "+baseRows);

		// (2) globalsingle=1 with an ADEQUATE candidate pool must pick F_sneaky.
		final File okOut=new File(dir, "ok.tsv");
		FamilyListBuilder.main(concat(common, new String[]{"out="+okOut.getAbsolutePath(), "ow=t",
			"globalsingle=1", "globalsinglecandidates=4"}));
		final ArrayList<String[]> okRows=parse(okOut);
		assertTrue(okRows.size()==3, "extended set size "+okRows.size()+" != 3: "+okRows);
		assertTrue(hasRow(okRows, "F_p1", 6) && hasRow(okRows, "F_p2", 6) && hasRow(okRows, "F_sneaky", 4),
			"extended set is not exactly {F_p1:6, F_p2:6, F_sneaky:4}: "+okRows);
		// (2b) Disjointness (gate criterion 1): the addition was never a per-phylum pick.
		assertTrue(!hasRow(baseRows, "F_sneaky", 4), "F_sneaky was already in the base set -- test fixture is wrong");

		// (3) The safety assert must actually FIRE with an INADEQUATE pool (T=12 from the junk
		// families > s1=0, since F_sneaky -- the true best pick, occ_single=4 -- never enters a
		// pool of size 3). assertDie halts the JVM, so this runs as a subprocess.
		final File badOut=new File(dir, "bad.tsv");
		final ProcessBuilder pb=new ProcessBuilder(
			System.getProperty("java.home")+"/bin/java", "-ea", "--add-modules", "jdk.incubator.vector",
			"-cp", System.getProperty("java.class.path"), "prot.FamilyListBuilder",
			"cluster="+cl, "reps="+rp, "taxpgm="+tp, "excluded="+ex,
			"ntop=1", "nbot=1", "minoccfrac=0.5", "minphylumgenomes=3",
			"out="+badOut.getAbsolutePath(), "ow=t", "globalsingle=1", "globalsinglecandidates=3");
		pb.redirectErrorStream(true);
		final Process proc=pb.start();
		final String out=new String(proc.getInputStream().readAllBytes());
		final int code=proc.waitFor();
		assertTrue(code!=0, "safety-assert subprocess exited 0 -- the crash-loud guard did not fire "
			+"when the candidate pool was too small to admit the true best pick");
		assertTrue(out.contains("safety bound violated") && out.contains("T=12") && out.contains("s1=0"),
			"safety-assert subprocess exited nonzero but not with the expected message: "+out);

		System.out.println("FamilyListBuilderPlus200Test PASS: globalsingle=0 backward-compat exact, "
			+"globalsingle=1 with an adequate pool correctly picks F_sneaky (disjoint from the base "
			+"set, occ_single correctly merged for the report), and the safety assert genuinely fires "
			+"(exit "+code+") when the candidate pool is too small to admit the true best candidate.");
	}

	private static boolean hasRow(ArrayList<String[]> rows, String rep, int occTotal){
		for(String[] r : rows){if(r[1].equals(rep) && Integer.parseInt(r[2])==occTotal){return true;}}
		return false;
	}

	private static ArrayList<String[]> parse(File f) throws Exception{
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		for(String line : Files.readAllLines(f.toPath())){
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			rows.add(line.split("\t"));
		}
		return rows;
	}

	private static void write(File f, String content) throws Exception{
		final FileWriter w=new FileWriter(f); w.write(content); w.close();
	}
	private static String p(File dir, String name){return new File(dir, name).getAbsolutePath();}
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
