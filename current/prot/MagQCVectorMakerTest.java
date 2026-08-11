package prot;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Self-contained test for {@link MagQCVectorMaker}'s ncRNA-subnet emission
 * (composite architecture, slice 1). Writes a tiny synthetic per-contig cache with
 * KNOWN ncRNA complements, runs the generator, and verifies:
 *
 * <ul>
 * <li><b>Global output unaffected:</b> the family/global/phylum training vector is
 * byte-identical whether or not the ncRNA subnet is also emitted.</li>
 * <li><b>Subset-relative-label invariant:</b> every subnet row's observed ncRNA
 * total never exceeds its native-complement target (a bin's target contigs are a
 * subset of the genome).</li>
 * <li><b>Target correctness:</b> each row's target equals a fixture organism's true
 * native ncRNA total.</li>
 * <li><b>The crux:</b> the same organism at different completeness yields the SAME
 * target — the denominator is organism-intrinsic, not bin-dependent (this is why
 * "50% is not trainable, 37.5% is").</li>
 * </ul>
 *
 * Run with assertions on: {@code java -ea prot.MagQCVectorMakerTest}
 *
 * @author Eru
 */
public class MagQCVectorMakerTest {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("magqcvm_test").toFile();
		dir.deleteOnExit();

		// --- fixture: 4 organisms, 2 per phylum, with hand-computable native ncRNA totals ---
		// cache cols: contigid tid domain len gc acgt cds mapped glenSum glenSq coding r16 r23 r5 rother trna famcounts
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

		// Native ncRNA totals (r16+r23+r5+rother+trna over each genome), hand-computed:
		// 100: (1,2,1,1,23)=28  200: (1,1,1,2,21)=26  300: (1,1,1,1,16)=20  400: (1,2,1,1,23)=28
		final HashMap<Integer,Integer> nativeTotal=new HashMap<Integer,Integer>();
		nativeTotal.put(100, 28); nativeTotal.put(200, 26); nativeTotal.put(300, 20); nativeTotal.put(400, 28);

		final String cache=p(dir, "cache.tsv"), size=p(dir, "sizemap.tsv");
		final String tax=p(dir, "taxpgm.tsv"), fam=p(dir, "familylist.tsv");

		// Run A: global only. Run B: global + ncRNA subnet. Same seed/args otherwise.
		final File gA=new File(dir, "gA.tsv"), gB=new File(dir, "gB.tsv");
		final File sn=new File(dir, "ncrna.tsv");
		final String[] common={"cache="+cache, "sizemap="+size, "taxpgm="+tax, "familylist="+fam,
			"n=200", "valn=0", "valfrac=0.5", "enc=two", "seed=1"};
		MagQCVectorMaker.main(concat(common, new String[]{"out="+gA.getAbsolutePath()}));
		MagQCVectorMaker.main(concat(common, new String[]{"out="+gB.getAbsolutePath(),
			"subnet=ncrna", "subnetout="+sn.getAbsolutePath()}));

		// (1) Global output byte-identical with the subnet on/off.
		final byte[] ba=Files.readAllBytes(gA.toPath()), bb=Files.readAllBytes(gB.toPath());
		assertTrue(java.util.Arrays.equals(ba, bb),
			"global vector changed when the subnet was enabled ("+ba.length+" vs "+bb.length+" bytes)");

		// Parse the subnet file.
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		String header=null;
		for(String line : Files.readAllLines(sn.toPath())){
			if(line.length()==0){continue;}
			if(header==null){header=line; continue;}
			rows.add(line.split("\t"));
		}
		// (2) Header: 5 obs + 3 phyla + 5 context = 13 inputs, 1 output.
		assertTrue("#dims\t13\t1\t0".equals(header), "bad subnet header: "+header);
		assertTrue(rows.size()==200, "expected 200 subnet rows, got "+rows.size());

		// (3)+(4): invariants + crux.
		final HashMap<Integer,java.util.HashSet<Integer>> obsByTarget=new HashMap<Integer,java.util.HashSet<Integer>>();
		for(String[] r : rows){
			assertTrue(r.length==14, "row width "+r.length+" != 14");
			int obsTotal=0;
			for(int i=0; i<5; i++){obsTotal+=Integer.parseInt(r[i]);}
			final int target=Integer.parseInt(r[13]);
			// (2-invariant) observed can never exceed the native complement.
			assertTrue(obsTotal<=target, "observed "+obsTotal+" > target "+target);
			// (3) target must be some fixture organism's true native total.
			assertTrue(nativeTotal.containsValue(target), "target "+target+" is not a fixture native total");
			java.util.HashSet<Integer> set=obsByTarget.get(target);
			if(set==null){obsByTarget.put(target, set=new java.util.HashSet<Integer>());}
			set.add(obsTotal);
		}
		// (4) THE CRUX: at least one target (organism) appears with >=2 distinct observed totals,
		// i.e. the same denominator is reused across different completeness levels.
		boolean cruxSeen=false;
		for(java.util.HashSet<Integer> set : obsByTarget.values()){if(set.size()>=2){cruxSeen=true; break;}}
		assertTrue(cruxSeen, "crux not demonstrated: no organism appeared at >=2 completeness levels");

		System.out.println("MagQCVectorMakerTest PASS: global byte-identical; "+rows.size()+
			" ncRNA-subnet rows; observed<=native on all; targets in "+obsByTarget.keySet()+
			"; crux (same denominator, varied observed) demonstrated.");
	}

	private static void write(File f, String content) throws Exception{
		FileWriter w=new FileWriter(f); w.write(content); w.close();
	}
	private static String p(File dir, String name){return new File(dir, name).getAbsolutePath();}
	private static String[] concat(String[] a, String[] b){
		String[] out=new String[a.length+b.length];
		System.arraycopy(a, 0, out, 0, a.length);
		System.arraycopy(b, 0, out, a.length, b.length);
		return out;
	}
	private static void assertTrue(boolean cond, String msg){
		if(!cond){throw new AssertionError(msg);}
	}
}
