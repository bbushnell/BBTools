package prot;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;

/**
 * Self-contained test for {@link AntiSetMiner}. Builds a 12-organism fixture with two
 * DESIGNED mutually-exclusive groups plus band-excluded distractors, and verifies the
 * greedy set-cover recovers exactly the designed structure:
 *
 * <ul>
 * <li>An exclusive PAIR: f0 (orgs 0-5) + f1 (orgs 6-11), union 100%.</li>
 * <li>An exclusive TRIPLE: f2 (0-3) + f3 (4-7) + f4 (8-11), union 100%.</li>
 * <li>f5 present in ALL organisms — above maxprev, must be excluded (it's a marker,
 * not an alternative).</li>
 * <li>f6 present in one organism — below minprev, must be excluded (noise).</li>
 * <li>The pooled subset file contains exactly the 5 grouped families; the matched
 * control is empty because the band holds nothing else (edge case exercised).</li>
 * </ul>
 *
 * Run with assertions on: {@code java -ea prot.AntiSetMinerTest}
 *
 * @author UMP45
 */
public class AntiSetMinerTest {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("antiset_test").toFile();
		//deleteOnExit() only removes an EMPTY dir; the fixture files inside would leak.
		//A shutdown hook recursively deletes the whole tree, even on assertion failure.
		Runtime.getRuntime().addShutdownHook(new Thread(()->deleteRecursive(dir)));

		//Per-org cache rows: tid domain bp gc acgt cds mapped glen glensq coding r16 r23 r5 rother trna famcounts
		final StringBuilder sb=new StringBuilder();
		sb.append("#tid\tdomain\tgenome_bp\tgc\tacgt\tcds\tmapped\tglen\tglensq\tcoding\tr16\tr23\tr5\trother\ttrna\tfamcounts\n");
		for(int o=0; o<12; o++){
			final ArrayList<Integer> fams=new ArrayList<Integer>();
			if(o<6){fams.add(0);}else{fams.add(1);}          //exclusive pair
			if(o<4){fams.add(2);}else if(o<8){fams.add(3);}else{fams.add(4);} //exclusive triple
			fams.add(5);                                      //universal (above band)
			if(o==0){fams.add(6);}                            //singleton (below band)
			final StringBuilder fc=new StringBuilder();
			for(int i=0; i<fams.size(); i++){if(i>0){fc.append(';');} fc.append(fams.get(i)).append(":1");}
			sb.append((1000+o)+"\tBacteria\t5000000\t2500000\t5000000\t5000\t2000\t4500000\t999\t4500000\t1\t1\t1\t0\t30\t"+fc+"\n");
		}
		write(new File(dir, "cache.tsv"), sb.toString());
		write(new File(dir, "familylist.tsv"), "rank\trep\n0\tf0\n1\tf1\n2\tf2\n3\tf3\n4\tf4\n5\tf5\n6\tf6\n7\tf7\n");

		final String outFile=p(dir, "groups.tsv"), prefix=p(dir, "anti");
		AntiSetMiner.main(new String[]{"cache="+p(dir, "cache.tsv"), "familylist="+p(dir, "familylist.tsv"),
			"out="+outFile, "subsetprefix="+prefix, "minprev=0.2", "maxprev=0.90",
			"uniontarget=0.98", "unionmin=0.95", "seed=1"});

		//Parse the groups report.
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		for(String line : Files.readAllLines(new File(outFile).toPath())){
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			rows.add(line.split("\t"));
		}
		assertTrue(rows.size()==2, "expected 2 groups, got "+rows.size());
		//Group 0 = the pair {0,1} (seeded from the highest-prevalence candidate).
		assertTrue(rows.get(0)[1].equals("2"), "group0 size: "+rows.get(0)[1]);
		assertTrue(memberRanks(rows.get(0)[4]).equals("0,1"), "group0 members: "+rows.get(0)[4]);
		assertTrue(rows.get(0)[2].startsWith("1.0"), "group0 union prev: "+rows.get(0)[2]);
		//Group 1 = the triple {2,3,4}.
		assertTrue(rows.get(1)[1].equals("3"), "group1 size: "+rows.get(1)[1]);
		assertTrue(memberRanks(rows.get(1)[4]).equals("2,3,4"), "group1 members: "+rows.get(1)[4]);
		assertTrue(rows.get(1)[2].startsWith("1.0"), "group1 union prev: "+rows.get(1)[2]);
		//f5 (universal) and f6 (singleton) appear in NO group.
		for(String[] r : rows){
			assertTrue(!memberRanks(r[4]).contains("5") && !memberRanks(r[4]).contains("6"),
				"band-excluded family leaked into group: "+r[4]);
		}
		//Pooled subset file: exactly the 5 grouped families, group order.
		final ArrayList<String> pool=new ArrayList<String>();
		for(String line : Files.readAllLines(new File(prefix+"_pool0.txt").toPath())){
			if(line.length()>0){pool.add(line);}
		}
		assertTrue(String.join(",", pool).equals("0,1,2,3,4"), "pool0: "+pool);
		//Control: band minus pool is empty here, so the control is empty (edge case).
		final long ctrlLines=Files.readAllLines(new File(prefix+"_ctrl0.txt").toPath())
			.stream().filter(s -> s.length()>0).count();
		assertTrue(ctrlLines==0, "ctrl0 should be empty, has "+ctrlLines);

		System.out.println("AntiSetMinerTest PASS: exclusive pair {0,1} and triple {2,3,4} recovered "
			+"exactly (union 100% each), universal + singleton families band-excluded, pooled subset "
			+"= the 5 grouped families, empty-band control edge case handled.");
	}

	/** Extracts "rank,rank,..." (sorted ascending) from a "rank:prev;rank:prev" member string. */
	private static String memberRanks(String members){
		final ArrayList<Integer> l=new ArrayList<Integer>();
		for(String m : members.split(";")){l.add(Integer.parseInt(m.substring(0, m.indexOf(':'))));}
		java.util.Collections.sort(l);
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<l.size(); i++){if(i>0){sb.append(',');} sb.append(l.get(i));}
		return sb.toString();
	}

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
	private static void assertTrue(boolean cond, String msg){
		if(!cond){throw new AssertionError(msg);}
	}
}
