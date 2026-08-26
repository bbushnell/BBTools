package prot;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;

/**
 * Self-contained test for {@link CacheBuilder}'s PRODUCER side -- the anticodon-tag parse and
 * the per-shred dimer-tracker reuse -- closing the gap UMP45 flagged: the three existing
 * differential tests ({@link MagQCVectorMakerTest}, {@link MagQCAggVectorTest},
 * {@link FamilyListBuilderPlus200Test}) only cover the emit/consumer side of the pipeline; none
 * of them exercises {@code CacheBuilder} itself.
 *
 * <ul>
 * <li><b>Anticodon parse</b> ({@link CacheBuilder#parseAnticodonCode}): a present, valid
 * ACGT-only tag decodes to the correct packed 2-bit-per-base code (0-63); an absent tag, an
 * ambiguous base (N) within the tag, and a truncated tag (fewer than 3 letters remaining) all
 * return -1 ("not recorded"), matching the ~31%-missing-in-real-data tolerance the class
 * javadoc documents -- verified directly on the byte-offset contract (attrA/attrLen into a
 * larger line buffer, not a standalone string), the same calling convention {@code loadGff}
 * actually uses.</li>
 * <li><b>Dimer-tracker clear-between-shreds correctness</b> (the reused {@code KmerTracker}
 * in {@code emit()}, {@code CacheBuilder.java:383,408}): runs the REAL {@code CacheBuilder}
 * end-to-end (not a reimplementation) over a synthetic 2-shred fixture with hand-computable,
 * maximally-distinguishing compositions -- shred A = "AAAA" (3x AA, dimer index 0b0000=0),
 * shred B = "CCCC" (3x CC, dimer index 0b0101=5) -- and asserts shred B's emitted dimers field
 * contains ONLY its own CC=3, with ZERO leaked AA count from shred A. If {@code clearAll()}
 * were missing or misplaced, shred B's row would show AA=3 leaked in from shred A's
 * accumulation -- this is exactly the failure mode the test is built to catch.</li>
 * </ul>
 *
 * Run with assertions on: {@code java -ea prot.CacheBuilderTest}
 *
 * @author Eru
 */
public class CacheBuilderTest {

	public static void main(String[] args) throws Exception{
		testAnticodonParse();
		testDimerClearBetweenShreds();
		System.out.println("CacheBuilderTest PASS: anticodon parse correct on present/absent/ambiguous/"
			+"truncated tags (byte-offset contract verified), and the reused per-shred dimer tracker "
			+"shows zero cross-shred leakage (clearAll() between shreds verified end-to-end via a real "
			+"CacheBuilder run, not a reimplementation).");
	}

	// ---------------------------------------------------------------------------------------
	// Part 1: parseAnticodonCode -- pure unit tests on the byte-offset contract.
	// ---------------------------------------------------------------------------------------

	private static void testAnticodonParse(){
		// A/C/G/T = 0/1/2/3 (dna.AminoAcid.baseToNumber). code=(b0<<4)|(b1<<2)|b2.

		// (a) Present, valid tag, offset into a larger buffer (prefix bytes before the attribute
		// column -- the real calling convention: attrA=lp.a(), not 0).
		{
			String prefix="chr1\tcallgenes\ttRNA\t100\t200\t.\t+\t.\t";
			String attrs="ID=1;ac=Gly;anticodon:AAC";
			byte[] line=(prefix+attrs).getBytes();
			int attrA=prefix.length(), attrLen=attrs.length();
			int code=CacheBuilder.parseAnticodonCode(line, attrA, attrLen);
			int expected=(0<<4)|(0<<2)|1; // A=0,A=0,C=1
			assertTrue(code==expected, "AAC should decode to "+expected+", got "+code);
		}

		// (b) A different real base combination, tag NOT at the start of the attribute column.
		{
			String attrs="ID=2;note=x;anticodon:TGA;extra=y";
			byte[] line=attrs.getBytes();
			int code=CacheBuilder.parseAnticodonCode(line, 0, attrs.length());
			int expected=(3<<4)|(2<<2)|0; // T=3,G=2,A=0
			assertTrue(code==expected, "TGA should decode to "+expected+", got "+code);
		}

		// (c) Absent tag entirely -> -1 ("not recorded", the ~31%-missing real-world case).
		{
			String attrs="ID=3;type=tRNA;note=no_anticodon_here";
			byte[] line=attrs.getBytes();
			int code=CacheBuilder.parseAnticodonCode(line, 0, attrs.length());
			assertTrue(code==-1, "absent tag should return -1, got "+code);
		}

		// (d) Tag present but with an ambiguous base (N) -> -1, not a silently-wrong code.
		{
			String attrs="ID=4;anticodon:NAC";
			byte[] line=attrs.getBytes();
			int code=CacheBuilder.parseAnticodonCode(line, 0, attrs.length());
			assertTrue(code==-1, "ambiguous base (N) in the tag should return -1, got "+code);
		}

		// (e) Tag present but truncated (fewer than 3 letters remain before the column ends)
		// -> the "i+tagLen+3<=attrLen" bound must reject it, not read past the buffer or return
		// a garbage code from whatever bytes happen to follow.
		{
			String attrs="ID=5;anticodon:AC"; // only 2 letters after the tag
			byte[] line=attrs.getBytes();
			int code=CacheBuilder.parseAnticodonCode(line, 0, attrs.length());
			assertTrue(code==-1, "truncated tag (2 letters) should return -1, got "+code);
		}

		// (f) All four DNA bases represented across codes, cross-checked against the documented
		// packing (2 bits/base via dna.AminoAcid.baseToNumber, A=0/C=1/G=2/T=3).
		{
			String[] codons={"AAA","CCC","GGG","TTT","ACG","TGC"};
			int[] expectedCodes={0, (1<<4)|(1<<2)|1, (2<<4)|(2<<2)|2, (3<<4)|(3<<2)|3,
				(0<<4)|(1<<2)|2, (3<<4)|(2<<2)|1};
			for(int i=0; i<codons.length; i++){
				String attrs="anticodon:"+codons[i];
				byte[] line=attrs.getBytes();
				int code=CacheBuilder.parseAnticodonCode(line, 0, attrs.length());
				assertTrue(code==expectedCodes[i], codons[i]+" should decode to "+expectedCodes[i]+", got "+code);
				assertTrue(code>=0 && code<=63, codons[i]+" code "+code+" out of the documented 0-63 range");
			}
		}
	}

	// ---------------------------------------------------------------------------------------
	// Part 2: dimer-tracker clear-between-shreds -- real CacheBuilder end-to-end run.
	// ---------------------------------------------------------------------------------------

	private static void testDimerClearBetweenShreds() throws Exception{
		final File dir=Files.createTempDirectory("cachebuilder_test").toFile();
		Runtime.getRuntime().addShutdownHook(new Thread(()->deleteRecursive(dir)));

		// Two shreds, same tid (both bacteria), maximally distinguishing compositions:
		// shred A = "AAAA" -> 3x AA (dimer index 0). shred B = "CCCC" -> 3x CC (dimer index 5).
		// If the reused KmerTracker isn't cleared between shreds, B's row leaks A's AA=3.
		final String shreds=">tid_1_testA\nAAAA\n>tid_1_testB\nCCCC\n";
		write(new File(dir, "shreds.fa"), shreds);

		// No genes needed -- dimers are computed directly from the shred sequence in emit(),
		// independent of the GFF/hits accumulator. Empty (comment-only) files satisfy the
		// required non-null args without asserting anything about gene content.
		write(new File(dir, "empty.gff"), "#empty\n");
		write(new File(dir, "empty.m8"), "");
		write(new File(dir, "empty.familylist.tsv"), "#rank\trep\n");

		// loadTidSet requires >=1 matching *.fna.gz per directory (even if unused) -- both dirs
		// must be non-empty or CacheBuilder's constructor-time assert fires before process() runs.
		final File bacteria4Dir=new File(dir, "bacteria4"); bacteria4Dir.mkdir();
		write(new File(bacteria4Dir, "tid_1_test.fna.gz"), ""); // tid=1, matches both shreds
		final File archaea4Dir=new File(dir, "archaea4"); archaea4Dir.mkdir();
		write(new File(archaea4Dir, "tid_999_dummy.fna.gz"), ""); // unrelated tid, just satisfies the non-empty gate

		final File out=new File(dir, "cache.tsv");
		CacheBuilder.main(new String[]{
			"shreds="+p(dir, "shreds.fa"), "gff="+p(dir, "empty.gff"), "hits="+p(dir, "empty.m8"),
			"familylist="+p(dir, "empty.familylist.tsv"),
			"archaea4="+archaea4Dir.getAbsolutePath(), "bacteria4="+bacteria4Dir.getAbsolutePath(),
			"out="+out.getAbsolutePath(), "topn=1", "ow=t"
		});

		final ArrayList<String[]> rows=new ArrayList<String[]>();
		for(String line : Files.readAllLines(out.toPath())){
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			rows.add(line.split("\t"));
		}
		assertTrue(rows.size()==2, "expected 2 shred rows, got "+rows.size());

		int[] dimersA=null, dimersB=null;
		for(String[] r : rows){
			// 19-field contract, 0-indexed: contig_id0 ... dimers18 (last field).
			final String contigId=r[0];
			final int[] dimers=parseDimers(r[r.length-1]);
			if(contigId.equals("tid_1_testA")){dimersA=dimers;}
			else if(contigId.equals("tid_1_testB")){dimersB=dimers;}
		}
		assertTrue(dimersA!=null, "shred A row not found");
		assertTrue(dimersB!=null, "shred B row not found");

		// Shred A: exactly 3 AA (index 0), zero everywhere else.
		assertTrue(dimersA[0]==3, "shred A: AA (index 0) should be 3, got "+dimersA[0]);
		for(int i=1; i<16; i++){
			assertTrue(dimersA[i]==0, "shred A: dimer index "+i+" should be 0, got "+dimersA[i]);
		}

		// Shred B: exactly 3 CC (index 5=0b0101), and CRITICALLY zero AA (index 0) -- if the
		// tracker weren't cleared after shred A, this would show the leaked 3.
		assertTrue(dimersB[5]==3, "shred B: CC (index 5) should be 3, got "+dimersB[5]);
		assertTrue(dimersB[0]==0, "shred B: AA (index 0) should be 0 (no leakage from shred A), got "+dimersB[0]
			+" -- if nonzero, the reused dimer tracker was NOT cleared between shreds");
		for(int i=0; i<16; i++){
			if(i==5){continue;}
			assertTrue(dimersB[i]==0, "shred B: dimer index "+i+" should be 0, got "+dimersB[i]);
		}
	}

	private static int[] parseDimers(String csv){
		String[] parts=csv.split(",");
		int[] out=new int[parts.length];
		for(int i=0; i<parts.length; i++){out[i]=Integer.parseInt(parts[i]);}
		return out;
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
	private static void assertTrue(boolean cond, String msg){
		if(!cond){throw new AssertionError(msg);}
	}
}
