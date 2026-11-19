package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import prok.NcrnaVoteDiagLogSink;
import prok.NcrnaVoteDiagnostics;
import prok.RrnaResourceIO;
import structures.IntList;

/** Optional generated-vote TSV cells preserve real zeroes while omitting absent values. */
public final class NcrnaVoteDiagLogSinkTest {
	public static void main(String[] args)throws Exception{
		final Path dir=Files.createTempDirectory("ncrna_vote_diag_");
		final Path out=dir.resolve("vote.tsv");
		final ByteStreamWriter writer=RrnaResourceIO.writer(out.toString());
		final NcrnaVoteDiagLogSink log=new NcrnaVoteDiagLogSink(writer);
		log.group(new NcrnaVoteDiagnostics.Group("fam\tone", "ctg\none", 1, 100, 2, 7, "legacy",
			"dropped", "WHY", NcrnaVoteDiagnostics.UNKNOWN, NcrnaVoteDiagnostics.UNKNOWN,
			0, 0, false, NcrnaVoteDiagnostics.UNKNOWN, NcrnaVoteDiagnostics.UNKNOWN,
			10, 20, 5, 3, 0, new int[]{10,20}, new long[]{NcrnaVoteDiagnostics.UNKNOWN,0x1aL},
			new int[]{-1,0}, new int[]{-1,7}, -1, -1, false, -1, -1, false, -1, -1, false));
		log.group(new NcrnaVoteDiagnostics.Group("fam2", "ctg2", 0, 90, 1, 8, "vote",
			"retained", "VOTE", -5, 5, 12.25, 0, true, 12, 0, 12, 0, 0, 0, 1,
			new int[0], new long[0], new int[0], new int[0], 3.5, 4.5, true, 0, 1.25,
			true, 0, 2.5, true));
		boolean failed=false;
		try{
			log.group(new NcrnaVoteDiagnostics.Group("bad", "ctg", 0, 90, 1, 9, "vote",
				"bad", "NONFINITE", -5, 5, Double.NaN, 1, true, 0, 1, 0, 1, 0, 0, 1,
				new int[0], new long[0], new int[0], new int[0], 1, 1, true, 0, 0, true, 0, 0, true));
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Present nonfinite vote diagnostics must fail before serialization");
		RrnaResourceIO.finish(writer);
		final java.util.List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
		check(lines.size()==2, "Expected two successful vote rows");
		check(lines.get(0).equals("ncrna_vote_diag_v2\tfam one\tctg one\t1\t100\t2\t7\tlegacy\tdropped\tWHY"
			+"\t\t\t\t\t\t\t10\t20\t5\t3\t0\t10,20\t,1a\t,0\t,7\t\t\t\t\t\t"),
			"Absent vote values must be empty cells/list elements, not legacy dot tokens: "+lines.get(0));
		check(lines.get(1).equals("ncrna_vote_diag_v2\tfam2\tctg2\t0\t90\t1\t8\tvote\tretained\tVOTE"
			+"\t-5\t5\t12.25\t0.0\t12\t0\t12\t0\t0\t0\t1\t\t\t\t\t3.5\t4.5\t0.0\t1.25\t0.0\t2.5"),
			"Present zero-valued vote diagnostics must remain numeric: "+lines.get(1));
		final ArrayList<NcrnaVoteDiagnostics.Group> groups=new ArrayList<NcrnaVoteDiagnostics.Group>();
		final NcrnaVoteDiagnostics diag=new NcrnaVoteDiagnostics(groups::add, "fam", "ctg", 0, 80, 1, 4,
			2, new int[]{10,20,30}, new long[]{0x10L,0x20L,0x30L});
		diag.core(0, 2, new long[]{0,1}, new int[]{0,1}, new double[]{10,14}, new double[]{20,28},
			12, 24, 3, 5, 100, 200, false);
		final IntList fallback=new IntList();
		fallback.add(10);
		fallback.add(30);
		diag.fallback(fallback, 0, 2, 110, 210, true);
		check(groups.size()==2, "Expected one core group and one fallback group");
		check(groups.get(0).rawSet && groups.get(0).weightSet && groups.get(0).rangeSet
			&& !groups.get(0).sdSet && groups.get(0).leftRange==4 && groups.get(0).rightRange==8,
			"Core groups carry raw weights/ranges but omit top10 SD");
		check(!groups.get(1).rawSet && !groups.get(1).weightSet && !groups.get(1).rangeSet
			&& !groups.get(1).sdSet && groups.get(1).parents[0]==0 && groups.get(1).parents[1]<0,
			"Fallback groups preserve parent ancestry while omitting nonexistent scalar predictions");
		System.out.println("PASS NcrnaVoteDiagLogSinkTest optional_cells versioned_schema finite_present_values");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
