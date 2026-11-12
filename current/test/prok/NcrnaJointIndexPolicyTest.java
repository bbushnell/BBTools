package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.TrnaKmerIndex;
import shared.Shared;

/** Literal cutoff fixtures over known model counts, independent of any caller
 * library. At k1 the three models share10,5,2postings with the query's three
 * distinct bases. This makes cap, absolute cutoff and both relative fractions
 * observable separately and together, including ceiling and keep-best behavior.
 * Whole-caller fixed/adaptive equivalence remains a separate integration gate.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaJointIndexPolicyTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaJointIndexPolicyTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				require(parser.parse(arg, a, b), "Unknown test parameter: "+arg);
			}
			multipleModels();
			singleModel();
			pp.outstream.println("NCRNA_JOINT_INDEX_POLICY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void multipleModels(){
		final byte[][] models={bytes("AAAAACCCCC"), bytes("AAAAC"), bytes("AC")};
		final byte[] query=bytes("AAAAAACCCCCGGGGG");
		final TrnaKmerIndex index=new TrnaKmerIndex(models, 1, true, 1, 0, 0, 1);
		expect(index.shortlist(query, 3), 0, 1, 2);
		require(index.lastSharedCount(0)==10 && index.lastSharedCount(1)==5
			&& index.lastSharedCount(2)==2 && index.lastUniqueQueryKmers()==3,
			"The literal query must deduplicate its three base types while retaining model posting multiplicity");
		expect(index.shortlist(query, 2), 0, 1);
		index.setPolicy(true, 6, 0, 0, 1);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(true, 1, .5f, 0, 1);
		expect(index.shortlist(query, 3), 0, 1);
		index.setPolicy(true, 1, .5001f, 0, 1);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(true, 1, 0, 1.5f, 1);
		expect(index.shortlist(query, 3), 0, 1);
		index.setPolicy(true, 1, 0, 2, 1);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(true, 4, .25f, 2, 1);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(false, 4, 100, 100, 3);
		expect(index.shortlist(query, 3), 0, 1);
		index.setPolicy(true, 3, 0, 0, 3);
		expect(index.shortlist(query, 3), 0, 1);
		index.setPolicy(true, 100, 0, 0, 3);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(false, 100, 0, 0, 100);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(true, 1, .6f, 0, 1);
		boolean rejected=false;
		try{index.setPolicy(false, 0, Float.NaN, 0, 0);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Nonfinite fractions must fail before mutating the index policy");
		index.setAdaptiveFloor(1);
		expect(index.shortlist(query, 3), 0);
		index.setPolicy(true, 1, 0, 0, 1);
		expect(index.shortlist(query, 3), 0, 1, 2);
	}

	/** Fixed0's no-hit bypass is deliberately not an adaptive-equivalence claim. */
	private static void singleModel(){
		final TrnaKmerIndex index=new TrnaKmerIndex(new byte[][]{bytes("AAAAACCCCC")}, 1, false, 0, 0, 0, 1);
		final byte[] hit=bytes("AC"), miss=bytes("GG");
		expect(index.shortlist(hit, 1), 0);
		expect(index.shortlist(miss, 1));
		index.setPolicy(true, 1, 0, 0, 1);
		expect(index.shortlist(hit, 1), 0);
		expect(index.shortlist(miss, 1));
		index.setPolicy(false, 0, 0, 0, 0);
		expect(index.shortlist(miss, 1), 0);
		index.setPolicy(true, 0, 0, 0, 0);
		expect(index.shortlist(miss, 1));
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void expect(final int[] actual, final int... expected){
		require(Arrays.equals(actual, expected), "Expected model IDs "+Arrays.toString(expected)+", observed "+Arrays.toString(actual));
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
