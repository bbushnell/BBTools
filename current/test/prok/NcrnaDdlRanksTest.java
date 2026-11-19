package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlQueryKeys;
import prok.NcrnaDdlRanks;
import shared.Shared;

/** Literal rank/tie oracles and a whole-window native-score composition check.
 * No genome accuracy, allocation-rate or performance result.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlRanksTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlRanksTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown rank fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Rank fixtures require assertions and own their literal inputs without files or sampling");
			DynamicDemiLog.setExponent(4);
			literalRanks();
			windowScores();
			maskedWords();
			canonicalParity();
			pp.outstream.println("NCRNA_DDL_RANKS_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Literal Oracles      ----------------*/
	/*--------------------------------------------------------------*/

	private static void literalRanks(){
		final float[] scores={.5f, .25f, .5f, 0f, .5f};
		final boolean[] acceptable={false, true, true, false, true};
		final int[] result={99, 99, 99, 99, 73};
		NcrnaDdlRanks.best(scores, acceptable, result);
		check(Arrays.equals(result, new int[]{2, 1, 3, 2, 73}),
			"Model2 is the first acceptable .5 tie; the full tie covers ranks1..3 and output tail is untouched");
		acceptable[0]=true;
		NcrnaDdlRanks.best(scores, acceptable, result);
		check(Arrays.equals(result, new int[]{1, 1, 3, 0, 73}), "Smaller acceptable model index wins an equal score");
		Arrays.fill(acceptable, false);
		acceptable[1]=true;
		NcrnaDdlRanks.best(scores, acceptable, result);
		check(Arrays.equals(result, new int[]{4, 4, 4, 1, 73}), "All three .5 models precede the single .25 model");
		acceptable[1]=false;
		NcrnaDdlRanks.best(scores, acceptable, result);
		check(Arrays.equals(result, new int[]{-1, -1, -1, -1, 73}), "No oracle success must remain absent, never rank1");
		acceptable[4]=true;
		NcrnaDdlRanks.best(new float[5], acceptable, result);
		check(Arrays.equals(result, new int[]{5, 1, 5, 4, 73}), "All-zero scores retain deterministic rank and complete ambiguity");
		for(final float invalid:new float[]{Float.NaN, Float.POSITIVE_INFINITY, -1f}){
			try{NcrnaDdlRanks.best(new float[]{invalid}, new boolean[]{true}, result);}
			catch(IllegalArgumentException expected){continue;}
			throw new AssertionError("Invalid score was accepted: "+invalid);
		}
	}

	/** A homopolymer supplies one canonical key. Two profiles match it, even when
	 * one has another nonmatching bucket; native WKID's minimum divisor is six.
	 */
	private static void windowScores(){
		final byte[] sequence="AAAAAAAAAAAAAAAAAAAA".getBytes(StandardCharsets.US_ASCII);
		final byte[] saved=sequence.clone();
		final NcrnaDdlQueryKeys raw=new NcrnaDdlQueryKeys(32, 5);
		raw.build(sequence);
		final int[] packed=raw.packedKeys();
		check(packed.length==1, "The separately tested raw codec supplies the homopolymer's one canonical key");
		final int bucket=packed[0]>>>16, key=packed[0]&65535;
		final char[][] keys=new char[4][32];
		keys[0][bucket]=keys[1][bucket]=(char)key;
		keys[1][(bucket+1)&31]=1;
		keys[2][bucket]=(char)(key==65535 ? 65534 : key+1);
		final NcrnaDdlProfile[] profiles=new NcrnaDdlProfile[4];
		for(int i=0; i<profiles.length; i++){
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(32);
			votes.add(keys[i]);
			profiles[i]=votes.finish();
		}
		final NcrnaDdlRanks ranks=new NcrnaDdlRanks(profiles, 5);
		final boolean[] acceptable={false, true, false, false};
		final int[] result=new int[4];
		try{ranks.best(NcrnaDdlRanks.SSU, acceptable, result);}
		catch(IllegalStateException expected){result[0]=17;}
		check(result[0]==17, "Unscored worker scratch cannot supply plausible ranks");
		for(int round=0; round<2; round++){
			ranks.score(sequence);
			for(int method=0; method<NcrnaDdlRanks.METHODS; method++){
				ranks.best(method, acceptable, result);
				check(Arrays.equals(result, new int[]{2, 1, 2, 1}), "Only models0/1 share the one nonzero query key");
				for(int model=0; model<profiles.length; model++){
					final float expected=model>1 ? 0 : method==NcrnaDdlRanks.SSU ? 1f/6f : 1f;
					check(Float.floatToIntBits(ranks.score(method, model))==Float.floatToIntBits(expected),
						"Exact native float score or one-bucket count differs: method="+method+", model="+model);
				}
			}
			ranks.score("NNNN".getBytes(StandardCharsets.US_ASCII));
			for(int method=0; method<NcrnaDdlRanks.METHODS; method++){
				ranks.best(method, acceptable, result);
				check(Arrays.equals(result, new int[]{2, 1, 4, 1}), "Replacing a query with no complete keys clears every previous score");
			}
		}
		check(Arrays.equals(sequence, saved), "Whole-window ranking borrows input bases read-only");
	}

	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}

	/** One central SNP in each complete k-mer span must disappear only in masked mode. */
	private static void maskedWords(){
		for(int k=7; k<=11; k++){
			final byte[] member=new byte[k];
			Arrays.fill(member, (byte)'A');
			final byte[] query=member.clone();
			query[k/2]='C';
			for(final int width:new int[]{256, 512, 1024, 2048}){
				for(final boolean middle:new boolean[]{false, true}){
					final DynamicDemiLog sketch=DynamicDemiLog.create(width, k, 12345, 0, false, false);
					new NcrnaDdlQueryKeys(width, k).build(member, sketch, false, middle);
					final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(width);
					votes.add(sketch.maxArray());
					final NcrnaDdlProfile[] profiles={votes.finish()};
					final NcrnaDdlRanks ranks=new NcrnaDdlRanks(profiles, k, false, middle);
					ranks.score(query);
					check(ranks.score(NcrnaDdlRanks.ALL_KEYS, 0)==(middle ? 1f : 0f)
						&& ranks.score(NcrnaDdlRanks.SSU, 0)==(middle ? 1f/6f : 0f),
						"A sole middle SNP scores zero plain and positive masked, for members and queries: k="+k+", width="+width);
					if(middle && (k&1)==0){
						final byte[] left=member.clone();
						left[k/2-1]='C';
						ranks.score(left);
						check(ranks.score(NcrnaDdlRanks.ALL_KEYS, 0)==0,
							"Even-k masking must clear the right central base, retaining the distinct left central SNP");
					}
					ranks.score(new byte[0]);
					check(ranks.score(NcrnaDdlRanks.ALL_KEYS, 0)==0 && ranks.score(NcrnaDdlRanks.SSU, 0)==0,
						"Replacing a forward/masked query clears both raw keys and native sketch");
				}
			}
			check(member[k/2]=='A' && query[k/2]=='C', "Word encoding must not mutate borrowed member/query bases");
		}
		final byte[] a="AAAAAAAAA".getBytes(StandardCharsets.US_ASCII), t="TTTTTTTTT".getBytes(StandardCharsets.US_ASCII);
		final DynamicDemiLog first=DynamicDemiLog.create(256, 9, 12345, 0, false);
		final DynamicDemiLog second=DynamicDemiLog.create(256, 9, 12345, 0, false);
		new NcrnaDdlQueryKeys(256, 9).build(a, first, false, false);
		new NcrnaDdlQueryKeys(256, 9).build(t, second, false, false);
		check(!Arrays.equals(first.maxArray(), second.maxArray()), "Forward-only words must retain the distinction between reverse complements");
	}

	/** The unchanged native loop independently controls the shared canonical codec. */
	private static void canonicalParity(){
		final byte[] sequence="ACGTCGTAGCTAGTGCATGACGTAGCGTACGATCGATCGNNNacgtcgatgctagcgtagctacgttgca".getBytes(StandardCharsets.US_ASCII);
		for(int k=7; k<=11; k++){
			for(final int width:new int[]{256, 512, 1024, 2048}){
				final DynamicDemiLog nativeSketch=DynamicDemiLog.create(width, k, 12345, 0, false);
				final DynamicDemiLog encoded=DynamicDemiLog.create(width, k, 12345, 0, false);
				nativeSketch.hash(sequence, null);
				final NcrnaDdlQueryKeys raw=new NcrnaDdlQueryKeys(width, k);
				raw.build(sequence, encoded, true, false);
				check(Arrays.equals(nativeSketch.maxArray(), encoded.maxArray()),
					"Canonical experiment words must exactly match the native sketch with ambiguity and lower-case bases");
			}
		}
		boolean rejected=false;
		try{new NcrnaDdlQueryKeys(256, 8).build(sequence, DynamicDemiLog.create(256, 8, 12345, 0, false), true, true);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "The new middle mask must not silently canonicalize even-k words");
	}
}
