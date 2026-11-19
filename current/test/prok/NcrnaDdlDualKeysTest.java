package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import map.IntHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlDualKeys;
import prok.NcrnaDdlQueryKeys;
import shared.Shared;
import shared.Tools;

/** Independent span enumeration, native winning-key parity and separate-half
 * controls for every proposed two-length geometry. No genome accuracy claim.
 * @author Brian, Raiden
 */
public final class NcrnaDdlDualKeysTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlDualKeysTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			boolean run=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("run")){run=Parse.parseBoolean(b);}
				else{check(parser.parse(arg, a, b), "Unknown dual-length fixture option: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			check(run && assertions && parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null,
				"Use run=t and -ea; the fixed synthetic fixture takes no external input or output");
			DynamicDemiLog.setExponent(4);
			check(!DynamicDemiLog.blacklistExists(), "Native winning-key controls require no process-global blacklist");
			final byte[] random=new byte[1024];
			final Random rng=new Random(20261113);
			for(int i=0; i<random.length; i++){random[i]=(byte)"ACGTNacgtn".charAt(rng.nextInt(10));}
			final byte[][] cases={bytes(""), bytes("AAAAAA"), bytes("ACGTACG"), bytes("ACGTACGTACG"),
				bytes("AAAAAAAAANAAAAAAAAA"), bytes("GGGGGGGGGGGNNNTTTTTTTTTTT"), random,
				bytes("AAAAAAAAAAAAAAAAAAAARYSWKMBDHVryswkmbdhvUu-CCCCCCCCCCCCCCCCCCCC")};
			int checked=0;
			for(final int shortK:new int[]{7,8,9}){
				for(final int longK:new int[]{10,11}){
					for(final int buckets:new int[]{2048,4096}){
						for(final boolean partitioned:new boolean[]{false,true}){
							final NcrnaDdlDualKeys query=new NcrnaDdlDualKeys(shortK, longK, buckets, partitioned);
							for(final byte[] bases:cases){
								verify(query, bases, shortK, longK, buckets, partitioned);
								checked++;
								query.build(bytes("N"));
								check(query.keyCount()==0 && empty(query.memberKeys()), "Reuse must clear all prior query/member keys");
							}
						}
					}
				}
			}
			NcrnaDdlPartitionKeysTest.run();
			pp.outstream.println("NCRNA_DDL_DUAL_KEYS_TEST_PASS\tcases="+checked);
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void verify(final NcrnaDdlDualKeys query, final byte[] bases, final int shortK,
			final int longK, final int buckets, final boolean partitioned){
		query.build(bases);
		final IntHashSet expected=new IntHashSet(256);
		final DynamicDemiLog nativeSketch=DynamicDemiLog.create(buckets, longK, 12345, 0, false, false);
		final char[] nativeHalves=new char[buckets];
		oracle(bases, shortK, false, buckets, partitioned, true, expected, nativeSketch, nativeHalves);
		oracle(bases, longK, true, buckets, partitioned, false, expected, nativeSketch, nativeHalves);
		check(query.buckets()==buckets && query.keyCount()==expected.size(), "Every geometry must conserve all distinct raw keys");
		final IntHashSet observed=new IntHashSet(256);
		final int[] prefix=query.packedArray();
		for(int i=0; i<query.keyCount(); i++){
			final int key=prefix[i];
			check(key==query.packedKey(i), "Borrowed storage must preserve the indexed live-prefix values");
			check(observed.add(key) && expected.contains(key),
				"Every live query key must be unique and come from an independently enumerated defined span");
		}
		check(Arrays.equals(query.memberKeys(), partitioned ? nativeHalves : nativeSketch.maxArray()),
			"Member votes must use native winning keys after both lengths compete, or independently within each half");
		if(partitioned){
			final byte[] normalized=bases.clone();
			for(int i=0; i<normalized.length; i++){if(!defined(normalized[i])){normalized[i]='N';}}
			final NcrnaDdlQueryKeys shortOnly=new NcrnaDdlQueryKeys(buckets/2, shortK);
			final NcrnaDdlQueryKeys longOnly=new NcrnaDdlQueryKeys(buckets/2, longK);
			shortOnly.build(normalized, null, false, false);
			longOnly.build(normalized, null, false, true);
			check(query.keyCount()==shortOnly.keyCount()+longOnly.keyCount(), "Partition halves cannot collide across lengths");
			for(int i=0; i<shortOnly.keyCount(); i++){
				check(expected.contains(shortOnly.packedKey(i)), "Short half must equal the untagged half-width plain codec");
			}
			for(int i=0; i<longOnly.keyCount(); i++){
				check(expected.contains(longOnly.packedKey(i)+((buckets/2)<<16)),
					"Weight zero must retain the untagged half-width masked long codec exactly");
			}
		}
		scores(query, expected, buckets, partitioned);
	}

	/** Full scalar bucket walk is independent of the inverted posting traversal. */
	private static void scores(final NcrnaDdlDualKeys query, final IntHashSet expected,
			final int buckets, final boolean partitioned){
		final char[] keys=query.memberKeys().clone();
		final char[][] profiles={keys, new char[buckets], keys.clone()};
		for(int bucket=0; bucket<buckets; bucket+=2){profiles[2][bucket]=0;}
		final DDLIndexSmall index=new DDLIndexSmall(profiles);
		final int[] scores=new int[profiles.length];
		final int[] shortScores=new int[profiles.length], longScores=new int[profiles.length];
		final int[][] halves={shortScores, longScores};
		query.count(index, scores);
		if(partitioned){query.countHalves(index, halves);}
		for(int model=0; model<profiles.length; model++){
			int shortMatches=0;
			int longMatches=0;
			for(int bucket=0; bucket<buckets; bucket++){
				final int key=profiles[model][bucket];
				if(key==0 || !expected.contains((bucket<<16)|key)){continue;}
				if(bucket<buckets/2){shortMatches++;}
				else{longMatches++;}
			}
			check(scores[model]==shortMatches+longMatches, "Inverted dual scores must equal a full scalar bucket comparison");
			if(partitioned){
				check(shortScores[model]==shortMatches && longScores[model]==longMatches,
					"Saved half-scores must preserve both weight-zero and weight-one controls exactly");
			}
		}
		query.build(new byte[0]);
		query.count(index, scores);
		for(final int score:scores){check(score==0, "An empty replacement must clear previous model scores");}
		if(partitioned){
			query.countHalves(index, halves);
			for(int model=0; model<profiles.length; model++){
				check(shortScores[model]==0 && longScores[model]==0, "Both saved half-score arrays must clear on query replacement");
			}
		}
		boolean rejected=false;
		try{query.count(new DDLIndexSmall(new char[profiles.length][1]), scores);}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "A shared index cannot silently accept a different collector bucket width");
		rejected=false;
		try{query.countHalves(index, partitioned ? new int[][]{shortScores, shortScores} : halves);}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "Half counts require partitioned geometry and distinct output rows");
	}

	/** Independent slow reference: rebuild each span and reject any undefined base.
	 * Native add(), rather than our key maximum, supplies the member control.
	 */
	private static void oracle(final byte[] bases, final int k, final boolean middle, final int buckets,
			final boolean partitioned, final boolean shortWord, final IntHashSet expected,
			final DynamicDemiLog combined, final char[] nativeHalves){
		assert(k>0 && buckets>=2048) : "Fixture spans and combined widths must match a declared dual-length geometry";
		final int width=partitioned ? buckets/2 : buckets;
		final int offset=partitioned && !shortWord ? width : 0;
		final DynamicDemiLog half=partitioned ? DynamicDemiLog.create(width, k, 12345, 0, false, false) : null;
		for(int start=0; start+k<=bases.length; start++){
			long word=0;
			boolean defined=true;
			for(int i=0; i<k; i++){
				final byte b=bases[start+i];
				if(!defined(b)){defined=false; break;}
				word=(word<<2)|((b>>1)&3);
			}
			if(!defined){continue;}
			if(middle){word&=~(3L<<(((k-1)/2)*2));}
			final long tagged=word^(partitioned ? 0L : shortWord ? 0x4000000000000000L : 0x8000000000000000L);
			final long hash=Tools.hash64shift(tagged^12345L);
			final int bucket=(int)(hash&(width-1))+offset;
			expected.add((bucket<<16)|NcrnaDdlQueryKeys.encodeKey(hash));
			if(partitioned){half.add(word);}
			else{combined.add(tagged);}
		}
		if(partitioned){System.arraycopy(half.maxArray(), 0, nativeHalves, offset, width);}
	}

	private static boolean defined(final byte b){
		return b=='A' || b=='C' || b=='G' || b=='T' || b=='a' || b=='c' || b=='g' || b=='t';
	}

	private static boolean empty(final char[] keys){
		assert(keys!=null) : "A member-key control must supply its actual dense array";
		for(final char key:keys){if(key!=0){return false;}}
		return true;
	}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
