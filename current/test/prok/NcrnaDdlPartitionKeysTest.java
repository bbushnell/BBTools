package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

import ddl.DDLIndexSmall;
import map.IntHashMap;
import map.IntHashSet;
import prok.NcrnaDdlPartitionKeys;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlQueryKeys;
import shared.Tools;

/** Independent complete-span enumeration and deliberate cross-length collisions.
 * Invoked by the existing dual-key fixture; no input files or biological claims.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaDdlPartitionKeysTest{

	static void run(){
		check(NcrnaDdlPartitionKeys.familyBucket(0, 1536)==0, "Unsigned range starts at zero");
		check(NcrnaDdlPartitionKeys.familyBucket(0x80000000L, 1536)==768, "Unsigned midpoint must reach the middle family bin");
		check(NcrnaDdlPartitionKeys.familyBucket(-1L, 1536)==1535, "Unsigned multiplication must not overflow an int");
		final byte[] random=new byte[1024];
		final Random rng=new Random(20261114);
		for(int i=0; i<random.length; i++){random[i]=(byte)"ACGTacgtNN".charAt(rng.nextInt(10));}
		final byte[][] cases={bytes(""), bytes("ACGTACG"), bytes("ACGTACGTAC"),
			bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("AAAAANAAAA"),
			bytes("AAAANAAAAARYSWKMBDHVUu-ACGTACGTACGT"), random};
		int checked=0;
		for(final int shortK:new int[]{0,7,8}){
			for(final int bins:new int[]{512,1024,2048}){
				final NcrnaDdlPartitionKeys codec=new NcrnaDdlPartitionKeys(shortK, bins);
				for(final byte[] sequence:cases){
					verify(codec, sequence, shortK, bins);
					checked++;
				}
			}
		}
		crossLengthCollision();
		votes();
		System.err.println("NCRNA_DDL_PARTITION_KEYS_TEST_PASS\tcases="+checked);
	}

	private static void verify(final NcrnaDdlPartitionKeys codec, final byte[] sequence, final int shortK, final int bins){
		final IntHashSet expected=new IntHashSet(256);
		final char[] member=new char[bins], consensus=new char[bins];
		if(shortK!=0){spans(sequence, shortK, true, shortK, bins, expected, member, consensus);}
		spans(sequence, 10, false, shortK, bins, expected, member, consensus);
		check(Arrays.equals(codec.memberKeys(sequence), member), "Member maxima must equal independent full-span enumeration");
		check(Arrays.equals(codec.consensusKeys(sequence), consensus), "Consensus maxima must stay entirely in the final quarter");
		codec.build(sequence);
		check(codec.keyCount()==expected.size(), "Every member/consensus lookup pair must be retained once");
		final IntHashSet observed=new IntHashSet(256);
		final int[] prefix=codec.packedArray();
		for(int i=0; i<codec.keyCount(); i++){
			final int pair=prefix[i];
			check(pair==codec.packedKey(i), "Borrowed storage must preserve the indexed live-prefix values");
			check(observed.add(pair) && expected.contains(pair), "Query keys must be unique complete bucket/key pairs");
		}
		final char[][] profiles={member, consensus, new char[bins]};
		final DDLIndexSmall index=new DDLIndexSmall(profiles);
		final int[] scores=new int[3];
		codec.count(index, scores);
		for(int model=0; model<profiles.length; model++){
			int count=0;
			for(int bucket=0; bucket<bins; bucket++){
				final int key=profiles[model][bucket];
				if(key!=0 && expected.contains((bucket<<16)|key)){count++;}
			}
			check(scores[model]==count, "Inverted scores must match a scalar walk for member-only, consensus-only and empty profiles");
		}
		codec.build(bytes("N"));
		codec.count(index, scores);
		check(codec.keyCount()==0 && scores[0]==0 && scores[1]==0 && scores[2]==0,
			"Query replacement must clear prior pairs and every model score");
		check(empty(codec.memberKeys(bytes("N"))) && empty(codec.consensusKeys(bytes("N"))),
			"Construction reuse must not retain a key across an undefined or empty span");
	}

	/** Slow reference rebuilds each span, including its masked base's validity. */
	private static void spans(final byte[] sequence, final int k, final boolean shortWord, final int shortK,
			final int bins, final IntHashSet expected, final char[] member, final char[] consensus){
		final int quarter=bins/4;
		for(int start=0; start+k<=sequence.length; start++){
			long word=0;
			boolean valid=true;
			for(int offset=0; offset<k; offset++){
				final int base=sequence[start+offset]&~32;
				if(base!='A' && base!='C' && base!='G' && base!='T'){
					valid=false;
					break;
				}
				word=(word<<2)|((base>>1)&3);
			}
			if(!valid){continue;}
			if(!shortWord){word&=~(3L<<8);}
			final long hash=Tools.hash64shift(word^12345L);
			final int bucket=shortK==0 ? (int)((hash&0xffffffffL)*(3*quarter)/0x100000000L)
				: shortWord ? (int)(hash&(quarter-1)) : quarter+(int)(hash&(2*quarter-1));
			put(bucket, hash, expected, member);
			final long consensusHash=shortK==0 ? hash
				: Tools.hash64shift(word^(shortWord ? 0x4000000000000000L : 0x8000000000000000L)^12345L);
			put(3*quarter+(int)(consensusHash&(quarter-1)), consensusHash, expected, consensus);
		}
	}

	private static void put(final int bucket, final long hash, final IntHashSet expected, final char[] maxima){
		check(bucket>=0 && bucket<maxima.length, "Independent reference must stay within the declared profile width");
		final char key=(char)NcrnaDdlQueryKeys.encodeKey(hash);
		expected.add((bucket<<16)|key);
		if(key>maxima[bucket]){maxima[bucket]=key;}
	}

	/** Find an actual compressed collision between differently tagged words.
	 * The finite search is fixture-only; the production loop never searches. */
	private static void crossLengthCollision(){
		final IntHashMap shortWords=new IntHashMap(16384);
		for(int word=0; word<(1<<14); word++){
			final long hash=Tools.hash64shift(word^0x4000000000000000L^12345L);
			shortWords.put(((int)(hash&127)<<16)|NcrnaDdlQueryKeys.encodeKey(hash), word);
		}
		for(int raw=0; raw<(1<<18); raw++){
			final long word=((long)(raw>>>8)<<10)|(raw&255);
			final long hash=Tools.hash64shift(word^0x8000000000000000L^12345L);
			final int pair=((int)(hash&127)<<16)|NcrnaDdlQueryKeys.encodeKey(hash);
			if(!shortWords.contains(pair)){continue;}
			final byte[] sequence=bytes(dna(shortWords.get(pair), 7)+"N"+dna(word, 10));
			final NcrnaDdlPartitionKeys codec=new NcrnaDdlPartitionKeys(7, 512);
			verify(codec, sequence, 7, 512);
			codec.build(sequence);
			final char[] profile=new char[512];
			profile[384+(pair>>>16)]=(char)pair;
			final int[] score=new int[1];
			codec.count(new DDLIndexSmall(new char[][]{profile}), score);
			check(score[0]==1, "A compressed shared-bin collision across lengths must score once, not twice");
			return;
		}
		throw new AssertionError("The bounded fixture must contain an actual short/long consensus collision");
	}

	private static void votes(){
		final char[] member=new char[512], consensus=new char[512];
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(512);
		member[3]=10;
		member[5]=40;
		votes.add(member);
		member[3]=20;
		votes.add(member);
		consensus[511]=65000;
		final NcrnaDdlProfile profile=votes.finishPartition(consensus);
		check(profile.members==2 && profile.keys()[3]==20 && profile.keys()[5]==40,
			"Real member plurality and highest-key ties must remain independent of consensus");
		check(profile.keys()[511]==65000 && profile.members==2,
			"The consensus quarter stores a maximum key without additional construction members");
		consensus[1]=3;
		try{votes.finishPartition(consensus);}
		catch(IllegalArgumentException expected){
			consensus[1]=0;
			member[511]=9;
			votes.add(member);
			try{votes.finishPartition(consensus);}
			catch(IllegalArgumentException alsoExpected){return;}
			throw new AssertionError("Member keys in the consensus quarter must fail");
		}
		throw new AssertionError("Consensus keys in a member partition must fail");
	}

	private static String dna(final long word, final int k){
		final char[] bases=new char[k];
		for(int i=0; i<k; i++){bases[k-1-i]="ACTG".charAt((int)((word>>>(2*i))&3));}
		return new String(bases);
	}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static boolean empty(final char[] keys){
		for(final char key:keys){if(key!=0){return false;}}
		return true;
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private NcrnaDdlPartitionKeysTest(){}
}
