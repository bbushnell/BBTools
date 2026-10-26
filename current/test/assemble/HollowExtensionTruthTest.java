package test.assemble;

import java.nio.charset.StandardCharsets;

/** Independent tiny cases for origin offsets, immutable prefixes and suffix-only accuracy.
 * @author Brian Bushnell, Fischl */
public final class HollowExtensionTruthTest {

	/** Exercises scorer boundaries before it is trusted on a genome-scale output. */
	public static void main(final String[] args){
		final byte[] ref=bases("AACCGGTTAA");
		final byte[] seed=bases("CTGG");//One mutated seed base; source positions2..5 are CCGG.
		final HollowExtensionTruth.Score score=new HollowExtensionTruth.Score();
		HollowExtensionTruth.score(bases("CTGGTT"), seed, ref, 2, 2, score);
		check(score.added==2 && score.bad==0 && score.beyond==0 && score.unknown==0);
		HollowExtensionTruth.score(bases("CTGGTA"), seed, ref, 2, 2, score);
		check(score.added==2 && score.bad==1);
		HollowExtensionTruth.score(seed, seed, ref, 2, 2, score);
		check(score.added==0 && score.bad==0);
		HollowExtensionTruth.score(bases("CTGGTN"), seed, ref, 2, 2, score);
		check(score.unknown==1 && score.bad==0);
		HollowExtensionTruth.score(bases("TTAAAC"), bases("TTAA"), ref, 6, 2, score);
		check(score.added==2 && score.beyond==2 && score.bad==0);
		boolean failed=false;
		try{HollowExtensionTruth.score(bases("CCGGTT"), seed, ref, 2, 2, score);}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed);
		failed=false;
		try{HollowExtensionTruth.score(bases("CTGGTTA"), seed, ref, 2, 2, score);}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed);
		System.out.println("HOLLOW_EXTENSION_TRUTH_TEST_PASS checks=7");
	}

	private static byte[] bases(final String sequence){return sequence.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition){
		if(!condition){throw new AssertionError("Coordinate scoring fixture failed.");}
	}
}
