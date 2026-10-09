package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import structures.ByteBuilder;
import tracker.EntropyTracker;

/** Native sequence-complexity inputs for the two original, untrimmed fusion tips. @author Fischl */
public final class FusionTipEntropy {

	/** Measures up to100 terminal bases, averaging only complete, defined50-base windows. */
	public float mean(final byte[] bases, final boolean rightEnd){
		assert(bases!=null) : "Tip entropy requires the original oriented contig sequence.";
		final int length=Math.min(TIP_BASES, bases.length);
		if(length<WINDOW_BASES){return UNAVAILABLE;}
		final int from=rightEnd ? bases.length-length : 0, to=from+length-1;
		int run=0;
		boolean available=false;
		for(int i=from; i<=to; i++){
			run=AminoAcid.isFullyDefined(bases[i]) ? run+1 : 0;
			if(run>=WINDOW_BASES){available=true; break;}
		}
		// EntropyTracker returns zero if every window has Ns; distinguish that
		// absence from a measured homopolymer, whose entropy is legitimately zero.
		if(!available){return UNAVAILABLE;}
		final float value=tracker.averageEntropy(bases, false, from, to);
		if(!Float.isFinite(value) || value<0 || value>1){
			throw new IllegalStateException("Native tip entropy is outside [0,1]: "+value);
		}
		return value;
	}

	/** Appends symmetric min/max inputs and returns the number of unavailable physical tips. */
	public int append(final ByteBuilder row, final byte[] source, final byte[] dest){
		assert(row!=null) : "Entropy values must extend an existing candidate record.";
		final float left=mean(source, true), right=mean(dest, false);
		appendValue(row, Math.min(left, right));
		appendValue(row, Math.max(left, right));
		return (left<0 ? 1 : 0)+(right<0 ? 1 : 0);
	}

	/** Writes the sentinel literally; measured values retain the vector format's nine decimals. */
	private static void appendValue(final ByteBuilder row, final float value){
		row.tab();
		if(value<0){row.append("-0.1");}else{row.append(value, 9);}
	}

	/** Tests exact endpoints, full sliding averages, missing data and reciprocal orientations. */
	static void selfTest(){
		final FusionTipEntropy entropy=new FusionTipEntropy();
		final byte[] poly=new byte[130], repeat=new byte[100], random=new byte[180];
		Arrays.fill(poly, (byte)'A');
		final Random rng=new Random(195);
		for(int i=0; i<repeat.length; i++){repeat[i]=(byte)"AC".charAt(i%2);}
		for(int i=0; i<random.length; i++){random[i]=(byte)"ACGT".charAt(rng.nextInt(4));}
		check(entropy.mean(poly, true)==0, "Measured homopolymers must not be unavailable.");
		check(entropy.mean(Arrays.copyOf(poly, 49), true)==UNAVAILABLE, "Short tips need the sentinel.");
		check(entropy.mean(Arrays.copyOf(poly, 50), true)==0, "One complete window must be measured.");
		final byte[] unknown=new byte[100];
		Arrays.fill(unknown, (byte)'N');
		check(entropy.mean(unknown, false)==UNAVAILABLE, "N-only tips must not become zero entropy.");
		System.arraycopy(poly, 0, unknown, 0, 50);
		check(entropy.mean(unknown, false)==0, "A valid window beside Ns must remain measurable.");
		check(entropy.mean(repeat, true)<entropy.mean(random, true), "Short repeats should have lower entropy.");
		final EntropyTracker direct=new EntropyTracker(ENTROPY_K, WINDOW_BASES, false);
		double sum=0;
		for(int i=random.length-TIP_BASES; i<=random.length-WINDOW_BASES; i++){
			sum+=direct.averageEntropy(random, false, i, i+WINDOW_BASES-1);
		}
		check(Math.abs(entropy.mean(random, true)-sum/51)<1e-6, "The terminal100bp mean must contain51 full windows.");
		final byte[] reversed=random.clone();
		AminoAcid.reverseComplementBasesInPlace(reversed);
		check(Math.abs(entropy.mean(random, true)-entropy.mean(reversed, false))<1e-6,
				"Reverse complement changed tip entropy.");
		final ByteBuilder forward=new ByteBuilder(), backward=new ByteBuilder();
		entropy.append(forward, random, repeat);
		final byte[] repeatReverse=repeat.clone();
		AminoAcid.reverseComplementBasesInPlace(repeatReverse);
		entropy.append(backward, repeatReverse, reversed);
		check(forward.toString().equals(backward.toString()), "Exchanging oriented join ends changed the inputs.");
		check(entropy.mean("ACGT".getBytes(StandardCharsets.US_ASCII), true)==UNAVAILABLE,
				"A tip shorter than the entropy K cannot supply a complete window.");
		System.err.println("FUSION_TIP_ENTROPY_TEST_PASS native_mean RC symmetry sentinel");
	}

	/** Keeps test failures loud even if assertions were disabled by a caller. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private final EntropyTracker tracker=new EntropyTracker(ENTROPY_K, WINDOW_BASES, false);
	public static final String VERSION="fusion_join_v2";
	public static final String COLUMNS="tip_entropy_min\ttip_entropy_max";
	public static final float UNAVAILABLE=-0.1f;
	public static final int ENTROPY_K=5, WINDOW_BASES=50, TIP_BASES=100;
}
