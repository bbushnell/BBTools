package synth;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import shared.Random;
import stream.Read;
import structures.ByteBuilder;

/** Deterministic regression for source-run context in the long-read error model.
 * @author Fischl
 */
public final class RandomReadsMGTest {

	public static void main(final String[] args){
		boolean assertions=false; assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run this regression with -ea.");}
		final RandomReadsMG generator=new RandomReadsMG(new String[]{"t=1","tree=f"});
		int literal=0, parity=0;
		for(final char b : new char[]{'A','C','G','T'}){
			final String run=""+b+b+b+b+b;
			check(generator,run,""+b+b+b+b+b+b+b+b+b,0,1,4); literal++;
			check(generator,run,""+b,0.01f,1,4); literal++;
			check(generator,run,""+b+b+b+b+b+b+b,0,0.3f,2); literal++;
			check(generator,""+b+b+'N'+b,""+b+b+b+'N'+b,0,1,1); literal++;
			final char other=b=='A' ? 'C' : 'A';
			check(generator,""+b+b+other+b+b,""+b+b+b+other+b+b+b,0,1,2); literal++;
		}
		// Independent archived hrate=0 semantics must retain exact bytes, event
		// totals, and subsequent RNG state, including ambiguous input bases.
		for(int seed=0; seed<128; seed++){
			final java.util.Random bases=new java.util.Random(seed);
			final byte[] input=new byte[257]; final byte[] alphabet=ascii("ACGTN");
			for(int i=0; i<input.length; i++){input[i]=alphabet[bases.nextInt(alphabet.length)];}
			final SeededRandom a=new SeededRandom(seed), b=new SeededRandom(seed);
			final Read actual=new Read(input.clone(),null,"parity",0), expected=new Read(input.clone(),null,"parity",0);
			final int changes=generator.mutateLongRead(actual,0.07f,0.08f,0.09f,0,a);
			final int oldChanges=legacy(expected,0.07f,0.08f,0.09f,b);
			assert(changes==oldChanges && Arrays.equals(actual.bases,expected.bases) && a.nextLong()==b.nextLong()) :
				"hrate=0 must preserve old output and random draw consumption; seed="+seed;
			parity++;
		}
		System.out.println("RANDOMREADSMG_HP_TEST_OK literal="+literal+" zero_boost_parity="+parity);
	}

	private static void check(final RandomReadsMG generator, final String input, final String expected,
		final float dRate, final float hRate, final int expectedChanges){
		assert(input.length()>0) : "Literal fixture must exercise a source run.";
		final Read r=new Read(ascii(input),null,"literal",0);
		final int changes=generator.mutateLongRead(r,0,0,dRate,hRate,new FixedRandom());
		assert(changes==expectedChanges && Arrays.equals(r.bases,ascii(expected))) :
			"Source-run bonus/reset regression: input="+input+" expected="+expected+" actual="+
			new String(r.bases,StandardCharsets.US_ASCII)+" changes="+changes+" expectedChanges="+expectedChanges;
	}

	/** Old model with its always-zero bonus, retained only as a parity oracle. */
	private static int legacy(final Read r, final float s, final float ins, final float del, final Random rng){
		assert(r.bases!=null) : "Parity oracle requires source bases.";
		final float error=s+ins+del, delProb=del/Math.max(0.000000000001f,ins+del);
		final ByteBuilder out=new ByteBuilder(); int changes=0;
		for(final byte base : r.bases){
			if(!AminoAcid.isFullyDefined(base)){out.append(base); continue;}
			final float f=rng.nextFloat();
			if(f>=error){out.append(base);}
			else if(f<s){out.append(AminoAcid.numberToBase[(AminoAcid.baseToNumber[base]+rng.nextInt3()+1)&3]); changes++;}
			else{
				if(rng.nextFloat()>=delProb){out.append(base); out.append(AminoAcid.numberToBase[rng.nextInt()&3]);}
				changes++;
			}
		}
		r.bases=out.toBytes(); return changes;
	}

	private static byte[] ascii(final String s){assert(s!=null) : "Literal fixture text is required."; return s.getBytes(StandardCharsets.US_ASCII);}

	private static final class FixedRandom implements Random {
		@Override public float nextFloat(){return 0.5f;}
		@Override public long nextLong(){throw new AssertionError("Literal fixture must only request float draws.");}
		@Override public long nextLong(long bound){throw new AssertionError("Unexpected bounded draw.");}
		@Override public double nextGaussian(){throw new AssertionError("Unexpected Gaussian draw.");}
		@Override public void setSeed(long seed){throw new AssertionError("Unexpected reseed.");}
	}
	private static final class SeededRandom implements Random {
		SeededRandom(final long seed){rng=new java.util.Random(seed);}
		@Override public long nextLong(){return rng.nextLong();}
		@Override public long nextLong(long bound){throw new AssertionError("Unexpected bounded draw.");}
		@Override public double nextGaussian(){return rng.nextGaussian();}
		@Override public void setSeed(long seed){rng.setSeed(seed);}
		final java.util.Random rng;
	}
}
