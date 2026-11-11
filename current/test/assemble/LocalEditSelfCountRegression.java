package test.assemble;

import assemble.HomopolymerIndelProposal;
import assemble.LocalEditCorrector;

import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Synthetic count-table regressions through the real LocalEditCorrector.
 * Checks clean self-counted abstention and substitution/insertion/deletion repairs
 * in both orientations. Assigned fixture depths are not biological truth evidence.
 * Run in a fresh JVM: main changes shared packing and quality flags without restoring them.
 * @author Fischl */
public final class LocalEditSelfCountRegression{

	/** Runs two clean and six repair cases; ignores arguments and keeps explicit oracle checks active without -ea. */
	public static void main(String[] args){
		Kmer.PACKED=true; Kmer.MASK_CORE=false; Read.CHANGE_QUALITY=false;
		for(boolean reverse : new boolean[]{false, true}){
			clean(reverse);
			for(int error=0; error<3; error++){positive(error, reverse);}
		}
		System.out.println("SELF_COUNT_REGRESSION_PASS clean=2 true_error_repairs=6 actual_caller=true quality=true");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds one actual query observation to assigned counts, then requires abstention and
	 * original array identities. The preceding forward held-out diagnostic is printed in
	 * both orientation runs but its result is deliberately not an acceptance requirement. */
	private static void clean(boolean reverse){
		final byte[] original=sequence(23), alternate=original.clone();
		alternate[P]='T';
		final Counts counts=new Counts();
		for(int s=0; s+K<=original.length; s++){counts.assign(original, s, s>=P-K+1 && s<=P ? 2 : 29);}
		for(int s=P-K+1; s<=P; s++){counts.assign(alternate, s, 21);}
		final Read heldout=new Read(original.clone(), quality(original.length), "heldout", 0, false);
		final int diagnostic=new LocalEditCorrector(K, counts).correctOne(heldout, true);
		System.out.println("HELDOUT_DIAGNOSTIC reverse="+reverse+" applied="+diagnostic+"; not an acceptance requirement");
		final byte[] input=reverse ? AminoAcid.reverseComplementBases(original) : original.clone();
		counts.addRead(input, 1);
		for(int s=P-K+1; s<=P; s++){
			require(counts.depth(original, s)==3 && counts.depth(alternate, s)==21, "SELF must add actual original-query words without increasing alternate support.");
		}
		final byte[] q=quality(input.length);
		final Read read=new Read(input, q, "self", 0, false);
		require(new LocalEditCorrector(K, counts).correctOne(read, true)==0, "Actual corrector must abstain on the self-counted clean repeat-variant fixture.");
		require(read.bases==input && read.quality==q, "SELF abstention must preserve sequence/quality array identities, not merely equal copies.");
	}

	/** Builds one known error (0 substitution, 1 extra base, 2 missing base), adds its
	 * self-count to 25 truth copies, and requires one correction plus exact oriented
	 * bases and literal quality splicing. No replacement quality is inferred from depth. */
	private static void positive(int error, boolean reverse){
		final byte[] truth=sequence(101+error), observed;
		if(error==0){observed=truth.clone(); observed[P]='T';}else if(error==1){
			observed=new byte[truth.length+1];
			System.arraycopy(truth, 0, observed, 0, P); observed[P]='T';
			System.arraycopy(truth, P, observed, P+1, truth.length-P);
		}else{
			observed=new byte[truth.length-1];
			System.arraycopy(truth, 0, observed, 0, P);
			System.arraycopy(truth, P+1, observed, P, truth.length-P-1);
		}
		final Counts counts=new Counts(); counts.addRead(truth, 25); counts.addRead(observed, 1);
		final byte[] q=quality(observed.length), expectedQ=new byte[truth.length];
		if(error==0){
			System.arraycopy(q, 0, expectedQ, 0, q.length); expectedQ[P]=0;
		}else if(error==1){
			System.arraycopy(q, 0, expectedQ, 0, P); System.arraycopy(q, P+1, expectedQ, P, q.length-P-1);
		}else{
			System.arraycopy(q, 0, expectedQ, 0, P); expectedQ[P]=0; System.arraycopy(q, P, expectedQ, P+1, q.length-P);
		}
		final Read read=new Read(reverse ? AminoAcid.reverseComplementBases(observed) : observed,
			reverse ? reversed(q) : q, "positive", 0, false);
		require(new LocalEditCorrector(K, counts).correctOne(read, true)==1, "Actual corrector must repair the genuine single-base positive control; a blanket no-op must fail.");
		require(Arrays.equals(read.bases, reverse ? AminoAcid.reverseComplementBases(truth) : truth), "S/I/D repair must restore exact truth in both orientations.");
		require(Arrays.equals(read.quality, reverse ? reversed(expectedQ) : expectedQ), "Retained Phred qualities must splice exactly; only new/replaced bases receive Q0.");
		System.out.println("POSITIVE_PASS error="+error+" reverse="+reverse);
	}

	/** Makes a seeded 220-base context with fixed endpoints and an A-C-G mutation neighborhood. */
	private static byte[] sequence(int seed){
		final Random random=new Random(seed);
		final byte[] out=new byte[220];
		for(int i=0; i<out.length; i++){out[i]=BASES[random.nextInt(4)];}
		out[0]=out[out.length-1]='A'; out[P-1]='A'; out[P]='C'; out[P+1]='G';
		return out;
	}
	/** Builds a position-varying Phred10–39 oracle so misplaced qualities are observable. */
	private static byte[] quality(int length){
		final byte[] q=new byte[length];
		for(int i=0; i<length; i++){q[i]=(byte)(10+i%30);}
		return q;
	}
	/** Returns a reversed copy for orientation-specific quality expectations. */
	private static byte[] reversed(byte[] in){
		final byte[] out=new byte[in.length];
		for(int i=0; i<in.length; i++){out[in.length-1-i]=in[i];}
		return out;
	}
	/** Keeps the finite oracle active independently of Java assertion mode. */
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Serial fixture-owned canonical table and rolling scratch, borrowed by each corrector. */
	private static final class Counts implements HomopolymerIndelProposal.CountLookup{
		/** Reused exact-K query word and its one-dimensional count table. */
		final Kmer key=new Kmer(K); final HashArrayU1D table=new HashArrayU1D(new int[]{2003}, key.k, K);
		/** Rebuilds one complete word from the supplied sequence window. */
		void word(byte[] seq, int s){
			assert(s>=0 && s+K<=seq.length) : "Exact fixture lookup requires a complete K62 window.";
			key.clearFast();
			for(int i=s; i<s+K; i++){key.addRight(seq[i]);}
		}
		/** Assigns an initially absent word by increments; collisions invalidate the setup. */
		void assign(byte[] seq, int s, int n){
			word(seq, s);
			require(table.getValue(key)<1, "Unexpected canonical-word collision invalidates independent count assignments.");
			for(int i=0; i<n; i++){table.increment(key);}
		}
		/** Adds every overlapping window once per requested copy, including existing words. */
		void addRead(byte[] seq, int copies){
			for(int s=0; s+K<=seq.length; s++){
				word(seq, s);
				for(int i=0; i<copies; i++){table.increment(key);}
			}
		}
		/** Returns a window's accumulated depth with absence normalized to zero. */
		int depth(byte[] seq, int s){word(seq, s); return Math.max(0, table.getValue(key));}
		/** Implements the corrector lookup directly; the backend's absent sentinel is retained. */
		@Override
		public int count(Kmer word){return table.getValue(word);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact word length and zero-based mutation coordinate in each 220-base fixture. */
	private static final int K=62, P=100;
	/** Seeded sequence alphabet in its fixed random-draw order. */
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
}
