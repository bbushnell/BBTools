package test.assemble;

import java.util.Arrays;
import java.util.Random;

import assemble.HollowSeedCounts;
import dna.AminoAcid;
import ukmer.Kmer;

/** Checks seed-only multiplicities against independently rebuilt windows.
 * @author Brian Bushnell, Fischl */
public final class HollowSeedCountsTest {

	/** Exercises partial words, repeat multiplicity, Ns, RC and forced primary collisions. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(7189);
		for(int k : new int[]{67, 96, 124, 128, 160, 224}){
			final HollowSeedCounts counts=new HollowSeedCounts(k);
			for(int trial=0; trial<12; trial++){
				final byte[] seed=new byte[k+90];
				for(int i=0; i<seed.length; i++){
					seed[i]=AminoAcid.numberToBase[trial%3==0 ? 0 : random.nextInt(4)];
				}
				if(trial%3==2){seed[k/2]='N';}
				final int length=seed.length-7;
				counts.reset(seed, length);
				for(int start=0; start<=seed.length-k; start++){
					final Kmer query=window(seed, start, k);
					if(query.len<k){continue;}
					final int expected=oracle(seed, length, query);
					check(counts.count(query)==expected, "Rebuilt-window multiplicity k="+k+" start="+start);
					final byte[] rc=Arrays.copyOfRange(seed, start, start+k);
					AminoAcid.reverseComplementBasesInPlace(rc);
					check(counts.count(window(rc, 0, k))==expected, "Canonical reverse-complement multiplicity");
					final byte[] changed=Arrays.copyOfRange(seed, start, start+k);
					changed[32]=changed[32]=='A' ? (byte)'C' : (byte)'A';
					check(counts.count(window(changed, 0, k))==expected, "Ignored interior does not change seed support");
				}
				counts.reset(seed, k-1);
				check(counts.count(0, 0)==0, "Short reset must erase prior seed state");
			}
		}
		collisions();
		System.err.println("HOLLOW_SEED_COUNTS_TEST_PASS checks="+checks);
	}

	/** Builds a fresh word instead of sharing the ledger's rolling state. */
	private static Kmer window(final byte[] bases, final int start, final int k){
		assert(start>=0 && start+k<=bases.length) : "The independent oracle must not read outside its explicit window.";
		final Kmer word=new Kmer(k, true);
		for(int i=start; i<start+k; i++){word.addRight(bases[i]);}
		return word;
	}

	/** Enumerates only original-prefix windows; extension bases must never count as self. */
	private static int oracle(final byte[] seed, final int length, final Kmer query){
		assert(length<=seed.length && query.len>=query.kbig) : "Oracle counts only complete keys inside the original prefix.";
		int count=0;
		for(int start=0; start<=length-query.kbig; start++){
			final Kmer word=window(seed, start, query.kbig);
			if(word.len>=word.kbig && word.xor()==query.xor() && word.xor2()==query.xor2()){count++;}
		}
		return count;
	}

	/** Repeats a forced primary collision across clears and map/list growth. */
	private static void collisions(){
		final HollowSeedCounts counts=new HollowSeedCounts(128);
		for(int round=0; round<3; round++){
			counts.clear();
			for(long key : new long[]{0, -1, Long.MIN_VALUE, Long.MAX_VALUE}){
				counts.add(key, 0);
				counts.add(key, -1);
				counts.add(key, Long.MIN_VALUE);
				counts.add(key, -1);
				check(counts.count(key, 0)==1, "Collision chain head retained");
				check(counts.count(key, -1)==2, "Secondary match incremented, not primary-only match");
				check(counts.count(key, Long.MIN_VALUE)==1, "Collision chain tail retained");
				check(counts.count(key, 5)==0, "Absent secondary not conflated");
			}
			for(int i=1; i<=500; i++){counts.add(i, i+1);}
			for(int i=1; i<=500; i++){check(counts.count(i, i+1)==1, "Ledger survives map/list resize");}
		}
		counts.clear();
		check(counts.count(0, -1)==0, "Clear removes collision chains");
	}

	/** Keeps failures loud even when somebody invokes the test without -ea. */
	private static void check(final boolean pass, final String message){
		checks++;
		if(!pass){throw new AssertionError(message);}
	}

	private static long checks;
}
