package test.ukmer;

import ukmer.Kmer;

import java.util.Arrays;
import java.util.Random;

/**
 * Rebuilding oracle for fixed-window insertions, deletions and restoration.
 * Covers packed/symmetric geometry, core masking and hash mixing on seeded and
 * patterned sequences. Expected windows are edited as independent numeric-base
 * arrays and rebuilt; lazy caches are warmed before in-place edits are compared.
 * Includes reverse complements, rolling lengths greater than K and invalid edits.
 * Requires assertions and a fresh JVM; checks is cumulative. Packing, core-mask
 * and mixing settings are restored even when a check throws.
 * @author Fischl
 */
public final class KmerIndelTest{

	/** Runs the deterministic layout/edit grid and invalid-input cases.
	 * @param args Unused command-line arguments
	 * @throws IllegalStateException If assertions are disabled */
	public static void main(String[] args){
		boolean enabled=false; assert(enabled=true);
		if(!enabled){throw new IllegalStateException("Run KmerIndelTest with -ea.");}
		final boolean packed=Kmer.PACKED, core=Kmer.MASK_CORE, mix=Kmer.FULL_MIX;
		try{
			for(boolean packing : new boolean[]{false, true}){
				for(boolean mask : new boolean[]{false, true}){
					for(boolean fullMix : new boolean[]{false, true}){
						Kmer.PACKED=packing; Kmer.MASK_CORE=mask; Kmer.FULL_MIX=fullMix;
						for(int requested : new int[]{1, 2, 5, 31, 32, 33, 62, 63, 64, 65, 94, 95, 96, 127}){
							final int k=new Kmer(requested).kbig;
							if(requested==94){
								check(k==(packing ? 94 : 93), "Requested94 geometry", requested, 0, 0);
								System.out.println("K94_GEOMETRY packed="+packing+" effectiveK="+k+" core="+mask+" fullMix="+fullMix);
							}
							final Random random=new Random(710013+k);
							for(int pattern=0; pattern<4; pattern++){
								final byte[] seq=new byte[k];
								for(int i=0; i<k; i++){seq[i]=(byte)(pattern==0 ? random.nextInt(4) : pattern==1 ? 3 : pattern==2 ? i%4 : (i<k/2 ? 0 : 3));}
								run(requested, seq);
							}
						}
					}
				}
			}
			invalid();
			System.out.println("KMER_INDEL_TEST_OK checks="+checks);
		}finally{Kmer.PACKED=packed; Kmer.MASK_CORE=core; Kmer.FULL_MIX=mix;}
	}

	/** Checks every position/incoming base against independent fixed-length windows.
	 * Insertions discard the last base; deletions append an incoming base. Also checks
	 * edit return values, restoration, reverse complements and preserved rolling length.
	 * @param requested Requested K, which symmetric geometry may round down
	 * @param seq Numeric bases 0-3, with length equal to the effective K */
	private static void run(final int requested, final byte[] seq){
		assert(seq.length==new Kmer(requested).kbig) : "Oracle sequence must use the effective packed/legacy K, including legacy rounding.";
		final int k=seq.length;
		final Kmer original=build(requested, seq);
		final byte[] reverse=new byte[k];
		for(int i=0; i<k; i++){reverse[i]=(byte)(3-seq[k-1-i]);}
		for(int p=0; p<k; p++){
			for(int x=0; x<4; x++){
				final Kmer key=new Kmer(original); warm(key);
				final byte[] deleted=new byte[k];
				for(int i=0; i<k-1; i++){deleted[i]=seq[i<p ? i : i+1];}
				deleted[k-1]=(byte)x;
				check(key.deleteBase(p, x)==seq[p], "Deleted-base return", requested, p, x);
				equal(key, build(requested, deleted), k, requested, p, x);
				check(key.insertBase(p, seq[p])==x, "Deletion restoration drops incoming base", requested, p, x);
				equal(key, original, k, requested, p, x);
				final byte[] inserted=new byte[k];
				for(int i=0; i<k; i++){inserted[i]=i==p ? (byte)x : seq[i<p ? i : i-1];}
				check(key.insertBase(p, x)==seq[k-1], "Insertion returns dropped final base", requested, p, x);
				equal(key, build(requested, inserted), k, requested, p, x);
				// Reuse the inserted window for all four alternatives, with warm caches.
				for(int y=0; y<4; y++){
					key.substituteBase(p, y); inserted[p]=(byte)y;
					equal(key, build(requested, inserted), k, requested, p, y);
				}
				check(key.deleteBase(p, seq[k-1])==3, "Insertion restoration returns final tested allele", requested, p, x);
				equal(key, original, k, requested, p, x);
				// Rolling into a full window may have len>K. Edits must not reset len.
				final byte[] rolled=seq.clone();
				System.arraycopy(seq, 1, rolled, 0, k-1); rolled[k-1]=(byte)x;
				key.addRightNumeric(x); warm(key);
				final long removed=key.deleteBase(p, 2); key.insertBase(p, removed);
				equal(key, build(requested, rolled), k+1, requested, p, x);
				final Kmer rc=new Kmer(original); rc.rcomp(); warm(rc);
				final byte[] reverseDeleted=new byte[k];
				for(int i=0; i<k-1; i++){reverseDeleted[i]=reverse[i<p ? i : i+1];}
				reverseDeleted[k-1]=(byte)x;
				check(rc.deleteBase(p, x)==reverse[p], "Deletion after rcomp", requested, p, x);
				equal(rc, build(requested, reverseDeleted), k, requested, p, x);
				rc.insertBase(p, reverse[p]); equal(rc, build(requested, reverse), k, requested, p, x);
				final byte[] reverseInserted=new byte[k];
				for(int i=0; i<k; i++){reverseInserted[i]=i==p ? (byte)x : reverse[i<p ? i : i-1];}
				check(rc.insertBase(p, x)==reverse[k-1], "Insertion after rcomp", requested, p, x);
				equal(rc, build(requested, reverseInserted), k, requested, p, x);
				rc.deleteBase(p, reverse[k-1]); equal(rc, build(requested, reverse), k, requested, p, x);
			}
		}
	}

	/** Rebuilds one complete numeric-base window and populates its lazy caches.
	 * @param requested Requested K under the active geometry settings
	 * @param seq Numeric bases whose length equals the effective K
	 * @return New independently populated key */
	private static Kmer build(final int requested, final byte[] seq){
		final Kmer key=new Kmer(requested);
		assert(seq.length==key.kbig) : "Independent rebuild must contain exactly K bases.";
		for(byte b : seq){key.addRightNumeric(b);}
		warm(key); return key;
	}

	/** Populates canonical-key, hash and palindrome caches before editing/comparison.
	 * @param key Complete key, possibly with rolling length greater than K */
	private static void warm(final Kmer key){
		assert(key.len()>=key.kbig) : "Warm all lazily cached values before an in-place edit.";
		key.key(); key.xor(); key.xor2(); key.corePalindrome();
	}

	/** Compares words/caches with an independent key and checks the expected logical length.
	 * @param actual Edited or restored key
	 * @param expected Independently built comparison key
	 * @param len Expected logical length, which can exceed the fixed window length
	 * @param k Requested K for diagnostics
	 * @param p Edited position for diagnostics
	 * @param b Incoming numeric base for diagnostics */
	private static void equal(final Kmer actual, final Kmer expected, final int len, final int k, final int p, final int b){
		check(Arrays.equals(actual.array1(), expected.array1()) && Arrays.equals(actual.array2(), expected.array2()), "Forward/RC words", k, p, b);
		check(Arrays.equals(actual.key(), expected.key()) && actual.xor()==expected.xor() && actual.xor2()==expected.xor2() && actual.corePalindrome()==expected.corePalindrome(), "Canonical key/hash cache", k, p, b);
		check(actual.len()==len && actual.verify(true), "Length and RC invariant", k, p, b);
	}

	/** Checks invalid-edit rejection and unchanged complete-key state.
	 * Incomplete-window cases check rejection and preserved logical length. */
	private static void invalid(){
		final Kmer key=build(31, new byte[31]);
		for(int op=0; op<2; op++){
			for(int bad=0; bad<4; bad++){
				final int p=bad==0 ? -1 : bad==1 ? 31 : 5;
				final long b=bad==2 ? -1 : bad==3 ? 4 : 0;
				boolean rejected=false;
				try{if(op==0){key.deleteBase(p, b);}else{key.insertBase(p, b);}}catch(AssertionError expected){rejected=true;}
				check(rejected, "Invalid edit rejected before mutation", 31, p, (int)b);
				equal(key, build(31, new byte[31]), 31, 31, p, (int)b);
			}
		}
		final Kmer shortKey=new Kmer(31); shortKey.addRightNumeric(0);
		boolean rejected=false; try{shortKey.deleteBase(0, 0);}catch(AssertionError expected){rejected=true;}
		check(rejected && shortKey.len()==1, "Incomplete deletion rejected", 31, 0, 0);
		rejected=false; try{shortKey.insertBase(0, 0);}catch(AssertionError expected){rejected=true;}
		check(rejected && shortKey.len()==1, "Incomplete insertion rejected", 31, 0, 0);
	}

	/** Counts one explicit comparison and throws immediately on failure, even under -da.
	 * @param ok Comparison result
	 * @param what Contract being checked
	 * @param k Requested K for diagnostics
	 * @param p Edited position for diagnostics
	 * @param b Incoming numeric base for diagnostics */
	private static void check(final boolean ok, final String what, final int k, final int p, final int b){
		checks++;
		if(!ok){throw new AssertionError(what+"; K="+k+", position="+p+", base="+b+", packed="+Kmer.PACKED+", core="+Kmer.MASK_CORE);}
	}

	/** Cumulative explicit comparisons; assertion-only preconditions are not included. */
	private static long checks;
}
