package test.assemble;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import assemble.HomopolymerIndelProposal;
import assemble.LocalEditEngine;
import assemble.LocalEditPatchCorrector;
import assemble.LocalSingleBaseEdit;
import dna.AminoAcid;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Integration controls for explicitly supplied original-coordinate patch regions.
 * Checks four k-mer lengths, both strands and optional qualities against known
 * sequence truth and the whole-read sequential correction path. Synthetic counts
 * use the production Kmer encoding; quality expectations share the sequential engine.
 * Covers read-wide budgets, retained input arrays and selected failure/reuse paths.
 * Does not test automatic grouping, population accuracy or performance.
 * @author Fischl */
public final class LocalEditPatchCorrectorTest{

	/** Runs the fixtures without arguments and restores the two modified Kmer flags.
	 * Requires isolated use: those flags are process-wide. Checks do not rely on -ea. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Deterministic patch integration test takes no inputs.");}
		final boolean packed=Kmer.PACKED, mask=Kmer.MASK_CORE;
		try{
			Kmer.PACKED=true; Kmer.MASK_CORE=false;
			for(int k:new int[]{31, 62, 95, 127}){for(boolean reverse:new boolean[]{false, true}){for(boolean quality:new boolean[]{false, true}){pairedRegions(k, reverse, quality);}}}
			failures();
			System.out.println("LOCAL_EDIT_PATCH_CORRECTOR_TEST_OK checks="+checks+"; supplied-region integration, no automatic grouping or benchmark claim.");
		}finally{Kmer.PACKED=packed; Kmer.MASK_CORE=mask;}
	}
	/** Tests two supplied regions containing four known edits, not paired-end reads.
	 * Exact-budget success, later budget rollback and initial-burden rejection use
	 * independent Read copies while reusing one patch-correction worker. */
	private static void pairedRegions(final int k, final boolean reverse, final boolean quality){
		final byte[] truth=sequence(20*k, 1603+k); final Oracle table=new Oracle(k, truth);
		final Read input=noisy(truth, k, quality); if(reverse){reverse(input);}
		final IntList regions=regions(input.length(), k, reverse);
		final Read expected=copy(input), actual=copy(input);
		final byte[] before=actual.bases, qbefore=actual.quality, saved=before.clone(), qsaved=qbefore==null ? null : qbefore.clone();
		final LocalEditPatchCorrector patch=new LocalEditPatchCorrector(k, table, 1, false, 8*k);
		check(new LocalEditEngine(k, table::count, 1).correct(expected, 8)==4, "The whole-read sequential fixture must recover all four known S/I/D errors.");
		check(patch.correct(actual, regions, 4)==4 && patch.stagedPatches==2, "Two independent local regions must stage all four repairs at the exact global limit.");
		check(Arrays.equals(actual.bases, reverse ? AminoAcid.reverseComplementBases(truth) : truth), "One final patch transaction must recover the full truth sequence.");
		check(Arrays.equals(actual.bases, expected.bases) && Arrays.equals(actual.quality, expected.quality), "Local-region integration must match whole-read sequential bases and qualities, including reverse/null-quality inputs.");
		check(Arrays.equals(before, saved) && Arrays.equals(qbefore, qsaved), "Local integration must not mutate original array aliases.");
		check(patch.substitutions==2 && patch.insertions==1 && patch.deletions==1 && patch.materializedBases==actual.length(), "Only committed local operations and the one final materialization enter output counters.");
		for(int limit:new int[]{2, 3}){
			final Read rejected=copy(input); final byte[] b=rejected.bases, q=rejected.quality;
			check(patch.correct(rejected, regions, limit)==0 && patch.budgetRejected && !patch.initialRejected && patch.stagedPatches==1, "Later local edits must trigger whole-read rollback, including zero remaining budget after the first patch.");
			check(rejected.bases==b && rejected.quality==q && patch.materializedBases==0 && patch.substitutions+patch.insertions+patch.deletions==0, "A read rejected after staging its first region must publish no partial sequence or edit counters.");
		}
		final Read initial=copy(input); final byte[] ib=initial.bases;
		check(patch.correct(initial, regions, 1)==0 && patch.initialRejected && initial.bases==ib, "The original whole-read burden preflight must run before local patches receive budgets.");
		check(patch.correct(copy(input), regions, 8)==4, "The same worker must be reusable after both preflight and later-budget rejection.");
	}
	/** Exercises bounded scratch, malformed regions, edits outside the mutable core
	 * and one tagged lookup failure after the first region has been staged. */
	private static void failures(){
		final int k=31, p=4*k, q=p+k+1; final byte[] truth=sequence(20*k, 1603+k); final Oracle table=new Oracle(k, truth);
		final LocalEditPatchCorrector patch=new LocalEditPatchCorrector(k, table, 1, false, 8*k);
		final Read input=noisy(truth, k, true); final byte[] b=input.bases, quals=input.quality;
		final IntList restricted=new IntList(); add(restricted, p-2*k, q+2*k, p-2, p+2);
		check(patch.correct(input, restricted, 8)==0 && patch.boundaryRejected==1 && input.bases==b && input.quality==quals, "Discovering an edit in padding must discard that local candidate, not change a neighboring region.");
		final IntList huge=new IntList(); add(huge, 0, input.length(), k, input.length()-k);
		check(patch.correct(input, huge, 8)==0 && patch.sizeRejected==1 && input.bases==b, "Scratch work is bounded by the explicit patch-size ceiling.");
		final IntList invalid=new IntList(); add(invalid, 10, 20, Integer.MIN_VALUE, Integer.MIN_VALUE+100);
		expect(()->patch.correct(input, invalid, 8), "Malformed core coordinates must fail before arithmetic can wrap into apparently valid flanks.");
		final IntList overlapping=regions(input.length(), k, false); add(overlapping, 0, 3*k, k, 2*k);
		expect(()->patch.correct(input, overlapping, 8), "Supplied copy regions must be in nonoverlapping original order.");
		table.fail=key(truth, 12*k-k/2, k);
		boolean failed=false; try{patch.correct(input, regions(input.length(), k, false), 8);}catch(IllegalStateException expected){failed=true;}
		check(failed && patch.stagedPatches==1, "Injected lookup failure must occur after the first patch was staged, not only in preflight.");
		check(input.bases==b && input.quality==quals && patch.materializedBases==0, "A later-region exception leaves the complete original read installed.");
		table.fail=null; check(patch.correct(input, regions(input.length(), k, false), 8)==4, "A worker recovers after abandoning a failed multi-region transaction.");
	}
	/** Returns two copy/core quadruples in increasing coordinates on the requested strand. */
	private static IntList regions(final int length, final int k, final boolean reverse){
		final IntList out=new IntList();
		if(reverse){region(out, length, k, 12*k, true); region(out, length, k, 4*k, true);
		}else{region(out, length, k, 4*k, false); region(out, length, k, 12*k, false);}
		return out;
	}
	/** Appends a prescribed half-open core and K-base flanks, reflecting all endpoints if needed. */
	private static void region(final IntList out, final int length, final int k, final int p, final boolean reverse){
		final int a=p-4, b=p+k+6, from=a-k, to=b+k;
		if(reverse){add(out, length-to, length-from, length-b, length-a);}else{add(out, from, to, a, b);}
	}
	/** Appends copy-start/end and core-start/end in the executor's four-field order. */
	private static void add(final IntList out, final int from, final int to, final int coreFrom, final int coreTo){out.add(from); out.add(to); out.add(coreFrom); out.add(coreTo);}
	/** Clones truth and introduces two substitutions, one deletion and one insertion.
	 * Optional initial qualities cycle from fifteen to thirty-nine; the editor assigns
	 * qualities to changed bases. Edits are applied from right to left. */
	private static Read noisy(final byte[] truth, final int k, final boolean quality){
		final byte[] qs=quality ? new byte[truth.length] : null; if(qs!=null){for(int i=0; i<qs.length; i++){qs[i]=(byte)(15+i%25);}}
		final Read r=new Read(truth.clone(), qs, "two-local-regions", 0, false); final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();
		editor.apply(r, LocalSingleBaseEdit.Operation.SUBSTITUTION, 13*k+1, other(r.bases[13*k+1]));
		editor.apply(r, LocalSingleBaseEdit.Operation.DELETION, 12*k);
		editor.apply(r, LocalSingleBaseEdit.Operation.INSERTION, 5*k+1, r.bases[5*k+1]);
		editor.apply(r, LocalSingleBaseEdit.Operation.SUBSTITUTION, 4*k, other(r.bases[4*k])); return r;
	}
	/** Chooses a different defined base for the uppercase fixture alphabet. */
	private static byte other(final byte b){return b=='A' ? (byte)'C' : (byte)'A';}
	/** Copies sequence/quality arrays and identifiers into a fresh unpaired Read. */
	private static Read copy(final Read r){return new Read(r.bases.clone(), r.quality==null ? null : r.quality.clone(), r.id, r.numericID, false);}
	/** Reverse-complements owned bases and reverses optional qualities in place. */
	private static void reverse(final Read r){AminoAcid.reverseComplementBasesInPlace(r.bases); if(r.quality!=null){for(int a=0, b=r.quality.length-1; a<b; a++, b--){final byte v=r.quality[a]; r.quality[a]=r.quality[b]; r.quality[b]=v;}}}
	/** Generates positive-length seeded A/C/G/T truth with terminal A bases. */
	private static byte[] sequence(final int length, final long seed){final Random random=new Random(seed); final byte[] b=new byte[length]; for(int i=0; i<length; i++){b[i]=ALPHABET[random.nextInt(4)];}b[0]=b[length-1]='A'; return b;}
	/** Encodes one truth window with the same Kmer representation used by lookup queries. */
	private static String key(final byte[] b, final int start, final int k){final Kmer word=new Kmer(k); for(int i=start; i<start+k; i++){word.addRight(b[i]);}return Arrays.toString(word.key());}
	/** Synthetic membership table: truth words have count twelve, absent words zero.
	 * Repeated truth words retain twelve rather than accumulating multiplicity. */
	private static final class Oracle implements HomopolymerIndelProposal.CountLookup{
		/** Records every complete truth window at a fixed supporting count. */
		Oracle(final int k, final byte[] truth){for(int i=0; i<=truth.length-k; i++){counts.put(key(truth, i, k), 12);}}
		/** Throws only for the selected failing key; otherwise returns the synthetic count. */
		@Override public int count(final Kmer word){final String s=Arrays.toString(word.key()); if(s.equals(fail)){throw new IllegalStateException("Injected later-region lookup failure.");}final Integer n=counts.get(s); return n==null ? 0 : n;}
		/** Fixture-owned encoded membership table. */
		final HashMap<String, Integer> counts=new HashMap<String, Integer>();
		/** Encoded failure-injection key, or null for normal lookup. */
		String fail;
	}
	/** Requires argument rejection; other exception types propagate. */
	private static void expect(final Runnable action, final String why){boolean failed=false; try{action.run();}catch(IllegalArgumentException expected){failed=true;}check(failed, why);}
	/** Enforces and counts a fixture condition regardless of JVM assertion settings. */
	private static void check(final boolean pass, final String why){checks++; if(!pass){throw new AssertionError(why);}}
	/** Cumulative checks in this JVM; main does not reset this field. */
	private static int checks;
	/** Defined alphabet used to generate sequence truth. */
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
