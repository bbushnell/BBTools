package test.assemble;

import java.util.Arrays;
import java.util.Random;

import assemble.HollowExtensionGuard;
import assemble.KmerWindowDepth;
import dna.AminoAcid;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayH1D;
import ukmer.HollowKmerLookup;
import ukmer.Kmer;

/** Checks rolling live-guard decisions against full diagnostic reconstruction.
 * @author Brian Bushnell, Fischl */
public final class HollowExtensionGuardTest {

	/** Covers sequential growth, resets, both orientations, absent context and buffer replacement. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(42122);
		for(int k : new int[]{67, 72, 96, 124, 128, 160, 192, 224, 310, 342}){
			for(int trial=0; trial<20; trial++){
				for(int mode=0; mode<3; mode++){
					final int seed=k+(mode==0 ? 0 : 5);
					final byte[] bases=new byte[seed+40];
					for(int i=0; i<bases.length; i++){bases[i]=(byte)AminoAcid.numberToBase[random.nextInt(4)];}
					if(mode==2){bases[seed-k-1]='N';}
					exercise(bases, seed, k, random);
					AminoAcid.reverseComplementBasesInPlace(bases);
					// Keep the N at the tested extra-context boundary in each orientation.
					if(mode==2){bases[bases.length-1-(seed-k-1)]='A'; bases[seed-k-1]='N';}
					exercise(bases, seed, k, random);
				}
			}
		}
		check(rejected>0 && retained>0, "Parity fixtures must exercise both decisions, not only empty-table retention.");
		invalidCalls();
		System.out.println("HOLLOW_EXTENSION_GUARD_TEST_PASS checks="+checks+" rejected="+rejected+" retained="+retained);
		final int before=checks;
		solidBranchRules();
		System.out.println("SOLID_BRANCH_PREDICATE_TEST_PASS checks="+(checks-before));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Exhausts support, symmetry, strict-ratio boundaries and overflow-safe old-mode parity. */
	private static void solidBranchRules(){
		for(int context : new int[]{0, 1, 2, 3, 22, 25, Integer.MAX_VALUE}){
			for(int chosen : new int[]{0, 1, 2, 3, 9, 12, 13, 23, Integer.MAX_VALUE}){
				for(int other : new int[]{0, 1, 2, 3, 9, 12, 13, 23, Integer.MAX_VALUE}){
					for(int support : new int[]{2, 3}){
						for(int ratio : new int[]{2, 4, Integer.MAX_VALUE}){
							final boolean old=context>=support && chosen==0 && other>=support;
							final boolean ambiguous=context>=support && chosen>=support && other>=support &&
									(long)chosen<(long)ratio*other && (long)other<(long)ratio*chosen;
							check(KmerWindowDepth.rejectSolidStep(context, chosen, other, support, ratio, false)==old,
									"Disabled Solid ambiguity must preserve the zero-choice contradiction predicate.");
							check(KmerWindowDepth.rejectSolidStep(context, chosen, other, support, ratio, true)==(old || ambiguous),
									"Both branches need support and neither may dominate; ratio multiplication must not overflow.");
						}
					}
				}
			}
		}
		check(!KmerWindowDepth.rejectSolidStep(25, 23, 1, 3, 4, true), "The supported correct checkpoint must remain traversable.");
		check(KmerWindowDepth.rejectSolidStep(22, 13, 9, 3, 4, true), "The supported ambiguous checkpoint must stop without choosing a replacement.");
		check(!KmerWindowDepth.rejectSolidStep(2, 13, 9, 3, 4, true), "Weak context cannot establish a Solid ambiguity.");
		check(!KmerWindowDepth.rejectSolidStep(22, 12, 3, 3, 4, true), "Exact fourfold dominance is not ambiguous.");
		check(KmerWindowDepth.rejectSolidStep(3, 0, 3, 3, 1, false), "Legacy off-mode must still accept its unused ratio of one.");
	}

	/** Loads every diagnostic context/neighbor first, then compares against an immutable table. */
	private static void exercise(final byte[] bases, final int seed, final int k, final Random random){
		assert(seed>=k && bases.length>seed) : "Every fixture has a valid nominal seed and at least one candidate.";
		final byte[] untouched=bases.clone();
		final AbstractKmerTableU[] tables=new AbstractKmerTableU[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArrayH1D(new int[]{1009});}
		final Kmer word=new Kmer(k, true), neighbor=new Kmer(k, true);
		for(int end=seed; end<bases.length; end++){
			for(int span : new int[]{k+1, k-1}){
				if(!defined(bases, end-span, end)){continue;}
				KmerWindowDepth.setHollowQuery(word, bases, end-span, span);
				put(tables, word, random.nextInt(80)+1);
				for(int base=0; base<4; base++){
					neighbor.setFrom(word); neighbor.addRightNumeric(base);
					put(tables, neighbor, random.nextInt(80)+1);
				}
			}
		}
		final HollowKmerLookup lookup=new HollowKmerLookup(tables, k, true);
		final HollowExtensionGuard rolling=new HollowExtensionGuard(lookup, k, 4, 3);
		final HollowExtensionGuard repeatedCalls=new HollowExtensionGuard(lookup, k, 4, 3);
		rolling.reset(bases, seed);
		for(int end=seed; end<bases.length; end++){
			final int candidate=AminoAcid.baseToNumber[bases[end]];
			final boolean expected=oracle(tables, bases, end, candidate, word, neighbor);
			if(expected){rejected++;}else{retained++;}
			check(rolling.reject(bases, end, candidate)==expected, "Rolling decision differs at K="+k+", end="+end);
			check(rolling.reject(bases, end, candidate)==expected, "Repeated candidate queries must not advance the endpoint.");
			// A larger backing buffer must not leave stale array references in the guard.
			check(rolling.reject(Arrays.copyOf(bases, bases.length+97), end, candidate)==expected, "Buffer replacement changed the decision.");
			repeatedCalls.reset(bases, end);
			check(repeatedCalls.reject(bases, end, candidate)==expected, "A new extension call must agree with uninterrupted traversal.");
		}
		check(Arrays.equals(bases, untouched), "The guard must never edit any seed or extension base.");
		// Reusing one worker for a shorter new read must clear its previous endpoint.
		rolling.reset(bases, seed);
		check(rolling.reject(bases, seed, AminoAcid.baseToNumber[bases[seed]])==
				oracle(tables, bases, seed, AminoAcid.baseToNumber[bases[seed]], word, neighbor), "Reset leaked state from the preceding read.");
	}

	/** Rebuilds both query words and applies the already validated diagnostic predicates. */
	private static boolean oracle(final AbstractKmerTableU[] tables, final byte[] bases, final int end,
			final int candidate, final Kmer word, final Kmer neighbor){
		assert(candidate>=0 && candidate<4) : "Oracle candidates are defined bases selected from the fixture path.";
		for(int span : new int[]{word.kbig+1, word.kbig-1}){
			if(!defined(bases, end-span, end)){continue;}
			KmerWindowDepth.setHollowQuery(word, bases, end-span, span);
			final int context=depth(tables, word);
			int chosen=0, other=0;
			for(int base=0; base<4; base++){
				neighbor.setFrom(word); neighbor.addRightNumeric(base);
				final int count=depth(tables, neighbor);
				if(base==candidate){chosen=count;}else{other=Math.max(other, count);}
			}
			if(KmerWindowDepth.rejectPathStep(context, chosen, other, 4, 3) ||
					KmerWindowDepth.ambiguousPathStep(context, chosen, other, 4, 3)){return true;}
		}
		return false;
	}

	/** Unsupported geometry and discontinuous use must fail loudly. */
	private static void invalidCalls(){
		final AbstractKmerTableU[] tables={new HashArrayH1D(new int[]{101})};
		final HollowKmerLookup lookup=new HollowKmerLookup(tables, 128, true);
		boolean failed=false;
		try{new HollowExtensionGuard(lookup, 66, 4, 3);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "K-1 needs two disjoint33bp contexts, so K<67 is unsupported.");
		final HollowExtensionGuard guard=new HollowExtensionGuard(lookup, 128, 4, 3);
		final byte[] bases=new byte[140];
		Arrays.fill(bases, (byte)'A');
		failed=false;
		try{guard.reject(bases, 128, 0);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Queries before reset must not use fabricated endpoint state.");
		guard.reset(bases, 128);
		failed=false;
		try{guard.reject(bases, 130, 0);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Skipping endpoints must require reset.");
		bases[127]='N';
		failed=false;
		try{guard.reset(bases, 128);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "An invalid nominal terminal window must not become count evidence.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Identifies unavailable observed context without encoding missing or N bases. */
	private static boolean defined(final byte[] bases, final int start, final int end){
		assert(end<=bases.length) : "Observed contexts cannot extend past the actual prefix.";
		if(start<0){return false;}
		for(int i=start; i<end; i++){if(AminoAcid.baseToNumber[bases[i]]<0){return false;}}
		return true;
	}

	/** Stores fixture counts with the original diagnostic representation. */
	private static void put(final AbstractKmerTableU[] tables, final Kmer word, final int count){
		assert(count>0) : "HashArrayH1D reserves zero for an absent key.";
		tables[word.mod(tables.length)].set(word, count);
	}

	/** Uses the original packed lookup rather than the endpoint helper being tested. */
	private static int depth(final AbstractKmerTableU[] tables, final Kmer word){
		assert(word.len>=word.kbig) : "Diagnostic counts need a complete table-sized query word.";
		return Math.max(0, tables[word.mod(tables.length)].getValue(word));
	}

	/** Test failures remain fatal if assertions were accidentally disabled. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}

	private static int checks, rejected, retained;
}
