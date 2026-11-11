package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import assemble.HomopolymerIndelProposal;
import assemble.LocalEditKmerProbe;
import ukmer.Kmer;

/**
 * Checks single-window candidate keys, depth attribution and lazy lookup budgets.
 * Candidate sequences are materialized independently as strings, then encoded with
 * Kmer for comparison with the probe's reused key. This tests edit construction,
 * not the correctness of Kmer's shared encoding or correction accuracy.
 * Uses synthetic counts and no real count table. The main method temporarily changes
 * global Kmer layout flags, so run in an isolated process, not beside other Kmer users.
 * Explicit checks throw AssertionError even when JVM assertions are disabled.
 * @author Fischl */
public final class LocalEditKmerProbeTest{

	/** Runs all four layout/core-mask combinations and restores their original settings.
	 * Arguments are ignored. Success prints only after all checks complete; counters are
	 * cumulative across repeated invocations in the same JVM. Failure checks run once,
	 * after the layout matrix, using its final settings.
	 */
	public static void main(final String[] args){
		final boolean packed=Kmer.PACKED, core=Kmer.MASK_CORE;
		try{
			for(boolean p:new boolean[]{false, true}){
				for(boolean c:new boolean[]{false, true}){
					Kmer.PACKED=p; Kmer.MASK_CORE=c;
					for(int k:new int[]{31, 62, 63}){
						for(int pos:new int[]{0, k/2, k-1}){
							for(int seed=0; seed<8; seed++){exact(k, pos, seed);}
						}
						invalid(k);
					}
				}
			}
			failures();
		}finally{Kmer.PACKED=packed; Kmer.MASK_CORE=core;}
		System.out.println("LOCAL_EDIT_PROBE_TEST_OK checks="+checks+" cases="+cases+"; exact 3+1+4 candidate keys/depths, lazy lookup budget, edges/undefined, K31/62/63 packed/legacy/core. Probe only, no correction accuracy claim.");
	}

	/** Checks one fully defined locus against ordered, independently edited strings.
	 * Includes first/middle/last positions within an interior window, all eight candidate
	 * depths, repeated lazy calls, direct indel-only query order and unchanged input bytes.
	 * @param k Effective word length
	 * @param relative Zero-based edit position within the original window
	 * @param seed Offset from the fixed random seed used to create the read
	 */
	private static void exact(final int k, final int relative, final int seed){
		final Random random=new Random(2026091005L+seed);
		final StringBuilder text=new StringBuilder(k+7);
		for(int i=0; i<k+7; i++){text.append("ACGT".charAt(random.nextInt(4)));}
		final String source=text.toString();
		final byte[] bases=source.getBytes(java.nio.charset.StandardCharsets.US_ASCII), before=bases.clone();
		final int w=3, p=w+relative, original="ACGT".indexOf(source.charAt(p));
		final AuditLookup lookup=new AuditLookup(k);
		for(int b=0; b<4; b++){
			if(b!=original){lookup.add(source.substring(w, p)+"ACGT".charAt(b)+source.substring(p+1, w+k));}
		}
		lookup.add(source.substring(w, p)+source.substring(p+1, w+k+1));
		for(int b=0; b<4; b++){lookup.add(source.substring(w, p)+"ACGT".charAt(b)+source.substring(p, w+k-1));}
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(k, lookup);
		probe.substitutions(bases, w, p);
		check(probe.queries==3 && lookup.cursor==3 && probe.substitutionsAvailable, "Substitution stage must not perform eager indel lookups.");
		check(probe.substitutionDepth[original]==-1, "The original allele is not one of the three substitution probes.");
		int ordinal=0;
		for(int b=0; b<4; b++){
			if(b!=original){check(probe.substitutionDepth[b]==depth(ordinal++), "Substitution depth attributed to wrong base.");}
		}
		probe.indels();
		check(probe.queries==8 && lookup.cursor==8 && probe.deletionAvailable && probe.insertionsAvailable, "Fully defined interior context must make exactly eight probes.");
		check(probe.deletionDepth==depth(3), "Deletion depth must correspond to the deleted-position candidate.");
		for(int b=0; b<4; b++){check(probe.insertionDepth[b]==depth(4+b), "Insertion depth attributed to wrong base.");}
		probe.indels(); check(probe.queries==8, "Repeated indel query must use the cached candidates.");
		check(Arrays.equals(bases, before), "Probes must never mutate input sequence."); cases++;
		final AuditLookup direct=new AuditLookup(k);
		for(int i=3; i<8; i++){direct.expected.add(lookup.expected.get(i));}
		final LocalEditKmerProbe indelOnly=new LocalEditKmerProbe(k, direct);
		indelOnly.indelsAt(bases, w, p);
		check(indelOnly.queries==5 && direct.cursor==5 && !indelOnly.substitutionsAvailable, "Direct indel stage must not repeat substitutions.");
		indelOnly.indels(); check(indelOnly.queries==5, "Direct indel results are cached too.");
	}

	/** Checks read-end and undefined-base availability with every queried key absent.
	 * Missing or retained undefined bases suppress queries; an insertion can discard
	 * the last undefined original-window base while deletion needs one extra called base.
	 */
	private static void invalid(final int k){
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){@Override public int count(Kmer key){return -1;}};
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(k, lookup);
		final byte[] bases=new byte[k+1]; Arrays.fill(bases, (byte)'A');
		probe.substitutions(Arrays.copyOf(bases, k), 0, k/2); probe.indels();
		check(probe.queries==7 && !probe.deletionAvailable && probe.deletionDepth==-1 && probe.insertionsAvailable, "Right read edge permits insertions, not a deletion needing an extra base.");
		bases[k-1]='N'; probe.substitutions(bases, 0, k/2); probe.indels();
		check(!probe.substitutionsAvailable && !probe.deletionAvailable && probe.insertionsAvailable && probe.queries==4, "Only insertion drops this terminal undefined base from its candidate window.");
		bases[k-1]='A'; bases[k]='N'; probe.substitutions(bases, 0, k/2); probe.indels();
		check(probe.queries==7 && !probe.deletionAvailable, "Undefined extra right base must invalidate deletion only.");
		bases[k]='A'; bases[1]='N'; probe.substitutions(bases, 0, k/2); probe.indels();
		check(probe.queries==0 && !probe.substitutionsAvailable && !probe.deletionAvailable && !probe.insertionsAvailable, "Undefined retained interior base must prevent all table queries.");
	}

	/** Expected normalized depth: ordinal zero represents an absent key, then distinct counts. */
	private static int depth(final int ordinal){return ordinal==0 ? 0 : ordinal+10;}

	/** Checks rejected count values, failed-locus invalidation and recovery on fresh input.
	 * The synthetic lookup returns -2 on its fourth call, then valid positive counts.
	 * An undefined requested allele must also invalidate previously cached indel state.
	 */
	private static void failures(){
		final byte[] bases=new byte[40]; Arrays.fill(bases, (byte)'A');
		final int[] calls={0};
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(31, new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(Kmer key){return ++calls[0]==4 ? -2 : 5;}
		});
		probe.substitutions(bases, 0, 15); boolean failed=false;
		try{probe.indels();}catch(IllegalStateException e){failed=true;}
		check(failed, "Invalid lookup depth must fail loudly.");
		failed=false;
		try{probe.indels();}catch(IllegalStateException e){failed=true;}
		check(failed && calls[0]==4, "A failed indel probe must not masquerade as cached success or retry partially.");
		probe.substitutions(bases, 0, 15); probe.indels();
		check(probe.queries==8, "Fresh input call must reset failure state.");
		bases[15]='N'; failed=false;
		try{probe.substitutions(bases, 0, 15);}catch(IllegalArgumentException e){failed=true;}
		check(failed, "N at requested position requires four alternatives and must not enter the three-substitution API.");
		failed=false;
		try{probe.indels();}catch(IllegalStateException e){failed=true;}
		check(failed, "An invalid new input must invalidate the previous cached locus.");
	}

	/** Ordered synthetic lookup; stores cloned encodings, never the borrowed probe key. */
	private static final class AuditLookup implements HomopolymerIndelProposal.CountLookup{
		/** Creates an empty expected-key sequence for one effective K. */
		AuditLookup(final int k_){k=k_;}
		/** Appends the canonical encoding of a materialized, exactly k-base candidate. */
		void add(final String s){
			check(s.length()==k, "Materialized candidate oracle must contain exactly k bases.");
			final Kmer key=new Kmer(k);
			for(int i=0; i<s.length(); i++){key.addRight(s.charAt(i));}
			expected.add(key.key().clone());
		}
		/** Checks the next complete encoded key; returns absent -1 first, then distinct depths. */
		@Override public int count(final Kmer key){
			check(cursor<expected.size() && Arrays.equals(key.key(), expected.get(cursor)), "Packed candidate differs from independently materialized edited sequence.");
			final int ordinal=cursor++; return ordinal==0 ? -1 : depth(ordinal);
		}
		/** Required materialized candidate length. */
		final int k;
		/** Number of expected keys already consumed. */
		int cursor;
		/** Ordered key snapshots; the direct-indel fixture shares their arrays without mutating them. */
		final ArrayList<long[]> expected=new ArrayList<long[]>();
	}

	/** Counts every check attempt and fails independently of JVM assertion settings. */
	private static void check(final boolean ok, final String why){checks++; if(!ok){throw new AssertionError(why);}}
	/** Lifetime check attempts and completed exact-locus primary checks, respectively. */
	private static long checks, cases;
}
