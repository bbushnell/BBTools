package test.assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import assemble.HomopolymerIndelProposal;
import assemble.LocalEditKmerProbe;
import ukmer.Kmer;

/**
 * Checks shifted probe groups against materialized edited reads and an exhaustive
 * legal-start search. Covers internal window boundaries, read ends with fewer
 * distinct windows, undefined-base vetoes and direct versus lazy indel probing.
 * Edit construction and group selection are independent of the production clamp
 * arithmetic; Kmer encoding remains shared. Counts are synthetic key hashes, not
 * evidence from a real read table, and this is not a correction-accuracy test.
 * PACKED is varied/restored, while MASK_CORE is inherited. Use an isolated process
 * because layout flags are global. Explicit checks remain active without -ea.
 * @author Fischl */
public final class LocalEditShiftedWindowProbeTest{

	/** Runs both layouts, three effective K values and requested window maxima of three/five.
	 * Each read is tested with no undefined base or one at either chosen flank position;
	 * every called locus uses three supplied window starts. Arguments are ignored.
	 * PACKED is restored on exit; counters accumulate across invocations in one JVM.
	 */
	public static void main(final String[] args){
		final boolean packed=Kmer.PACKED;
		try{
			for(boolean layout:new boolean[]{false, true}){
				Kmer.PACKED=layout;
				for(int k:new int[]{31, 62, 63}){
					for(int n:new int[]{3, 5}){
						final StringBuilder text=new StringBuilder();
						final Random random=new Random(2026091103L+k);
						for(int i=0; i<k+8; i++){text.append("ACGT".charAt(random.nextInt(4)));}
						for(int undefined:new int[]{-1, 3, k+6}){
							final StringBuilder source=new StringBuilder(text);
							if(undefined>=0){source.setCharAt(undefined, 'N');}
							for(int p=0; p<source.length(); p++){
								if(p!=undefined){
									for(int w:new int[]{Math.min(p, source.length()-k), Math.max(0, p-k+1), Math.max(0, Math.min(p-k/2, source.length()-k))}){
										exact(source.toString(), k, n, w, p);
									}
								}
							}
						}
					}
				}
			}
		}finally{Kmer.PACKED=packed;}
		System.out.println("LOCAL_EDIT_SHIFTED_WINDOW_TEST_OK cases="+cases+" checks="+checks+"; complete edited-window oracle, shifted groups, fewer distinct windows at ends, N veto, lazy budgets.");
	}

	/** Materializes each family's chosen windows, then checks the probe's ordered evidence.
	 * Unavailable windows contribute -1 to the minimum without creating lookup entries.
	 * Also checks direct/lazy indel equality, deletion/insertion availability and unchanged
	 * input. The unchanged substitution allele remains unprobed at -1.
	 * @param source Original read containing A/C/G/T and at most the fixture's one N
	 * @param k Effective K for the current layout
	 * @param n Requested odd maximum adjacent windows per edit family
	 * @param w Supplied in-bounds original window start, used as the preferred group center
	 * @param p Called original-read position inside that window
	 */
	private static void exact(final String source, final int k, final int n, final int w, final int p){
		final Audit lookup=new Audit(k);
		final int original="ACGT".indexOf(source.charAt(p));
		final int[] expectedS={-1, -1, -1, -1}, expectedI={-1, -1, -1, -1}; int expectedD=-1;
		final int[] startsS=starts(source.length(), k, n, w, p, 0), startsD=starts(source.length()-1, k, n, w, p, 1), startsI=starts(source.length()+1, k, n, w, p, 2);
		for(int j=0; j<startsS.length; j++){
			for(int b=0; b<4; b++){
				if(b!=original){
					final String edited=source.substring(0, p)+"ACGT".charAt(b)+source.substring(p+1);
					final int depth=lookup.add(edited.substring(startsS[j], startsS[j]+k));
					expectedS[b]=j==0 ? depth : Math.min(expectedS[b], depth);
				}
			}
		}
		final int substitutionQueries=lookup.keys.size();
		for(int j=0; j<startsD.length; j++){
			final String edited=source.substring(0, p)+source.substring(p+1);
			final int depth=lookup.add(edited.substring(startsD[j], startsD[j]+k));
			expectedD=j==0 ? depth : Math.min(expectedD, depth);
		}
		for(int j=0; j<startsI.length; j++){
			for(int b=0; b<4; b++){
				final String edited=source.substring(0, p)+"ACGT".charAt(b)+source.substring(p);
				final int depth=lookup.add(edited.substring(startsI[j], startsI[j]+k));
				expectedI[b]=j==0 ? depth : Math.min(expectedI[b], depth);
			}
		}
		final byte[] bases=source.getBytes(StandardCharsets.US_ASCII), before=bases.clone();
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(k, lookup, n);
		probe.substitutions(bases, w, p);
		check(probe.queries==substitutionQueries, "S stage must query only defined substitution windows.");
		check(Arrays.equals(expectedS, probe.substitutionDepth), "Shifted S minima differ from materialized oracle.");
		probe.indels();
		check(probe.deletionDepth==expectedD && Arrays.equals(expectedI, probe.insertionDepth), "D/I windows or minimum attribution differ from the materialized oracle.");
		check(probe.queries==lookup.keys.size() && lookup.cursor==lookup.keys.size() && probe.queries<=8*n, "Missing/repeated/wasted window lookups.");
		check(probe.deletionAvailable==(expectedD>=0) && probe.insertionsAvailable==(expectedI[0]>=0), "Availability must reflect the complete chosen group, not just one member.");
		probe.indels(); check(probe.queries==lookup.keys.size(), "Lazy indels must remain cached.");
		check(Arrays.equals(before, bases), "Window shifting must not edit the input.");
		final Audit direct=new Audit(k);
		direct.keys.addAll(lookup.keys.subList(substitutionQueries, lookup.keys.size()));
		final LocalEditKmerProbe only=new LocalEditKmerProbe(k, direct, n);
		only.indelsAt(bases, w, p);
		check(only.queries==direct.keys.size() && direct.cursor==direct.keys.size(), "Direct indel stage repeats S probes or misses D/I windows.");
		check(only.deletionDepth==expectedD && Arrays.equals(only.insertionDepth, expectedI), "Direct and lazy indel paths must agree.");
		cases++;
	}

	/** Exhaustively finds legal spanning starts and chooses the nearest complete group.
	 * Uses fewer starts only when fewer exist. Deletion must straddle the junction;
	 * substitution/insertion must contain the changed base. Ties retain the earlier group.
	 * Distance compares the group's first start with w-requested/2, even for a smaller group.
	 * This enumeration intentionally avoids the production clamp expression.
	 * @param editedLength Length after the virtual edit
	 * @param k Complete window length
	 * @param requested Requested odd maximum window count
	 * @param w Original preferred center start
	 * @param p Changed-base index or deletion junction in the edited read
	 * @param family Zero for substitution, one for deletion, two for insertion
	 * @return Ordered starts for the nearest group, or an empty array when none is legal
	 */
	private static int[] starts(final int editedLength, final int k, final int requested, final int w, final int p, final int family){
		final ArrayList<Integer> valid=new ArrayList<Integer>();
		for(int s=0; s+k<=editedLength; s++){
			if(s<=p && p<s+k && (family!=1 || s<p)){valid.add(s);}
		}
		final int size=Math.min(requested, valid.size());
		if(size==0){return new int[0];}
		int best=-1, bestDistance=Integer.MAX_VALUE;
		for(int j=0; j+size<=valid.size(); j++){
			final int distance=Math.abs(valid.get(j)-(w-requested/2));
			if(distance<bestDistance){best=j; bestDistance=distance;}
		}
		check(best>=0, "A nonempty legal-start set must contain a nearest group.");
		final int[] result=new int[size];
		for(int j=0; j<size; j++){result[j]=valid.get(best+j);}
		return result;
	}

	/** Ordered defined-key oracle with deterministic synthetic counts and no real table. */
	private static final class Audit implements HomopolymerIndelProposal.CountLookup{
		/** Creates an empty query sequence for one effective K. */
		Audit(final int k_){k=k_;}
		/** Returns -1 for an N-containing window without adding a query; otherwise snapshots its key. */
		int add(final String s){
			if(s.indexOf('N')>=0){return -1;}
			final Kmer key=new Kmer(k);
			for(int j=0; j<s.length(); j++){key.addRight(s.charAt(j));}
			final long[] value=key.key().clone(); keys.add(value); return depth(value);
		}
		/** Checks exact key/order and returns its synthetic depth without retaining the borrowed key. */
		@Override public int count(final Kmer key){
			check(cursor<keys.size() && Arrays.equals(keys.get(cursor), key.key()), "Queried key differs from the independent edited-read window.");
			return depth(keys.get(cursor++));
		}
		/** Maps an encoded key to a repeatable count in 0..30, including zero. */
		private static int depth(final long[] key){return (Arrays.hashCode(key)&Integer.MAX_VALUE)%31;}
		/** Effective encoded word length. */
		final int k;
		/** Next expected query index. */
		int cursor;
		/** Key snapshots; the direct-indel oracle shares these arrays without changing them. */
		final ArrayList<long[]> keys=new ArrayList<long[]>();
	}

	/** Counts each check attempt and throws independently of JVM assertion settings. */
	private static void check(final boolean ok, final String why){checks++; if(!ok){throw new AssertionError(why);}}
	/** Completed locus cases and total check attempts, cumulative across main invocations. */
	private static long cases, checks;
}
