package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import ukmer.Kmer;

/** Independent materialized-string oracle for adjacent-window minima.
 * @author Fischl */
public final class LocalEditWindowProbeTest {
	public static void main(final String[] args){
		final boolean packed=Kmer.PACKED;
		try{
			for(boolean layout:new boolean[]{true,false}){
				Kmer.PACKED=layout;
				for(int k:new int[]{31,62,63}){for(int n:new int[]{3,5}){for(int seed=0;seed<12;seed++){exact(k,n,seed);}}}
			}
			edges();
		}finally{Kmer.PACKED=packed;}
		System.out.println("LOCAL_EDIT_WINDOW_PROBE_TEST_OK checks="+checks+"; materialized 8*N keys/minima, lazy indels, edges, undefined bases, invalid N.");
	}
	private static void exact(final int k,final int n,final int seed){
		final Random random=new Random(2026091101L+seed);final StringBuilder text=new StringBuilder();
		for(int i=0;i<k+20;i++){text.append("ACGT".charAt(random.nextInt(4)));}
		final String raw=text.toString();final byte[] bases=raw.getBytes(StandardCharsets.US_ASCII),before=bases.clone();
		final int w=8,p=w+k/2,original="ACGT".indexOf(raw.charAt(p));
		final Audit lookup=new Audit(k);final int[] expectedS={-1,-1,-1,-1},expectedI={-1,-1,-1,-1};int expectedD=-1;
		for(int shift=-n/2;shift<=n/2;shift++){
			for(int b=0;b<4;b++){if(b!=original){
				final String edited=raw.substring(0,p)+"ACGT".charAt(b)+raw.substring(p+1);
				final int depth=lookup.add(edited.substring(w+shift,w+shift+k));
				expectedS[b]=shift==-n/2 ? depth : Math.min(expectedS[b],depth);
			}}
		}
		for(int shift=-n/2;shift<=n/2;shift++){
			final String deleted=raw.substring(0,p)+raw.substring(p+1);
			final int depth=lookup.add(deleted.substring(w+shift,w+shift+k));
			expectedD=shift==-n/2 ? depth : Math.min(expectedD,depth);
		}
		for(int shift=-n/2;shift<=n/2;shift++){
			for(int b=0;b<4;b++){
				final String inserted=raw.substring(0,p)+"ACGT".charAt(b)+raw.substring(p);
				final int d=lookup.add(inserted.substring(w+shift,w+shift+k));
				expectedI[b]=shift==-n/2 ? d : Math.min(expectedI[b],d);
			}
		}
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(k,lookup,n);
		probe.substitutions(bases,w,p);
		check(probe.queries==3*n && lookup.cursor==3*n,"Substitution stage must make exactly 3*N lookups, without eager indels.");
		check(Arrays.equals(expectedS,probe.substitutionDepth),"Substitution minimum assigned to wrong alternative.");
		probe.indels();check(probe.queries==8*n && lookup.cursor==8*n,"Fully defined interior probe must enumerate exactly 8*N kmers.");
		check(probe.deletionDepth==expectedD && Arrays.equals(expectedI,probe.insertionDepth),"Indel minimum differs from the materialized oracle.");
		probe.indels();check(probe.queries==8*n,"Repeated lazy indels must be cached.");
		check(Arrays.equals(before,bases),"Candidate enumeration must not mutate original bases.");
	}
	private static void edges(){
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){@Override public int count(final Kmer key){return 7;}};
		final byte[] b=new byte[40];Arrays.fill(b,(byte)'A');
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(31,lookup,3);
		probe.substitutions(b,0,15);probe.indels();
		check(probe.substitutionsAvailable && probe.insertionsAvailable && probe.deletionAvailable && probe.queries==24,"Shift the group right when the supplied center is at the read's left edge.");
		probe.substitutions(b,4,4);probe.indels();
		check(probe.substitutionsAvailable && probe.insertionsAvailable && probe.deletionAvailable && probe.queries==24,"An internal candidate must not disappear merely because it lies at the original window's edge.");
		b[3]='N';probe.substitutions(b,4,19);probe.indels();
		check(!probe.substitutionsAvailable && !probe.insertionsAvailable && !probe.deletionAvailable,"A single undefined flank replicate vetoes each candidate.");
		for(int n:new int[]{0,2,33}){
			boolean failed=false;try{new LocalEditKmerProbe(31,lookup,n);}catch(IllegalArgumentException e){failed=true;}
			check(failed,"N must be positive, odd and no larger than K.");
		}
	}
	private static final class Audit implements HomopolymerIndelProposal.CountLookup {
		Audit(final int k_){k=k_;}
		int add(final String s){
			final Kmer key=new Kmer(k);for(int i=0;i<s.length();i++){key.addRight(s.charAt(i));}
			expected.add(key.key().clone());return depth(expected.size()-1);
		}
		@Override public int count(final Kmer key){
			check(cursor<expected.size() && Arrays.equals(expected.get(cursor),key.key()),"Adjacent probe key differs from an independently materialized edited read.");
			return depth(cursor++);
		}
		private static int depth(final int ordinal){return (ordinal*17+13)%29;}
		final int k;int cursor;final ArrayList<long[]> expected=new ArrayList<long[]>();
	}
	private static void check(final boolean ok,final String why){checks++;if(!ok){throw new AssertionError(why);}}
	private static long checks;
}
