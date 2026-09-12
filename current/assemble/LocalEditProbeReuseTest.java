package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import ukmer.Kmer;

/** Cache/reuse regression against fresh probes, complementing materialized-window oracles.
 * @author Fischl */
public final class LocalEditProbeReuseTest {
	public static void main(String[] args){
		final boolean packed=Kmer.PACKED;
		try{
			for(boolean layout:new boolean[]{false,true}){
				Kmer.PACKED=layout;
				for(int requested:new int[]{31,62,63,94}){for(int n:new int[]{1,3,5}){
					final int k=new Kmer(requested).kbig;
					final byte[] bases=new byte[3*k+5];final Random random=new Random(551+k);
					for(int i=0;i<bases.length;i++){bases[i]=(byte)"ACGT".charAt(random.nextInt(4));}
					final Trace trace=new Trace();final LocalEditKmerProbe cached=new LocalEditKmerProbe(k,trace,n);
					for(int pass=0;pass<2;pass++){
						if(pass==1){bases[k+2]='N';bases[bases.length-1]='N';}
						cached.beginRead(bases);
						for(int p=0;p<bases.length;p++){if(bases[p]!='N'){
							final int w=Math.max(0,Math.min(p-k/2,bases.length-k));
							compare(cached,trace,k,n,bases,w,p,false);
						}}
						for(int p=bases.length-1;p>=0;p--){if(bases[p]!='N'){
							final int w=Math.max(0,Math.min(p-k/2,bases.length-k));
							compare(cached,trace,k,n,bases,w,p,true);
						}}
					}
					budget(k,n);
				}}
			}
			System.out.println("LOCAL_EDIT_PROBE_REUSE_TEST_OK checks="+checks);
		}finally{Kmer.PACKED=packed;}
	}
	private static void compare(final LocalEditKmerProbe cached,final Trace trace,final int k,final int n,final byte[] bases,final int w,final int p,final boolean direct){
		final Trace expectedTrace=new Trace();final LocalEditKmerProbe fresh=new LocalEditKmerProbe(k,expectedTrace,n);
		trace.keys.clear();final byte[] before=bases.clone();
		if(direct){fresh.indelsAt(bases,w,p);cached.indelsAt(bases,w,p);}
		else{fresh.substitutions(bases,w,p);fresh.indels();cached.substitutions(bases,w,p);cached.indels();}
		check(cached.queries==fresh.queries && trace.keys.size()==expectedTrace.keys.size(),"Cached query budget differs from fresh probe.");
		for(int i=0;i<trace.keys.size();i++){check(Arrays.equals(trace.keys.get(i),expectedTrace.keys.get(i)),"Cached key/query order differs from fresh probe.");}
		check(Arrays.equals(cached.substitutionDepth,fresh.substitutionDepth) && Arrays.equals(cached.insertionDepth,fresh.insertionDepth) && cached.deletionDepth==fresh.deletionDepth,"Cached candidate depths differ.");
		check(Arrays.equals(before,bases),"Window reuse must never mutate the input sequence.");
	}
	private static void budget(final int k,final int n){
		final byte[] bases=new byte[k+10];Arrays.fill(bases,(byte)'A');
		final Trace trace=new Trace();final LocalEditKmerProbe probe=new LocalEditKmerProbe(k,trace,n);
		probe.beginRead(bases);
		for(int p=k/2;p<k/2+5;p++){probe.substitutions(bases,3,p);probe.indelsAt(bases,3,p);}
		check(probe.windowBuilds()==1,"One immutable original window group should need one full build across positions/families.");
		probe.substitutions(bases,4,k/2);
		check(probe.windowBuilds()==1 && probe.windowRolls()>0,"Adjacent groups must roll, not rebuild.");
		bases[k/2]='C';probe.beginRead(bases);
		compare(probe,trace,k,n,bases,3,k/2,false);
		check(probe.windowBuilds()==2,"Beginning a new read scope must invalidate even same-array cached content.");
		trace.failAt=0;boolean failed=false;
		try{probe.substitutions(bases,3,k/2);}catch(IllegalStateException expected){failed=true;}
		check(failed,"Throwing lookup fixture must execute.");trace.failAt=-1;
		compare(probe,trace,k,n,bases,3,k/2,false);
		trace.failAt=0;failed=false;
		try{probe.indelsAt(bases,3,k/2);}catch(IllegalStateException expected){failed=true;}
		check(failed,"Deletion lookup failure fixture must execute.");trace.failAt=-1;
		compare(probe,trace,k,n,bases,3,k/2,false);
		// The direct single-window last window cannot delete; fail during insertion.
		if(n==1){
			trace.failAt=0;failed=false;
			try{probe.indelsAt(bases,10,k/2);}catch(IllegalStateException expected){failed=true;}
			check(failed,"Insertion lookup failure fixture must execute.");trace.failAt=-1;
			compare(probe,trace,k,n,bases,10,k/2,true);
		}
	}
	private static final class Trace implements HomopolymerIndelProposal.CountLookup {
		@Override public int count(final Kmer key){
			if(failAt==0){throw new IllegalStateException("Intentional lookup failure.");}if(failAt>0){failAt--;}
			final long[] value=key.key().clone();keys.add(value);
			return Arrays.hashCode(value)&1023;
		}
		final ArrayList<long[]> keys=new ArrayList<long[]>();int failAt=-1;
	}
	private static void check(final boolean ok,final String why){checks++;if(!ok){throw new AssertionError(why);}}
	private static long checks;
}
