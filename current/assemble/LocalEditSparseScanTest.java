package assemble;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import ukmer.Kmer;

/** Independent dense-profile oracle for sampling, refinement and query reuse.
 * String maps are test-only; no boxed collections in the production scan.
 * @author Fischl */
public final class LocalEditSparseScanTest {
	public static void main(final String[] args){
		final boolean oldMask=Kmer.MASK_CORE,oldPacked=Kmer.PACKED;
		try{
			Kmer.MASK_CORE=false;Kmer.PACKED=true;
			for(int k:new int[]{31,62}){
				for(int phase=0;phase<8;phase++){
					for(int width:new int[]{1,2,3,4,5,6,7,k-1,k,k+1,3*k}){
						geometry(k,16+phase,width,false,false);
					}
					geometry(k,phase,k,true,false);geometry(k,16+phase,k,false,true);
				}
				geometry(k,0,k,false,false);geometry(k,8*k-k,k,false,false);
				for(boolean reverse:new boolean[]{false,true}){correction(k,reverse);}
			}
			for(int k:new int[]{31,62,95,127}){for(int phase=0;phase<8;phase++){sevenBaseIsland(k,phase);}}
			fallbacks();invalidBoundary();
			System.out.println("LOCAL_EDIT_SPARSE_SCAN_TEST_OK checks="+checks);
		}finally{Kmer.MASK_CORE=oldMask;Kmer.PACKED=oldPacked;}
	}
	/** Brian's explicit K-low/7-high/K-low example, including an island entirely
	 * between samples and long kmers spanning more than two packed words. */
	private static void sevenBaseIsland(final int k,final int phase){
		final byte[] bases=sequence(7*k+31,421+k);final int n=bases.length-k+1;
		final int left=16+phase,island=left+k,right=island+7;
		final int[] values=new int[n];Arrays.fill(values,12);
		Arrays.fill(values,left,island,1);Arrays.fill(values,right,right+k,2);
		final ProfileOracle oracle=new ProfileOracle(k,bases,values);
		final LocalEditDepthProfile sparse=new LocalEditDepthProfile(k,oracle);sparse.fill(bases,8);
		final LocalEditTroughLocator locator=new LocalEditTroughLocator(k,3);locator.reset(bases,sparse.depths,true);
		for(int start:new int[]{left,right}){
			check(locator.next(),"Both K-wide low runs must be detected even when the seven-position high island is unsampled.");
			check(locator.depthStart==start && locator.depthEnd==start+k,"Refinement must return the two exact cliffs, not one merged 2K+7 region.");
			check(sparse.depths.get(start-1)==12 && sparse.depths.get(start+k)==12,"Each error region has actual measured high flanks.");
			for(int j=start;j<start+k;j++){check(sparse.depths.get(j)==values[j] && oracle.visits[j]==1,"Every low-region depth must be quantified once, including unsampled interiors.");}
		}
		check(!locator.next(),"The high island must not create a third candidate or merge the two errors.");
		check(sparse.depths.get(island)==12 && sparse.depths.get(right-1)==12,"Both sides of an unsampled high island must be measured.");
		for(int j=0;j<n;j++){check(oracle.visits[j]<=1,"Boundary refinement must reuse previously measured windows.");}
	}
	private static void geometry(final int k,final int start,final int width,final boolean withN,final boolean neighbor){
		final byte[] bases=sequence(9*k-1,421+k);final int n=bases.length-k+1;
		if(withN){bases[Math.min(bases.length-1,start+k/2)]='N';}
		final int[] values=new int[n];Arrays.fill(values,12);
		for(int i=start;i<Math.min(n,start+width);i++){values[i]=i%3;}
		if(neighbor){for(int i=start+width+1;i<Math.min(n,start+2*width+1);i++){values[i]=1;}}
		final ProfileOracle oracle=new ProfileOracle(k,bases,values);
		final LocalEditDepthProfile dense=new LocalEditDepthProfile(k,oracle),sparse=new LocalEditDepthProfile(k,oracle);
		dense.fill(bases,1);oracle.resetVisits();sparse.fill(bases,8);
		final boolean[] expectedKnown=new boolean[n];
		for(int i=0;i<n;i++){expectedKnown[i]=oracle.invalid[i] || i%8==0 || i==n-1;}
		for(int i=0;i<n;){
			if(oracle.invalid[i] || values[i]>=3){i++;continue;}
			final int a=i;boolean sampled=false;
			while(i<n && !oracle.invalid[i] && values[i]<3){sampled|=(i%8==0 || i==n-1);i++;}
			if(sampled){for(int p=Math.max(0,a-1);p<=Math.min(n-1,i);p++){expectedKnown[p]=true;}}
		}
		long expectedQueries=0;
		for(int i=0;i<n;i++){
			final int expected=oracle.invalid[i] ? LocalEditDepthProfile.INVALID_N : expectedKnown[i] ? values[i] : LocalEditDepthProfile.UNKNOWN;
			check(sparse.depths.get(i)==expected,"Sampled components must match independently derived exact boundaries, depths and flanks.");
			check(dense.depths.get(i)==(oracle.invalid[i] ? 0 : values[i]),"Dense mode preserves original N-zero and exact-count semantics.");
			check(oracle.visits[i]==(expectedKnown[i] && !oracle.invalid[i] ? 1 : 0),"Each queried window is looked up once; unknown and invalid windows never query.");
			if(expectedKnown[i] && !oracle.invalid[i]){expectedQueries++;}
		}
		check(sparse.queries==expectedQueries && sparse.queries<=dense.queries,"Profile counters must equal actual queries, never exceed dense fill.");
		final LocalEditTroughLocator locator=new LocalEditTroughLocator(k,3);locator.reset(bases,sparse.depths,true);
		while(locator.next()){
			check(sparse.depths.get(locator.depthStart-1)>=3 && sparse.depths.get(locator.depthEnd)>=3,"Returned sparse troughs need actual supported flanks.");
			for(int i=locator.depthStart;i<locator.depthEnd;i++){check(sparse.depths.get(i)>=0 && sparse.depths.get(i)<3,"Unknown/N values may not enter correction candidates.");}
		}
		oracle.resetVisits();sparse.fill(bases,8);
		check(sparse.queries==expectedQueries,"Reusing a worker profile cannot retain stale depth/query state.");
	}
	private static void correction(final int k,final boolean reverse){
		final byte[] truth=sequence(16*k,100+k);
		final HashMap<String,Integer> table=new HashMap<String,Integer>();final Kmer key=new Kmer(k);
		for(byte b:truth){key.addRight(b);if(key.len()>=k){table.put(Arrays.toString(key.key()),12);}}
		final LocalEditEngine.CountLookup lookup=new LocalEditEngine.CountLookup(){
			@Override public int count(final Kmer candidate){final Integer d=table.get(Arrays.toString(candidate.key()));return d==null ? 0 : d;}
		};
		String noisy=new String(truth,java.nio.charset.StandardCharsets.US_ASCII);
		noisy=noisy.substring(0,11*k)+noisy.substring(11*k+1);
		noisy=noisy.substring(0,7*k)+"A"+noisy.substring(7*k);
		final byte[] mutated=noisy.getBytes(java.nio.charset.StandardCharsets.US_ASCII);mutated[3*k]=mutated[3*k]=='A' ? (byte)'C' : (byte)'A';
		final byte[] input=reverse ? AminoAcid.reverseComplementBases(mutated) : mutated;
		for(boolean pairs:new boolean[]{false,true}){
			for(int cap:new int[]{1,8}){
				final Read a=read(input),b=read(input);final byte[] original=b.bases,quality=b.quality;
				final LocalEditEngine dense=new LocalEditEngine(k,lookup,1,false,1),sparse=new LocalEditEngine(k,lookup,1,false,8);
				final int ea=dense.correct(a,cap,pairs),eb=sparse.correct(b,cap,pairs);
				check(ea==eb && Arrays.equals(a.bases,b.bases) && Arrays.equals(a.quality,b.quality),"Dense/sparse must agree for separated single-base errors, both strands and caps.");
				if(cap==8){check(eb==3 && Arrays.equals(b.bases,reverse ? AminoAcid.reverseComplementBases(truth) : truth),"S/I/D correction must recover independent original truth.");}
				else{check(eb==0 && b.bases==original && b.quality==quality,"Over-budget rejection preserves original arrays.");}
				if(pairs){check(dense.profileQueries==sparse.profileQueries,"Pair witnesses force a dense profile even when stride8 was requested.");}
				else{check(sparse.profileQueries<dense.profileQueries,"Sparse correction must actually reduce table accesses on the controlled read.");}
			}
		}
	}
	private static void fallbacks(){
		for(int k:new int[]{5,8,31}){
			final LocalEditDepthProfile p=new LocalEditDepthProfile(k,new HomopolymerIndelProposal.CountLookup(){@Override public int count(final Kmer key){return 12;}});
			p.fill(sequence(k==31 ? k+4 : 10*k,7),8);
			check(!p.sparse && p.stride==1 && p.queries==p.depths.size,"Small K and short profiles must fall back to dense sampling.");
		}
	}
	private static void invalidBoundary(){
		final byte[] bases=sequence(93,99);final structures.IntList counts=new structures.IntList();
		for(int i=0;i<63;i++){counts.add(LocalEditDepthProfile.UNKNOWN);}counts.set(32,0);
		final LocalEditTroughLocator locator=new LocalEditTroughLocator(31,3);locator.reset(bases,counts,true);
		boolean threw=false;try{locator.next();}catch(IllegalStateException expected){threw=true;}
		check(threw,"An unexpanded sparse low region cannot treat UNKNOWN as a supported flank.");
		locator.reset(bases,counts);threw=false;try{locator.next();}catch(IllegalArgumentException expected){threw=true;}
		check(threw,"Dense consumers reject sparse sentinels instead of silently changing semantics.");
	}
	private static Read read(final byte[] bases){final byte[] q=new byte[bases.length];Arrays.fill(q,(byte)30);return new Read(bases.clone(),q,"sparse",0,false);}
	private static byte[] sequence(final int length,final long seed){final Random r=new Random(seed);final byte[] b=new byte[length],a={'A','C','G','T'};for(int i=0;i<length;i++){b[i]=a[r.nextInt(4)];}return b;}
	private static final class ProfileOracle implements HomopolymerIndelProposal.CountLookup {
		ProfileOracle(final int k,final byte[] bases,final int[] depths_){
			depths=depths_;visits=new int[depths.length];invalid=new boolean[depths.length];
			for(int start=0;start<depths.length;start++){
				final Kmer key=new Kmer(k);
				for(int j=start;j<start+k;j++){if(bases[j]=='N'){invalid[start]=true;}key.addRight(bases[j]);}
				if(!invalid[start]){check(starts.put(Arrays.toString(key.key()),start)==null,"Oracle fixture requires unique canonical kmers to count visits per position.");}
			}
		}
		@Override public int count(final Kmer key){final Integer start=starts.get(Arrays.toString(key.key()));if(start==null){throw new AssertionError("Refinement reconstructed a kmer absent from the original input.");}visits[start]++;return depths[start];}
		void resetVisits(){Arrays.fill(visits,0);}
		final HashMap<String,Integer> starts=new HashMap<String,Integer>();
		final int[] depths,visits;final boolean[] invalid;
	}
	private static void check(final boolean pass,final String why){checks++;if(!pass){throw new AssertionError(why);}}
	private static long checks;
}
