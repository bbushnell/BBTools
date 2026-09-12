package bloom;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import assemble.LocalEditEngine;
import dna.AminoAcid;
import stream.Read;
import ukmer.Kmer;

/** Small real-count-min-table correction/geometry regression and load diagnostic.
 * Test-only exact string keys deliberately do not approximate the Bloom oracle.
 * The sweep reports adverse results; it does not assert noisy tables are safe.
 * @author Fischl */
public final class BloomLocalEditTest {
	public static void main(final String[] args){
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		try{
			if(args.length==0){
				for(int k:new int[]{31,62,63}){for(boolean dual:new boolean[]{false,true}){run(k,2000003,dual,6,true);}}
				System.out.println("BLOOM_LOCAL_EDIT_TEST_OK checks="+checks);
			}else{
				if(args.length!=2){throw new IllegalArgumentException("Load diagnostic: K cells, or no arguments for unit gates.");}
				run(Integer.parseInt(args[0]),Long.parseLong(args[1]),true,96,false);
			}
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
	}
	private static void run(final int k,final long cells,final boolean dual,final int samples,final boolean strict){
		KCountArray7MTA.setSeed(20260911);BloomFilter.OVERRIDE_CELLS=cells;
		final BloomFilter filter=new BloomFilter(null,null,null,k,k,4,3,1,true,false,false,0.01f,dual);
		final Exact counts=new Exact();final ArrayList<String> truth=new ArrayList<String>(),queries=new ArrayList<String>();
		final Random random=new Random(2026091102L);
		for(int sample=0;sample<samples;sample++){
			final StringBuilder s=new StringBuilder();for(int i=0;i<4*k;i++){s.append("ACGT".charAt(random.nextInt(4)));}
			final int p=2*k;
			if((sample&1)==0){for(int i=p-3;i<p+3;i++){s.setCharAt(i,'A');}}
			final String ref=s.toString();add(filter,counts,ref,12);
			// Three error classes and three untouched controls, including supported
			// real +/-3 alleles in a higher-depth reference background.
			final int type=sample%6;
			String q=ref,expected=ref;
			if(type==0){q=ref.substring(0,p)+(ref.charAt(p)=='A' ? 'C' : 'A')+ref.substring(p+1);}
			else if(type==1){q=ref.substring(0,p)+ref.substring(p+1);}
			else if(type==2){q=ref.substring(0,p)+ref.charAt(p)+ref.substring(p);}
			else if(type==4){q=ref.substring(0,p)+"ACG"+ref.substring(p);expected=q;add(filter,counts,q,6);}
			else if(type==5){q=ref.substring(0,p)+ref.substring(p+3);expected=q;add(filter,counts,q,6);}
			add(filter,counts,q,1);queries.add(q);truth.add(expected);
		}
		// Same deterministic key population at every capacity. Never insert a
		// second copy merely because the query is evaluated at multiple N values.
		if(!strict){for(int i=0;i<10000;i++){
			final StringBuilder s=new StringBuilder();for(int j=0;j<k;j++){s.append("ACGT".charAt(random.nextInt(4)));}
			add(filter,counts,s.toString(),i%5==0 ? 8 : 1);
		}}
		System.out.println("TABLE k="+k+" cells="+filter.filter.cells+" dual="+dual+" unique_exact="+counts.map.size()+
				" used1="+filter.filter.usedFraction(1)+" used2="+filter.filter.usedFraction(2)+" used4="+filter.filter.usedFraction(4));
		for(int n:new int[]{1,3,5}){
			final LocalEditEngine bloom=BloomFilterCorrectorWrapper.makeLocalEditEngine(filter,n);
			final LocalEditEngine exact=new LocalEditEngine(k,counts,n);
			int repaired=0,wrong=0,abstained=0,controlsChanged=0,exactRepaired=0,disagree=0;
			for(int i=0;i<queries.size();i++){
				final String q=queries.get(i),expected=truth.get(i);final boolean control=q.equals(expected);
				final Read read=read(q),reference=read(q);final byte[] old=read.bases,oldQ=read.quality,before=old.clone(),beforeQ=oldQ.clone();
				final int edits=bloom.correct(read,1);exact.correct(reference,1);
				final boolean correct=Arrays.equals(read.bases,expected.getBytes(StandardCharsets.US_ASCII));
				if(control){if(edits>0){controlsChanged++;}}
				else if(edits==0){abstained++;}else if(correct){repaired++;}else{wrong++;}
				if(!control && Arrays.equals(reference.bases,expected.getBytes(StandardCharsets.US_ASCII))){exactRepaired++;}
				if(!Arrays.equals(read.bases,reference.bases)){disagree++;}
				check(Arrays.equals(old,before) && Arrays.equals(oldQ,beforeQ),"Correction must not mutate the original arrays.");
				check(read.quality.length==read.bases.length && read.id.equals("fixture"),"Output quality/header contract broken.");
				final Read reverse=read(new String(AminoAcid.reverseComplementBases(before),StandardCharsets.US_ASCII));
				for(int j=0;j<beforeQ.length;j++){reverse.quality[j]=beforeQ[beforeQ.length-1-j];}
				final int reverseEdits=bloom.correct(reverse,1);
				check(reverseEdits==edits && Arrays.equals(AminoAcid.reverseComplementBases(read.bases),reverse.bases),"Bloom correction must commute with reverse complementation.");
				for(int j=0;j<read.quality.length;j++){check(read.quality[j]==reverse.quality[read.quality.length-1-j],"Retained/edit quality placement must commute with reverse complementation.");}
			}
			System.out.println("RESULT k="+k+" n="+n+" repaired="+repaired+" wrong="+wrong+" abstained="+abstained+
					" controls_changed="+controlsChanged+" exact_repaired="+exactRepaired+" disagreed="+disagree+
					" profile_queries="+bloom.profileQueries+" probe_queries="+bloom.probeQueries+" verification_queries="+bloom.verificationQueries);
			if(strict){check(repaired==samples/2 && wrong==0 && controlsChanged==0 && disagree==0,"Low-load Bloom oracle must repair all three error classes and preserve all three controls.");}
		}
	}
	private static Read read(final String s){
		final byte[] q=new byte[s.length()];for(int i=0;i<q.length;i++){q[i]=(byte)(i%41);}
		return new Read(s.getBytes(StandardCharsets.US_ASCII),q,"fixture",0,false);
	}
	private static void add(final BloomFilter filter,final Exact exact,final String s,final int copies){
		assert(copies>0) : "The fixture loader adds known positive multiplicities to both independent count stores.";
		final Kmer key=new Kmer(filter.k);
		for(int i=0;i<s.length();i++){
			key.addRight(s.charAt(i));if(key.len()<filter.k){continue;}
			final String identity=Arrays.toString(key.key());final Integer previous=exact.map.get(identity);
			exact.map.put(identity,(previous==null ? 0 : previous)+copies);
			final long first=filter.k<=31 ? filter.toKey(key.array1()[0],key.array2()[0]) : key.xor();
			if(filter.dualHash){filter.filter.increment(first,filter.k<=31 ? KCountArray.secondaryHash(first) : key.xor2(),copies);}
			else{for(int j=0;j<copies;j++){filter.filter.increment(first);}}
		}
	}
	private static final class Exact implements LocalEditEngine.CountLookup {
		@Override public int count(final Kmer key){final Integer v=map.get(Arrays.toString(key.key()));return v==null ? 0 : v;}
		final HashMap<String,Integer> map=new HashMap<String,Integer>();
	}
	private static void check(final boolean ok,final String why){checks++;if(!ok){throw new AssertionError(why);}}
	private static long checks;
}
