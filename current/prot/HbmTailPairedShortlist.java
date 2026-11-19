package prot;

import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;

/** Reusable DDL top-N plus covering5 top-N shortlist with DDL-first primitive de-duplication. */
final class HbmTailPairedShortlist {

	static HbmTailPairedShortlist load(String profilesList,String sidecarPath,String consensusPath,
			String kmersPath,int buckets,String mode,int k,long seed,boolean indexed,String namespace,
			int ddlTop,int coverTop) throws Exception{
		require(buckets>=128 && buckets<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(buckets)==1 &&
			k>=1 && k<=12 && (mode.equals("literal") || mode.equals("nonempty")),
			"Invalid DDL sketch recipe");
		final DDLFamilySearch.Index ddl=DDLFamilySearch.loadForSearch(
			Paths.get(profilesList), buckets, mode, k, seed, indexed, false, false);
		final FamilyShortlistSidecar sidecar=FamilyShortlistSidecar.load(sidecarPath, consensusPath, kmersPath);
		return new HbmTailPairedShortlist(ddl, sidecar, namespace, buckets, k, seed, ddlTop, coverTop);
	}

	HbmTailPairedShortlist(DDLFamilySearch.Index ddl_,FamilyShortlistSidecar sidecar_,String namespace_,
			int buckets_,int k_,long seed_,int ddlTop,int coverTop){
		ddl=ddl_; sidecar=sidecar_; namespace=namespace_; buckets=buckets_; k=k_; seed=seed_;
		require(sidecar.coveringK==5, "Paired shortlist requires covering5 sidecar");
		require(ddlTop>0 && coverTop>0 && ddlTop+coverTop<=1000,
			"Expected positive top counts and a bounded union buffer");
		ddlLimit=Math.min(ddlTop, ddl.families());
		coverLimit=Math.min(coverTop, sidecar.nFamilies);
		require(ddlLimit>0 && coverLimit>0, "Paired shortlist requires nonempty rosters");
		ddlToSidecar=mapDdlToSidecar(ddl, sidecar, namespace);
	}

	Scratch makeScratch(){
		return new Scratch(ddl, ddlLimit, sidecar, coverLimit, buckets);
	}

	int limit(){
		return ddlLimit+coverLimit;
	}

	void validateRoster(String[] ids){
		require(Arrays.equals(ids, sidecar.repIds), "Consensus and paired sidecar rosters must match exactly");
	}

	int rank(byte[] rawAscii,ProteinSequence query,Scratch scratch){
		DDLFamilyProfiles.sketch(rawAscii, k, seed, scratch.ddlQuery);
		ddl.score(scratch.ddlQuery, scratch.ddlWorkspace);
		scratch.ddlFound=ddl.top(scratch.ddlWorkspace, scratch.ddlBest);
		ProteinSearcher.scoreFamiliesKmers(query, sidecar, scratch.coverScratch);
		scratch.coverFound=topCovering(scratch.coverScratch.sk, scratch.coverScratch.topIdx);
		scratch.queryRegisters=nonzero(scratch.ddlQuery);
		scratch.ddlCandidates=scratch.ddlWorkspace.size;
		scratch.ddlPostings=scratch.ddlWorkspace.postings;
		scratch.coverQueryRegisters=scratch.coverScratch.qSize;
		scratch.coverCandidates=positive(scratch.coverScratch.sk);
		scratch.unionFound=union(scratch.ddlBest, scratch.ddlFound,
			scratch.coverScratch.topIdx, scratch.coverFound, scratch);
		return scratch.unionFound;
	}

	/** Selects positive scores in descending order, breaking ties by ascending roster index.
	 * The output buffer must be nonempty; entries beyond the returned size are unchanged.
	 */
	static int topCovering(int[] scores,int[] out){
		int size=0;
		for(int r=0; r<scores.length; r++){
			final int score=scores[r];
			if(score<=0 || (size==out.length && score<=scores[out[size-1]])){continue;}
			int at=Math.min(size, out.length-1);
			while(at>0 && score>scores[out[at-1]]){at--;}
			if(size<out.length){size++;}
			for(int j=size-1; j>at; j--){out[j]=out[j-1];}
			out[at]=r;
		}
		return size;
	}

	private int union(int[] ddlBest,int ddlFound,int[] coverBest,int coverFound,Scratch scratch){
		if(scratch.epoch==Integer.MAX_VALUE){Arrays.fill(scratch.seen, 0); scratch.epoch=1;}
		final int epoch=scratch.epoch++;
		int size=0;
		for(int i=0; i<ddlFound; i++){
			final int family=ddlToSidecar[ddlBest[i]];
			scratch.seen[family]=epoch; scratch.union[size++]=family;
		}
		for(int i=0; i<coverFound; i++){
			final int family=coverBest[i];
			if(scratch.seen[family]!=epoch){scratch.seen[family]=epoch; scratch.union[size++]=family;}
		}
		return size;
	}

	private static int[] mapDdlToSidecar(DDLFamilySearch.Index ddl,FamilyShortlistSidecar sidecar,String namespace){
		final HashMap<String,Integer> sidecarRanks=new HashMap<String,Integer>();
		for(int i=0; i<sidecar.repIds.length; i++){
			require(sidecarRanks.put(sidecar.repIds[i], Integer.valueOf(i))==null,
				"Duplicate sidecar family ID: "+sidecar.repIds[i]);
		}
		final int[] map=new int[ddl.names.length];
		final boolean[] seen=new boolean[sidecar.nFamilies];
		for(int i=0; i<ddl.names.length; i++){
			final String family=familyId(ddl.names[i], namespace);
			final Integer rank=sidecarRanks.get(family);
			require(rank!=null, "DDL family absent from sidecar: "+ddl.names[i]);
			require(!seen[rank.intValue()], "Duplicate DDL family after namespace mapping: "+ddl.names[i]);
			seen[rank.intValue()]=true; map[i]=rank.intValue();
		}
		for(int i=0; i<seen.length; i++){require(seen[i], "Sidecar family absent from DDL profiles: "+sidecar.repIds[i]);}
		return map;
	}

	private static String familyId(String raw,String namespace){
		if(namespace!=null && namespace.length()>0 && !namespace.equals(".")){
			final String prefix=namespace+":";
			require(raw.startsWith(prefix), "DDL family lacks namespace "+namespace+": "+raw);
			return raw.substring(prefix.length());
		}
		return raw;
	}

	private static int nonzero(char[] values){
		int count=0; for(char value : values){if(value>0){count++;}} return count;
	}

	private static int positive(int[] values){
		int count=0; for(int value : values){if(value>0){count++;}} return count;
	}

	private static void require(boolean ok,String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	static final class Scratch {
		Scratch(DDLFamilySearch.Index ddl,int ddlTop,FamilyShortlistSidecar sidecar,int coverTop,int buckets){
			ddlQuery=new char[buckets];
			ddlWorkspace=new DDLFamilySearch.Workspace(ddl.families());
			ddlBest=new int[ddlTop];
			coverScratch=new ProteinSearcher.ShortlistScratch(sidecar.dims, sidecar.nFamilies, coverTop);
			seen=new int[sidecar.nFamilies];
			union=new int[ddlTop+coverTop];
		}

		long checksum(){
			long sum=17;
			for(int i=0; i<ddlFound; i++){
				final int family=ddlBest[i];
				sum=31*sum+family;
				sum=31*sum+ddlWorkspace.scores[family];
			}
			sum=31*sum-7;
			for(int i=0; i<coverFound; i++){
				final int family=coverScratch.topIdx[i];
				sum=31*sum+family;
				sum=31*sum+coverScratch.sk[family];
			}
			sum=31*sum-11;
			for(int i=0; i<unionFound; i++){sum=31*sum+union[i];}
			return sum;
		}

		final char[] ddlQuery;
		final DDLFamilySearch.Workspace ddlWorkspace;
		final int[] ddlBest, seen, union;
		final ProteinSearcher.ShortlistScratch coverScratch;
		int ddlFound, coverFound, queryRegisters, ddlCandidates, coverQueryRegisters, coverCandidates, unionFound;
		long ddlPostings;
		int epoch=1;
	}

	final DDLFamilySearch.Index ddl;
	final FamilyShortlistSidecar sidecar;
	final String namespace;
	final int buckets, k, ddlLimit, coverLimit;
	final long seed;
	private final int[] ddlToSidecar;
}
