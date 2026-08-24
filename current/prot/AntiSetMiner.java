package prot;

import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Mines ANTI-NETWORKS: groups of protein families with complementary, near-mutually-
 * exclusive presence — "an organism must have A, B, or C to live, but not all of them"
 * (Brian, 2026-08-11). Members are individually BELOW the universal-marker prevalence
 * floor, but each group's UNION covers ~all organisms, so the group functions as a
 * disjunctive marker built from sub-threshold genes. Because members live in different
 * organisms (and different loci), their retention under fragmentation is ~independent —
 * no shared operon fate, unlike positively-correlated modules (which lost to random).
 * A native complement of ~1 per group also makes co-occurring members strong
 * contamination evidence: exclusive alternatives almost never co-occur natively.
 *
 * <p>Algorithm: greedy SET COVER with an exclusivity penalty. Universe = organisms;
 * each candidate family covers the organisms carrying it. A group grows by repeatedly
 * adding the family maximizing (newly covered organisms) - lambda*(overlap with already
 * covered), until the union prevalence target is reached. The penalty is what pushes
 * toward mutual exclusivity rather than mere complementarity. Groups are mined
 * repeatedly with used families removed, yielding many disjoint "one-of" groups.
 *
 * <p>Outputs: a groups report (members + union prevalence + overlap stats), pooled
 * subset files compatible with MagQCVectorMaker subnet=famset (each pools
 * groupsperset groups so the native total is ~groupsperset — trainable), and
 * size-matched random-control subset files drawn from the SAME prevalence band
 * (candidates minus the pool's members), so the control isolates the grouping's value
 * from the band's value.
 *
 * <p>Phylum note: phylogenetically complementary families (the Proteobacteria version
 * of essential function F vs the Bacillota version) are deliberately INCLUDED — a
 * group only reaches ~98% union prevalence if it genuinely tiles the corpus, and
 * "one of the per-clade implementations" is a valid anti-network. A within-phylum-
 * residual variant (true isozyme alternatives) is a planned v2; raw first, and the
 * random controls judge both.
 *
 * <p>Idiom rewrite (2026-08): I/O runs on {@link ByteFile}/{@link ByteStreamWriter}
 * (byte[] lines, zero per-line String allocation), field extraction on a reused
 * {@link LineParser1} (SIMD delimiter scan, byte-range parsing), collections on
 * primitive {@link IntList}/boolean[] (no boxing). {@code loadPresence}'s file scan
 * and {@code mine}'s per-candidate best-scan are both thread-per-core (see
 * {@link FillThread}, {@link ScanThread}) — the candidate scan's parallel merge is
 * ORDER-PRESERVING (shard 0..T-1, strict-greater-wins) so it reproduces the sequential
 * scan's tie-break exactly; organism-index assignment in loadPresence is order-free by
 * construction (the algorithm only ever reads bitset SET membership, never organism
 * identity), so its threaded fill needs no such care. Algorithm and output format are
 * UNCHANGED — see {@code AntiSetMinerTest} for the pinned behavioral contract.
 *
 * Usage: java prot.AntiSetMiner cache=perorg_cache.tsv familylist=familylist.tsv
 *   out=anti_groups.tsv subsetprefix=anti [domain=bacteria minprev=0.02 maxprev=0.90
 *   uniontarget=0.98 unionmin=0.95 lambda=1.0 mingain=0.005 maxmembers=30
 *   maxgroups=1000 groupsperset=50 controls=t seed=1 patience=20]
 *
 * <p>{@code patience} (added 2026-08-22, default 20 preserves prior behavior exactly): the
 * consecutive-rejected-seeds limit before {@code mine()} gives up. A seed that grows a group
 * failing to reach {@code unionmin} is marked used and never retried regardless of patience --
 * patience only controls how many DIFFERENT (lower-prevalence) seeds get a chance before the
 * miner stops looking for more groups. Tune this before lowering unionMin if groups aren't
 * forming: a low patience can exhaust the budget on high-prevalence seeds that can't reach the
 * union target, before ever trying a seed/path that could.
 *
 * @author UMP45, Eru
 */
public class AntiSetMiner implements Accumulator<AntiSetMiner.FillThread> {

	public static void main(String[] args){
		AntiSetMiner x=new AntiSetMiner(args);
		x.process();
	}

	public AntiSetMiner(String[] args){
		for(String arg : args){
			int eq=arg.indexOf('=');
			if(eq<0){continue;}
			String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("cache")){cacheFile=b;}
			else if(a.equals("familylist")){familyFile=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("subsetprefix")){subsetPrefix=b;}
			else if(a.equals("domain")){domain=b.toLowerCase();}
			else if(a.equals("minprev")){minPrev=Double.parseDouble(b);}
			else if(a.equals("maxprev")){maxPrev=Double.parseDouble(b);}
			else if(a.equals("uniontarget")){unionTarget=Double.parseDouble(b);}
			else if(a.equals("unionmin")){unionMin=Double.parseDouble(b);}
			else if(a.equals("lambda")){lambda=Double.parseDouble(b);}
			else if(a.equals("mingain")){minGain=Double.parseDouble(b);}
			else if(a.equals("maxmembers")){maxMembers=Integer.parseInt(b);}
			else if(a.equals("maxgroups")){maxGroups=Integer.parseInt(b);}
			else if(a.equals("groupsperset")){groupsPerSet=Integer.parseInt(b);}
			else if(a.equals("controls")){controls=parseBool(b);}
			else if(a.equals("exclude")){excludeFile=b;}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("patience")){patience=Integer.parseInt(b);}
			else{System.err.println("Warning: unknown arg "+arg);}
		}
		if(cacheFile==null || familyFile==null || out==null || subsetPrefix==null){
			throw new RuntimeException("Required: cache= familylist= out= subsetprefix=");
		}
	}

	void process(){
		loadFamilyCount();
		loadPresence();
		mine();
		write();
	}

	private static boolean parseBool(String s){
		return s==null || s.equals("t") || s.equals("true") || s.equals("1") || s.equals("yes");
	}

	/** Counts family rows (id line count, header excluded). Byte-level scan, no String per line. */
	private void loadFamilyCount(){
		final ByteFile bf=ByteFile.makeByteFile(familyFile, true);
		bf.nextLine();//header
		int c=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){c++;}
		}
		bf.close();
		numFam=c;
	}

	/** Loads per-organism presence bitsets from the per-org cache (famcounts col 15).
	 *  Two passes: (1) count qualifying (domain-matched) rows cheaply, sizing presence[][]/words
	 *  up front — no whole-file buffering; (2) thread-per-core fill, each thread claiming
	 *  organism indices from a shared counter and buffering (org,rank) pairs thread-locally;
	 *  a single-threaded accumulate() ORs them into the shared bitsets (no concurrent writers).
	 *  Organism-index assignment is order-free: the algorithm only ever reads bitset SET
	 *  membership (popcount/andPopcount), never which physical row a bit corresponds to. */
	private void loadPresence(){
		final byte[] domainPrefix=(domain==null) ? null :
			domain.substring(0, Math.min(4, domain.length())).getBytes();

		int n=0;
		{
			final ByteFile bf=ByteFile.makeByteFile(cacheFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()<16){continue;}
				if(domainPrefix!=null && !prefixMatchIgnoreCase(lp, 1, domainPrefix)){continue;}
				n++;
			}
			bf.close();
		}
		nOrgs=n;
		words=(nOrgs+63)>>6;
		presence=new long[numFam][];
		prevalence=new int[numFam];

		final ByteFile bf2=ByteFile.makeByteFile(cacheFile, true);
		final AtomicInteger orgCounter=new AtomicInteger(0);
		final int threads=Tools.max(1, Shared.threads());
		final ArrayList<FillThread> alpt=new ArrayList<FillThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new FillThread(bf2, orgCounter, domainPrefix, numFam));
		}
		final boolean success=ThreadWaiter.startAndWait(alpt, this);
		bf2.close();
		if(!success){throw new RuntimeException("loadPresence: a fill thread failed.");}
		assert(orgCounter.get()==nOrgs) : "Row count changed between passes: "+orgCounter.get()+" vs "+nOrgs;
		assert(nOrgs>0) : "No organisms loaded (domain filter too strict, or empty/header-only cache)?";

		for(int f=0; f<numFam; f++){
			if(presence[f]!=null){prevalence[f]=popcount(presence[f]);}
		}
		System.err.println("orgs="+nOrgs+" ("+domain+"), numFam="+numFam);
	}

	@Override
	public void accumulate(FillThread t){
		for(int i=0; i<t.pairOrg.size(); i++){
			final int o=t.pairOrg.get(i), r=t.pairRank.get(i);
			if(presence[r]==null){presence[r]=new long[words];}
			presence[r][o>>6]|=(1L<<(o&63));
		}
		errorState|=!t.success;
	}
	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return rwlock;}

	/** True if line[a, a+prefixLower.length) case-insensitively equals prefixLower (already lowercase). */
	private static boolean prefixMatchIgnoreCase(LineParser1 lp, int term, byte[] prefixLower){
		final int len=lp.length(term);
		if(len<prefixLower.length){return false;}
		final byte[] line=lp.line();
		final int a=lp.a();
		for(int i=0; i<prefixLower.length; i++){
			byte c=line[a+i];
			if(c>='A' && c<='Z'){c+=32;}
			if(c!=prefixLower[i]){return false;}
		}
		return true;
	}

	/** Parses "rank:count;rank:count;..." within [a,b) of line, adding (org,rank) pairs to the
	 *  thread-local output lists. Count is unused (matches original semantics — only presence
	 *  matters, not copy number). An out-of-range rank crashes loud (a familylist/cache version
	 *  mismatch should be visible, not silently under-counted). */
	private static void parseFamCounts(byte[] line, int a, int b, int org, int numFam,
			IntList outOrg, IntList outRank){
		int start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==';'){
				if(i>start){
					int colon=-1;
					for(int j=start; j<i; j++){if(line[j]==':'){colon=j; break;}}
					assert(colon>start) : "Malformed famcounts field (missing ':' in a rank:count pair).";
					final int r=Parse.parseInt(line, start, colon);
					assert(r<numFam) : "rank "+r+" out of bounds (numFam="+numFam+"); familylist/cache mismatch?";
					outOrg.add(org); outRank.add(r);
				}
				start=i+1;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Mining            ----------------*/
	/*--------------------------------------------------------------*/

	/** Greedy set-cover-with-exclusivity mining loop. */
	private void mine(){
		final int lo=(int)Math.ceil(minPrev*nOrgs), hi=(int)Math.floor(maxPrev*nOrgs);
		final boolean[] candidate=new boolean[numFam];
		int nCand=0;
		for(int f=0; f<numFam; f++){
			if(presence[f]!=null && prevalence[f]>=lo && prevalence[f]<=hi){candidate[f]=true; nCand++;}
		}
		//exclude= removes families (e.g. an already-trained pool) from BOTH mining and the
		//control band, so an incremental run is disjoint-by-construction from its parent.
		if(excludeFile!=null){
			int nEx=0;
			final ByteFile bf=ByteFile.makeByteFile(excludeFile, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				int a=0, b=line.length;
				while(a<b && line[a]<=' '){a++;}
				while(b>a && line[b-1]<=' '){b--;}
				if(b>a){
					final int f=Parse.parseInt(line, a, b);
					if(f<numFam && candidate[f]){candidate[f]=false; nEx++;}
				}
			}
			bf.close();
			nCand-=nEx;
			System.err.println("excluded "+nEx+" families via "+excludeFile);
		}
		System.err.println("candidates in prevalence band ["+minPrev+","+maxPrev+"]: "+nCand);
		candidateBand=candidate.clone();//preserved for control sampling

		final int scanThreads=Tools.max(1, Tools.min(Shared.threads(), (numFam+63)/64));
		final boolean[] used=new boolean[numFam];
		final boolean[] inGroup=new boolean[numFam];//S4: O(1) membership, reused across groups
		final IntList touched=new IntList();//indices set in inGroup this group, for cheap clearing
		final long[] covered=new long[words];
		final int minGainOrgs=Math.max(1, (int)Math.ceil(minGain*nOrgs));
		int consecutiveFails=0;
		while(groups.size()<maxGroups && consecutiveFails<patience){
			//Seed: highest-prevalence unused candidate.
			int seedFam=-1;
			for(int f=0; f<numFam; f++){
				if(candidate[f] && !used[f] && (seedFam<0 || prevalence[f]>prevalence[seedFam])){seedFam=f;}
			}
			if(seedFam<0){break;}
			java.util.Arrays.fill(covered, 0);
			touched.clear();
			final IntList members=new IntList();
			or(covered, presence[seedFam]);
			members.add(seedFam);
			inGroup[seedFam]=true; touched.add(seedFam);
			long overlapSum=0;
			while(members.size()<maxMembers){
				final int cov=popcount(covered);
				if(cov>=unionTarget*nOrgs){break;}
				final int[] best=findBestCandidate(candidate, used, inGroup, presence, covered, prevalence,
					lambda, scanThreads, minGainOrgs);
				if(best[0]<0){break;}
				or(covered, presence[best[0]]);
				members.add(best[0]);
				inGroup[best[0]]=true; touched.add(best[0]);
				overlapSum+=best[2];
			}
			for(int i=0; i<touched.size(); i++){inGroup[touched.get(i)]=false;}

			final double unionPrev=popcount(covered)/(double)nOrgs;
			assert(unionPrev>=0 && unionPrev<=1) : unionPrev;
			if(members.size()>=2 && unionPrev>=unionMin){
				assert(members.size()<=maxMembers);
				final Group g=new Group();
				g.members=members;
				g.unionPrev=unionPrev;
				g.meanOverlapFrac=overlapSum/(double)Math.max(1, popcount(covered));
				groups.add(g);
				for(int i=0; i<members.size(); i++){used[members.get(i)]=true;}
				consecutiveFails=0;
			}else{
				//Seed can't anchor an accepted group; retire it so the miner moves on.
				used[seedFam]=true;
				consecutiveFails++;
			}
		}
		System.err.println("mined "+groups.size()+" groups");
	}

	/** Parallel best-next-candidate scan (S5): partitions [0,numFam) into contiguous shards;
	 *  each {@link ScanThread} finds its local best via the SAME strict-greater-wins rule
	 *  scanning f ascending within its shard. A single-threaded ordered merge (shard 0..T-1,
	 *  same strict '&gt;' rule) reproduces the sequential scan's tie-break EXACTLY: for any
	 *  tie between shards, the smaller-f shard is processed first and is never overwritten
	 *  by an equal-scoring later shard — identical to what a straight f=0..numFam-1 scan
	 *  would pick. Returns [bestF, bestNew, bestOv] (bestF=-1 if nothing qualifies). */
	private static int[] findBestCandidate(boolean[] candidate, boolean[] used, boolean[] inGroup,
			long[][] presence, long[] covered, int[] prevalence, double lambda, int threads,
			int minGainOrgs){
		final int chunk=(candidate.length+threads-1)/threads;
		final ScanThread[] arr=new ScanThread[threads];
		for(int i=0; i<threads; i++){
			final int lo=i*chunk, hi=Math.min(candidate.length, lo+chunk);
			arr[i]=new ScanThread(lo, hi, candidate, used, inGroup, presence, covered, prevalence,
				lambda, minGainOrgs);
		}
		ThreadWaiter.startAndWait(java.util.Arrays.asList(arr));
		int bestF=-1; double bestScore=0; int bestNew=0, bestOv=0;
		for(ScanThread t : arr){
			if(t.localBestF>=0 && t.localBestScore>bestScore){
				bestF=t.localBestF; bestScore=t.localBestScore; bestNew=t.localBestNew; bestOv=t.localBestOv;
			}
		}
		return new int[]{bestF, bestNew, bestOv};
	}

	static final class ScanThread extends Thread {
		ScanThread(int lo_, int hi_, boolean[] candidate_, boolean[] used_, boolean[] inGroup_,
				long[][] presence_, long[] covered_, int[] prevalence_, double lambda_, int minGainOrgs_){
			lo=lo_; hi=hi_; candidate=candidate_; used=used_; inGroup=inGroup_;
			presence=presence_; covered=covered_; prevalence=prevalence_; lambda=lambda_;
			minGainOrgs=minGainOrgs_;
		}
		@Override
		public void run(){
			int bf=-1; double bs=0; int bn=0, bo=0;
			for(int f=lo; f<hi; f++){
				if(!candidate[f] || used[f] || inGroup[f]){continue;}
				final int ov=andPopcount(presence[f], covered);
				final int nw=prevalence[f]-ov;
				final double score=nw-lambda*ov;
				if(nw>=minGainOrgs && score>bs){bf=f; bs=score; bn=nw; bo=ov;}
			}
			localBestF=bf; localBestScore=bs; localBestNew=bn; localBestOv=bo;
		}
		final int lo, hi;
		final boolean[] candidate, used, inGroup;
		final long[][] presence;
		final long[] covered;
		final int[] prevalence;
		final double lambda;
		final int minGainOrgs;
		int localBestF=-1; double localBestScore=0; int localBestNew=0, localBestOv=0;
	}

	static final class FillThread extends Thread {
		FillThread(ByteFile bf_, AtomicInteger orgCounter_, byte[] domainPrefix_, int numFam_){
			bf=bf_; orgCounter=orgCounter_; domainPrefix=domainPrefix_; numFam=numFam_;
		}
		@Override
		public void run(){
			final LineParser1 lp=new LineParser1((byte)'\t');
			ListNum<byte[]> ln=bf.nextList();
			while(ln!=null && ln.size()>0){
				for(byte[] line : ln){
					if(line.length==0 || line[0]=='#'){continue;}
					lp.set(line);
					if(lp.terms()<16){continue;}
					if(domainPrefix!=null && !prefixMatchIgnoreCase(lp, 1, domainPrefix)){continue;}
					final int o=orgCounter.getAndIncrement();
					lp.setBounds(15);
					parseFamCounts(lp.line(), lp.a(), lp.b(), o, numFam, pairOrg, pairRank);
				}
				ln=bf.nextList();
			}
			success=true;
		}
		final ByteFile bf;
		final AtomicInteger orgCounter;
		final byte[] domainPrefix;
		final int numFam;
		final IntList pairOrg=new IntList();
		final IntList pairRank=new IntList();
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output             ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		//Groups report. write() is bounded (<=maxGroups*maxMembers records, typically <=30000),
		//so appendSlow's String.format-equivalent precision is used throughout rather than
		//ByteBuilder's fast append(double,decimals) - which special-cases exact whole numbers
		//(1.0 -> "1", not "1.0000") and would silently corrupt the fixed-decimal-places contract.
		final ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, true);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(4096);
		bb.append("#group\tsize\tunion_prev\tmean_overlap_frac\tmembers(rank:prev%)").nl();
		bsw.print(bb); bb.clear();
		int totalFams=0;
		for(int i=0; i<groups.size(); i++){
			final Group g=groups.get(i);
			bb.append(i).tab().append(g.members.size()).tab();
			bb.appendSlow(g.unionPrev, 4).tab();
			bb.appendSlow(g.meanOverlapFrac, 4).tab();
			for(int j=0; j<g.members.size(); j++){
				final int f=g.members.get(j);
				if(j>0){bb.semi();}
				bb.append(f).colon();
				bb.appendSlow(100.0*prevalence[f]/nOrgs, 1);
			}
			bb.nl();
			bsw.print(bb); bb.clear();
			totalFams+=g.members.size();
		}
		bsw.poisonAndWait();
		System.err.println("wrote "+groups.size()+" groups ("+totalFams+" families) to "+out);

		//Pooled subset files (groupsPerSet groups each) + size-matched random controls
		//drawn from the candidate band minus the pool's own members.
		final Random rnd=new Random(seed);
		int setIdx=0;
		for(int start=0; start<groups.size(); start+=groupsPerSet, setIdx++){
			final int end=Math.min(groups.size(), start+groupsPerSet);
			final IntList pool=new IntList();
			final boolean[] inPool=new boolean[numFam];
			for(int i=start; i<end; i++){
				final IntList mem=groups.get(i).members;
				for(int j=0; j<mem.size(); j++){
					final int f=mem.get(j);
					pool.add(f);
					inPool[f]=true;
				}
			}
			final String sf=subsetPrefix+"_pool"+setIdx+".txt";
			writeRanks(sf, pool);
			System.err.println("pool"+setIdx+": groups "+start+"-"+(end-1)+", "+pool.size()+" families -> "+sf);
			if(controls){
				final IntList bandRest=new IntList();
				for(int f=0; f<numFam; f++){
					if(candidateBand[f] && !inPool[f]){bandRest.add(f);}
				}
				final IntList ctrl=new IntList();
				int bandSize=bandRest.size();
				for(int k=0; k<pool.size() && bandSize>0; k++){
					final int idx=rnd.nextInt(bandSize);
					ctrl.add(bandRest.get(idx));
					bandSize--;
					bandRest.set(idx, bandRest.get(bandSize));//swap-to-end-and-shrink, no element shift
				}
				final String cf=subsetPrefix+"_ctrl"+setIdx+".txt";
				writeRanks(cf, ctrl);
				System.err.println("ctrl"+setIdx+": "+ctrl.size()+" band-matched random families -> "+cf);
			}
		}
	}

	private static void writeRanks(String file, IntList ranks){
		final ByteStreamWriter bsw=new ByteStreamWriter(file, true, false, true);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(Math.max(64, ranks.size()*4));
		for(int i=0; i<ranks.size(); i++){bb.append(ranks.get(i)).nl();}
		bsw.print(bb);
		bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/

	static final class Group {
		IntList members;
		double unionPrev;
		double meanOverlapFrac;
	}

	private static int popcount(long[] a){
		int c=0; for(long x : a){c+=Long.bitCount(x);} return c;
	}
	private static int andPopcount(long[] a, long[] b){
		int c=0; for(int i=0; i<a.length; i++){c+=Long.bitCount(a[i]&b[i]);} return c;
	}
	private static void or(long[] dst, long[] src){
		for(int i=0; i<dst.length; i++){dst[i]|=src[i];}
	}

	/*--------------------------------------------------------------*/

	private String cacheFile, familyFile, out, subsetPrefix, excludeFile;
	private String domain="bacteria";
	private double minPrev=0.02, maxPrev=0.90, unionTarget=0.98, unionMin=0.95;
	private double lambda=1.0, minGain=0.005;
	private int maxMembers=30, maxGroups=1000, groupsPerSet=50, patience=20;
	private boolean controls=true;
	private long seed=1;

	private int numFam, nOrgs, words;
	private long[][] presence;
	private int[] prevalence;
	private boolean[] candidateBand;
	private final ArrayList<Group> groups=new ArrayList<Group>();

	private boolean errorState=false;
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
}
