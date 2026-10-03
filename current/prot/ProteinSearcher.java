package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;
import map.LongHashSet;
import map.LongLongListHashMap;
import structures.IntList;
import structures.LongList;

/**
 * In-memory, blastp-style protein similarity search (Phase-1 MVP).
 *
 * <p><b>This is the load-bearing API.</b> {@link #search(List, List)} takes
 * query and database protein sequences as in-memory objects and returns a list
 * of {@link ProteinHit} objects — no disk round-trip. It is callable from other
 * BBTools code (e.g. so a binner could later assess an in-memory bin against a
 * marker set).</p>
 *
 * <p>Algorithm (basic version, correctness-first):</p>
 * <ul>
 * <li><b>K-mer seed:</b> exact (or amino8-reduced) k-mer overlap selects
 * candidate targets for each query — a target is a candidate when it shares at
 * least {@link #minSeedHits} distinct k-mers with the query.
 * <li><b>Extend/score:</b> each candidate is aligned with BLOSUM62 affine-gap
 * local Smith-Waterman ({@link AAAligner}), yielding the single best HSP per
 * query-target pair.
 * <li><b>Statistics:</b> rigorous gapped-BLOSUM62 bitscore, plus an approximate
 * Karlin-Altschul E-value (no edge-length correction; flagged approximate).
 * <li><b>Filter + sort:</b> E-value cutoff and optional min raw score / min
 * pident; {@code maxTargetSeqs} caps distinct targets by best-HSP score; output
 * is sorted by the frozen total order (query, E-value, bitscore, target, tstart,
 * qstart) so results are deterministic.
 * </ul>
 *
 * <p>Multithreading (added 2026-09-02, PATH_TO_PRODUCTION_v1 A2): {@link #threads} parallelizes
 * over queries with a bounded worker pool; output is byte-identical to the single-threaded path
 * for any thread count (see {@link #threads}'s javadoc for why).</p>
 *
 * <p>Deferred (reported, not built): rigorous edge-corrected E-values, an indexed
 * fast-search mode (the current index only prunes candidates, not full alignment),
 * multi-HSP per pair, and blastx / six-frame translated search.</p>
 *
 * @author Eru
 */
public final class ProteinSearcher {

	/** K-mer length for seeding. */
	public int k=5;
	/** Minimum distinct shared k-mers for a target to be a candidate. */
	public int minSeedHits=1;
	/** E-value significance cutoff (frozen default 10). */
	public double evalueCutoff=10.0;
	/** Minimum raw score to report (0 = disabled). */
	public int minRawScore=0;
	/** Minimum percent identity to report (0 = disabled). */
	public double minPident=0;
	/**
	 * Minimum BIDIRECTIONAL coverage to report (0 = disabled): the alignment must cover at
	 * least this fraction of BOTH the query and the target length, i.e.
	 * {@code alnLen/qLen >= minCoverage AND alnLen/tLen >= minCoverage} -- equivalently
	 * {@code alnLen >= minCoverage*max(qLen,tLen)}. This is mmseqs' {@code --cov-mode 0}
	 * (the default when only {@code -c} is given, no explicit --cov-mode). Matching this
	 * matters beyond mirroring mmseqs specifically: without it, ProteinSearcher accepts
	 * short, statistically weak local alignments (a shared domain/motif) as if they were
	 * true family membership. Verified on real data (Eru+UMP45, 2026-08-29, ecoli.fa vs
	 * consensus_reps_v4b.fasta): with minCoverage=0, 70% of hits (6340/9112) had under 20%
	 * mutual coverage, several with e-values of 1.8-6.7 (statistically expected by chance),
	 * while only 13% (1198/9112) passed a true bidirectional >=0.7 test -- consistent with
	 * mmseqs' independently-measured ~952 real gene-family assignments on the same input.
	 * Default 0 preserves prior behavior for callers that haven't opted in; family-assignment
	 * callers (feeding MAG-QC vectors) MUST set 0.7 to match the training distribution.
	 */
	public double minCoverage=0;
	/** Cap on distinct targets per query, ranked by best-HSP score. */
	public int maxTargetSeqs=Integer.MAX_VALUE;
	/** Use the amino8 reduced alphabet for seeding (more sensitive). */
	public boolean reducedSeed=false;
	/**
	 * Worker threads for {@link #search(List, List)} (default 1: single-threaded, no thread
	 * overhead, and the original unthreaded code path). Parallelism is over QUERIES only
	 * (embarrassingly parallel -- each query's hit set depends only on that query and the
	 * shared, read-only target index). Bounded: exactly this many worker threads are created
	 * regardless of query count (never one thread per query). Deterministic: each query writes
	 * its own hit list to its own array slot (no shared mutable state between queries other than
	 * the read-only k-mer index and AAAligner's own thread-local scratch), and the final sort by
	 * {@link #TOTAL_ORDER} is applied to the complete set exactly as in the single-threaded path
	 * -- so output is byte-identical for any thread count. Added 2026-09-02 (PATH_TO_PRODUCTION_v1
	 * A2).
	 */
	public int threads=1;
	/**
	 * Shortlist size (Brian, 2026-09-03, relayed by UMP45: "integrate the hybrid shortlist into
	 * the assigner" -- the assigner is the current bottleneck at 10 genes/s per 32 threads).
	 * Default 0: current full behavior, UNCHANGED and byte-identical -- every seed-passing
	 * candidate is aligned (see {@link #searchOneQuery}). {@code shortlist>0}: {@link #sidecar}
	 * must be set; instead of the k-mer seed filter, every query is scored against every family
	 * by F4 (z(s_k) + 2*z(s_c) - d_L, the production-winning fusion from HybridShortlistAssay's
	 * evaluation, {@code results/hybrid_shortlist_v1.md}) using the sidecar's precomputed
	 * composition/covering-set/length data, and only the top {@code shortlist} families by F4
	 * (via the deterministic scalar {@link #topN(float[], int, int[])}) are aligned. The L-hard length
	 * window is OFF by default (soft d_L only) pending Brian's cutoff decision -- see
	 * {@code results/hybrid_shortlist_v1.md} §2 ("soft length consistently beats hard filtering").
	 */
	public int shortlist=0;
	/** Aligner used by {@link #searchOneQueryShortlist}: "aaaligner" retains the legacy default;
	 *  "d55" uses the canonical {@link GlocalAminoLinear} recurrence. The historical "glocal"
	 *  and "blosum" names both select {@link AAAligner#alignGlocal}; they do not invoke the
	 *  archived GlocalAmino/GlocalAminoBlosum prototypes. Only consulted when {@link #shortlist}
	 *  is positive; the unshortlisted legacy search still uses AAAligner directly. */
	public String aligner="aaaligner";
	/**
	 * K-mer term normalization inside F4 (UMP45, 2026-09-09, plans/SHORTLIST_KMER_NORMALIZATION_ASSAY_v1.md;
	 * Brian's idea "kmer matches/(ref length + 100) or similar ... to compensate for coincidental matches").
	 * Only consulted when {@link #shortlist}&gt;0. {@code "none"} (default): today's path, BYTE-IDENTICAL
	 * (the same arithmetic on the same ints). {@code "sub"}: the phase-1 winner on the held-out harness
	 * ({@code results/shortlist_kmer_normalization_v1.md}): {@code s_k} is replaced by
	 * {@code s_k - rho_q*|C_f|}, where {@code |C_f|} is family f's covering-set size ({@link
	 * FamilyShortlistSidecar#coveringSetSize}) and {@code rho_q = sum_f s_k / sum_f |C_f|} is this
	 * query's own coincidental-hit rate over all families -- the expected coincidental count is
	 * subtracted before z-scoring. {@code "cs300"}: Brian's ratio form, {@code s_k/(|C_f|+300)} (the
	 * phase-1 top-1 maximizer, carried per Yoimiya 16:55Z so the trade-off is measured on the production
	 * path). Any other value crashes loud in {@link #search}.
	 */
	public String kmerNorm=KMERNORM_NONE;
	public static final String KMERNORM_NONE="none", KMERNORM_SUB="sub", KMERNORM_CS300="cs300";
	private static boolean validKmerNorm(final String v){
		return v.equals(KMERNORM_NONE) || v.equals(KMERNORM_SUB) || v.equals(KMERNORM_CS300);
	}
	/** Precomputed per-family shortlist data ({@link FamilyShortlistSidecarBuilder}'s output,
	 *  loaded via {@link FamilyShortlistSidecar#load}). Required when {@link #shortlist}&gt;0. */
	public FamilyShortlistSidecar sidecar=null;
	/**
	 * A3 (2026-09-04, plans/PER_FAMILY_THRESHOLDS_v1.md): apply the v3 sidecar's per-family
	 * thresholds. Only consulted when {@link #shortlist}&gt;0. Default false: current behavior,
	 * BYTE-IDENTICAL (this flag adds no new code to the false path at all -- see {@link
	 * #searchOneQueryShortlist}). {@code true}: stage 1 (length-ratio window {@code r_lo<=r<=r_hi}),
	 * stage 2 (dimer-cosine floor {@code s_c>=minDimer}), and stage 3 (covering-kmer floor
	 * {@code s_k>=minKmers}) form a survivor mask applied BEFORE {@link #topN(float[], int, int[])} -- {@link
	 * #shortlist} candidates are selected from the mask, never padded with a rejected family when
	 * fewer survive (sec 1: "fewer only when fewer survive"); stage 4 (identity/score-ratio/raw-
	 * score/coverage, sec 1 4a-4d) applies each family's own APPLIED threshold (already max/min-
	 * combined with the global floor at sidecar-build time, sec 2) inside the hit filter, alongside
	 * (not instead of) {@link #minRawScore}/{@link #minPident}/{@link #minCoverage}. Requires a v3
	 * sidecar ({@link FamilyShortlistSidecar#schemaVersion}&gt;=2); crashes loud on a v2 sidecar
	 * (sec 3: "crash-loud on a sidecar v2 when triage=t").
	 */
	public boolean triage=false;
	/**
	 * Stage 4's coverage definition when {@link #triage} is true (sec 1's own open question --
	 * "Brian picks the definition with the cutoffs"; kept behind this flag "so today's behaviour
	 * remains reproducible" per sec 3, literally). {@code "column"} (default): {@code
	 * alnLen/qLen}, {@code alnLen/tLen} -- matches {@link #minCoverage}'s existing formula exactly
	 * (alnLen counts gap columns). {@code "span"}: {@code (qStop-qStart+1)/qLen},
	 * {@code (tStop-tStart+1)/tLen} -- the aligned SPAN per side, no gap-column over-credit (sec 1's
	 * NOTE on the current filter's known over-credit issue).
	 */
	public String covDef="column";
	/** Per-stage triage counters (A3, sec 3: "Per-stage counters... summed into .meta"), summed
	 *  across every query and every worker thread. Only incremented when {@link #triage} is true;
	 *  read by {@link ProteinSearch#writeSidecar} for the run manifest. */
	public final java.util.concurrent.atomic.AtomicLong triageStage1Survivors=new java.util.concurrent.atomic.AtomicLong(),
		triageStage2Survivors=new java.util.concurrent.atomic.AtomicLong(),
		triageStage3Survivors=new java.util.concurrent.atomic.AtomicLong(),
		triageAligned=new java.util.concurrent.atomic.AtomicLong(),
		triageHits=new java.util.concurrent.atomic.AtomicLong();
	/**
	 * True: the emitted E-values omit the KA edge-length correction and are
	 * therefore approximate (bitscores remain rigorous). Fixed true for the MVP.
	 */
	public static final boolean EVALUE_APPROXIMATE=true;

	/** Maps extended-alphabet codes (0-19) to amino8 reduced groups (0-6). */
	private static final int[] EXT_TO_AMINO8=buildReducedMap();

	private static int[] buildReducedMap(){
		final int[] map=new int[Blosum62.X_CODE+1];
		java.util.Arrays.fill(map, -1);
		for(char c='A'; c<='Z'; c++){
			final int ext=AminoAcid.acidToNumber[c];
			if(ext>=0 && ext<=19){map[ext]=AminoAcid.acidToNumber8[c];}
		}
		return map;
	}

	/**
	 * Searches every query against the database and returns all passing hits,
	 * sorted in the frozen deterministic total order.
	 *
	 * @param queries In-memory query protein sequences.
	 * @param targets In-memory database protein sequences.
	 * @return Hits passing the significance filters, deterministically sorted.
	 */
	public List<ProteinHit> search(final List<ProteinSequence> queries,
			final List<ProteinSequence> targets){
		if(queries==null || targets==null){
			throw new RuntimeException("Null query or target list.");
		}
		checkDuplicateIds(queries, "query");
		checkDuplicateIds(targets, "database");

		long totalDbResiduesSum=0;
		for(ProteinSequence t : targets){totalDbResiduesSum+=t.length();}
		final long totalDbResidues=totalDbResiduesSum;

		if(shortlist>0 && sidecar==null){
			throw new RuntimeException("shortlist>0 requires sidecar to be set (FamilyShortlistSidecar.load(...)).");
		}
		if(shortlist>0){
			if(!validKmerNorm(kmerNorm)){
				throw new RuntimeException("Unknown kmernorm: '"+kmerNorm+"' (expected "+KMERNORM_NONE+", "+KMERNORM_SUB+" or "+KMERNORM_CS300+").");
			}
			if(sidecar.nFamilies!=targets.size()){
				throw new RuntimeException("sidecar family count ("+sidecar.nFamilies+") != targets.size() ("+
					targets.size()+") -- sidecar was built for a different target list than the one being searched.");
			}
			//Count equality alone is not enough (UMP45's review, 2026-09-03): a different
			//same-sized target list would silently mis-score every family, since sidecar row r's
			//composition/covering-set/length data is only valid for targets.get(r). Cheap, done
			//once per search() call, not per query: every target id must match the sidecar's rep_id
			//at the same index.
			for(int i=0; i<targets.size(); i++){
				final String targetId=targets.get(i).id, sidecarId=sidecar.repIds[i];
				if(!targetId.equals(sidecarId)){
					throw new RuntimeException("Target identity mismatch at index "+i+": targets has '"+
						targetId+"' but sidecar expects '"+sidecarId+"' -- the target list does not match "+
						"the family order the sidecar was built from.");
				}
			}
		}
		//Build the target k-mer index: kmer -> distinct target indices. Shared read-only across
		//all queries and all worker threads -- never mutated after this point. Skipped entirely in
		//shortlist mode (shortlist>0): that path uses the sidecar's covering-set index instead, and
		//this seed index would otherwise be pure wasted work on the assigner's hot path.
		final LongLongListHashMap index=(shortlist>0 ? null : buildIndex(targets));

		final int nQ=queries.size();
		final int nThreads=Math.max(1, threads);
		@SuppressWarnings("unchecked")
		final List<ProteinHit>[] perQuery=new List[nQ];

		if(nThreads==1 || nQ<=1){
			//Single-threaded path (also used for a trivial query count): identical to the
			//pre-threading implementation, byte-for-byte -- kept separate so the zero-thread-
			//overhead case never touches Thread/AtomicInteger machinery.
			final int[] hitCount=new int[targets.size()];
			final IntList touched=new IntList();
			final ShortlistScratch shortlistScratch=(shortlist>0 ? new ShortlistScratch(sidecar.dims, targets.size(), shortlist) : null);
			for(int qi=0; qi<nQ; qi++){
				perQuery[qi]=(shortlist>0 ?
					searchOneQueryShortlist(queries.get(qi), targets, totalDbResidues, shortlistScratch) :
					searchOneQuery(queries.get(qi), targets, totalDbResidues, index, hitCount, touched));
			}
		}else{
			//Bounded worker pool (exactly nThreads threads, regardless of query count), work-
			//stealing over query indices via a shared counter. Each thread gets its OWN
			//hitCount[]/touched scratch (these are per-query working state, not shared state --
			//giving each thread its own copy is the only change needed for correctness, since a
			//single query's hit set never depends on any other query). Each query writes to its
			//own perQuery[qi] slot, so there is no write contention; Thread.join() below
			//establishes the happens-before needed to safely read every slot afterward.
			final java.util.concurrent.atomic.AtomicInteger next=new java.util.concurrent.atomic.AtomicInteger(0);
			final Thread[] workers=new Thread[nThreads];
			final Throwable[] errors=new Throwable[nThreads];
			for(int w=0; w<nThreads; w++){
				final int wi=w;
				workers[w]=new Thread(){
					@Override
					public void run(){
						try{
							final int[] hitCount=new int[targets.size()];
							final IntList touched=new IntList();
							final ShortlistScratch shortlistScratch=(shortlist>0 ?
								new ShortlistScratch(sidecar.dims, targets.size(), shortlist) : null);
							int qi;
							while((qi=next.getAndIncrement())<nQ){
								perQuery[qi]=(shortlist>0 ?
									searchOneQueryShortlist(queries.get(qi), targets, totalDbResidues, shortlistScratch) :
									searchOneQuery(queries.get(qi), targets, totalDbResidues, index, hitCount, touched));
							}
						}catch(Throwable e){
							errors[wi]=e;
						}
					}
				};
				workers[w].start();
			}
			for(Thread t : workers){
				try{t.join();}catch(InterruptedException e){throw new RuntimeException("Interrupted waiting for search worker", e);}
			}
			for(Throwable e : errors){
				if(e!=null){throw new RuntimeException("ProteinSearcher worker thread failed", e);}
			}
		}

		final ArrayList<ProteinHit> all=new ArrayList<ProteinHit>();
		for(List<ProteinHit> qhits : perQuery){all.addAll(qhits);}
		Collections.sort(all, TOTAL_ORDER);
		return all;
	}

	/**
	 * Computes one query's passing hits: seed, extend/score every candidate target, apply the
	 * significance filters, and cull to {@link #maxTargetSeqs}. Factored out of {@link
	 * #search(List, List)} so the identical logic runs whether called from the single-threaded
	 * path or a worker thread; {@code hitCount}/{@code touched} are caller-provided scratch
	 * (must be sized to {@code targets.size()} / start empty on first use) so each thread reuses
	 * its own arrays across queries without any cross-thread state.
	 *
	 * @param q The query.
	 * @param targets Full target list (read-only).
	 * @param totalDbResidues Sum of all target lengths, for the E-value search space.
	 * @param index Shared read-only target k-mer index.
	 * @param hitCount Per-target seed-hit tally scratch, reused across calls on the same thread.
	 * @param touched Indices touched during seeding, reused across calls on the same thread.
	 * @return This query's passing hits (not yet in the frozen total order).
	 */
	private List<ProteinHit> searchOneQuery(final ProteinSequence q, final List<ProteinSequence> targets,
			final long totalDbResidues, final LongLongListHashMap index,
			final int[] hitCount, final IntList touched){
		final double searchSpace=(double)q.length()*(double)totalDbResidues;

		//Seed: tally shared k-mers per target.
		for(int ti=0; ti<touched.size; ti++){hitCount[touched.get(ti)]=0;}
		touched.clear();
		final LongHashSet qkmers=kmerSet(q.enc);
		if(qkmers.isEmpty()){
			//Query too short for a k-mer: fall back to aligning against all targets. Every
			//touched index must be recorded so the NEXT call's reset step actually zeroes it --
			//otherwise hitCount[] stays latched at minSeedHits for every target forever,
			//permanently defeating the seed filter for the rest of this thread's queries.
			for(int i=0; i<targets.size(); i++){hitCount[i]=minSeedHits; touched.add(i);}
		}else{
			final long[] kmerArray=qkmers.toArray();
			for(long km : kmerArray){
				final LongList list=index.get(km);
				if(list!=null){
					for(int j=0; j<list.size; j++){
						final int idx=(int)list.get(j);
						if(hitCount[idx]==0){touched.add(idx);}
						hitCount[idx]++;
					}
				}
			}
		}

		//Extend + score every candidate target.
		final ArrayList<ProteinHit> qhits=new ArrayList<ProteinHit>();
		final int nTargets=targets.size();
		for(int i=0; i<nTargets; i++){
			if(qkmers.isEmpty() ? true : hitCount[i]>=minSeedHits){
				final ProteinSequence t=targets.get(i);
				final AAAlignment aln=AAAligner.align(q.enc, t.enc);
				if(aln==null){continue;}
				if(aln.rawScore<minRawScore){continue;}
				final double pid=aln.pident();
				if(pid<minPident){continue;}
				if(minCoverage>0){
					//Bidirectional (mmseqs cov-mode 0): BOTH sequences must be covered.
					final double qCov=aln.length/(double)q.length();
					final double tCov=aln.length/(double)t.length();
					if(qCov<minCoverage || tCov<minCoverage){continue;}
				}
				final double e=aln.evalue(searchSpace);
				if(e>evalueCutoff){continue;}
				qhits.add(new ProteinHit(q.id, t.id, aln, e, EVALUE_APPROXIMATE));
			}
		}

		//maxTargetSeqs: keep top-N distinct targets by best-HSP score (bitscore).
		if(maxTargetSeqs<qhits.size()){
			Collections.sort(qhits, BY_SCORE_DESC);
			while(qhits.size()>maxTargetSeqs){qhits.remove(qhits.size()-1);}
		}
		return qhits;
	}

	/** Per-thread reusable scratch for {@link #searchOneQueryShortlist} -- zero heap allocation
	 *  per query in steady state (matches {@code hitCount}/{@code touched}'s role for the
	 *  unshortlisted path). Sized once per thread from the sidecar's dims/family count.
	 *  Package-private (2026-09-03, plans/PER_FAMILY_THRESHOLDS_v1.md WP-A A1): {@link
	 *  FamilyCalibrationAssay} needs its own thread-local scratch to call {@link #scoreFamiliesF4}
	 *  without a second implementation of the scratch layout. */
	public static final class ShortlistScratch{
		ShortlistScratch(final int dims, final int nFamilies, final int shortlistN){
			counts=new double[dims]; centered=new double[dims]; qVec=new float[dims];
			compScores=new float[nFamilies]; f4=new float[nFamilies]; sk=new int[nFamilies]; tk=new double[nFamilies];
			topIdx=new int[Math.min(shortlistN, nFamilies)];
		}
		final double[] counts, centered; final float[] qVec, compScores, f4; final int[] sk, topIdx;
		/** The normalized k-mer term per family ({@code kmernorm=sub}); untouched under {@code none}. */
		final double[] tk;
		/** Qsize for the CURRENT query: the count of distinct query covering-set k-mer words, set by {@link
		 *  ProteinSearcher#scoreFamiliesF4} each call from {@code qWords.size()}. It is the k-mer query-density
		 *  denominator (D105 {@code s_k/Qsize}) read by {@link ProteinSearcher#assignFamily}; 0 for a query with no
		 *  covering words. Per-query mutable (like the other buffers' contents), captured once during the single F4
		 *  pass so the assigner never rebuilds the word set. Added 2026-09-14 (UMP45, Patch 2). */
		int qSize;
		/** Caller-owned distinct-word scratch shared by F4 and direct-family k-mer scoring. */
		final FamilyKmerScratch familyKmers=new FamilyKmerScratch();
		/** Reused D55/D60 Option-C metric holder for the shared {@link ProteinSearcher#d55Metrics}
		 *  computation (the calibration emitter and the final assigner both fill it) -- one per thread,
		 *  zero per-candidate allocation on the hot path. Untouched by the AAAligner (non-d55) path. */
		final D55Metrics d55m=new D55Metrics();
		/** Reused, query-local construction evidence; cleared by assignFamily before
		 *  every fixed-construction call and untouched by production calls. */
		public final ConstructionTrace constructionTrace=new ConstructionTrace();
		/** Metrics of the most recently aligned fixed-construction candidate. Filled
		 *  only by alignConstructionCandidate; kept here to avoid a per-alignment object. */
		int constructionIdx,constructionRawScore;
		float constructionR,constructionIdentity,constructionMutualCovQ,
			constructionMutualCovT,constructionMutualOverlap;
		/** Candidate, traceback, and acceptance counts from the latest assignment. */
		public int lastAlignedCandidates,lastTracebackAlignments,lastAcceptedCandidates;
	}

	/** Per-thread caller-owned scratch for direct-family k-mer scoring. */
	static final class FamilyKmerScratch{
		final LongHashSet set=new LongHashSet(256);
		final LongList words=new LongList(256);
	}

	/** Reusable result holder for one query against one family's covering set. */
	static final class FamilyKmerMetrics{
		int sK, qSize;
		/** {@code sK/qSize}; NaN when the query has no valid covering-set word. */
		double qDensity;
	}

	/**
	 * Scores one query directly against one family's covering set without computing composition or
	 * F4 for every family.  Query-word enumeration and the loaded covering index are the same ones
	 * used by {@link #scoreFamiliesF4}, so {@code sK}/{@code qSize} have live assignment semantics.
	 * The caller owns and reuses both scratch and result.
	 */
	static void scoreFamilyKmers(final ProteinSequence q, final FamilyShortlistSidecar sc,
			final int familyRank, final FamilyKmerScratch scratch, final FamilyKmerMetrics out){
		if(q==null || sc==null || scratch==null || out==null){
			throw new IllegalArgumentException("scoreFamilyKmers requires non-null query, sidecar, scratch, and result.");
		}
		if(familyRank<0 || familyRank>=sc.nFamilies){
			throw new IllegalArgumentException("scoreFamilyKmers family rank "+familyRank+
				" out of [0,"+sc.nFamilies+").");
		}
		ReducedAlphabetSeedAssay.fillKmerSet(q.enc,sc.coveringAlphabet,sc.coveringK,
			ReducedAlphabetSeedAssay.Boundary.NONE,scratch.set,scratch.words);
		int shared=0;
		for(int i=0; i<scratch.words.size; i++){
			final LongList hit=sc.coveringIndex.get(scratch.words.get(i));
			if(hit==null){continue;}
			for(int j=0; j<hit.size; j++){if(hit.get(j)==familyRank){shared++;}}
		}
		out.sK=shared;
		out.qSize=scratch.words.size;
		out.qDensity=(out.qSize<1 ? Double.NaN : out.sK/(double)out.qSize);
		assert(out.sK<=out.qSize) : "shared distinct k-mers exceed query distinct k-mers for rank "+
			familyRank+": sK="+out.sK+" qSize="+out.qSize+
			" (FamilyShortlistSidecar covering rows must contain each family rank at most once per word)";
	}

	/**
	 * Computes the stage 1-3 per-family shortlist scoring (composition cosine {@code s.compScores},
	 * covering-set k-mer count {@code s.sk}, and the fused F4 = z(s_k)+2*z(s_c)-d_L in {@code s.f4})
	 * for query {@code q} against every family in {@code sc}. This is the ONLY implementation of
	 * F4 -- both {@link #searchOneQueryShortlist} (production) and {@code FamilyCalibrationAssay}
	 * (calibration/evaluation, plans/PER_FAMILY_THRESHOLDS_v1.md WP-A A1) call this method, so a
	 * calibration run's stage 1-3 metrics are guaranteed to match production exactly -- extracted
	 * 2026-09-03 from what was previously inline in {@link #searchOneQueryShortlist}, with NO
	 * logic change (verified by the A1 byte-identical differential gate on {@code ProteinSearch}
	 * output before/after this extraction).
	 *
	 * @param q The query.
	 * @param sc The sidecar (family composition/covering-set/length data).
	 * @param s Caller-provided, thread-local scratch (reused across calls on the same thread);
	 *        {@code s.compScores}, {@code s.sk}, and {@code s.f4} are filled for every family.
	 */
	static void scoreFamiliesF4(final ProteinSequence q, final FamilyShortlistSidecar sc, final ShortlistScratch s){
		scoreFamiliesF4(q, sc, s, KMERNORM_NONE);
	}

	/**
	 * As {@link #scoreFamiliesF4(ProteinSequence, FamilyShortlistSidecar, ShortlistScratch)} with the
	 * k-mer term chosen by {@code kmerNorm} ({@link #kmerNorm}): {@code "none"} runs the original
	 * arithmetic unchanged; {@code "sub"} z-scores {@code s_k - rho_q*|C_f|} instead of {@code s_k}
	 * (the composition and length terms are untouched). The 3-arg overload is the {@code none} path,
	 * so {@code FamilyCalibrationAssay} and every pre-existing caller are byte-identical.
	 */
	static void scoreFamiliesF4(final ProteinSequence q, final FamilyShortlistSidecar sc, final ShortlistScratch s,
			final String kmerNorm){
		final int nFamilies=sc.nFamilies;
		if(!validKmerNorm(kmerNorm)){
			throw new RuntimeException("Unknown kmernorm: '"+kmerNorm+"' (expected "+KMERNORM_NONE+", "+KMERNORM_SUB+" or "+KMERNORM_CS300+").");
		}
		final boolean sub=!kmerNorm.equals(KMERNORM_NONE);//true for sub AND cs300: both substitute the k-mer term
		final boolean cs300=kmerNorm.equals(KMERNORM_CS300);

		//Composition: same recipe as the sidecar (aa20 2-mer, centered vs sc.background) -- MUST
		//use the sidecar's own recorded rawpseudo/variant/background, never recomputed, or a live
		//query's vector would not be comparable to the family vectors (see FamilyShortlistSidecar's
		//javadoc and FamilyShortlistSidecarBuilder's class javadoc for why).
		Arrays.fill(s.counts, 0);
		CompositionProfileAssay.addCounts(q.enc, sc.compositionAlphabet, sc.compositionK, s.counts);
		//toFreq() allocates its own return array (a shared CompositionProfileAssay utility used by
		//several callers with this signature; not worth widening just for this one caller) -- every
		//other per-query buffer below is thread-local reused scratch, never reallocated.
		final double[] freq=CompositionProfileAssay.toFreq(s.counts, sc.dims, sc.rawPseudocount);
		CompositionProfileAssay.applyVariant(freq, sc.background, sc.variant, s.centered);
		CompositionProfileAssay.normalize(s.centered);
		for(int d=0; d<sc.dims; d++){s.qVec[d]=(float)s.centered[d];}
		//Fixed-order SCALAR dot (mag-qc GATE-1 #8, Yoimiya Option A, 2026-09-11), NOT Vector.dotRows: the SIMD
		//reduction is not bit-reproducible across JIT states, and s.compScores (composition cosine s_c) is the ONLY
		//float-nondeterminism source in F4 -- s.sk is integer, topN is deterministic -- so the SIMD dot flipped the
		//top-50 shortlist RANK ORDER (not membership) across fresh JVMs, breaking hash-pinnability. Keep this scalar
		//operation in the canonical MAG-QC source rather than depending on an unsealed BBTools helper: the accepted
		//v5 load-gate evidence exposed that Vector.dotRowsScalar had never entered its declared clean BBTools commit.
		//The loop is left-to-right JLS float summation, decoupled from Shared.SIMD, so the exact SIMD aligner and its
		//in-run parity stay on. Scalar shortlist cost is negligible (measured).
		dotRowsScalar(s.qVec, sc.compositionMatrix, sc.dims, nFamilies, s.compScores);

		//Covering-set k-mer hits: distinct query k-mers found in each family's covering set.
		//The set/list are caller-owned per-thread scratch.  Direct-family calibration calls the
		//same ReducedAlphabetSeedAssay.fillKmerSet seam, preventing an alternate packing/reset path.
		Arrays.fill(s.sk, 0);
		final FamilyKmerScratch kscratch=s.familyKmers;
		ReducedAlphabetSeedAssay.fillKmerSet(q.enc,sc.coveringAlphabet,sc.coveringK,
			ReducedAlphabetSeedAssay.Boundary.NONE,kscratch.set,kscratch.words);
		//Qsize = distinct query covering-set words (D105 s_k/Qsize denominator), captured once here from the same set
		//scoreFamiliesF4 already builds, so assignFamily need not rebuild it (Patch 2, UMP45 2026-09-14). ADDITIVE:
		//writes a NEW scratch field only -- f4/sk/compScores are byte-identical, so the A1 differential still holds
		//(please re-run it to confirm). The list contains exactly the words newly inserted into the set.
		s.qSize=kscratch.words.size;
		for(int i=0; i<kscratch.words.size; i++){
			final long word=kscratch.words.get(i);
			final LongList hit=sc.coveringIndex.get(word);
			if(hit==null){continue;}
			for(int j=0; j<hit.size; j++){s.sk[(int)hit.get(j)]++;}
		}

		//z-scores of s_k/s_c across all families (same formula as HybridShortlistAssay's F1/F4).
		double sumSk=0, sumSkSq=0, sumSc=0, sumScSq=0;
		for(int r=0; r<nFamilies; r++){
			sumSk+=s.sk[r]; sumSkSq+=(double)s.sk[r]*s.sk[r];
			sumSc+=s.compScores[r]; sumScSq+=(double)s.compScores[r]*s.compScores[r];
		}
		if(sub){
			//kmernorm=sub: t_f = s_k - rho_q*|C_f|, rho_q = (sum over all families of s_k)/(sum of |C_f|);
			//kmernorm=cs300: t_f = s_k/(|C_f|+300) (HybridShortlistAssay methods 20-25 / f4n_sub, f4n_cs300;
			//results/shortlist_kmer_normalization_v1.md). The k-mer moments are then those of t, not s_k;
			//s_c and d_L are untouched.
			//Zero-background case (all covering sets empty, coveringSetSizeSum==0 -- the loader has always
			//allowed it): no k-mer can hit anything, so sumSk==0 too; rho_q is defined as 0 and every t_f is 0
			//(sub) or 0/(0+300)=0 (cs300), i.e. the k-mer term is a constant and F4 reduces to composition+length,
			//exactly what kmernorm=none computes when s_k==0 everywhere. Not a crash: legacy behavior preserved.
			assert(sc.coveringSetSizeSum>0 || sumSk==0) : "no covering k-mers loaded but s_k hits were counted: sumSk="+sumSk+
				" (scoreFamiliesF4 counts hits only from sc.coveringIndex, which is built from the same rows as coveringSetSize)";
			final double rhoQ=(sc.coveringSetSizeSum>0 ? sumSk/sc.coveringSetSizeSum : 0);
			assert(rhoQ>=0 && rhoQ<=1) : "rho_q="+rhoQ+" must be a rate in [0,1]: a query's shared k-mers with a family cannot exceed "+
				"that family's covering-set size (each set k-mer counts at most once per query); sumSk="+sumSk+" sizes="+sc.coveringSetSizeSum;
			sumSk=0; sumSkSq=0;
			for(int r=0; r<nFamilies; r++){
				final double t=(cs300 ? s.sk[r]/(sc.coveringSetSize[r]+300.0) : s.sk[r]-rhoQ*sc.coveringSetSize[r]);
				s.tk[r]=t; sumSk+=t; sumSkSq+=t*t;
			}
		}
		final double meanSk=sumSk/nFamilies, meanSc=sumSc/nFamilies;
		final double stdSk=Math.sqrt(Math.max(0, sumSkSq/nFamilies-meanSk*meanSk));
		final double stdSc=Math.sqrt(Math.max(0, sumScSq/nFamilies-meanSc*meanSc));

		final int qLen=q.length();
		final double lnQLen=Math.log(qLen);
		for(int r=0; r<nFamilies; r++){
			final double dL=Math.abs(lnQLen-Math.log(sc.medLen[r]))/Math.max(sc.sdLnLen[r], 0.05);
			final double kTerm=(sub ? s.tk[r] : s.sk[r]);//none: int->double promotion, exactly as before
			final double zk=(stdSk>1e-12 ? (kTerm-meanSk)/stdSk : 0);
			final double zc=(stdSc>1e-12 ? (s.compScores[r]-meanSc)/stdSc : 0);
			s.f4[r]=(float)(zk+2*zc-dL);
		}
	}

	/**
	 * Deterministic row-wise dot products for shortlist composition scores.
	 * Each matrix row begins at {@code row*stride}; only {@code q.length}
	 * values participate, accumulated in strict ascending-index order.
	 */
	static void dotRowsScalar(final float[] q, final float[] matrix, final int stride,
			final int nRows, final float[] out){
		if(q==null || matrix==null || out==null || stride<q.length || stride<1 ||
				nRows<0 || out.length<nRows || (nRows>0 && (long)nRows*stride>matrix.length)){
			throw new IllegalArgumentException("Invalid scalar row-dot arguments: q="+
				(q==null ? -1 : q.length)+" matrix="+(matrix==null ? -1 : matrix.length)+
				" stride="+stride+" rows="+nRows+" out="+(out==null ? -1 : out.length));
		}
		for(int row=0; row<nRows; row++){
			final int base=row*stride;
			float sum=0;
			for(int d=0; d<q.length; d++){sum+=q[d]*matrix[base+d];}
			out[row]=sum;
		}
	}

	/** {@code selfScore(q) = sum_i Blosum62.score(q_i, q_i)} -- the score ratio's denominator
	 *  (plans/PER_FAMILY_THRESHOLDS_v1.md sec 1, metric 4b; same definition Ady used on
	 *  ady/two-pass-scoring dd7879c, recomputed here rather than imported per the design).
	 *  Package-private (2026-09-04, WP-A A3): promoted from a private copy in {@link
	 *  FamilyCalibrationAssay} so there is exactly one implementation, matching A1's own "one
	 *  implementation of F4" principle for shared scoring math. */
	static int selfScore(final ProteinSequence q){
		int sum=0;
		for(final byte b : q.enc){sum+=Blosum62.score(b, b);}
		return sum;
	}

	/**
	 * D61 zeroed self-score of an encoded span: {@code Σ max(0, Blosum62.score(x,x))} over the inclusive
	 * index range [from,to]. X's BLOSUM62 diagonal is -1 (the only nonpositive diagonal), so it
	 * contributes 0 -- the R-denominator convention (D58/D61), NOT the plain {@link #selfScore}.
	 * Promoted 2026-09-14 (UMP45) from a private copy in {@link FamilyCalibrationAssay} so there is
	 * exactly ONE implementation shared by the calibration emitter and the final assigner (same
	 * principle as {@link #selfScore}/{@link #scoreFamiliesF4}).
	 */
	static int zeroedSelfScore(final byte[] enc, final int from, final int to){
		assert(from>=0 && to<enc.length && from<=to)
			: "bad span ["+from+","+to+"] len="+enc.length+" (D61 zeroedSelfScore invariant, promoted from FamilyCalibrationAssay)";
		int sum=0;
		for(int k=from; k<=to; k++){
			final int d=Blosum62.score(enc[k], enc[k]);
			if(d>0){sum+=d;}//X (d=-1) contributes 0
		}
		return sum;
	}

	/**
	 * Reusable holder for the D55/D60 Option-C metrics computed by {@link #d55Metrics} -- the SINGLE
	 * implementation shared by {@link FamilyCalibrationAssay}'s row emitter and the final assigner, so a
	 * calibration run's metrics are guaranteed identical to the assigner's decision inputs. Caller reuses
	 * one instance per thread (see {@link ShortlistScratch#d55m}) -- zero per-candidate allocation.
	 */
	static final class D55Metrics{
		/** Raw BLOSUM alignment score A (the D107 score-filter axis). */
		int A;
		/** D61 zeroed self-scores of the query span and the reference span (R's denominator arms). */
		int Sq, St;
		/** Symmetric normalized score {@code R = 2A/(Sq+St)}, unclamped (0 when Sq+St==0). D54 arm;
		 *  the D111 score-all RANKING metric (unweighted, Yoimiya 2026-09-14). */
		double R;
		/** Span coverages and their min (D55 span overlap). Under q-global geometry covq==1. */
		double covq, covt, overlap;
		/** Non-gap aligned-residue coverages. Unlike {@link #overlap}, these do not
		 *  count a residue aligned to a gap. Construction D119 uses their minimum. */
		double mutualCovQ, mutualCovT, mutualOverlap;
		/** Overlap-weighted arm {@code R*overlap^0.25} -- HISTORICAL (D107 dropped the 0.25 exponent
		 *  from acceptance); retained ONLY to keep calibration output byte-identical. Not a new-decision input. */
		double Rov;
	}

	/**
	 * Computes the D55/D60 Option-C metrics for one (query, family) D55 alignment into {@code out} -- the
	 * ONE implementation used by both {@link FamilyCalibrationAssay#emitRowD55} (calibration rows) and the
	 * final assigner (acceptance/ranking), so the numbers cannot diverge (the "one implementation of F4"
	 * principle, plans/PER_FAMILY_THRESHOLDS_v1.md A1, extended to the D55 metrics 2026-09-14 by UMP45).
	 * Reproduces exactly what {@code emitRowD55} computed inline; {@link D55MetricsSharedTest} pins byte
	 * parity against an independent oracle. {@code aln} is a D55 alignment ({@link GlocalAminoLinear#align},
	 * never null for nonempty input); q-global geometry (full query span) is asserted so a geometry drift
	 * crashes loud here rather than silently changing overlap.
	 * @param q Query. @param t The family consensus/target. @param aln The D55 alignment. @param out Reused holder.
	 */
	static void d55Metrics(final ProteinSequence q, final ProteinSequence t, final AAAlignment aln, final D55Metrics out){
		final int Lq=q.length(), Lt=t.length();
		assert(aln.qStart==0 && aln.qStop==Lq-1)
			: "D55 q-global expects full query span; qStart="+aln.qStart+" qStop="+aln.qStop+" Lq="+Lq+" q='"+q.id+"'";
		final int A=aln.rawScore;
		final int Sq=zeroedSelfScore(q.enc, 0, Lq-1);//query span = full query (q-global)
		final int St=zeroedSelfScore(t.enc, aln.tStart, aln.tStop);//reference span [tStart,tStop]
		final int denom=Sq+St;
		final double R=(denom>0) ? (2.0*A/denom) : 0.0;
		final double covq=(aln.qStop-aln.qStart+1)/(double)Lq;
		final double covt=(aln.tStop-aln.tStart+1)/(double)Lt;
		final double overlap=Math.min(covq, covt);
		final int paired=aln.identities+aln.mismatches;
		final double mutualCovQ=paired/(double)Lq, mutualCovT=paired/(double)Lt;
		final double mutualOverlap=Math.min(mutualCovQ,mutualCovT);
		final double Rov=R*Math.pow(overlap, 0.25);
		out.A=A; out.Sq=Sq; out.St=St; out.R=R; out.covq=covq; out.covt=covt; out.overlap=overlap;
		out.mutualCovQ=mutualCovQ; out.mutualCovT=mutualCovT; out.mutualOverlap=mutualOverlap; out.Rov=Rov;
	}

	/** Fixed, membership-independent acceptance profile used only while rebuilding
	 *  consensuses. The profile text is parsed to float32 once and its canonical
	 *  float strings are hashed, so every inclusive comparison and every receipt
	 *  names the exact values used. */
	public static final class FixedConstructionProfile{
		public static final int PROFILE_COUNT=3;
		public final float minLengthRatio,maxLengthRatio,minR,minMutualOverlap;
		public final float[] minIdentity;
		public final String[] minIdentityText;
		public final String minLengthRatioText,maxLengthRatioText,minRText,minMutualOverlapText,sha80;

		public FixedConstructionProfile(final String minLengthRatio_, final String maxLengthRatio_,
				final String minR_, final String minMutualOverlap_, final String identityCsv){
			minLengthRatio=parseFiniteFloat32(minLengthRatio_,"min_length_ratio",0,Float.POSITIVE_INFINITY);
			maxLengthRatio=parseFiniteFloat32(maxLengthRatio_,"max_length_ratio",minLengthRatio,Float.POSITIVE_INFINITY);
			minR=parseFiniteFloat32(minR_,"min_r",-Float.MAX_VALUE,Float.MAX_VALUE);
			minMutualOverlap=parseFiniteFloat32(minMutualOverlap_,"min_mutual_overlap",0,1);
			final String[] split=identityCsv.split(",",-1);
			if(split.length!=PROFILE_COUNT){
				throw new IllegalArgumentException("identity_profiles must contain exactly "+
					PROFILE_COUNT+" comma-separated float32 values");
			}
			minIdentity=new float[PROFILE_COUNT]; minIdentityText=new String[PROFILE_COUNT];
			for(int i=0; i<PROFILE_COUNT; i++){
				minIdentity[i]=parseFiniteFloat32(split[i],"identity_profiles["+i+"]",0,100);
				if(i>0 && minIdentity[i]<=minIdentity[i-1]){throw new IllegalArgumentException("identity_profiles must be strictly increasing");}
				minIdentityText[i]=Float.toString(minIdentity[i]);
			}
			minLengthRatioText=Float.toString(minLengthRatio); maxLengthRatioText=Float.toString(maxLengthRatio);
			minRText=Float.toString(minR); minMutualOverlapText=Float.toString(minMutualOverlap);
			final String canonical="min_length_ratio\t"+minLengthRatioText+'\n'+
				"max_length_ratio\t"+maxLengthRatioText+'\n'+"min_r\t"+minRText+'\n'+
				"min_mutual_overlap\t"+minMutualOverlapText+'\n'+"identity_profiles\t"+
				minIdentityText[0]+','+minIdentityText[1]+','+minIdentityText[2]+'\n'+
				"raw_score_gate\tPASS_ALL_INT\n"+
				"kmer_count_gate\tPASS_ALL_INT\n"+
				"kmer_density_gate\tPASS_ALL_DOUBLE\n"+
				"hbm_path_gate\tPASS_ALL_FLOAT\n";
			sha80=DigestSuffix.bytes(canonical.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		}

		private static float parseFiniteFloat32(final String text, final String label, final float min, final float max){
			final float value;
			try{value=Float.parseFloat(text);}catch(final NumberFormatException e){throw new IllegalArgumentException(label+" is not a float32: "+text,e);}
			if(!Float.isFinite(value) || value<min || value>max){throw new IllegalArgumentException(label+" out of range ["+min+","+max+"]: "+text);}
			final String canonical=Float.toString(value);
			if(Float.floatToIntBits(Float.parseFloat(canonical))!=Float.floatToIntBits(value)){
				throw new IllegalArgumentException(label+" failed float32 round-trip: "+text+" -> "+canonical);
			}
			return value;
		}
	}

	/** Reused construction-only evidence from the most recent {@link #assignFamily}
	 *  call on a scratch object. Slots 0/1 are the best and runner-up aligned
	 *  candidates; slots 2..4 are the exact winners at the three identity floors. */
	public static final class ConstructionTrace{
		public static final int BEST=0,RUNNER_UP=1,PROFILE_OFFSET=2,SLOT_COUNT=5;
		public final int[] familyIdx=new int[SLOT_COUNT],rawScore=new int[SLOT_COUNT],coreLength=new int[SLOT_COUNT];
		public final int[] profileSurvivors=new int[FixedConstructionProfile.PROFILE_COUNT];
		public final float[] R=new float[SLOT_COUNT],identity=new float[SLOT_COUNT],mutualCovQ=new float[SLOT_COUNT],
			mutualCovT=new float[SLOT_COUNT],mutualOverlap=new float[SLOT_COUNT];
		public int lengthSurvivors,alignedCandidates,rSurvivors,mutualSurvivors;
		ConstructionTrace(){clear();}
		void clear(){
			Arrays.fill(familyIdx,-1); Arrays.fill(R,Float.NEGATIVE_INFINITY);
			Arrays.fill(rawScore,0); Arrays.fill(coreLength,0); Arrays.fill(identity,Float.NaN);
			Arrays.fill(mutualCovQ,Float.NaN); Arrays.fill(mutualCovT,Float.NaN); Arrays.fill(mutualOverlap,Float.NaN);
			Arrays.fill(profileSurvivors,0);
			lengthSurvivors=alignedCandidates=rSurvivors=mutualSurvivors=0;
		}
		void set(final int slot, final int idx, final int A, final int l50, final float r, final float id,
				final float covQ, final float covT, final float mutual){
			familyIdx[slot]=idx; rawScore[slot]=A; coreLength[slot]=l50; R[slot]=r; identity[slot]=id;
			mutualCovQ[slot]=covQ; mutualCovT[slot]=covT; mutualOverlap[slot]=mutual;
		}
		void copy(final int from, final int to){
			familyIdx[to]=familyIdx[from]; rawScore[to]=rawScore[from]; coreLength[to]=coreLength[from];
			R[to]=R[from]; identity[to]=identity[from]; mutualCovQ[to]=mutualCovQ[from];
			mutualCovT[to]=mutualCovT[from]; mutualOverlap[to]=mutualOverlap[from];
		}
	}

	/** Result of exhaustive fixed-construction scoring for shortlist-recall measurement.
	 *  The exhaustive winner and runner-up are the two highest-R families that pass
	 *  the strict construction profile across every length-eligible family. Shortlist
	 *  ranks are one-based; zero means absent. */
	public static final class ConstructionRecallResult{
		public final AssignReason reason;
		public final int shortlistSize,lengthSurvivors,alignedCandidates,acceptedCandidates;
		public final int winnerIdx,runnerUpIdx,winnerShortlistRank,runnerUpShortlistRank;
		public final int winnerRawScore,runnerUpRawScore;
		public final float winnerR,runnerUpR,winnerIdentity,runnerUpIdentity;
		public final float winnerMutualOverlap,runnerUpMutualOverlap;
		private ConstructionRecallResult(final AssignReason reason_, final int shortlistSize_,
				final int lengthSurvivors_, final int alignedCandidates_, final int acceptedCandidates_,
				final int winnerIdx_, final int runnerUpIdx_, final int winnerShortlistRank_,
				final int runnerUpShortlistRank_, final int winnerRawScore_, final int runnerUpRawScore_,
				final float winnerR_, final float runnerUpR_, final float winnerIdentity_,
				final float runnerUpIdentity_, final float winnerMutualOverlap_,
				final float runnerUpMutualOverlap_){
			reason=reason_; shortlistSize=shortlistSize_; lengthSurvivors=lengthSurvivors_;
			alignedCandidates=alignedCandidates_; acceptedCandidates=acceptedCandidates_;
			winnerIdx=winnerIdx_; runnerUpIdx=runnerUpIdx_; winnerShortlistRank=winnerShortlistRank_;
			runnerUpShortlistRank=runnerUpShortlistRank_; winnerRawScore=winnerRawScore_;
			runnerUpRawScore=runnerUpRawScore_; winnerR=winnerR_; runnerUpR=runnerUpR_;
			winnerIdentity=winnerIdentity_; runnerUpIdentity=runnerUpIdentity_;
			winnerMutualOverlap=winnerMutualOverlap_; runnerUpMutualOverlap=runnerUpMutualOverlap_;
		}
		public boolean hasWinner(){return winnerIdx>=0;}
		public boolean hasRunnerUp(){return runnerUpIdx>=0;}
		public boolean winnerWithin(final int n){return winnerShortlistRank>0 && winnerShortlistRank<=n;}
		public boolean runnerUpWithin(final int n){return runnerUpShortlistRank>0 && runnerUpShortlistRank<=n;}
	}

	/** Endpoint policy for {@link #assignFamily} (P3 selects the winner between them, PATH_TO_FINISHED_PRODUCT_v1 §D). */
	public enum AssignPolicy{
		/** First clearing candidate in F4 rank order (the highest-F4 acceptance-passer). */
		FIRST_CLEAR,
		/** Continue X candidates beyond the last strict accepted R advance, then choose greatest R. */
		BOUNDED_LOOKAHEAD,
		/** All acceptance-passers, choose the greatest normalized score R (Yoimiya 2026-09-14, unweighted; tie: R desc, then rep_id asc). */
		SCORE_ALL
	}

	/** Stable receipt/header value for the bounded-lookahead stopping rule. */
	public static final String BOUNDED_LOOKAHEAD_SEMANTICS=
		"after_last_strict_accepted_normalized_d55_R_advance_v1";

	/** Validates the shared policy/lookahead contract used by every assignment entry point. */
	static void validateLookahead(final AssignPolicy policy, final int lookahead){
		if(policy==null){throw new IllegalArgumentException("Assignment policy must not be null.");}
		if(policy==AssignPolicy.BOUNDED_LOOKAHEAD){
			if(lookahead<0){
				throw new IllegalArgumentException(
					"Bounded lookahead must be nonnegative: "+lookahead);
			}
		}else if(lookahead!=0){
			throw new IllegalArgumentException("lookahead="+lookahead+
				" requires BOUNDED_LOOKAHEAD policy.");
		}
	}

	/**
	 * Updates a bounded-lookahead stop boundary after one accepted candidate.
	 * A candidate resets the boundary only when its normalized D55 R is strictly
	 * greater than the prior running best. Equal-R rep-id tie improvements can
	 * change the winner, but do not buy more alignments.
	 *
	 * @param priorStopExclusive Existing exclusive stop boundary.
	 * @param zeroBasedRank Candidate's zero-based shortlist rank.
	 * @param cap Number of shortlisted candidates available.
	 * @param lookahead Number of candidates to evaluate after a strict R advance.
	 * @param candidateR Accepted candidate's normalized D55 R.
	 * @param priorBestR Best accepted R before this candidate.
	 * @return The unchanged boundary, or {@code min(cap,rank+1+lookahead)} after a strict advance.
	 */
	static int boundedStopExclusive(final int priorStopExclusive,
			final int zeroBasedRank, final int cap, final int lookahead,
			final double candidateR, final double priorBestR){
		if(cap<1 || zeroBasedRank<0 || zeroBasedRank>=cap){
			throw new IllegalArgumentException("Invalid bounded-lookahead rank/cap: rank="+
				zeroBasedRank+" cap="+cap);
		}
		if(priorStopExclusive<1 || priorStopExclusive>cap){
			throw new IllegalArgumentException("Invalid bounded-lookahead stop boundary "+
				priorStopExclusive+" for cap="+cap);
		}
		if(lookahead<0){throw new IllegalArgumentException(
			"Bounded lookahead must be nonnegative: "+lookahead);}
		if(candidateR<=priorBestR){return priorStopExclusive;}
		return (int)Math.min((long)cap,(long)zeroBasedRank+1L+(long)lookahead);
	}

	/** Why canonical assignment returned what it did (an assigned tracked family, a bait intercept,
	 *  or one of the typed rejections). {@code INVALID_SEQUENCE} is produced by the canonical input
	 *  adapters before {@link #assignFamily}, because an invalid raw record cannot become a
	 *  {@link ProteinSequence}; it still belongs to this shared result vocabulary.
	 *  {@code NO_VALID_KMER} added D117 (Brian, 2026-09-15): for alphabet=20/k=5/boundary=NONE, a called
	 *  protein with no valid covering-set 5-mer (empty-Q, {@code Qsize==0}) is discarded from family
	 *  assignment with this typed reason, contributes NO family count, and processing continues -- see
	 *  {@link #assignFamily}'s empty-Q guard. */
	public enum AssignReason{ ASSIGNED, BAIT_INTERCEPTED, MASKED_INTERCEPTED,
		NO_PRESELECT_SURVIVOR, NO_ACCEPTANCE_PASS, NO_VALID_KMER, INVALID_SEQUENCE }

	/** The canonical assigner's zero-or-one result (§E: "returns zero or one family per protein"). */
	public static final class FamilyAssignment{
		/** Winning tracked, masked, or bait family rank, or -1 when no candidate won. */
		public final int familyIdx;
		/** Winning tracked, masked, or bait family rep_id, or null when no candidate won. */
		public final String repId;
		public final AssignReason reason;
		/** Winning candidate metrics for assigned, masked-intercepted, and bait-intercepted outcomes. */
		public final int rawScore, kmerCount; public final double identity, overlap, R, kmerQdensity;
		/** HBM scores for an accepted winner; zero sentinels for typed rejections. */
		public final float hbmRaw,hbmPathRelative;
		private FamilyAssignment(final int familyIdx, final String repId, final AssignReason reason,
				final int rawScore, final double identity, final double overlap, final double R,
				final int kmerCount, final double kmerQdensity, final float hbmRaw,
				final float hbmPathRelative){
			this.familyIdx=familyIdx; this.repId=repId; this.reason=reason; this.rawScore=rawScore;
			this.identity=identity; this.overlap=overlap; this.R=R; this.kmerCount=kmerCount; this.kmerQdensity=kmerQdensity;
			this.hbmRaw=hbmRaw; this.hbmPathRelative=hbmPathRelative;
		}
		static FamilyAssignment rejected(final AssignReason reason){
			return new FamilyAssignment(-1, null, reason, 0, 0, 0, 0, 0, 0, 0, 0);
		}
		static FamilyAssignment assigned(final int idx, final String repId, final int rawScore, final double identity,
				final double overlap, final double R, final int kmerCount, final double kmerQdensity,
				final float hbmRaw, final float hbmPathRelative){
			return new FamilyAssignment(idx, repId, AssignReason.ASSIGNED, rawScore, identity,
				overlap, R, kmerCount, kmerQdensity,hbmRaw,hbmPathRelative);
		}
		static FamilyAssignment assigned(final int idx, final String repId,
				final int rawScore, final double identity, final double overlap,
				final double R, final int kmerCount, final double kmerQdensity,
				final float hbmRaw){
			return assigned(idx,repId,rawScore,identity,overlap,R,kmerCount,
				kmerQdensity,hbmRaw,0);
		}
		static FamilyAssignment baitIntercepted(final int idx, final String repId, final int rawScore, final double identity,
				final double overlap, final double R, final int kmerCount, final double kmerQdensity,
				final float hbmRaw, final float hbmPathRelative){
			return new FamilyAssignment(idx, repId, AssignReason.BAIT_INTERCEPTED,
				rawScore, identity, overlap, R, kmerCount, kmerQdensity,hbmRaw,
				hbmPathRelative);
		}
		static FamilyAssignment baitIntercepted(final int idx, final String repId,
				final int rawScore, final double identity, final double overlap,
				final double R, final int kmerCount, final double kmerQdensity,
				final float hbmRaw){
			return baitIntercepted(idx,repId,rawScore,identity,overlap,R,kmerCount,
				kmerQdensity,hbmRaw,0);
		}
		static FamilyAssignment maskedIntercepted(final int idx, final String repId,
				final int rawScore, final double identity, final double overlap,
				final double R, final int kmerCount, final double kmerQdensity,
				final float hbmRaw, final float hbmPathRelative){
			return new FamilyAssignment(idx, repId, AssignReason.MASKED_INTERCEPTED,
				rawScore, identity, overlap, R, kmerCount, kmerQdensity,hbmRaw,
				hbmPathRelative);
		}
		/** True only for a countable tracked/network assignment. */
		public boolean isAssigned(){return reason==AssignReason.ASSIGNED;}
		public boolean isBaitIntercepted(){return reason==AssignReason.BAIT_INTERCEPTED;}
		public boolean isMaskedIntercepted(){return reason==AssignReason.MASKED_INTERCEPTED;}
	}

	/**
	 * A trusted binding of Eru's {@link FamilyAcceptanceConfig} to the sidecar + target consensus list the assigner runs
	 * against. Built ONLY through {@link #load} (private constructor), which loads all three from the live files, hash-binds
	 * them, and OWNS them privately -- so a caller cannot inject a mismatched/mutated resource or reorder the targets, and
	 * cannot mutate anything after validation (there is no external alias). {@link #assignFamily} consumes the binding and
	 * does no per-query roster/resource validation of its own. (Yoimiya's consolidated review, 2026-09-14: a validate-only
	 * check over caller-supplied objects could not prove residue bytes, could skip the sidecar-file binding, and left a
	 * mutate-after-validate window -- all closed by loading and owning the objects here.)
	 */
	public static final class AssignmentBinding{
		private final FamilyShortlistSidecar sidecar;
		private final List<ProteinSequence> targets;
		private final FamilyAcceptanceConfig cfg;
		private final HbmBundleLoader.Loaded hbm;
		private final FixedConstructionProfile construction;
		private final FamilyCoreCoordinates core;
		private AssignmentBinding(final FamilyShortlistSidecar sidecar, final List<ProteinSequence> targets,
				final FamilyAcceptanceConfig cfg, final HbmBundleLoader.Loaded hbm,
				final FixedConstructionProfile construction_, final FamilyCoreCoordinates core_){
			this.sidecar=sidecar; this.targets=targets; this.cfg=cfg; this.hbm=hbm; construction=construction_; core=core_;
			if((cfg==null)!=(construction_!=null)){throw new IllegalArgumentException("AssignmentBinding must contain exactly one acceptance profile");}
			if(construction_!=null && core_==null){throw new IllegalArgumentException("Construction binding requires core coordinates");}
			if(cfg!=null && cfg.isSchema7()!=(core_!=null)){
				throw new IllegalArgumentException("Only schema-7 production bindings require core coordinates");
			}
		}
		/** Allocates one worker's scratch without exposing the bound mutable resources. */
		public ShortlistScratch newScratch(){
			return new ShortlistScratch(sidecar.dims, sidecar.nFamilies, Math.min(50, sidecar.nFamilies));
		}

		/** Immutable content identity of the loaded threshold profile: the artifact's own data-row
		 *  SHA-256, recomputed and verified against the artifact header at load ({@link
		 *  FamilyAcceptanceConfig#thresholdsSha256}). For report provenance -- returns a String
		 *  (immutable), never the private cfg arrays. */
		public String profileThresholdsSha80(){return construction==null ? cfg.thresholdsSha256 : construction.sha80;}
		/** Complete sealed schema-7 artifact identity, or {@code NA} for legacy/construction profiles. */
		public String profileArtifactSha80(){return construction==null ? cfg.profileArtifactSha80 : "NA";}
		public String hbmBundleSha80(){return construction==null ? cfg.hbmBundleSha256 : "NA";}
		public String hbmSemanticProvenanceSha80(){return construction==null ? cfg.hbmSemanticProvenanceSha256 : "NA";}
		public boolean isConstruction(){return construction!=null;}
		public int coveringK(){return sidecar.coveringK;}
		public String repId(final int rank){return sidecar.repIds[rank];}
		public int familyCount(){return sidecar.nFamilies;}
		public String coreCoordinatesSha80(){return core==null ? "NA" : core.artifactSha80;}
		/** Legacy source-compatibility alias; the returned text is always a 20-character suffix. */
		@Deprecated public String profileThresholdsSha256(){return profileThresholdsSha80();}
		/** Legacy source-compatibility alias; the returned text is always a 20-character suffix. */
		@Deprecated public String hbmBundleSha256(){return hbmBundleSha80();}
		/** Legacy source-compatibility alias; the returned text is always a 20-character suffix. */
		@Deprecated public String hbmSemanticProvenanceSha256(){return hbmSemanticProvenanceSha80();}
		/**
		 * The ONLY way to build a binding. Loads the config, sidecar, and target consensus list from the live files and
		 * owns them privately. Every resource is hash-bound: {@link FamilyAcceptanceConfig#load} re-hashes the
		 * roster/consensus/covering-set/sidecar FILES against the artifact (and checks the caller-supplied,
		 * reviewed-provenance aligner/BLOSUM/kmer-definition), so tampered target residues (a changed consensus) or a
		 * swapped sidecar fail HERE; the sidecar OBJECT and the targets are then loaded from those SAME verified files.
		 *
		 * @param artifactPath The strict schema-4 {@code d55_d107_hbm_family_thresholds} artifact.
		 * @param familyListPath Live roster {@code #rank\trep_id\tocc_total} (defines the authoritative rank order).
		 * @param consensusPath Live consensus FASTA (targets are read from THIS file, so their residues cannot be a
		 *        caller-injected variant -- a different one changes the hash and fails config load).
		 * @param roleManifestPath Live rank/rep-ID/role/network-feature manifest, hash-bound by the threshold profile.
		 * @param coveringSetsPath Live covering-set k-mer file.
		 * @param sidecarPath Live {@link FamilyShortlistSidecar} TSV (the sidecar OBJECT is loaded from this exact file).
		 * @param trustedAlignerSha256 Reviewed provenance hash of the live D55 aligner source (NOT the artifact's own claim).
		 * @param trustedBlosum62Sha256 Reviewed provenance hash of the live BLOSUM62 source.
		 * @param trustedKmerDefinition Reviewed canonical kmer-definition contract value.
		 * @param hbmBundlePath Live MQHB bundle hash-bound by the profile.
		 * @param hbmProvenancePath Live 16-row semantic-provenance file hash-bound by the profile.
		 * @return A validated, self-owned binding.
		 */
		public static AssignmentBinding load(final String artifactPath, final String familyListPath,
				final String consensusPath, final String roleManifestPath, final String coveringSetsPath, final String sidecarPath,
				final String trustedAlignerSha256, final String trustedBlosum62Sha256, final String trustedKmerDefinition,
				final String hbmBundlePath, final String hbmProvenancePath){
			return load(artifactPath,familyListPath,consensusPath,roleManifestPath,
				coveringSetsPath,sidecarPath,trustedAlignerSha256,
				trustedBlosum62Sha256,trustedKmerDefinition,hbmBundlePath,
				hbmProvenancePath,1);
		}

		public static AssignmentBinding load(final String artifactPath, final String familyListPath,
				final String consensusPath, final String roleManifestPath, final String coveringSetsPath, final String sidecarPath,
				final String trustedAlignerSha256, final String trustedBlosum62Sha256, final String trustedKmerDefinition,
				final String hbmBundlePath, final String hbmProvenancePath,
				final int hbmMinCount){
			//1. Config: re-hashes roster/consensus/coveringsets/SIDECAR files vs the artifact + verifies the trusted
			//aligner/BLOSUM/kmer-definition provenance (Eru's loader). Tampered consensus residues or a swapped sidecar fail here.
			final FamilyAcceptanceConfig cfg=FamilyAcceptanceConfig.load(artifactPath, familyListPath, consensusPath, roleManifestPath,
				"D55", trustedAlignerSha256, trustedBlosum62Sha256, coveringSetsPath, sidecarPath, trustedKmerDefinition,
				hbmBundlePath,hbmProvenancePath);
			//2. Sidecar OBJECT from the SAME sidecar file the config just hash-verified (expectedAligner=null: the sidecar's
			//own v3 per-family thresholds are NEVER used by the assigner -- only cfg.* governs acceptance -- so no aligner binding).
			final FamilyShortlistSidecar sidecar=FamilyShortlistSidecar.load(sidecarPath, consensusPath, coveringSetsPath);
			//3. Targets from the SAME consensus file the config just hash-verified: residues cannot be a caller-injected/mutated
			//variant (that would have changed the consensus hash and failed step 1).
			final List<ProteinSequence> loaded=ProteinSearch.readFasta(consensusPath);
			//4. In-memory cross-checks (defense in depth; all three came from hash-verified files).
			final int nF=cfg.nFamilies;
			if(sidecar.nFamilies!=nF || loaded.size()!=nF){
				throw new RuntimeException("AssignmentBinding family-count mismatch: cfg="+nF+" sidecar="+sidecar.nFamilies+" targets="+loaded.size());
			}
			for(int i=0; i<nF; i++){
				final String tid=loaded.get(i).id, sid=sidecar.repIds[i], cid=cfg.repId(i);
				if(!tid.equals(sid) || !sid.equals(cid)){
					throw new RuntimeException("AssignmentBinding roster mismatch at index "+i+": target='"+tid+"' sidecar='"+sid+"' cfg='"+cid+"'.");
				}
			}
			if(!"D55".equals(cfg.aligner)){throw new RuntimeException("AssignmentBinding requires cfg.aligner==\"D55\"; got '"+cfg.aligner+"'.");}
			if(!cfg.consensusRefSha256.equals(sidecar.consensusSha256)){
				throw new RuntimeException("AssignmentBinding consensus mismatch: cfg="+cfg.consensusRefSha256+" sidecar="+sidecar.consensusSha256+".");
			}
			if(!cfg.coveringSetsSha256.equals(sidecar.kmerSetsSha256)){
				throw new RuntimeException("AssignmentBinding covering-set mismatch: cfg="+cfg.coveringSetsSha256+" sidecar="+sidecar.kmerSetsSha256+".");
			}
			if(cfg.coveringK!=sidecar.coveringK){
				throw new RuntimeException("AssignmentBinding covering_k mismatch: cfg="+cfg.coveringK+" sidecar="+sidecar.coveringK+".");
			}
			//Explicit header-alphabet bind (hash alone does not imply the header alphabet -- Yoimiya). sidecar.coveringAlphabet
			//is a ReducedAlphabetSeedAssay.Alphabet; its .name field round-trips the header string for the production
			//uppercase-letters alphabet (ReducedAlphabet.parse maps a non-preset letters spec to itself).
			if(!cfg.coveringAlphabet.equals(sidecar.coveringAlphabet.name)){
				throw new RuntimeException("AssignmentBinding covering_alphabet mismatch: cfg='"+cfg.coveringAlphabet+
					"' sidecar='"+sidecar.coveringAlphabet.name+"'.");
			}
			if(!"NONE".equals(cfg.kmerBoundary)){
				throw new RuntimeException("AssignmentBinding requires cfg.kmerBoundary==\"NONE\" "+
					"(Boundary.NONE, the mode scoreFamiliesF4 uses); got '"+
					cfg.kmerBoundary+"'.");
			}
			//5. Load and own the HBM bundle from the exact roster, consensus residues, and runtime-semantic
			//provenance already bound by cfg. The loader clones consensus bytes and exposes no graph alias.
			final ArrayList<String> roster=new ArrayList<String>(nF);
			final HashMap<String,byte[]> consensusByRep=new HashMap<String,byte[]>(nF*2);
			for(final ProteinSequence target : loaded){roster.add(target.id); consensusByRep.put(target.id,target.enc);}
			final byte[][] semanticProvenance=HbmBundleLoader.loadSemanticProvenance(hbmProvenancePath);
			final byte[][] trustedRuntimeHashes=Arrays.copyOf(semanticProvenance,HbmBundleFormat.RUNTIME_SEMANTIC_COUNT);
			final HbmBundleLoader.Loaded hbm;
			try{
				hbm=HbmBundleLoader.load(java.nio.file.Paths.get(hbmBundlePath),roster,
					new HbmBundleLoader.ConsensusProvider(){
						@Override public byte[] consensusFor(final String repId){return consensusByRep.get(repId);}
					},trustedRuntimeHashes,hbmMinCount);
			}catch(java.io.IOException e){throw new RuntimeException("Failed to load verified HBM bundle: "+hbmBundlePath,e);}
			if(hbm.familyCount()!=nF){throw new RuntimeException("AssignmentBinding HBM family-count mismatch: hbm="+hbm.familyCount()+" cfg="+nF);}
			for(int i=0; i<nF; i++){
				if(!hbm.repId(i).equals(cfg.repId(i))){throw new RuntimeException("AssignmentBinding HBM roster mismatch at rank "+i+
					": hbm='"+hbm.repId(i)+"' cfg='"+cfg.repId(i)+"'.");}
			}
			//6. Concurrent-edit detection (Yoimiya, 2026-09-14): cfg.load hashed the files, then the object loaders read them AGAIN;
			//re-hash now and compare to cfg's recorded hashes, so a file rewritten DURING this factory's loads is caught. A
			//final hash comparison to cfg is sufficient for ordinary concurrent-edit detection (no broader framework needed).
			verifyUnchanged(familyListPath, cfg.familyListSha256, "family list");
			verifyUnchanged(consensusPath, cfg.consensusRefSha256, "consensus");
			verifyUnchanged(roleManifestPath, cfg.roleManifestSha256, "role manifest");
			verifyUnchanged(coveringSetsPath, cfg.coveringSetsSha256, "covering sets");
			verifyUnchanged(sidecarPath, cfg.shortlistSidecarSha256, "sidecar");
			verifyUnchanged(hbmBundlePath,cfg.hbmBundleSha256,"HBM bundle");
			verifyUnchanged(hbmProvenancePath,cfg.hbmSemanticProvenanceSha256,"HBM semantic provenance");
			//7. Own the targets privately: an unmodifiable copy, so the returned binding exposes no mutable resource-graph alias.
			final List<ProteinSequence> owned=Collections.unmodifiableList(new ArrayList<ProteinSequence>(loaded));
			return new AssignmentBinding(sidecar, owned, cfg,hbm,null,null);
		}

		/**
		 * Loads the canonical schema-7 production binding, including the frozen
		 * per-family core coordinates required by the mutual-overlap gate.
		 */
		public static AssignmentBinding loadSchema7(final String artifactPath,
				final String expectedArtifactSha80,
				final String rosterPath, final String consensusPath,
				final String roleManifestPath, final String coreCoordinatesPath,
				final String coveringSetsPath, final String sidecarPath,
				final String hbmBundlePath, final String hbmProvenancePath){
			return loadSchema7(artifactPath,expectedArtifactSha80,rosterPath,
				consensusPath,roleManifestPath,coreCoordinatesPath,coveringSetsPath,
				sidecarPath,hbmBundlePath,hbmProvenancePath,1);
		}

		public static AssignmentBinding loadSchema7(final String artifactPath,
				final String expectedArtifactSha80,
				final String rosterPath, final String consensusPath,
				final String roleManifestPath, final String coreCoordinatesPath,
				final String coveringSetsPath, final String sidecarPath,
				final String hbmBundlePath, final String hbmProvenancePath,
				final int hbmMinCount){
			final FamilyAcceptanceConfig cfg=FamilyAcceptanceConfig.loadSchema7(
				artifactPath,expectedArtifactSha80,rosterPath,consensusPath,roleManifestPath,
				coreCoordinatesPath,coveringSetsPath,sidecarPath,hbmBundlePath,
				hbmProvenancePath);
			final FamilyShortlistSidecar sidecar=FamilyShortlistSidecar.load(
				sidecarPath,consensusPath,coveringSetsPath);
			final List<ProteinSequence> loaded=ProteinSearch.readFasta(consensusPath);
			final int nF=cfg.nFamilies;
			if(sidecar.nFamilies!=nF || loaded.size()!=nF){
				throw new RuntimeException("Schema-7 binding family-count mismatch: cfg="+
					nF+" sidecar="+sidecar.nFamilies+" targets="+loaded.size());
			}
			final int[] consensusLengths=new int[nF];
			for(int rank=0; rank<nF; rank++){
				final String target=loaded.get(rank).id;
				if(!target.equals(sidecar.repIds[rank]) || !target.equals(cfg.repId(rank))){
					throw new RuntimeException("Schema-7 binding roster mismatch at rank "+rank+".");
				}
				consensusLengths[rank]=loaded.get(rank).length();
			}
			if(!cfg.consensusRefSha256.equals(sidecar.consensusSha256) ||
					!cfg.coveringSetsSha256.equals(sidecar.kmerSetsSha256) ||
					cfg.coveringK!=sidecar.coveringK ||
					!cfg.coveringAlphabet.equals(sidecar.coveringAlphabet.name)){
				throw new RuntimeException("Schema-7 sidecar/config resource mismatch.");
			}
			final FamilyCoreCoordinates core=FamilyCoreCoordinates.load(
				coreCoordinatesPath,rosterPath,consensusPath,sidecar.repIds,
				consensusLengths,cfg.coreFamilyListSha80);
			if(!core.artifactSha80.equals(cfg.coreCoordinatesSha80)){
				throw new RuntimeException("Schema-7 core-coordinate artifact hash mismatch.");
			}
			final ArrayList<String> roster=new ArrayList<String>(nF);
			final HashMap<String,byte[]> consensusByRep=new HashMap<String,byte[]>(nF*2);
			for(final ProteinSequence target : loaded){
				roster.add(target.id); consensusByRep.put(target.id,target.enc);
			}
			final byte[][] semantic=HbmBundleLoader.loadSemanticProvenance(hbmProvenancePath);
			final byte[][] trustedRuntime=Arrays.copyOf(semantic,
				HbmBundleFormat.RUNTIME_SEMANTIC_COUNT);
			final HbmBundleLoader.Loaded hbm;
			try{
				hbm=HbmBundleLoader.load(java.nio.file.Paths.get(hbmBundlePath),roster,
					new HbmBundleLoader.ConsensusProvider(){
						@Override public byte[] consensusFor(final String repId){
							return consensusByRep.get(repId);
						}
					},trustedRuntime,hbmMinCount);
			}catch(java.io.IOException e){
				throw new RuntimeException("Failed to load schema-7 HBM bundle: "+
					hbmBundlePath,e);
			}
			for(int rank=0; rank<nF; rank++){
				if(!hbm.repId(rank).equals(cfg.repId(rank))){
					throw new RuntimeException("Schema-7 HBM roster mismatch at rank "+rank+".");
				}
			}
			if(!cfg.hbmScoreContract.equals(hbm.scoreContract())){
				throw new RuntimeException("Schema-7 HBM score-contract mismatch: config='"+
					cfg.hbmScoreContract+"' runtime='"+hbm.scoreContract()+"'.");
			}
			verifyUnchanged(rosterPath,cfg.familyListSha256,"schema-7 roster",true);
			verifyUnchanged(consensusPath,cfg.consensusRefSha256,"schema-7 consensus",true);
			verifyUnchanged(roleManifestPath,cfg.roleManifestSha256,"schema-7 role manifest",true);
			verifyUnchanged(coreCoordinatesPath,cfg.coreCoordinatesSha80,"schema-7 core coordinates",true);
			verifyUnchanged(coveringSetsPath,cfg.coveringSetsSha256,"schema-7 covering sets",true);
			verifyUnchanged(sidecarPath,cfg.shortlistSidecarSha256,"schema-7 sidecar",true);
			verifyUnchanged(hbmBundlePath,cfg.hbmBundleSha256,"schema-7 HBM bundle");
			verifyUnchanged(hbmProvenancePath,cfg.hbmSemanticProvenanceSha256,
				"schema-7 HBM provenance",true);
			final List<ProteinSequence> owned=Collections.unmodifiableList(
				new ArrayList<ProteinSequence>(loaded));
			return new AssignmentBinding(sidecar,owned,cfg,hbm,null,core);
		}

		/** Loads the same consensus/shortlist resource graph as production assignment,
		 *  but binds a fixed membership-independent construction profile and no HBM.
		 *  The roster, sidecar and consensus must have identical rank order. */
		public static AssignmentBinding loadConstruction(final String familyListPath,
				final String consensusPath, final String coveringSetsPath, final String sidecarPath,
				final String coreCoordinatesPath,
				final FixedConstructionProfile profile){
			if(profile==null){throw new IllegalArgumentException("Construction profile must not be null");}
			final String rosterBefore=DigestSuffix.file(familyListPath),consensusBefore=DigestSuffix.file(consensusPath),
				coveringBefore=DigestSuffix.file(coveringSetsPath),sidecarBefore=DigestSuffix.file(sidecarPath);
			final HbmRosterLoader roster=HbmRosterLoader.load(familyListPath);
			final FamilyShortlistSidecar sidecar=FamilyShortlistSidecar.load(sidecarPath,consensusPath,coveringSetsPath);
			final List<ProteinSequence> loaded=ProteinSearch.readFasta(consensusPath);
			if(roster.size()!=sidecar.nFamilies || loaded.size()!=sidecar.nFamilies){
				throw new RuntimeException("Construction binding family-count mismatch: roster="+roster.size()+
					" sidecar="+sidecar.nFamilies+" targets="+loaded.size());
			}
			for(int i=0; i<sidecar.nFamilies; i++){
				final String rid=roster.repIdByRank[i],sid=sidecar.repIds[i],tid=loaded.get(i).id;
				if(!rid.equals(sid) || !sid.equals(tid)){throw new RuntimeException("Construction binding roster mismatch at rank "+i+
					": roster='"+rid+"' sidecar='"+sid+"' consensus='"+tid+"'.");}
				requireAsciiToken(rid,"construction rep_id at rank "+i);
			}
			final int[] consensusLengths=new int[loaded.size()];
			for(int i=0; i<loaded.size(); i++){consensusLengths[i]=loaded.get(i).length();}
			final String semanticFamilySha=roster.sourceRosterSha80==null?rosterBefore:roster.sourceRosterSha80;
			final FamilyCoreCoordinates core=FamilyCoreCoordinates.load(coreCoordinatesPath,familyListPath,
				consensusPath,roster.repIdByRank,consensusLengths,semanticFamilySha);
			verifyUnchanged(familyListPath,rosterBefore,"construction family list");
			verifyUnchanged(consensusPath,consensusBefore,"construction consensus");
			verifyUnchanged(coveringSetsPath,coveringBefore,"construction covering sets");
			verifyUnchanged(sidecarPath,sidecarBefore,"construction sidecar");
			final List<ProteinSequence> owned=Collections.unmodifiableList(new ArrayList<ProteinSequence>(loaded));
			return new AssignmentBinding(sidecar,owned,null,null,profile,core);
		}

		private static void requireAsciiToken(final String value, final String label){
			if(value.isEmpty()){throw new RuntimeException(label+" is empty");}
			for(int i=0; i<value.length(); i++){
				final char c=value.charAt(i);
				if(c>127 || c==' ' || c=='\t' || c=='\r' || c=='\n'){
					throw new RuntimeException(label+" must be an ASCII non-whitespace token: '"+value+"'");
				}
			}
		}
		private static void verifyUnchanged(final String path, final String expectedSha256, final String label){
			verifyUnchanged(path, expectedSha256, label, false);
		}
		/** Text transport compression preserves content pins; HBM remains a stored-byte pin. */
		private static void verifyUnchanged(final String path, final String expectedSha256, final String label, final boolean text){
			final String live=text ? MagQCTextResource.sha80(path) : DigestSuffix.file(path);
			if(!live.equals(expectedSha256)){
				throw new RuntimeException("AssignmentBinding: "+label+" file '"+path+"' changed during binding load (live "+
					live+" != config-recorded "+expectedSha256+") -- a concurrent edit; refusing to bind stale/rewritten resources.");
			}
		}
	}

	/**
	 * Canonical zero-or-one family assignment for one query (§E). Stage order (pinned; see the patch doc):
	 * F4 score (kmerNorm=none, the fixed canonical baseline) -> empty-Q guard -> preselect (length + k-mer count + k-mer
	 * query density) BEFORE topN -> topN (at most min(50,nFamilies), D50) -> D55-align each survivor -> post-align
	 * acceptance (raw score AND identity AND overlap) -> policy. Returns the one assigned family or a typed rejection;
	 * never more than one, never a silent drop. Uses the D55 aligner (this searcher's {@code aligner} mode must be
	 * {@code "d55"}) and the shared {@link #d55Metrics} so the decision inputs are IDENTICAL to the calibration metrics
	 * (Patch 1). Cutoffs come from {@code binding.cfg} (Eru's accepted config); this method invents no VALUES. The
	 * {@link AssignmentBinding} was loaded+owned by {@link AssignmentBinding#load}, so its resource/roster guarantees
	 * cannot be bypassed by a caller.
	 *
	 * @param binding A loaded config+sidecar+targets binding.
	 * @param q Query protein.
	 * @param scratch Caller-provided thread-local scratch (reused; carries {@link ShortlistScratch#d55m}/{@link ShortlistScratch#qSize}).
	 * @param policy FIRST_CLEAR or SCORE_ALL (non-null); bounded lookahead uses the overload with X.
	 * @return The assigned family, or a typed rejection.
	 */
	public FamilyAssignment assignFamily(final AssignmentBinding binding, final ProteinSequence q,
			final ShortlistScratch scratch, final AssignPolicy policy){
		return assignFamily(binding,q,scratch,policy,0);
	}

	/**
	 * Canonical assignment with explicit bounded-lookahead distance. For
	 * {@link AssignPolicy#BOUNDED_LOOKAHEAD}, each accepted candidate whose normalized
	 * D55 R strictly exceeds the running best resets the exclusive stop boundary to
	 * {@code min(cap,rank+lookahead)}, where rank is one-based. Equal-R rep-id tie
	 * improvements can change the winner but do not reset the boundary. The method
	 * therefore evaluates X candidates after the last strict R advance and chooses
	 * the greatest R, breaking ties by ascending rep_id. If no family passes, it
	 * evaluates the complete capped shortlist. Other policies require
	 * {@code lookahead==0}.
	 */
	public FamilyAssignment assignFamily(final AssignmentBinding binding, final ProteinSequence q,
			final ShortlistScratch scratch, final AssignPolicy policy, final int lookahead){
		final FamilyShortlistSidecar sidecar=binding.sidecar;
		final List<ProteinSequence> targets=binding.targets;
		final FamilyAcceptanceConfig cfg=binding.cfg;
		final FixedConstructionProfile construction=binding.construction;
		final int nF=sidecar.nFamilies;

		//Per-query guards as RUNTIME throws (fire even under -da) -- a null policy, wrong aligner mode, or mis-sized
		//scratch must never silently produce a wrong assignment.
		if(policy==null){throw new RuntimeException("assignFamily: policy must not be null.");}
		validateLookahead(policy,lookahead);
		if(!aligner.equals("d55")){throw new RuntimeException("assignFamily requires this searcher's aligner mode==\"d55\" "+
			"(the D55 aligner of record; alignWith dispatches on it); got '"+aligner+"'.");}
		if(scratch.counts.length!=sidecar.dims || scratch.centered.length!=sidecar.dims || scratch.qVec.length!=sidecar.dims){
			throw new RuntimeException("assignFamily scratch composition-dim mismatch: scratch dims={"+scratch.counts.length+
				","+scratch.centered.length+","+scratch.qVec.length+"} but sidecar.dims="+sidecar.dims+".");
		}
		if(scratch.compScores.length!=nF || scratch.f4.length!=nF || scratch.sk.length!=nF){
			throw new RuntimeException("assignFamily scratch family-array length mismatch vs nFamilies="+nF+": compScores="+
				scratch.compScores.length+" f4="+scratch.f4.length+" sk="+scratch.sk.length+".");
		}
		final int maxCandidates=Math.min(50, nF);//D50 (Brian, verbatim, DECISIONS_v1 D50): "Repeat until you hit 50 OR one
		//is accepted" -- the canonical assigner visits AT MOST 50 candidates in rank order (METHODS_v1: "maximum 50").
		if(scratch.topIdx.length<1 || scratch.topIdx.length>maxCandidates){
			throw new RuntimeException("assignFamily scratch topIdx.length must be in [1,min(50,nFamilies)="+maxCandidates+
				"], got "+scratch.topIdx.length+" (D50: visit at most 50 candidates in rank order).");
		}
		if(construction!=null && policy!=AssignPolicy.SCORE_ALL){
			throw new RuntimeException("Fixed consensus construction requires SCORE_ALL; got "+policy+".");
		}
		if(construction==null && !cfg.hbmEnabled){throw new RuntimeException("Production assignFamily requires HBM gating.");}
		if(construction==null && cfg.isSchema7()!=(binding.core!=null)){
			throw new RuntimeException("Schema-7 assignment requires its bound core coordinates.");
		}
		if(construction!=null){scratch.constructionTrace.clear();}
		scratch.lastAlignedCandidates=scratch.lastTracebackAlignments=
			scratch.lastAcceptedCandidates=0;

		//Stage 1: F4 score (fills scratch.f4/sk/compScores + scratch.qSize). kmerNorm PINNED to the canonical baseline
		//(KMERNORM_NONE); acceptance is on RAW s_k / s_k/Qsize (kmerNorm-independent), and sub/cs300 are separate declared
		//comparison runs, never the production default (Yoimiya, 2026-09-14).
		scoreFamiliesF4(q, sidecar, scratch, KMERNORM_NONE);

		//Empty-Q (Qsize==0 -> density s_k/Qsize would be 0/0): D117 (Brian, 2026-09-15) -- for alphabet=20,
		//k=5, boundary=NONE, a called protein with no valid covering-set 5-mer is discarded from family
		//assignment as a typed rejection, contributes NO family count, and the assembly continues. This
		//is NOT a guessed density value and NOT a silent skip: it is the same typed-rejection mechanism as
		//NO_PRESELECT_SURVIVOR/NO_ACCEPTANCE_PASS below, so every caller (FamilyAssignCLI, GeneHitsAssigner)
		//already has a code path for "rejected, not assigned" and needs no new fallback logic. Training and
		//inference share this exact behavior (both call this same method). Formerly crashed loud pending
		//this decision (root_acceptance_v7.md gate 2); see AssignFamilyTest for the regression.
		if(scratch.qSize==0){
			return FamilyAssignment.rejected(AssignReason.NO_VALID_KMER);
		}
		final int Lq=q.length(), Qsize=scratch.qSize;

		//Stage 2: cheap preselect mask BEFORE topN -- absolute Lq window (D107, NOT a length RATIO) AND s_k count AND
		//s_k/Qsize query-density. A family failing any can never be the assignment, so mask it (f4=-inf). Cutoffs via
		//Eru's config ACCESSORS (its arrays are private). Old sidecar per-family thresholds are never read here.
		int survivors=0;
		for(int r=0; r<nF; r++){
			final int sk=scratch.sk[r];
			final boolean lenOk;
			final boolean countOk,qdenOk;
			if(construction==null){
				lenOk=(Lq>=cfg.lqLo(r) && Lq<=cfg.lqHi(r));
				countOk=(sk>=cfg.minKmerCount(r));
				qdenOk=(((double)sk/Qsize)>=cfg.minKmerQdensity(r));
			}else{
				lenOk=constructionLengthPass(binding,Lq,r);
				countOk=qdenOk=true;//D119: kmers retrieve/rank; they are not construction acceptance gates.
			}
			if(lenOk && countOk && qdenOk){survivors++;}
			else{scratch.f4[r]=Float.NEGATIVE_INFINITY;}
		}
		if(construction!=null){scratch.constructionTrace.lengthSurvivors=survivors;}
		if(survivors==0){return FamilyAssignment.rejected(AssignReason.NO_PRESELECT_SURVIVOR);}
		final int cap=Math.min(survivors, scratch.topIdx.length);
		topN(scratch.f4, cap, scratch.topIdx);//topIdx[0..cap) = survivors, F4 desc; F4 ties -> index asc

		//Stages 4-5-6: align each survivor (D55), post-align D107 acceptance, then policy.
		int bestIdx=-1; double bestR=Double.NEGATIVE_INFINITY; String bestRep=null;
		int bestA=0, bestSk=0; double bestId=0, bestOv=0, bestQden=0;
		float bestHbmRaw=0,bestHbmPath=0;
		int stopExclusive=cap;
		for(int i=0; i<cap; i++){
			final int idx=scratch.topIdx[i];
			final ProteinSequence t=targets.get(idx);
			scratch.lastAlignedCandidates++;
			final AAAlignment aln=(construction==null ? alignWith(q.enc,t.enc) : null);
			if(construction!=null){
				alignConstructionCandidate(binding,q,idx,scratch);
				continue;
			}
			if(aln==null){
				throw new RuntimeException("assignFamily: D55 alignment unexpectedly null for query='"+q.id+"' vs target='"+
					t.id+"' -- GlocalAminoLinear.align must never return null for nonempty input.");
			}
			d55Metrics(q, t, aln, scratch.d55m);
			final int A=scratch.d55m.A;
			final double id=aln.pident(),wholeOverlap=scratch.d55m.overlap,
				R=scratch.d55m.R;
			final int sk=scratch.sk[idx];
			final double qden=(double)sk/Qsize;
			if(cfg.isSchema7()){
				final boolean lengthPass=Lq>=cfg.lqLo(idx) && Lq<=cfg.lqHi(idx);
				if(!FamilyAcceptanceGate.passesBeforeHbm(cfg.acceptanceEnabled(idx),
						lengthPass,A,id,R,0,
						sk,qden,cfg.minRawScore(idx),cfg.minId(idx),cfg.minR(idx),
						FamilyAcceptanceGate.PASS_ALL_DOUBLE,cfg.minKmerCount(idx),
						cfg.minKmerQdensity(idx))){
					if(i+1>=stopExclusive){break;}else{continue;}
				}
			}else if(A<cfg.minRawScore(idx) || id<cfg.minId(idx) ||
					wholeOverlap<cfg.minOverlap(idx)){
				if(i+1>=stopExclusive){break;}else{continue;}
			}
			//D118: raw HBM is the final inclusive GATE only. Traceback is generated only for the rare
			//candidate that already passed every length/kmer/D55 base gate, then the original D55 R remains
			//the cross-family winner metric. The repeated D55 call is deliberate: path recording is avoided
			//for base-filter failures, as required by the final-conditional-scorer design.
			final AAAlignment pathAln=GlocalAminoLinear.align(q.enc,t.enc,true);
			scratch.lastTracebackAlignments++;
			if(pathAln==null || pathAln.match==null){throw new RuntimeException("assignFamily: path-recorded D55 alignment missing for query='"+
				q.id+"' vs target='"+t.id+"'.");}
			if(!sameAlignment(aln,pathAln)){throw new RuntimeException("assignFamily: D55 path recording changed alignment metrics for query='"+
				q.id+"' vs target='"+t.id+"'; refusing inconsistent base/HBM scoring.");}
			final double acceptedOverlap;
			if(cfg.isSchema7()){
				constructionCoreMetrics(pathAln,q.length(),binding.core.coreStart[idx],
					binding.core.coreEnd[idx],binding.core.coreLength[idx],scratch.d55m);
				acceptedOverlap=scratch.d55m.mutualOverlap;
				if(acceptedOverlap<cfg.minCoreMutualOverlap(idx)){
					if(i+1>=stopExclusive){break;}else{continue;}
				}
			}else{acceptedOverlap=wholeOverlap;}
			final float[] hbmScore=binding.hbm.score(idx,q.enc,pathAln,true);
			if(hbmScore==null || hbmScore.length!=2 || !Float.isFinite(hbmScore[0]) ||
					(cfg.isSchema7() && !Float.isFinite(hbmScore[1]))){
				throw new RuntimeException("assignFamily: invalid HBM score for query='"+q.id+"' vs target='"+t.id+"'.");
			}
			final float hbmRaw=hbmScore[0],hbmPath=hbmScore[1];
			if(cfg.isSchema7()){
				final boolean lengthPass=Lq>=cfg.lqLo(idx) && Lq<=cfg.lqHi(idx);
				if(!FamilyAcceptanceGate.passes(cfg.acceptanceEnabled(idx),lengthPass,
					A,id,R,acceptedOverlap,sk,qden,hbmPath,cfg.minRawScore(idx),
					cfg.minId(idx),cfg.minR(idx),cfg.minCoreMutualOverlap(idx),
					cfg.minKmerCount(idx),cfg.minKmerQdensity(idx),
					cfg.minHbmPathRelative(idx))){
					if(i+1>=stopExclusive){break;}else{continue;}
				}
			}else if(hbmRaw<cfg.minHbmRaw99(idx)){
				if(i+1>=stopExclusive){break;}else{continue;}
			}
			scratch.lastAcceptedCandidates++;
			if(policy==AssignPolicy.FIRST_CLEAR){
				return acceptedCandidate(cfg,idx,sidecar.repIds[idx],A,id,
					acceptedOverlap,R,sk,qden,hbmRaw,hbmPath);
			}
			//Non-greedy policies keep greatest R; deterministic tie-break R desc then rep_id ascending.
			final String rep=sidecar.repIds[idx];
			final double priorBestR=bestR;
			if(policy==AssignPolicy.BOUNDED_LOOKAHEAD){
				stopExclusive=boundedStopExclusive(stopExclusive,i,cap,lookahead,R,priorBestR);
			}
			if(R>priorBestR || (R==priorBestR && (bestRep==null || rep.compareTo(bestRep)<0))){
				bestR=R; bestIdx=idx; bestRep=rep; bestA=A; bestId=id;
				bestOv=acceptedOverlap; bestSk=sk; bestQden=qden;
				bestHbmRaw=hbmRaw; bestHbmPath=hbmPath;
			}
			if(i+1>=stopExclusive){break;}
		}
		if(construction!=null){
			final ConstructionTrace trace=scratch.constructionTrace;
			final int strictSlot=ConstructionTrace.PROFILE_OFFSET+FixedConstructionProfile.PROFILE_COUNT-1;
			final int idx=trace.familyIdx[strictSlot];
			if(idx<0){return FamilyAssignment.rejected(AssignReason.NO_ACCEPTANCE_PASS);}
			return FamilyAssignment.assigned(idx,sidecar.repIds[idx],trace.rawScore[strictSlot],
				trace.identity[strictSlot],trace.mutualOverlap[strictSlot],trace.R[strictSlot],
				scratch.sk[idx],(double)scratch.sk[idx]/Qsize,0);
		}
		if(bestIdx<0){return FamilyAssignment.rejected(AssignReason.NO_ACCEPTANCE_PASS);}
		return acceptedCandidate(cfg,bestIdx,bestRep,bestA,bestId,bestOv,bestR,
			bestSk,bestQden,bestHbmRaw,bestHbmPath);
	}

	/** Measures recall of the production F4 top-50 shortlist against an exhaustive
	 *  fixed-construction search. Both arms share the exact construction length,
	 *  D55/core-overlap, acceptance, R ordering, and tie-break implementations used
	 *  by {@link #assignFamily}; the only difference is candidate extent. */
	public ConstructionRecallResult measureConstructionRecall(final AssignmentBinding binding,
			final ProteinSequence q, final ShortlistScratch scratch){
		if(binding==null || q==null || scratch==null){
			throw new IllegalArgumentException("measureConstructionRecall requires non-null binding, query, and scratch");
		}
		if(binding.construction==null){throw new IllegalArgumentException("measureConstructionRecall requires a fixed-construction binding");}
		if(!aligner.equals("d55")){throw new IllegalArgumentException("measureConstructionRecall requires aligner=d55, got "+aligner);}
		final FamilyShortlistSidecar sidecar=binding.sidecar;
		final int nF=sidecar.nFamilies,Lq=q.length();
		if(scratch.compScores.length!=nF || scratch.f4.length!=nF || scratch.sk.length!=nF || scratch.topIdx.length>50){
			throw new IllegalArgumentException("measureConstructionRecall received incompatible shortlist scratch");
		}
		final ConstructionTrace trace=scratch.constructionTrace;
		trace.clear();
		scoreFamiliesF4(q,sidecar,scratch,KMERNORM_NONE);
		if(scratch.qSize==0){
			return emptyConstructionRecall(AssignReason.NO_VALID_KMER,0,0,0);
		}
		int survivors=0;
		for(int r=0;r<nF;r++){
			if(constructionLengthPass(binding,Lq,r)){survivors++;}
			else{scratch.f4[r]=Float.NEGATIVE_INFINITY;}
		}
		trace.lengthSurvivors=survivors;
		if(survivors==0){return emptyConstructionRecall(AssignReason.NO_PRESELECT_SURVIVOR,0,0,0);}
		final int cap=Math.min(survivors,scratch.topIdx.length);
		topN(scratch.f4,cap,scratch.topIdx);

		int accepted=0,winner=-1,runner=-1,winnerA=0,runnerA=0;
		float winnerR=Float.NEGATIVE_INFINITY,runnerR=Float.NEGATIVE_INFINITY;
		float winnerId=Float.NaN,runnerId=Float.NaN,winnerMutual=Float.NaN,runnerMutual=Float.NaN;
		final int strict=FixedConstructionProfile.PROFILE_COUNT-1;
		for(int r=0;r<nF;r++){
			if(!constructionLengthPass(binding,Lq,r)){continue;}
			alignConstructionCandidate(binding,q,r,scratch);
			final float R=scratch.constructionR,id=scratch.constructionIdentity,mutual=scratch.constructionMutualOverlap;
			if(!constructionProfilePass(binding.construction,strict,R,id,mutual)){continue;}
			accepted++;
			if(betterConstruction(r,R,winner,winnerR,sidecar.repIds)){
				runner=winner;runnerA=winnerA;runnerR=winnerR;runnerId=winnerId;runnerMutual=winnerMutual;
				winner=r;winnerA=scratch.constructionRawScore;winnerR=R;winnerId=id;winnerMutual=mutual;
			}else if(betterConstruction(r,R,runner,runnerR,sidecar.repIds)){
				runner=r;runnerA=scratch.constructionRawScore;runnerR=R;runnerId=id;runnerMutual=mutual;
			}
		}
		if(trace.alignedCandidates!=survivors){
			throw new RuntimeException("Exhaustive construction aligned "+trace.alignedCandidates+
				" candidates but length gate admitted "+survivors+" for query '"+q.id+"'");
		}
		final int strictSlot=ConstructionTrace.PROFILE_OFFSET+strict;
		if(trace.familyIdx[strictSlot]!=winner){
			throw new RuntimeException("Exhaustive strict winner diverged from canonical construction trace for query '"+q.id+"'");
		}
		if(winner<0){return emptyConstructionRecall(AssignReason.NO_ACCEPTANCE_PASS,cap,survivors,trace.alignedCandidates);}
		final int winnerRank=shortlistRank(scratch.topIdx,cap,winner),runnerRank=shortlistRank(scratch.topIdx,cap,runner);
		return new ConstructionRecallResult(AssignReason.ASSIGNED,cap,survivors,trace.alignedCandidates,accepted,
			winner,runner,winnerRank,runnerRank,winnerA,runnerA,winnerR,runnerR,winnerId,runnerId,
			winnerMutual,runnerMutual);
	}

	private static ConstructionRecallResult emptyConstructionRecall(final AssignReason reason,
			final int shortlistSize, final int lengthSurvivors, final int alignedCandidates){
		return new ConstructionRecallResult(reason,shortlistSize,lengthSurvivors,alignedCandidates,0,
			-1,-1,0,0,0,0,Float.NEGATIVE_INFINITY,Float.NEGATIVE_INFINITY,
			Float.NaN,Float.NaN,Float.NaN,Float.NaN);
	}

	private static int shortlistRank(final int[] ranks, final int size, final int target){
		if(target<0){return 0;}
		for(int i=0;i<size;i++){if(ranks[i]==target){return i+1;}}
		return 0;
	}

	private static boolean constructionLengthPass(final AssignmentBinding binding, final int queryLength,
			final int familyIdx){
		final float ratio=queryLength/(float)binding.core.coreLength[familyIdx];
		return ratio>=binding.construction.minLengthRatio && ratio<=binding.construction.maxLengthRatio;
	}

	private static void alignConstructionCandidate(final AssignmentBinding binding, final ProteinSequence q,
			final int idx, final ShortlistScratch scratch){
		final ProteinSequence t=binding.targets.get(idx);
		final AAAlignment aln=GlocalAminoLinear.align(q.enc,t.enc,true);
		if(aln==null || aln.match==null){throw new RuntimeException("Fixed-construction D55 path missing for query='"+
			q.id+"' target='"+t.id+"'");}
		d55Metrics(q,t,aln,scratch.d55m);
		constructionCoreMetrics(aln,q.length(),binding.core.coreStart[idx],binding.core.coreEnd[idx],
			binding.core.coreLength[idx],scratch.d55m);
		final float R=(float)scratch.d55m.R,id=(float)aln.pident(),covQ=(float)scratch.d55m.mutualCovQ,
			covT=(float)scratch.d55m.mutualCovT,mutual=(float)scratch.d55m.mutualOverlap;
		scratch.constructionIdx=idx;scratch.constructionRawScore=scratch.d55m.A;
		scratch.constructionR=R;scratch.constructionIdentity=id;scratch.constructionMutualCovQ=covQ;
		scratch.constructionMutualCovT=covT;scratch.constructionMutualOverlap=mutual;
		considerConstructionCandidate(scratch.constructionTrace,binding.construction,binding.sidecar.repIds,
			idx,scratch.d55m.A,binding.core.coreLength[idx],R,id,covQ,covT,mutual);
	}

	private static boolean constructionProfilePass(final FixedConstructionProfile profile, final int profileIdx,
			final float R, final float identity, final float mutual){
		return FamilyAcceptanceGate.passes(true,true,0,identity,R,mutual,0,0,0,
			FamilyAcceptanceGate.PASS_ALL_INT,profile.minIdentity[profileIdx],
			profile.minR,profile.minMutualOverlap,FamilyAcceptanceGate.PASS_ALL_INT,
			FamilyAcceptanceGate.PASS_ALL_DOUBLE,FamilyAcceptanceGate.PASS_ALL_FLOAT);
	}

	/** Updates the reusable construction trace without allocating per candidate. */
	private static void considerConstructionCandidate(final ConstructionTrace trace,
			final FixedConstructionProfile profile, final String[] repIds, final int idx, final int A,
			final int coreLength, final float R, final float identity, final float covQ,
			final float covT, final float mutual){
		trace.alignedCandidates++;
		if(betterConstruction(idx,R,trace.familyIdx[ConstructionTrace.BEST],trace.R[ConstructionTrace.BEST],repIds)){
			trace.copy(ConstructionTrace.BEST,ConstructionTrace.RUNNER_UP);
			trace.set(ConstructionTrace.BEST,idx,A,coreLength,R,identity,covQ,covT,mutual);
		}else if(betterConstruction(idx,R,trace.familyIdx[ConstructionTrace.RUNNER_UP],trace.R[ConstructionTrace.RUNNER_UP],repIds)){
			trace.set(ConstructionTrace.RUNNER_UP,idx,A,coreLength,R,identity,covQ,covT,mutual);
		}
		if(R<profile.minR){return;}
		trace.rSurvivors++;
		if(mutual<profile.minMutualOverlap){return;}
		trace.mutualSurvivors++;
		for(int p=0; p<FixedConstructionProfile.PROFILE_COUNT; p++){
			if(!constructionProfilePass(profile,p,R,identity,mutual)){continue;}
			trace.profileSurvivors[p]++;
			final int slot=ConstructionTrace.PROFILE_OFFSET+p;
			if(betterConstruction(idx,R,trace.familyIdx[slot],trace.R[slot],repIds)){
				trace.set(slot,idx,A,coreLength,R,identity,covQ,covT,mutual);
			}
		}
	}

	/** Construction coverage is measured only inside the frozen per-family core.
	 *  Whole-alignment identity and normalized D55 R remain unchanged. */
	static void constructionCoreMetrics(final AAAlignment aln, final int queryLength,
			final int coreStart, final int coreEnd, final int coreLength, final D55Metrics out){
		if(aln==null || aln.match==null){throw new RuntimeException("Construction core coverage requires a recorded D55 path");}
		if(queryLength<1 || coreStart<0 || coreEnd<coreStart || coreLength!=coreEnd-coreStart+1){
			throw new IllegalArgumentException("Invalid construction core dimensions");
		}
		int q=aln.qStart,ref=aln.tStart,pairedInside=0;
		for(final byte op : aln.match){
			if(op=='m'){
				if(ref>=coreStart && ref<=coreEnd){pairedInside++;}
				q++; ref++;
			}else if(op=='D'){
				ref++;
			}else if(op=='I'){
				q++;
			}else{throw new RuntimeException("Unknown D55 path operation: "+(char)op);}
			if(q>queryLength || ref>aln.tStop+1){throw new RuntimeException("D55 construction path exceeds declared spans");}
		}
		if(q!=queryLength || q!=aln.qStop+1 || ref!=aln.tStop+1 || pairedInside>queryLength || pairedInside>coreLength){
			throw new RuntimeException("D55 construction path does not consume declared spans or exceeds core bounds");
		}
		out.mutualCovQ=pairedInside/(double)queryLength;
		out.mutualCovT=pairedInside/(double)coreLength;
		out.mutualOverlap=pairedInside/(double)Math.max(queryLength,coreLength);
	}

	/** Construction uses canonical float32 R; rep IDs are ASCII-validated at binding load,
	 *  so String order is bytewise unsigned order for exact-R ties. */
	private static boolean betterConstruction(final int idx, final float R, final int oldIdx,
			final float oldR, final String[] repIds){
		return oldIdx<0 || R>oldR || (R==oldR && repIds[idx].compareTo(repIds[oldIdx])<0);
	}

	/** Path recording must not change the deterministic D55 optimum used by the base gate and HBM. */
	private static boolean sameAlignment(final AAAlignment a, final AAAlignment b){
		return a.rawScore==b.rawScore && a.qStart==b.qStart && a.qStop==b.qStop && a.tStart==b.tStart && a.tStop==b.tStop
			&& a.identities==b.identities && a.mismatches==b.mismatches && a.gapOpens==b.gapOpens && a.length==b.length;
	}

	/** Java-8-safe scalar top-N: score descending, index ascending, NaN last. */
	private static void topN(final float[] scores, final int n, final int[] outIdx){
		if(scores==null || outIdx==null || n<0 || n>scores.length || n>outIdx.length){
			throw new IllegalArgumentException("Invalid topN arguments: scores="+(scores==null ? -1 : scores.length)+
				" n="+n+" outIdx="+(outIdx==null ? -1 : outIdx.length));
		}
		if(n==0){return;}
		int size=0;
		for(int idx=0; idx<scores.length; idx++){
			if(size==n){
				final int worst=outIdx[size-1];
				if(!betterTop(scores[idx],idx,scores[worst],worst)){continue;}
			}
			int pos=(size<n ? size : size-1);
			while(pos>0 && betterTop(scores[idx],idx,scores[outIdx[pos-1]],outIdx[pos-1])){pos--;}
			if(size<n){size++;}
			for(int j=size-1; j>pos; j--){outIdx[j]=outIdx[j-1];}
			outIdx[pos]=idx;
		}
	}
	private static boolean betterTop(final float scoreA, final int indexA, final float scoreB, final int indexB){
		final boolean nanA=Float.isNaN(scoreA),nanB=Float.isNaN(scoreB);
		return nanB ? !nanA : (!nanA && (scoreA>scoreB || (scoreA==scoreB && indexA<indexB)));
	}

	/** Maps one acceptance-winning rank to a countable assignment or a typed non-counting intercept. */
	private static FamilyAssignment acceptedCandidate(final FamilyAcceptanceConfig cfg, final int idx,
			final String repId, final int rawScore, final double identity, final double overlap, final double R,
			final int kmerCount, final double kmerQdensity, final float hbmRaw,
			final float hbmPathRelative){
		if(cfg.isBait(idx)){
			return FamilyAssignment.baitIntercepted(idx,repId,rawScore,identity,
				overlap,R,kmerCount,kmerQdensity,hbmRaw,hbmPathRelative);
		}
		if(cfg.isTracked(idx) && !cfg.isNetworkFeature(idx)){
			return FamilyAssignment.maskedIntercepted(idx,repId,rawScore,identity,
				overlap,R,kmerCount,kmerQdensity,hbmRaw,hbmPathRelative);
		}
		if(!cfg.isTracked(idx)){
			throw new RuntimeException("Accepted family has no valid schema-7 role at rank "+idx+" rep_id='"+repId+"'.");
		}
		return FamilyAssignment.assigned(idx,repId,rawScore,identity,overlap,R,
			kmerCount,kmerQdensity,hbmRaw,hbmPathRelative);
	}

	/**
	 * A3 (2026-09-04, plans/PER_FAMILY_THRESHOLDS_v1.md sec 3): stages 1-3 of the triage --
	 * length-ratio window, dimer-cosine floor, covering-kmer floor -- as a survivor mask over
	 * {@code s.f4}, applied BEFORE {@link #topN(float[], int, int[])} so a rejected family can never be padded
	 * into the shortlist when fewer than {@link #shortlist} families survive (sec 1: "fewer only
	 * when fewer survive"). Rejected families get {@code s.f4[r]=Float.NEGATIVE_INFINITY} (always
	 * worse than any real, bounded F4 score) rather than NaN -- {@link #topN(float[], int, int[])}'s contract
	 * sorts NaN last but still fills the requested count with it if not enough real scores exist,
	 * which would silently pad with rejects; requesting exactly the survivor count with rejects at
	 * -Infinity gives exactly the survivors, correctly ranked, with no padding.
	 * @param q The query.
	 * @param sc The sidecar (must be schema_version&gt;=2 -- crashes loud otherwise, sec 3).
	 * @param s Caller-provided scratch, already filled by {@link #scoreFamiliesF4}.
	 * @return The number of families surviving stage 3 (capped at {@code s.topIdx.length}) --
	 *         the count to request from {@link #topN(float[], int, int[])}.
	 */
	private int applyTriageMask(final ProteinSequence q, final FamilyShortlistSidecar sc, final ShortlistScratch s){
		if(sc.schemaVersion<2){
			throw new RuntimeException("triage=t requires a sidecar v3 (schema_version>=2, built with "+
				"thresholds=) -- this sidecar is schema_version="+sc.schemaVersion+" (v2, no per-family thresholds).");
		}
		final int nFamilies=sc.nFamilies;
		final double lQ=q.length();
		long s1=0, s2=0, s3=0;
		for(int r=0; r<nFamilies; r++){
			final double ratio=lQ/sc.medLen[r];
			if(ratio<sc.rLo[r] || ratio>sc.rHi[r]){s.f4[r]=Float.NEGATIVE_INFINITY; continue;}
			s1++;
			if(s.compScores[r]<sc.minDimer[r]){s.f4[r]=Float.NEGATIVE_INFINITY; continue;}
			s2++;
			if(s.sk[r]<sc.minKmers[r]){s.f4[r]=Float.NEGATIVE_INFINITY; continue;}
			s3++;
		}
		triageStage1Survivors.addAndGet(s1); triageStage2Survivors.addAndGet(s2); triageStage3Survivors.addAndGet(s3);
		return (int)Math.min((long)s.topIdx.length, s3);
	}

	/**
	 * Shortlist-mode query scoring (Brian, 2026-09-03: "integrate the hybrid shortlist"):
	 * scores the query against every family by F4 = z(s_k) + 2*z(s_c) - d_L using {@link #sidecar}'s
	 * precomputed data, selects the top {@link #shortlist} families via the scalar
	 * {@link #topN(float[], int, int[])} (deterministic tie-break), and aligns only
	 * those -- replacing the unshortlisted path's "align every seed-passing candidate" behavior.
	 * The L-hard length window is intentionally NOT applied here (soft d_L only), matching the
	 * production recommendation in {@code results/hybrid_shortlist_v1.md} §2/§5.
	 *
	 * @param q The query.
	 * @param targets Full target list (read-only); MUST be in the same order as {@link #sidecar}'s
	 *        families (verified once, up front, in {@link #search(List, List)} by both family-count
	 *        AND per-index id equality against {@code sidecar.repIds} -- UMP45's review, 2026-09-03).
	 * @param totalDbResidues Sum of all target lengths, for the E-value search space.
	 * @param s Caller-provided, thread-local scratch (reused across calls on the same thread).
	 * @return This query's passing hits (not yet in the frozen total order).
	 */
	private List<ProteinHit> searchOneQueryShortlist(final ProteinSequence q, final List<ProteinSequence> targets,
			final long totalDbResidues, final ShortlistScratch s){
		final double searchSpace=(double)q.length()*(double)totalDbResidues;
		final FamilyShortlistSidecar sc=sidecar;

		scoreFamiliesF4(q, sc, s, kmerNorm);
		//triage=f: cap MUST equal s.topIdx.length exactly as before A3 (byte-identical gate,
		//plans/PER_FAMILY_THRESHOLDS_v1.md sec 3) -- none of the triage-only code below runs.
		final int cap=(triage ? applyTriageMask(q, sc, s) : s.topIdx.length);
		topN(s.f4, cap, s.topIdx);
		if(triage){triageAligned.addAndGet(cap);}
		//selfScore(q) is the same for every candidate family -- computed once per query, not
		//per candidate. Unused (0) when triage is off.
		final int selfScoreQ=(triage ? selfScore(q) : 0);

		final ArrayList<ProteinHit> qhits=new ArrayList<ProteinHit>();
		for(int ti=0; ti<cap; ti++){//NOT a for-each over the whole array: s.topIdx has stale
			//entries beyond `cap` whenever triage reduced the survivor count below its length.
			final int idx=s.topIdx[ti];
			final ProteinSequence t=targets.get(idx);
			final AAAlignment aln=alignWith(q.enc, t.enc);
			if(aln==null){continue;}
			if(aln.rawScore<minRawScore){continue;}
			final double pid=aln.pident();
			if(pid<minPident){continue;}
			if(minCoverage>0){
				final double qCov=aln.length/(double)q.length();
				final double tCov=aln.length/(double)t.length();
				if(qCov<minCoverage || tCov<minCoverage){continue;}
			}
			if(triage){
				//Stage 4 (sec 1, 4a-4d): each family's own APPLIED threshold (already max/min-
				//combined with the global floor at sidecar-build time, sec 2), alongside (not
				//instead of) minRawScore/minPident/minCoverage above.
				if(pid<sc.minId[idx]){continue;}
				final double ratio=aln.rawScore/(double)selfScoreQ;
				if(ratio<sc.minRatio[idx]){continue;}
				if(aln.rawScore<sc.minScore[idx]){continue;}
				final double covQ, covT;
				if(covDef.equals("span")){
					covQ=(aln.qStop-aln.qStart+1)/(double)q.length(); covT=(aln.tStop-aln.tStart+1)/(double)t.length();
				}else{
					covQ=aln.length/(double)q.length(); covT=aln.length/(double)t.length();
				}
				if(covQ<sc.minCovQ[idx] || covT<sc.minCovT[idx]){continue;}
			}
			final double e=aln.evalue(searchSpace);
			if(e>evalueCutoff){continue;}
			if(triage){triageHits.incrementAndGet();}
			qhits.add(new ProteinHit(q.id, t.id, aln, e, EVALUE_APPROXIMATE));
		}
		if(maxTargetSeqs<qhits.size()){
			Collections.sort(qhits, BY_SCORE_DESC);
			while(qhits.size()>maxTargetSeqs){qhits.remove(qhits.size()-1);}
		}
		return qhits;
	}

	/** Dispatches the shortlist search without changing historical aliases. The D55 recurrence
	 *  is the canonical production choice. The archived prototypes used different APIs/scoring
	 *  contracts and must not be substituted for the legacy aliases as a supposed missing import.
	 *  The unshortlisted legacy path calls {@link AAAligner#align} directly. */
	private AAAlignment alignWith(final byte[] q, final byte[] t){
		if(aligner.equals("aaaligner")){return AAAligner.align(q, t, false);}
		//d55: the canonical D55/D56 detailed aligner (query-global/reference-local BLOSUM62, LINEAR gap-4,
		//no gap-open) -- the production alignment of record (PATH_TO_FINISHED_PRODUCT_v1 §1.4/§2, D55).
		//Same entry signature as AAAligner.align; never null for nonempty input (GlocalAminoLinear asserts
		//m>0 && n>0). Added 2026-09-14 (UMP45) for the canonical-assigner integration.
		if(aligner.equals("d55")){return GlocalAminoLinear.align(q, t, false);}
		//Preserve these historical aliases for old callers; they are not the canonical D55 path.
		if(aligner.equals("glocal") || aligner.equals("blosum")){return AAAligner.alignGlocal(q, t, false);}
		throw new RuntimeException("Unknown aligner: "+aligner+" (expected aaaligner/d55/glocal/blosum)");
	}

	/**
	 * Convenience: single-pair search returning the best HSP or null.
	 * @param q Query sequence.
	 * @param t Target sequence.
	 * @param searchSpace Effective search space for the E-value.
	 * @return Best hit, or null if none scores positively.
	 */
	public ProteinHit searchPair(final ProteinSequence q, final ProteinSequence t,
			final double searchSpace){
		final AAAlignment aln=AAAligner.align(q.enc, t.enc);
		if(aln==null){return null;}
		return new ProteinHit(q.id, t.id, aln, aln.evalue(searchSpace), EVALUE_APPROXIMATE);
	}

	/** Builds the target k-mer index (kmer packed as a long -> target indices). */
	private LongLongListHashMap buildIndex(final List<ProteinSequence> targets){
		final LongLongListHashMap index=new LongLongListHashMap();
		for(int i=0; i<targets.size(); i++){
			final long[] kmerArray=kmerSet(targets.get(i).enc).toArray();
			for(long km : kmerArray){
				index.put(km, (long)i);//distinct per target, so no duplicate index per kmer
			}
		}
		return index;
	}

	/**
	 * Extracts the set of distinct k-mers from an encoded sequence. K-mers
	 * containing X (ambiguous) are skipped. Uses the amino8 reduced alphabet
	 * when {@link #reducedSeed} is set.
	 * @param enc Encoded residues.
	 * @return Set of packed k-mers.
	 */
	private LongHashSet kmerSet(final byte[] enc){
		final LongHashSet set=new LongHashSet();
		final int bits=reducedSeed ? 3 : 5;
		assert(bits*k<=62) : "k too large for packed k-mer: k="+k+", bits="+bits;
		if(enc.length<k){return set;}
		final long mask=(1L<<(bits*k))-1;
		long kmer=0;
		int valid=0;
		for(int i=0; i<enc.length; i++){
			final byte e=enc[i];
			int code;
			if(reducedSeed){
				code=(e>=0 && e<EXT_TO_AMINO8.length) ? EXT_TO_AMINO8[e] : -1;
			}else{
				code=Blosum62.isStandard(e) ? e : -1;
			}
			if(code<0){kmer=0; valid=0; continue;}//reset on X/ambiguous residue
			kmer=((kmer<<bits)|code)&mask;
			valid++;
			if(valid>=k){set.add(kmer);}
		}
		return set;
	}

	/** Fails loudly on duplicate identifiers within one input (frozen contract). */
	private static void checkDuplicateIds(final List<ProteinSequence> seqs, final String which){
		final HashSet<String> seen=new HashSet<String>();
		for(ProteinSequence s : seqs){
			if(!seen.add(s.id)){
				throw new RuntimeException("Duplicate "+which+" identifier: '"+s.id+"'.");
			}
		}
	}

	/** Orders hits by best-HSP score descending (for maxTargetSeqs culling). */
	private static final Comparator<ProteinHit> BY_SCORE_DESC=new Comparator<ProteinHit>(){
		@Override
		public int compare(ProteinHit a, ProteinHit b){
			if(a.bitscore!=b.bitscore){return a.bitscore>b.bitscore ? -1 : 1;}
			return a.target.compareTo(b.target);
		}
	};

	/**
	 * The frozen total output order: query asc, E-value asc, bitscore desc,
	 * target asc, tstart asc, qstart asc. Deterministic (a total order).
	 */
	private static final Comparator<ProteinHit> TOTAL_ORDER=new Comparator<ProteinHit>(){
		@Override
		public int compare(ProteinHit a, ProteinHit b){
			int c=a.query.compareTo(b.query);
			if(c!=0){return c;}
			if(a.evalue!=b.evalue){return a.evalue<b.evalue ? -1 : 1;}
			if(a.bitscore!=b.bitscore){return a.bitscore>b.bitscore ? -1 : 1;}
			c=a.target.compareTo(b.target);
			if(c!=0){return c;}
			if(a.tstart!=b.tstart){return a.tstart<b.tstart ? -1 : 1;}
			return Integer.compare(a.qstart, b.qstart);
		}
	};
}
