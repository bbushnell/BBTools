package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Bounded membership pilot, not the legacy highest-eligible-family assignment.
 * Selects the best HBM among F4's top200 (opt-in top100), then tests ONLY its consensus with BLOSUM.
 * Optional patience20 accepts only an eligible current winner after20 nonimprovements.
 * No fallback, paired-terminal requirement or reciprocal-coverage gate is used.
 * @author Sayu
 */
public final class HbmWinnerCalibration {
	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, HbmWinnerCalibration.class, false);
			try{runForTest(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Small, completely pinned input only; there is deliberately no full-corpus switch. */
	public static void runForTest(String[] args) throws Exception{
		final HashMap<String,String> o=HmmComparisonData.options(args);
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList("out", "expected", "t", "threads", "profile", "hbmscoreonly", "shortlist", "patience"));
		for(String key : INPUTS){allowed.add(key); allowed.add(key+"sha80");}
		for(String key : o.keySet()){require(allowed.contains(key), "Unknown winner-pilot option: "+key);}
		final int expected=Integer.parseInt(HmmComparisonData.required(o, "expected"));
		require(expected>0 && expected<=MAX_QUERIES, "Pilot requires1..4096 queries; full-corpus execution is not released");
		require(!o.containsKey("t") || !o.containsKey("threads") || o.get("t").equals(o.get("threads")), "t and threads must agree");
		final int threads=Integer.parseInt(o.getOrDefault("t", o.getOrDefault("threads", "1")));
		require(threads>=1 && threads<=8, "Bounded pilot allows1..8 workers");
		final String profileOption=o.getOrDefault("profile", "f");
		require(profileOption.equals("t") || profileOption.equals("f"), "profile must be t or f");
		final boolean profile=profileOption.equals("t");
		final String scoreOption=o.getOrDefault("hbmscoreonly", "f");
		require(scoreOption.equals("t") || scoreOption.equals("f"), "hbmscoreonly must be t or f");
		final boolean scoreOnly=scoreOption.equals("t");
		final int shortlist=Integer.parseInt(o.getOrDefault("shortlist", "200"));
		require(shortlist==100 || shortlist==SHORTLIST, "Bounded shortlist comparison allows only100 or200 candidates");
		final int patience=Integer.parseInt(o.getOrDefault("patience", "0"));
		require(patience==0 || patience==20, "Bounded early-stop comparison allows patience=0 or20");
		require(patience==0 || shortlist==SHORTLIST, "Patience20 is independent of the top100 experiment; use shortlist=200");
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Pilot output directory must be fresh: "+out);
		for(String key : INPUTS){
			require(Files.isRegularFile(Paths.get(HmmComparisonData.required(o, key))), "Pilot input must be a regular file: "+key);
			DigestSuffix.requireSuffix(HmmComparisonData.required(o, key+"sha80"), key+" input pin");
		}
		Shared.setThreads(threads); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
		final long setupStart=System.nanoTime();
		verify(o);
		final ProteinSequence[] queries=readQueries(o.get("in"), expected);
		final HbmTailAssign.Library library=new HbmTailAssign.Library(o, HbmTailProfileRefiner.TailPolicy.LEGACY_CONSTRUCTION);
		require(library.ids.length>=SHORTLIST, "Both comparison arms require the same roster of at least200 models");
		final Engine engine=new Engine(library.ids, library.refs, library.profiles);
		final double setupSeconds=seconds(setupStart);
		System.err.println("HBM_WINNER_PILOT_SETUP queries="+queries.length+" families="+library.ids.length+" seconds="+setupSeconds);
		final long searchStart=System.nanoTime();
		final Result[] results=new Result[queries.length];
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final int active=Math.min(threads, queries.length);
		try{
			for(int i=0; i<active; i++){
				final ProcessThread worker=new ProcessThread(i, active, queries, results, engine, library.sidecar, failure, profile, scoreOnly, shortlist, patience);
				worker.start(); workers.add(worker);
			}
		}catch(Throwable t){failure.compareAndSet(null, t);}
		finally{ThreadWaiter.waitForThreadsToFinish(workers);}
		if(failure.get()!=null){throw new RuntimeException("Winner pilot worker failed; no PASS published", failure.get());}
		final double searchSeconds=seconds(searchStart);
		for(Result result : results){
			require(result!=null && result.scored<=shortlist && result.blosumScored>=1 && result.blosumScored<=result.scored,
				"Every query must retain actual bounded HBM and BLOSUM work counts");
			require(patience==0 ? result.scored==shortlist && result.blosumScored==1
				: result.scored>=patience+1 && (result.scored==shortlist || result.metrics.accepted()),
				"Only an eligible current winner may stop early after20 nonimproving candidates");
		}
		verify(o);//Publication follows both input verifications and all successful workers.
		Files.createDirectory(out);
		writeResults(out, queries, results, library.ids, o, active, shortlist, setupSeconds, searchSeconds);
		if(profile){writeProfile(out, workers, searchSeconds);}
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_WINNER_CALIBRATION_PILOT_PASS"); HmmComparisonData.close(pass);
		System.err.println("HBM_WINNER_PILOT_PASS queries="+results.length+" search_seconds="+searchSeconds);
	}

	private static void verify(HashMap<String,String> o){
		HbmTailAssign.verifyReference(o);
		HbmProfilePilot.requireHash(HmmComparisonData.required(o, "in"), HmmComparisonData.required(o, "insha80"));
	}
	/** Shared production/test seam: the reduced arm is the exact dense-rank F4 prefix. */
	public static void selectCandidates(float[] scores,int[] out){
		require(scores!=null && scores.length>=SHORTLIST && out!=null && (out.length==100 || out.length==SHORTLIST),
			"Both bounded shortlist arms require the same200-or-more-family score roster");
		HbmCompetitiveAssign.top(scores,out);
	}

	/** Native FASTA reader; the bounded panel fits in memory and keeps output order deterministic. */
	private static ProteinSequence[] readQueries(String path,int expected){
		assert(expected>0 && expected<=MAX_QUERIES) : "Preloading is safe only for the explicitly bounded pilot";
		final ArrayList<ProteinSequence> queries=new ArrayList<ProteinSequence>(expected);
		final HashSet<String> ids=new HashSet<String>();
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, false), null, true, expected+1L);
		input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(queries.size()<expected, "Input exceeds declared panel count; no implicit truncation is allowed");
					require(read.mate==null && read.bases!=null, "Pilot expects unpaired amino-acid records");
					final ProteinSequence query=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
					require(ids.add(query.id) && query.length()==read.bases.length, "Duplicate ID or stripped terminal stop: "+query.id);
					validateTail(query.enc); queries.add(query);
				}
			}
		}finally{input.close();}
		require(!input.errorState() && queries.size()==expected, "Pilot FASTA read/count failed: expected="+expected+" observed="+queries.size());
		return queries.toArray(new ProteinSequence[queries.size()]);
	}

	/** Immutable model roster; each thread owns its candidate marks and consensus cache. */
	public static final class Engine {
		public Engine(String[] ids_,byte[][] refs_,HbmPositionModel[] profiles_){
			require(ids_!=null && refs_!=null && profiles_!=null && ids_.length>0 && ids_.length==refs_.length
				&& ids_.length==profiles_.length, "Pilot model roster dimensions must agree");
			ids=ids_.clone(); refs=new byte[refs_.length][]; profiles=profiles_.clone();
			final HashSet<String> seen=new HashSet<String>();
			for(int i=0; i<ids.length; i++){
				require(ids[i]!=null && ids[i].matches("[!-~]+") && seen.add(ids[i]), "Family IDs must be unique printable ASCII tokens");
				validateTail(refs_[i]); refs[i]=refs_[i].clone();
				require(profiles[i]!=null && profiles[i].length()==refs[i].length, "HBM dimensions must match its consensus at rank="+i);
				profiles[i].requireConsensus(refs[i]);
			}
		}

		/** HBM score zero participates in ties, never an extra membership rejection gate. */
		public Result select(byte[] query,int[] candidates,Scratch work){
			final long hbmStart=work!=null && work.profile ? System.nanoTime() : 0;
			validateTail(query);
			require(candidates!=null && candidates.length>0 && candidates.length<=SHORTLIST && work!=null
				&& work.seen.length==ids.length, "Pilot shortlist/scratch dimensions differ");
			if(work.patience>0){
				require(candidates.length==SHORTLIST, "Patience20 requires the full200-candidate shortlist");
				return selectWithPatience(query, candidates, work, hbmStart);
			}
			if(++work.generation==0){Arrays.fill(work.seen, 0); work.generation=1;}
			int winner=-1, best=-1, ties=0, nonpositive=0;
			for(int rank : candidates){
				require(rank>=0 && rank<ids.length && work.seen[rank]!=work.generation, "Duplicate or out-of-range candidate: "+rank);
				work.seen[rank]=work.generation;
				final int score;
				if(work.scoreOnly){score=profiles[rank].scoreTail(query, 200);}
				else{
					final HbmPositionModel.TailResult h=profiles[rank].alignTail(query, 200, false);
					score=h==null ? 0 : h.score;
				}
				if(score==0){nonpositive++;}
				if(score>best){best=score; winner=rank; ties=1;}
				else if(score==best){ties++; if(ids[rank].compareTo(ids[winner])<0){winner=rank;}}
			}
			assert(winner>=0 && ties>=1) : "A nonempty HBM shortlist always chooses one winner, even at score zero";
			final long blosumStart=work.profile ? System.nanoTime() : 0;
			if(work.profile){work.hbmNanos+=blosumStart-hbmStart;}
			// At most one consensus scorer per encountered winner per worker; no per-candidate BLOSUM work.
			if(work.owner!=this){Arrays.fill(work.consensus, null); work.owner=this;}
			HbmPositionModel scorer=work.consensus[winner];
			if(scorer==null){scorer=work.consensus[winner]=HbmPositionModel.fromConsensus(refs[winner]);}
			final HbmPositionModel.TailResult alignment=scorer.alignTail(query, 200, true);
			final Result result=new Result(winner, best, ties, candidates.length, nonpositive, measure(query, refs[winner], alignment));
			if(work.profile){work.blosumNanos+=System.nanoTime()-blosumStart;}
			return result;
		}

		/** Independent opt-in: never fall back from a better but ineligible HBM winner. */
		private Result selectWithPatience(byte[] query,int[] candidates,Scratch work,long hbmStart){
			assert(work.patience==20 && candidates.length==SHORTLIST) : "This pilot isolates patience20 from shortlist reduction";
			if(++work.generation==0){Arrays.fill(work.seen, 0); work.generation=1;}
			// Validate even candidates beyond a possible early stop; malformed input must not hide there.
			for(int rank : candidates){
				require(rank>=0 && rank<ids.length && work.seen[rank]!=work.generation, "Duplicate or out-of-range candidate: "+rank);
				work.seen[rank]=work.generation;
			}
			int winner=-1, best=-1, ties=0, nonpositive=0, scored=0, stale=0, checkedWinner=-1, blosumScored=0;
			Metrics metrics=null;
			long segmentStart=hbmStart;
			for(int rank : candidates){
				final int score;
				if(work.scoreOnly){score=profiles[rank].scoreTail(query, 200);}
				else{
					final HbmPositionModel.TailResult h=profiles[rank].alignTail(query, 200, false);
					score=h==null ? 0 : h.score;
				}
				scored++;
				if(score==0){nonpositive++;}
				if(score>best){best=score; winner=rank; ties=1; stale=0;}
				else{
					stale++;//Equal scores never reset patience, even when the ASCII winner changes.
					if(score==best){ties++; if(ids[rank].compareTo(ids[winner])<0){winner=rank;}}
				}
				if(stale>=work.patience){
					if(checkedWinner!=winner){
						final long start=work.profile ? System.nanoTime() : 0;
						if(work.profile){work.hbmNanos+=start-segmentStart;}
						metrics=consensusMetrics(query, winner, work); checkedWinner=winner; blosumScored++;
						if(work.profile){segmentStart=System.nanoTime(); work.blosumNanos+=segmentStart-start;}
					}
					if(metrics.accepted()){break;}
				}
			}
			assert(winner>=0 && ties>=1 && scored>=work.patience+1) : "First candidate establishes best;20 subsequent nonimprovements precede early acceptance";
			final long finalStart=work.profile ? System.nanoTime() : 0;
			if(work.profile){work.hbmNanos+=finalStart-segmentStart;}
			if(checkedWinner!=winner){
				metrics=consensusMetrics(query, winner, work); blosumScored++;
				if(work.profile){work.blosumNanos+=System.nanoTime()-finalStart;}
			}
			return new Result(winner, best, ties, scored, nonpositive, metrics, blosumScored);
		}
		private Metrics consensusMetrics(byte[] query,int winner,Scratch work){
			assert(winner>=0 && winner<refs.length) : "Only the current HBM winner receives a consensus eligibility check";
			if(work.owner!=this){Arrays.fill(work.consensus, null); work.owner=this;}
			HbmPositionModel scorer=work.consensus[winner];
			if(scorer==null){scorer=work.consensus[winner]=HbmPositionModel.fromConsensus(refs[winner]);}
			return measure(query, refs[winner], scorer.alignTail(query, 200, true));
		}
		private final String[] ids;
		private final byte[][] refs;
		private final HbmPositionModel[] profiles;
	}

	public static final class Scratch {
		public Scratch(int families){
			this(families, false);
		}
		/** Opt-in wall timers; no per-candidate clock calls or shared counters. */
		public Scratch(int families,boolean profile_){
			this(families, profile_, false);
		}
		public Scratch(int families,boolean profile_,boolean scoreOnly_){
			this(families, profile_, scoreOnly_, 0);
		}
		public Scratch(int families,boolean profile_,boolean scoreOnly_,int patience_){
			require(families>0, "Candidate generation marks need a nonempty roster");
			require(patience_==0 || patience_==20, "Only default0 or independent patience20 is supported");
			profile=profile_;
			scoreOnly=scoreOnly_;
			patience=patience_;
			seen=new int[families]; consensus=new HbmPositionModel[families];
		}
		public long hbmNanos(){return hbmNanos;}
		public long blosumNanos(){return blosumNanos;}
		private final boolean profile, scoreOnly;
		private final int patience;
		private long hbmNanos, blosumNanos;
		private final int[] seen;
		private int generation;
		private final HbmPositionModel[] consensus;
		private Engine owner;
	}

	/** Measures the BLOSUM traceback, not the HBM traceback. Coordinates are zero-based inclusive. */
	public static Metrics measure(byte[] query,byte[] ref,HbmPositionModel.TailResult a){
		validateTail(query); validateTail(ref);
		// Same ordinary diagonal sum as ProteinSearcher.selfScore(ProteinSequence), without a decoded copy.
		int self=0; for(byte aa : query){self+=Blosum62.score(aa, aa);}
		require(self>0 && Float.isFinite(self), "A/Sq requires a positive full-query ordinary BLOSUM self score");
		if(a==null){return new Metrics(null, self, 0, 0, 0, 0, 0, 15);}
		require(a.path!=null && a.path.length>0 && a.score>0 && a.queryStart>=0 && a.start>=0
			&& a.queryStart<query.length && a.start<ref.length && a.queryEnd==query.length-1 && a.end==ref.length-1,
			"Winner trace must consume both actual right endpoints; terminal gaps are allowed");
		int q=a.queryStart, r=a.start, matches=0, paired=0, score=0;
		for(byte op : a.path){
			if(op=='m'){
				require(q<query.length && r<ref.length, "Paired trace column exceeds a sequence endpoint");
				if(query[q]==ref[r]){matches++;} paired++;
				score+=(int)(HbmPositionModel.SCALE*Blosum62.score(query[q++], ref[r++]));
			}else if(op=='I'){require(q<query.length, "Insertion exceeds query endpoint"); q++; score-=HbmPositionModel.GAP;}
			else if(op=='D'){require(r<ref.length, "Deletion exceeds consensus endpoint"); r++; score-=HbmPositionModel.GAP;}
			else{throw new IllegalArgumentException("Unknown BLOSUM trace operation: "+op);}
		}
		require(q==query.length && r==ref.length && score==a.score, "Trace consumption/score differs from bottom-right BLOSUM alignment");
		final int qSpan=a.queryEnd-a.queryStart+1, rSpan=a.end-a.start+1;
		final int mask=failureMask(qSpan, rSpan, matches, a.path.length, a.score, self);
		return new Metrics(a, self, qSpan, rSpan, matches, paired, a.path.length, mask);
	}

	/** Rational cutoffs avoid rounding near30% identity and0.118 normalized score. */
	public static int failureMask(int querySpan,int refSpan,int matches,int columns,int score64,float self){
		require(querySpan>=0 && refSpan>=0 && matches>=0 && columns>=matches && score64>=0 && self>0 && Float.isFinite(self),
			"Gate metrics must be finite, nonnegative and dimensionally consistent");
		int mask=0;
		if(querySpan<100){mask|=1;}
		if(refSpan<100){mask|=2;}
		if(columns==0 || 10L*matches<3L*columns){mask|=4;}
		if(1000.0*score64<118.0*HbmPositionModel.SCALE*self){mask|=8;}
		return mask;
	}

	public static final class Metrics {
		Metrics(HbmPositionModel.TailResult alignment_,float self_,int querySpan_,int refSpan_,int matches_,int paired_,int columns_,int failures_){
			alignment=alignment_; self=self_; querySpan=querySpan_; refSpan=refSpan_; matches=matches_; paired=paired_; columns=columns_; failures=failures_;
		}
		public final HbmPositionModel.TailResult alignment;
		public final float self;
		public final int querySpan, refSpan, matches, paired, columns, failures;
		public boolean accepted(){return failures==0;}
	}
	public static final class Result {
		Result(int winner_,int hbmScore64_,int ties_,int scored_,int nonpositive_,Metrics metrics_){
			this(winner_, hbmScore64_, ties_, scored_, nonpositive_, metrics_, 1);
		}
		Result(int winner_,int hbmScore64_,int ties_,int scored_,int nonpositive_,Metrics metrics_,int blosumScored_){
			winner=winner_; hbmScore64=hbmScore64_; ties=ties_; scored=scored_; nonpositive=nonpositive_; metrics=metrics_;
			blosumScored=blosumScored_;
		}
		public final int winner, hbmScore64, ties, scored, nonpositive, blosumScored;
		public final Metrics metrics;
	}

	private static final class ProcessThread extends Thread {
		ProcessThread(int id_,int stride_,ProteinSequence[] queries_,Result[] results_,Engine engine_,FamilyShortlistSidecar sidecar_,AtomicReference<Throwable> failure_,boolean profile_,boolean scoreOnly_,int shortlistSize_,int patience_){
			id=id_; stride=stride_; queries=queries_; results=results_; engine=engine_; sidecar=sidecar_; failure=failure_;
			profile=profile_;
			scoreOnly=scoreOnly_; shortlistSize=shortlistSize_; patience=patience_;
			setDaemon(true);
		}
		@Override public void run(){
			final long start=profile ? System.nanoTime() : 0;
			try{
				assert(id>=0 && id<stride) : "Disjoint worker strides preserve one result per input index";
				final ProteinSearcher.ShortlistScratch shortlist=new ProteinSearcher.ShortlistScratch(sidecar.dims, engine.ids.length, shortlistSize);
				work=new Scratch(engine.ids.length, profile, scoreOnly, patience);
				for(int i=id; i<queries.length && failure.get()==null; i+=stride){
					final long f4Start=profile ? System.nanoTime() : 0;
					ProteinSearcher.scoreFamiliesF4(queries[i], sidecar, shortlist);
					final long topStart=profile ? System.nanoTime() : 0;
					if(profile){f4Nanos+=topStart-f4Start;}
					selectCandidates(shortlist.f4, shortlist.topIdx);
					if(profile){topNanos+=System.nanoTime()-topStart;}
					results[i]=engine.select(queries[i].enc, shortlist.topIdx, work);
					if(++count%128==0){System.err.println("HBM_WINNER_PILOT_PROGRESS worker="+id+" queries="+count);}
				}
			}catch(Throwable t){failure.compareAndSet(null, t);}
			finally{if(profile){elapsedNanos=System.nanoTime()-start;}}
		}
		final boolean profile, scoreOnly;
		final int shortlistSize, patience;
		Scratch work;
		long f4Nanos, topNanos, elapsedNanos;
		int count;
		final int id, stride;
		final ProteinSequence[] queries;
		final Result[] results;
		final Engine engine;
		final FamilyShortlistSidecar sidecar;
		final AtomicReference<Throwable> failure;
	}

	/** Per-worker elapsed phases overlap across threads; their sum is NOT process wall or CPU time. */
	private static void writeProfile(Path out,ArrayList<ProcessThread> workers,double searchSeconds){
		assert(!workers.isEmpty()) : "A successful nonempty pilot must have completed workers";
		final ByteStreamWriter writer=HmmComparisonData.writer(out.resolve("profile.tsv").toString());
		final ByteBuilder row=new ByteBuilder();
		writer.println("worker\tqueries\tf4_seconds\ttop_seconds\thbm_seconds\tblosum_metrics_seconds\tother_seconds\tworker_seconds\tsearch_wall_seconds");
		try{
			for(ProcessThread w : workers){
				final long other=w.elapsedNanos-w.f4Nanos-w.topNanos-w.work.hbmNanos-w.work.blosumNanos;
				require(other>=0 && w.count>0, "Disjoint profile phases must fit inside worker elapsed time");
				writer.println(row.clear().append(w.id).tab().append(w.count).tab().append(w.f4Nanos*1e-9, 9)
					.tab().append(w.topNanos*1e-9, 9).tab().append(w.work.hbmNanos*1e-9, 9)
					.tab().append(w.work.blosumNanos*1e-9, 9).tab().append(other*1e-9, 9)
					.tab().append(w.elapsedNanos*1e-9, 9).tab().append(searchSeconds, 9));
			}
		}finally{HmmComparisonData.close(writer);}
	}

	private static void writeResults(Path out,ProteinSequence[] queries,Result[] results,String[] ids,HashMap<String,String> options,
			int threads,int shortlist,double setupSeconds,double searchSeconds){
		assert(queries.length==results.length) : "Publication must conserve one row per original query";
		final ByteStreamWriter calls=HmmComparisonData.writer(out.resolve("assignments.tsv").toString());
		final ByteBuilder row=new ByteBuilder(); long accepted=0, nonpositive=0, tied=0, hbmScored=0, blosumScored=0;
		calls.println(ASSIGNMENT_HEADER);
		try{
			for(int i=0; i<results.length; i++){
				final Result result=results[i]; final Metrics m=result.metrics;
				if(m.accepted()){accepted++;} if(result.ties>1){tied++;} nonpositive+=result.nonpositive;
				hbmScored+=result.scored; blosumScored+=result.blosumScored;
				calls.println(formatRow(row, queries[i].id, result, ids));
			}
		}finally{HmmComparisonData.close(calls);}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("queries\tassigned\trejected\tfamilies\tshortlist\thbm_scored\thbm_nonpositive\ttied_winners\tblosum_scored\tworkers\tsetup_seconds\tsearch_seconds");
		summary.println(row.clear().append(results.length).tab().append(accepted).tab().append(results.length-accepted).tab().append(ids.length)
			.tab().append(shortlist).tab().append(hbmScored).tab().append(nonpositive).tab().append(tied)
			.tab().append(blosumScored).tab().append(threads).tab().append(setupSeconds, 9).tab().append(searchSeconds, 9));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		for(String key : INPUTS){recipe.println(key+"\t"+options.get(key)); recipe.println(key+"_sha80\t"+options.get(key+"sha80"));}
		// Leave legacy/default recipe bytes unchanged; the opt-in kernel is explicit provenance.
		if(options.getOrDefault("hbmscoreonly", "f").equals("t")){recipe.println("hbm_kernel\tscoreTail_rolling_row_exact_endpoint");}
		if(Integer.parseInt(options.getOrDefault("patience", "0"))==20){
			recipe.println("early_stop\t20_nonimproving_candidates_current_winner_passes_gates_ties_do_not_reset");
		}
		recipe.println("policy\thbm_winner_then_blosum_no_fallback\nshortlist\tF4_top"+shortlist+"_dense_rank_ties\nhbm_ranking\tscore64_then_ascii_family_id_including_zero\nprofile\tlogodds_beta0.01_clip-4_11_gap4\nblosum\tBLOSUM62_linear_gap4_bottom_right_zero_floor_left_trim\nminimum_actual_span_each\t100\nidentity\tmatches_over_all_path_columns_ge_0.30\nnormalized_score\tA_over_full_query_ordinary_self_score_ge_0.118\nterminal_gaps\tallowed_with_penalty\nfailure_mask\t1_query_span_2_consensus_span_4_identity_8_normalized_score\nfull_run\tNOT_RELEASED");
		HmmComparisonData.close(recipe);
	}
	/** Shared byte-for-byte pilot/streaming representation, including rejected winners. */
	public static ByteBuilder formatRow(ByteBuilder row,String query,Result result,String[] ids){
		assert(result.winner>=0 && result.winner<ids.length) : "Engine winner must index the same roster used in output";
		final Metrics m=result.metrics; final HbmPositionModel.TailResult a=m.alignment;
		return row.clear().append(query).tab().append(m.accepted() ? "ASSIGNED" : "REJECTED").tab().append(result.winner)
			.tab().append(ids[result.winner]).tab().append(result.hbmScore64).tab().append(result.ties).tab().append(result.scored).tab().append(result.nonpositive)
			.tab().append(a==null ? 0 : a.score).tab().append(m.self, 0).tab().append(a==null ? 0 : a.score/(HbmPositionModel.SCALE*(double)m.self), 12)
			.tab().append(a==null ? -1 : a.queryStart).tab().append(a==null ? -1 : a.queryEnd).tab().append(a==null ? -1 : a.start).tab().append(a==null ? -1 : a.end)
			.tab().append(m.querySpan).tab().append(m.refSpan).tab().append(m.paired).tab().append(m.matches).tab().append(m.columns).tab().append(m.failures);
	}
	public static final String ASSIGNMENT_HEADER="query\tstatus\twinner_rank\tfamily_id\thbm_score64\thbm_ties\thbm_scored\thbm_nonpositive\tblosum_score64\tquery_self_score\tA_over_Sq\tquery_start\tquery_end\tconsensus_start\tconsensus_end\tquery_span\tconsensus_span\tpaired\tmatches\tcolumns\tfailure_mask";
	public static void validateTail(byte[] sequence){
		require(sequence!=null && sequence.length>=100 && sequence.length<=200, "Pilot tails must contain100..200 residues");
		for(byte aa : sequence){require(aa>=0 && aa<20, "Pilot tail contains ambiguous or stop residue: "+aa);}
	}
	private static double seconds(long start){return (System.nanoTime()-start)*1e-9;}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String[] INPUTS={"in", "models", "consensus", "sidecar", "kmers", "background", "hbm", "provenance", "packrecipe"};
	private static final int MAX_QUERIES=4096, SHORTLIST=200;
	private HbmWinnerCalibration(){}
}
