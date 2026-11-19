package prot;

import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/** One bounded worker pool scores complete labeled panels through a loaded classifier.
 * Worker counts and optional calibration maxima are merged only after join.
 * @author Brian, Keqing
 */
public final class FamilyPanelRunner implements Accumulator<FamilyPanelRunner.ProcessThread>{

	/** Resident evaluation with optional hit telemetry and frozen classifier gates.
	 * Copies the panel roster; descriptors are immutable. Input pins remain the caller's
	 * responsibility, as does protein-reader configuration (Shared.AMINO_IN=true),
	 * as for existing panel evaluation. No calibration updates occur. */
	public FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, PrintStream outstream_, String label_, int minQueryAa_,
			int maxQueryAa_, boolean hitTelemetry_){
		this(classifier_,panels_==null ? null : panels_.clone(),threads_,inputThreads_,null,outstream_,label_,null,
			minQueryAa_,maxQueryAa_,false,false,hitTelemetry_);
	}

	/** Resident GC bindings are shared read-only; their roster is copied and verified each round.
	 * Explicit bindings may also accompany an ordinary classifier to measure baseline coverage. */
	public FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, PrintStream outstream_, String label_, int minQueryAa_,
			int maxQueryAa_, boolean hitTelemetry_, FamilyQueryGc[] queryGc_){
		this(classifier_, panels_, threads_, inputThreads_, null, outstream_, label_, null,
			minQueryAa_, maxQueryAa_, false, false, hitTelemetry_, 1, queryGc_);
	}

	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_,Panel[] panels_,int threads_,int inputThreads_,
			FamilyCalibrationUpdates updates_,PrintStream outstream_,String label_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, null);
	}

	/** Optional length boundaries add compact counters, never additional sequence readers or scorers. */
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_, int threads_, int inputThreads_,
			FamilyCalibrationUpdates updates_, PrintStream outstream_, String label_, int[] lengthEdges_){
		this(classifier_,panels_,threads_,inputThreads_,updates_,outstream_,label_,lengthEdges_,0,Integer.MAX_VALUE);
	}
	/** Query bounds validate the classifier's normalized length; they never select or crop input records. */
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, lengthEdges_,
			minQueryAa_, maxQueryAa_, false);
	}
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_, boolean admissionReport_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, lengthEdges_,
			minQueryAa_, maxQueryAa_, admissionReport_, false);
	}
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_, boolean admissionReport_,
			boolean stageCpu_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, lengthEdges_,
			minQueryAa_, maxQueryAa_, admissionReport_, stageCpu_, false);
	}
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_, boolean admissionReport_,
			boolean stageCpu_, boolean hitTelemetry_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, lengthEdges_,
			minQueryAa_, maxQueryAa_, admissionReport_, stageCpu_, hitTelemetry_, 1);
	}
	/** Calibrate against an explicit scaled view while retaining the immutable raw
	 * parent for update publication. Existing evaluation-only views are untouched
	 * by the default constructor path. */
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_, boolean admissionReport_,
			boolean stageCpu_, boolean hitTelemetry_, double updateMarginScale_){
		this(classifier_, panels_, threads_, inputThreads_, updates_, outstream_, label_, lengthEdges_,
			minQueryAa_, maxQueryAa_, admissionReport_, stageCpu_, hitTelemetry_, updateMarginScale_, null);
	}
	/** The same immutable query bindings survive calibration rounds and margin views. */
	FamilyPanelRunner(FamilyTailCascade.LoadedClassifier classifier_, Panel[] panels_,
			int threads_, int inputThreads_, FamilyCalibrationUpdates updates_, PrintStream outstream_,
			String label_, int[] lengthEdges_, int minQueryAa_, int maxQueryAa_, boolean admissionReport_,
			boolean stageCpu_, boolean hitTelemetry_, double updateMarginScale_, FamilyQueryGc[] queryGc_){
		require(Double.isFinite(updateMarginScale_) && updateMarginScale_>0 && (updateMarginScale_==1 || updates_!=null),
			"A calibration margin requires a positive finite scale and raw-parent update accumulator");
		updateMarginScale=updateMarginScale_;
		require(minQueryAa_>=0 && maxQueryAa_>=minQueryAa_,"Panel query bounds must be nonnegative and ordered");
		minQueryAa=minQueryAa_;maxQueryAa=maxQueryAa_;
		require(classifier_!=null && classifier_.familyCount()>0 && panels_!=null && panels_.length>0,
			"Panel execution requires a loaded live classifier and nonempty panel roster");
		require(threads_>0 && inputThreads_>0 && outstream_!=null && label_!=null,
			"Panel execution needs positive worker counts and a progress destination");
		classifier=updateMarginScale==1 ? classifier_ : classifier_.withEvaluationMarginScale(updateMarginScale);
		panels=panels_.clone();threads=threads_;inputThreads=inputThreads_;
		updates=updates_;outstream=outstream_;label=label_;
		require(!hitTelemetry_ || updates==null,"Hit telemetry evaluates frozen gates without calibration updates");
		for(Panel panel:panels){require(panel!=null,"Panel descriptors must be nonnull");}
		queryGc=queryGc_==null ? null : queryGc_.clone();
		require(!classifier.usesFamilyGcPenalty() || queryGc!=null,
			"GC ranking requires explicit measured/neutral bindings for every panel");
		if(queryGc!=null){
			require(queryGc.length==panels.length,"Query GC and scoring panel rosters must have equal length");
			for(int p=0; p<panels.length; p++){
				require(queryGc[p]!=null,"Every scoring panel needs its explicit measured/neutral GC binding");
				queryGc[p].requirePanel(panels[p].path, panels[p].expectedQueries);
			}
		}
		require(updates==null || updates.boundTo(classifier_.config()),
			"Calibration maxima must belong to the exact classifier roster and acceptance policy being scored");
		counts=new PanelCounts[panels.length];
		for(int i=0;i<counts.length;i++){counts[i]=new PanelCounts(classifier.familyCount());}
		lengthEdges=lengthEdges_==null ? null : lengthEdges_.clone();
		if(lengthEdges!=null){
			require(lengthEdges.length>0, "Length reporting requires at least one boundary");
			for(int i=0; i<lengthEdges.length; i++){
				require(lengthEdges[i]>0 && (i==0 || lengthEdges[i]>lengthEdges[i-1]),
					"Length boundaries must increase strictly from a positive first boundary");
			}
		}
		lengthCounts=newLengthCounts();
		admissionReport=admissionReport_; stageCpu=stageCpu_;
		admissionCounts=newAdmissionCounts();
		hitTelemetry=hitTelemetry_;
		hitHistograms=newHitHistograms();
	}

	/** The caller pins inputs and publishes output; this method only scores and reduces. */
	public void process(){
		require(!started,"A panel runner represents one evaluation; use a fresh runner for each round");started=true;
		verifyQueryGc();
		final Streamer[] inputs=new Streamer[panels.length];
		for(int i=0;i<inputs.length;i++){
			final FileFormat format=FileFormat.testInput(panels[i].path.toString(),FileFormat.FASTA,null,true,true);
			require(format.fasta() || format.fastq(),"Protein panels must be FASTA or FASTQ: "+panels[i].path);
			inputs[i]=StreamerFactory.makeStreamer(format,null,true,-1,false,true,inputThreads);
		}
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0;i<threads;i++){workers.add(new ProcessThread(i,inputs));}
		try{
			for(Streamer input:inputs){input.start();}
			if(!ThreadWaiter.startAndWait(workers,this)){errorState=true;}
		}catch(Throwable problem){
			problem.printStackTrace(outstream);
			KillSwitch.kill("Family panel execution could not start or join its scoring workers");
		}
		for(Streamer input:inputs){errorState|=ReadWrite.closeStream(input);}
		if(failure.get()!=null){
			throw new IllegalStateException("Family panel execution failed; no completed report is available",
				failure.get());
		}
		require(!errorState,"Family panel input or worker failed");
		long total=0;
		for(int i=0;i<counts.length;i++){
			counts[i].validate(panels[i].name);total=Math.addExact(total,counts[i].queries);
			if(queryGc!=null){
				final PanelCounts c=counts[i];
				require(Math.addExact(c.measuredGcQueries, c.neutralGcQueries)==c.queries &&
					Math.addExact(c.measuredGcAccepts, c.neutralGcAccepts)==c.accepted &&
					c.measuredGcAccepts<=c.measuredGcQueries && c.neutralGcAccepts<=c.neutralGcQueries,
					"GC coverage must partition every scored query and accepted winner");
			}
			require(panels[i].expectedQueries<0 || counts[i].queries==panels[i].expectedQueries,
				"Panel "+panels[i].name+" count differs from its declared population: "+
				counts[i].queries+" != "+panels[i].expectedQueries);
		}
		require(total==progress.get(),"Worker totals disagree with processed input batches");
		validateLengths();
		validateHitHistograms();
		if(acceptingRanks!=null){
			for(int p=0;p<panels.length;p++){acceptingRanks[p].validate(counts[p].accepted);}
		}
		if(admissionCounts!=null){
			for(int p=0;p<panels.length;p++){
				require(admissionCounts[p][0][0]==counts[p].queries,
					"Admission callbacks must cover every panel query exactly once");
				require(admissionCounts[p][0][5]==counts[p].alignmentAttempts,"Admission and worker alignment attempts must agree");
				if(lengthEdges!=null){
					long queries=0,pool=0,dp=0;
					for(int b=1;b<admissionCounts[p].length;b++){
						queries=Math.addExact(queries,admissionCounts[p][b][0]);pool=Math.addExact(pool,admissionCounts[p][b][1]);
						dp=Math.addExact(dp,admissionCounts[p][b][5]);
					}
					require(queries==counts[p].queries && pool==admissionCounts[p][0][1] && dp==admissionCounts[p][0][5],
						"Admission length bins must partition all queries, candidate visits and DP attempts");
				}
			}
		}
		if(stageCpu){
			stageCpuCounts=new FamilyTailCascade.StageCpuCounts();
			for(ProcessThread worker:workers){stageCpuCounts.add(worker.stageCpuCountsT);}
		}
		verifyQueryGc();complete=true;
	}
	/** Metadata verification is outside workers; a changed input invalidates the entire round. */
	private void verifyQueryGc(){
		if(queryGc==null){return;}
		try{FamilyQueryGc.verifyAll(queryGc);}
		catch(java.io.IOException problem){
			throw new IllegalStateException("Could not verify retained query GC inputs", problem);
		}
	}
	/** Coverage requires explicit GC bindings, including for a coefficient-zero baseline. */
	public long gcQueries(int panel, boolean measured){
		final PanelCounts c=completedPanel(panel);require(queryGc!=null,"GC coverage was not requested for this evaluation");
		return measured ? c.measuredGcQueries : c.neutralGcQueries;
	}
	public long gcAccepts(int panel, boolean measured){
		final PanelCounts c=completedPanel(panel);require(queryGc!=null,"GC coverage was not requested for this evaluation");
		return measured ? c.measuredGcAccepts : c.neutralGcAccepts;
	}
	/** Small final coverage report; no per-query intermediate is emitted. */
	public void writeGcCoverage(Path path){
		require(complete && queryGc!=null,"GC coverage requires a completed explicitly bound evaluation");
		final ByteBuilder out=new ByteBuilder("panel\tclass\tgc_coverage\tqueries\taccepted\n");
		for(int p=0; p<panels.length; p++){
			for(int measured=1; measured>=0; measured--){
				out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative").tab()
					.append(measured==1 ? "measured" : "neutral").tab().append(gcQueries(p, measured==1))
					.tab().append(gcAccepts(p, measured==1)).nl();
			}
		}
		write(path, out);
	}

	/** Detached seven-bin winner rows, indexed by live family ordinal; bin6 means6+. */
	public long[] copyWinnerBins(final int panel){return completedHistogram(panel).copyWinnerBins();}
	/** Detached uncensored raw-hit rows for a negative panel, including derived zero bins. */
	public long[] copyRawBins(final int panel){return completedHistogram(panel).copyRawBins();}
	/** Actual completed valid queries in this panel, never multiplied by the family count. */
	public long completedQueries(final int panel){return completedPanel(panel).queries;}
	public long acceptedQueries(final int panel){return completedPanel(panel).accepted;}
	/** Candidate alignment-stage visits, including unavailable calibration; requires admission capture. */
	public long dpAttempts(final int panel){return completedAdmission(panel)[5];}
	/** Actual alignment invocations, including calls that returned no usable measurement. */
	public long alignmentCalls(final int panel){
		final PanelCounts c=completedPanel(panel);
		return c.alignmentAttempts-c.alignmentUnavailable;
	}
	public long alignmentUnavailable(final int panel){return completedPanel(panel).alignmentUnavailable;}
	public long alignmentFailures(final int panel){return completedPanel(panel).alignmentFailures;}
	/** Sum of uncapped admitted family candidates across queries, not accepted-query count. */
	public long admittedCandidates(final int panel){return completedAdmission(panel)[1];}
	private long[] completedAdmission(final int panel){
		completedPanel(panel);
		require(admissionCounts!=null,"Candidate and DP work require admission capture enabled during this evaluation");
		return admissionCounts[panel][0];
	}

	/** Primitive payload of ALL retained histogram arrays during worker reduction.
	 * Each positive panel needs64F bytes (7F longs +2F ints); each negative needs120F
	 * (14F longs +2F ints). Multiply by workers+one aggregate, not by one active panel.
	 * Excludes preexisting counters/index/scratch, object/array headers and reference
	 * arrays. There are (workers+1)*panels collector objects and roster reference slots,
	 * with three primitive arrays per positive collector and four per negative one.
	 * Completion validation temporarily copies one7F-long row table at a time;
	 * each caller-held detached export adds another56F bytes plus its array header. */
	public long histogramPrimitiveBytes(){
		if(!hitTelemetry){return 0;}
		long perWorker=0;
		for(Panel panel:panels){
			perWorker=Math.addExact(perWorker,Math.multiplyExact((long)classifier.familyCount(),panel.positive ? 64L : 120L));
		}
		return Math.multiplyExact(perWorker,(long)threads+1);
	}
	private FamilyHitHistogram completedHistogram(final int panel){
		completedPanel(panel);
		require(hitHistograms!=null,"Hit telemetry requires capture enabled");
		return hitHistograms[panel];
	}
	private PanelCounts completedPanel(final int panel){
		require(complete,"Panel counts require a completed evaluation");
		require(panel>=0 && panel<panels.length,"Panel index is outside the evaluated roster");
		return counts[panel];
	}
	private FamilyHitHistogram[] newHitHistograms(){
		if(!hitTelemetry){return null;}
		final FamilyHitHistogram[] histograms=new FamilyHitHistogram[panels.length];
		for(int p=0;p<panels.length;p++){histograms[p]=new FamilyHitHistogram(classifier.familyCount(),!panels[p].positive);}
		return histograms;
	}
	private void validateHitHistograms(){
		if(hitHistograms==null){return;}
		for(int p=0;p<panels.length;p++){
			final FamilyHitHistogram histogram=hitHistograms[p];
			require(histogram.completedQueries()==counts[p].queries && histogram.acceptedQueries()==counts[p].accepted,
				"Hit telemetry denominators disagree with terminal observations");
			validateHitRows(histogram.copyWinnerBins(),p,false);
			if(!panels[p].positive){validateHitRows(histogram.copyRawBins(),p,true);}
		}
	}
	private void validateHitRows(final long[] bins,final int panel,final boolean raw){
		long winners=0;
		for(int f=0;f<classifier.familyCount();f++){
			long total=0;final int offset=f*FamilyHitHistogram.BINS;
			for(int b=0;b<FamilyHitHistogram.BINS;b++){
				require(bins[offset+b]>=0,"Hit telemetry cannot contain negative counts");
				total=Math.addExact(total,bins[offset+b]);
			}
			require(total==(raw ? counts[panel].queries : counts[panel].winners[f]),
				"Hit row does not conserve its raw-query or terminal-winner denominator");
			if(!raw){
				require(bins[offset]==0,"A terminal winner must have a positive raw hit count");
				winners=Math.addExact(winners,total);
			}
		}
		if(!raw){require(winners==counts[panel].accepted,"Winner bins do not conserve accepted queries");}
	}

	/** Untrimmed admitted-pool counts share the existing scoring pass and optional length bins. */
	void writeAdmission(Path path){
		require(complete && admissionCounts!=null,
			"Admission reporting requires a completed evaluation with telemetry enabled");
		final ByteBuilder out=new ByteBuilder(
			"panel\tclass\tscope\tlength_min_inclusive\tlength_max_exclusive\tqueries"+
			"\tadmitted_sum\tadmitted_per_query\tadmitted_max\tfloor_min\tfloor_max"+
			"\tdp_attempts_sum\tdp_attempts_per_query\n");
		for(int p=0;p<panels.length;p++){
			for(int b=0;b<admissionCounts[p].length;b++){
				final long[] row=admissionCounts[p][b];
				out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative")
					.tab().append(b==0 ? "all" : "length_bin").tab().append(b<=1 ? 0 : lengthEdges[b-2]).tab();
				if(b==0 || b==lengthEdges.length+1){out.append("NA");}else{out.append(lengthEdges[b-1]);}
				out.tab().append(row[0]).tab().append(row[1]).tab();
				if(row[0]==0){out.append("NA\tNA\tNA\tNA\t0\tNA\n");}
				else{out.appendSlow(row[1]/(double)row[0]).tab().append(row[2]).tab().append(row[3]).tab().append(row[4])
					.tab().append(row[5]).tab().appendSlow(row[5]/(double)row[0]).nl();}
			}
		}
		write(path,out);
	}
	/** Optional terminal-rank telemetry reuses ordinary scoring, never full candidate tracing. */
	public void enableAcceptingRankStats(){
		require(!started && updates==null && acceptingRanks==null,
			"Accepting ranks require one fresh frozen evaluation, never calibration updates or repeated enablement");
		acceptingRanks=newAcceptingRankHistograms();
	}
	private AcceptingRankHistogram[] newAcceptingRankHistograms(){
		final int cap=Math.min(classifier.config().top,classifier.familyCount());
		final AcceptingRankHistogram[] result=new AcceptingRankHistogram[panels.length];
		for(int p=0;p<result.length;p++){result[p]=new AcceptingRankHistogram(cap);}
		return result;
	}
	/** Counts and quantiles cover ranked accepts only; shortcut accepts retain explicit NA. */
	public void writeAcceptingRanks(Path path){
		require(complete && acceptingRanks!=null,"Accepting ranks require a completed telemetry evaluation");
		final ByteBuilder out=new ByteBuilder("#family_accepting_rank_summary\t1\n"+
			"panel\tclass\tqueries\taccepted\tranked_accepts\tunranked_accepts\tmean\tmedian\trank1_fraction\tp90\n");
		for(int p=0;p<panels.length;p++){
			final AcceptingRankHistogram h=acceptingRanks[p];h.validate(counts[p].accepted);
			out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative")
				.tab().append(counts[p].queries).tab().append(counts[p].accepted)
				.tab().append(h.ranked()).tab().append(h.unavailable()).tab();
			if(h.ranked()==0){out.append("NA\tNA\tNA\tNA");}
			else{
				out.appendSlow(h.mean()).tab().appendSlow(h.median()).tab()
					.appendSlow(h.rank1Fraction()).tab().append(h.p90());
			}
			out.nl();
		}
		write(path,out);
	}

	/** Bounded primitive histogram: bins[i] represents rank i+1, never rank zero.
	 * Median averages the middle pair; p90 is nearest-rank ceil(0.9*N).
	 * Unavailable shortcut ranks are counted separately and excluded from all four statistics. */
	public static final class AcceptingRankHistogram{
		public AcceptingRankHistogram(int maximumRank){
			require(maximumRank>0,"A rank histogram must have a positive strict-shortlist cap");
			bins=new long[maximumRank];
		}
		public void observe(int rank){
			require(rank==-1 || rank>=1 && rank<=bins.length,
				"An accepting rank must be within the strict shortlist, or -1 for an unranked shortcut");
			if(rank==-1){unavailable=Math.addExact(unavailable,1);}
			else{bins[rank-1]=Math.addExact(bins[rank-1],1);ranked=Math.addExact(ranked,1);}
		}
		public void mergeFrom(AcceptingRankHistogram other){
			require(other!=null && other!=this && other.bins.length==bins.length,
				"Rank reduction requires a distinct worker histogram with the same strict-shortlist cap");
			for(int i=0;i<bins.length;i++){bins[i]=Math.addExact(bins[i],other.bins[i]);}
			ranked=Math.addExact(ranked,other.ranked);unavailable=Math.addExact(unavailable,other.unavailable);
		}
		public void validate(long accepted){
			long total=0;
			for(long n:bins){require(n>=0,"Rank bins must be nonnegative counts");total=Math.addExact(total,n);}
			require(total==ranked && unavailable>=0 && Math.addExact(ranked,unavailable)==accepted,
				"Ranked and unavailable accepts must partition the independently counted terminal accepts");
		}
		public long ranked(){return ranked;}
		public long unavailable(){return unavailable;}
		public long countAt(int rank){
			require(rank>=1 && rank<=bins.length,"Histogram rank lookup is one-based within the strict shortlist");
			return bins[rank-1];
		}
		public double mean(){
			require(ranked>0,"No ranked accepts: report NA instead of a numeric mean");
			double sum=0;for(int i=0;i<bins.length;i++){sum+=(i+1.0)*bins[i];}
			return sum/ranked;
		}
		public double median(){
			require(ranked>0,"No ranked accepts: report NA instead of a numeric median");
			return (orderStatistic((ranked-1)/2+1)+(double)orderStatistic(ranked/2+1))/2;
		}
		public double rank1Fraction(){
			require(ranked>0,"No ranked accepts: report NA instead of a numeric fraction");
			return bins[0]/(double)ranked;
		}
		public int p90(){
			require(ranked>0,"No ranked accepts: report NA instead of a numeric percentile");
			return orderStatistic(ranked-ranked/10);
		}
		private int orderStatistic(long ordinal){
			require(ordinal>0 && ordinal<=ranked,"A rank order statistic must index a measured accepted query");
			long count=0;
			for(int i=0;i<bins.length;i++){count=Math.addExact(count,bins[i]);if(count>=ordinal){return i+1;}}
			throw new IllegalStateException("Accepting-rank histogram lost its recorded observations");
		}
		private final long[] bins;
		private long ranked,unavailable;
	}

	/** Final evaluation only: each compared pair is an actual candidate DP, not a unique query.
	 * Input-too-short and unavailable candidates are excluded; failed-span DPs are included. */
	public void enableOverhangDiagnostics(){
		require(!started && updates==null && classifier.config().alignmentMetric.mode==
			FamilyAlignmentMetric.Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY,
			"Overhang diagnostics require a fresh frozen overhang evaluation, never calibration updates");
		overhangDiagnostics=true;
	}
	void writeOverhangDiagnostics(Path path){
		require(complete && overhangDiagnostics,"Overhang reporting requires a completed diagnostic evaluation");
		final ByteBuilder out=new ByteBuilder(
			"panel\tclass\tqueries\tcandidate_attempts\tcompared_dp_pairs\tmaximal_trim_positive"+
			"\ttrim_positive_overhang_nonpositive\tspan_qualified_pairs\tspan_failed_pairs"+
			"\tidentity_path_eligible\tidentity_path_removed\n");
		for(int p=0;p<panels.length;p++){
			final PanelCounts c=counts[p];
			out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative").tab().append(c.queries)
				.tab().append(c.alignmentAttempts).tab().append(c.overhangCompared).tab().append(c.overhangTrimPositive)
				.tab().append(c.overhangTrimPositiveNonpositive).tab().append(c.overhangMeasured)
				.tab().append(c.overhangSpanFailed).tab().append(c.overhangIdentityEligible)
				.tab().append(c.overhangIdentityRemoved).nl();
		}
		write(path,out);
	}
	void writeAlignments(Path path){
		require(complete,"Alignment reporting requires a completed evaluation");
		final ByteBuilder out=new ByteBuilder(
			"panel\tclass\tqueries\talignment_attempts\talignment_unavailable"+
			"\talignment_calls\talignment_failed\talignment_measured\talignment_calls_per_query\n");
		for(int p=0;p<panels.length;p++){
			final PanelCounts c=counts[p];final long calls=alignmentCalls(p);
			out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative").tab().append(c.queries)
				.tab().append(c.alignmentAttempts).tab().append(c.alignmentUnavailable).tab().append(calls)
				.tab().append(c.alignmentFailures).tab().append(calls-c.alignmentFailures)
				.tab().appendSlow(calls/(double)c.queries).nl();
		}
		write(path,out);
	}
	void writeStageCpu(Path path){
		require(complete && stageCpu && stageCpuCounts!=null,"Stage CPU reporting requires a completed timed evaluation");
		final long search=stageCpuCounts.indexNanos+stageCpuCounts.shortlistNanos+stageCpuCounts.alignmentNanos;
		final long other=stageCpuCounts.workerNanos-search;
		final String[] stages={"index","shortlist","alignment","other_worker"};
		final long[] nanos={stageCpuCounts.indexNanos,stageCpuCounts.shortlistNanos,stageCpuCounts.alignmentNanos,other};
		final ByteBuilder out=new ByteBuilder("stage\tnanos\tseconds\tpercent_worker_cpu\tpercent_search_cpu\n");
		for(int i=0;i<stages.length;i++){
			out.append(stages[i]).tab().append(nanos[i]).tab().appendSlow(nanos[i]*1e-9).tab()
				.appendSlow(stageCpuCounts.workerNanos==0 ? 0 : 100.0*nanos[i]/stageCpuCounts.workerNanos).tab()
				.appendSlow(i==3 || search==0 ? 0 : 100.0*nanos[i]/search).nl();
		}
		out.append("total\t").append(stageCpuCounts.workerNanos).tab().appendSlow(stageCpuCounts.workerNanos*1e-9)
			.append("\t100\t100\n");
		write(path,out);
	}
	/** Row0 is the whole panel; remaining rows mirror length bins when requested. */
	private long[][][] newAdmissionCounts(){
		return admissionReport ? new long[panels.length]
			[lengthEdges==null ? 1 : lengthEdges.length+2][6] : null;
	}
	private static void observeAdmission(long[] row,int floor,int maximum,int pool,int dp){
		assert(row!=null && row.length==6 && floor>=0 && maximum>=floor && pool>=0 && dp>=0 && dp<=pool) :
			"Admission telemetry must describe the raw pool and its actual downstream DP attempts";
		row[3]=row[0]==0 ? floor : Math.min(row[3],floor);row[4]=Math.max(row[4],maximum);
		row[0]++;row[1]+=pool;row[2]=Math.max(row[2],pool);
		row[5]+=dp;
	}

	/** Bins are [lower,upper); underflow and overflow retain every input query, including empty bins. */
	void writeLengths(Path path){
		require(complete && lengthCounts!=null, "Length summaries require a completed evaluation with explicit boundaries");
		final ByteBuilder out=new ByteBuilder(
			"panel\tclass\tlength_min_inclusive\tlength_max_exclusive\tqueries\taccepted"+
			"\trejected\tresidues\tTPR\tFNR\tFPR\tfalse_accepts_per_million");
		final String[] stages=stages();
		for(String stage:stages){out.tab().append(stage);}out.nl();
		for(int p=0; p<panels.length; p++){
			for(int b=0; b<lengthCounts[p].length; b++){
				final long[] row=lengthCounts[p][b];
				final long queries=row[0], accepted=row[3]+row[4]+row[5];
				out.append(panels[p].name).tab().append(panels[p].positive ? "positive" : "negative")
					.tab().append(b==0 ? 0 : lengthEdges[b-1]).tab();
				if(b==lengthEdges.length){out.append("NA");}else{out.append(lengthEdges[b]);}
				out.tab().append(queries).tab().append(accepted).tab().append(queries-accepted).tab().append(row[1]);
				if(queries==0){out.append("\tNA\tNA\tNA\tNA");}
				else{
					final double rate=accepted/(double)queries;
					if(panels[p].positive){out.tab().appendSlow(rate).tab().appendSlow(1-rate).append("\tNA\tNA");}
					else{out.append("\tNA\tNA\t").appendSlow(rate).tab().appendSlow(rate*1000000);}
				}
				for(int s=0; s<stages.length; s++){out.tab().append(row[s+2]);}out.nl();
			}
		}
		write(path, out);
	}

	/** Per-panel/bin query, residue and terminal counts; no per-family array is multiplied by bin count. */
	private long[][][] newLengthCounts(){
		assert(panels!=null) : "Length counters must share the runner's already validated panel roster";
		return lengthEdges==null ? null : new long[panels.length][lengthEdges.length+1][METRIC_STAGES.length+2];
	}

	/** Every length bin must partition the ordinary panel counts exactly, with no query dropping. */
	private void validateLengths(){
		if(lengthCounts==null){return;}
		for(int p=0; p<panels.length; p++){
			final long[] totals=new long[METRIC_STAGES.length+2];
			for(long[] row:lengthCounts[p]){
				long queries=0;
				for(int s=0; s<METRIC_STAGES.length; s++){queries=Math.addExact(queries, row[s+2]);}
				require(queries==row[0], "Length-bin terminal counts disagree with its query count");
				for(int i=0; i<row.length; i++){totals[i]=Math.addExact(totals[i], row[i]);}
			}
			require(totals[0]==counts[p].queries && totals[1]==counts[p].residues,
				"Length reporting lost or duplicated panel queries or residues: "+panels[p].name);
			for(int s=0; s<METRIC_STAGES.length; s++){
				require(totals[s+2]==counts[p].stages[s], "Length reporting changed a panel terminal stage: "+METRIC_STAGES[s]);
			}
		}
	}

	/** Uses observed denominators; only positive panels have TPR/FNR fields. */
	void writeSummary(Path path){
		require(complete,"Panel results may only be published after every worker has finished");
		final ByteBuilder out=new ByteBuilder(
			"panel\tclass\tqueries\taccepted\trejected\tTPR\tFNR\tFPR\tfalse_accepts_per_million");
		final String[] stages=stages();
		for(String stage:stages){out.tab().append(stage);}out.nl();
		for(int i=0;i<panels.length;i++){
			final Panel p=panels[i];final PanelCounts c=counts[i];final double rate=c.accepted/(double)c.queries;
			out.append(p.name).tab().append(p.positive ? "positive" : "negative").tab().append(c.queries)
				.tab().append(c.accepted).tab().append(c.queries-c.accepted).tab();
			if(p.positive){out.appendSlow(rate).tab().appendSlow(1-rate).append("\tNA\tNA");}
			else{out.append("NA\tNA\t").appendSlow(rate).tab().appendSlow(rate*1000000);}
			for(int s=0; s<stages.length; s++){out.tab().append(c.stages[s]);}out.nl();
		}
		write(path,out);
	}

	/** Emits terminal winners for all callable families, including zero-count rows. */
	void writeFamilies(Path path){
		require(complete,"Per-family results require a completed evaluation");
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		try{
			writer.println("panel\troster_rank\tfamily_id\tnegative_queries\tfalse_accepts\tfalse_accepts_per_million");
			final ByteBuilder line=new ByteBuilder();
			for(int p=0;p<panels.length;p++){
				if(panels[p].positive){continue;}
				for(int f=0;f<classifier.familyCount();f++){
					final long n=counts[p].winners[f];
					line.clear().append(panels[p].name).tab().append(classifier.familyRankAt(f)).tab().append(classifier.familyIdAt(f))
						.tab().append(counts[p].queries).tab().append(n).tab().appendSlow(n*1000000.0/counts[p].queries).nl();
					writer.print(line);
				}
			}
		}finally{require(!writer.poisonAndWait(),"Per-family false-accept writer failed");}
	}

	@Override public void accumulate(ProcessThread worker){
		assert(worker.countsT.length==counts.length) : "Each worker must use the same ordered panel roster";
		if(hitHistograms!=null && worker.success){
			for(int p=0;p<panels.length;p++){
				final FamilyHitHistogram source=worker.hitHistogramsT[p];
				require(source.completedQueries()==worker.countsT[p].queries &&
					source.acceptedQueries()==worker.countsT[p].accepted,
					"Completed worker histogram must cover every terminal observation exactly once");
				hitHistograms[p].mergeFrom(source);
			}
		}
		for(int i=0;i<counts.length;i++){counts[i].add(worker.countsT[i]);}
		if(acceptingRanks!=null && worker.success){
			for(int p=0;p<panels.length;p++){
				worker.acceptingRanksT[p].validate(worker.countsT[p].accepted);
				acceptingRanks[p].mergeFrom(worker.acceptingRanksT[p]);
			}
		}
		if(lengthCounts!=null){
			for(int p=0; p<panels.length; p++){
				for(int b=0; b<lengthCounts[p].length; b++){
					for(int c=0; c<lengthCounts[p][b].length; c++){
						lengthCounts[p][b][c]=Math.addExact(lengthCounts[p][b][c], worker.lengthCountsT[p][b][c]);
					}
				}
			}
		}
		if(updates!=null){updates.add(worker.updatesT);}
		if(admissionCounts!=null){
			for(int p=0;p<panels.length;p++){
				for(int b=0;b<admissionCounts[p].length;b++){
					final long[] target=admissionCounts[p][b],source=worker.admissionCountsT[p][b];
					if(source[0]==0){continue;}
					target[3]=target[0]==0 ? source[3] : Math.min(target[3],source[3]);
					target[0]=Math.addExact(target[0],source[0]);target[1]=Math.addExact(target[1],source[1]);
					target[2]=Math.max(target[2],source[2]);target[4]=Math.max(target[4],source[4]);
					target[5]=Math.addExact(target[5],source[5]);
				}
			}
		}
		errorState|=!worker.success;
	}
	@Override public boolean success(){return !errorState && failure.get()==null;}
	@Override public ReadWriteLock rwlock(){return lock;}

	/** Streamers distribute batches; scratch and observations never cross worker threads. */
	final class ProcessThread extends Thread implements FamilyTailCascade.AssignmentObserver,
			FamilyTailCascade.RawHitObserver{
		ProcessThread(int id,Streamer[] inputs_){
			super("FamilyPanel-"+id);inputs=inputs_;
			countsT=new PanelCounts[panels.length];
			for(int i=0;i<countsT.length;i++){countsT[i]=new PanelCounts(classifier.familyCount());}
			updatesT=updates==null ? null : updates.newAccumulator();
			lengthCountsT=newLengthCounts();
			admissionCountsT=newAdmissionCounts();
			hitHistogramsT=newHitHistograms();
			acceptingRanksT=acceptingRanks==null ? null : newAcceptingRankHistograms();
		}
		@Override public void run(){
			try{
				final FamilyTailCascade.Worker scorer=classifier.newWorker(stageCpu);
				if(overhangDiagnostics){scorer.enableOverhangDiagnostics();}
				for(int panel=0;panel<inputs.length;panel++){
					final long attempts=scorer.alignmentAttempts(), unavailable=scorer.alignmentUnavailable(),
						failed=scorer.alignmentFailures();
					final long compared=scorer.overhangCompared(),trimPositive=scorer.overhangTrimPositive(),
						crossed=scorer.overhangTrimPositiveNonpositive(),measured=scorer.overhangMeasured(),
						spanFailed=scorer.overhangSpanFailed(),
						identityEligible=scorer.overhangIdentityEligible(),identityRemoved=scorer.overhangIdentityRemoved();
					for(ListNum<Read> batch=inputs[panel].nextList();batch!=null;batch=inputs[panel].nextList()){
						try{
							if(failure.get()==null){
								currentCounts=countsT[panel];currentPanel=panel;negative=!panels[panel].positive;process(batch,scorer);
							}
						}catch(Throwable problem){failure.compareAndSet(null,problem);}
					}
					countsT[panel].alignmentAttempts=scorer.alignmentAttempts()-attempts;
					countsT[panel].alignmentUnavailable=scorer.alignmentUnavailable()-unavailable;
					countsT[panel].alignmentFailures=scorer.alignmentFailures()-failed;
					countsT[panel].overhangCompared=scorer.overhangCompared()-compared;
					countsT[panel].overhangTrimPositive=scorer.overhangTrimPositive()-trimPositive;
					countsT[panel].overhangTrimPositiveNonpositive=scorer.overhangTrimPositiveNonpositive()-crossed;
					countsT[panel].overhangMeasured=scorer.overhangMeasured()-measured;
					countsT[panel].overhangSpanFailed=scorer.overhangSpanFailed()-spanFailed;
					countsT[panel].overhangIdentityEligible=scorer.overhangIdentityEligible()-identityEligible;
					countsT[panel].overhangIdentityRemoved=scorer.overhangIdentityRemoved()-identityRemoved;
				}
				stageCpuCountsT=scorer.stageCpuCounts();
				success=failure.get()==null;
			}catch(Throwable problem){
				problem.printStackTrace(outstream);
				KillSwitch.kill("Family panel execution encountered an unrecoverable input/worker failure");
			}
		}
		private void process(ListNum<Read> batch,FamilyTailCascade.Worker scorer){
			assert(batch!=null) : "A null batch is EOF and cannot contribute panel records";
			for(Read read:batch.list){
				require(read!=null && read.bases!=null && read.mate==null,"Family panels require unpaired protein sequences");
				if(panels[currentPanel].verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL){
					require(read.id!=null && !FamilyTailCascade.firstToken(read.id).isEmpty(),
						"Structural panels require nonempty query identities");
				}
				final String id=FamilyTailCascade.firstToken(read.id);
				if(queryGc==null){
					if(hitHistogramsT==null){scorer.assignObserve(id, read.bases, this);}
					else{scorer.assignObserve(id, read.bases, this, this);}
				}else{
					final FamilyQueryGc gc=queryGc[currentPanel];gc.checkQuery(read.numericID, id);
					currentGcMeasured=gc.available(read.numericID);
					scorer.assignObserve(id, read.bases, currentGcMeasured, gc.fraction(read.numericID), this,
						hitHistogramsT==null ? null : this);
				}
			}
			final long total=progress.addAndGet(batch.size()),next=nextProgress.get();
			if(total>=next && nextProgress.compareAndSet(next,total+1000000)){
				outstream.println(label+" progress: "+total+" queries scored.");
			}
		}
		@Override public void accept(final int[] touched,final int[] scores,final int size){
			require(hitHistogramsT!=null && currentCounts!=null,"Raw observation needs an active telemetry panel");
			hitHistogramsT[currentPanel].beginQuery(touched,scores,size);
		}
		@Override public void admission(int length,int floor,int maximum,int pool,int dp){
			if(admissionCountsT==null){return;}
			require(currentCounts!=null && length>=0 && floor>=0 && maximum>=floor &&
					pool>=0 && pool<=classifier.familyCount() && dp>=0 && dp<=pool,
				"Admission callback must describe one query in the currently claimed panel and the loaded roster");
			observeAdmission(admissionCountsT[currentPanel][0],floor,maximum,pool,dp);
			if(lengthEdges!=null){
				int bin=0;while(bin<lengthEdges.length && length>=lengthEdges[bin]){bin++;}
				observeAdmission(admissionCountsT[currentPanel][bin+1],floor,maximum,pool,dp);
			}
		}
		@Override public void acceptingRank(int family,int rank){
			if(acceptingRanksT==null){return;}
			require(currentCounts!=null && family>=0 && family<classifier.familyCount() && acceptingFamily==-1,
				"Each terminal winner must have exactly one rank callback in the currently claimed panel");
			acceptingFamily=family;
			acceptingRanksT[currentPanel].observe(rank);
		}
		@Override public void accept(int length,String terminal,int family,int ddl,float identity,float blosum){
			accept(length,terminal,family,ddl,identity,blosum,FamilyAlignmentMetric.State.MEASURED);
		}
		@Override public void accept(int length,String terminal,int family,int ddl,float identity,float blosum,
				FamilyAlignmentMetric.State alignmentState){
			assert(currentCounts!=null) : "The synchronous observer counts against the currently claimed panel batch";
			if(length<minQueryAa || length>maxQueryAa){
				throw new IllegalArgumentException("Panel "+panels[currentPanel].name+": normalized query length "+length+
					" is outside "+minQueryAa+".."+maxQueryAa+"; evaluation refuses to filter or truncate records");
			}
			final PanelCounts count=currentCounts;final int stage=stage(terminal);
			if(acceptingRanksT!=null){
				require(acceptingFamily==family,
					"The rank callback must identify this query's terminal winning family, or be absent for rejection");
				acceptingFamily=-1;
			}
			count.queries++;count.residues+=length;count.stages[stage]++;
			if(queryGc!=null){
				if(currentGcMeasured){count.measuredGcQueries++;if(family>=0){count.measuredGcAccepts++;}}
				else{count.neutralGcQueries++;if(family>=0){count.neutralGcAccepts++;}}
			}
			if(lengthCountsT!=null){
				require(length>=0, "The synchronous assignment observer must report a normalized nonnegative query length");
				int bin=0;
				while(bin<lengthEdges.length && length>=lengthEdges[bin]){bin++;}
				final long[] row=lengthCountsT[currentPanel][bin];
				row[0]++;row[1]+=length;row[stage+2]++;
			}
			if(family>=0){
				require(stage>=1 && stage<=3 && family<count.winners.length,
					"Accepted query has inconsistent terminal family/stage");
				count.accepted++;count.winners[family]++;
				if(negative && updatesT!=null){
					if(stage==1 || alignmentState==FamilyAlignmentMetric.State.MEASURED){
						updatesT.observeAtMargin(family,stage,ddl,identity,blosum,alignmentState,updateMarginScale);
					}else{
						throw new IllegalArgumentException("Accepted family lacks a measured alignment: "+terminal);
					}
				}
			}else{require(stage!=1 && stage!=2 && stage!=3,"Accepting stage has no winner");}
			if(hitHistogramsT!=null){hitHistogramsT[currentPanel].finishQuery(family,ddl);}
		}
		final Streamer[] inputs;
		final PanelCounts[] countsT;
		final FamilyCalibrationUpdates updatesT;
		final long[][][] lengthCountsT;
		final long[][][] admissionCountsT;
		final FamilyHitHistogram[] hitHistogramsT;
		final AcceptingRankHistogram[] acceptingRanksT;
		int acceptingFamily=-1;
		FamilyTailCascade.StageCpuCounts stageCpuCountsT;
		int currentPanel;
		PanelCounts currentCounts;
		boolean negative,success,currentGcMeasured;
	}

	public static final class Panel{
		public Panel(String name_,String path_,String pin_,boolean positive_){
			this(name_,path_,pin_,positive_,FamilyCalibration.VerificationMode.SHA80,-1);
		}
		/** Explicit expected population is checked after workers join; the inherited pin is not recomputed. */
		public Panel(String name_,String path_,String pin_,boolean positive_,
				FamilyCalibration.VerificationMode mode,long expectedQueries_){
			require(mode!=null && expectedQueries_!=0 && expectedQueries_>=-1,
				"Panel populations must be absent or positive");
			verificationMode=mode;expectedQueries=expectedQueries_;
			require(name_!=null && name_.matches("[A-Za-z0-9_.-]+"),"Panel labels must be nonempty tokens");
			name=name_;path=Paths.get(path_).toAbsolutePath().normalize();pin=pin_;positive=positive_;
			require(Files.isRegularFile(path),"Missing sequence panel: "+path);
			require(pin==null || pin.matches("[0-9a-f]{20}") ||
					mode==FamilyCalibration.VerificationMode.STRUCTURAL && pin.equals("unmeasured"),
				"Panel fingerprint requires20hex, or unmeasured with explicit structural verification");
		}
		void verify(){
			if(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL){
				require(Files.isRegularFile(path) && !Files.isSymbolicLink(path),"Missing/nonregular structural panel: "+name);
			}else if(pin!=null){require(pin.equals(DigestSuffix.file(path.toString())),"Panel fingerprint differs: "+name);}
		}
		final FamilyCalibration.VerificationMode verificationMode;
		final long expectedQueries;
		public final String name,pin;
		public final Path path;
		public final boolean positive;
	}

	static final class PanelCounts{
		PanelCounts(int families){require(families>0,"Panel counters require a live roster");winners=new long[families];}
		void add(PanelCounts other){
			assert(winners.length==other.winners.length) : "Worker reduction must retain the classifier roster";
			queries=Math.addExact(queries,other.queries);
			accepted=Math.addExact(accepted,other.accepted);
			residues=Math.addExact(residues,other.residues);
			alignmentAttempts=Math.addExact(alignmentAttempts,other.alignmentAttempts);
			alignmentUnavailable=Math.addExact(alignmentUnavailable,other.alignmentUnavailable);
			alignmentFailures=Math.addExact(alignmentFailures,other.alignmentFailures);
			overhangCompared=Math.addExact(overhangCompared,other.overhangCompared);
			overhangTrimPositive=Math.addExact(overhangTrimPositive,other.overhangTrimPositive);
			overhangTrimPositiveNonpositive=Math.addExact(overhangTrimPositiveNonpositive,other.overhangTrimPositiveNonpositive);
			overhangMeasured=Math.addExact(overhangMeasured,other.overhangMeasured);
			overhangSpanFailed=Math.addExact(overhangSpanFailed,other.overhangSpanFailed);
			overhangIdentityEligible=Math.addExact(overhangIdentityEligible,other.overhangIdentityEligible);
			overhangIdentityRemoved=Math.addExact(overhangIdentityRemoved,other.overhangIdentityRemoved);
			measuredGcQueries=Math.addExact(measuredGcQueries, other.measuredGcQueries);
			neutralGcQueries=Math.addExact(neutralGcQueries, other.neutralGcQueries);
			measuredGcAccepts=Math.addExact(measuredGcAccepts, other.measuredGcAccepts);
			neutralGcAccepts=Math.addExact(neutralGcAccepts, other.neutralGcAccepts);
			for(int i=0;i<winners.length;i++){winners[i]=Math.addExact(winners[i],other.winners[i]);}
			for(int i=0;i<stages.length;i++){stages[i]=Math.addExact(stages[i],other.stages[i]);}
		}
		void validate(String name){
			require(queries>0 && accepted>=0 && accepted<=queries,"Empty or inconsistent panel: "+name);
			require(alignmentUnavailable>=0 && alignmentFailures>=0 && alignmentAttempts>=alignmentUnavailable &&
				alignmentFailures<=alignmentAttempts-alignmentUnavailable,"Inconsistent alignment work counts: "+name);
			long hits=0,total=0;
			require(overhangCompared>=overhangTrimPositive && overhangTrimPositive>=overhangTrimPositiveNonpositive &&
				overhangTrimPositiveNonpositive>=0 && overhangCompared<=alignmentAttempts-alignmentUnavailable &&
				overhangCompared==overhangMeasured+overhangSpanFailed && overhangSpanFailed>=0 &&
				overhangMeasured>=overhangIdentityEligible &&
				overhangIdentityEligible>=overhangIdentityRemoved && overhangIdentityRemoved>=0,
				"Inconsistent overhang candidate-pair diagnostic denominators: "+name);
			for(long n:winners){require(n>=0,"Negative family accept count");hits=Math.addExact(hits,n);}
			for(long n:stages){require(n>=0,"Negative terminal count");total=Math.addExact(total,n);}
			require(hits==accepted && total==queries && stages[1]+stages[2]+stages[3]==accepted,
				"Per-family, terminal and panel counts disagree: "+name);
		}
		long queries,accepted,residues;
		long measuredGcQueries,neutralGcQueries,measuredGcAccepts,neutralGcAccepts;
		long alignmentAttempts,alignmentUnavailable,alignmentFailures;
		long overhangCompared,overhangTrimPositive,overhangTrimPositiveNonpositive;
		long overhangMeasured,overhangSpanFailed,overhangIdentityEligible,overhangIdentityRemoved;
		final long[] winners,stages=new long[METRIC_STAGES.length];
	}

	static int stage(String terminal){
		for(int i=0;i<METRIC_STAGES.length;i++){if(METRIC_STAGES[i].equals(terminal)){return i;}}
		throw new IllegalStateException("Unknown classifier terminal stage: "+terminal);
	}
	private String[] stages(){
		return classifier.config().alignmentMetric.same(FamilyAlignmentMetric.LEGACY) ? STAGES : METRIC_STAGES;
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static void write(Path path,ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		writer.print(text);require(!writer.poisonAndWait(),"Panel summary writer failed");
	}
	private final FamilyTailCascade.LoadedClassifier classifier;
	private final Panel[] panels;
	private final FamilyQueryGc[] queryGc;
	private final int threads,inputThreads;
	private final double updateMarginScale;
	private final int minQueryAa,maxQueryAa;
	private final PrintStream outstream;
	private final String label;
	final PanelCounts[] counts;
	private final int[] lengthEdges;
	private final long[][][] lengthCounts;
	private final boolean admissionReport;
	private final boolean stageCpu;
	private final boolean hitTelemetry;
	private final FamilyHitHistogram[] hitHistograms;
	private AcceptingRankHistogram[] acceptingRanks;
	private final long[][][] admissionCounts;
	private FamilyTailCascade.StageCpuCounts stageCpuCounts;
	final FamilyCalibrationUpdates updates;
	private boolean errorState,started,complete;
	private boolean overhangDiagnostics;
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private final AtomicLong progress=new AtomicLong(),nextProgress=new AtomicLong(1000000);
	private final ReadWriteLock lock=new ReentrantReadWriteLock();
	static final String[] STAGES={"ddl_reject","ddl_accept","identity_accept","blosum_accept",
		"identity_reject","blosum_reject",
		"no_hit","scalar_reject"};
	static final String[] METRIC_STAGES={"ddl_reject","ddl_accept","identity_accept","blosum_accept",
		"identity_reject","blosum_reject","no_hit","scalar_reject","alignment_unavailable",
		"input_too_short","no_alignment","span_too_short"};
}
