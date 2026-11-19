package prok;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.atomic.AtomicBoolean;

import dna.GeneticCode;
import map.LongHashMap;
import prot.FamilyTailCascade;
import prot.FamilyMatchInformation;
import prot.FamilyMatchInformation.View;
import prok.StopTailCandidateScan.Candidate;
import prok.StopTailCandidateScan.Scan;
import stream.Read;

/** Shared classifier workers and immutable per-contig stop evidence.
 * One already loaded classifier serves the whole caller invocation. Each eligible
 * stop is searched once per applicable family set, independently of alternative gene starts.
 * This layer changes no score, pruning rule or classifier acceptance threshold.
 * @author Keqing
 */
public final class StopTailCascadeSupport implements AutoCloseable{

	/** The caller supplies its pinned classifier and closes this pool after all caller workers join.
	 * A single search worker still uses the pool, so concurrent callers cannot exceed that limit. */
	public StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier, final int threads_){
		this(classifier, null, null, threads_);
	}
	public StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,final int threads,
			final FamilyMatchInformation information){
		this(classifier,null,null,threads,false,false,null,40,false,false,StopTailCandidateScan.TAIL_LENGTH,information);
	}

	/** Optional short bins use separately calibrated classifiers and never inherit the
	 * full-length margin policy. The intervening 125..149 bin is deliberately absent. */
	public StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_){
		this(classifier, short100, short150, threads_, false);
	}

	/** Explicit research binding; the supplied classifier must be calibrated for exactly100-aa queries. */
	public static StopTailCascadeSupport fixed100(final FamilyTailCascade.LoadedClassifier classifier, final int threads){
		return new StopTailCascadeSupport(classifier, null, null, threads, true);
	}
	public static StopTailCascadeSupport fixed100(final FamilyTailCascade.LoadedClassifier classifier,final int threads,
			final FamilyMatchInformation information){
		return new StopTailCascadeSupport(classifier,null,null,threads,true,false,null,40,false,false,100,information);
	}

	/** Explicit independent-family100..149 experiment;150..199 and200 retain their old rosters. */
	public static StopTailCascadeSupport independent100To149(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier small, final FamilyTailCascade.LoadedClassifier medium, final int threads){
		if(small==null || medium==null){
			throw new IllegalArgumentException("Independent100..149 dispatch requires both short classifiers");
		}
		return new StopTailCascadeSupport(full, small, medium, threads, false, true);
	}

	/** Adds an independent40..99 family set while preserving the100..149,150..199 and200 policies. */
	public static StopTailCascadeSupport independent40To99(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier small,
			final FamilyTailCascade.LoadedClassifier medium, final int threads){
		return independent40To99(full, tiny, small, medium, threads, 40);
	}

	public static StopTailCascadeSupport independent40To99(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier small,
			final FamilyTailCascade.LoadedClassifier medium, final int threads, final int minimumShortLength){
		return independent40To99(full, tiny, small, medium, threads, minimumShortLength, false);
	}

	/** Explicit query-boundary experiment; family gates and worker routing otherwise stay unchanged. */
	public static StopTailCascadeSupport independent40To99(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier small,
			final FamilyTailCascade.LoadedClassifier medium, final int threads, final int minimumShortLength,
			final boolean firstStart){
		if(tiny==null || small==null || medium==null){
			throw new IllegalArgumentException("40..99 dispatch requires all three short classifiers");
		}
		return new StopTailCascadeSupport(full, small, medium, threads, false, true, tiny, minimumShortLength, firstStart);
	}

	/** Search short-family suffixes as well as the applicable long-family view; retain one binary stop decision. */
	public static StopTailCascadeSupport allFamilySets(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier small,
			final FamilyTailCascade.LoadedClassifier medium, final int threads){
		return allFamilySets(full,tiny,small,medium,threads,true);
	}
	/** The explicit false policy retains stop-to-stop windows for comparison controls. */
	public static StopTailCascadeSupport allFamilySets(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier small,
			final FamilyTailCascade.LoadedClassifier medium, final int threads, final boolean firstStart){
		return new StopTailCascadeSupport(full,small,medium,threads,false,true,tiny,40,firstStart,true);
	}

	/** Explicit cap100 comparison: one exact100 long-family gate plus the unchanged
	 * independent short-family search. Larger intervals are clipped to100 aa. */
	public static StopTailCascadeSupport capped100(final FamilyTailCascade.LoadedClassifier exact,
			final FamilyTailCascade.LoadedClassifier tiny,final int threads,final boolean firstStart){
		return new StopTailCascadeSupport(exact,null,null,threads,false,true,tiny,40,firstStart,true,100);
	}
	public static StopTailCascadeSupport capped100(final FamilyTailCascade.LoadedClassifier exact,
			final FamilyTailCascade.LoadedClassifier tiny,final int threads,final boolean firstStart,
			final FamilyMatchInformation information){
		return new StopTailCascadeSupport(exact,null,null,threads,false,true,tiny,40,firstStart,true,100,information);
	}

	/** Explicit cap120 with independently calibrated100..119 and exact120 gates. */
	public static StopTailCascadeSupport capped120(final FamilyTailCascade.LoadedClassifier exact,
			final FamilyTailCascade.LoadedClassifier tiny, final FamilyTailCascade.LoadedClassifier lower,
			final int threads, final boolean firstStart){
		return new StopTailCascadeSupport(exact,lower,null,threads,false,true,tiny,40,firstStart,true,120);
	}

	/** Explicit cap150 comparison: one long-family roster with100..149 and exact150 gates,
	 * plus the unchanged independent short-family search. No200-aa gate is reused. */
	public static StopTailCascadeSupport capped150(final FamilyTailCascade.LoadedClassifier exact,
			final FamilyTailCascade.LoadedClassifier tiny,final FamilyTailCascade.LoadedClassifier lower,final int threads,
			final boolean firstStart){
		return new StopTailCascadeSupport(exact,lower,null,threads,false,true,tiny,40,firstStart,true,150);
	}

	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_, final boolean fixed100_){
		this(classifier, short100, short150, threads_, fixed100_, false);
	}
	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_, final boolean fixed100_, final boolean independent100_){
		this(classifier, short100, short150, threads_, fixed100_, independent100_, null, 40);
	}
	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_, final boolean fixed100_, final boolean independent100_,
			final FamilyTailCascade.LoadedClassifier short40,
			final int minimumShortLength_){
		this(classifier, short100, short150, threads_, fixed100_, independent100_, short40, minimumShortLength_, false);
	}
	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_, final boolean fixed100_, final boolean independent100_,
			final FamilyTailCascade.LoadedClassifier short40,
			final int minimumShortLength_, final boolean firstStart_){
		this(classifier, short100, short150, threads_, fixed100_, independent100_, short40,
			minimumShortLength_, firstStart_, false);
	}
	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100, final FamilyTailCascade.LoadedClassifier short150,
			final int threads_, final boolean fixed100_, final boolean independent100_,
			final FamilyTailCascade.LoadedClassifier short40,
			final int minimumShortLength_, final boolean firstStart_, final boolean allFamilySets_){
		this(classifier, short100, short150, threads_, fixed100_, independent100_, short40,
			minimumShortLength_, firstStart_, allFamilySets_,
			fixed100_ ? 100 : StopTailCandidateScan.TAIL_LENGTH);
	}
	private StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100,final FamilyTailCascade.LoadedClassifier short150,
			final int threads_,final boolean fixed100_,final boolean independent100_,
			final FamilyTailCascade.LoadedClassifier short40,
			final int minimumShortLength_,final boolean firstStart_,final boolean allFamilySets_,final int maximumTailLength_){
		this(classifier,short100,short150,threads_,fixed100_,independent100_,short40,
			minimumShortLength_,firstStart_,allFamilySets_,maximumTailLength_,null);
	}
	/** Resource loader binds evidence before any worker or search is created. */
	StopTailCascadeSupport(final FamilyTailCascade.LoadedClassifier classifier,
			final FamilyTailCascade.LoadedClassifier short100,final FamilyTailCascade.LoadedClassifier short150,
			final int threads_,final boolean fixed100_,final boolean independent100_,
			final FamilyTailCascade.LoadedClassifier short40,
			final int minimumShortLength_,final boolean firstStart_,final boolean allFamilySets_,final int maximumTailLength_,
			final FamilyMatchInformation information){
		final boolean capped100=maximumTailLength_==100 && !fixed100_;
		final boolean cappedLower=maximumTailLength_==120 || maximumTailLength_==150;
		if((capped100 || cappedLower) && (fixed100_ || !allFamilySets_ || !independent100_ ||
			short150!=null || capped100 && short100!=null) ||
			!(capped100 || cappedLower) && maximumTailLength_!=(fixed100_ ? 100 : StopTailCandidateScan.TAIL_LENGTH)){
			throw new IllegalArgumentException("Capped composite callers require explicit all-family-set search; "
				+"other callers retain their original cap");
		}
		if(allFamilySets_ && (fixed100_ || !independent100_ || minimumShortLength_!=40 || short40==null ||
			(capped100 ? short100!=null || short150!=null : short100==null || !cappedLower && short150==null))){
			throw new IllegalArgumentException("All-family-set search requires the explicitly calibrated bins "
				+"for its cap and minimum40");
		}
		if(allFamilySets_ && (short40.usesLengthScaledShortAdmission() || short40.usesFamilyMeanAdmission())){
			throw new IllegalArgumentException("All-family-set and length-based admission comparisons must remain separate");
		}
		allFamilySets=allFamilySets_;
		maximumTailLength=maximumTailLength_;
		if(minimumShortLength_<40 || minimumShortLength_>99){
			throw new IllegalArgumentException("Short-bin scanner minimum must be in [40,99]");
		}
		if(firstStart_ && (fixed100_ || short40==null || !independent100_)){
			throw new IllegalArgumentException("First-start queries require an explicit independent short-family binding");
		}
		firstStart=firstStart_;
		minimumShortLength=minimumShortLength_;
		fixed100=fixed100_;
		independent100=independent100_;
		if(classifier==null || threads_<1){
			throw new IllegalArgumentException("A loaded classifier and positive search worker count are required");
		}
		if(!independent100 || allFamilySets){requireSameRoster(classifier, short100);}
		requireSameRoster(classifier, short150);
		bin40=short40!=null; bin100=short100!=null; bin150=short150!=null;
		gcBound=classifier.usesFamilyGcPenalty() || short40!=null && short40.usesFamilyGcPenalty() ||
			short100!=null && short100.usesFamilyGcPenalty() || short150!=null && short150.usesFamilyGcPenalty();
		if(gcBound && (!allFamilySets || maximumTailLength!=100 || short100!=null || short150!=null)){
			throw new IllegalArgumentException("Nucleotide GC caller support requires exact LONG100/SHORT40..99 "
				+"all-set routing");
		}
		alignmentOnly=classifier.hasAlignmentOnlyAcceptance() &&
			(short100==null || short100.hasAlignmentOnlyAcceptance())
			&& (short150==null || short150.hasAlignmentOnlyAcceptance())
			&& (short40==null || short40.hasAlignmentOnlyAcceptance());
		final int[] rawScales=new int[4];
		final View[] rawViews=new View[4];
		final String[] rawBindings=new String[4];
		int rawCount=0;
		rawBindings[rawCount]="long"; rawScales[rawCount]=classifier.alignmentMatrix().scale;
		rawViews[rawCount]=View.LONG; rawCount++;
		if(short40!=null){
			rawBindings[rawCount]="short40"; rawScales[rawCount]=short40.alignmentMatrix().scale;
			rawViews[rawCount]=View.SHORT; rawCount++;
		}
		if(short100!=null){
			rawBindings[rawCount]="short100"; rawScales[rawCount]=short100.alignmentMatrix().scale;
			rawViews[rawCount]=independent100_ && !allFamilySets_ ? View.SHORT : View.LONG; rawCount++;
		}
		if(short150!=null){
			rawBindings[rawCount]="short150"; rawScales[rawCount]=short150.alignmentMatrix().scale;
			rawViews[rawCount]=View.LONG; rawCount++;
		}
		rawScoreBindingCount=rawCount;
		rawScoreMatrixScales=rawScales;
		rawScoreViews=rawViews;
		rawScoreBindings=rawBindings;
		informationBound=information!=null;
		if(informationBound){
			classifier.requireMatchInformation();
			if(short100!=null){short100.requireMatchInformation();}
			if(short150!=null){short150.requireMatchInformation();}
			if(short40!=null){short40.requireMatchInformation();}
		}
		threads=threads_;
		final FamilyMatchInformation.Binding longInfo=information==null ? null : information.bind(View.LONG);
		final FamilyMatchInformation.Binding shortInfo=information==null ? null : information.bind(View.SHORT);
		scratch=ThreadLocal.withInitial(()->new FamilyTailCascade.Worker[]{classifier.newWorker(longInfo),
			short100==null ? null : short100.newWorker(independent100_ && !allFamilySets_ ? shortInfo : longInfo),
			short150==null ? null : short150.newWorker(longInfo),
			short40==null ? null : short40.newWorker(shortInfo)});
		final AtomicInteger names=new AtomicInteger();
		pool=Executors.newFixedThreadPool(threads, task->{
			final Thread thread=new Thread(task, "CallGenes-TailCascade-"+names.incrementAndGet());
			thread.setDaemon(true);
			return thread;
		});
	}

	/** Check once before enabling a graded caller; shortcuts have no calibrated alignment margin. */
	public void requireAlignmentMargins(){
		if(allFamilySets){
			throw new IllegalArgumentException("All-family-set support uses one binary bonus; "
				+"a combined graded-margin policy is not defined");
		}
		if(!alignmentOnly){
			throw new IllegalArgumentException("Margin bonuses require alignment-only acceptance in every loaded query bin; "
				+"disable DDL/cosine shortcuts in the calibrated package");
		}
	}
	/** The caller invokes this only for an enabled nonzero information bonus. */
	public void requireMatchInformation(){
		if(!informationBound){
			throw new IllegalArgumentException("Information bonus requires frozen background and view normalizers "
				+"bound before classification");
		}
	}
	/** Formula experiments retain primary-view precedence and need measured alignment evidence. */
	public void requireFormulaEvidence(final prot.FamilyBonusFormula formula){
		requireFormulaEvidence(formula, false);
	}
	/** Formula experiments retain primary-view precedence and need measured alignment evidence. */
	public void requireFormulaEvidence(final prot.FamilyBonusFormula formula,final boolean rawScoreI0Override){
		requireFormulaInputs(formula);
		for(int i=0; i<rawScoreBindingCount; i++){
			formula.requireRawScoreI0(rawScoreMatrixScales[i], rawScoreViews[i],
				rawScoreBindings[i], rawScoreI0Override);
		}
	}

	/** Return the effective immutable formula before recipe emission or worker publication. */
	public prot.FamilyBonusFormula bindFormulaEvidence(final prot.FamilyBonusFormula formula,
			final boolean rawScoreI0Override){
		requireFormulaInputs(formula);
		prot.FamilyBonusFormula bound=formula;
		for(int i=0; i<rawScoreBindingCount; i++){
			bound=bound.bindRawScoreI0(rawScoreMatrixScales[i],rawScoreViews[i],
				rawScoreBindings[i],rawScoreI0Override);
		}
		return bound;
	}

	/** Validate alignment evidence and provenance required by the requested formula. */
	private void requireFormulaInputs(final prot.FamilyBonusFormula formula){
		if(formula==null || !alignmentOnly){
			throw new IllegalArgumentException("Formula scoring requires alignment-only acceptance in every loaded family set");
		}
		if(formula.requiresInformation()){requireMatchInformation();}
		if(formula.requiresFamilySetView() && !allFamilySets){
			throw new IllegalArgumentException("Per-view formula scoring requires explicit primary family-set provenance");
		}
	}

	/** These packages are alternative gates for the same frozen families, not different databases. */
	private static void requireSameRoster(final FamilyTailCascade.LoadedClassifier full,
			final FamilyTailCascade.LoadedClassifier shorter){
		if(shorter==null){return;}
		if(full.familyCount()!=shorter.familyCount()){
			throw new IllegalArgumentException("Short and full classifier family counts differ");
		}
		for(int i=0; i<full.familyCount(); i++){
			if(full.familyRankAt(i)!=shorter.familyRankAt(i) || !full.familyIdAt(i).equals(shorter.familyIdAt(i))){
				throw new IllegalArgumentException("Short and full classifier rosters differ at ordinal "+i);
			}
		}
	}

	/** Enumerate only lengths for which this invocation has an explicitly loaded classifier. */
	public Scan scan(final String contig, final byte[] bases, final GeneticCode code){
		if(gcBound){return StopTailCandidateScan.scanCappedWithGc(contig,bases,code,minimumShortLength,firstStart);}
		if(maximumTailLength!=StopTailCandidateScan.TAIL_LENGTH && !fixed100){
			return StopTailCandidateScan.scanCapped(contig,bases,code,minimumShortLength,maximumTailLength,firstStart);
		}
		return fixed100 ? StopTailCandidateScan.scanFixed100(contig, bases, code)
			: StopTailCandidateScan.scanBins(contig, bases, code, bin40, bin100, independent100,
				bin150, minimumShortLength, firstStart);
	}
	public boolean usesFirstStart(){return firstStart;}
	public boolean usesAllFamilySets(){return allFamilySets;}
	/** Configured amino-acid cap, including empty scans; immutable after the pool closes. */
	public int maximumTailLength(){return maximumTailLength;}
	private int bin(final int length){
		if(fixed100){
			if(length!=100){throw new IllegalArgumentException("Fixed100 classifier requires exactly100-aa queries");}
			return 0;
		}
		if(length==maximumTailLength){return 0;}
		if(bin40 && length>=minimumShortLength && length<=99){return 3;}
		if(bin100 && length>=100 && length<=(independent100 ? 149 : 124)){return 1;}
		if(bin150 && length>=150 && length<=199){return 2;}
		throw new IllegalArgumentException("No calibrated stop-tail classifier is bound for query length "+length);
	}

	/** May be called by multiple contig workers. Retains decisions and coordinates, not query arrays.
	 * Returned evidence is reused for every start and every bonus arm without another search.
	 * Input scan/code/query arrays must remain unchanged until this call completes. */
	public ContigEvidence classify(final Scan scan){
		if(scan==null || closed){
			throw new IllegalStateException("Stop-tail classification requires a scan and an open pool");
		}
		if(scan.maximumTailLength!=maximumTailLength || scan.firstStart!=firstStart){
			throw new IllegalArgumentException("Scan cap or start boundary differs from classifier query policy");
		}
		final int count=scan.candidates.size();
		final LongHashMap index=new LongHashMap(Math.max(16, count));
		final long[] keys=new long[count];
		final int[] upstream=new int[count];
		final int[] queryLengths=new int[count];
		final int[] boundaries=new int[count];
		final FamilyTailCascade.Assignment[] decisions=new FamilyTailCascade.Assignment[count];
		final FamilyTailCascade.Assignment[] shortDecisions=allFamilySets ? new FamilyTailCascade.Assignment[count] : null;
		final FamilyTailCascade.Assignment[] longDecisions=allFamilySets ? new FamilyTailCascade.Assignment[count] : null;
		for(int i=0; i<count; i++){
			final Candidate candidate=scan.candidates.get(i);
			bin(candidate.tail.length);// Fail before worker dispatch instead of transferring a 200-aa cutoff.
			if(candidate.queryCodons>candidate.upstreamCodons
					|| (!firstStart && candidate.queryCodons!=candidate.upstreamCodons) ||
				candidate.tail.length!=Math.min(candidate.queryCodons, scan.maximumTailLength) || !index.put(candidate.key(), i)){
				throw new IllegalArgumentException("Stop-tail scan contains a wrong-length query or duplicate key in "+scan.contig);
			}
			keys[i]=candidate.key(); upstream[i]=candidate.upstreamCodons;
			queryLengths[i]=candidate.tail.length;
			boundaries[i]=candidate.intervalBoundary;
		}
		if(count>0){assign(scan,decisions,shortDecisions,longDecisions);}
		int accepted=0;
		for(int i=0; i<count; i++){
			if(decisions[i]==null){
				throw new IllegalStateException("Missing stop-tail decision for "+scan.contig+" key="+keys[i]);
			}
			if(allFamilySets){
				if(shortDecisions[i]==null || (queryLengths[i]>=100)!=(longDecisions[i]!=null)){
					throw new IllegalStateException("All-family-set evidence must retain every applicable search at "
						+scan.contig+" key="+keys[i]);
				}
				final FamilyTailCascade.Assignment longResult=longDecisions[i],shortResult=shortDecisions[i];
				final FamilyTailCascade.Assignment primary=longResult==null ? shortResult :
					longResult.assigned() || !shortResult.assigned() ? longResult : shortResult;
				if(decisions[i]!=primary){
					throw new IllegalStateException("Stop decision must preserve the declared long-first acceptance precedence");
				}
			}
			if(decisions[i].assigned()){accepted++;}
		}
		return new ContigEvidence(scan.contig, scan.length, scan.code, index, keys, upstream,
			queryLengths, boundaries, decisions, accepted, scan.minimumTailLength,
			scan.lengthBins, scan.maximumTailLength, firstStart, scan.noStartStops,
			shortDecisions, longDecisions);
	}

	/** Classify once for repeated native/bonus arms. Own one private nucleotide copy
	 * so changing or recycling the input Read cannot attach cached evidence to other bases.
	 * The copy is retained only for this contig; no query arrays survive classification. */
	public Replay prepareReplay(final Read input, final GeneticCode code){
		if(input==null || input.bases==null || input.id==null || code==null){
			throw new IllegalArgumentException("Tail replay requires a named nucleotide contig and explicit code");
		}
		final Read copy=new Read(input.bases.clone(), null, input.id, input.numericID, false);
		final Scan scan=scan(copy.id, copy.bases, code);
		return new Replay(this, copy, classify(scan), scan.under100, scan.from100To199);
	}

	/** A single contig and its immutable decisions, reused sequentially across caller arms.
	 * Calls use the ordinary GeneCaller path, which restores nucleotide orientation.
	 * The classifier pool may be closed after preparation; replay does not use it. */
	public static final class Replay{
		private Replay(final StopTailCascadeSupport owner_, final Read input_, final ContigEvidence evidence_,
				final long under100_, final long from100To199_){
			owner=owner_; input=input_; evidence=evidence_; under100=under100_; from100To199=from100To199_;
		}
		/** Each caller must be thread-confined and already configured for native calls
		 * or this exact classifier. Bonuses/protection belong to the caller, not the cache. */
		public ArrayList<Orf> call(final GeneCaller caller, final GeneModel model){
			return call(caller,model,-1,0);
		}
		/** Diagnostic intervention on one accepted stop, not a new classifier decision.
		 * The caller retains its ordinary protection and all other stop bonuses. A later
		 * ordinary call has no override; success establishes only recovery at this tested bonus. */
		public ArrayList<Orf> callWithStopBonus(final GeneCaller caller, final GeneModel model,
				final int stopLeft, final int strand, final float bonus){
			final int ordinal=evidence.ordinalAtStop(stopLeft,strand);
			if(ordinal<0 || !evidence.assignment(ordinal).assigned()){
				throw new IllegalArgumentException("Targeted replay requires a stop accepted by its frozen family evidence");
			}
			return call(caller,model,ordinal,bonus);
		}
		private ArrayList<Orf> call(final GeneCaller caller, final GeneModel model, final int target, final float bonus){
			if(caller==null || model==null){throw new IllegalArgumentException("Replay requires a caller and fixed PGM");}
			if(!inUse.compareAndSet(false, true)){throw new IllegalStateException("A contig cannot be replayed concurrently");}
			try{return target<0 ? caller.callGenesWithPreparedTail(input,model,this)
				: caller.callGenesWithTargetedTailBonus(input,model,this,target,bonus);}
			finally{inUse.set(false);}
		}
		public int queryCount(){return evidence.size();}
		public int acceptedCount(){return evidence.accepted;}
		/** Immutable stop coordinates and decisions for observation; no query or input arrays escape. */
		public ContigEvidence evidence(){return evidence;}
		final StopTailCascadeSupport owner;
		final ContigEvidence evidence;
		final long under100, from100To199;
		private final Read input;
		private final AtomicBoolean inUse=new AtomicBoolean();
	}

	/** Reuses StopTailSupport's coarse chunk-claim pattern; futures surface worker failures to the caller. */
	private void assign(final Scan scan, final FamilyTailCascade.Assignment[] decisions,
			final FamilyTailCascade.Assignment[] shortDecisions,final FamilyTailCascade.Assignment[] longDecisions){
		assert(decisions.length==scan.candidates.size() && decisions.length>0) :
			"Each scanned stop needs one stable output slot";
		final int chunks=1+(decisions.length-1)/BATCH;
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Callable<Void> task=()->{
			try{
				final FamilyTailCascade.Worker[] workers=scratch.get();
				final byte[] suffix=allFamilySets ? shortSuffix.get() : null;
				for(int chunk=next.getAndIncrement(); chunk<chunks && failure.get()==null; chunk=next.getAndIncrement()){
					if(Thread.currentThread().isInterrupted()){throw new IllegalStateException("Interrupted classifier worker");}
					final int first=chunk*BATCH, end=first+Math.min(BATCH, decisions.length-first);
					for(int i=first; i<end; i++){
						final Candidate candidate=scan.candidates.get(i);
						try{
							// Truth diagnostics are disabled in the callable API; contig ID suffices without a temporary ID per stop.
							if(gcBound && !candidate.gcAvailable){
								throw new IllegalArgumentException("GC caller cannot reuse a scan without measured nucleotide tails");
							}
							final FamilyTailCascade.Worker worker=workers[bin(candidate.tail.length)];
							final FamilyTailCascade.Assignment original=gcBound ?
								worker.assign(scan.contig,candidate.tail,true,candidate.tailGc()) :
								worker.assign(scan.contig,candidate.tail);
							if(!allFamilySets){decisions[i]=original;}
							else if(candidate.tail.length<100){shortDecisions[i]=decisions[i]=original;}
							else{
								System.arraycopy(candidate.tail,candidate.tail.length-99,suffix,0,99);
								final FamilyTailCascade.Assignment shorter=gcBound ?
									workers[3].assign(scan.contig,suffix,true,candidate.shortTailGc()) :
									workers[3].assign(scan.contig,suffix);
								longDecisions[i]=original;shortDecisions[i]=shorter;
								// Always complete both searches, then keep A37 precedence when both accept.
								decisions[i]=original.assigned() || !shorter.assigned() ? original : shorter;
							}
						}catch(RuntimeException | Error e){
							throw new IllegalStateException("Classifier failed at "+scan.contig+" stop="+candidate.stopLeft+
								" strand="+candidate.strand, e);
						}
					}
				}
				return null;
			}catch(RuntimeException | Error e){failure.compareAndSet(null, e); throw e;}
		};
		final ArrayList<Callable<Void>> tasks=new ArrayList<Callable<Void>>(Math.min(threads, chunks));
		for(int i=0; i<Math.min(threads, chunks); i++){tasks.add(task);}
		try{
			for(Future<Void> future : pool.invokeAll(tasks)){future.get();}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Interrupted while waiting for stop-tail classification", e);
		}catch(ExecutionException e){throw new IllegalStateException("Stop-tail classifier worker failed", e.getCause());}
		if(failure.get()!=null){throw new IllegalStateException("Incomplete stop-tail classification", failure.get());}
	}

	/** Close after callers have finished. A concurrent new request fails rather than returning partial evidence. */
	@Override
	public void close(){
		closed=true;
		pool.shutdown();
		try{
			while(!pool.awaitTermination(1, TimeUnit.DAYS)){}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Interrupted while closing stop-tail classifier workers", e);
		}
	}

	/** Coordinates and decisions are scoped to this contig and exact code, including custom code objects.
	 * No mutable map or array is exposed. A missing key means no eligible query; an unassigned
	 * nonnull decision means that an eligible query was searched and rejected. */
	public static final class ContigEvidence{
		private ContigEvidence(final String contig_, final int length_, final GeneticCode code_, final LongHashMap index_,
				final long[] keys_, final int[] upstream_, final int[] queryLengths_, final int[] boundaries_,
				final FamilyTailCascade.Assignment[] decisions_, final int accepted_,
				final int minimum_, final int bins_, final int maximum_, final boolean firstStart_, final long noStartStops_,
				final FamilyTailCascade.Assignment[] shortDecisions_,final FamilyTailCascade.Assignment[] longDecisions_){
			contig=contig_; length=length_; code=code_; index=index_; keys=keys_; upstream=upstream_;
			decisions=decisions_; accepted=accepted_;
			minimumTailLength=minimum_; maximumTailLength=maximum_; lengthBins=bins_;
			queryLengths=queryLengths_; firstStart=firstStart_;
			boundaries=boundaries_; noStartStops=noStartStops_;
			shortDecisions=shortDecisions_;longDecisions=longDecisions_;
			if((shortDecisions==null)!=(longDecisions==null) || shortDecisions!=null &&
				(shortDecisions.length!=decisions.length || longDecisions.length!=decisions.length)){
				throw new IllegalArgumentException("Per-set evidence must cover the same unique-stop roster");
			}
		}
		public boolean searchesAllFamilySets(){return shortDecisions!=null;}
		/** Null means that set was not searched, never an observed rejection. */
		public FamilyTailCascade.Assignment familySetAssignment(final int ordinal,final boolean shortSet){
			if(!searchesAllFamilySets()){
				throw new IllegalStateException("Per-set evidence requires explicit all-family-set search");
			}
			return shortSet ? shortDecisions[ordinal] : longDecisions[ordinal];
		}
		public int familySetQueryLength(final int ordinal,final boolean shortSet){
			if(familySetAssignment(ordinal,shortSet)==null){
				throw new IllegalStateException("An unsearched set has no measured query length");
			}
			return shortSet ? Math.min(99,queryLengths[ordinal]) : queryLengths[ordinal];
		}
		public String primaryFamilySet(final int ordinal){
			if(!searchesAllFamilySets()){
				throw new IllegalStateException("Primary-set provenance requires explicit all-family-set search");
			}
			return !decisions[ordinal].assigned() ? "none" : decisions[ordinal]==shortDecisions[ordinal] ? "short" : "long";
		}
		/** Actual winning family-set provenance, independent of peptide/trimmed length.
		 * Null means rejected or unavailable under legacy single-set evidence; callers
		 * requiring per-view normalization must reject unavailable provenance. */
		public View primaryView(final int ordinal){
			if(ordinal<0 || ordinal>=decisions.length){
				throw new IllegalArgumentException("Primary view requires an existing stop ordinal");
			}
			if(!searchesAllFamilySets() || !decisions[ordinal].assigned()){return null;}
			if(decisions[ordinal]==shortDecisions[ordinal]){return View.SHORT;}
			if(decisions[ordinal]!=longDecisions[ordinal]){
				throw new IllegalStateException("Accepted primary decision must be one of its recorded family-set results");
			}
			return View.LONG;
		}
		public boolean emitsLength(final int length){
			if(length<minimumTailLength || length>maximumTailLength){return false;}
			return (lengthBins&StopTailCandidateScan.lengthBit(length))!=0;
		}
		public int queryLength(final int ordinal){return queryLengths[ordinal];}
		public int intervalBoundary(final int ordinal){return boundaries[ordinal];}
		public FamilyTailCascade.Assignment atStop(final int stopLeft, final int strand){
			final int ordinal=ordinalAtStop(stopLeft, strand);
			return ordinal<0 ? null : decisions[ordinal];
		}
		/** Stable index for per-contig accounting; a negative value means no eligible query. */
		public int ordinalAtStop(final int stopLeft, final int strand){
			return index.get(StopTailCandidateScan.stopKey(stopLeft, strand));
		}
		public int size(){return decisions.length;}
		public long key(final int ordinal){return keys[ordinal];}
		public int upstreamCodons(final int ordinal){return upstream[ordinal];}
		public FamilyTailCascade.Assignment assignment(final int ordinal){return decisions[ordinal];}
		public final String contig;
		public final int length, accepted;
		public final int minimumTailLength, maximumTailLength, lengthBins;
		public final boolean firstStart;
		public final long noStartStops;
		public final GeneticCode code;
		private final LongHashMap index;
		private final long[] keys;
		private final int[] upstream, queryLengths, boundaries;
		private final FamilyTailCascade.Assignment[] decisions;
		private final FamilyTailCascade.Assignment[] shortDecisions,longDecisions;
	}

	private final int threads, minimumShortLength, maximumTailLength;
	private final ExecutorService pool;
	private final ThreadLocal<FamilyTailCascade.Worker[]> scratch;
	private final boolean bin40, bin100, bin150, fixed100, independent100, alignmentOnly, firstStart,allFamilySets;
	private final boolean informationBound;
	private final boolean gcBound;
	private final int rawScoreBindingCount;
	private final int[] rawScoreMatrixScales;
	private final View[] rawScoreViews;
	private final String[] rawScoreBindings;
	private final ThreadLocal<byte[]> shortSuffix=ThreadLocal.withInitial(()->new byte[99]);
	private volatile boolean closed;
	private static final int BATCH=32;
}
