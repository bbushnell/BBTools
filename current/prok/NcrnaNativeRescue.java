package prok;

import java.util.Arrays;

import structures.ByteBuilder;

/** Worker-owned state for the explicit rejected shared-window experiment.
 * Actual attempts survive all ordinary verifier stages. Only primary,
 * span-eligible alignment scores supply the gate. Candidate counts and actual
 * nested request counts are separate quantities.
 * @author Ganyu
 */
public final class NcrnaNativeRescue{

	public NcrnaNativeRescue(final int models, final int count_, final double minimumScore_){
		require(models>0 && count_>0 && count_<=models && Double.isFinite(minimumScore_),
			"Native rescue needs a real model axis, bounded total candidate budget and finite score gate");
		count=count_;
		minimumScore=minimumScore_;
		attempted=new boolean[models];
		nativeOrder=new int[models];
		selected=new int[count];
	}

	/** Reset once before the complete primary/edge sequence, never between stages. */
	public void begin(){
		require(!active, "Native rescue windows cannot nest or omit their request-accounting completion");
		Arrays.fill(attempted, false);
		primaryPresent=false;
		primaryBest=0;
		reached=false;
		rescuing=false;
		gateChecked=false;
		gateEligible=false;
		selectedCount=0;
		windowDirect=0;
		active=true;
	}

	/** Called at kernel entry, after length skips; repeated ordinary-stage attempts
	 * are legal, but a new native rescue must never repeat an attempted model. */
	public void aligned(final int model){
		require(active && model>=0 && model<attempted.length,
			"Actual model attempts require the active window and original library identity");
		if(rescuing){
			require(!attempted[model], "Native rescue must exclude every earlier actual model attempt");
			windowDirect++;
			require(windowDirect<=selectedCount, "Direct rescue kernels cannot exceed the offered total model budget");
		}
		attempted[model]=true;
	}

	/** A primary score can be negative or zero; absence is an independent flag. */
	public void primaryScore(final double score){
		require(active && !rescuing && Double.isFinite(score),
			"Only a real finite span-eligible primary alignment may update the rescue gate");
		if(!primaryPresent || score>primaryBest){primaryBest=score;}
		primaryPresent=true;
	}

	/** Visited ranks include prealignment length skips, matching saved sizing. */
	public void completedStage(final int visited, final int threshold){
		require(active && !rescuing && visited>=0 && visited<=NcrnaScorePatience.CEILING
			&& threshold>0 && threshold<=NcrnaScorePatience.CEILING,
			"Rescue eligibility uses actual primary/edge visited ranks and the declared positive threshold");
		reached|=visited>=threshold;
	}

	/** Invoke only after every ordinary path returned no accepted proposal. */
	public boolean eligible(){
		require(active && !rescuing && !gateChecked, "Evaluate the native-rescue gate once after the normal paths reject");
		gateChecked=true;
		rejectedWindows++;
		if(!reached){return false;}
		thresholdWindows++;
		if(!primaryPresent){missingScoreWindows++;return false;}
		if(primaryBest<minimumScore){return false;}
		gateWindows++;
		gateEligible=true;
		return true;
	}

	/** The supplied prefix is the complete native-eligible order, either native
	 * ranking or the explicitly configured rescue profile, never the primary cap.
	 * All members share one budget, including candidates later length-skipped. */
	public int select(final int[] order, final int size){
		require(active && !rescuing && gateEligible && order!=null && size>=0 && size<=order.length,
			"Native rescue selection needs one complete actual native-order prefix");
		selectedCount=0;
		for(int rank=0; rank<size && selectedCount<count; rank++){
			final int model=order[rank];
			require(model>=0 && model<attempted.length, "Native ranking must use the bound original model axis");
			if(attempted[model]){continue;}
			for(int previous=0; previous<selectedCount; previous++){
				require(selected[previous]!=model, "Native ranking cannot offer a duplicate rescue candidate");
			}
			selected[selectedCount++]=model;
		}
		offeredModels+=selectedCount;
		rescuing=true;
		return selectedCount;
	}

	public int selected(final int index){
		require(rescuing && index>=0 && index<selectedCount, "Only the active selected candidate prefix is readable");
		return selected[index];
	}

	/** total and rescue are actual enclosing logical-request deltas, including
	 * downstream refinement. Neither is inferred from the candidate budget. */
	public void finish(final long total, final long rescue, final boolean rescued){
		require(active && total>=rescue && rescue>=windowDirect && (!rescued || rescuing)
			&& (rescuing || rescue==0),
			"Actual main and rescue request deltas must conserve every nested request in one window");
		windows++;
		mainRequests+=total-rescue;
		rescueRequests+=rescue;
		directRequests+=windowDirect;
		if(rescued){rescuedWindows++;}
		active=false;
		rescuing=false;
	}

	public boolean primaryPresent(){return primaryPresent;}
	public double primaryBest(){
		require(primaryPresent, "Missing primary score must not be interpreted as numeric zero");
		return primaryBest;
	}
	public long mainRequests(){return mainRequests;}
	public long rescueRequests(){return rescueRequests;}
	public long directRequests(){return directRequests;}
	public void addTo(final long[] totals){
		require(!active && totals!=null && totals.length==10, "Merge only completed input counters on the declared accounting axis");
		totals[0]+=windows;
		totals[1]+=rejectedWindows;
		totals[2]+=thresholdWindows;
		totals[3]+=missingScoreWindows;
		totals[4]+=gateWindows;
		totals[5]+=offeredModels;
		totals[6]+=directRequests;
		totals[7]+=rescuedWindows;
		totals[8]+=mainRequests;
		totals[9]+=rescueRequests;
	}

	public ByteBuilder append(final ByteBuilder out){
		require(!active && out!=null, "Publish input totals only after its final scheduled window completes");
		return out.append(count).tab().append(Double.toString(minimumScore)).tab().append(windows)
			.tab().append(rejectedWindows).tab().append(thresholdWindows).tab().append(missingScoreWindows)
			.tab().append(gateWindows).tab().append(offeredModels).tab().append(directRequests)
			.tab().append(rescuedWindows).tab().append(mainRequests).tab().append(rescueRequests)
			.tab().append(mainRequests+rescueRequests);
	}

	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalStateException(reason);}
	}
	public static final String HEADER="candidate_budget\tminimum_primary_score\twindows\tnormal_rejections"
		+"\tthreshold_rejections\tmissing_primary_score\tgate_passed_windows\toffered_models"
		+"\tdirect_rescue_requests\trescued_window_proposals\tmain_requests\trescue_requests\ttotal_requests";
	public final int count;
	public final double minimumScore;
	final int[] nativeOrder;
	final int[] positions=new int[7];
	private final boolean[] attempted;
	private final int[] selected;
	private int selectedCount, windowDirect;
	private boolean active, reached, primaryPresent, rescuing, gateChecked, gateEligible;
	private double primaryBest;
	private long windows, rejectedWindows, thresholdWindows, missingScoreWindows, gateWindows;
	private long offeredModels, directRequests, rescuedWindows, mainRequests, rescueRequests;
}
