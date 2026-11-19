package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;

/** Fits partition weights from immutable, already counted window/model matches.
 * No query hashing, sequence reading or alignment occurs during a trial.
 * Borrowed inputs must not change for the lifetime of this single-owner engine.
 * @author Ganyu
 */
public final class NcrnaPartitionWeightSearch{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaPartitionWeightSearch(final int models_, final char[][] counts_,
			final int[][] parents_, final long[][][] acceptable_, final int[] strata_,
			final int[] sides_, final int[] bins_){
		this(models_, 3, counts_, parents_, acceptable_, strata_, sides_, bins_);
	}

	public NcrnaPartitionWeightSearch(final int models_, final int terms_, final char[][] counts_,
			final int[][] parents_, final long[][][] acceptable_, final int[] strata_,
			final int[] sides_, final int[] bins_){
		this(models_, terms_, counts_, parents_, acceptable_, strata_, sides_, bins_, DRAWS);
	}

	public NcrnaPartitionWeightSearch(final int models_, final int terms_, final char[][] counts_,
			final int[][] parents_, final long[][][] acceptable_, final int[] strata_,
			final int[] sides_, final int[] bins_, final int randomDraws_){
		this(models_, terms_, counts_, parents_, acceptable_, strata_, sides_, bins_, randomDraws_, terms_<=2 ? 0 : 1);
	}

	public NcrnaPartitionWeightSearch(final int models_, final int terms_, final char[][] counts_,
			final int[][] parents_, final long[][][] acceptable_, final int[] strata_,
			final int[] sides_, final int[] bins_, final int randomDraws_, final int referenceTerm_){
		require(randomDraws_>0, "The predeclared search must contain a positive number of random draws");
		randomDraws=randomDraws_;
		require(referenceTerm_>=0 && referenceTerm_<terms_, "The fixed scale term must belong to this layout's actual score terms");
		referenceTerm=referenceTerm_;
		require(models_>0 && models_<=65536 && counts_!=null && parents_!=null && acceptable_!=null
			&& strata_!=null && sides_!=null && bins_!=null && terms_>=1 && terms_<=4,
			"Weight search requires one to four actual score terms and complete counted windows and parent metadata");
		models=models_;
		terms=terms_;
		counts=counts_;
		parents=parents_;
		acceptable=acceptable_;
		strata=strata_;
		sides=sides_;
		bins=bins_;
		require(counts.length==parents.length && counts.length==acceptable.length
			&& strata.length==sides.length && strata.length==bins.length,
			"Window and parent arrays must retain their original aligned populations");
		for(int p=0; p<strata.length; p++){
			require(strata[p]>=0 && strata[p]<3 && bins[p]>=0 && bins[p]<5
				&& (strata[p]!=0 || sides[p]==0 || sides[p]==1), "Every not-seen parent needs its frozen assembly-group half");
		}
		active=new boolean[counts.length];
		for(int w=0; w<counts.length; w++){
			require(counts[w]!=null && counts[w].length==terms*models && parents[w]!=null
				&& acceptable[w]!=null && parents[w].length==acceptable[w].length,
				"Each counted window needs all score terms per model and one acceptable bitset per parent edge");
			for(int e=0; e<parents[w].length; e++){
				final int p=parents[w][e];
				require(p>=0 && p<strata.length && acceptable[w][e]!=null
					&& acceptable[w][e].length==(models+63)/64, "Window edges must bind valid original parents and complete model bitsets");
				if((models&63)!=0){
					require((acceptable[w][e][models>>>6]>>> (models&63))==0, "Acceptable bits cannot name models outside the original roster");
				}
				active[w]|=strata[p]==0;
			}
		}
		scores=new double[models];
		ranks=new int[strata.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Seed and bounded refinement are predeclared; the comparator sees fit only.
	 * Fit -1 is a final all-not-seen refit, never an independent judgment. */
	public Fit fit(final int side){
		require(side>=-1 && side<=1, "Fit side is A=0, B=1 or all-data refit=-1");
		final int metric=side<0 ? 2 : side;
		final ArrayList<Trial> trials=new ArrayList<Trial>(initial());
		require(trials.get(0).metrics[metric].parents>0, "Cannot fit an empty not-seen population");
		final Comparator<Trial> comparator=(a, b)->compare(a, b, metric);
		for(int stage=0; stage<STEPS.length; stage++){
			final ArrayList<Trial> ordered=new ArrayList<Trial>(trials);
			Collections.sort(ordered, comparator);
			final int centers=Math.min(5, ordered.size());
			for(int center=0; center<centers; center++){
				final Trial best=ordered.get(center);
				refine(trials, best.weights, best.weights.clone(), 0, STEPS[stage], "REFINE"+(stage+1));
			}
		}
		Trial winner=trials.get(0);
		for(final Trial trial:trials){if(comparator.compare(trial, winner)<0){winner=trial;}}
		final double[] range=new double[2*terms];
		for(int term=0; term<terms; term++){range[2*term]=Double.POSITIVE_INFINITY;}
		int near=0;
		for(final Trial trial:trials){
			if(near(trial.metrics[metric].cappedSum, winner.metrics[metric].cappedSum)){
				near++;
				final double[] normalized=trial.normalized();
				for(int term=0; term<terms; term++){
					range[2*term]=Math.min(range[2*term], normalized[term]);
					range[2*term+1]=Math.max(range[2*term+1], normalized[term]);
				}
			}
		}
		assert(near>0) : "The optimum must lie in its own evaluated one-percent neighborhood";
		return new Fit(side, trials, winner, range, near,
			near(trials.get(0).metrics[metric].cappedSum, winner.metrics[metric].cappedSum));
	}

	/** The unit and seeded random points have the same weights in every fit.
	 * Their retained metrics already cover both halves and all not-seen parents.
	 * Fit-local histories borrow these read-only results; refinements stay local. */
	private ArrayList<Trial> initial(){
		if(initialTrials!=null){return initialTrials;}
		final ArrayList<Trial> prepared=new ArrayList<Trial>(terms==1 ? 1 : randomDraws+1);
		final double[] unit=new double[terms];
		Arrays.fill(unit, 1);
		append(prepared, unit, "UNIT");
		//Membership only: output order remains the original random draw order.
		//These arrays are never mutated after insertion; refinement uses copies.
		final TreeSet<double[]> seen=new TreeSet<double[]>(NcrnaPartitionWeightSearch::compareWeights);
		seen.add(unit);
		final Random random=new Random(SEED);
		for(int i=0; terms>1 && i<randomDraws; i++){
			final double[] weights=unit.clone();
			for(int term=0; term<terms; term++){
				if(term!=referenceTerm){weights[term]=StrictMath.exp(LOG_MIN+LOG_SPAN*random.nextDouble());}
			}
			if(seen.add(weights)){append(prepared, weights, "RANDOM");}
		}
		assert(prepared.size()>0 && (terms==1 || prepared.size()>1) && prepared.size()<=randomDraws+1)
			: "The fixed seed supplies a bounded shared initial trial population";
		initialTrials=prepared;
		return initialTrials;
	}

	/** Reports both halves, but fit() never uses judge metrics for selection.
	 * Missing ranks cost the cap and remain explicitly unavailable. */
	public Trial evaluate(final double shortWeight, final double consensusWeight){
		require(terms==3 && referenceTerm==1, "The two-free-weight overload is only for the historical three-term geometry");
		return evaluate(new double[]{shortWeight, 1, consensusWeight});
	}

	public Trial evaluate(final double[] weights){
		return evaluate(weights, false);
	}

	private Trial evaluate(final double[] weights, final boolean allParents){
		require(weights!=null && weights.length==terms && weights[referenceTerm]==1,
			"Weights must bind every actual score term with the declared scale reference fixed to one");
		for(final double weight:weights){
			require(Double.isFinite(weight) && weight>=MIN && weight<=MAX,
				"Relative partition weights must be finite within the declared [.1,10] domain");
		}
		Arrays.fill(ranks, -1);
		for(int w=0; w<counts.length; w++){
			if(!allParents && !active[w]){continue;}
			final char[] row=counts[w];
			for(int m=0; m<models; m++){
				final int offset=terms*m;
				double sum=0;
				for(int term=0; term<terms; term++){sum+=weights[term]*row[offset+term];}
				scores[m]=sum;
			}
			for(int e=0; e<parents[w].length; e++){
				final int p=parents[w][e];
				if(!allParents && strata[p]!=0){continue;}
				final int rank=firstRank(scores, acceptable[w][e]);
				if(rank>0 && (ranks[p]<0 || rank<ranks[p])){ranks[p]=rank;}
			}
		}
		final Metrics[] metrics={new Metrics(), new Metrics(), new Metrics()};
		for(int p=0; p<ranks.length; p++){
			if(strata[p]!=0){continue;}
			metrics[sides[p]].add(ranks[p], bins[p]==0);
			metrics[2].add(ranks[p], bins[p]==0);
		}
		return new Trial(weights, metrics);
	}

	/** Unit-weight parity may inspect all original parents, not just fitted loci. */
	public int[] ranks(final double shortWeight, final double consensusWeight){
		require(terms==3 && referenceTerm==1, "Two free weights identify only the historical three-term geometry");
		evaluate(new double[]{shortWeight, 1, consensusWeight}, true);
		return ranks.clone();
	}
	public int[] ranks(final double[] weights){
		evaluate(weights, true);
		return ranks.clone();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Complete descending-score / ascending-model-ordinal ranking, including zeros. */
	public static int firstRank(final double[] values, final long[] allowed){
		require(values!=null && values.length>0 && allowed!=null && allowed.length==(values.length+63)/64,
			"Ranking requires a complete score vector and original acceptable-model bits");
		int selected=-1;
		for(int m=0; m<values.length; m++){
			assert(Double.isFinite(values[m]) && values[m]>=0) : "Weighted counts must remain finite and nonnegative";
			if((allowed[m>>>6]&(1L<<(m&63)))!=0 && (selected<0 || values[m]>values[selected])){selected=m;}
		}
		if(selected<0){return -1;}
		int rank=1;
		for(int m=0; m<values.length; m++){
			if(values[m]>values[selected] || values[m]==values[selected] && m<selected){rank++;}
		}
		assert(rank>0 && rank<=values.length) : "First acceptable model rank must lie within the original model population";
		return rank;
	}

	private void refine(final ArrayList<Trial> trials, final double[] center, final double[] point,
			final int term, final double step, final String stage){
		assert(term>=0 && term<=terms) : "Refinement visits each finite-dimensional coordinate exactly once";
		if(term==terms){
			add(trials, point, stage);
			return;
		}
		if(term==referenceTerm){
			refine(trials, center, point, term+1, step, stage);
			return;
		}
		for(int offset=-1; offset<=1; offset++){
			point[term]=clamp(center[term]*StrictMath.exp(offset*step));
			refine(trials, center, point, term+1, step, stage);
		}
	}
	private void add(final ArrayList<Trial> trials, final double[] weights, final String stage){
		assert(trials!=null && stage!=null) : "Every predeclared trial requires its fit-local history and phase";
		for(final Trial trial:trials){if(Arrays.equals(trial.weights, weights)){return;}}
		append(trials, weights, stage);
	}

	private void append(final ArrayList<Trial> trials, final double[] weights, final String stage){
		assert(trials!=null && stage!=null) : "A unique point retains its original history order and search phase";
		final Trial trial=evaluate(weights);
		trial.stage=stage;
		trials.add(trial);
	}

	/** Lexical membership is exact for the finite positive weights accepted by evaluate. */
	private static int compareWeights(final double[] a, final double[] b){
		assert(a.length==b.length) : "One initial population has a single declared score dimension";
		for(int term=0; term<a.length; term++){
			final int order=Double.compare(a[term], b[term]);
			if(order!=0){return order;}
		}
		return 0;
	}

	private static int compare(final Trial a, final Trial b, final int metric){
		assert(metric>=0 && metric<=2) : "Selection may inspect only its declared fit population";
		int order=Long.compare(a.metrics[metric].cappedSum, b.metrics[metric].cappedSum);
		if(order==0){order=Double.compare(distance(a), distance(b));}
		for(int term=0; order==0 && term<a.weights.length; term++){order=Double.compare(a.weights[term], b.weights[term]);}
		return order;
	}
	private static double distance(final Trial t){
		double sum=0;
		for(final double weight:t.weights){
			final double value=StrictMath.log(weight);
			sum+=value*value;
		}
		return sum;
	}
	private static double clamp(final double v){return Math.max(MIN, Math.min(MAX, v));}
	private static boolean near(final long value, final long best){return value<=best+best/100;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Metrics{
		private void add(final int rank, final boolean low){
			assert(rank==-1 || rank>0) : "An unavailable rank is exactly -1, never zero";
			parents++;
			final int cost=rank<0 ? CAP : Math.min(rank, CAP);
			cappedSum+=cost;
			if(rank<0){unavailable++;}
			else{
				availableRankSum+=rank;
				if(rank<=10){within10++;}
				if(rank<=50){within50++;}
			}
			if(low){
				lowParents++;
				lowCappedSum+=cost;
				if(rank>0){
					lowAvailableRankSum+=rank;
					if(rank<=10){lowWithin10++;}
					if(rank<=50){lowWithin50++;}
				}
			}
		}
		public long cappedSum, availableRankSum, lowCappedSum, lowAvailableRankSum;
		public int parents, unavailable, within10, within50, lowParents, lowWithin10, lowWithin50;
	}
	/** Read-only result, including its borrowed weights/metrics arrays. */
	public static final class Trial{
		private Trial(final double[] w, final Metrics[] m){
			weights=w.clone();
			shortWeight=w[0];
			consensusWeight=w.length>2 ? w[2] : w[w.length-1];
			metrics=m;
		}
		public double[] normalized(){
			final double[] result=weights.clone();
			double sum=0;
			for(final double weight:weights){sum+=weight;}
			for(int i=0; i<result.length; i++){result[i]/=sum;}
			return result;
		}
		public final double[] weights;
		public final double shortWeight, consensusWeight;
		public final Metrics[] metrics;
		public String stage;
	}
	public static final class Fit{
		private Fit(final int s, final ArrayList<Trial> t, final Trial w, final double[] r, final int n, final boolean u){
			side=s;
			trials=t;
			winner=w;
			range=r;
			nearTrials=n;
			unitNear=u;
		}
		public final int side, nearTrials;
		public final ArrayList<Trial> trials;
		public final Trial winner;
		/** Per-term normalized min/max over evaluated fit trials within 1%. */
		public final double[] range;
		public final boolean unitNear;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int models, terms;
	public final int randomDraws;
	public final int referenceTerm;
	private final char[][] counts;
	private final int[][] parents;
	private final long[][][] acceptable;
	private final int[] strata, sides, bins;
	private final boolean[] active;
	private final double[] scores;
	private final int[] ranks;
	private ArrayList<Trial> initialTrials;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int CAP=384, DRAWS=200;
	public static final long SEED=20261113;
	private static final double MIN=.1, MAX=10, LOG_MIN=StrictMath.log(MIN), LOG_SPAN=StrictMath.log(MAX/MIN);
	private static final double[] STEPS={StrictMath.log(2)/2, StrictMath.log(2)/10};
}
