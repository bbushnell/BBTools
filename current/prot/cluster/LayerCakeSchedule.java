package prot.cluster;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

/** Temperature and stopping control for completed structural passes.
 * Fission implementations finish all recursive waves before returning; fusion
 * implementations publish exactly one frozen set of disjoint reciprocal pairs.
 * This driver never reads negatives or performs alignments itself.
 * @author Keqing
 */
public final class LayerCakeSchedule{

	public LayerCakeSchedule(final Settings settings_){
		if(settings_==null){throw new IllegalArgumentException("Layer-cake scheduling needs explicit settings");}
		settings=settings_;
	}

	/** Runs opt-in structural work. Operations must not mutate their frozen inputs.
	 * A positive maxCycles is an explicit study bound, reported as CYCLE_LIMIT
	 * rather than convergence; zero imposes no cycle cap.
	 */
	public Result run(final Family[] initial,final Operations operations,final Listener listener){
		if(initial==null || settings.mode!=Mode.NONE && operations==null){
			throw new IllegalArgumentException("Structural scheduling requires families and the enabled pass implementation");
		}
		Family[] current=initial.clone();
		Summary summary=summarize(current,settings.familyBudget,settings.inputTails);
		final ArrayList<PhaseSummary> history=new ArrayList<PhaseSummary>();
		if(settings.mode==Mode.NONE){return new Result(current,summary,history,Stop.DISABLED,0);}
		if(summary.assignedTails==0){return new Result(current,summary,history,Stop.EMPTY,0);}
		if(reached(summary)){return new Result(current,summary,history,Stop.COVERAGE,0);}
		double temperature=settings.initialTemperature;
		for(int cycle=1; ; cycle++){
			long changes=0;
			for(Phase phase:PHASES){
				if(!settings.mode.includes(phase)){continue;}
				final long start=System.nanoTime();
				final PassResult pass=operations.run(phase,current,temperature,cycle);
				if(pass==null){throw new IllegalStateException("An enabled structural pass returned no result");}
				final Summary next=summarize(pass.families,settings.familyBudget,settings.inputTails);
				validateTransition(summary,next,phase,pass.acceptedOperations);
				final PhaseSummary event=new PhaseSummary(cycle,phase,temperature,pass.acceptedOperations,
					summary,next,System.nanoTime()-start);
				history.add(event);
				if(listener!=null){listener.completed(event);}
				changes+=pass.acceptedOperations;
				current=pass.families; summary=next;
				if(reached(summary)){return new Result(current,summary,history,Stop.COVERAGE,cycle);}
				temperature=Math.max(settings.minimumTemperature,temperature*settings.cooling);
			}
			if(changes==0){return new Result(current,summary,history,Stop.UNCHANGED,cycle);}
			if(settings.maxCycles>0 && cycle>=settings.maxCycles){
				return new Result(current,summary,history,Stop.CYCLE_LIMIT,cycle);
			}
			if(cycle==Integer.MAX_VALUE){throw new IllegalStateException("Structural cycle counter exhausted before convergence");}
		}
	}

	/** Checks population and binary-operation algebra, in addition to each operator's ordinal proof. */
	private static void validateTransition(final Summary before,final Summary after,final Phase phase,final long changes){
		final long expected=phase==Phase.FISSION ? before.families+changes : before.families-changes;
		if(after.assignedTails!=before.assignedTails || expected!=after.families ||
			phase==Phase.FUSION && changes>before.families/2){
			throw new IllegalStateException("A structural pass must conserve original tails and count binary changes; fusion pairs must be disjoint");
		}
	}

	/** Uses the original-input denominator, including any initially unassigned tails. */
	private boolean reached(final Summary summary){
		return settings.inputTails>0 && summary.topMembers>=Math.ceil(settings.coverageTarget*settings.inputTails);
	}

	/** Exact top-budget membership sum in one pass, with only budget-sized primitive scratch. */
	public static Summary summarize(final Family[] families,final int budget,final long inputTails){
		if(families==null || budget<1 || inputTails<0){throw new IllegalArgumentException("Coverage needs a roster, positive budget and nonnegative input census");}
		final int[] heap=new int[Math.min(budget,families.length)];
		int used=0;
		long assigned=0,top=0;
		for(Family family:families){
			if(family==null || family.size()<1){throw new IllegalArgumentException("Structural populations contain only nonempty families");}
			final int size=family.size(); assigned=Math.addExact(assigned,size);
			if(used<heap.length){
				int child=used++;
				while(child>0){
					final int parent=(child-1)>>>1;
					if(heap[parent]<=size){break;}
					heap[child]=heap[parent]; child=parent;
				}
				heap[child]=size; top+=size;
			}else if(size>heap[0]){
				top+=size-heap[0];
				int parent=0;
				while(parent<(used>>>1)){
					int child=(parent<<1)+1;
					if(child+1<used && heap[child+1]<heap[child]){child++;}
					if(heap[child]>=size){break;}
					heap[parent]=heap[child]; parent=child;
				}
				heap[parent]=size;
			}
		}
		if(assigned>inputTails){throw new IllegalStateException("Structural member count exceeds the original input denominator");}
		assert(top<=assigned && used==Math.min(budget,families.length)) : "The coverage heap retains exactly the largest available budget of whole families";
		return new Summary(families.length,assigned,top,used==budget ? heap[0] : 0);
	}

	/** Prices signed local moves at an interior rank of the current active family count, capped by the model budget. */
	public static int pricedRank(final int activeFamilies,final int maximumBudget,final double pricedFraction){
		if(activeFamilies<1 || maximumBudget<1 || !Double.isFinite(pricedFraction) || pricedFraction<=0 || pricedFraction>1){
			throw new IllegalArgumentException("A priced rank needs at least one active family, a positive budget and 0<fraction<=1");
		}
		final BigDecimal count=BigDecimal.valueOf(pricedFraction).multiply(BigDecimal.valueOf(activeFamilies))
			.setScale(0,RoundingMode.FLOOR);
		return Math.min(maximumBudget,Math.max(1,count.intValue()));
	}

	public interface Operations{
		PassResult run(Phase phase,Family[] frozenFamilies,double temperature,int cycle);
	}
	public interface Listener{
		void completed(PhaseSummary summary);
	}
	public enum Phase{FISSION,FUSION}
	public enum Stop{DISABLED,EMPTY,COVERAGE,UNCHANGED,CYCLE_LIMIT}
	public enum Mode{
		NONE,FISSION,FUSION,LAYERCAKE;
		boolean includes(final Phase phase){return this==LAYERCAKE || this==FISSION && phase==Phase.FISSION || this==FUSION && phase==Phase.FUSION;}
	}

	public static final class Settings{
		public Settings(final Mode mode_,final double initialTemperature_,final double minimumTemperature_,
				final double cooling_,final int familyBudget_,final double coverageTarget_,final long inputTails_,final int maxCycles_){
			if(mode_==null || !finiteNonnegative(initialTemperature_) || !finiteNonnegative(minimumTemperature_) ||
				minimumTemperature_>initialTemperature_ || !Double.isFinite(cooling_) || cooling_<=0 || cooling_>=1 ||
				familyBudget_<1 || !Double.isFinite(coverageTarget_) || coverageTarget_<=0 || coverageTarget_>1 ||
				inputTails_<0 || maxCycles_<0){
				throw new IllegalArgumentException("Invalid schedule: require finite temperatures, 0<cooling<1, positive budget/target and nonnegative census/cycle cap");
			}
			mode=mode_; initialTemperature=initialTemperature_; minimumTemperature=minimumTemperature_;
			cooling=cooling_; familyBudget=familyBudget_; coverageTarget=coverageTarget_; inputTails=inputTails_; maxCycles=maxCycles_;
		}
		public final Mode mode;
		public final double initialTemperature,minimumTemperature,cooling,coverageTarget;
		public final int familyBudget,maxCycles;
		public final long inputTails;
	}

	public static final class PassResult{
		public PassResult(final Family[] families_,final long acceptedOperations_){
			if(families_==null || acceptedOperations_<0){throw new IllegalArgumentException("A completed pass needs a population and nonnegative accepted-operation count");}
			families=families_; acceptedOperations=acceptedOperations_;
		}
		public final Family[] families;
		public final long acceptedOperations;
	}
	public static final class Summary{
		private Summary(final int families_,final long assigned_,final long top_,final int slotCost_){
			families=families_; assignedTails=assigned_; topMembers=top_; slotCost=slotCost_;
		}
		public final int families;
		/** Current budget-th largest member count; zero means the model budget has unused slots. */
		public final int slotCost;
		public final long assignedTails,topMembers;
	}
	public static final class PhaseSummary{
		private PhaseSummary(final int cycle_,final Phase phase_,final double temperature_,final long changes_,
				final Summary before_,final Summary after_,final long elapsed_){
			cycle=cycle_; phase=phase_; temperature=temperature_; acceptedOperations=changes_;
			before=before_; after=after_; elapsedNanos=elapsed_;
		}
		public final int cycle;
		public final Phase phase;
		public final double temperature;
		public final long acceptedOperations,elapsedNanos;
		public final Summary before,after;
	}
	public static final class Result{
		private Result(final Family[] families_,final Summary summary_,final ArrayList<PhaseSummary> history_,final Stop stop_,final int cycles_){
			families=families_; summary=summary_; history=history_.toArray(new PhaseSummary[history_.size()]); stop=stop_; cycles=cycles_;
		}
		public final Family[] families;
		public final Summary summary;
		public final PhaseSummary[] history;
		public final Stop stop;
		public final int cycles;
	}
	private static boolean finiteNonnegative(final double value){return Double.isFinite(value) && value>=0;}
	public final Settings settings;
	private static final Phase[] PHASES={Phase.FISSION,Phase.FUSION};
}
