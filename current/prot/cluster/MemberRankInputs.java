package prot.cluster;

import prot.ProteinSequence;

/** Immutable ordinal-bound length/KD cache for one clustering input population.
 * Computes native member KD once, before workers start. Retains only primitive
 * arrays (13 bytes per input before array headers), plus a borrowed population
 * reference for identity checks, not another peptide copy. Callers keep that
 * population unchanged while using the cache.
 * Unknown KD has an explicit validity bit; requesting it as measured data fails.
 * @author Brian Bushnell, Keqing
 */
public final class MemberRankInputs{

	public MemberRankInputs(final ProteinSequence[] sequences){
		this(sequences,DEFAULT_HYDRO_WEIGHT);
	}

	/** Coefficient is frozen with the population binding before any worker starts. */
	public MemberRankInputs(final ProteinSequence[] sequences,final double hydroWeight_){
		hydroWeight=checkedHydroWeight(hydroWeight_);
		if(sequences==null){throw new IllegalArgumentException("Member ranking inputs require the original ordered protein population");}
		population=sequences;
		lengths=new int[sequences.length]; means=new double[sequences.length]; defined=new boolean[sequences.length];
		final CandidateFilter.Context scalar=new CandidateFilter.Context();
		for(int i=0; i<sequences.length; i++){
			final ProteinSequence sequence=sequences[i];
			if(sequence==null || sequence.length()<1){
				throw new IllegalArgumentException("Missing original protein at member ordinal "+i);
			}
			lengths[i]=sequence.length();
			scalar.reset(sequence.enc);
			defined[i]=scalar.defined(); means[i]=scalar.hydrophobicity();
			assert(Double.isFinite(means[i]) && (!defined[i] || means[i]>=-4.5 && means[i]<=4.5)) :
				"CandidateFilter.Context must retain finite native KD and explicit unavailable state";
		}
	}

	/** Reuses immutable observations when initial seed selection/hydration used the default.
	 * Bind the returned view before attaching families; existing bindings remain unchanged.
	 */
	public MemberRankInputs withHydroWeight(final double value){
		final double weight=checkedHydroWeight(value);
		return weight==hydroWeight ? this : new MemberRankInputs(this,weight);
	}

	private MemberRankInputs(final MemberRankInputs source,final double weight){
		assert(source!=null) : "Coefficient views must reuse the original frozen observations, not recalculate member KD";
		population=source.population;lengths=source.lengths;means=source.means;defined=source.defined;
		hydroWeight=checkedHydroWeight(weight);
	}

	/** Supported experiment range: KD contributes at most one shared-hit score unit. */
	public static double checkedHydroWeight(final double value){
		final boolean valid=Double.isFinite(value) && value>=0 && value<=1;
		assert(valid) : "memberhydroweight must lie in [0,1]: nonnegative KD similarity preserves the zero lower bound; the experiment caps its contribution at one hit; got "+value;
		if(!valid){throw new IllegalArgumentException("memberhydroweight must be finite and in [0,1]: "+value);}
		return value==0 ? 0 : value;
	}

	public static final double DEFAULT_HYDRO_WEIGHT=0.2;
	public final double hydroWeight;

	public int size(){return lengths.length;}
	/** Borrowed population identity prevents a same-size unrelated cache from binding to an engine. */
	public boolean belongsTo(final ProteinSequence[] sequences){return population==sequences;}
	public int length(final int ordinal){checkOrdinal(ordinal); return lengths[ordinal];}
	public boolean defined(final int ordinal){checkOrdinal(ordinal); return defined[ordinal];}

	/** Missing measurements cannot become neutral KD or consensus-derived replacements. */
	public double meanKd(final int ordinal){
		checkOrdinal(ordinal);
		if(!defined[ordinal]){
			throw new IllegalArgumentException("No measured member KD at input ordinal "+ordinal+"; explicit unavailable handling is required");
		}
		return means[ordinal];
	}

	/** Called only by the family prefix publisher; this cache does not decide membership ownership. */
	public void addTo(final MemberRankStatistics statistics, final int ordinal){
		if(statistics==null){throw new IllegalArgumentException("Member publication requires its owning statistics accumulator");}
		final double kd=meanKd(ordinal);
		statistics.add(lengths[ordinal],kd);
	}

	private void checkOrdinal(final int ordinal){
		if(ordinal<0 || ordinal>=lengths.length){
			throw new IllegalArgumentException("Member ordinal "+ordinal+" is outside the cached input population of "+lengths.length);
		}
	}
	private final int[] lengths;
	private final ProteinSequence[] population;
	private final double[] means;
	private final boolean[] defined;
}
