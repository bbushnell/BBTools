package prot.cluster;

/** Shared local split/fusion gain; this is not calibrated positive retention.
 * Evidence is immutable and must describe the same frozen comparison context.
 * No sequences, negative panels, mutable family state or hidden thresholds are read.
 * @author Keqing
 */
public final class LocalFamilyObjective{

	public LocalFamilyObjective(final Parameters parameters_){
		if(parameters_==null){throw new IllegalArgumentException("Local objective needs explicit parameters");}
		parameters=parameters_;
	}

	/** Computes member-weighted quality/separation/size gain minus one budget-slot cost.
	 * Parent and children must describe the same conserved original-member population.
	 * Callers obtain all three evidence objects under identical scoring and neighbor-context rules.
	 * slotCost is the current budget-th family size, frozen for the whole structural pass.
	 */
	public double splitGain(final Evidence parent,final Evidence left,final Evidence right,final double temperature,final int slotCost){
		if(parent==null || left==null || right==null || (long)left.members+right.members!=parent.members){
			throw new IllegalArgumentException("A local split comparison requires one parent and two children conserving member count");
		}
		if(!Double.isFinite(temperature) || temperature<0 || slotCost<0){
			throw new IllegalArgumentException("Local objective temperature and frozen budget-slot cost must be nonnegative and finite");
		}
		final double a=left.members/(double)parent.members,b=right.members/(double)parent.members;
		final double cohesion=a*left.cohesion+b*right.cohesion-parent.cohesion;
		final double separation=a*(left.cohesion-left.neighbor)+b*(right.cohesion-right.neighbor)-
			(parent.cohesion-parent.neighbor);
		final double balance=1-Math.abs((long)left.members-right.members)/(double)parent.members;
		assert(balance>0 && balance<=1) : "Two positive conserved child counts define a balance in (0,1]";
		final double pressure=temperature*parameters.sizeForm.pressure(parent.members,parameters.largeBoundary)*balance;
		final double gain=parent.members*(parameters.cohesionScale*cohesion+
			parameters.separationScale*separation+parameters.sizeScale*pressure)-slotCost;
		if(!Double.isFinite(gain)){
			throw new IllegalArgumentException("Local objective overflow: reduce temperature or objective scales");
		}
		return gain;
	}

	/** Exact inverse gain for the same evidence, temperature, frozen slot cost and parameters.
	 * Regenerating a different parent model or changing evidence/context is a different proposal.
	 */
	public double fusionGain(final Evidence left,final Evidence right,final Evidence merged,final double temperature,final int slotCost){
		return -splitGain(merged,left,right,temperature,slotCost);
	}

	/** Both operation directions use this single strict positive-margin rule. */
	public boolean accepts(final double gain){
		if(!Double.isFinite(gain)){throw new IllegalArgumentException("Acceptance requires a finite measured local gain");}
		return gain>parameters.margin;
	}

	/** Complete normalized evidence; missing measurements must not be replaced by zero.
	 * A known absence of any eligible neighbor may use zero only under the caller's explicit context rule.
	 */
	public static final class Evidence{
		public Evidence(final int members_,final double cohesion_,final double neighbor_){
			if(members_<1 || !unit(cohesion_) || !unit(neighbor_)){
				throw new IllegalArgumentException("Local evidence needs positive member count and finite normalized quality in [0,1]");
			}
			members=members_; cohesion=cohesion_; neighbor=neighbor_;
		}
		public final int members;
		public final double cohesion,neighbor;
	}

	/** Explicit study parameters; constructing this object does not enable any program mode. */
	public static final class Parameters{
		public Parameters(final int largeBoundary_,final double margin_,final double cohesionScale_,
				final double separationScale_,final double sizeScale_,final SizeForm sizeForm_){
			if(largeBoundary_<1 || !Double.isFinite(margin_) || margin_<=0 ||
				!nonnegative(cohesionScale_) || !nonnegative(separationScale_) || !nonnegative(sizeScale_) || sizeForm_==null){
				throw new IllegalArgumentException("Local objective needs positive boundary/margin, finite nonnegative scales and a size form");
			}
			largeBoundary=largeBoundary_; margin=margin_; cohesionScale=cohesionScale_;
			separationScale=separationScale_; sizeScale=sizeScale_; sizeForm=sizeForm_;
		}
		public final int largeBoundary;
		public final double margin,cohesionScale,separationScale,sizeScale;
		public final SizeForm sizeForm;
	}

	/** Dimensionless pressure above an explicit member-count boundary; zero below it. */
	public enum SizeForm{
		LOG2,LINEAR,QUADRATIC;

		public double pressure(final int members,final int boundary){
			if(members<1 || boundary<1){throw new IllegalArgumentException("Size pressure requires positive members and boundary");}
			if(members<=boundary){return 0;}
			final double ratio=members/(double)boundary;
			switch(this){
				case LOG2: return Math.log(ratio)/LOG_TWO;
				case LINEAR: return ratio-1;
				case QUADRATIC: return ratio*ratio-1;
				default: throw new AssertionError("Every declared size-pressure form must have an explicit formula");
			}
		}
	}

	private static boolean unit(final double value){return Double.isFinite(value) && value>=0 && value<=1;}
	private static boolean nonnegative(final double value){return Double.isFinite(value) && value>=0;}
	public final Parameters parameters;
	private static final double LOG_TWO=Math.log(2);
}
