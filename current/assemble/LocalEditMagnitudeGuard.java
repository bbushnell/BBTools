package assemble;

/** Stateless magnitude-based abstention predicate over supplied local-edit depths.
 * Invoke only after resolving candidate-family ambiguity; this helper does not
 * perform candidate discovery, resolve ambiguity or edit a read. Factor0 disables it.
 * Negative-one probe slots are unavailable/unchanged, distinct from measured zero.
 * Arrays are borrowed for the call and never modified; keep them stable while used.
 * @author Fischl */
final class LocalEditMagnitudeGuard{
	/** Vetoes only when all eight actual edits are measured, seven have depth at least
	 * two below original, and exactly one has depth at least original times factor.
	 * Factor0 or original depth below2 returns false before validating probe arrays.
	 * Otherwise both arrays must have four slots and all supplied depths must be at
	 * least -1. The unchanged substitution is unavailable, leaving at most eight
	 * actual edits; that upper bound is assertion-checked. Missing evidence abstains
	 * from this veto, so false does not mean a candidate is safe to accept.
	 * @param original Nonnegative original-window depth
	 * @param factor Zero to disable, otherwise an integer at least2
	 * @param substitutions Four A/C/G/T substitution depths; unchanged base has -1
	 * @param deletion Deletion depth or -1 when unavailable
	 * @param insertions Four A/C/G/T insertion depths, with -1 for unavailable probes
	 */
	static boolean veto(final int original, final int factor, final int[] substitutions,
		final int deletion, final int[] insertions){
		if(original<0 || factor<0 || factor==1){throw new IllegalArgumentException("Depth is nonnegative; magnitude factor is0(disabled) or>=2.");}
		if(factor==0 || original<2){return false;}
		if(substitutions==null || insertions==null || substitutions.length!=4 || insertions.length!=4){throw new IllegalArgumentException("Magnitude guard requires4substitution slots and4insertion slots.");}
		int available=0, low=0, high=0;
		for(int slot=0; slot<9; slot++){
			final int depth=slot<4 ? substitutions[slot] : slot==4 ? deletion : insertions[slot-5];
			if(depth<-1){throw new IllegalArgumentException("Probe depth below unavailable sentinel.");}
			if(depth<0){continue;}
			available++;
			if((long)original-depth>1){low++;}
			if(depth/(long)original>=factor){high++;}
		}
		assert(available<=8) : "LocalEditKmerProbe excludes unchanged substitution, leaving at most3+1+4actual edits.";
		return available==8 && low==7 && high==1;
	}
}
