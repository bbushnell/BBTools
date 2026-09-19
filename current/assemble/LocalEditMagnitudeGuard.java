package assemble;

/** Conservative abstention predicate used by Tadpole fixindels.
 * Invoke only after resolving candidate-family ambiguity. Factor0 is disabled.
 * Negative-one probe slots are unavailable/unchanged, never zero evidence.
 * @author Fischl */
final class LocalEditMagnitudeGuard{
	static boolean veto(final int original,final int factor,final int[] substitutions,
		final int deletion,final int[] insertions){
		if(original<0 || factor<0 || factor==1){throw new IllegalArgumentException("Depth is nonnegative; magnitude factor is0(disabled) or>=2.");}
		if(factor==0 || original<2){return false;}
		if(substitutions==null || insertions==null || substitutions.length!=4 || insertions.length!=4){throw new IllegalArgumentException("Magnitude guard requires4substitution slots and4insertion slots.");}
		int available=0,low=0,high=0;
		for(int slot=0;slot<9;slot++){
			final int depth=slot<4?substitutions[slot]:slot==4?deletion:insertions[slot-5];
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
