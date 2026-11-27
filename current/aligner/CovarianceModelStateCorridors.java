package aligner;

/** Immutable per-consensus residual bounds for CM state endpoint placements.
 * Four signed offsets adjust raw glocal start-low/high and end-low/high anchors.
 * Unobserved states retain the declared uniform radius. An inverted projected
 * interval is empty and requests fallback; it is never repaired by sorting.
 * The caller owns and reuses the four-int projection buffer.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelStateCorridors{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CovarianceModelStateCorridors(final CovarianceModel model_, final short[] offsets_, final boolean[] observed_){
		require(model_!=null && offsets_!=null && observed_!=null && observed_.length==model_.states()
			&& (long)observed_.length*4==offsets_.length,
			"Learned state corridors require four signed offsets and one support flag per state of their exact CM");
		model=model_;
		offsets=offsets_.clone();
		observed=observed_.clone();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Projection          ----------------*/
	/*--------------------------------------------------------------*/

	/** Projects one raw state anchor in place. The result may be an empty interval;
	 * false tells the scorer to represent no cells, rather than widen it silently.
	 * State support is joint across its four bounds, never inferred from zeroes.
	 */
	public boolean project(final int state, final int length, final int[] bounds, final int uniformRadius){
		require(state>=0 && state<observed.length && length>=0 && bounds!=null && bounds.length>=4 && uniformRadius>=0,
			"State projection needs an existing state, actual window length and reusable four-bound buffer");
		require(bounds[0]>=0 && bounds[0]<=bounds[1] && bounds[1]<=length
			&& bounds[2]>=0 && bounds[2]<=bounds[3] && bounds[3]<=length,
			"Residual bounds apply to raw ordered glocal intervals, before any uniform expansion");
		return projectOrdered(state, length, bounds, uniformRadius);
	}

	/** Raw anchors translated from a different context may lie outside the search
	 * window. Add residuals before intersecting that window: clipping an anchor
	 * first would change its signed learned offset or invent an edge cell.
	 */
	public boolean projectUnclipped(final int state, final int length, final int[] bounds, final int uniformRadius){
		require(state>=0 && state<observed.length && length>=0 && bounds!=null && bounds.length>=4 && uniformRadius>=0,
			"Context projection needs an existing state, actual search length and reusable four-bound buffer");
		require(bounds[0]<=bounds[1] && bounds[2]<=bounds[3],
			"Translated raw anchors may cross a search edge but must retain their original interval order");
		return projectOrdered(state, length, bounds, uniformRadius);
	}

	/** Both public contracts validate their raw coordinate domain before this fill. */
	private boolean projectOrdered(final int state, final int length, final int[] bounds, final int uniformRadius){
		assert(bounds[0]<=bounds[1] && bounds[2]<=bounds[3]) :
			"Signed residuals adjust distinct ordered raw endpoints before the actual search-domain intersection";
		final int base=state*4;
		final boolean learned=observed[state];
		final long startLow=Math.max(0L, (long)bounds[0]+(learned ? offsets[base] : -uniformRadius));
		final long startHigh=Math.min(length, (long)bounds[1]+(learned ? offsets[base+1] : uniformRadius));
		final long endLow=Math.max(0L, (long)bounds[2]+(learned ? offsets[base+2] : -uniformRadius));
		final long endHigh=Math.min(length, (long)bounds[3]+(learned ? offsets[base+3] : uniformRadius));
		if(startLow>startHigh || endLow>endHigh){
			// Clamping BOTH sides to an edge would invent a cell when the whole
			// learned interval lies outside the window. Preserve an empty set.
			bounds[0]=1;
			bounds[1]=0;
			bounds[2]=1;
			bounds[3]=0;
			return false;
		}
		bounds[0]=(int)startLow;
		bounds[1]=(int)startHigh;
		bounds[2]=(int)endLow;
		bounds[3]=(int)endHigh;
		return true;
	}

	/** Allocation-free membership of the same triangular/state-length domain used
	 * by placement layouts. Call project first, even when it returns false; the
	 * resulting inverted bounds then reject every cell naturally.
	 */
	public boolean contains(final int state, final int length, final int[] projected, final int start, final int end){
		require(state>=0 && state<observed.length && length>=0 && projected!=null && projected.length>=4,
			"Containment uses this model's state and the caller's already projected endpoint bounds");
		return start>=0 && start<=end && end<=length && start>=projected[0] && start<=projected[1]
			&& end>=projected[2] && end<=projected[3] && (model.type[state]!=CovarianceModel.E || start==end);
	}

	/** Whether this state was actually visited by retained training trajectories. */
	public boolean observed(final int state){
		require(state>=0 && state<observed.length, "Support lookup must index this learned resource's exact state axis");
		return observed[state];
	}

	/** One signed endpoint offset; only observed states have learned offsets. */
	public int offset(final int state, final int bound){
		require(state>=0 && state<observed.length && observed[state] && bound>=0 && bound<4,
			"Only observed states have four learned signed endpoint offsets");
		return offsets[state*4+bound];
	}

	public boolean belongsTo(final CovarianceModel expected){return model==expected;}
	public int states(){return observed.length;}

	/** Encloses the observed residuals of one member's supported profiles. The
	 * supplied bounds already carry their selected margin; never expand it again.
	 * This pool supplies no per-consensus support and does not alter any input.
	 */
	public static CovarianceModelStateCorridors pool(final CovarianceModel model,
			final CovarianceModelStateCorridors[] profiles){
		require(model!=null && profiles!=null && profiles.length>0,
			"Member pooling requires an explicit CM and its complete existing consensus axis");
		final short[] bounds=new short[Math.multiplyExact(model.states(), 4)];
		final boolean[] seen=new boolean[model.states()];
		boolean any=false;
		for(final CovarianceModelStateCorridors profile:profiles){
			if(profile==null){continue;}
			require(profile.belongsTo(model), "Residual envelopes cannot pool different member CMs, even with equal state dimensions");
			for(int state=0; state<seen.length; state++){
				if(!profile.observed[state]){continue;}
				final int base=state*4;
				for(int bound=0; bound<4; bound++){
					final short value=profile.offsets[base+bound];
					if(!seen[state] || ((bound&1)==0 ? value<bounds[base+bound] : value>bounds[base+bound])){
						bounds[base+bound]=value;
					}
				}
				seen[state]=true;
				any=true;
			}
		}
		return any ? new CovarianceModelStateCorridors(model, bounds, seen) : null;
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final CovarianceModel model;
	private final short[] offsets;
	private final boolean[] observed;
}
