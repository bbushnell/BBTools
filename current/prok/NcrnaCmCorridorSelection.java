package prok;

import map.IntHashMap;
import structures.IntList;
import structures.LongList;

/** Selects whole training trajectories for one consensus and one source fold.
 * Input call ordinals are strictly increasing; negative departure means the
 * call has no usable finite anchored path. At most floor(N*trimPercent/100) greatest
 * departures are removed, with lower call ordinal first on ties. Trimming
 * stops before either retained support minimum would be broken.
 * Scratch is primitive, reusable and confined to one fitting thread.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmCorridorSelection{

	public NcrnaCmCorridorSelection(final int minCalls_, final int minGroups_, final int folds_){
		this(minCalls_, minGroups_, folds_, 1);
	}

	public NcrnaCmCorridorSelection(final int minCalls_, final int minGroups_, final int folds_, final int trimPercent_){
		require(minCalls_>0 && minGroups_>0 && folds_>1 && trimPercent_>=0 && trimPercent_<100,
			"Whole-path fitting requires positive retained-call/group support, source-group folds and trimpercent=0..99");
		minCalls=minCalls_;
		minGroups=minGroups_;
		folds=folds_;
		trimPercent=trimPercent_;
	}

	/** Sets one previously unused bit per retained call, preserving other fits.
	 * excludedFold=-1 fits all training sources; otherwise group modulo folds
	 * excludes an entire resolved assembly/alias group before support or trimming.
	 * Returns retained call support; zero selects uniform fallback for this fit.
	 */
	public int select(final IntList calls, final int[] callGroups, final int[] departures,
			final int excludedFold, final byte[] retained, final int bit){
		require(calls!=null && callGroups!=null && departures!=null && retained!=null
			&& callGroups.length==departures.length && retained.length==departures.length
			&& excludedFold>=-1 && excludedFold<folds && bit>0 && bit<=128 && (bit&(bit-1))==0,
			"Each fit needs one call axis, one valid source fold and one distinct retained-mask bit");
		groups.clear();
		ranked.clear();
		int previous=-1;
		for(int i=0; i<calls.size; i++){
			final int call=calls.get(i);
			require(call>previous && call<retained.length && callGroups[call]>=0 && (retained[call]&bit)==0,
				"Model call ordinals must be unique/increasing and each fit may assign its mask bit only once");
			previous=call;
			if(departures[call]<0 || excludedFold>=0 && callGroups[call]%folds==excludedFold){continue;}
			groups.increment(callGroups[call]);
			// Both packed fields are nonnegative, so signed long ordering is the
			// declared departure order; complementing the ordinal resolves ties.
			ranked.add(((long)departures[call]<<32)|(Integer.MAX_VALUE-call));
		}
		int count=ranked.size, distinct=groups.size();
		if(count<minCalls || distinct<minGroups){return 0;}
		for(int i=0; i<ranked.size; i++){
			final int call=Integer.MAX_VALUE-(int)ranked.get(i);
			retained[call]|=(byte)bit;
		}
		ranked.sort();
		final int budget=(int)((long)count*trimPercent/100);
		for(int dropped=0; dropped<budget; dropped++){
			final int call=Integer.MAX_VALUE-(int)ranked.get(ranked.size-1-dropped);
			final int group=callGroups[call], multiplicity=groups.get(group);
			assert(multiplicity>0) : "Every trimming candidate was counted in this fit's resolved-group support";
			if(count-1<minCalls || distinct-(multiplicity==1 ? 1 : 0)<minGroups){break;}
			retained[call]&=(byte)~bit;
			groups.increment(group, -1);
			count--;
			if(multiplicity==1){distinct--;}
		}
		assert(count>=minCalls && distinct>=minGroups && ranked.size-count<=budget) :
			"The learned corridor's retained trajectories must preserve both support minima and the declared whole-path budget";
		return count;
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	private final int minCalls, minGroups, folds, trimPercent;
	private final IntHashMap groups=new IntHashMap();
	private final LongList ranked=new LongList();
}
