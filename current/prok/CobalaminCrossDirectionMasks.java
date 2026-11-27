package prok;

/** Six emitted-to-scoring directions, excluding the three own-target diagonals.
 * Bits are A->B, A->C, B->A, B->C, C->A, C->B, where A/B/C are the native
 * main/second/third member axes. Own support is never removed by a cross mask.
 * This class owns no resources and never changes its caller's arrays.
 * @author Brian Bushnell, Diona
 */
public final class CobalaminCrossDirectionMasks{

	private CobalaminCrossDirectionMasks(){}

	/*--------------------------------------------------------------*/
	/*----------------          Direction Masks     ----------------*/
	/*--------------------------------------------------------------*/

	/** One-time CLI parser: null/all enables every cross direction; none disables
	 * all of them. Otherwise use distinct RFAM:RFAM pairs separated by commas. */
	public static int parse(final String text){
		if(text==null || text.equals("all")){return 63;}
		if(text.equals("none")){return 0;}
		int allowed=0;
		for(final String pair:text.split(",", -1)){
			final int colon=pair.indexOf(':');
			require(colon>0 && colon<pair.length()-1 && pair.indexOf(':', colon+1)<0,
				"Cross directions require all, none, or comma-separated RFAM:RFAM pairs: "+pair);
			final int emitted=member(pair.substring(0, colon));
			final int target=member(pair.substring(colon+1));
			final int bit=1<<direction(emitted, target);
			require((allowed&bit)==0, "A cross direction must not occur twice: "+pair);
			allowed|=bit;
		}
		assert(allowed>0 && allowed<=63) : "An explicit nonempty pair list must enable at least one of the six cross directions";
		return allowed;
	}
	/** Returns three target bits for this emitted row, always excluding own.
	 * Resources can cache these bits once and use one test per live target. */
	public static int targets(final int allowed, final int emitted){
		require(allowed>=0 && allowed<=63 && emitted>=0 && emitted<3,
			"Target bits require six allowed-direction bits and an emitted member in [0,2]");
		final int row=(allowed>>>(emitted*2))&3;
		final int lowerMask=(1<<emitted)-1;
		final int result=(row&lowerMask)|((row>>>emitted)<<(emitted+1));
		assert((result&~7)==0 && (result&(1<<emitted))==0)
			: "Cached cross-target bits must preserve the separate unconditional own-scoring path";
		return result;
	}
	/** Returns a bit index from zero through five, not a shifted bit value. */
	public static int direction(final int emitted, final int target){
		require(emitted>=0 && emitted<3 && target>=0 && target<3 && emitted!=target,
			"A cross direction needs two distinct native member axes in [0,2]");
		return emitted*2+target-(target>emitted ? 1 : 0);
	}
	/** Retains any own pass, or a pass on at least one allowed cross direction.
	 * passTargets uses target bits 0..2; allowed uses direction bits 0..5. */
	public static boolean retained(final int passTargets, final int emitted, final int allowed){
		require(passTargets>=0 && passTargets<=7 && emitted>=0 && emitted<3 && allowed>=0 && allowed<=63,
			"Retention needs three target-pass bits, a valid emitted member and six allowed-direction bits");
		if((passTargets&(1<<emitted))!=0){return true;}
		// Remove the own-target bit. The other two pass bits now align with
		// this emitted member's consecutive pair in the six-direction word.
		final int lowerMask=(1<<emitted)-1;
		final int crossPass=(passTargets&lowerMask)|((passTargets>>>(emitted+1))<<emitted);
		assert((crossPass&~3)==0) : "Deleting the own bit must leave exactly the two columns of one emitted-member row";
		return (crossPass&(allowed>>>(emitted*2)))!=0;
	}
	/** Identity at the ceiling bypasses CM tests. Below it, own support survives
	 * regardless of the cross-fallback floor; only cross support is floor-gated.
	 * Equal bounds describe an empty cross-only interval, not an invalid mask. */
	public static boolean retained(final int passTargets, final int emitted, final int allowed,
			final float identity, final float fallbackFloor, final float ceiling){
		// Validate the discrete axes even when identity alone would accept.
		final boolean supported=retained(passTargets, emitted, allowed);
		require(Float.isFinite(identity) && identity>=0 && identity<=1 && Float.isFinite(fallbackFloor)
			&& Float.isFinite(ceiling) && fallbackFloor>=0 && fallbackFloor<=ceiling && ceiling<=1,
			"Identity and ordered fallback/ceiling bounds must be finite in [0,1]");
		return identity>=ceiling || (passTargets&(1<<emitted))!=0 || identity>=fallbackFloor && supported;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Cost Kernel      ----------------*/
	/*--------------------------------------------------------------*/

	/** Predicts cross work after an own failure. Skip own/disallowed targets,
	 * pay ranking once only if a target is visited, and stop on the first finite
	 * score-minus-own-GA pass. NaN is unavailable; negative infinity is impossible;
	 * both fail. These supplied service clocks are not measurements made here. */
	public static double crossCpuNanos(final float[] scores, final float[] ga, final double[] cutoffs,
			final double[] service, final double ranking, final int emitted, final int[] priority, final int allowed){
		require(scores!=null && ga!=null && cutoffs!=null && service!=null && scores.length==3
			&& ga.length==3 && cutoffs.length==3 && service.length==3,
			"Cross costs require three matching score, GA, cutoff and service axes");
		require(emitted>=0 && emitted<3 && allowed>=0 && allowed<=63,
			"Cross costs require an emitted member in [0,2] and an allowed mask in [0,63]");
		require(Double.isFinite(ranking) && ranking>=0, "Ranking cost must be finite and nonnegative");
		for(int target=0; target<3; target++){
			require(validScore(scores[target]) && Float.isFinite(ga[target]) && Double.isFinite(cutoffs[target]),
				"Scores may be finite, unavailable or impossible; GA and cutoffs must be finite");
			require(Double.isFinite(service[target]) && service[target]>=0, "Target service costs must be finite and nonnegative");
		}
		require(priority!=null && priority.length==3, "Cross priority must contain all three native members");
		int seen=0;
		for(final int target:priority){
			require(target>=0 && target<3 && (seen&(1<<target))==0, "Cross priority must be a permutation of native members [0,2]");
			seen|=1<<target;
		}
		assert(seen==7) : "Every native member must be represented before own/disallowed directions are filtered";
		require(!passes(scores[emitted], ga[emitted], cutoffs[emitted]),
			"Cross cost prediction requires the own target to have failed already");
		double cpu=0;
		boolean visited=false;
		for(final int target:priority){
			if(target==emitted || (allowed&(1<<direction(emitted, target)))==0){continue;}
			if(!visited){
				cpu=ranking;
				visited=true;
			}
			cpu+=service[target];
			require(Double.isFinite(cpu), "Summed cross service and ranking costs must not overflow to an unusable prediction");
			if(passes(scores[target], ga[target], cutoffs[target])){break;}
		}
		return cpu;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private static int member(final String name){
		for(int i=0; i<MEMBERS.length; i++){if(MEMBERS[i].equals(name)){return i;}}
		throw new IllegalArgumentException("Unknown cross-direction RFAM member: "+name);
	}
	private static boolean validScore(final float score){
		return Float.isFinite(score) || Float.isNaN(score) || score==Float.NEGATIVE_INFINITY;
	}
	/** Matches CobalaminCmAcceptance: promote before subtraction, preserving the
	 * native finite score-minus-GA comparison at the cutoff boundary. */
	private static boolean passes(final float score, final float ga, final double cutoff){
		return Float.isFinite(score) && (double)score-ga>=cutoff;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
}
