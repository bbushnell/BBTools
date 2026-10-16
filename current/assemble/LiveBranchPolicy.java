package assemble;

/**
 * Origin-independent screening and features for conservative live branch rescue.
 * This class does not alter the graph or decide whether a read-derived path is
 * true. Native accepts bypass inference; only grey-zone junctions need a model.
 * The caller must separately enforce actual-predecessor, seed and ownership
 * guards. No field depends on the partial contig or the thread reaching a node.
 * @author Fischl
 */
final class LiveBranchPolicy {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Utility class; all scratch belongs to the calling worker. */
	private LiveBranchPolicy(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Screens a node using both immutable neighborhoods, not one directed edge.
	 * A true return from the native branch predicate is supplied separately for
	 * each side. Swapping sides and their flags must leave the result unchanged.
	 * The initial half-depth competitor screen is experimental, not a default
	 * change to Tadpole. Equality at the screen boundary remains a hard stop.
	 */
	static int level(final int count, final int[] left, final int[] right,
			final boolean leftBranch, final boolean rightBranch,
			final int minExtend, final double maxCompetitorFraction){
		validateCounts(left); validateCounts(right);
		if(count<0 || minExtend<1 || !Double.isFinite(maxCompetitorFraction) ||
				maxCompetitorFraction<=0 || maxCompetitorFraction>1){
			throw new IllegalArgumentException("Branch screening requires nonnegative depth, positive mce and a fraction in (0,1].");
		}
		// Do not retroactively veto the native heuristic's accepted decisions.
		if(!leftBranch && !rightBranch){return NATIVE_ACCEPT;}
		if(count<minExtend || hardSide(left, minExtend, maxCompetitorFraction) ||
				hardSide(right, minExtend, maxCompetitorFraction)){return HARD_STOP;}
		return GREY_ZONE;
	}

	/** Missing preferred edges, equal maxima and strong competitors cannot be rescued. */
	private static boolean hardSide(final int[] counts, final int minExtend, final double fraction){
		assert(counts.length==4 && minExtend>0) : "Each directed DNA neighborhood has four counts; mce protects the chosen edge.";
		int highest=0, second=0;
		for(int count : counts){
			if(count>highest){second=highest; highest=count;}
			else if(count>second){second=count;}
		}
		return highest<minExtend || second==highest || second>=highest*fraction;
	}

	/**
	 * Returns the only permitted preferred base, or -1 for a tie/absent side.
	 * Sorting features never changes these raw nucleotide identities. Arriving
	 * through any other predecessor must remain a hidden-branch stop in the caller.
	 */
	static int preferredBase(final int[] counts){
		validateCounts(counts);
		int highest=0, best=-1;
		boolean tied=false;
		for(int i=0; i<4; i++){
			if(counts[i]>highest){highest=counts[i]; best=i; tied=false;}
			else if(counts[i]==highest){tied=true;}
		}
		return tied ? -1 : best;
	}

	/**
	 * Writes current depth, two ordered sorted side blocks, K and global means.
	 * The raw arrays are never changed. Side order is chosen AFTER normalization,
	 * so equal float blocks are identical even if large integer depths round alike.
	 * This is a new 12-input schema, not the 27-input completed-contig postpass.
	 */
	static void fill(final float[] out, final int k, final int count,
			final int[] left, final int[] right, final double uniqueMean,
			final double observationMean){
		validateCounts(left); validateCounts(right);
		if(out==null || out.length!=INPUTS || k<2 || count<0){
			throw new IllegalArgumentException("Live branch features require 12 inputs, K>=2 and nonnegative node depth.");
		}
		out[0]=depth(count);
		for(int i=0; i<4; i++){
			out[1+i]=depth(left[i]); out[5+i]=depth(right[i]);
		}
		sortFour(out, 1); sortFour(out, 5);
		for(int i=0; i<4; i++){
			if(out[1+i]<out[5+i]){break;}
			if(out[1+i]>out[5+i]){
				for(int j=0; j<4; j++){
					final float temp=out[1+j]; out[1+j]=out[5+j]; out[5+j]=temp;
				}
				break;
			}
		}
		out[9]=(float)(0.125*Math.log(k)*INV_LOG_2);
		out[10]=depth(uniqueMean); out[11]=depth(observationMean);
	}

	/** Log1p keeps zero distinct from singleton evidence and does not clip high depth. */
	private static float depth(final double value){
		if(!Double.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid live branch depth: "+value);}
		return (float)(0.125*Math.log1p(value)*INV_LOG_2);
	}

	/** Descending in-place sort of exactly one four-base feature block. */
	private static void sortFour(final float[] values, final int offset){
		assert(offset==1 || offset==5) : "Only the two four-base blocks may be sorted; node/global features keep fixed positions.";
		for(int i=offset+1; i<offset+4; i++){
			final float value=values[i];
			int j=i;
			while(j>offset && values[j-1]<value){values[j]=values[j-1]; j--;}
			values[j]=value;
		}
	}

	/** Rejects malformed evidence instead of silently converting corruption to a stop. */
	private static void validateCounts(final int[] counts){
		if(counts==null || counts.length!=4){throw new IllegalArgumentException("Expected four base-indexed extension counts.");}
		for(int count : counts){if(count<0){throw new IllegalArgumentException("Negative extension depth: "+count);}}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int NATIVE_ACCEPT=0, HARD_STOP=1, GREY_ZONE=2;
	static final int INPUTS=12;
	static final String VERSION="livebranch12-v1";
	private static final double INV_LOG_2=1/Math.log(2);
}
