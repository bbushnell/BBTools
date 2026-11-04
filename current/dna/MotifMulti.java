package dna;

import java.util.Arrays;

/**
 * A composite {@link Motif} that ORs several sub-motifs: it matches (exactly or extended)
 * at a position if ANY constituent sub-motif matches there, and reports the maximum
 * per-sub normalized {@link #matchStrength}, with zero as the lower bound. All sub-motifs are assumed to share the
 * length and center of sub[0]. Because it cannot know which sub produced a strength,
 * {@link #normalize} is an identity cast. The child array is borrowed and mutable;
 * thread safety also depends on every child's matching implementation.
 * @author Brian Bushnell
 */
public class MotifMulti extends Motif{

	/** Borrows the child array. The caller must supply a nonempty array of nonnull children
	 * with compatible length/center geometry; these requirements are not explicitly validated.
	 * Geometry is copied only from the first child; later array edits do not refresh it. */
	public MotifMulti(String name_, Motif...args){
		super(name_, args[0].length, args[0].center);
		commonLetters=Arrays.toString(args);
		sub=args;
	}

	/** Returns at the first exact child match at center a; child exceptions propagate. */
	@Override
	public boolean matchesExactly(byte[] source, int a){
		for(int i=0; i<sub.length; i++){
			Motif m=sub[i];
			if(m.matchesExactly(source, a)){
				return true;
			}
		}
		return false;
	}

	/** Returns at the first extended child match at center a; child exceptions propagate. */
	@Override
	public boolean matchesExtended(byte[] source, int a){
		for(int i=0; i<sub.length; i++){
			Motif m=sub[i];
			if(m.matchesExtended(source, a)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns the strength value unchanged as normalization is not applicable for composite motifs.
	 * Normalization requires knowledge of the specific sub-motif that generated the strength value.
	 * @param strength The strength value to normalize
	 * @return The strength value cast to float without modification
	 */
	@Override
	public float normalize(double strength){
		return (float)strength;
//		throw new RuntimeException("MotifMulti can't normalize without knowing the submotif.");
	}

	/**
	 * Calculates the maximum normalized match strength among all sub-motifs at the specified position.
	 * Starts from zero and applies the local comparison helper to each normalized score.
	 * A negative score cannot lower a finite maximum; NaN behavior follows that helper's
	 * ordinary comparisons rather than a special NaN policy.
	 *
	 * @param source The sequence to evaluate
	 * @param a Candidate center coordinate
	 * @return Maximum normalized strength with a zero floor when all child scores are finite
	 */
	@Override
	public float matchStrength(byte[] source, int a){
		float max=0;
		for(int i=0; i<sub.length; i++){
			Motif m=sub[i];
			float temp=m.matchStrength(source, a);
			temp=m.normalize(temp);
			max=max(max, temp);
		}
		return max;
	}

	/**
	 * Delegates the model-size convention to the first child.
	 * @return First child's numBases value; its unit depends on that child type
	 */
	@Override
	public int numBases(){
		return sub[0].numBases();
	}

	/** Borrowed mutable child array; it and its elements must remain nonnull. */
	public final Motif[] sub;

}
