package dna;
import java.util.Arrays;


public class MotifMulti extends Motif {
	
	public MotifMulti(String name_, Motif...args){
		super(name_, args[0].length, args[0].center);
		commonLetters=Arrays.toString(args);
		sub=args;
	}
	
	
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
	 * Evaluates each constituent motif and returns the highest normalized strength value.
	 *
	 * @param source The sequence to evaluate
	 * @param a Starting position in the source sequence
	 * @return Maximum normalized match strength among all sub-motifs
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
	 * Returns the number of bases used by the motifs.
	 * Since all sub-motifs have the same length, returns the base count of the first sub-motif.
	 * @return Number of bases in the motif patterns
	 */
	@Override
	public int numBases() {
		return sub[0].numBases();
	}
	
	public final Motif[] sub;
	
}
