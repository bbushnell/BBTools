package dna;
import java.util.Arrays;


/**
 * Composite motif that represents multiple alternative motifs using OR logic.
 * Returns positive matches if any of its constituent sub-motifs match the target sequence.
 * All sub-motifs must have the same length and center position.
 * @author Brian Bushnell
 */
public class MotifMulti extends Motif {
	
	/**
	 * Constructs a composite motif from multiple alternative motifs.
	 * All motifs must have the same length and center position.
	 * @param name_ Name identifier for this composite motif
	 * @param args Variable number of constituent motifs to combine
	 */
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
	
	@Override
	public float normalize(double strength){
		return (float)strength;
//		throw new RuntimeException("MotifMulti can't normalize without knowing the submotif.");
	}
	
	
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

	@Override
	public int numBases() {
		return sub[0].numBases();
	}
	
	/** Array of constituent motifs that comprise this composite motif */
	public final Motif[] sub;
	
}
