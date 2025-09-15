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
	
	public final Motif[] sub;
	
}
