package test.prot;

import prot.HbmTailRobustness;

import java.util.Arrays;

/** Boundary fixtures for short-query perturbations and unchanged long-query spans. @author Keqing */
public final class HbmTailRobustnessTest{
	public static void main(String[] args){
		int tested=0;
		for(int length : new int[]{30, 40, 41, 500}){
			final byte[] original=new byte[length];
			for(int i=0; i<length; i++){original[i]=(byte)(i%20);}
			final byte[] saved=original.clone();
			final int[] lengths={length, length, length, length, length, length-10, length-20,
				Math.min(500, length+10), Math.min(500, length+20)};
			for(int v=0; v<9; v++){
				final byte[] out=HbmTailRobustness.variant(original, v);
				check(Arrays.equals(original, saved), "A perturbation changed the source query");
				check(out!=original && out.length==lengths[v], "Wrong independent output length for variant "+v+" at length "+length);
				if(v==0){check(Arrays.equals(out, original), "Original variant changed sequence bytes");}
				else if(v<=3){
					final int changed=Math.min(length, v==1 ? 10 : v==2 ? 20 : 40);
					for(int i=0; i<length; i++){
						check((out[i]!=original[i])==(i>=length-changed), "Replacement changed the wrong residue span");
					}
				}else if(v==4){
					final int[] before=new int[20], after=new int[20];
					for(int i=0; i<length; i++){
						before[original[i]]++; after[out[i]]++;
						if(i<length-20){check(out[i]==original[i], "Shuffle changed the preserved prefix");}
					}
					check(Arrays.equals(before, after), "Shuffle changed residue composition");
				}else if(v<=6){
					for(int i=0; i<out.length; i++){check(out[i]==original[i], "Trimming changed retained residues");}
				}else{
					final int added=v==7 ? 10 : 20, dropped=Math.max(0, length+added-500);
					for(int i=0; i<length-dropped; i++){check(out[i]==original[i+dropped], "Appending moved the retained suffix incorrectly");}
				}
				for(byte aa : out){check(aa>=0 && aa<20, "Perturbation introduced an invalid residue");}
				tested++;
			}
		}
		System.err.println("HBM_TAIL_ROBUSTNESS_TEST_PASS cases="+tested+" short_queries=true legacy_spans=true source_immutable=true");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private HbmTailRobustnessTest(){}
}
