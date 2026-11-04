package test.dna;

import dna.Gene;

/** Literal genome-build parsing checks; does not load reference or chain resources.
 * Both assertion modes must preserve numeric suffixes: callers perform build remapping.
 * @author Shinobu
 */
public final class GeneBuildTest{

	/** Runs ordinary controls or the two-character hg-prefix regression. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Require case=controls|hg");}
		final String[] inputs;
		final int[] expected;
		if(args[0].equals("case=controls")){
			inputs=new String[]{"36", "build36", "BUILD=37", "b=36", "B37", "=36"};
			expected=new int[]{36, 36, 37, 36, 37, 36};
		}else if(args[0].equals("case=hg")){
			inputs=new String[]{"hg18", "HG19", "hG=38"};
			expected=new int[]{18, 19, 38};
		}else{throw new IllegalArgumentException(args[0]);}
		int failures=0;
		for(int i=0; i<inputs.length; i++){
			try{
				final int actual=Gene.toBuild(inputs[i]);
				System.out.println(inputs[i]+" => "+actual+" expected="+expected[i]);
				if(actual!=expected[i]){failures++;}
			}catch(RuntimeException | AssertionError e){
				failures++;
				System.out.println(inputs[i]+" => "+e.getClass().getSimpleName()+" expected="+expected[i]);
			}
		}
		if(failures!=0){throw new AssertionError("Gene.toBuild must remove complete prefixes and retain numeric suffixes: failures="+failures);}
		System.out.println("PASS "+args[0]);
	}
}
