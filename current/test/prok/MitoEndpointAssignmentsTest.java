package test.prok;


import java.util.Arrays;
import java.util.Random;
import idaligner.QuantumAligner;
import prok.TrnaKmerIndex;

/** Checks bounded native-index assignment, identity selection, ties and no-hit state.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointAssignmentsTest {
	public static void main(String[] args){
		final byte[] gene=new byte[256],alphabet={'A','C','G','T'};final Random random=new Random(20261024);
		for(int i=0;i<gene.length;i++){gene[i]=alphabet[random.nextInt(4)];}
		final byte[][] models=new byte[10][];
		for(int m=0;m<models.length;m++){models[m]=gene.clone();if(m<8){for(int i=m;i<gene.length;i+=29){models[m][i]=models[m][i]=='A'?(byte)'C':(byte)'A';}}}
		final TrnaKmerIndex index=new TrnaKmerIndex(models,9,false,1,0,0,1);
		final MitoEndpointAssignments.Selection s=MitoEndpointAssignments.assign(gene,models,index);
		check(s.tried==8,"Construction mapping must align only top8 even when10models share9mers");
		check(s.model==8 && s.identity==1f,"Maximum identity must win; tied exact models retain shared-index rank order");
		final byte[] none=new byte[256];Arrays.fill(none,(byte)'N');
		final MitoEndpointAssignments.Selection absent=MitoEndpointAssignments.assign(none,models,index);
		check(absent.model==-1 && absent.tried==0 && Float.isNaN(absent.identity),"No shared9mer remains explicit, never a fabricated assigned model");
		final byte[][] shifted=new byte[8][];final float[] scores=new float[8];final int[][] positions=new int[8][2];
		for(int i=0;i<shifted.length;i++){
			shifted[i]=new byte[gene.length+96];for(int j=0;j<shifted[i].length;j++){shifted[i][j]=alphabet[random.nextInt(4)];}
			System.arraycopy(gene,0,shifted[i],7+9*i,gene.length);
			for(int j=19;j<gene.length;j+=31){shifted[i][7+9*i+j]=(byte)'N';}
			scores[i]=QuantumAligner.alignStatic(gene,shifted[i],positions[i]);
		}
		final int[] reused={123456,-77};
		for(int i=shifted.length-1;i>=0;i--){
			final float value=MitoEndpointAssignments.alignIndependent(gene,shifted[i],reused);
			check(value==scores[i] && Arrays.equals(reused,positions[i]),"Dirty hints and reversed work order must reproduce fresh native reference coordinates and identities");
		}
		System.out.println("PASS MitoEndpointAssignmentsTest native_k9_top8_max_identity_stable_ties_no_hit_independent_hints");
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
