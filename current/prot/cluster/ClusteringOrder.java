package prot.cluster;

import prot.ProteinSequence;

/** Reorders an already shuffled eligible-input permutation without changing input ordinals.
 * Stable ordering retains the seeded permutation for equal hydrophobicity or bins.
 * Undefined hydrophobicity sorts after defined values, using explicit availability.
 * @author Keqing
 */
public final class ClusteringOrder {

	public enum Mode {SHUFFLE, HYDRO, HYDROBIN}

	/** Uses native residue hydrophobicities and primitive scratch allocated only for sorted modes. */
	public static void apply(int[] order,ProteinSequence[] sequences,Mode mode){
		if(order==null || sequences==null || mode==null){throw new IllegalArgumentException("Ordering needs inputs, a permutation and a mode");}
		if(mode==Mode.SHUFFLE || order.length<2){return;}
		final double[] keys=new double[sequences.length];
		final boolean[] defined=new boolean[sequences.length];
		final CandidateFilter.Context context=new CandidateFilter.Context();
		for(int input : order){
			if(input<0 || input>=sequences.length || sequences[input]==null){
				throw new IllegalArgumentException("Eligible input ordinal is missing or out of bounds: "+input);
			}
			context.reset(sequences[input].enc); defined[input]=context.defined();
			if(defined[input]){
				final double hydro=context.hydrophobicity();
				keys[input]=mode==Mode.HYDRO ? hydro : Math.floor(hydro/BIN_WIDTH);
				assert(Double.isFinite(keys[input])) : "Native finite hydrophobicities and fixed positive bin width must produce finite ordering keys";
			}
		}
		// Stable bottom-up merge over primitive ordinals: no per-protein wrapper objects.
		int[] source=order,target=new int[order.length];
		for(int width=1; width<order.length; width=width>order.length/2 ? order.length : width<<1){
			for(int left=0; left<order.length; ){
				final int middle=left+Math.min(width,order.length-left);
				final int right=middle+Math.min(width,order.length-middle);
				int a=left,b=middle,out=left;
				while(a<middle && b<right){
					if(compare(source[a],source[b],keys,defined)<=0){target[out++]=source[a++];}
					else{target[out++]=source[b++];}
				}
				while(a<middle){target[out++]=source[a++];}
				while(b<right){target[out++]=source[b++];}
				assert(out==right) : "Every merge must preserve the complete eligible-input permutation";
				left=right;
			}
			final int[] swap=source; source=target; target=swap;
		}
		if(source!=order){System.arraycopy(source,0,order,0,order.length);}
	}

	private static int compare(int a,int b,double[] keys,boolean[] defined){
		if(defined[a]!=defined[b]){return defined[a] ? -1 : 1;}
		return defined[a] ? Double.compare(keys[a],keys[b]) : 0;
	}

	public static final double BIN_WIDTH=0.05;
}
