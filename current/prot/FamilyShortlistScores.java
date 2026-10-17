package prot;

/** Global F4 normalization and shortlist selection for independently stored packs.
 * Inputs must follow a stable pack order and each pack's local active-index order.
 * Raw composition scores must use the same background and transform. Callers
 * validate that reference compatibility before using this arithmetic.
 * Matches ProteinSearcher.scoreFamiliesF4's kmernorm=none formula; normalizing
 * separately inside each pack would change the meaning of both feature weights.
 * All output and input arrays belong to the caller; no query scratch is allocated.
 * @author Keqing
 */
public final class FamilyShortlistScores {

	/** Computes one query's moments across every family, then scores each pack.
	 * Outputs must be distinct from all composition input arrays.
	 */
	public static void score(final int[][] hits, final float[][] composition,
			final double[][] medianLength, final double[][] sdLogLength,
			final int queryLength, final float[][] out){
		require(queryLength>0, "Global F4 requires a nonempty query for its log-length term");
		final int count=validate(hits, composition, medianLength, sdLogLength, out);
		double sumK=0, squareK=0, sumC=0, squareC=0;
		int finiteC=0, nanC=0;
		for(int p=0; p<hits.length; p++){
			for(int r=0; r<hits[p].length; r++){
				final int k=hits[p][r]; final float c=composition[p][r];
				if(k<0 || Float.isInfinite(c)){
					throw new IllegalArgumentException("Invalid raw shortlist feature at pack "+p+", rank "+r+": hits="+k+", composition="+c);
				}
				sumK+=k; squareK+=(double)k*k;
				if(Float.isNaN(c)){nanC++;}
				else{finiteC++; sumC+=c; squareC+=(double)c*c;}
			}
		}
		if(nanC>0 && finiteC>0){throw new IllegalArgumentException("Mixed finite and NaN raw composition scores");}
		final double meanK=sumK/count, meanC=(finiteC>0 ? sumC/finiteC : Double.NaN);
		final double sdK=Math.sqrt(Math.max(0, squareK/count-meanK*meanK));
		final double sdC=Math.sqrt(Math.max(0, (finiteC>0 ? squareC/finiteC-meanC*meanC : Double.NaN)));
		final double logQuery=Math.log(queryLength);
		for(int p=0; p<hits.length; p++){
			for(int r=0; r<hits[p].length; r++){
				final double length=medianLength[p][r], spread=sdLogLength[p][r];
				if(!(length>0) || !Double.isFinite(length) || spread<0 || !Double.isFinite(spread)){
					throw new IllegalArgumentException("Invalid shortlist length statistics at pack "+p+", rank "+r);
				}
				final double distance=Math.abs(logQuery-Math.log(length))/Math.max(spread, 0.05);
				final double zk=sdK>1e-12 ? (hits[p][r]-meanK)/sdK : 0;
				final double zc=sdC>1e-12 ? (composition[p][r]-meanC)/sdC : 0;
				out[p][r]=(float)(zk+2*zc-distance);
				if(!Float.isFinite(out[p][r])){
					throw new IllegalArgumentException("Nonfinite global F4 at pack "+p+", rank "+r);
				}
			}
		}
	}

	/** Selects descending F4, breaking ties by canonical pack order then local rank.
	 * Both output arrays have the requested shortlist length; every slot is filled.
	 */
	public static void top(final float[][] scores, final int[] packs, final int[] ranks){
		require(scores!=null && scores.length>0 && packs!=null && ranks!=null
			&& packs!=ranks && packs.length>0 && packs.length==ranks.length,
			"Global shortlist needs scores and distinct, equally sized pack/rank outputs");
		long count=0;
		for(float[] row : scores){require(row!=null, "Missing pack score array"); count+=row.length;}
		require(count>=packs.length, "Requested shortlist exceeds the combined family count");
		int size=0;
		for(int p=0; p<scores.length; p++){
			for(int r=0; r<scores[p].length; r++){
				final float score=scores[p][r];
				if(!Float.isFinite(score)){throw new IllegalArgumentException("Nonfinite F4 at pack "+p+", rank "+r);}
				if(size==packs.length && score<=scores[packs[size-1]][ranks[size-1]]){continue;}
				int at=Math.min(size, packs.length-1);
				while(at>0 && score>scores[packs[at-1]][ranks[at-1]]){at--;}
				if(size<packs.length){size++;}
				for(int i=size-1; i>at; i--){packs[i]=packs[i-1]; ranks[i]=ranks[i-1];}
				packs[at]=p; ranks[at]=r;
			}
		}
		assert(size==packs.length) : "Validated global family count must fill every shortlist slot";
	}

	private static int validate(final int[][] hits, final float[][] composition,
			final double[][] medianLength, final double[][] sdLogLength, final float[][] out){
		require(hits!=null && hits.length>0 && composition!=null && medianLength!=null
			&& sdLogLength!=null && out!=null, "Global F4 requires all pack feature arrays");
		require(composition.length==hits.length && medianLength.length==hits.length
			&& sdLogLength.length==hits.length && out.length==hits.length,
			"Global F4 feature pack counts differ");
		long count=0;
		for(int p=0; p<hits.length; p++){
			require(hits[p]!=null && composition[p]!=null && medianLength[p]!=null
				&& sdLogLength[p]!=null && out[p]!=null, "Missing global F4 feature row");
			final int n=hits[p].length;
			if(n<1 || composition[p].length!=n || medianLength[p].length!=n
					|| sdLogLength[p].length!=n || out[p].length!=n){
				throw new IllegalArgumentException("Global F4 family dimensions differ or are empty at pack "+p);
			}
			for(int q=0; q<hits.length; q++){
				require(out[p]!=composition[q], "F4 outputs must not overwrite raw composition inputs needed for later packs");
				if(q<p){require(out[p]!=out[q], "Each pack requires a distinct F4 output array");}
			}
			count+=n;
		}
		require(count<=Integer.MAX_VALUE, "Combined family count exceeds the supported integer index range");
		return (int)count;
	}

	private static void require(final boolean ok, final String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private FamilyShortlistScores(){}
}
