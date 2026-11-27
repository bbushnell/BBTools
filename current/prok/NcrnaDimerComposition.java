package prok;

import java.util.Arrays;

/** Reusable, strand-oriented nucleotide composition for RNA ranking features.
 * Bin order is AA,AC,AG,AT,CA,...,TT. U is T; ambiguity breaks adjacency.
 * Missing observations have explicit presence, never an invented zero similarity.
 * @author Brian, Ganyu
 */
public final class NcrnaDimerComposition{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Scan one half-open interval. A preceding base outside the interval never
	 * contributes a pair. Counts and totals replace the previous observation. */
	public NcrnaDimerComposition scan(final byte[] bases, final int from, final int to){
		require(bases!=null && from>=0 && to>=from && to<=bases.length,
			"Composition requires a valid half-open interval in the actual oriented sequence");
		clear();
		int previous=-1;
		for(int i=from; i<to; i++){
			if((bases[i]|32)=='u'){sourceUracils++;}
			final int base=code(bases[i]);
			if(base>=0){
				knownBases++;
				if(base==1 || base==2){gcBases++;}
				if(previous>=0){
					counts[(previous<<2)|base]++;
					pairs++;
				}
			}
			previous=base;
		}
		assert(pairs<=Math.max(0, knownBases-1) && gcBases<=knownBases)
			: "Only adjacent known bases form dimers; GC is a subset of the known-base denominator";
		return this;
	}

	/** Reverse-complement the 16 bins without scanning or copying a contig again.
	 * The destination is separate worker-owned storage; counts are not pooled. */
	public void reverseComplementTo(final NcrnaDimerComposition destination){
		require(destination!=null && destination!=this,
			"Reverse-complement composition requires separate reusable destination storage");
		destination.clear();
		for(int pair=0; pair<BINS; pair++){
			final int reverse=((3-(pair&3))<<2)|(3-(pair>>>2));
			destination.counts[reverse]=counts[pair];
		}
		destination.pairs=pairs;
		destination.knownBases=knownBases;
		destination.gcBases=gcBases;
		destination.sourceUracils=sourceUracils;
	}

	/** Writes relative dimer frequencies. False means unavailable; the cleared
	 * buffer is only scratch and must not be serialized as a measured vector. */
	public boolean frequencies(final double[] destination){
		require(destination!=null && destination.length==BINS,
			"Every nucleotide dimer vector has exactly sixteen ordered coordinates");
		if(pairs==0){
			Arrays.fill(destination, 0);
			return false;
		}
		final double scale=1.0/pairs;
		for(int pair=0; pair<BINS; pair++){destination[pair]=counts[pair]*scale;}
		return true;
	}

	/** Cosine of nonnegative finite frequency/mean-frequency vectors. NaN is an
	 * in-memory missing value when either norm is zero, not a measured cosine.
	 * Callers must preserve presence and write NA for an unavailable feature. */
	public static double cosine(final double[] first, final double[] second){
		require(first!=null && second!=null && first.length==BINS && second.length==BINS,
			"Cosine compares complete sixteen-coordinate nucleotide vectors");
		double dot=0, a=0, b=0;
		for(int i=0; i<BINS; i++){
			final double x=first[i], y=second[i];
			require(Double.isFinite(x) && Double.isFinite(y) && x>=0 && y>=0 && x<=1 && y<=1,
				"Frequency and equal-member mean coordinates must be finite fractions, not raw counts or missing sentinels");
			dot+=x*y;
			a+=x*x;
			b+=y*y;
		}
		if(a==0 || b==0){return Double.NaN;}
		final double result=dot/Math.sqrt(a*b);
		require(Double.isFinite(result) && result<=1+1e-12,
			"Nonnegative finite vector cosine cannot exceed unity beyond rounding error");
		return Math.min(1, result);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/

	public boolean hasGc(){return knownBases>0;}
	public boolean hasDimers(){return pairs>0;}
	public double gcFraction(){return knownBases==0 ? Double.NaN : gcBases/(double)knownBases;}
	public long knownBases(){return knownBases;}
	public long gcBases(){return gcBases;}
	public long pairs(){return pairs;}
	/** Literal U count in the original scanned interval, retained as source audit
	 * metadata through a bin permutation; not a count in an emitted RC sequence. */
	public long sourceUracils(){return sourceUracils;}
	public long count(final int pair){
		require(pair>=0 && pair<BINS, "Dimer indices use the complete ordered AA..TT axis");
		return counts[pair];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private void clear(){
		Arrays.fill(counts, 0);
		knownBases=0;
		gcBases=0;
		pairs=0;
		sourceUracils=0;
	}
	private static int code(final byte base){
		switch(base|32){
			case 'a': return 0;
			case 'c': return 1;
			case 'g': return 2;
			case 't':
			case 'u': return 3;
			default: return -1;
		}
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final long[] counts=new long[BINS];
	private long knownBases, gcBases, pairs, sourceUracils;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int BINS=16;
}
