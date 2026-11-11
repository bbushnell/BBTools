package assemble;

import dna.AminoAcid;
import kmer.Walker;
import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;
import ukmer.WalkerU;

/**
 * Shared, truth-blind feature calculation for same-K boundary decisions.
 * Endpoint neighborhoods retain nucleotide identities until the actual inward
 * edge has been selected. Only the normalized model representation is sorted.
 * @author Fischl
 */
public final class SameKJoinFeatures{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Fills the caller-owned 27-input vector without changing either raw neighborhood.
	 * Partner/strand reversal exchanges the endpoint blocks before their lexicographic
	 * ordering; the final sorted vector is invariant when all oriented inputs agree.
	 * Each raw neighborhood is current/right-A,C,G,T/left-A,C,G,T. The actual inward
	 * count comes from that nucleotide-indexed input, independently of sorted output.
	 * Means are nonnegative finite depths; position counts are positive window counts.
	 * @param out Exactly 27 writable entries; may be partly overwritten if validation fails
	 * @param k Word length in bases, at least two
	 * @param source Source endpoint's nine nucleotide-indexed counts
	 * @param dest Destination endpoint's nine nucleotide-indexed counts
	 * @param sourceInwardBase Source's actual inward nucleotide index, zero through three
	 * @param destInwardBase Destination's actual inward nucleotide index, zero through three
	 * @param sourceMean Mean depth across the source arm's positions
	 * @param destMean Mean depth across the destination arm's positions
	 * @param sourcePositions Number of source arm windows contributing to its mean
	 * @param destPositions Number of destination arm windows contributing to its mean
	 * @param uniqueMean Population depth averaged over positive unique keys
	 * @param observationMean Population depth weighted by count observations
	 */
	public static void fill(final float[] out, final int k, final int[] source, final int[] dest,
			final int sourceInwardBase, final int destInwardBase, final double sourceMean,
			final double destMean, final long sourcePositions, final long destPositions,
			final double uniqueMean, final double observationMean){
		if(out==null || out.length!=INPUTS || k<2){
			throw new IllegalArgumentException("Same-K features require exactly 27 inputs and K>=2.");
		}
		endpoint(out, 0, source, true, sourceInwardBase, sourceMean, sourcePositions);
		endpoint(out, BLOCK, dest, false, destInwardBase, destMean, destPositions);
		// Ordering after normalization makes ties harmless: equal blocks are identical.
		for(int i=0; i<BLOCK; i++){
			if(out[i]<out[BLOCK+i]){break;}
			if(out[i]>out[BLOCK+i]){
				for(int j=0; j<BLOCK; j++){
					final float temp=out[j]; out[j]=out[BLOCK+j]; out[BLOCK+j]=temp;
				}
				break;
			}
		}
		out[24]=(float)(0.125*Math.log(k)*INV_LOG_2);
		out[25]=depth(uniqueMean); out[26]=depth(observationMean);
	}

	/** Encodes current, outward4, inward4, actual inward, mean and position count. */
	private static void endpoint(final float[] out, final int offset, final int[] counts,
			final boolean outwardRight, final int inwardBase, final double mean, final long positions){
		if(counts==null || counts.length!=9 || inwardBase<0 || inwardBase>3 || positions<1){
			throw new IllegalArgumentException("An endpoint requires nine counts, a defined inward base and positive length.");
		}
		assert(offset==0 || offset==BLOCK) : "Two endpoint blocks precede the three global features.";
		for(int count : counts){if(count<0){throw new IllegalArgumentException("Negative neighborhood depth: "+count);}}
		final int outward=outwardRight ? 1 : 5, inward=outwardRight ? 5 : 1;
		out[offset]=depth(counts[0]);
		for(int i=0; i<4; i++){
			out[offset+1+i]=depth(counts[outward+i]);
			out[offset+5+i]=depth(counts[inward+i]);
		}
		sortFour(out, offset+1); sortFour(out, offset+5);
		out[offset+9]=depth(counts[inward+inwardBase]);
		out[offset+10]=depth(mean);
		out[offset+11]=(float)(0.0625*Math.log(positions)*INV_LOG_2);
	}

	/** Tiny in-place descending insertion sort; no per-record sorting scratch. */
	private static void sortFour(final float[] values, final int start){
		assert(start>=0 && start+4<=values.length) : "Only one four-extension block may be sorted.";
		for(int i=start+1; i<start+4; i++){
			final float value=values[i];
			int j=i;
			while(j>start && values[j-1]<value){values[j]=values[j-1]; j--;}
			values[j]=value;
		}
	}

	/** Encodes a finite nonnegative depth as log2(1+value)/8.
	 * Preserves absence versus singleton evidence without clipping high depth. */
	public static float depth(final double value){
		if(!Double.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid feature depth: "+value);}
		return (float)(0.125*Math.log1p(value)*INV_LOG_2);
	}

	/** Queries current/right4/left4 counts in the requested contig orientation.
	 * Mutates caller-owned keys and out scratch, leaves contig bases unchanged and
	 * uses the count provider only for lookups. Absent-key sentinels become zero;
	 * extension slots retain A,C,G,T order. */
	static void neighborhood(final TadpoleGraph.Counts counts, final Contig contig,
			final boolean reverse, final int start, final int k, final Kmer key,
			final Kmer neighbor, final int[] out){
		assert(out.length==9 && key.kbig==k && neighbor.kbig==k) : "Neighborhood scratch must match the live table K.";
		assert(start>=0 && start+k<=contig.length()) : "A complete endpoint Kmer must fit its contig.";
		key.clear();
		for(int i=0; i<k; i++){
			final byte base=SameKJoinCollector.baseAt(contig, reverse, start+i);
			if(AminoAcid.baseToNumber[base]<0){throw new IllegalArgumentException("Undefined base in an endpoint Kmer.");}
			key.addRight(base);
		}
		out[0]=Math.max(0, counts.count(key));
		for(int base=0; base<4; base++){
			neighbor.setFrom(key); neighbor.addRightNumeric(base);
			out[1+base]=Math.max(0, counts.count(neighbor));
			neighbor.setFrom(key); neighbor.addLeftNumeric(base);
			out[5+base]=Math.max(0, counts.count(neighbor));
		}
	}

	/** Visits the existing explicit table once; never builds another table or rereads reads. */
	static Moments moments(final Tadpole tad){
		if(!(tad instanceof Tadpole2)){return momentsByWalker(tad);}
		final Moments result=new Moments();
		final KmerTableSetU tables=((Tadpole2)tad).tables;
		for(int i=0; i<tables.ways; i++){addLongCounts(result, tables.getTable(i));}
		return result;
	}

	/**
	 * Preserves WalkerU1D's ascending main-cell then victim-list order exactly,
	 * including floating-point addition order, without reconstructing unused keys.
	 * Positive values are authoritative occupancy (HashArrayU1D.fillHistogram and
	 * WalkerU1D.next); this is intentionally not applied to the short-key tables.
	 */
	public static void addLongCounts(final Moments result, final HashArrayU1D table){
		assert(result!=null && table!=null) : "Moment accumulation requires a live explicit table and its caller-owned result.";
		for(int count : table.values()){if(count>0){result.add(count);}}
		for(KmerNodeU node : table.victims().toList()){
			final int count=node.count();
			if(count>0){result.add(count);}
		}
	}

	/** Retains the original key-aware traversal as the short-K path and diagnostic oracle. */
	static Moments momentsByWalker(final Tadpole tad){
		final Moments result=new Moments();
		if(tad instanceof Tadpole1){
			final Walker walker=((Tadpole1)tad).tables.walk();
			while(walker.next()){if(walker.value()>0){result.add(walker.value());}}
		}else if(tad instanceof Tadpole2){
			final WalkerU walker=((Tadpole2)tad).tables.walk();
			while(walker.next()){if(walker.value()>0){result.add(walker.value());}}
		}else{throw new IllegalArgumentException("Exact depth moments require a known explicit Tadpole table.");}
		// Only positive live counts belong to the population. Keep that filter explicit
		// on both walker paths, independently of each walker's empty-cell handling.
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact population, with floating-point sums to avoid integer-square overflow. */
	public static final class Moments{
		/** Public constructor for relocated tests that compare long-K table traversal. */
		public Moments(){}
		/** Adds one positive live count, not a depth-histogram bin or observation. */
		public void add(final int count){
			if(count<1){throw new IllegalArgumentException("Depth moments count only live positive entries.");}
			unique=Math.addExact(unique, 1);
			observations+=count; squares+=(double)count*count;
			assert(Double.isFinite(squares)) : "Finite exact-table counts must produce representable double moments.";
		}
		/** Sum of counts divided by positive-key count; zero for an empty population. */
		public double uniqueMean(){return unique==0 ? 0 : observations/unique;}
		/** Sum of squared counts divided by total observations; zero when empty. */
		public double observationMean(){return observations==0 ? 0 : squares/observations;}
		/** Number of positive live keys accumulated, with checked integer overflow. */
		public long unique;
		/** Count sum and squared-count sum, accumulated as doubles in traversal order. */
		public double observations, squares;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Two 12-column endpoint blocks followed by K and two population means. */
	static final int INPUTS=27, BLOCK=12;
	/** Fixed exported/inference schema identity; normalization changes require migration. */
	public static final String VERSION="samek27-v1";
	/** Converts natural logarithms to base two before the feature-specific scale. */
	private static final double INV_LOG_2=1/Math.log(2);
}
