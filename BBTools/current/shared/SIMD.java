package shared;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import ml.Cell;

/**
 * Holds SIMD methods.
 * @author Brian Bushnell
 * @date Sep 12, 2023?
 *
 */
public final class SIMD{

	// Example from https://medium.com/@Styp/java-18-vector-api-do-we-get-free-speed-up-c4510eda50d2
	@SuppressWarnings("restriction")
	private static final VectorSpecies<Float> FSPECIES=FloatVector.SPECIES_256;// FloatVector.SPECIES_PREFERRED; //This needs to be final or performance drops.
	private static final int FWIDTH=FSPECIES.length();
//	private static final int boundMask=~(FWIDTH-1);

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Byte> BSPECIES=ByteVector.SPECIES_256;
	private static final int BWIDTH=BSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Integer> ISPECIES=IntVector.SPECIES_256;
	private static final int IWIDTH=ISPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Short> SSPECIES=ShortVector.SPECIES_256;
	private static final int SWIDTH=SSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Double> DSPECIES=DoubleVector.SPECIES_256;
	private static final int DWIDTH=DSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Long> LSPECIES=LongVector.SPECIES_256;
	private static final int LWIDTH=LSPECIES.length();

	@SuppressWarnings("restriction")
	/**
	 * Vectorized version of "c+=a[i]*b[i]" where a and b are equal-length arrays.
	 * @param a A vector to multiply.
	 * @param b A vector to multiply.
	 * @return Sum of products of vector elements.
	 */
	static final float fma(final float[] a, final float[] b){
		assert (a.length==b.length);

		// Note: FSPECIES=FloatVector.SPECIES_256 and FWIDTH=8
		final int limit=FSPECIES.loopBound(a.length);

		FloatVector sum=FloatVector.zero(FSPECIES);
		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
			sum=va.fma(vb, sum);
		}
		float c=sum.reduceLanes(VectorOperators.ADD);
		for(; i<a.length; i++){// Residual scalar loop
			c+=a[i]*b[i];
		}
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Vectorized version of "c+=a[i]*b[bSet[i]]" where a and bSet are equal-length arrays, and bSet stores indices of b, in ascending contiguous blocks of 8.
	 * @param a A vector to multiply.
	 * @param b A vector to multiply.
	 * @return Sum of products of vector elements.
	 */
	static final float fmaSparse(final float[] a, final float[] b, int[] bSet){
		assert (a.length==bSet.length);
		assert (a.length<b.length);// Otherwise should do normal fma

		// Note: FSPECIES=FloatVector.SPECIES_256 and FWIDTH=8
		final int limit=FSPECIES.loopBound(bSet.length);
//		assert(FWIDTH==8);
//		int elements=0;

		FloatVector sum=FloatVector.zero(FSPECIES);
		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			int idx=bSet[i];
//			elements+=FWIDTH;
//			assert(idx%8==0) : idx+", "+i+", "+Arrays.toString(bSet);
//			assert(bSet[i+1]==idx+1);
//			assert(bSet[i+7]==idx+7);
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, idx);
			sum=va.fma(vb, sum);
		}
		float c=sum.reduceLanes(VectorOperators.ADD);
		for(; i<bSet.length; i++){// Residual scalar loop
//			elements++;
			c+=a[i]*b[bSet[i]];
		}

//		float c2=0;
//		for (int j=0; j<bSet.length; j++) {//Verification loop
//			c2+=a[j]*b[bSet[j]];
//		}
//		assert(Tools.absdif(c, c2)<0.0001f) : c+", "+c2;
//		assert(elements==bSet.length);

		return c;
	}

	// Isla
	public static final float absDif(long[] a, long[] b, float inva, float invb){
		assert (a.length==b.length);

		final int length=a.length;
		final int limit=LSPECIES.loopBound(length);

		FloatVector sumVec=FloatVector.zero(FSPECIES);
		int i=0;

		// SIMD loop for aligned portion
		for(; i<limit; i+=LWIDTH){
			LongVector va=LongVector.fromArray(LSPECIES, a, i);
			LongVector vb=LongVector.fromArray(LSPECIES, b, i);

			// Convert longs to floats
			FloatVector fa=(FloatVector)va.convertShape(VectorOperators.L2F, FSPECIES, 0);
			FloatVector fb=(FloatVector)vb.convertShape(VectorOperators.L2F, FSPECIES, 0);

			// Apply scaling factors
			fa=fa.mul(inva);
			fb=fb.mul(invb);

			// Calculate absolute difference and accumulate
			FloatVector diff=fa.sub(fb).abs();
			sumVec=sumVec.add(diff);
		}

		// For residual elements, just use scalar loop
		float sum=sumVec.reduceLanes(VectorOperators.ADD);
		for(; i<length; i++){
			float ai=a[i]*inva;
			float bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}

		return sum;
	}

	// Isla
	// Unfortunately, this dumps core, is very slow, and gives the wrong answer
	public static float absDifComp(long[] a, long[] b, int k, int[] gcmap){
		final int length=a.length;

		// Calculate GC bucket sums - this can't be easily vectorized
		final float[] aSums=new float[k+1];
		final float[] bSums=new float[k+1];

		for(int i=0; i<length; i++){
			int gc=gcmap[i];
			aSums[gc]+=a[i];
			bSums[gc]+=b[i];
		}

		final float inv=1f/(k+1);

		// Compute normalization factors
		for(int i=0; i<=k; i++){
			aSums[i]=inv/Math.max(aSums[i], 1);
			bSums[i]=inv/Math.max(bSums[i], 1);
		}

		// Process by GC content - one efficient way to vectorize with different scaling factors
		FloatVector sumVec=FloatVector.zero(FSPECIES);

		// Process elements grouped by their GC content
		for(int gc=0; gc<=k; gc++){
			float aFactor=aSums[gc];
			float bFactor=bSums[gc];

			// Find all indices with this GC content in chunks for efficient processing
			int currentIndex=0;
			while(currentIndex<length){
				// Find start of a chunk with this GC
				while(currentIndex<length && gcmap[currentIndex]!=gc){ currentIndex++; }

				// If we found a starting point
				if(currentIndex<length){
					int chunkStart=currentIndex;

					// Find end of the chunk
					while(currentIndex<length && gcmap[currentIndex]==gc){ currentIndex++; }

					int chunkEnd=currentIndex;
					int chunkSize=chunkEnd-chunkStart;

					// Process this chunk with SIMD
					if(chunkSize>=LWIDTH){
						int limit=chunkStart+(chunkSize/LWIDTH)*LWIDTH;

						for(int i=chunkStart; i<limit; i+=LWIDTH){
							LongVector va=LongVector.fromArray(LSPECIES, a, i);
							LongVector vb=LongVector.fromArray(LSPECIES, b, i);

							// Convert to float
							FloatVector fa=(FloatVector)va.convertShape(VectorOperators.L2F, FSPECIES, 0);
							FloatVector fb=(FloatVector)vb.convertShape(VectorOperators.L2F, FSPECIES, 0);

							// Apply GC-specific scaling
							fa=fa.mul(aFactor);
							fb=fb.mul(bFactor);

							// Calculate absolute difference
							sumVec=sumVec.add(fa.sub(fb).abs());
						}

						// Handle remainder
						for(int i=limit; i<chunkEnd; i++){
							float aComp=a[i]*aFactor;
							float bComp=b[i]*bFactor;
							sumVec=sumVec.add(Math.abs(aComp-bComp));
						}
					}else{
						// Small chunk - handle with scalar code
						for(int i=chunkStart; i<chunkEnd; i++){
							float aComp=a[i]*aFactor;
							float bComp=b[i]*bFactor;
							sumVec=sumVec.add(Math.abs(aComp-bComp));
						}
					}
				}
			}
		}

		float sum=sumVec.reduceLanes(VectorOperators.ADD);
		return Tools.mid(0, 1, (Float.isFinite(sum) && sum>0 ? sum : 0));
	}

	/**
	 * This is here to keep all the vector operations in a single loop, to prevent going in and out of SIMD mode too often... hopefully. ~20% measured speed increase compared to calling fma() for
	 * ScoreSequence.
	 */
	@SuppressWarnings("restriction")
	public static void feedForward(final Cell[] layer, final float[] b){
		assert (false)
			:"This was giving incorrect results for nets made made with simd=f and vice versa.  Needs validation.";
		final int limit=FSPECIES.loopBound(b.length);

		for(int cnum=0; cnum<layer.length; cnum++){
			final Cell cell=layer[cnum];
			final float[] a=cell.weights;
			FloatVector sum=FloatVector.zero(FSPECIES);
			for(int i=0; i<limit; i+=FWIDTH){// SIMD loop
				FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
				FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
				sum=va.fma(vb, sum);
			}
			cell.sum=sum.reduceLanes(VectorOperators.ADD);
		}

		for(int cnum=0; cnum<layer.length; cnum++){
			final Cell cell=layer[cnum];
			final float[] a=cell.weights;
			float residual=cell.bias;
			for(int i=limit+FWIDTH; i<a.length; i++){// Residual scalar loop
				residual+=a[i]*b[i];
			}
			cell.sum+=residual;
			final float v=(float)cell.activation(cell.sum);
			cell.setValue(v);
		}
	}

	/**
	 * This is here to keep all the vector operations in a single loop, to prevent going in and out of SIMD mode too often... hopefully. ~20% measured speed increase compared to calling fma() for Train.
	 */
	public static void backPropFma(Cell[] layer, float[] a, float[][] bb){
		final int limit=FSPECIES.loopBound(a.length);

		for(int cnum=0; cnum<layer.length; cnum++){
			Cell cell=layer[cnum];
			float[] b=bb[cnum];
			FloatVector sum=FloatVector.zero(FSPECIES);
			for(int i=0; i<limit; i+=FWIDTH){// SIMD loop
				FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
				FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
				sum=va.fma(vb, sum);
			}
			cell.eTotalOverOut=sum.reduceLanes(VectorOperators.ADD);
		}

		if(limit+FWIDTH>=a.length){ return; }// Shortcut when length is divisible by 8.

		for(int cnum=0; cnum<layer.length; cnum++){
			Cell cell=layer[cnum];
			float[] b=bb[cnum];
			float residual=0;
			for(int i=limit+FWIDTH; i<a.length; i++){// Residual scalar loop
				residual+=a[i]*b[i];
			}
			cell.eTotalOverOut+=residual;
		}
	}

	/**
	 * Performs "a+=incr" where a and incr are equal-length arrays.
	 * @param a A vector to increment.
	 * @param b Increment amount.
	 */
	@SuppressWarnings("restriction")
	static final void add(final float[] a, final float[] b){
		// final int width=SPECIES.length();
		final int limit=FSPECIES.loopBound(a.length);
//		final int limit=a.length&boundMask;

		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
			FloatVector sum=va.add(vb);
			sum.intoArray(a, i);
		}
		for(; i<a.length; i++){// Residual scalar loop; TODO: replace with vector mask
			a[i]+=b[i];
		}
	}

	/**
	 * Performs "a+=b*mult" where a and b are equal-length arrays.
	 * @param a A vector to increment.
	 * @param b Increment amount.
	 * @param mult Increment multiplier.
	 */
	@SuppressWarnings("restriction")
	static final void addProduct(final float[] a, final float[] b, final float mult){
		// final int width=SPECIES.length();
		final int limit=FSPECIES.loopBound(a.length);
//		final int limit=a.length&boundMask;

		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
			FloatVector sum=va.add(vb.mul(mult));
			sum.intoArray(a, i);
		}
		for(; i<a.length; i++){// Residual scalar loop; TODO: replace with vector mask
			a[i]+=b[i]*mult;
		}
	}

	@SuppressWarnings("restriction")
	static final void addProductSparse(final float[] a, final float[] b, final int[] bSet,
		final float mult){
		// final int width=SPECIES.length();
		final int limit=FSPECIES.loopBound(bSet.length);
//		final int limit=a.length&boundMask;

		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			int idx=bSet[i];
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, idx);
			FloatVector sum=va.add(vb.mul(mult));
			sum.intoArray(a, i);
		}
		for(; i<bSet.length; i++){// Residual scalar loop; TODO: replace with vector mask
			a[i]+=b[bSet[i]]*mult;
		}
	}

	// a is dest
	@SuppressWarnings("restriction")
	static final void copy(final float[] a, final float[] b){
		final int length=Tools.min(a.length, b.length);
		// final int width=SPECIES.length();
		final int limit=FSPECIES.loopBound(length);
//		final int limit=a.length&boundMask;

		int i=0;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
			vb.intoArray(a, i);
		}
		for(; i<length; i++){// Residual scalar loop; TODO: replace with vector mask
			a[i]=b[i];
		}
	}

	/** Returns number of matches */
	@SuppressWarnings("restriction")
	static final int countMatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		final int length=b2-a2+1;
		final int limit0=BSPECIES.loopBound(length);
		final int limit=a2+limit0;

		int i=a1, j=a2;
		int matches=0;
		for(; j<limit; i+=BWIDTH, j+=BWIDTH){// SIMD loop
			ByteVector v1=ByteVector.fromArray(BSPECIES, s1, i);
			ByteVector v2=ByteVector.fromArray(BSPECIES, s2, j);
			VectorMask<Byte> x=v1.eq(v2);
			matches+=x.trueCount();// This might be slow, or might not
		}
		for(; j<=b2; i++, j++){
			final byte x=s1[i], y=s2[j];
			final int m=((x==y) ? 1 : 0);
			matches+=m;
		}
		return matches;
	}

	/** Returns index of symbol */
	@SuppressWarnings("restriction")
	static final int find(final byte[] a, final byte symbol, final int from, final int to){// 15% Slower than scalar code, at least for ByteFile1
		final int length=to-from;// Intentionally exclusive
		final int limit0=BSPECIES.loopBound(length);
		final int limit=from+limit0;

		int pos=from;
		for(; pos<limit; pos+=BWIDTH){// SIMD loop
			ByteVector v=ByteVector.fromArray(BSPECIES, a, pos);
			VectorMask<Byte> x=v.eq(symbol);
			int t=x.firstTrue();
			if(t<BWIDTH){ return pos+t; }
//			if(x.anyTrue()) {break;}
		}
		while(pos<to && a[pos]!=symbol){ pos++; }
		return pos;
	}

	@SuppressWarnings("restriction")
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final float sum(final float[] a, final int from, final int to){
		final int length=to-from+1;// Intentionally inclusive
		final int limit0=FSPECIES.loopBound(length);
		final int limit=from+limit0;

		FloatVector sum=FloatVector.zero(FSPECIES);
		int i=from;
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			sum=sum.add(va);
		}
		float c=sum.reduceLanes(VectorOperators.ADD);
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final long sum(final long[] a, final int from, final int to){
		final int length=to-from+1;
		final int limit0=LSPECIES.loopBound(length);
		final int limit=from+limit0;

		LongVector sum=LongVector.zero(LSPECIES);
		int i=from;
		for(; i<limit; i+=LWIDTH){// SIMD loop
			LongVector va=LongVector.fromArray(LSPECIES, a, i);
			sum=sum.add(va);
		}
		long c=sum.reduceLanes(VectorOperators.ADD);
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final long sum(final int[] a, final int from, final int to){// Tested as 1.5x scalar speed
		final int length=to-from+1;
		final int limit0=ISPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		long c=0;
		for(; i<limit; i+=IWIDTH){// SIMD loop
			IntVector va=IntVector.fromArray(ISPECIES, a, i);
			c+=va.reduceLanesToLong(VectorOperators.ADD);// This is probably slow
		}
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final long sum(final byte[] a, final int from, final int to){// Tested as 4x scalar speed
		// TODO: Test speed.
		final int length=to-from+1;
		final int limit0=BSPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		long c=0;
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, a, i);
			c+=va.reduceLanesToLong(VectorOperators.ADD);
		}
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final double sum(final double[] a, final int from, final int to){
		final int length=to-from+1;
		final int limit0=DSPECIES.loopBound(length);
		final int limit=from+limit0;

		DoubleVector sum=DoubleVector.zero(DSPECIES);
		int i=from;
		for(; i<limit; i+=DWIDTH){// SIMD loop
			DoubleVector va=DoubleVector.fromArray(DSPECIES, a, i);
			sum=sum.add(va);
		}
		double c=sum.reduceLanes(VectorOperators.ADD);
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
		return c;
	}

//	public static float absDifFloat(float[] a, float[] b) {
//		assert(a.length==b.length);
//		float sum=0;
//		for(int i=0; i<a.length; i++){
//			sum+=Math.abs(a[i]-b[i]);
//		}
//		return (float)sum;
//	}

	/**
	 * Calculates the sum of the absolute differences between corresponding elements of two float arrays.
	 *
	 * @param a the first float array
	 * @param b the second float array
	 * @return the sum of the absolute differences between corresponding elements of the two arrays
	 * @throws IllegalArgumentException if the lengths of the arrays do not match
	 */
	public static float absDifFloat(float[] a, float[] b){
		if(a.length!=b.length){
			throw new IllegalArgumentException("Arrays must have the same length");
		}

		final int length=a.length;
		final int limit0=FSPECIES.loopBound(length);
		final int limit=limit0;

		FloatVector sumVec=FloatVector.zero(FSPECIES);
		int i=0;
		for(; i<limit; i+=FWIDTH){ // SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i);
			FloatVector diff=va.sub(vb).abs();
			sumVec=sumVec.add(diff);
		}

		// Scalar residual loop
//        float sum=sumVec.reduceLanes(VectorOperators.ADD);
//        for (; i<length; i++) { // Residual scalar loop
//            sum+=Math.abs(a[i]-b[i]);
//        }
//        return sum;

		// Handle the residual elements using lanewise masking
		// Lightly tested and seems to work
		if(i<length){
			VectorMask<Float> mask=FSPECIES.indexInRange(i, length);
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i, mask);
			FloatVector vb=FloatVector.fromArray(FSPECIES, b, i, mask);
			FloatVector diff=va.sub(vb).abs();
			sumVec=sumVec.add(diff, mask);
		}
		float sum=sumVec.reduceLanes(VectorOperators.ADD);
		return sum;
	}

	// Isla
	public static float cosineSimilarity(int[] a, int[] b, float inva, float invb){
		assert (a.length==b.length);

		int length=a.length;
		int upperBound=ISPECIES.loopBound(length);

		// Accumulation vectors
		FloatVector dotProductVec=FloatVector.zero(FSPECIES);
		FloatVector normVec1Vec=FloatVector.zero(FSPECIES);
		FloatVector normVec2Vec=FloatVector.zero(FSPECIES);

		int i=0;
		for(; i<upperBound; i+=IWIDTH){
			IntVector va=IntVector.fromArray(ISPECIES, a, i);
			IntVector vb=IntVector.fromArray(ISPECIES, b, i);

			FloatVector fa=(FloatVector)va.convertShape(VectorOperators.I2F, FSPECIES, 0);
			FloatVector fb=(FloatVector)vb.convertShape(VectorOperators.I2F, FSPECIES, 0);
			fa=fa.mul(inva);
			fb=fb.mul(invb);

			// Accumulate in vector space
			dotProductVec=dotProductVec.add(fa.mul(fb));
			normVec1Vec=normVec1Vec.add(fa.mul(fa));
			normVec2Vec=normVec2Vec.add(fb.mul(fb));
		}

		// Reduce once at the end
		float dotProduct=dotProductVec.reduceLanes(VectorOperators.ADD);
		float normVec1=normVec1Vec.reduceLanes(VectorOperators.ADD);
		float normVec2=normVec2Vec.reduceLanes(VectorOperators.ADD);

		// Handle remaining elements
		for(; i<length; i++){
			float ai=a[i]*inva;
			float bi=b[i]*invb;
			dotProduct+=ai*bi;
			normVec1+=ai*ai;
			normVec2+=bi*bi;
		}

		normVec1=Math.max(normVec1, 1e-15f);
		normVec2=Math.max(normVec2, 1e-15f);

		return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
	}

	@SuppressWarnings("restriction")
	/**
	 * Finds the maximum value.
	 * @param a A vector.
	 * @return The max.
	 */
	static final int max(final int[] a, final int from, final int to){// Tested as 5x scalar speed
		final int length=to-from+1;
		final int limit0=ISPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		IntVector max=IntVector.broadcast(ISPECIES, a[from]);
		for(; i<limit; i+=IWIDTH){// SIMD loop
			IntVector va=IntVector.fromArray(ISPECIES, a, i);
			max=max.max(va);
		}
		int c=max.reduceLanes(VectorOperators.MAX);
		for(; i<=to; i++){// Residual scalar loop
			final int x=a[i];
			c=(x>c ? x : c);
		}
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Finds the maximum value.
	 * @param a A vector.
	 * @return The max.
	 */
	static final long max(final long[] a, final int from, final int to){
		final int length=to-from+1;
		final int limit0=LSPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		LongVector max=LongVector.broadcast(LSPECIES, a[from]);
		for(; i<limit; i+=LWIDTH){// SIMD loop
			LongVector va=LongVector.fromArray(LSPECIES, a, i);
			max=max.max(va);
		}
		long c=max.reduceLanes(VectorOperators.MAX);
		for(; i<=to; i++){// Residual scalar loop
			final long x=a[i];
			c=(x>c ? x : c);
		}
		return c;
	}

	@SuppressWarnings("restriction")
	/**
	 * Finds the maximum value.
	 * @param a A vector.
	 * @return The max.
	 */
	static final float max(final float[] a, final int from, final int to){
		final int length=to-from+1;
		final int limit0=FSPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		FloatVector max=FloatVector.broadcast(FSPECIES, a[from]);
		for(; i<limit; i+=FWIDTH){// SIMD loop
			FloatVector va=FloatVector.fromArray(FSPECIES, a, i);
			max=max.max(va);
		}
		float c=max.reduceLanes(VectorOperators.MAX);
		for(; i<=to; i++){// Residual scalar loop
			final float x=a[i];
			c=(x>c ? x : c);
		}
		return c;
	}

}
