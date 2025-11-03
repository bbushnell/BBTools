package shared;

import dna.AminoAcid;
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
import structures.IntList;

/**
 * Holds SIMD methods.
 * @author Brian Bushnell
 * @date Sep 12, 2023?
 *
 */
final class SIMD{
	
	private static int maxVectorLength=ByteVector.SPECIES_PREFERRED.vectorBitSize();
	public static final int maxVectorLength() {return maxVectorLength;}

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
	static final float absDif(long[] a, long[] b, float inva, float invb){
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
	static float absDifComp(long[] a, long[] b, int k, int[] gcmap){
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
	static void feedForward(final Cell[] layer, final float[] b){
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
	static void backPropFma(Cell[] layer, float[] a, float[][] bb){
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

	//	static float absDifFloat(float[] a, float[] b) {
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
	static float absDifFloat(float[] a, float[] b){
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
	static float cosineSimilarity(int[] a, int[] b, float inva, float invb){
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

	/**
	 * Find positions of the given symbol in a byte array using SIMD.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param symbol Character to find
	 * @param positions IntList to store newline positions
	 * @return Number of symbols found, including pre-existing ones
	 */
	@SuppressWarnings("restriction")
	static final int findSymbols(final byte[] buffer, final int from, 
		final int to, final byte symbol, final IntList positions){
		final int limit=BSPECIES.loopBound(to-from);
		final ByteVector newlineVec=ByteVector.broadcast(BSPECIES, symbol);

		int i=from;

		// SIMD loop
		for(; i<from+limit; i+=BWIDTH){
			ByteVector vec=ByteVector.fromArray(BSPECIES, buffer, i);
			VectorMask<Byte> mask=vec.eq(newlineVec);
			//			if(!mask.anyTrue()) {continue;}//Hopefully common case - Not faster, maybe 1% slower
			// Convert mask to long bitmask
			long bits=mask.toLong();

			// Extract set bit positions using bit manipulation
			//Brian version
			//			for(int lane=0; bits!=0; lane++){//Looks strange but lane needs to be incremented
			//				int zeros=Long.numberOfTrailingZeros(bits);
			//				lane+=zeros;
			//				positions.add(lane+i);
			//				bits>>>=(zeros+1);
			//			}
			while(bits!=0){//Isla version - 5% faster
				int lane=Long.numberOfTrailingZeros(bits);
				positions.add(i+lane);
				bits&=(bits-1); // Clear lowest set bit
			}
		}

		// Residual scalar loop
		for(; i<to; i++){
			if(buffer[i]==symbol){positions.add(i);}
		}

		return positions.size();
	}
	
	/**
	 * Find the last symbol in buffer by scanning backwards using SIMD.
	 * @param buffer Buffer to scan
	 * @param limit Scan backwards from this position (exclusive)
	 * @return Position of last newline, or -1 if none found
	 */
	@SuppressWarnings("restriction")
	static int findLastSymbol(byte[] buffer, int limit, final byte symbol){
		final ByteVector newlineVec=ByteVector.broadcast(BSPECIES, symbol);
		
		// Start from the last aligned chunk and work backwards
		int i=((limit-1)/BWIDTH)*BWIDTH; // Round down to vector boundary
		
		// SIMD loop - work backwards
		for(; i>=0; i-=BWIDTH){
			ByteVector vec=ByteVector.fromArray(BSPECIES, buffer, i);
			VectorMask<Byte> mask=vec.eq(newlineVec);
			long bits=mask.toLong();
			
			if(bits!=0){
				// Found newline(s) in this chunk - find the highest bit
				int lane=63-Long.numberOfLeadingZeros(bits); // Highest set bit
				return i+lane;
			}
		}
		
		// Residual scalar loop for beginning of buffer
		for(i+=BWIDTH-1; i>=0; i--){
			if(buffer[i]==symbol){return i;}
		}
		
		return -1;
	}

	@SuppressWarnings("restriction")
	/**
	 * Adds a constant delta to all bytes in an array using SIMD.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 */
	static final void add(final byte[] array, final byte delta){
		if(array==null){return;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final ByteVector vdelta=ByteVector.broadcast(BSPECIES, delta);
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, array, i);
			ByteVector vresult=va.add(vdelta);
			vresult.intoArray(array, i);
		}
		for(; i<length; i++){// Residual scalar loop
			array[i]+=delta;
		}
	}

	@SuppressWarnings("restriction")
	/**
	 * Adds delta to all bytes, caps at minimum value, and returns the minimum encountered.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 * @param cap The minimum allowed value after addition.
	 * @return The minimum value encountered after addition (before capping).
	 */
	static final byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES, delta);
		ByteVector vcap=ByteVector.broadcast(BSPECIES, (byte)cap);
		ByteVector vmin=ByteVector.broadcast(BSPECIES, (byte)127);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, array, i);
			ByteVector vresult=va.add(vdelta);
			vmin=vmin.min(vresult);
			ByteVector vcapped=vresult.max(vcap);
			vcapped.intoArray(array, i);
		}

		int min=vmin.reduceLanes(VectorOperators.MIN);

		for(; i<length; i++){// Residual scalar loop
			int b=array[i]+delta;
			min=Math.min(min, b);
			array[i]=(byte)Math.max(cap, b);
		}

		return (byte)min;
	}

	@SuppressWarnings("restriction")
	/**
	 * Applies quality offset delta, zeros quality for N bases, caps others at 2.
	 * @param quals Quality array to modify in-place.
	 * @param bases Sequence bases array.
	 * @param delta The offset to add to quality scores.
	 */
	static final void applyQualOffset(final byte[] quals, final byte[] bases, final int delta){
		if(quals==null){return;}

		final int length=quals.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES, (byte)delta);
		ByteVector vn=ByteVector.broadcast(BSPECIES, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES, (byte)2);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vquals=ByteVector.fromArray(BSPECIES, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES, bases, i);

			// Add delta
			ByteVector vresult=vquals.add(vdelta);

			// Cap at 2 for non-N bases
			vresult=vresult.max(vcap2);

			// Create mask: where bases == 'N'
			VectorMask<Byte> maskN=vbases.eq(vn);

			// Blend: if N then 0, else capped result
			vresult=vresult.blend(vzero, maskN);

			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){// Residual scalar loop
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	@SuppressWarnings("restriction")
	/**
	 * Converts U to T and u to t in a byte array using SIMD.
	 * @param bases The base array to modify in-place.
	 */
	static final void uToT(final byte[] bases){
		if(bases==null){return;}

		final int length=bases.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vU=ByteVector.broadcast(BSPECIES, (byte)'U');
		ByteVector vu=ByteVector.broadcast(BSPECIES, (byte)'u');
		ByteVector vT=ByteVector.broadcast(BSPECIES, (byte)'T');
		ByteVector vt=ByteVector.broadcast(BSPECIES, (byte)'t');

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vbases=ByteVector.fromArray(BSPECIES, bases, i);

			// Create masks for U and u
			VectorMask<Byte> maskU=vbases.eq(vU);
			VectorMask<Byte> masku=vbases.eq(vu);

			// Replace U with T, u with t
			vbases=vbases.blend(vT, maskU);
			vbases=vbases.blend(vt, masku);

			vbases.intoArray(bases, i);
		}

		for(; i<length; i++){// Residual scalar loop
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}

	@SuppressWarnings("restriction")
	/**
	 * Converts lowercase letters to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean lowerCaseToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final byte a='a', N='N';

		ByteVector va=ByteVector.broadcast(BSPECIES, a);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);

			// Create mask: where b > a (lowercase)
			VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va);

			// Replace with N where lowercase
			ByteVector vresult=vb.blend(vN, maskLower);

			vresult.intoArray(array, i);
		}

		for(; i<length; i++){// Residual scalar loop
//			byte b=array[i];
//			b=(byte)(b<a ? b : N);
//			array[i]=b;
			array[i]=AminoAcid.lowerCaseToNocall[array[i]];
		}

		return true;
	}

	@SuppressWarnings("restriction")
	/**
	 * Converts dot, dash, and X to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean dotDashXToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final byte dot='.', dash='-', X='X', N='N';

		ByteVector vdot=ByteVector.broadcast(BSPECIES, dot);
		ByteVector vdash=ByteVector.broadcast(BSPECIES, dash);
		ByteVector vX=ByteVector.broadcast(BSPECIES, X);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);

			// Create masks for each character
			VectorMask<Byte> maskDot=vb.eq(vdot);
			VectorMask<Byte> maskDash=vb.eq(vdash);
			VectorMask<Byte> maskX=vb.eq(vX);

			// Combine masks: any of the three
			VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX);

			// Replace with N where mask is true
			ByteVector vresult=vb.blend(vN, maskAny);

			vresult.intoArray(array, i);
		}

		for(; i<length; i++){// Residual scalar loop
//			byte b=array[i];
//			b=(b==dot || b==dash || b==X) ? N : b;
//			array[i]=b;
			array[i]=AminoAcid.dotDashXToNocall[array[i]];
		}

		return true;
	}
	
	@SuppressWarnings("restriction")
	/**
	 * Checks if array contains common amino acids (E or L).
	 * @param array The byte array to check.
	 * @return true if E or L found (likely protein).
	 */
	static final boolean isProtein(final byte[] array){
		if(array==null){return false;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte E='E', L='L';
		
		ByteVector vE=ByteVector.broadcast(BSPECIES, E);
		ByteVector vL=ByteVector.broadcast(BSPECIES, L);
		
		VectorMask<Byte> foundProtein=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);
			
			VectorMask<Byte> isE=vb.eq(vE);
			VectorMask<Byte> isL=vb.eq(vL);
			
			foundProtein=isE.or(isL).or(foundProtein);
		}
		
		boolean protein=foundProtein.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
			byte b=array[i];
			boolean nuc=AminoAcid.baseToNumberExtended[b]>=0;
			boolean amino=AminoAcid.acidToNumberExtended[b]>=0;
//			protein|=(b==E || b==L);
			protein|=(amino && !nuc);
		}
		
		return protein;
	}
	
	/* */
	
	@SuppressWarnings("restriction")
	/**
	 * Converts to uppercase using bitmask. Returns false if non-letter found.
	 * @param array The byte array to modify in-place.
	 * @return false if any byte is outside A-Z range after conversion.
	 */
	static final boolean toUpperCase(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			vb.intoArray(array, i);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> valid=validLow.and(validHigh);
			
			invalid=valid.not().or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
//			byte b=(byte)(array[i] & mask);
//			array[i]=b;
//			success&=(b>=A && b<=Z);
			array[i]=AminoAcid.toUpperCase[array[i]];
		}
		
		return success;
	}

	@SuppressWarnings("restriction")
	/**
	 * Checks if all bytes are letters (case-insensitive check via mask).
	 * @param array The byte array to check.
	 * @return false if any non-letter found.
	 */
	static final boolean allLetters(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> valid=validLow.and(validHigh);
			
			invalid=valid.not().or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// TODO: Change to lookup table
			int b=(array[i] & mask);
			success&=(b>=A && b<=Z);
		}
		
		return success;
	}

	@SuppressWarnings("restriction")
	/**
	 * Converts IUPAC ambiguity codes to N, preserves A/C/G/T/U (case-insensitive).
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean iupacToN(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', C='C', G='G', T='T', U='U', N='N';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vC=ByteVector.broadcast(BSPECIES, C);
		ByteVector vG=ByteVector.broadcast(BSPECIES, G);
		ByteVector vT=ByteVector.broadcast(BSPECIES, T);
		ByteVector vU=ByteVector.broadcast(BSPECIES, U);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask for comparison
			ByteVector vb=vb0.and(vmask);
			
			// Check if A, C, G, T, or U
			VectorMask<Byte> isA=vb.eq(vA);
			VectorMask<Byte> isC=vb.eq(vC);
			VectorMask<Byte> isG=vb.eq(vG);
			VectorMask<Byte> isT=vb.eq(vT);
			VectorMask<Byte> isU=vb.eq(vU);
			VectorMask<Byte> isValid=isA.or(isC).or(isG).or(isT).or(isU);
			
			// Keep original if valid, replace with N if not
			ByteVector vresult=vb0.blend(vN, isValid.not());
			
			vresult.intoArray(array, i);
		}
		
		for(; i<length; i++){// Residual scalar loop
//			final int b0=array[i];
//			int b=(b0 & mask);
//			b=(b==A || b==C || b==G || b==T || b==U) ? b0 : N;
//			array[i]=(byte)b;
			array[i]=AminoAcid.baseToACGTN[array[i]];
		}
		
		return true;
	}

	@SuppressWarnings("restriction")
	/**
	 * Checks if all letters and no E or L (nucleotide validation).
	 * @param array The byte array to check.
	 * @return true if valid nucleotide sequence.
	 */
	static final boolean isNucleotide(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte E='E', L='L';
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		ByteVector vE=ByteVector.broadcast(BSPECIES, E);
		ByteVector vL=ByteVector.broadcast(BSPECIES, L);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> isLetter=validLow.and(validHigh);
			
			// Check if E or L
			VectorMask<Byte> isE=vb.eq(vE);
			VectorMask<Byte> isL=vb.eq(vL);
			VectorMask<Byte> isProtein=isE.or(isL);
			
			// Invalid if not letter OR is protein marker
			invalid=isLetter.not().or(isProtein).or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
//			int b=(array[i] & mask);
//			success&=((b>=A && b<=Z) && (b!=E && b!=L));
			success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);
		}
		
		return success;
	}

	/**
	 * Benchmark uToT vs uToT2
	 * Usage: java SIMD <iterations>
	 */
	static void main(String[] args){
		int iterations=args.length>0 ? Integer.parseInt(args[0]) : 1000;
		int numArrays=100;
		int arrayLength=1000;
		
		// Generate base arrays with random characters
		byte[][] baseArrays=new byte[numArrays][arrayLength];
		java.util.Random rand=new java.util.Random(42);
		for(int i=0; i<numArrays; i++){
			for(int j=0; j<arrayLength; j++){
				// Mix of various characters including U/u
				int r=rand.nextInt(10);
				if(r<2){baseArrays[i][j]='U';}
				else if(r<4){baseArrays[i][j]='u';}
				else if(r<5){baseArrays[i][j]='A';}
				else if(r<6){baseArrays[i][j]='C';}
				else if(r<7){baseArrays[i][j]='G';}
				else if(r<8){baseArrays[i][j]='T';}
				else if(r<9){baseArrays[i][j]='N';}
				else{baseArrays[i][j]=(byte)('a'+rand.nextInt(26));}
			}
		}
		long start1=0, start2=0, time1=0, time2=0;

		{
			// Working copies
			byte[][] workArrays1=new byte[numArrays][arrayLength];

			// Warmup
			for(int w=0; w<100; w++){
				for(int i=0; i<numArrays; i++){
					System.arraycopy(baseArrays[i], 0, workArrays1[i], 0, arrayLength);
				}
				for(int i=0; i<numArrays; i++){
					uToT(workArrays1[i]);
				}
			}

			// Benchmark uToT (original)
			start1=System.nanoTime();
			for(int iter=0; iter<iterations; iter++){
				for(int i=0; i<numArrays; i++){
					System.arraycopy(baseArrays[i], 0, workArrays1[i], 0, arrayLength);
				}
				for(int i=0; i<numArrays; i++){
					uToT(workArrays1[i]);
				}
			}
			time1=System.nanoTime()-start1;
		}
//		{
//
//			byte[][] workArrays2=new byte[numArrays][arrayLength];
//			// Warmup uToT2
//			for(int w=0; w<100; w++){
//				for(int i=0; i<numArrays; i++){
//					System.arraycopy(baseArrays[i], 0, workArrays2[i], 0, arrayLength);
//				}
//				for(int i=0; i<numArrays; i++){
//					uToT2(workArrays2[i]);
//				}
//			}
//
//			// Benchmark uToT2 (new)
//			start2=System.nanoTime();
//			for(int iter=0; iter<iterations; iter++){
//				for(int i=0; i<numArrays; i++){
//					System.arraycopy(baseArrays[i], 0, workArrays2[i], 0, arrayLength);
//				}
//				for(int i=0; i<numArrays; i++){
//					uToT2(workArrays2[i]);
//				}
//			}
//			time2=System.nanoTime()-start2;
//		}
//		// Benchmark uToT2 (new)
//		long start3=System.nanoTime();
//		for(int iter=0; iter<iterations; iter++){
//			for(int i=0; i<numArrays; i++){
//				System.arraycopy(baseArrays[i], 0, workArrays2[i], 0, arrayLength);
//			}
//			for(int i=0; i<numArrays; i++){
//				uToT3(workArrays2[i]);
//			}
//		}
//		long time3=System.nanoTime()-start3;
//		
//		// Benchmark uToT2 (new)
//		long start4=System.nanoTime();
//		for(int iter=0; iter<iterations; iter++){
//			for(int i=0; i<numArrays; i++){
//				System.arraycopy(baseArrays[i], 0, workArrays2[i], 0, arrayLength);
//			}
//			for(int i=0; i<numArrays; i++){
//				uToT4(workArrays2[i]);
//			}
//		}
//		long time4=System.nanoTime()-start4;
//		
//		// Benchmark uToT2 (new)
//		long start5=System.nanoTime();
//		for(int iter=0; iter<iterations; iter++){
//			for(int i=0; i<numArrays; i++){
//				System.arraycopy(baseArrays[i], 0, workArrays2[i], 0, arrayLength);
//			}
//			for(int i=0; i<numArrays; i++){
//				uToT5(workArrays2[i]);
//			}
//		}
//		long time5=System.nanoTime()-start5;
		
		// Verify results match
//		boolean match=true;
//		for(int i=0; i<numArrays; i++){
//			for(int j=0; j<arrayLength; j++){
//				if(workArrays1[i][j]!=workArrays2[i][j]){
//					match=false;
//					break;
//				}
//			}
//		}
		
		System.out.println("Arrays processed: "+numArrays);
		System.out.println("Iterations: "+iterations);
//		System.out.println("Results match: "+match);
		System.out.println();
		System.out.println("uToT  time: "+(time1/1_000_000)+" ms");
		System.out.println("uToT2 time: "+(time2/1_000_000)+" ms");
//		System.out.println("uToT3 time: "+(time3/1_000_000)+" ms");
//		System.out.println("uToT4 time: "+(time4/1_000_000)+" ms");
//		System.out.println("uToT5 time: "+(time5/1_000_000)+" ms");
		System.out.println("Speedup: "+String.format("%.2f", (double)time1/time2)+"x");
	}
	
//	static void uToT3(byte[] bases){//Fast
//		for(int i=0; i<bases.length; i++){
//			bases[i]=AminoAcid.uToT[bases[i]];
//		}
//	}
//	static void uToT4(byte[] bases){//Slow
//		final byte u='u', T='T';
//		for(int i=0; i<bases.length; i++){
//			byte b=bases[i];
//			bases[i]=((b|32)!=u) ? b : (byte)(T|(b&32));
//		}
//	}
//	static void uToT5(byte[] bases){//Slow
//		for(int i=0; i<bases.length; i++){// Residual scalar loop
//			if(bases[i]=='U'){bases[i]='T';}
//			else if(bases[i]=='u'){bases[i]='t';}
//		}
//	}

}
