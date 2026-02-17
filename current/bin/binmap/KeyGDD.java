package bin.binmap;

import bin.Bin;
import shared.Tools;

/**
 * Quantization utility for binning genomic data by GC content and depth.
 * Converts continuous genomic measurements into discrete quantization levels
 * for efficient comparison and binning operations.
 *
 * @author Brian Bushnell
 * @date February 4, 2025
 */
class KeyGDD extends Key {

	/** Creates an uninitialized KeyGDD with default quantization levels */
	public KeyGDD() {}

	/**
	 * Sets this KeyGDD's values based on a Bin's characteristics.
	 * @param a The Bin to extract GC content and depth values from
	 * @return This KeyGDD instance for method chaining
	 */
	public KeyGDD set(Bin a) {
		return setValue(a.gc(), a.depth(0), a.depth(1));
	}

	/**
	 * Directly sets the quantized levels for GC content and coverage.
	 *
	 * @param gcLevel_ Quantized GC content level
	 * @param covLevel_ Primary quantized coverage level
	 * @param covLevel2_ Secondary quantized coverage level
	 * @return This KeyGDD instance for method chaining
	 */
	public KeyGDD setLevel(int gcLevel_, int covLevel_, int covLevel2_) {
		gcLevel=gcLevel_;
		dim2=covLevel_;
		dim3=covLevel2_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=maxDepthLevel) : dim2+", "+maxDepthLevel;
		assert(dim3>=0 && dim3<=maxDepthLevel) : dim3+", "+maxDepthLevel;
		return this;
	}

	/**
	 * Sets KeyGDD values using continuous measurements.
	 * Automatically quantizes the input values into discrete levels.
	 *
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @param cov Primary coverage/depth value
	 * @param cov2 Secondary coverage/depth value
	 * @return This KeyGDD instance for method chaining
	 */
	public KeyGDD setValue(float gc, float cov, float cov2) {
		assert(gc>=0 && gc<=1) : gc;
		assert(cov>=0) : cov;
		assert(cov2>=0) : cov;
		return setLevel(quantizeGC(gc), quantizeDepth(cov), quantizeDepth(cov2));
	}

	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim2, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim2, range, maxGrid, maxDepthRatio);
	}

	@Override
	public int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<2 ? 0 : lowerBoundDepth(a.depth(1), dim3, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<2 ? 0 : upperBoundDepth(a.depth(1), dim3, range, maxGrid, maxDepthRatio);
	}

	/** Returns a string representation showing the three quantization levels.
	 * @return String in format "(gcLevel,covLevel,covLevel2)" */
	public String toString() {
		return "("+gcLevel+","+dim2+","+dim3+")";
	}

}
