package bin.binmap;

import bin.Bin;

/**
 * Quantization utility for binning genomic data by GC content and Depth.
 * Uses a 3D key (GC, Depth) to shatter large clusters in the initial binning phase.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 16, 2026
 */
class KeyGD extends Key {

	/** Creates an uninitialized Key with default quantization levels */
	public KeyGD() {}

	/**
	 * Sets this Key's values based on a Bin's characteristics.
	 * Uses GC and primary Depth.
	 * @param a The Bin to extract values from
	 * @return This Key instance for method chaining
	 */
	public KeyGD set(Bin a) {
		return setValue(a.gc(), a.depth(0));
	}

	/**
	 * Directly sets the quantized levels for GC and Coverage.
	 *
	 * @param dim1_ Quantized GC content level
	 * @param depthLevel_ Quantized coverage level
	 * @return This Key instance for method chaining
	 */
	public KeyGD setLevel(int gcLevel_, int depthLevel_) {
		gcLevel=gcLevel_;
		dim2=depthLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=maxDepthLevel) : dim2+", "+maxDepthLevel;
		return this;
	}

	/**
	 * Sets Key values using continuous measurements.
	 * Automatically quantizes the input values into discrete levels.
	 *
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @param depth Primary coverage/depth value
	 * @return This Key instance for method chaining
	 */
	public KeyGD setValue(float gc, float depth) {
		assert(gc>=0 && gc<=1) : gc;
		assert(depth>=0) : depth;
		return setLevel(quantizeGC(gc), quantizeDepth(depth));
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
		return 0;
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return 0;
	}

	/** Returns a string representation showing the three quantization levels.
	 * @return String in format "(dim1,hhLevel,depthLevel)" */
	public String toString() {
		return "("+gcLevel+","+dim2+")";
	}
}