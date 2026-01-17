package bin;

/**
 * Quantization utility for binning genomic data by GC content, HH (Homo/Het), and Depth.
 * Uses a 3D key (GC, HH, Depth) to shatter large clusters in the initial binning phase.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 16, 2026
 */
class KeyGHD extends Key {

	/** Creates an uninitialized Key with default quantization levels */
	public KeyGHD() {}

	/**
	 * Sets this Key's values based on a Bin's characteristics.
	 * Uses GC, HH, and primary Depth.
	 * @param a The Bin to extract values from
	 * @return This Key instance for method chaining
	 */
	public KeyGHD set(Bin a) {
		//Assuming Bin.hh is accessible directly or via accessor; using field access based on prompt instructions.
		return setValue(a.gc(), a.hh, a.depth(0));
	}

	/**
	 * Directly sets the quantized levels for GC, HH, and Coverage.
	 *
	 * @param dim1_ Quantized GC content level
	 * @param hhLevel_ Quantized HH (homo/het) level
	 * @param depthLevel_ Quantized coverage level
	 * @return This Key instance for method chaining
	 */
	public KeyGHD setLevel(int gcLevel_, int hhLevel_, int depthLevel_) {
		gcLevel=gcLevel_;
		dim2=hhLevel_;
		dim3=depthLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)hhLevelMult);
		assert(dim3>=0 && dim3<=maxDepthLevel) : dim3+", "+maxDepthLevel;
		return this;
	}

	/**
	 * Sets Key values using continuous measurements.
	 * Automatically quantizes the input values into discrete levels.
	 *
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @param hh HH ratio as a fraction (0.0 to 1.0)
	 * @param depth Primary coverage/depth value
	 * @return This Key instance for method chaining
	 */
	public KeyGHD setValue(float gc, float hh, float depth) {
		assert(gc>=0 && gc<=1) : gc;
		assert(hh>=0 && hh<=1) : hh;
		assert(depth>=0) : depth;
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeDepth(depth));
	}

	@Override
	public Key setValue(float gc, float hh, float depth, float d4){
		assert(gc>=0 && gc<=1) : gc;
		assert(hh>=0 && hh<=1) : hh;
		assert(depth>=0) : depth;
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeDepth(depth));
	}

	@Override
	public KeyGHD clone() {return (KeyGHD)(super.clone());}

	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundHH(a.hh, dim2, range, minGrid, maxGCDif);
	}

	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundHH(a.hh, dim2, range, maxGrid, maxGCDif);
	}

	@Override
	public int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim3, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim3, range, maxGrid, maxDepthRatio);
	}

	/** Returns a string representation showing the three quantization levels.
	 * @return String in format "(dim1,hhLevel,depthLevel)" */
	public String toString() {
		return "("+gcLevel+","+dim2+","+dim3+")";
	}

	/** Quantized GC content level */
	int gcLevel() {return gcLevel;}
	/** Quantized HH (Homo/Het) level */
	int hhLevel() {return dim2;}
	/** Quantized coverage/depth level */
	int depthLevel() {return dim3;}

}