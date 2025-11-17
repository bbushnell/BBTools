package bin;

import shared.Tools;

/**
 * Quantization utility for binning genomic data by GC content and depth.
 * Converts continuous genomic measurements into discrete quantization levels
 * for efficient comparison and binning operations.
 *
 * @author Brian Bushnell
 * @date 2014
 */
class Key implements Cloneable {
	
	/**
	 * Parses command-line arguments for Key quantization parameters.
	 * Handles gcwidth, gcmult, depthwidth, and depthmult parameter settings.
	 *
	 * @param arg The full argument string (unused)
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	static boolean parse(String arg, String a, String b) {
		
		if(a.equalsIgnoreCase("gcwidth")){
			float f=Float.parseFloat(b);
			setGCWidth(f);
		}else if(a.equalsIgnoreCase("gcmult")){
			float f=Float.parseFloat(b);
			setGCMult(f);
		}else if(a.equalsIgnoreCase("depthwidth")){
			float f=Float.parseFloat(b);
			setDepthWidth(f);
		}else if(a.equalsIgnoreCase("depthmult")){
			float f=Float.parseFloat(b);
			setDepthMult(f);
		}else {
			return false;
		}
		return true;
	}
	
	/**
	 * Creates a Key with specified GC content and coverage values.
	 * Quantizes the continuous values into discrete levels.
	 *
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @param cov Primary coverage/depth value
	 * @param cov2 Secondary coverage/depth value
	 */
	public Key(float gc, float cov, float cov2) {
		setValue(gc, cov, cov2);
	}
	
	/** Creates an uninitialized Key with default quantization levels */
	public Key() {}

	/**
	 * Sets this Key's values based on a Bin's characteristics.
	 * @param a The Bin to extract GC content and depth values from
	 * @return This Key instance for method chaining
	 */
	public Key set(Bin a) {
		return setValue(a.gc(), a.depth(0), a.depth(1));
	}
	
	/**
	 * Directly sets the quantized levels for GC content and coverage.
	 *
	 * @param gcLevel_ Quantized GC content level
	 * @param covLevel_ Primary quantized coverage level
	 * @param covLevel2_ Secondary quantized coverage level
	 * @return This Key instance for method chaining
	 */
	public Key setLevel(int gcLevel_, int covLevel_, int covLevel2_) {
		gcLevel=gcLevel_;
		covLevel=covLevel_;
		covLevel2=covLevel2_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(covLevel>=0 && covLevel<=maxDepthLevel) : covLevel+", "+maxDepthLevel;
		assert(covLevel2>=0 && covLevel2<=maxDepthLevel) : covLevel2+", "+maxDepthLevel;
		return this;
	}
	
	/**
	 * Sets Key values using continuous measurements.
	 * Automatically quantizes the input values into discrete levels.
	 *
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @param cov Primary coverage/depth value
	 * @param cov2 Secondary coverage/depth value
	 * @return This Key instance for method chaining
	 */
	public Key setValue(float gc, float cov, float cov2) {
		assert(gc>=0 && gc<=1) : gc;
		assert(cov>=0) : cov;
		assert(cov2>=0) : cov;
		return setLevel(quantizeGC(gc), quantizeDepth(cov), quantizeDepth(cov2));
	}
	
	@Override
	public boolean equals(Object other) {
		return equals((Key)other);
	}
	
	/**
	 * Compares this Key with another Key for equality.
	 * Keys are equal if all three quantization levels match.
	 * @param b The Key to compare with
	 * @return true if both Keys have identical quantization levels
	 */
	public boolean equals(Key b) {
		return gcLevel==b.gcLevel && covLevel==b.covLevel && covLevel2==b.covLevel2;
	}
	
	@Override
	public int hashCode() {
		return covLevel+(covLevel2<<10)+(gcLevel<<20);
	}
	
	@Override
	public Key clone() {
		try {
			return (Key)(super.clone());
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	//This is probably faster but not as simple to explain or adjust
//	public static int quantizeDepth(float depth) {
//		float yf=depth*depth*16;
//		long y=(long)yf;
//		int zeros=Long.numberOfLeadingZeros(y);
//		int level=Tools.min(maxDepthLevel, 64-zeros);
//		return level;
//	}

	/**
	 * Quantizes a depth/coverage value into a discrete level.
	 * Uses logarithmic scaling with low-depth offset for better resolution
	 * at low coverage values.
	 *
	 * @param depth The depth/coverage value to quantize
	 * @return Quantized depth level as an integer
	 */
	public static int quantizeDepth(float depth) {
		depth=Tools.min(depth, maxDepth);
		float yf=((float)(Tools.log2(depth+0.0625f)+4));
		int level=(int)(yf*depthLevelMult);
		return level;
	}

	/**
	 * Quantizes a GC content fraction into a discrete level.
	 * Uses linear scaling across the 0.0 to 1.0 GC range.
	 * @param gc GC content as a fraction (0.0 to 1.0)
	 * @return Quantized GC level as an integer
	 */
	public static int quantizeGC(float gc) {
		return (int)(Tools.mid(0,gc,1)*gcLevelMult);
	}
	
	/** Returns a string representation showing the three quantization levels.
	 * @return String in format "(gcLevel,covLevel,covLevel2)" */
	public String toString() {
		return "("+gcLevel+","+covLevel+","+covLevel2+")";
	}
	
	/** Sets GC quantization multiplier (inverse of width).
	 * @param f Multiplier value (must be >= 2) */
	static final void setGCMult(float f) {
		assert(f>=2);
		setGCWidth(1/f);
	}
	
	/**
	 * Sets the GC quantization width parameter.
	 * Smaller values provide finer GC resolution.
	 * @param f Width value (must be > 0 and <= 0.5)
	 */
	static final void setGCWidth(float f) {
		assert(f>0 && f<=0.5f);
		gcLevelWidth=f;
		gcLevelMult=1f/gcLevelWidth;
	}
	
	/** Sets depth quantization width (inverse of multiplier).
	 * @param f Width value (must be > 0) */
	static final void setDepthWidth(float f) {
		assert(f>0);
		setDepthMult(1/f);
	}
	
	/**
	 * Sets the depth quantization multiplier parameter.
	 * Higher values provide finer depth resolution.
	 * @param f Multiplier value (must be > 0)
	 */
	static final void setDepthMult(float f) {
		assert(f>0);
		depthLevelMult=f;
		maxDepthLevel=quantizeDepth(maxDepth);
		assert(maxDepthLevel>0) : "maxDepthLevel="+maxDepthLevel+", depthLevelMult="+depthLevelMult+
			", maxDepth="+maxDepth+", yf="+(Tools.log2(maxDepth+0.0625f)+4);
	}
	
	/** Quantized GC content level */
	int gcLevel;
	/** Primary quantized coverage/depth level */
	int covLevel;
	/** Secondary quantized coverage/depth level */
	int covLevel2;
	
	//gcwidth=0.01, depthwidth=0.25 seems faster, more sensitive, and more specific (halving both of them).
	/** Maximum depth value for quantization calculations */
	private static final float maxDepth=1000000;
	/** Multiplier for depth quantization resolution */
	private static float depthLevelMult=2f;
	/** Width parameter for GC content quantization */
	private static float gcLevelWidth=0.02f;
	/** Multiplier for GC content quantization (inverse of width) */
	private static float gcLevelMult=1f/gcLevelWidth;
	/** Maximum quantized depth level corresponding to maxDepth */
	private static int maxDepthLevel=quantizeDepth(maxDepth);
	
}
