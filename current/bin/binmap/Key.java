package bin.binmap;

import bin.Bin;
import shared.Tools;

/**
 * Quantization utility for binning genomic data by dimensions such as GC content,
 * HH (Homo/Het), CAGA, and Depth.
 *
 * @author Brian Bushnell
 * @date January 16, 2026
 */
public abstract class Key implements Cloneable {

	public static Key makeKey() {
		if(defaultType==GHDType) {
			return new KeyGHD();
		}else if(defaultType==GDDType) {
			return new KeyGDD();
		}else if(defaultType==GHCType) {
			return new KeyGHC();
		}else if(defaultType==GHCDType) {
			return new KeyGHCD();
		}else if(defaultType==GHDDType) {
			return new KeyGHDD();
		}else if(defaultType==GHCDDType) {
			return new KeyGHCDD();
		}else if(defaultType==GHDDDType) {
			return new KeyGHDDD();
		}else if(defaultType==ZEROType) {
			return new Key0();
		}else if(defaultType==GType) {
			return new KeyG();
		}else if(defaultType==GHType) {
			return new KeyGH();
		}else if(defaultType==GCType) {
			return new KeyGC();
		}else if(defaultType==GDType) {
			return new KeyGD();
		}else {

		}
		throw new RuntimeException("Unknown key type "+defaultType);
	}

	/**
	 * Sets this Key's values based on a Bin's characteristics.
	 * @param a The Bin to extract values from
	 * @return This Key instance for method chaining
	 */
	public abstract Key set(Bin a);

	/**
	 * Directly sets the quantized dimensions.
	 * @return This Key instance for method chaining
	 */
	public Key setLevel(int gc, int d2, int d3) {
		gcLevel=gc;
		dim2=d2;
		dim3=d3;
		return this;
	}

	/**
	 * Directly sets the quantized dimensions.
	 * @return This Key instance for method chaining
	 */
	public Key setLevel(int gc, int d2, int d3, int d4) {
		gcLevel=gc;
		dim2=d2;
		dim3=d3;
		dim4=d4;
		return this;
	}

	/**
	 * Directly sets the quantized dimensions.
	 * @return This Key instance for method chaining
	 */
	public Key setLevel(int gc, int d2, int d3, int d4, int d5) {
		gcLevel=gc;
		dim2=d2;
		dim3=d3;
		dim4=d4;
		dim5=d5;
		return this;
	}

//	/**
//	 * Sets Key values using continuous measurements.
//	 * Automatically quantizes the input values into discrete levels.
//	 * @return This Key instance for method chaining
//	 */
//	public abstract Key setValue(float gc, float d2, float d3);
//
//	/**
//	 * Sets Key values using continuous measurements.
//	 * Automatically quantizes the input values into discrete levels.
//	 * @return This Key instance for method chaining
//	 */
//	public abstract Key setValue(float gc, float d2, float d3, float d4);

	@Override
	public final boolean equals(Object other) {
		return equals((Key)other);
	}

	/**
	 * Compares this Key with another Key for equality.
	 * Keys are equal if all three quantization levels match.
	 * @param b The Key to compare with
	 * @return true if both Keys have identical quantization levels
	 */
	public final boolean equals(Key b) {
		return gcLevel==b.gcLevel && dim2==b.dim2 && dim3==b.dim3 && dim4==b.dim4 && dim5==b.dim5;
	}

	@Override
	public final int hashCode() {
		//Shifted to avoid collisions. 
		//Depth is usually 0-60 (6 bits). HH is 0-50 (6 bits). GC is 0-50 (6 bits).
		//Allocating 7 bits each is safe and fast.
		return gcLevel+(dim2<<7)+(dim3<<14)+(dim4<<21)+Integer.rotateLeft(dim5, 28);
	}

	@Override
	public final Key clone() {
		try {
			return (Key)(super.clone());
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	public int lowerBoundDim1(float gc, int range, int minGrid, float maxGCDif) {
		return Tools.max(minGrid, gcLevel-range, quantizeGC(gc-maxGCDif));
	}
	
	public int upperBoundDim1(float gc, int range, int maxGrid, float maxGCDif) {
		return Tools.min(maxGrid, gcLevel+range, quantizeGC(gc+maxGCDif));
	}

	public abstract int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio);
	public abstract int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio);

	public abstract int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio);
	public abstract int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio);

	public int lowerBoundDim4(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {return 0;}
	public int upperBoundDim4(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {return 0;}

	public int lowerBoundDim5(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {return 0;}
	public int upperBoundDim5(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {return 0;}

	public static int lowerBoundDepth(float val, int level, int range, int minGrid, float maxDepthRatio) {
		// Depth Physics: Use division/multiplication
		int quant=quantizeDepth(val/maxDepthRatio);
//		System.err.println("lower: v="+val+", q="+quant+", mng="+minGrid+
//			", l="+level+", r="+range+", l-r="+(level-range)+" -> max="+
//			Tools.max(minGrid, level-range, quant));
		// Combine: Grid Bound vs (Key Index - Range) vs (Physics Lower Bound)
		return Tools.max(minGrid, level-range, quant);
	}

	public static int upperBoundDepth(float val, int level, int range, int maxGrid, float maxDepthRatio) {
		int quant=quantizeDepth(val*maxDepthRatio);
//		System.err.println("upper: v="+val+", q="+quant+", mxg="+maxGrid+
//			", l="+level+", r="+range+", l+r="+(level+range)+" -> min="+
//			Tools.min(maxGrid, level+range, quant));
		return Tools.min(maxGrid, level+range, quant);
	}

	public static int lowerBoundHH(float val, int level, int range, int minGrid, float maxGCDif) {
		// HH Physics: Use subtraction (Linear)
		// Note: Reusing maxGCDif for HH tolerance for now?
		int quant=quantizeHH(val-maxGCDif*hhDifMult); 
		return Tools.max(minGrid, level-range, quant);
	}

	public static int upperBoundHH(float val, int level, int range, int maxGrid, float maxGCDif) {
		int quant=quantizeHH(val+maxGCDif*hhDifMult);
		return Tools.min(maxGrid, level+range, quant);
	}

	public static int lowerBoundCAGA(float val, int level, int range, int minGrid, float maxGCDif) {
		int quant=quantizeCAGA(val-maxGCDif*cagaDifMult); 
		return Tools.max(minGrid, level-range, quant);
	}

	public static int upperBoundCAGA(float val, int level, int range, int maxGrid, float maxGCDif) {
		int quant=quantizeCAGA(val+maxGCDif*cagaDifMult);
		return Tools.min(maxGrid, level+range, quant);
	}

	/**
	 * Parses command-line arguments for Key quantization parameters.
	 * Handles gcwidth, hhwidth, depthwidth, and their multipliers.
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
		}else if(a.equalsIgnoreCase("hhwidth")){
			float f=Float.parseFloat(b);
			setHHWidth(f);
		}else if(a.equalsIgnoreCase("hhmult")){
			float f=Float.parseFloat(b);
			setHHMult(f);
		}else if(a.equalsIgnoreCase("hhdifmult")){
			float f=Float.parseFloat(b);
			hhDifMult=f;
		}else if(a.equalsIgnoreCase("cagawidth")){
			float f=Float.parseFloat(b);
			setCagaWidth(f);
		}else if(a.equalsIgnoreCase("cagamult")){
			float f=Float.parseFloat(b);
			setCagaMult(f);
		}else if(a.equalsIgnoreCase("cagadifmult")){
			float f=Float.parseFloat(b);
			cagaDifMult=f;
		}else if(a.equalsIgnoreCase("depthwidth")){
			float f=Float.parseFloat(b);
			setDepthWidth(f);
		}else if(a.equalsIgnoreCase("depthmult")){
			float f=Float.parseFloat(b);
			setDepthMult(f);
		}else if(a.equalsIgnoreCase("dimensions")){
			int d=Integer.parseInt(b);
			assert(d>0 && d<=4);
			dimensions=d;
		}else if(a.equalsIgnoreCase("key") || a.equalsIgnoreCase("keytype")){
			if(b.equalsIgnoreCase("auto")){
				defaultType=-1;
			}else if(b.equalsIgnoreCase("gdd")){
				defaultType=GDDType;
			}else if(b.equalsIgnoreCase("ghd")){
				defaultType=GHDType;
			}else if(b.equalsIgnoreCase("ghc")){
				defaultType=GHCType;
			}else if(b.equalsIgnoreCase("ghcd")){
				defaultType=GHCDType;
			}else if(b.equalsIgnoreCase("ghdd")){
				defaultType=GHDDType;
			}else if(b.equalsIgnoreCase("ghcdd")){
				defaultType=GHCDDType;
			}else if(b.equalsIgnoreCase("ghddd")){
				defaultType=GHDDDType;
			}else if(b.equalsIgnoreCase("zero") || b.equalsIgnoreCase("0")){
				defaultType=ZEROType;
			}else if(b.equalsIgnoreCase("g")){
				defaultType=GType;
			}else if(b.equalsIgnoreCase("gh")){
				defaultType=GHType;
			}else if(b.equalsIgnoreCase("gc")){
				defaultType=GCType;
			}else if(b.equalsIgnoreCase("gd")){
				defaultType=GDType;
			}else{
				throw new RuntimeException("Unknown key type "+arg);
			}
			setType=(defaultType>=0);
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Quantizes a depth/coverage value into a discrete level.
	 * Uses logarithmic scaling.
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

	/**
	 * Quantizes an HH ratio fraction into a discrete level.
	 * Uses linear scaling across the 0.0 to 1.0 range.
	 * @param hh HH ratio as a fraction (0.0 to 1.0)
	 * @return Quantized HH level as an integer
	 */
	public static int quantizeHH(float hh) {
		return (int)(Tools.mid(0,hh,1)*hhLevelMult);
	}

	/**
	 * Quantizes an CAGA ratio fraction into a discrete level.
	 * Uses linear scaling across the 0.0 to 1.0 range.
	 * @param caga CAGA ratio as a fraction (0.0 to 1.0)
	 * @return Quantized CAGA level as an integer
	 */
	public static int quantizeCAGA(float caga) {
		return (int)(Tools.mid(0,caga,1)*cagaLevelMult);
	}

	/** Sets GC quantization multiplier (inverse of width). */
	static final void setGCMult(float f) {
		assert(f>=2);
		setGCWidth(1/f);
	}

	/** Sets the GC quantization width parameter. */
	static final void setGCWidth(float f) {
		assert(f>0 && f<=0.5f);
		gcLevelWidth=f;
		gcLevelMult=1f/gcLevelWidth;
	}

	/** Sets HH quantization multiplier (inverse of width). */
	static final void setHHMult(float f) {
		assert(f>=2);
		setHHWidth(1/f);
	}

	/** Sets the HH quantization width parameter. */
	static final void setHHWidth(float f) {
		assert(f>0 && f<=0.5f);
		hhLevelWidth=f;
		hhLevelMult=1f/hhLevelWidth;
	}

	/** Sets CAGA quantization multiplier (inverse of width). */
	static final void setCagaMult(float f) {
		assert(f>=2);
		setCagaWidth(1/f);
	}

	/** Sets the CAGA quantization width parameter. */
	static final void setCagaWidth(float f) {
		assert(f>0 && f<=0.5f);
		cagaLevelWidth=f;
		cagaLevelMult=1f/cagaLevelWidth;
	}

	/** Sets depth quantization width (inverse of multiplier). */
	static final void setDepthWidth(float f) {
		assert(f>0);
		setDepthMult(1/f);
	}

	/** Sets the depth quantization multiplier parameter. */
	static final void setDepthMult(float f) {
		assert(f>0);
		depthLevelMult=f;
		maxDepthLevel=quantizeDepth(maxDepth);
	}
	
	static int setType(int samples, long contigs) {
		if(setType) {return defaultType;}
		return defaultType=pickType(samples, contigs);
	}
	
	static int pickType(int samples, long contigs) {
		if(samples<1) {return GHCType;}
		else if(contigs<200000) {return GHDType;}
		
		if(samples<2) {
			return contigs<1000000 ? GHCDType : GHDDType;
		}else if(samples<3) {
			return contigs<1000000 ? GHDDType : GHCDDType;
		}else {
			return contigs<1000000 ? GHDDType : GHDDDType;
		}
	}

	int gcLevel=0;
	int dim2=0;
	int dim3=0;
	int dim4=0;
	int dim5=0;

	private static final int GDDType=0, GHDType=1, GHCType=2, GHCDType=3, GHDDType=4, GHDDDType=5, GHCDDType=6, ZEROType=7, GType=8, GHType=9, GCType=10, GDType=11;
	protected static int dimensions=3;
	protected static int defaultType=GHDDType;
	private static boolean setType=false;

	/** Maximum depth value for quantization calculations */
	protected static final float maxDepth=1000000;

	/** Multiplier for depth quantization resolution */
	protected static float depthLevelMult=2f;
	protected static int maxDepthLevel=quantizeDepth(maxDepth);

	/** Width parameter for GC content quantization */
	protected static float gcLevelWidth=0.02f;
	/** Multiplier for GC content quantization */
	protected static float gcLevelMult=1f/gcLevelWidth;

	/** Width parameter for HH content quantization */
	protected static float hhLevelWidth=0.02f;
	/** Multiplier for HH content quantization */
	protected static float hhLevelMult=1f/hhLevelWidth;
	private static float hhDifMult=1f;

	/** Width parameter for CAGA content quantization */
	protected static float cagaLevelWidth=0.02f;
	/** Multiplier for CAGA content quantization */
	protected static float cagaLevelMult=1f/cagaLevelWidth;
	private static float cagaDifMult=1f;
}
