package bin.binmap;

import bin.Bin;

/**
 * Quantization utility for binning genomic data by GC, HH, and two Depth dimensions.
 * A 4D key useful for multi-sample binning where differential coverage is key.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 16, 2026
 */
class KeyGHDD extends Key {

	public KeyGHDD() {}

	@Override
	public KeyGHDD set(Bin a) {
		return setValue(a.gc(), a.hh, a.depth(0), a.depth(1));
	}

	public KeyGHDD setLevel(int gcLevel_, int hhLevel_, int depthLevel_, int depthLevel2_) {
		gcLevel=gcLevel_;
		dim2=hhLevel_;
		dim3=depthLevel_;
		dim4=depthLevel2_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)hhLevelMult);
		assert(dim3>=0 && dim3<=maxDepthLevel);
		assert(dim4>=0 && dim4<=maxDepthLevel);
		return this;
	}
	
	public KeyGHDD setValue(float gc, float hh, float depth, float depth2) {
		assert(gc>=0 && gc<=1) : gc;
		assert(hh>=0 && hh<=1) : hh;
		assert(depth>=0) : depth;
		assert(depth2>=0) : depth2;
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeDepth(depth), quantizeDepth(depth2));
	}

	// Dim 2 = HH
	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundHH(a.hh, dim2, range, minGrid, maxGCDif);
	}

	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundHH(a.hh, dim2, range, maxGrid, maxGCDif);
	}

	// Dim 3 = Depth 1
	@Override
	public int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim3, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim3, range, maxGrid, maxDepthRatio);
	}

	// Dim 4 = Depth 2
	@Override
	public int lowerBoundDim4(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<2 ? 0 : lowerBoundDepth(a.depth(1), dim4, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim4(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<2 ? 0 : upperBoundDepth(a.depth(1), dim4, range, maxGrid, maxDepthRatio);
	}

	public String toString() {
		return "("+gcLevel+","+dim2+","+dim3+","+dim4+")";
	}
}