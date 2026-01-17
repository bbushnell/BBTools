package bin;

/**
 * Quantization utility for binning by GC, HH, CAGA, and Depth.
 * A 4D key that maximizes sequence composition resolution plus abundance.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 16, 2026
 */
class KeyGHCD extends Key {

	public KeyGHCD() {}

	@Override
	public KeyGHCD set(Bin a) {
		return setValue(a.gc(), a.hh, a.caga, a.depth(0));
	}

	public KeyGHCD setLevel(int gcLevel_, int hhLevel_, int cagaLevel_, int depthLevel_) {
		gcLevel=gcLevel_;
		dim2=hhLevel_;
		dim3=cagaLevel_;
		dim4=depthLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)hhLevelMult);
		assert(dim3>=0 && dim3<=(int)cagaLevelMult);
		assert(dim4>=0 && dim4<=maxDepthLevel);
		return this;
	}

	@Override
	public KeyGHCD setValue(float gc, float hh, float caga, float depth) {
		assert(gc>=0 && gc<=1) : gc;
		assert(hh>=0 && hh<=1) : hh;
		assert(caga>=0 && caga<=1) : caga;
		assert(depth>=0) : depth;
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeCAGA(caga), quantizeDepth(depth));
	}

	@Override
	public Key setValue(float gc, float hh, float caga) {
		return setValue(gc, hh, caga, 0);
	}

	@Override
	public KeyGHCD clone() {return (KeyGHCD)(super.clone());}

	// Dim 2 = HH
	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundHH(a.hh, dim2, range, minGrid, maxGCDif);
	}

	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundHH(a.hh, dim2, range, maxGrid, maxGCDif);
	}

	// Dim 3 = CAGA
	@Override
	public int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundCAGA(a.caga, dim3, range, minGrid, maxGCDif);
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundCAGA(a.caga, dim3, range, maxGrid, maxGCDif);
	}

	// Dim 4 = Depth
	@Override
	public int lowerBoundDim4(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim4, range, minGrid, maxDepthRatio);
	}

	@Override
	public int upperBoundDim4(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim4, range, maxGrid, maxDepthRatio);
	}

	public String toString() {
		return "("+gcLevel+","+dim2+","+dim3+","+dim4+")";
	}
}