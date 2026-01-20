package bin;

import java.util.Arrays;

/**
 * Quantization utility for binning genomic data by GC content, HH (Homo/Het), and CAGA.
 * Uses a 3D key (GC, HH, CAGA) based purely on sequence composition, ignoring depth.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 16, 2026
 */
class KeyGHC extends Key {

	public KeyGHC() {}

	@Override
	public KeyGHC set(Bin a) {
		assert(Float.isFinite(a.gc())) : a;
		assert(Float.isFinite(a.hh)) : a;
		assert(Float.isFinite(a.caga)) : "\n"+a.caga+"\n"+Arrays.toString(a.dimers)+
			"\n"+a+"\n"+(a.toSeq());
		return setValue(a.gc(), a.hh, a.caga);
	}

	public KeyGHC setLevel(int gcLevel_, int hhLevel_, int cagaLevel_) {
		gcLevel=gcLevel_;
		dim2=hhLevel_;
		dim3=cagaLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)hhLevelMult);
		assert(dim3>=0 && dim3<=(int)cagaLevelMult);
		return this;
	}
	
	public KeyGHC setValue(float gc, float hh, float caga) {
		assert(gc>=0 && gc<=1) : gc;
		assert(hh>=0 && hh<=1) : hh;
		assert(caga>=0 && caga<=1) : caga;
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeCAGA(caga));
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

	// Dim 3 = CAGA
	@Override
	public int lowerBoundDim3(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundCAGA(a.caga, dim3, range, minGrid, maxGCDif);
	}

	@Override
	public int upperBoundDim3(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundCAGA(a.caga, dim3, range, maxGrid, maxGCDif);
	}

	public String toString() {
		return "("+gcLevel+","+dim2+","+dim3+")";
	}
}