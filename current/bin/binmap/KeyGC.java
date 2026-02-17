package bin.binmap;

import bin.Bin;

class KeyGC extends Key {
	public KeyGC() {}
	
	@Override
	public KeyGC set(Bin a) {
		return setValue(a.gc(), a.caga);
	}
	
	public KeyGC setLevel(int gcLevel_, int cagaLevel_) {
		gcLevel=gcLevel_;
		dim2=cagaLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)cagaLevelMult);
		return this;
	}
	
	public KeyGC setValue(float gc, float caga) {
		return setLevel(quantizeGC(gc), quantizeCAGA(caga));
	}

	// Dim 2 is CAGA here, not HH
	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundCAGA(a.caga, dim2, range, minGrid, maxGCDif);
	}
	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundCAGA(a.caga, dim2, range, maxGrid, maxGCDif);
	}
	
	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return 0;}

	@Override
	public String toString() {return "("+gcLevel+","+dim2+")";}
}