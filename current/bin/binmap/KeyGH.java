package bin.binmap;

import bin.Bin;

class KeyGH extends Key {
	public KeyGH() {}
	
	@Override
	public KeyGH set(Bin a) {
		return setValue(a.gc(), a.hh);
	}
	
	public KeyGH setLevel(int gcLevel_, int hhLevel_) {
		gcLevel=gcLevel_;
		dim2=hhLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(dim2>=0 && dim2<=(int)hhLevelMult);
		return this;
	}
	
	public KeyGH setValue(float gc, float hh) {
		return setLevel(quantizeGC(gc), quantizeHH(hh));
	}

	@Override
	public int lowerBoundDim2(Bin a, int range, int minGrid, float maxGCDif, float maxDepthRatio) {
		return lowerBoundHH(a.hh, dim2, range, minGrid, maxGCDif);
	}
	@Override
	public int upperBoundDim2(Bin a, int range, int maxGrid, float maxGCDif, float maxDepthRatio) {
		return upperBoundHH(a.hh, dim2, range, maxGrid, maxGCDif);
	}
	
	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return 0;}

	@Override
	public String toString() {return "("+gcLevel+","+dim2+")";}
}