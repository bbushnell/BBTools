package bin;

class KeyGHDDD extends Key {
	public KeyGHDDD() {}
	
	@Override
	public KeyGHDDD set(Bin a) {
		return setValue(a.gc(), a.hh, a.depth(0), a.depth(1), a.depth(2));
	}
	
	@Override
	public KeyGHDDD setLevel(int gc, int hh, int dp1, int dp2) {
		throw new RuntimeException("All 5 dims must be set.");
	}
	
	public KeyGHDDD setLevel(int gc, int hh, int dp1, int dp2, int dp3) {
		gcLevel=gc; dim2=hh; dim3=dp1; dim4=dp2; dim5=dp3;
		return this;
	}
	
	public KeyGHDDD setValue(float gc, float hh, float d1, float d2, float d3) {
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeDepth(d1), quantizeDepth(d2), quantizeDepth(d3));
	}
    
    // Dim 2 = HH
	@Override public int lowerBoundDim2(Bin a, int r, int min, float mg, float mdr) {return lowerBoundHH(a.hh, dim2, r, min, mg);}
	@Override public int upperBoundDim2(Bin a, int r, int max, float mg, float mdr) {return upperBoundHH(a.hh, dim2, r, max, mg);}

    // Dim 3 = Depth 0
	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim3, r, min, mdr);}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim3, r, max, mdr);}

    // Dim 4 = Depth 1
	@Override public int lowerBoundDim4(Bin a, int r, int min, float mg, float mdr) {return a.numDepths()<2 ? 0 : lowerBoundDepth(a.depth(1), dim4, r, min, mdr);}
	@Override public int upperBoundDim4(Bin a, int r, int max, float mg, float mdr) {return a.numDepths()<2 ? 0 : upperBoundDepth(a.depth(1), dim4, r, max, mdr);}

    // Dim 5 = Depth 2
	@Override public int lowerBoundDim5(Bin a, int r, int min, float mg, float mdr) {return a.numDepths()<3 ? 0 : lowerBoundDepth(a.depth(2), dim5, r, min, mdr);}
	@Override public int upperBoundDim5(Bin a, int r, int max, float mg, float mdr) {return a.numDepths()<3 ? 0 : upperBoundDepth(a.depth(2), dim5, r, max, mdr);}

	@Override public String toString() {return "("+gcLevel+","+dim2+","+dim3+","+dim4+","+dim5+")";}
}