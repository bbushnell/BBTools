package bin;

class KeyGHCDD extends Key {
	public KeyGHCDD() {}
	
	@Override
	public KeyGHCDD set(Bin a) {
		return setValue(a.gc(), a.hh, a.caga, a.depth(0), a.depth(1));
	}
	
	@Override
	public KeyGHDDD setLevel(int gc, int hh, int caga, int dp1) {
		throw new RuntimeException("All 5 dims must be set.");
	}
	
	public KeyGHCDD setLevel(int gc, int hh, int caga, int dp1, int dp2) {
		gcLevel=gc; dim2=hh; dim3=caga; dim4=dp1; dim5=dp2;
		return this;
	}
	
	public KeyGHCDD setValue(float gc, float hh, float caga, float d1, float d2) {
		return setLevel(quantizeGC(gc), quantizeHH(hh), quantizeCAGA(caga), quantizeDepth(d1), quantizeDepth(d2));
	}
    
    // Dim 2 = HH
	@Override public int lowerBoundDim2(Bin a, int r, int min, float mg, float mdr) {return lowerBoundHH(a.hh, dim2, r, min, mg);}
	@Override public int upperBoundDim2(Bin a, int r, int max, float mg, float mdr) {return upperBoundHH(a.hh, dim2, r, max, mg);}

    // Dim 3 = CAGA
	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return lowerBoundCAGA(a.caga, dim3, r, min, mg);}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return upperBoundCAGA(a.caga, dim3, r, max, mg);}

    // Dim 4 = Depth 0
	@Override public int lowerBoundDim4(Bin a, int r, int min, float mg, float mdr) {return a.numDepths()<1 ? 0 : lowerBoundDepth(a.depth(0), dim4, r, min, mdr);}
	@Override public int upperBoundDim4(Bin a, int r, int max, float mg, float mdr) {return a.numDepths()<1 ? 0 : upperBoundDepth(a.depth(0), dim4, r, max, mdr);}

    // Dim 5 = Depth 1
	@Override public int lowerBoundDim5(Bin a, int r, int min, float mg, float mdr) {return a.numDepths()<2 ? 0 : lowerBoundDepth(a.depth(1), dim5, r, min, mdr);}
	@Override public int upperBoundDim5(Bin a, int r, int max, float mg, float mdr) {return a.numDepths()<2 ? 0 : upperBoundDepth(a.depth(1), dim5, r, max, mdr);}

	@Override public String toString() {return "("+gcLevel+","+dim2+","+dim3+","+dim4+","+dim5+")";}
}