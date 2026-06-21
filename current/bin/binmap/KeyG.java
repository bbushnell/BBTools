package bin.binmap;

import bin.Bin;

class KeyG extends Key {
	public KeyG() {}

	@Override
	//claim: KeyG = GC only (dim1); dim2/dim3 bounds return 0. Simplest non-null variant.
	public KeyG set(Bin a) {
		return setValue(a.gc());
	}

	public KeyG setLevel(int gcLevel_) {
		gcLevel=gcLevel_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		return this;
	}

	public KeyG setValue(float gc) {
		return setLevel(quantizeGC(gc));
	}

	@Override public int lowerBoundDim2(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim2(Bin a, int r, int max, float mg, float mdr) {return 0;}

	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return 0;}

	@Override
	public String toString() {return "("+gcLevel+")";}
}