package bin;

class Key0 extends Key {
	public Key0() {}
	
	@Override
	public Key0 set(Bin a) {
		return this; // No state to set
	}
	
	@Override
	public int lowerBoundDim1(float gc, int range, int minGrid, float maxGCDif) {return 0;}
	
	@Override
	public int upperBoundDim1(float gc, int range, int maxGrid, float maxGCDif) {return 0;}
	
	@Override public int lowerBoundDim2(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim2(Bin a, int r, int max, float mg, float mdr) {return 0;}
	
	@Override public int lowerBoundDim3(Bin a, int r, int min, float mg, float mdr) {return 0;}
	@Override public int upperBoundDim3(Bin a, int r, int max, float mg, float mdr) {return 0;}
	
	@Override
	public String toString() {return "()";}
}