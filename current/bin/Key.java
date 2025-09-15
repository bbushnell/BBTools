package bin;

import shared.Tools;

class Key implements Cloneable {
	
	static boolean parse(String arg, String a, String b) {
		
		if(a.equalsIgnoreCase("gcwidth")){
			float f=Float.parseFloat(b);
			setGCWidth(f);
		}else if(a.equalsIgnoreCase("gcmult")){
			float f=Float.parseFloat(b);
			setGCMult(f);
		}else if(a.equalsIgnoreCase("depthwidth")){
			float f=Float.parseFloat(b);
			setDepthWidth(f);
		}else if(a.equalsIgnoreCase("depthmult")){
			float f=Float.parseFloat(b);
			setDepthMult(f);
		}else {
			return false;
		}
		return true;
	}
	
	public Key(float gc, float cov, float cov2) {
		setValue(gc, cov, cov2);
	}
	
	public Key() {}

	public Key set(Bin a) {
		return setValue(a.gc(), a.depth(0), a.depth(1));
	}
	
	public Key setLevel(int gcLevel_, int covLevel_, int covLevel2_) {
		gcLevel=gcLevel_;
		covLevel=covLevel_;
		covLevel2=covLevel2_;
		assert(gcLevel>=0 && gcLevel<=(int)gcLevelMult);
		assert(covLevel>=0 && covLevel<=maxDepthLevel) : covLevel+", "+maxDepthLevel;
		assert(covLevel2>=0 && covLevel2<=maxDepthLevel) : covLevel2+", "+maxDepthLevel;
		return this;
	}
	
	public Key setValue(float gc, float cov, float cov2) {
		assert(gc>=0 && gc<=1) : gc;
		assert(cov>=0) : cov;
		assert(cov2>=0) : cov;
		return setLevel(quantizeGC(gc), quantizeDepth(cov), quantizeDepth(cov2));
	}
	
	@Override
	public boolean equals(Object other) {
		return equals((Key)other);
	}
	
	public boolean equals(Key b) {
		return gcLevel==b.gcLevel && covLevel==b.covLevel && covLevel2==b.covLevel2;
	}
	
	@Override
	public int hashCode() {
		return covLevel+(covLevel2<<10)+(gcLevel<<20);
	}
	
	@Override
	public Key clone() {
		try {
			return (Key)(super.clone());
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	//This is probably faster but not as simple to explain or adjust
//	public static int quantizeDepth(float depth) {
//		float yf=depth*depth*16;
//		long y=(long)yf;
//		int zeros=Long.numberOfLeadingZeros(y);
//		int level=Tools.min(maxDepthLevel, 64-zeros);
//		return level;
//	}

	public static int quantizeDepth(float depth) {
		depth=Tools.min(depth, maxDepth);
		float yf=((float)(Tools.log2(depth+0.0625f)+4));
		int level=(int)(yf*depthLevelMult);
		return level;
	}

	public static int quantizeGC(float gc) {
		return (int)(Tools.mid(0,gc,1)*gcLevelMult);
	}
	
	public String toString() {
		return "("+gcLevel+","+covLevel+","+covLevel2+")";
	}
	
	static final void setGCMult(float f) {
		assert(f>=2);
		setGCWidth(1/f);
	}
	
	static final void setGCWidth(float f) {
		assert(f>0 && f<=0.5f);
		gcLevelWidth=f;
		gcLevelMult=1f/gcLevelWidth;
	}
	
	static final void setDepthWidth(float f) {
		assert(f>0);
		setDepthMult(1/f);
	}
	
	static final void setDepthMult(float f) {
		assert(f>0);
		depthLevelMult=f;
		maxDepthLevel=quantizeDepth(maxDepth);
		assert(maxDepthLevel>0) : "maxDepthLevel="+maxDepthLevel+", depthLevelMult="+depthLevelMult+
			", maxDepth="+maxDepth+", yf="+(Tools.log2(maxDepth+0.0625f)+4);
	}
	
	int gcLevel;
	int covLevel;
	int covLevel2;
	
	//gcwidth=0.01, depthwidth=0.25 seems faster, more sensitive, and more specific (halving both of them).
	private static final float maxDepth=1000000;
	private static float depthLevelMult=2f;
	private static float gcLevelWidth=0.02f;
	private static float gcLevelMult=1f/gcLevelWidth;
	private static int maxDepthLevel=quantizeDepth(maxDepth);
	
}
