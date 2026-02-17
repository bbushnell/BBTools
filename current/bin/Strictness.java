package bin;

import shared.Tools;

public class Strictness{
	
	public static float parseStrictness(String arg, String a, String b) {
		return toStrictness(arg, a, b);
	}
	
	static float toStrictness(String arg, String a, String b) {
		if(a.equalsIgnoreCase("strictness") || a.equalsIgnoreCase("stringency")) {
			int x=findStrictness(b);
			return x>=0 ? (float)values[x] : Float.parseFloat(b);
		}
		int x=findStrictness(arg);
//		assert(x>=0) : arg+", "+Arrays.toString(names);
		return x>=0 ? (float)values[x] : -1;
	}
	
	static int findStrictness(String term) {
		if(term==null) {return -1;}
		term=term.toLowerCase();
		int x=Tools.find(term, shortNames);
		if(x>=0) {return x;}
		return Tools.find(term, names);
	}

	public static final String[] shortNames={
		"is", "zs", "ys", //inf, zeta, yatta
		"xs", "hs", "us", "vs", //extra, hyper, ultra, very
		"s", "n", "l", //strict, normal, loose
		"vl", "ul", "hl", "xl", 
		"yl", "zl", "il"
	};
	public static final String[] names;
	public static final double[] values={
		0.25, 0.4, 0.5,
		0.6, 0.65, 0.7, 0.8,
		0.9, 1.0, 1.125,
		1.25, 1.375, 1.44, 1.5,
		2.0, 2.5, 5.0
	};
	
	static {
		names=new String[shortNames.length];
		for(int i=0; i<shortNames.length; i++) {//Order is important here!
			names[i]=shortNames[i].replace("s", "strict").replace("l", "loose").replace("n", "normal");
		}
	}
	
}
