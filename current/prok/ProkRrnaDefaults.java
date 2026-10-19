package prok;

/** Independent default methods; explicit rrna17 profiles override this policy.
 * @author Brian Bushnell, Raiden */
public final class ProkRrnaDefaults {

	static final Method DEFAULT_16S=Method.RRNA17;
	static final Method DEFAULT_23S=Method.RRNA17;
	static final Method DEFAULT_5S=Method.RRNA17;

	static int mask(boolean call16, boolean call23, boolean call5){
		return mask(DEFAULT_16S, DEFAULT_23S, DEFAULT_5S, call16, call23, call5);
	}
	public static int mask(Method s16, Method s23, Method s5, boolean call16, boolean call23, boolean call5){
		if(s16==null || s23==null || s5==null){throw new IllegalArgumentException("Every prokaryotic rRNA family requires an explicit default method");}
		return (call16 && s16==Method.RRNA17 ? S16 : 0)
			| (call23 && s23==Method.RRNA17 ? S23 : 0)
			| (call5 && s5==Method.RRNA17 ? S5 : 0);
	}
	static int bit(int type){
		if(type==ProkObject.r16S){return S16;}
		if(type==ProkObject.r23S){return S23;}
		if(type==ProkObject.r5S){return S5;}
		throw new IllegalArgumentException("Default method selection is limited to prokaryotic16S/23S/5S, type="+type);
	}
	public enum Method {PGM, RRNA17}
	public static final int S16=1, S23=2, S5=4;
}
