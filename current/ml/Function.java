package ml;

import java.util.Arrays;
import java.util.Random;

import shared.Tools;

public abstract class Function {
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	public abstract double activate(double x);
	
	public abstract double derivativeX(double x);
	
	public abstract double derivativeFX(double fx);
	
	public abstract double derivativeXFX(double x, double fx);
	
	public abstract int type();
	
	public abstract String name();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	static final int toType(String b) {
		return toType(b, true);
	}

	static final int toType(String b, boolean assertValid) {
			int type;
			if(Tools.startsWithLetter(b)) {
				type=Tools.findIC(b, TYPES);
				if(type<0) {type=Tools.findIC(b, TYPES_LONG);}
			}else{
				type=Integer.parseInt(b);
				throw new RuntimeException("Numbers are not allowed for defining types: "+b);
			}
			assert(!assertValid || type>=0 && type<TYPES.length) : type;
	//		System.err.println(b+" -> "+type);
			return type;
		}

	public static synchronized final void normalizeTypeRates() {
		assert(TYPE_RATES_CUM==null);
		double sum=shared.Vector.sum(TYPE_RATES);
		assert(sum>=0) : sum;
		
		if(sum<=0) {
			TYPE_RATES_CUM=null;
			return;
		}
		if(Tools.absdif(sum, 1)>0.000001){
			final double mult=1.0/sum;
			for(int i=0; i<TYPE_RATES.length; i++) {
				double r=TYPE_RATES[i];
				assert(r>=0) : i+": "+r;
				TYPE_RATES[i]=(float)(r*mult);
			}
		}
		
		TYPE_RATES_CUM=new float[TYPE_RATES.length];
		double c=0;
		for(int i=0; i<TYPE_RATES.length; i++) {
			double r=TYPE_RATES[i];
			assert(r>=0) : i+": "+r;
			c+=r;
			TYPE_RATES_CUM[i]=(float)c;
		}
		assert(Tools.absdif(c, 1)<0.00001);
		TYPE_RATES_CUM[TYPE_RATES_CUM.length-1]=1;
	}
	
	/*--------------------------------------------------------------*/
	
	public static final Function getFunction(int type) {
		return functions[type];
	}
	
	private static final Function[] makeFunctions() {
		assert(functions==null);
		Function[] array=new Function[TYPES.length];
		array[SIG]=Sigmoid.instance;
		array[TANH]=Tanh.instance;
		array[RSLOG]=RSLog.instance;
		array[MSIG]=MSig.instance;
		array[SWISH]=Swish.instance;
		array[ESIG]=ExtendedSigmoid.instance;
		array[EMSIG]=ExtendedMSig.instance;
		array[BELL]=Bell.instance;
		for(int i=0; i<array.length; i++) {
			Function f=array[i];
			assert(f!=null) : i+", "+TYPES[i]+", "+f;
			assert(f.type()==i) : i+", "+TYPES[i]+", "+f;
			assert(f.name().equals(TYPES[i])) : i+", "+f;
		}
		return array;
	}
	
	/*--------------------------------------------------------------*/
	
	static final Function randomFunction(Random randy) {
		final int type=randomType(randy, TYPE_RATES_CUM);
		return functions[type];
	}
	
	static final int randomType(Random randy, float[] cumRate) {
		float f=randy.nextFloat();
		for(int i=0; i<cumRate.length; i++) {
			if(cumRate[i]>=f) {return i;}
		}
		assert(false) : f+", "+Arrays.toString(cumRate);
		return cumRate.length-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public static final int SIG=0, TANH=1, RSLOG=2, MSIG=3, SWISH=4, ESIG=5, EMSIG=6, BELL=7;
	
	static final String[] TYPES=new String[] {"SIG", "TANH", "RSLOG", "MSIG", "SWISH", "ESIG", "EMSIG", "BELL"};
	
	static final String[] TYPES_LONG=new String[] {"SIGMOID", "HYPERBOLICTANGENT", 
	"ROTATIONALLYSYMMETRICLOGARITHM", "MIRROREDSIGMOID", "SWISH",
	"EXTENDEDSIGMOID", "EXTENDEDMIRROREDSIGMOID", "GAUSSIAN"};
	
	private static final Function[] functions=makeFunctions();
	
	public static final float[] TYPE_RATES=new float[TYPES.length];
	//tanh=.4 sig=.6 msig=.02 rslog=.02 swish=0
	
	public static float[] TYPE_RATES_CUM=null;

	static {
		TYPE_RATES[TANH]=0.4f;
		TYPE_RATES[SIG]=0.6f;
		TYPE_RATES[MSIG]=0.02f;
		TYPE_RATES[RSLOG]=0.02f;
	}
}
