package ml;

import java.util.Arrays;
import shared.Random;

import shared.Tools;

/**
 * Abstract base class for defining activation functions in neural network architectures.
 * Provides a standardized interface for implementing various activation function
 * transformations and their derivatives in machine learning models.
 * Supports multiple built-in activation function types including Sigmoid, Tanh,
 * Swish, RSLog, MSig, ExtendedSigmoid, ExtendedMSig, and Bell functions.
 *
 * @author Brian Bushnell
 * @date December 16, 2013
 */
public abstract class Function {
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Computes the activation function output for the given input value.
	 * @param x Input value to transform
	 * @return Activation function result
	 */
	public abstract double activate(double x);
	
	/**
	 * Calculates the derivative of the activation function with respect to input x.
	 * Used during backpropagation to compute gradients.
	 * @param x Input value at which to calculate derivative
	 * @return Derivative value at x
	 */
	public abstract double derivativeX(double x);
	
	/**
	 * Calculates the derivative given the function's output value.
	 * More efficient than derivativeX when the function output is already known.
	 * @param fx The pre-computed activation function output
	 * @return Derivative value corresponding to fx
	 */
	public abstract double derivativeFX(double fx);
	
	/**
	 * Calculates the derivative using both input and output values.
	 * Provides flexibility for functions that can optimize derivative computation
	 * when both input and output are available.
	 *
	 * @param x Input value
	 * @param fx Pre-computed activation function output
	 * @return Derivative value
	 */
	public abstract double derivativeXFX(double x, double fx);
	
	/** Returns the numeric type identifier for this activation function.
	 * @return Function type constant (SIG, TANH, RSLOG, etc.) */
	public abstract int type();
	
	/** Returns the string name identifier for this activation function.
	 * @return Function name string */
	public abstract String name();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts a string representation to a function type integer with validation.
	 * @param b String representation of function type
	 * @return Function type integer constant
	 */
	static final int toType(String b) {
		return toType(b, true);
	}

	/**
	 * Converts string representation to function type integer with optional validation.
	 * Searches both short and long type name arrays for matching function names.
	 *
	 * @param b String representation of function type
	 * @param assertValid Whether to assert that the type is valid
	 * @return Function type integer constant
	 */
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

	/**
	 * Normalizes the TYPE_RATES array to sum to 1.0 and computes cumulative rates.
	 * Creates the TYPE_RATES_CUM array for weighted random function selection.
	 * Must be called before using randomFunction() method.
	 * @return void
	 */
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
	
	/**
	 * Returns the singleton instance of the specified function type.
	 * @param type Function type constant
	 * @return Function instance for the specified type
	 */
	public static final Function getFunction(int type) {
		return functions[type];
	}
	
	/**
	 * Creates and initializes the array of function instances.
	 * Populates the functions array with singleton instances of all activation functions.
	 * Validates that each function has correct type and name assignments.
	 * @return Array of Function instances indexed by type
	 */
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
	
	/**
	 * Selects a random activation function based on weighted probabilities.
	 * Uses the TYPE_RATES_CUM array for weighted selection.
	 * @param randy Random number generator
	 * @return Randomly selected Function instance
	 */
	static final Function randomFunction(Random randy) {
		final int type=randomType(randy, TYPE_RATES_CUM);
		return functions[type];
	}
	
	/**
	 * Selects a random type index based on cumulative probability distribution.
	 * @param randy Random number generator
	 * @param cumRate Cumulative probability array
	 * @return Selected type index
	 */
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
	
	/** Array of singleton Function instances indexed by type */
	private static final Function[] functions=makeFunctions();
	
	/** Probability weights for random function selection */
	public static final float[] TYPE_RATES=new float[TYPES.length];
	//tanh=.4 sig=.6 msig=.02 rslog=.02 swish=0
	
	/** Cumulative probability distribution for weighted random selection */
	public static float[] TYPE_RATES_CUM=null;

	static {
		TYPE_RATES[TANH]=0.4f;
		TYPE_RATES[SIG]=0.6f;
		TYPE_RATES[MSIG]=0.02f;
		TYPE_RATES[RSLOG]=0.02f;
	}
}
