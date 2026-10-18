package test.ml;

import ml.VectorUtils;

import java.lang.reflect.Field;
import java.util.ArrayList;

/** Constructor-only checks for output fraction normalization; no files are opened.
 * @author Brian Bushnell, Shinobu */
public final class VectorUtilsFractionTest{

	/** Runs one named case in a fresh JVM with assertions enabled. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		final String name=args[0];
		if(name.equals("large_equal")){
			check("a:3e38,b:3e38", new float[]{0.5f, 0.5f});
		}else if(name.equals("positive_infinity")){
			reject("a:Infinity,b:1", "finite");
		}else if(name.equals("negative_infinity")){
			reject("a:-Infinity,b:1", "finite");
		}else if(name.equals("nan")){
			reject("a:NaN,b:1", null);
		}else if(name.equals("ratio")){
			check("a:1,b:3", new float[]{0.25f, 0.75f});
		}else if(name.equals("unspecified")){
			check("a,b", new float[]{0.5f, 0.5f});
		}else if(name.equals("mixed")){
			check("a:0.25,b", new float[]{0.25f, 0.75f});
		}else if(name.equals("negative_unspecified")){
			check("a:-2,b:0.25", new float[]{0.75f, 0.25f});
		}else if(name.equals("zero_total")){
			reject("a:0,b:0", "sum to 0");
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		System.err.println("VectorUtilsFractionTest PASS "+name);
	}

	/** Reads normalized constructor state without starting process(), readers, writers or workers. */
	private static void check(final String specification, final float[] expected) throws Exception{
		final VectorUtils tool=new VectorUtils(new String[]{"in=unused.tsv", "out="+specification});
		final Field field=VectorUtils.class.getDeclaredField("outFracs");
		field.setAccessible(true);
		final ArrayList<?> actual=(ArrayList<?>)field.get(tool);
		if(actual.size()!=expected.length){throw new AssertionError("Unexpected fraction count: "+actual.size());}
		for(int i=0; i<expected.length; i++){
			final float value=((Float)actual.get(i)).floatValue();
			if(Float.floatToIntBits(value)!=Float.floatToIntBits(expected[i])){
				throw new AssertionError("Fraction "+i+" expected "+expected[i]+" but observed "+value);
			}
		}
	}

	/** Requires constructor rejection; checks a stable diagnostic fragment when specified. */
	private static void reject(final String specification, final String diagnostic){
		try{
			new VectorUtils(new String[]{"in=unused.tsv", "out="+specification});
		}catch(AssertionError e){
			if(diagnostic!=null && (e.getMessage()==null || !e.getMessage().contains(diagnostic))){throw e;}
			return;
		}
		throw new AssertionError("Accepted invalid output fractions: "+specification);
	}

	private VectorUtilsFractionTest(){}
}
