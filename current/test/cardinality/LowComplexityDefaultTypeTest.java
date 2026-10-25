package test.cardinality;

import java.lang.reflect.Field;

import cardinality.DDLCalibrationDriver;
import cardinality.LowComplexityCalibrationDriver;
import parse.Parser;

/** Checks lowcomplexcalibrate default type normalization before CF initialization.
 * @author Sucrose
 */
public final class LowComplexityDefaultTypeTest{

	/** Runs the captured baseline or corrected effective-type contract. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation: "+args[0]);
		final String initialType=Parser.loglogType;
		try{
			checkCase(fixed, null, fixed);
			checkCase(true, "type=DDL", true);
		}finally{
			Parser.loglogType=initialType;
		}
	}

	/** Verifies whether the effective type reaches DDL's mantissa-specific V3 columns. */
	private static void checkCase(final boolean expectMeanM, final String typeArg,
			final boolean fixed) throws Exception{
		Parser.loglogType="DDL";
		if(typeArg==null){LowComplexityCalibrationDriver.main(ARGS);}
		else{
			final String[] args=new String[ARGS.length+1];
			System.arraycopy(ARGS, 0, args, 0, ARGS.length);
			args[ARGS.length]=typeArg;
			LowComplexityCalibrationDriver.main(args);
		}
		final boolean hasMeanM=contains(v3ColumnNames(), "MeanM_cf");
		check(hasMeanM==expectMeanM,
			"Expected MeanM_cf presence "+expectMeanM+" for "+typeArg+
			"; fixed="+fixed+" observed="+hasMeanM+" Parser.loglogType="+Parser.loglogType);
		System.out.println("PASS\t"+(typeArg==null ? "default" : typeArg)+"\tmeanm="+hasMeanM);
	}

	/** Returns DDLCalibrationDriver's active V3 column names. */
	private static String[] v3ColumnNames() throws Exception{
		final Field f=DDLCalibrationDriver.class.getDeclaredField("V3_COL_NAMES");
		f.setAccessible(true);
		return (String[])f.get(null);
	}

	/** Returns true if the array contains the exact string. */
	private static boolean contains(final String[] array, final String value){
		for(String s : array){
			if(value.equals(s)){return true;}
		}
		return false;
	}

	/** Keeps oracle checks active under both assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String[] ARGS={
		"card=1", "ddls=1", "iter=1", "buckets=64",
		"k=1", "t=1", "reportfrac=1", "minrand=f"
	};
}
