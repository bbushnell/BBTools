package test.json;

import java.util.Arrays;

import json.JsonObject;
import json.JsonParser;

/** Checks strict JSON validation and non-finite numeric serialization policies.
 * @author Jean */
public final class JsonPolicyTest{

	/** Runs all policy cases without depending on JVM assertions.
	 * @param args No arguments are accepted */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		strictNumbers();
		strictLiterals();
		strictStringEscapes();
		lenientParseControls();
		nonFiniteScalars();
		nonFiniteArrays();
		rawLiteral();
		System.out.println("JSON_POLICY_RESULT cases="+cases+" failures="+failures);
		if(failures!=0){throw new AssertionError("JSON policy failures="+failures);}
	}

	/** Standard validation accepts only JSON number grammar. */
	private static void strictNumbers(){
		final String[] valid={"0", "-0", "10", "0.1", "-12.34", "1e2", "1E-2", "1.25e+2"};
		for(final String text : valid){
			run("strict-number-valid-"+text, new Check(){public void check(){
				require(new JsonParser(text).validate(), "standard number rejected: "+text);
			}});
		}
		final String[] invalid={"+1", "01", "00", "1.", ".5", "1.e2", "1e", "1e+", "-"};
		for(final String text : invalid){
			run("strict-number-invalid-"+text, new Check(){public void check(){
				require(!new JsonParser(text).validate(), "non-standard number accepted: "+text);
			}});
		}
	}

	/** Standard validation accepts only lowercase JSON literals. */
	private static void strictLiterals(){
		final String[] valid={"true", "false", "null"};
		for(final String text : valid){
			run("strict-literal-valid-"+text, new Check(){public void check(){
				require(new JsonParser(text).validate(), "lowercase literal rejected: "+text);
			}});
		}
		final String[] invalid={"TRUE", "False", "t", "F", "NULL", "Null"};
		for(final String text : invalid){
			run("strict-literal-invalid-"+text, new Check(){public void check(){
				require(!new JsonParser(text).validate(), "non-standard literal accepted: "+text);
			}});
		}
	}

	/** Validation accepts the JSON escape set and rejects historical unknown escapes. */
	private static void strictStringEscapes(){
		final String[] valid={"\\\"", "\\\\", "\\/", "\\b", "\\f", "\\n", "\\r", "\\t",
				"\\u0041"};
		for(final String escape : valid){
			run("strict-escape-valid-"+escapeName(escape), new Check(){public void check(){
				require(new JsonParser(quote(escape)).validate(), "standard escape rejected");
			}});
		}
		final String[] invalid={"\\x", "\\'", "\\v", "\\0", "\\U0041"};
		for(final String escape : invalid){
			run("strict-escape-invalid-"+escapeName(escape), new Check(){public void check(){
				require(!new JsonParser(quote(escape)).validate(), "unknown escape accepted");
			}});
		}
		run("strict-escape-invalid-bad-unicode", new Check(){public void check(){
			require(!new JsonParser(quote("\\u12G4")).validate(), "bad Unicode escape accepted");
		}});
	}

	/** Parsing remains permissive for older callers while validate is strict. */
	private static void lenientParseControls(){
		run("lenient-legacy-literals", new Check(){public void check(){
			final JsonParser parser=new JsonParser("[T,F,Null,t,f]");
			final Object[] array=parser.parseJsonArray();
			require(!parser.errorState(), "permissive literal parse failed");
			require(Arrays.equals(array, new Object[]{true, false, null, true, false}),
					"permissive literals changed");
		}});
		run("lenient-legacy-numbers", new Check(){public void check(){
			final JsonParser parser=new JsonParser("[+1,01,1.,.5]");
			final Object[] array=parser.parseJsonArray();
			require(!parser.errorState(), "permissive number parse failed");
			require(Arrays.equals(array, new Object[]{1L, 1L, 1.0, 0.5}),
					"permissive numbers changed");
		}});
		run("lenient-unknown-escape", new Check(){public void check(){
			final JsonParser parser=new JsonParser(quote("\\x"));
			final JsonObject wrapper=new JsonParser("{\"v\":"+quote("\\x")+"}").parseJsonObject();
			require(!parser.validate(), "strict validation should reject the control escape");
			require("x".equals(wrapper.getString("v")), "permissive string escape changed");
		}});
		run("strict-validation-does-not-leak", new Check(){public void check(){
			final JsonParser parser=new JsonParser("[T]");
			require(!parser.validate(), "strict validation should reject legacy true");
			final Object[] array=parser.parseJsonArray("[T]");
			require(!parser.errorState(), "parse-with-input must reset validation error");
			require(Arrays.equals(array, new Object[]{true}), "legacy true parse changed");
		}});
	}

	/** Double, Float and formatted JsonLiteral output must map non-finite values to null. */
	private static void nonFiniteScalars(){
		checkNonFiniteValue("double-nan", Double.valueOf(Double.NaN));
		checkNonFiniteValue("double-pos-inf", Double.valueOf(Double.POSITIVE_INFINITY));
		checkNonFiniteValue("double-neg-inf", Double.valueOf(Double.NEGATIVE_INFINITY));
		checkNonFiniteValue("float-nan", Float.valueOf(Float.NaN));
		checkNonFiniteValue("float-pos-inf", Float.valueOf(Float.POSITIVE_INFINITY));
		checkNonFiniteValue("float-neg-inf", Float.valueOf(Float.NEGATIVE_INFINITY));
		run("restricted-nonfinite", new Check(){public void check(){
			try{
				JsonObject.setDecimals(3);
				checkNullObject("restricted-double-nan", new JsonObject("v", Double.NaN).toString());
				checkNullObject("restricted-float-inf",
						new JsonObject("v", Float.POSITIVE_INFINITY).toString());
			}finally{JsonObject.setDecimals(-1);}
		}});
		run("numeric-literal-nonfinite", new Check(){public void check(){
			final JsonObject o=new JsonObject();
			o.addLiteral("nan", Double.NaN, 2);
			o.addLiteral("posInf", Double.POSITIVE_INFINITY, 2);
			o.addLiteral("negInf", Double.NEGATIVE_INFINITY, 2);
			checkValidJson("numeric literal non-finite", o.toString());
			require(o.toString().contains("\"nan\": null"), "NaN literal must become null");
			require(o.toString().contains("\"posInf\": null"), "+Infinity literal must become null");
			require(o.toString().contains("\"negInf\": null"), "-Infinity literal must become null");
		}});
		run("numeric-literal-negative-decimals", new Check(){public void check(){
			final JsonObject o=new JsonObject();
			try{
				o.addLiteral("nan", Double.NaN, -1);
				throw new AssertionError("negative decimals bypassed validation");
			}catch(IllegalArgumentException expected){/* expected */}
		}});
	}

	/** Primitive float/double arrays must also route non-finite values through null. */
	private static void nonFiniteArrays(){
		run("double-array-nonfinite", new Check(){public void check(){
			checkNullArray(new JsonObject("v", new double[]{Double.NaN, Double.POSITIVE_INFINITY,
					Double.NEGATIVE_INFINITY, 1.5}).toString());
		}});
		run("float-array-nonfinite", new Check(){public void check(){
			checkNullArray(new JsonObject("v", new float[]{Float.NaN, Float.POSITIVE_INFINITY,
					Float.NEGATIVE_INFINITY, 1.5f}).toString());
		}});
	}

	/** Raw string literals stay a trusted unmodified escape hatch. */
	private static void rawLiteral(){
		run("raw-literal-nonfinite", new Check(){public void check(){
			final JsonObject o=new JsonObject();
			o.addLiteral("nan", "NaN");
			final String text=o.toString();
			require(text.contains("\"nan\": NaN"), "raw literal should be emitted verbatim: "+text);
			require(!new JsonParser(text).validate(), "raw NaN should not be made valid");
		}});
	}

	/** Checks one scalar member that should serialize as null.
	 * @param label Case label
	 * @param value Non-finite value */
	private static void checkNonFiniteValue(final String label, final Object value){
		run(label, new Check(){public void check(){
			checkNullObject(label, new JsonObject("v", value).toString());
		}});
	}

	/** Requires valid JSON with no non-finite tokens.
	 * @param label Case label
	 * @param text Serialized JSON */
	private static void checkValidJson(final String label, final String text){
		require(text.indexOf("NaN")<0, label+" emitted NaN: "+text);
		require(text.indexOf("Infinity")<0, label+" emitted Infinity: "+text);
		require(new JsonParser(text).validate(), label+" emitted invalid JSON: "+text);
	}

	/** Checks an object whose v member should parse back as null.
	 * @param label Case label
	 * @param text Serialized JSON object */
	private static void checkNullObject(final String label, final String text){
		checkValidJson(label, text);
		require(new JsonParser(text).parseJsonObject().containsKey("v"), label+" missing key");
		require(text.contains("\"v\": null"), label+" did not emit null: "+text);
	}

	/** Checks an object with an array value whose first three values are null.
	 * @param text Serialized JSON object */
	private static void checkNullArray(final String text){
		checkValidJson("primitive array non-finite", text);
		final Object[] array=new JsonParser(text).parseJsonObject().getArray("v");
		require(array.length==4, "primitive array length changed");
		require(array[0]==null && array[1]==null && array[2]==null, "non-finite values changed");
		require(Double.valueOf(1.5).equals(array[3]), "finite array value changed: "+array[3]);
	}

	/** Wraps a JSON string body in double quotes.
	 * @param body Raw JSON string body
	 * @return Quoted JSON string */
	private static String quote(final String body){return "\""+body+"\"";}

	/** Makes an escape printable as a case label.
	 * @param escape Escape text
	 * @return Label-safe text */
	private static String escapeName(final String escape){
		return escape.replace('\\', '_').replace('"', 'q').replace('\'', 's');
	}

	/** Runs one independent case and retains failures so baseline reports stay complete.
	 * @param label Case name
	 * @param check Case implementation */
	private static void run(final String label, final Check check){
		cases++;
		try{check.check();}
		catch(Throwable t){failures++; System.out.println("FAIL "+label+": "+t);}
	}

	/** Enforces a literal expectation with assertions enabled or disabled.
	 * @param condition Required condition
	 * @param detail Failure detail */
	private static void require(final boolean condition, final String detail){
		if(!condition){throw new AssertionError(detail);}
	}

	/** A finite independent test case. */
	private interface Check{
		/** Executes the case or throws on failure. */
		void check();
	}

	/** Prevents construction of this standalone fixture. */
	private JsonPolicyTest(){}
	/** Number of cases executed. */
	private static int cases;
	/** Number of failed cases. */
	private static int failures;
}
