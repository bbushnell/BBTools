package test.json;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import json.JsonObject;
import json.JsonParser;

/** Checks JSON text, object mutation, and structural validation against literal oracles.
 * @author Shinobu */
public final class JsonRegressionTest{

	/** Runs all finite cases without depending on JVM assertions.
	 * @param args No arguments are accepted */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		unicode();
		encodingValidation();
		objects();
		validation();
		controls();
		formatting();
		System.out.println("JSON_RESULT cases="+cases+" failures="+failures);
		if(failures!=0){throw new AssertionError("JSON regression failures="+failures);}
	}

	/** Distinguishes genuine replacement characters from malformed UTF-8 bytes. */
	private static void encodingValidation(){
		final byte[][] malformed={
				{'[', '"', (byte)0xc3, '"', ']'},
				{'[', '"', (byte)0xc0, (byte)0xaf, '"', ']'},
				{'[', '"', (byte)0xed, (byte)0xa0, (byte)0x80, '"', ']'}};
		for(final byte[] bytes : malformed){
			run("malformed-utf8", new Check(){public void check(){
				require(!new JsonParser(bytes).validate(), "malformed UTF-8 must set validation error");
			}});
		}
		run("literal-replacement-character", new Check(){public void check(){
			final JsonParser parser=new JsonParser("[\"\ufffd\"]");
			require(parser.validate(), "a genuinely encoded U+FFFD is valid");
			require("\ufffd".equals(parser.parseJsonArray()[0]), "validation must preserve cursor");
		}});
		run("mixed-raw-and-escaped-unicode", new Check(){public void check(){
			final Object[] values=new JsonParser("[\"\u00e9\\u4e2d\ud83d\ude00\"]").parseJsonArray();
			require("\u00e9\u4e2d\ud83d\ude00".equals(values[0]), "raw UTF-8 spans around escapes");
		}});
	}

	/** Checks numeric formatting beyond the byte builder's safe long/precision range. */
	private static void formatting(){
		run("large-finite-number", new Check(){public void check(){
			try{
				JsonObject.setDecimals(3);
				final JsonObject o=new JsonObject("v", 1e30);
				final JsonObject parsed=new JsonParser(o.toString()).parseJsonObject();
				require(parsed.getDouble("v")==1e30, "large finite value must not saturate to Long.MAX_VALUE");
			}finally{JsonObject.setDecimals(-1);}
		}});
		run("high-precision-number", new Check(){public void check(){
			try{
				JsonObject.setDecimals(30);
				require(new JsonObject("v", 1.25).toString().contains("1.250000000000000000000000000000"),
						"high precision must not overflow the fractional accumulator");
			}finally{JsonObject.setDecimals(-1);}
		}});
	}

	/** Exercises raw UTF-8 and escaped UTF-16 independently of serialization. */
	private static void unicode(){
		final String[] values={"caf\u00e9", "\u03bb\u4e2d", "\ud83d\ude00", "a\nb\t\"\\"};
		final String[] escaped={"caf\\u00e9", "\\u03bb\\u4e2d", "\\ud83d\\ude00", "a\\nb\\t\\\"\\\\"};
		for(int i=0; i<values.length; i++){
			final String expected=values[i], token=escaped[i];
			run("escaped-value-"+i, new Check(){public void check(){
				require(expected.equals(new JsonParser("{\"v\":\""+token+"\"}").parseJsonObject().getString("v")), "decoded value");
			}});
			run("escaped-key-"+i, new Check(){public void check(){
				require("yes".equals(new JsonParser("{\""+token+"\":\"yes\"}").parseJsonObject().getString(expected)), "decoded key");
			}});
			if(i<3){
				run("utf8-value-"+i, new Check(){public void check(){
					final byte[] text=("{\"v\":\""+expected+"\"}").getBytes(StandardCharsets.UTF_8);
					require(expected.equals(new JsonParser(text).parseJsonObject().getString("v")), "UTF-8 value");
				}});
			}
			run("serialized-value-"+i, new Check(){public void check(){
				final String text=new JsonObject("v", expected).toString();
				require(text.contains("\""+token+"\""), "ASCII-safe JSON quoting: "+text);
			}});
		}
	}

	/** Checks last-write-wins and unique renaming across both internal maps. */
	private static void objects(){
		run("scalar-to-object", new Check(){public void check(){
			final JsonObject o=new JsonObject("k", "old"), child=new JsonObject("n", 3L);
			o.add("k", child);
			require(o.getJson("k")==child && o.omapSize()==0 && o.jmapSize()==1, "one nested value");
		}});
		run("object-to-scalar", new Check(){public void check(){
			final JsonObject o=new JsonObject("k", new JsonObject("n", 3L));
			o.add("k", "new");
			require("new".equals(o.getString("k")) && o.getJson("k")==null && o.jmapSize()==0, "one scalar value");
		}});
		run("object-to-literal", new Check(){public void check(){
			final JsonObject o=new JsonObject("k", new JsonObject("n", 3L));
			o.addLiteral("k", "7");
			require(o.jmapSize()==0 && o.omapSize()==1 && o.toString().contains("\"k\": 7"), "literal replacement");
		}});
		run("rename-across-maps", new Check(){public void check(){
			final JsonObject o=new JsonObject("k", 1L), child=new JsonObject("n", 2L);
			o.addAndRename("k", child);
			o.addAndRename("k", "third");
			require(o.getLong("k")==1L && o.getJson("k 2")==child && "third".equals(o.getString("k 3")), "cross-map suffixes");
		}});
		run("null-child", new Check(){public void check(){
			final JsonObject o=new JsonObject();
			o.add("k", (JsonObject)null);
			require(o.containsKey("k") && o.toString().contains("\"k\": null"), "JSON null value");
		}});
		run("duplicate-input-key", new Check(){public void check(){
			final JsonObject o=new JsonParser("{\"k\":1,\"k\":{\"n\":2}}").parseJsonObject();
			require(o.omapSize()==0 && o.getJson("k").getLong("n")==2L, "last duplicate wins across types");
		}});
	}

	/** Checks structural rules that do not depend on optional token extensions. */
	private static void validation(){
		final String[] valid={"{}", "[]", " \n{\"x\":1}\t", "[true,false,null,{\"a\":[1,2]}]",
				"null", "true", "-1.25e2", "\"string\""};
		for(final String text : valid){
			run("valid-"+text, new Check(){public void check(){
				require(new JsonParser(text).validate(), "valid document");
			}});
		}
		final String[] invalid={"{\"a\" 1}", "{\"a\":1}garbage", "[1 2]", "[1,]", "{\"a\":1,}",
				"{\"a\":\"one\" \"b\":\"two\"}", "{\"a\":\"line\nfeed\"}", "[", "{\"a\":1",
				"", " \t\n", "[\"\\uXY12\"]", "[\"unterminated]"};
		for(final String text : invalid){
			run("invalid-"+text.replace('\n', ' '), new Check(){public void check(){
				require(!new JsonParser(text).validate(), "malformed document must fail validation");
			}});
		}
		run("repeat-validation", new Check(){public void check(){
			final JsonParser p=new JsonParser("{\"a\":1}");
			require(p.validate() && p.validate(), "validation must not depend on a prior validation call");
		}});
	}

	/** Keeps ordinary nested values, primitive arrays and literal formatting covered. */
	private static void controls(){
		run("nested-values", new Check(){public void check(){
			final JsonObject o=new JsonParser("{\"v\":[1,2.5,true,null,{\"n\":\"ok\"}]}").parseJsonObject();
			final Object[] a=o.getArray("v");
			require(a.length==5 && a[0].equals(1L) && a[1].equals(2.5) && a[2].equals(true) && a[3]==null
					&& "ok".equals(((JsonObject)a[4]).getString("n")), "nested types");
		}});
		run("integer-array", new Check(){public void check(){
			final JsonObject o=new JsonObject("v", new int[]{1, -2, 3});
			require(o.toString().contains("[1, -2, 3]"), "primitive array serialization");
		}});
		run("empty-array", new Check(){public void check(){
			require("[]".equals(JsonObject.toString(new Object[0])), "empty array");
		}});
		run("ordinary-array", new Check(){public void check(){
			require(Arrays.equals(new JsonParser("[1,2,3]").parseJsonArray(), new Object[]{1L, 2L, 3L}), "array integers");
		}});
		run("finite-literal", new Check(){public void check(){
			final JsonObject o=new JsonObject();
			o.addLiteral("v", 3.125, 2);
			require(o.toString().contains("\"v\": 3.13"), "formatted numeric literal");
		}});
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
	private JsonRegressionTest(){}
	/** Number of cases executed. */
	private static int cases;
	/** Number of failed cases. */
	private static int failures;
}
