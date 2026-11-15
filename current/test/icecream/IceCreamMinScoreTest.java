package test.icecream;

import java.nio.charset.StandardCharsets;

import aligner.AlignmentResult;
import icecream.IceCreamAligner;
import shared.Shared;

/** Tests full-alignment threshold width with a literal two-base matching alignment.
 * @author Brian Bushnell, Sucrose
 */
public final class IceCreamMinScoreTest{

	/** Requires the real native backend by default; backend=java explicitly checks Java only.
	 * @param args Optional backend=both or backend=java
	 */
	public static void main(String[] args){
		boolean nativeRequired=true;
		for(String arg : args){
			if(arg.equals("backend=java")){nativeRequired=false;}
			else if(arg.equals("backend=both")){nativeRequired=true;}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		require(IceCreamAligner.class.desiredAssertionStatus(), "Production assertions must be enabled");
		final boolean originalJNI=Shared.USE_JNI;
		int failures=0;
		try{
			System.out.println("backend\tmin_score\texpected\tobserved\tresult");
			Shared.USE_JNI=false;
			final IceCreamAligner java=IceCreamAligner.makeAligner(32);
			require(java.getClass().getName().equals("icecream.IceCreamAlignerJava"),
					"Java control must use the production Java implementation");
			failures+=runCases("java", java);
			if(nativeRequired){
				Shared.USE_JNI=true;
				final IceCreamAligner jni=IceCreamAligner.makeAligner(32);
				require(jni.getClass().getName().equals("icecream.IceCreamAlignerJNI"),
						"Native cases require JNI; silently falling back to Java would not test threshold narrowing");
				failures+=runCases("jni", jni);
			}
		}finally{
			Shared.USE_JNI=originalJNI;
		}
		require(failures==0, failures+" full-alignment threshold observations differ from the literal oracle");
	}

	/** Includes both short boundaries, both adjacent out-of-range values and an ordinary control. */
	private static int runCases(String backend, IceCreamAligner aligner){
		int failures=0;
		failures+=observe(backend, aligner, -32769, "2:2:1:1.0");
		failures+=observe(backend, aligner, -32768, "2:2:1:1.0");
		failures+=observe(backend, aligner, -1, "2:2:1:1.0");
		failures+=observe(backend, aligner, 32767, "null");
		failures+=observe(backend, aligner, 32768, "null");
		return failures;
	}

	/** Literal score2, query endpoint2 and reference endpoint1 follow from the exact AC/AC match. */
	private static int observe(String backend, IceCreamAligner aligner, int minScore, String expected){
		final byte[] query="AC".getBytes(StandardCharsets.US_ASCII);
		final byte[] reference="AC".getBytes(StandardCharsets.US_ASCII);
		final AlignmentResult result=aligner.alignForward(query, reference, 0, 1, minScore, 0.9f);
		final String observed=result==null ? "null" : result.maxScore+":"+result.maxQpos+":"+result.maxRpos+":"+result.ratio;
		final boolean pass=expected.equals(observed);
		System.out.println(backend+"\t"+minScore+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Keeps fixture/backend identity and oracle checks active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Static-only fixture. */
	private IceCreamMinScoreTest(){}
}
