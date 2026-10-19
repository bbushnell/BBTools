package test.prok;

import java.nio.charset.StandardCharsets;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;

/** Actual-alignment checks that each trim gate preserves an accepted input span.
 * @author Raiden */
public final class NcrnaTrimGuardTest {
	public static void main(String[] args){
		NcrnaScavengerSoftMaskTest.main(new String[0]);
		final byte[] consensus=bytes("ACGTTGCAAGTCGATCGTACGATGC");
		final byte[] exact=bytes("TTTTTTTTTTACGTTGCAAGTCGATCGTACGATGCGGGGGGGGGG");
		final byte[] mismatch=exact.clone();mismatch[20]='A';
		check(consensus,exact,10,.75f,1f,true,"equal-quality trim");
		check(consensus,exact,consensus.length+1,.75f,1f,false,"minimum length");
		check(consensus,mismatch,10,1f,.50f,false,"family identity gate");
		check(consensus,mismatch,10,.50f,1f,false,"identity degradation");
		check(consensus,exact,10,.75f,.90f,true,"improved identity");
		boolean rejected=false;
		try{check(consensus,exact,10,.75f,Float.NaN,false,"missing identity");}
		catch(IllegalArgumentException expected){rejected=true;}
		if(!rejected){throw new AssertionError("Missing verifier identity must fail before endpoint mutation");}
		System.out.println("PASS NcrnaTrimGuardTest length identity degradation equal improvement invalid");
	}
	private static void check(byte[] consensus,byte[] bases,int minLen,float pass,float previous,boolean expected,String label){
		assert(consensus.length>0 && bases.length>36):"Fixture has one internal model and an explicit pre-trim span";
		final NcrnaScavenger s=new NcrnaScavenger(new byte[][]{consensus},null,null,null,17,minLen,10,
			7,1,false,0f,0f,0f,1,0f,1f,pass,.50f);
		final Orf o=new Orf(label,8,36,0,0,bases,false,ProkObject.RNA);
		final boolean trimmed=s.trimToAlignmentExtent(o,bases,0,8,36,previous);
		if(trimmed!=expected){throw new AssertionError(label+": expected trim="+expected+" observed="+trimmed);}
		if(expected){
			if(o.start!=10 || o.stop!=10+consensus.length-1){throw new AssertionError(label+": exact internal consensus bounds lost");}
		}else if(o.start!=8 || o.stop!=36){throw new AssertionError(label+": rejected trim changed the verified span");}
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
}
