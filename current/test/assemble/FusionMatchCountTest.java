package test.assemble;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import assemble.FusionJoinDiagnostic;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/** Exact-count controls, a cost pilot and a fixed overflow boundary for the real match helper.
 * Reflection preserves its package-private API; no counter state is supplied or changed.
 * Reference entries are explicit forward/reverse pairs, with immutable aliases in the
 * larger case. This tests helper counts/last coordinates, not reference loading,
 * biological label accuracy, the indexed labeler, or whole-program performance.
 * A separate tiny consumer case checks seedless labeler saturation against literal counts.
 * @author Brian Bushnell, Sucrose
 */
public final class FusionMatchCountTest{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepts one fixed pilot, boundary or consumer case plus launcher tokens; no custom size. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FusionMatchCountTest.class, false);
		try{
			boolean selected=false, boundary=false, consumer=false;
			for(String arg : pp.args){
				if(Parser.isJavaFlag(arg)){continue;}
				if(!selected && (arg.equalsIgnoreCase("case=pilot") || arg.equalsIgnoreCase("case=boundary") ||
						arg.equalsIgnoreCase("case=consumer"))){
					selected=true; boundary=arg.equalsIgnoreCase("case=boundary"); consumer=arg.equalsIgnoreCase("case=consumer");
				}else{throw new IllegalArgumentException("Expected one case=pilot, case=boundary or case=consumer option: "+arg);}
			}
			if(!selected){throw new IllegalArgumentException("Required case=pilot, case=boundary or case=consumer.");}
			final FusionMatchCountTest test=new FusionMatchCountTest();
			if(consumer){test.consumer();}else{test.run(boundary);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves the actual method and result fields once; count may be int or long. */
	private FusionMatchCountTest() throws ReflectiveOperationException{
		find=FusionJoinDiagnostic.class.getDeclaredMethod("find", ArrayList.class, byte[].class);
		find.setAccessible(true);
		final Class<?> hit=Class.forName("assemble.FusionJoinDiagnostic$Hit");
		count=hit.getDeclaredField("count"); sequence=hit.getDeclaredField("sequence"); pos=hit.getDeclaredField("pos");
		count.setAccessible(true); sequence.setAccessible(true); pos.setAccessible(true);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs five controls plus the pilot, or one boundary call; each invocation uses fresh state. */
	private void run(final boolean boundary) throws ReflectiveOperationException{
		System.out.println("case\treference_sequences\tquery\texpected_count\tobserved_count\texpected_sequence\tobserved_sequence\texpected_start\tobserved_start\tpassed");
		int failed=0;
		if(!boundary){
			if(!observe("absent", pair("CCCCC", "GGGGG"), "AAA", 0, -1, -1)){failed++;}
			if(!observe("forward-unique", pair("CCAAACC", "GGTTTGG"), "AAA", 1, 0, 2)){failed++;}
			if(!observe("reverse-unique", pair("GGTTTGG", "CCAAACC"), "AAA", 1, 1, 2)){failed++;}
			if(!observe("overlapping", pair("AAAA", "TTTT"), "AAA", 2, 0, 1)){failed++;}
			if(!observe("palindromic-pair", pair("CATG", "CATG"), "AT", 2, 1, 1)){failed++;}
		}
		final char[] bases=new char[65538];
		Arrays.fill(bases, 'A'); final String forward=new String(bases);
		Arrays.fill(bases, 'T'); final String reverse=new String(bases);
		final int pairs=boundary ? BOUNDARY_PAIRS : PILOT_PAIRS;
		final ArrayList<String> refs=new ArrayList<String>(2*pairs);
		for(int i=0; i<pairs; i++){refs.add(forward); refs.add(reverse);}
		// Every A string has 65536 overlapping AAA starts; the paired T string has none.
		final long expected=pairs*65536L;
		assert(expected==(boundary ? 2147483648L : 67108864L)) : "Fixed cases must retain their independent exact-count oracles.";
		if(!observe(boundary ? "paired-alias-boundary" : "paired-alias-pilot", refs, "AAA", expected, 2*(pairs-1), 65535)){failed++;}
		check(refs.size()==2*pairs, "Exhaustive lookup must preserve the caller's reference-list size.");
		for(int i=0; i<pairs; i++){
			check(refs.get(2*i)==forward && refs.get(2*i+1)==reverse,
					"Exhaustive lookup must retain borrowed forward/reverse String identities at pair "+i);
		}
		final int cases=boundary ? 1 : 6;
		System.out.println("FUSION_MATCH_COUNT_TEST case="+(boundary ? "boundary" : "pilot")+" assertions="+
				FusionMatchCountTest.class.desiredAssertionStatus()+" cases="+cases+" failed="+failed);
		System.out.flush();
		if(failed!=0){throw new AssertionError("Fusion match-count contract failed: "+failed+" of "+cases);}
	}

	/** Invokes production once, prints independent expectations and observations, and returns parity. */
	private boolean observe(final String name, final ArrayList<String> refs, final String query,
			final long expected, final int expectedSequence, final int expectedStart) throws ReflectiveOperationException{
		assert(!query.isEmpty()) : "FusionJoinDiagnostic.find requires a nonempty literal query.";
		final Object hit;
		try{hit=find.invoke(null, refs, query.getBytes(StandardCharsets.US_ASCII));}
		catch(InvocationTargetException e){throw new IllegalStateException("Match helper threw in "+name, e.getCause());}
		final Object value=count.get(hit);
		check(value instanceof Number, "Hit.count must retain a numeric field type for before/after comparison.");
		final long observed=((Number)value).longValue();
		final int observedSequence=sequence.getInt(hit), observedStart=pos.getInt(hit);
		final boolean passed=observed==expected && observedSequence==expectedSequence && observedStart==expectedStart;
		System.out.println(name+"\t"+refs.size()+"\t"+query+"\t"+expected+"\t"+observed+"\t"+
				expectedSequence+"\t"+observedSequence+"\t"+expectedStart+"\t"+observedStart+"\t"+passed);
		System.out.flush();
		return passed;
	}

	/** Checks three literal overlapping starts and the seedless labeler's cap of two.
	 * Saturated labeler coordinates are intentionally not part of this contract. */
	private void consumer() throws ReflectiveOperationException{
		final String forward="AAAAA", reverse="TTTTT";
		final ArrayList<String> refs=pair(forward, reverse);
		final byte[] query="AAA".getBytes(StandardCharsets.US_ASCII);
		final ArrayList<byte[]> patterns=new ArrayList<byte[]>();
		patterns.add(query);
		final Method findAll=Class.forName("assemble.FusionJoinLabeler").getDeclaredMethod("findAll", ArrayList.class, ArrayList.class);
		findAll.setAccessible(true);
		final Object exact;
		final Object[] labels;
		try{
			exact=find.invoke(null, refs, query);
			labels=(Object[])findAll.invoke(null, refs, patterns);
		}catch(InvocationTargetException e){throw new IllegalStateException("Match consumer threw in seedless-repeat", e.getCause());}
		check(labels.length==1, "One literal query must produce one labeler result.");
		final Object exactValue=count.get(exact), labelValue=count.get(labels[0]);
		check(exactValue instanceof Number && labelValue instanceof Number, "Both hit counts must remain numeric.");
		final long exactCount=((Number)exactValue).longValue(), labelCount=((Number)labelValue).longValue();
		final int exactSequence=sequence.getInt(exact), exactStart=pos.getInt(exact);
		check(refs.size()==2 && refs.get(0)==forward && refs.get(1)==reverse,
				"Match consumers must preserve borrowed reference entries.");
		check(patterns.size()==1 && patterns.get(0)==query && Arrays.equals(query, new byte[]{'A', 'A', 'A'}),
				"Match consumers must preserve the borrowed literal query.");
		// AAA occurs at starts 0, 1 and 2 of AAAAA, and nowhere in TTTTT.
		final boolean passed=exactCount==3 && labelCount==2 && exactSequence==0 && exactStart==2;
		System.out.println("case\texpected_exact\tobserved_exact\texpected_labeler\tobserved_labeler\texpected_sequence\tobserved_sequence\texpected_start\tobserved_start\tpassed");
		System.out.println("seedless-repeat\t3\t"+exactCount+"\t2\t"+labelCount+"\t0\t"+exactSequence+"\t2\t"+exactStart+"\t"+passed);
		System.out.println("FUSION_MATCH_COUNT_TEST case=consumer assertions="+
				FusionMatchCountTest.class.desiredAssertionStatus()+" cases=1 failed="+(passed ? 0 : 1));
		System.out.flush();
		if(!passed){throw new AssertionError("Fusion match-count contract failed: 1 of 1");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal orientation pairs; no computed reverse-complement oracle is shared with production. */
	private static ArrayList<String> pair(final String forward, final String reverse){
		return new ArrayList<String>(Arrays.asList(forward, reverse));
	}

	/** Setup/ownership failures stay loud with assertions disabled. */
	private static void check(final boolean pass, final String why){
		if(!pass){throw new IllegalStateException("Match-count fixture setup: "+why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Actual exhaustive lookup; invoked once per case without replacing its implementation. */
	private final Method find;
	/** Read-only accessors for numeric hit count and last int coordinates; never inject result state. */
	private final Field count, sequence, pos;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed small pilot population; repeated immutable references avoid bulk input allocation. */
	private static final int PILOT_PAIRS=1024;
	/** Exact 2^31-hit population; its invocation is separate from the cost pilot. */
	private static final int BOUNDARY_PAIRS=32768;
}
