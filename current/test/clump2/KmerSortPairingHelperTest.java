package test.clump2;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import clump2.KmerSort;
import stream.Read;

/** Finite no-I/O oracle for KmerSort pairing helpers.
 * @author Jean
 */
public final class KmerSortPairingHelperTest{

	/** Records all eight helper cases, then checks the requested baseline/fixed phase. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected phase: baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown phase: "+args[0]);

		System.out.println("case\tsize\tpaired\tretained\tbaseline_ok\tfixed_ok");
		boolean baselineOk=true, fixedOk=true;
		final Result[] results=new Result[] {
				checkEmpty("name-empty", NAME_SORT),
				checkPair("name-pair", NAME_SORT, "pair/1", "pair/2", 1),
				checkTrailing("name-trailing", NAME_SORT, true),
				checkOdd("name-odd", NAME_SORT, true),
				checkEmpty("id-empty", ID_SORT),
				checkPair("id-pair", ID_SORT, "pair1", "pair2", 10),
				checkTrailing("id-trailing", ID_SORT, false),
				checkOdd("id-odd", ID_SORT, false)
		};
		for(Result result : results){
			System.out.println(result);
			baselineOk&=result.baselineOk;
			fixedOk&=result.fixedOk;
		}
		check(fixed ? fixedOk : baselineOk,
				"Unexpected "+args[0]+" KmerSort pairing-helper behavior");
	}

	/** Ensures an empty helper input remains empty in both phases. */
	private static Result checkEmpty(String label, Method method) throws Exception{
		return evaluate(label, invoke(method), new Read[0], new Read[0], null, null);
	}

	/** Verifies that a matched pair is retained and linked in both phases. */
	private static Result checkPair(String label, Method method, String id1, String id2,
			long numericID) throws Exception{
		final Read pair1=read(id1, numericID, 0);
		final Read pair2=read(id2, numericID, 1);
		final Read[] expected=new Read[] {pair1};
		return evaluate(label, invoke(method, pair1, pair2), expected, expected, pair1, pair2);
	}

	/** Verifies that baseline drops, and fixed code retains, one trailing singleton. */
	private static Result checkTrailing(String label, Method method, boolean names)
			throws Exception{
		final Read pair1=read(names ? "alpha/1" : "alpha1", 2, 0);
		final Read pair2=read(names ? "alpha/2" : "alpha2", 2, 1);
		final Read singleton=read(names ? "omega/1" : "omega1", 3, 0);
		return evaluate(label, invoke(method, pair1, pair2, singleton),
				new Read[] {pair1}, new Read[] {pair1, singleton}, pair1, pair2);
	}

	/** Verifies that baseline drops, and fixed code retains, the last odd singleton. */
	private static Result checkOdd(String label, Method method, boolean names)
			throws Exception{
		final Read one=read(names ? "apple" : "one", 4, 0);
		final Read two=read(names ? "banana" : "two", 5, 0);
		final Read three=read(names ? "carrot" : "three", 6, 0);
		return evaluate(label, invoke(method, one, two, three),
				new Read[] {one, two}, new Read[] {one, two, three}, null, null);
	}

	/** Computes exact baseline and fixed expectations without throwing immediately. */
	private static Result evaluate(String label, ArrayList<Read> out,
			Read[] baselineExpected, Read[] fixedExpected, Read first, Read mate){
		final boolean paired=linkedPair(out, first, mate);
		final int retained=retainedSingletons(out);
		final boolean baselineOk=matches(out, baselineExpected, first, mate);
		final boolean fixedOk=matches(out, fixedExpected, first, mate);
		return new Result(label, out.size(), paired, retained, baselineOk, fixedOk);
	}

	/** Checks exact top-level read identity and mate state. */
	private static boolean matches(ArrayList<Read> out, Read[] expected,
			Read first, Read mate){
		if(out.size()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			final Read read=out.get(i);
			if(read!=expected[i]){return false;}
			if(read==first && mate!=null){
				if(!linkedPair(out, first, mate)){return false;}
			}else if(read.mate!=null){return false;}
		}
		return true;
	}

	/** True only when the first top-level read owns its literal expected mate. */
	private static boolean linkedPair(ArrayList<Read> out, Read first, Read mate){
		if(out.isEmpty() || first==null || mate==null || out.get(0)!=first){return false;}
		return first.mate==mate && mate.mate==first
				&& first.pairnum()==0 && mate.pairnum()==1;
	}

	/** Counts top-level singleton reads retained by the helper. */
	private static int retainedSingletons(ArrayList<Read> out){
		int singletons=0;
		for(Read read : out){
			if(read.mate==null){singletons++;}
		}
		return singletons;
	}

	/** Invokes one static KmerSort helper. */
	@SuppressWarnings("unchecked")
	private static ArrayList<Read> invoke(Method method, Read... reads) throws Exception{
		final ArrayList<Read> list=new ArrayList<Read>(reads.length);
		for(Read read : reads){list.add(read);}
		return (ArrayList<Read>)method.invoke(null, list, Boolean.TRUE);
	}

	/** Makes a tiny read with deterministic pair state. */
	private static Read read(String id, long numericID, int pairnum){
		final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null, id, numericID);
		read.setPairnum(pairnum);
		return read;
	}

	/** Lightweight assertion independent of JVM assertion flags. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Finds an otherwise package-private static helper. */
	private static Method method(String name){
		try{
			final Method method=KmerSort.class.getDeclaredMethod(name, ArrayList.class,
					boolean.class);
			method.setAccessible(true);
			return method;
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Recorded result for one helper case. */
	private static final class Result{

		Result(String label_, int size_, boolean paired_, int retained_,
				boolean baselineOk_, boolean fixedOk_){
			label=label_;
			size=size_;
			paired=paired_;
			retained=retained_;
			baselineOk=baselineOk_;
			fixedOk=fixedOk_;
		}

		@Override
		public String toString(){
			return label+'\t'+size+'\t'+paired+'\t'+retained+'\t'+baselineOk+'\t'+fixedOk;
		}

		final String label;
		final int size;
		final boolean paired;
		final int retained;
		final boolean baselineOk;
		final boolean fixedOk;
	}

	private static final Method NAME_SORT=method("nameSort");
	private static final Method ID_SORT=method("idSort");

}
