package test.pacbio;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import pacbio.RemoveAdapters3;
import stream.Read;

/** Captures the existing paired-input split assertion without starting worker or writer threads.
 * This is a baseline guard observation, not a paired-CLI termination test.
 * @author Brian Bushnell, Sucrose
 */
public final class RemoveAdapters3PairedSplitTest{

	/** Invokes only the private list-splitting helper and checks its input ownership.
	 * @param args No arguments
	 * @throws Exception If the private helper seam changes or throws an unexpected failure
	 */
	public static void main(String[] args) throws Exception{
		require(args.length==0, "This bounded baseline probe takes no arguments");
		require(RemoveAdapters3.class.desiredAssertionStatus(), "The production single-end assertion must remain enabled");
		final Class<?> worker=Class.forName("pacbio.RemoveAdapters3$ProcessThread");
		final Method split=worker.getDeclaredMethod("split", ArrayList.class);
		split.setAccessible(true);
		final Read left=new Read("AAAA".getBytes(StandardCharsets.US_ASCII), null, "pair/1", 0);
		final Read right=new Read("TTTT".getBytes(StandardCharsets.US_ASCII), null, "pair/2", 0);
		final ArrayList<Read> input=new ArrayList<Read>();
		input.add(left);
		final Object result=split.invoke(null, input);
		require(result instanceof ArrayList<?>, "Unpaired split must return a list");
		final ArrayList<?> unpaired=(ArrayList<?>)result;
		require(unpaired.size()==1 && unpaired.get(0)==left,
				"An unmarked unpaired read must retain its object identity");
		left.mate=right;
		right.mate=left;
		left.setPairnum(0);
		right.setPairnum(1);
		boolean asserted=false;
		try{split.invoke(null, input);}
		catch(InvocationTargetException failure){
			if(!(failure.getCause() instanceof AssertionError)){throw failure;}
			asserted=true;
		}
		require(asserted, "Baseline RemoveAdapters3 split must expose its paired-input assertion");
		require(input.size()==1 && input.get(0)==left && left.mate==right && right.mate==left,
				"The rejected baseline list and reciprocal mate links must remain intact");
		System.out.println("case\texpected\tobserved\tresult");
		System.out.println("unpaired-unmarked\toriginal_read\toriginal_read\tPASS");
		System.out.println("paired-split-guard\tAssertionError\tAssertionError\tCAPTURED");
	}

	/** Keeps the oracle visible even if the fixture itself has assertions disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private RemoveAdapters3PairedSplitTest(){}
}
