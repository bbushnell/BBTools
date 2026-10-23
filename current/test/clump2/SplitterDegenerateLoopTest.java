package test.clump2;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Set;

import clump2.Clump;
import stream.Read;

/**
 * Finite oracle for a degenerate correlated split in Splitter's work-list loop.
 * @author Sucrose
 */
public final class SplitterDegenerateLoopTest{

	/** Runs baseline or fixed expectations for the pending-sibling loop case. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final GuardedLoopList queue=new GuardedLoopList(
				makeDegenerateCorrelatedClump(), makeSplitSibling());
		if(fixed){
			checkFixedVisit(queue);
		}else{
			checkBaselineSkip(queue);
		}
		System.out.println("case\tpendingVisited\tdegenerateSelected\tout2\tout4\tout24");
		System.out.println("degenerate-pending\t"+queue.pendingVisited+'\t'+
				queue.degenerateSelected+'\t'+queue.out2+'\t'+
				queue.out4+'\t'+queue.out24);
	}

	/** Checks the pre-fix path that selects the appended copy before the sibling. */
	private static void checkBaselineSkip(final GuardedLoopList queue) throws Exception{
		try{
			splitOnPivot(queue);
			throw new AssertionError("Expected the appended degenerate copy to be selected");
		}catch(DegenerateSelectedException expected){
			check(!queue.pendingVisited, "Pending sibling should not have been visited");
			check(queue.size()==3, "Expected one degenerate append, got "+queue.size());
			check(queue.rawGet(0)==null, "Expected the first slot to be cleared");
			check(queue.rawGet(1)==queue.pendingSibling, "Expected pending sibling at index 1");
			check(queue.rawGet(2).size()==DEGENERATE_READS,
					"Expected every degenerate read in the appended copy");
			queue.degenerateSelected=true;
		}
	}

	/** Checks the repaired path visits and splits the pending sibling. */
	private static void checkFixedVisit(final GuardedLoopList queue) throws Exception{
		final ArrayList<Clump> out=splitOnPivot(queue);
		for(Clump clump : out){
			if(clump.size()==2){queue.out2++;}
			else if(clump.size()==4){queue.out4++;}
			else if(clump.size()==DEGENERATE_READS){queue.out24++;}
			else{throw new AssertionError("Unexpected clump size "+clump.size());}
		}
		checkReturnedReadsOnce(queue, out);
		check(queue.pendingVisited, "Expected the pending sibling to be visited");
		check(queue.out2==1, "Expected one 2-read sibling split, got "+queue.out2);
		check(queue.out4==1, "Expected one 4-read sibling split, got "+queue.out4);
		check(queue.out24==1, "Expected one degenerate 24-read original, got "+queue.out24);
	}

	/** Checks that fixed-mode output returns every input read identity once. */
	private static void checkReturnedReadsOnce(final GuardedLoopList queue,
			final ArrayList<Clump> out){
		final Set<Read> expected=Collections.newSetFromMap(
				new IdentityHashMap<Read, Boolean>());
		final Set<Read> observed=Collections.newSetFromMap(
				new IdentityHashMap<Read, Boolean>());
		expected.addAll(queue.expectedReads);

		for(Clump clump : out){
			for(Read read : clump){
				check(expected.contains(read), "Returned an unexpected read: "+read);
				check(observed.add(read), "Returned a read identity twice: "+read);
			}
		}
		check(observed.size()==expected.size(), "Expected "+expected.size()+
				" returned reads, got "+observed.size());
	}

	/** Creates a clump whose correlated pivots cover all reads. */
	private static Clump makeDegenerateCorrelatedClump() throws Exception{
		final Clump clump=Clump.makeClump(KMER);
		int id=0;
		for(int i=0; i<20; i++){addRead(clump, id++, -1, VAR_A, VAR_B);}
		for(int i=0; i<2; i++){addRead(clump, id++, VAR_B>>2, VAR_A);}
		for(int i=0; i<2; i++){addRead(clump, id++, VAR_A>>2, VAR_B);}
		return clump;
	}

	/** Creates a sibling with a literal 4/2 single-pivot split. */
	private static Clump makeSplitSibling() throws Exception{
		final Clump clump=Clump.makeClump(KMER);
		int id=100;
		for(int i=0; i<4; i++){addRead(clump, id++, -1);}
		for(int i=0; i<2; i++){addRead(clump, id++, VAR_A>>2, VAR_A);}
		return clump;
	}

	/** Adds a low-quality read with optional high-quality consensus support. */
	private static void addRead(final Clump clump, final int id,
			final int highQualityRefPosition, final int... vars) throws Exception{
		final byte[] bases=basesFor(vars).getBytes(StandardCharsets.US_ASCII);
		final byte[] qualities=new byte[READ_LENGTH];
		Arrays.fill(qualities, (byte)2);
		if(highQualityRefPosition>=0){qualities[highQualityRefPosition]=(byte)40;}

		final Read read=new Read(bases, qualities, "splitter-degenerate-"+id, id);
		makeKey(read);
		clump.add(read);
	}

	/** Creates and attaches the package-private ReadKey used by Splitter. */
	private static void makeKey(final Read read) throws Exception{
		final Class<?> readKey=Class.forName("clump2.ReadKey");
		final Method makeKey=readKey.getDeclaredMethod(
				"makeKey", Read.class, boolean.class);
		final Method set=readKey.getDeclaredMethod(
				"set", long.class, int.class, boolean.class);
		makeKey.setAccessible(true);
		set.setAccessible(true);
		set.invoke(makeKey.invoke(null, read, true), KMER, READ_POS, false);
	}

	/** Calls package-private Splitter.splitOnPivot. */
	@SuppressWarnings("unchecked")
	private static ArrayList<Clump> splitOnPivot(final ArrayList<Clump> list)
			throws Exception{
		final Class<?> splitter=Class.forName("clump2.Splitter");
		final Method splitOnPivot=splitter.getDeclaredMethod("splitOnPivot",
				ArrayList.class);
		splitOnPivot.setAccessible(true);
		try{
			return (ArrayList<Clump>)splitOnPivot.invoke(null, list);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof DegenerateSelectedException){
				throw (DegenerateSelectedException)cause;
			}
			throw e;
		}
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Creates a read that contains the supplied variants at their clump columns. */
	private static String basesFor(final int... vars){
		final char[] bases=new char[READ_LENGTH];
		Arrays.fill(bases, 'T');
		for(int var : vars){
			bases[var>>2]=BASES[var&0x3];
		}
		return new String(bases);
	}

	/** ArrayList guard that makes a repeating degenerate append finite. */
	private static final class GuardedLoopList extends ArrayList<Clump>{

		GuardedLoopList(final Clump degenerate, final Clump pendingSibling_){
			add(degenerate);
			add(pendingSibling=pendingSibling_);
			expectedReads=new ArrayList<Read>(degenerate.size()+pendingSibling.size());
			expectedReads.addAll(degenerate);
			expectedReads.addAll(pendingSibling);
		}

		@Override
		public Clump get(int index){
			check(++selections<=MAX_SELECTIONS,
					"Exceeded the guarded loop selection cap: "+selections);
			if(index==2 && !pendingVisited){
				throw new DegenerateSelectedException();
			}
			return super.get(index);
		}

		@Override
		public Clump set(int index, Clump element){
			if(index==1 && element==null && super.get(index)==pendingSibling){
				pendingVisited=true;
			}
			return super.set(index, element);
		}

		Clump rawGet(final int index){return super.get(index);}

		final Clump pendingSibling;
		final ArrayList<Read> expectedReads;
		boolean pendingVisited;
		boolean degenerateSelected;
		int out2;
		int out4;
		int out24;
		private int selections;

		private static final int MAX_SELECTIONS=16;
		private static final long serialVersionUID=1L;
	}

	/** Marker for selecting a copy that would otherwise repeat indefinitely. */
	private static final class DegenerateSelectedException extends RuntimeException{
		private static final long serialVersionUID=1L;
	}

	private SplitterDegenerateLoopTest(){}

	private static final long KMER=3L;
	private static final int VAR_A=(10<<2)|0;
	private static final int VAR_B=(20<<2)|1;
	private static final int DEGENERATE_READS=24;
	private static final int READ_LENGTH=40;
	private static final int READ_POS=0;
	private static final char[] BASES={'A', 'C', 'G', 'T'};
}
