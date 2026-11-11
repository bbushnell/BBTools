package test.assemble;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;

import assemble.Edge;

/** Finite, boundary and overflowing depth cases for dormant Edge.merge.
 * Invokes the package-private merge helper by reflection so the baseline can
 * exercise current Edge bytes without widening the production API.
 * @author Jean
 */
public final class EdgeMergeDepthTest{

	/** Runs every row before reporting failures so both overflow branches are visible. */
	public static void main(final String[] args){
		if(args.length>0){throw new IllegalArgumentException("No arguments are supported.");}
		final ArrayList<String> failures=new ArrayList<String>();
		final StringBuilder report=new StringBuilder();
		for(final Case c : CASES){checkCase(c, report, failures);}
		System.out.print(report.toString());
		if(!failures.isEmpty()){
			final StringBuilder sb=new StringBuilder("Edge.merge depth failures=").append(failures.size());
			for(final String failure : failures){sb.append('\n').append(failure);}
			throw new AssertionError(sb.toString());
		}
		System.out.println("EDGE_MERGE_DEPTH_TEST_OK checks="+CASES.length);
	}

	/** Merges one same-topology pair and verifies depth, payload selection and path-stat retention. */
	private static void checkCase(final Case c, final StringBuilder report, final ArrayList<String> failures){
		final Edge oldEdge=edge(c.oldDepth, c.oldLength, c.oldBases, c.oldOverlap, c.oldSourceTrim,
				c.oldDestTrim, c.oldPathMin, c.oldPathMax, c.oldPathMean);
		final Edge incoming=edge(c.incomingDepth, c.incomingLength, c.incomingBases, c.incomingOverlap,
				c.incomingSourceTrim, c.incomingDestTrim, c.incomingPathMin, c.incomingPathMax,
				c.incomingPathMean);
		final int expectedLength=c.incomingWins ? c.incomingLength : c.oldLength;
		final byte[] expectedBases=c.incomingWins ? c.incomingBases : c.oldBases;
		final int expectedOverlap=c.incomingWins ? c.incomingOverlap : c.oldOverlap;
		final int expectedSourceTrim=c.incomingWins ? c.incomingSourceTrim : c.oldSourceTrim;
		final int expectedDestTrim=c.incomingWins ? c.incomingDestTrim : c.oldDestTrim;
		Throwable thrown=null;
		try{
			MERGE.invoke(oldEdge, incoming);
		}catch(InvocationTargetException e){
			thrown=e.getCause();
		}catch(IllegalAccessException e){
			thrown=e;
		}catch(IllegalArgumentException e){
			thrown=e;
		}
		final boolean depthOk=thrown==null && oldEdge.depth==c.expectedDepth;
		final boolean payloadOk=thrown==null && oldEdge.length==expectedLength && oldEdge.bases==expectedBases &&
				oldEdge.overlap==expectedOverlap && oldEdge.sourceTrim==expectedSourceTrim &&
				oldEdge.destTrim==expectedDestTrim;
		final boolean topologyOk=thrown==null && oldEdge.origin==ORIGIN && oldEdge.destination==DESTINATION &&
				oldEdge.orientation==ORIENTATION;
		final boolean pathStatsOk=thrown==null && oldEdge.pathMinDepth==c.oldPathMin &&
				oldEdge.pathMaxDepth==c.oldPathMax &&
				Double.compare(oldEdge.pathMeanDepth, c.oldPathMean)==0;
		appendReport(report, c, oldEdge, depthOk, payloadOk, topologyOk, pathStatsOk, thrown);
		check(thrown==null, failures, c.name+" threw "+thrown);
		check(depthOk, failures, c.name+" depth mismatch: expected="+c.expectedDepth+
				", actual="+oldEdge.depth+", mathematical="+c.expectedLong);
		check(payloadOk, failures, c.name+" payload mismatch: expected incoming="+c.incomingWins+
				", length="+oldEdge.length+", overlap="+oldEdge.overlap+
				", sourceTrim="+oldEdge.sourceTrim+", destTrim="+oldEdge.destTrim);
		check(topologyOk, failures, c.name+" topology changed: origin="+oldEdge.origin+
				", destination="+oldEdge.destination+", orientation="+oldEdge.orientation);
		check(pathStatsOk, failures, c.name+" path statistics changed: min="+oldEdge.pathMinDepth+
				", max="+oldEdge.pathMaxDepth+", mean="+oldEdge.pathMeanDepth);
	}

	/** Creates a same-topology edge with caller-selected mutable fields. */
	private static Edge edge(final int depth, final int length, final byte[] bases, final int overlap,
			final int sourceTrim, final int destTrim, final int pathMin, final int pathMax,
			final double pathMean){
		final Edge edge=new Edge(ORIGIN, DESTINATION, length, ORIENTATION, depth, bases, overlap,
				sourceTrim, destTrim);
		edge.pathMinDepth=pathMin;
		edge.pathMaxDepth=pathMax;
		edge.pathMeanDepth=pathMean;
		return edge;
	}

	/** Appends one observed Edge.merge row for baseline comparison. */
	private static void appendReport(final StringBuilder report, final Case c, final Edge actual,
			final boolean depthOk, final boolean payloadOk, final boolean topologyOk,
			final boolean pathStatsOk, final Throwable thrown){
		report.append("case=").append(c.name);
		report.append(" old=").append(c.oldDepth);
		report.append(" incoming=").append(c.incomingDepth);
		report.append(" expectedLong=").append(c.expectedLong);
		report.append(" expectedDepth=").append(c.expectedDepth);
		report.append(" actualInt=").append(actual.depth);
		report.append(" incomingWins=").append(c.incomingWins);
		report.append(" depthOk=").append(depthOk);
		report.append(" payloadOk=").append(payloadOk);
		report.append(" topologyOk=").append(topologyOk);
		report.append(" pathStatsOk=").append(pathStatsOk);
		if(thrown!=null){report.append(" thrown=").append(thrown.getClass().getName());}
		report.append('\n');
	}

	/** Explicit check helper that is active even when assertions are disabled. */
	private static void check(final boolean condition, final ArrayList<String> failures, final String message){
		if(!condition){failures.add(message);}
	}

	/** Finds the package-private Edge.merge helper without changing Edge's API. */
	private static Method mergeMethod(){
		try{
			final Method method=Edge.class.getDeclaredMethod("merge", Edge.class);
			method.setAccessible(true);
			return method;
		}catch(NoSuchMethodException e){
			throw new ExceptionInInitializerError(e);
		}catch(SecurityException e){
			throw new ExceptionInInitializerError(e);
		}
	}

	/** Literal test-case depths and distinct old/incoming payload fields. */
	private static final class Case{
		Case(final String name_, final int oldDepth_, final int incomingDepth_,
				final long expectedLong_, final int expectedDepth_, final boolean incomingWins_){
			name=name_;
			oldDepth=oldDepth_;
			incomingDepth=incomingDepth_;
			expectedLong=expectedLong_;
			expectedDepth=expectedDepth_;
			incomingWins=incomingWins_;
			oldLength=OLD_BASES.length;
			incomingLength=INCOMING_BASES.length;
			oldBases=OLD_BASES;
			incomingBases=INCOMING_BASES;
			oldOverlap=OLD_OVERLAP;
			incomingOverlap=INCOMING_OVERLAP;
			oldSourceTrim=OLD_SOURCE_TRIM;
			incomingSourceTrim=INCOMING_SOURCE_TRIM;
			oldDestTrim=OLD_DEST_TRIM;
			incomingDestTrim=INCOMING_DEST_TRIM;
			oldPathMin=OLD_PATH_MIN;
			incomingPathMin=INCOMING_PATH_MIN;
			oldPathMax=OLD_PATH_MAX;
			incomingPathMax=INCOMING_PATH_MAX;
			oldPathMean=OLD_PATH_MEAN;
			incomingPathMean=INCOMING_PATH_MEAN;
		}

		final String name;
		final int oldDepth;
		final int incomingDepth;
		final long expectedLong;
		final int expectedDepth;
		final boolean incomingWins;
		final int oldLength;
		final int incomingLength;
		final byte[] oldBases;
		final byte[] incomingBases;
		final int oldOverlap;
		final int incomingOverlap;
		final int oldSourceTrim;
		final int incomingSourceTrim;
		final int oldDestTrim;
		final int incomingDestTrim;
		final int oldPathMin;
		final int incomingPathMin;
		final int oldPathMax;
		final int incomingPathMax;
		final double oldPathMean;
		final double incomingPathMean;
	}

	private static final int ORIGIN=17;
	private static final int DESTINATION=29;
	private static final int ORIENTATION=3;
	private static final Method MERGE=mergeMethod();
	private static final byte[] OLD_BASES=new byte[]{'A', 'C', 'G', 'T'};
	private static final byte[] INCOMING_BASES=new byte[]{'T', 'G', 'C', 'A', 'A'};
	private static final int OLD_OVERLAP=10;
	private static final int INCOMING_OVERLAP=30;
	private static final int OLD_SOURCE_TRIM=50;
	private static final int INCOMING_SOURCE_TRIM=70;
	private static final int OLD_DEST_TRIM=90;
	private static final int INCOMING_DEST_TRIM=110;
	private static final int OLD_PATH_MIN=130;
	private static final int INCOMING_PATH_MIN=150;
	private static final int OLD_PATH_MAX=170;
	private static final int INCOMING_PATH_MAX=190;
	private static final double OLD_PATH_MEAN=150.25;
	private static final double INCOMING_PATH_MEAN=170.25;
	private static final Case[] CASES={
			new Case("old-greater-finite", 5, 3, 8L, 8, false),
			new Case("equal-finite", 9, 9, 18L, 18, false),
			new Case("incoming-greater-finite", 3, 5, 8L, 8, true),
			new Case("exact-boundary", Integer.MAX_VALUE-7, 7, 2147483647L, Integer.MAX_VALUE, false),
			new Case("old-greater-overflow", Integer.MAX_VALUE-7, 11, 2147483651L,
					Integer.MAX_VALUE, false),
			new Case("incoming-greater-overflow", 11, Integer.MAX_VALUE-7, 2147483651L,
					Integer.MAX_VALUE, true)
	};

}
