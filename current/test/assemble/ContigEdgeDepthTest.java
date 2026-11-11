package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;

import assemble.Contig;
import assemble.Edge;

/** Finite duplicate-edge depth cases for Contig's retained left/right edge merge paths.
 * Same destination/orientation edges with different payloads bypass the duplicate-payload
 * assertion and merge their supporting depths. The proposed overflow contract
 * follows HashArray1D/HashArrayU1D positive count saturation.
 * @author Jean
 */
public final class ContigEdgeDepthTest{

	/** Runs left and right merge cases that exercise ordinary, boundary and overflowing sums. */
	public static void main(final String[] args){if(args.length>0){runTiming(args); return;}
		final Case[] cases={
				new Case(3, 5, 8L, 8),
				new Case(1, Integer.MAX_VALUE-1, Integer.MAX_VALUE, Integer.MAX_VALUE),
				new Case(1500000000, 1500000000, 3000000000L, Integer.MAX_VALUE),
				new Case(1, Integer.MAX_VALUE, 2147483648L, Integer.MAX_VALUE)
		};
		int checks=0;
		final ArrayList<String> failures=new ArrayList<String>();
		final StringBuilder report=new StringBuilder();
		for(final Case c : cases){
			checkLeft(c, report, failures);
			checkRight(c, report, failures);
			checks+=2;
		}
		System.out.print(report.toString());
		if(!failures.isEmpty()){
			final StringBuilder sb=new StringBuilder("Contig edge-depth failures=").append(failures.size());
			for(final String failure : failures){sb.append('\n').append(failure);}
			throw new AssertionError(sb.toString());
		}
		System.out.println("CONTIG_EDGE_DEPTH_TEST_OK checks="+checks);
	}

	/** Runs the optional non-overflow retained-edge timing mode. */
	private static void runTiming(final String[] args){
		final TimingConfig config=TimingConfig.parse(args);
		final java.lang.management.ThreadMXBean bean=cpuClock();
		final TimingState state=new TimingState(config);
		final long warmupChecksum=state.run(config.warmup);
		state.checkEndpointState(config.warmup);
		final long expectedWarmup=state.expectedChecksum(config.warmup);
		if(warmupChecksum!=expectedWarmup){
			throw new AssertionError("Warmup checksum mismatch: observed="+warmupChecksum+
					", expected="+expectedWarmup);
		}
		final long startWall=System.nanoTime(), startCpu=currentCpuTime(bean);
		final long checksum=state.run(config.iterations);
		final long elapsedWall=System.nanoTime()-startWall, elapsedCpu=currentCpuTime(bean)-startCpu;
		state.checkEndpointState(config.iterations);
		final long expected=state.expectedChecksum(config.iterations);
		if(checksum!=expected){
			throw new AssertionError("Checksum mismatch: observed="+checksum+", expected="+expected);
		}
		System.out.println("CONTIG_EDGE_DEPTH_TIMING_OK iterations="+config.iterations+
				" warmup="+config.warmup+" checksum="+checksum+" expected="+expected+
				" warmupChecksum="+warmupChecksum+" expectedWarmup="+expectedWarmup+
				" wallNs="+elapsedWall+" mainThreadCpuNs="+elapsedCpu);
	}

	/** Adds two merge-compatible left edges and checks depth, count and payload retention. */
	private static void checkLeft(final Case c, final StringBuilder report, final ArrayList<String> failures){
		final Contig contig=new Contig(new byte[]{'A'}, 0);
		final Edge oldEdge=edge(c.oldDepth, 0, payload('A'));
		final Edge incoming=edge(c.incomingDepth, 0, payload('C'));
		contig.addLeftEdge(oldEdge);
		contig.addLeftEdge(incoming);
		final int count=contig.leftEdgeCount();
		final Edge retained=contig.getLeftEdge(1, 0);
		final int actual=(retained==null ? Integer.MIN_VALUE : retained.depth);
		final boolean countOk=count==1, inPlaceOk=retained==oldEdge;
		final boolean payloadOk=retained!=null && retained.bases==incoming.bases && Arrays.equals(retained.bases, payload('C'));
		appendReport(report, "left", c, actual, count, countOk, inPlaceOk, payloadOk);
		check(countOk, failures, "Left duplicate should merge into one retained edge; count="+count);
		check(inPlaceOk, failures, "Left merge should update the retained edge in place.");
		check(payloadOk, failures,
				"Left merge should retain the incoming payload when incoming depth is at least old depth.");
		if(retained!=null){checkDepth("left", c, actual, failures);}
	}

	/** Adds two merge-compatible right edges and checks depth, count and payload retention. */
	private static void checkRight(final Case c, final StringBuilder report, final ArrayList<String> failures){
		final Contig contig=new Contig(new byte[]{'A'}, 0);
		final Edge oldEdge=edge(c.oldDepth, 1, payload('A'));
		final Edge incoming=edge(c.incomingDepth, 1, payload('C'));
		contig.addRightEdge(oldEdge);
		contig.addRightEdge(incoming);
		final int count=contig.rightEdgeCount();
		final Edge retained=contig.getRightEdge(1, 1);
		final int actual=(retained==null ? Integer.MIN_VALUE : retained.depth);
		final boolean countOk=count==1, inPlaceOk=retained==oldEdge;
		final boolean payloadOk=retained!=null && retained.bases==incoming.bases && Arrays.equals(retained.bases, payload('C'));
		appendReport(report, "right", c, actual, count, countOk, inPlaceOk, payloadOk);
		check(countOk, failures, "Right duplicate should merge into one retained edge; count="+count);
		check(inPlaceOk, failures, "Right merge should update the retained edge in place.");
		check(payloadOk, failures,
				"Right merge should retain the incoming payload when incoming depth is at least old depth.");
		if(retained!=null){checkDepth("right", c, actual, failures);}
	}

	/** Creates a merge-compatible edge with caller-selected depth and payload. */
	private static Edge edge(final int depth, final int orientation, final byte[] bases){
		return new Edge(0, 1, bases.length, orientation, depth, bases);
	}

	/** Creates a same-length payload with caller-selected contents. */
	private static byte[] payload(final char c){
		final byte[] bases=new byte[7];
		Arrays.fill(bases, (byte)c);
		return bases;
	}

	/** Returns a main-thread CPU clock or throws before the timing mode can claim success. */
	private static java.lang.management.ThreadMXBean cpuClock(){
		final java.lang.management.ThreadMXBean bean=java.lang.management.ManagementFactory.getThreadMXBean();
		if(!bean.isCurrentThreadCpuTimeSupported()){
			throw new UnsupportedOperationException("Current-thread CPU time is not supported.");
		}
		if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
		final long time=bean.getCurrentThreadCpuTime();
		if(time<0){throw new IllegalStateException("Current-thread CPU time is not enabled.");}
		return bean;
	}

	/** Returns nonnegative main-thread CPU nanoseconds. */
	private static long currentCpuTime(final java.lang.management.ThreadMXBean bean){
		final long time=bean.getCurrentThreadCpuTime();
		if(time<0){throw new IllegalStateException("Current-thread CPU time became unavailable.");}
		return time;
	}

	/** Appends one observed retained-edge merge row for baseline comparison. */
	private static void appendReport(final StringBuilder report, final String side, final Case c,
			final int actual, final int count, final boolean countOk, final boolean inPlaceOk, final boolean payloadOk){
		report.append("side=").append(side).append(" old=").append(c.oldDepth);
		report.append(" incoming=").append(c.incomingDepth);
		report.append(" expectedLong=").append(c.expectedLong);
		report.append(" expectedDepth=").append(c.expectedDepth);
		report.append(" actualInt=").append(actual);
		report.append(" retainedCount=").append(count);
		report.append(" countOk=").append(countOk);
		report.append(" inPlaceOk=").append(inPlaceOk);
		report.append(" payloadOk=").append(payloadOk).append('\n');
	}

	/** Checks that a retained edge reports the saturated expected sum, not a wrapped int. */
	private static void checkDepth(final String side, final Case c, final int actual, final ArrayList<String> failures){
		check(actual==c.expectedDepth, failures, side+" depth sum mismatch: old="+c.oldDepth+
				", incoming="+c.incomingDepth+", expectedLong="+c.expectedLong+
				", expectedDepth="+c.expectedDepth+", actualInt="+actual);
	}

	/** Explicit check helper that is active even when assertions are disabled. */
	private static void check(final boolean condition, final ArrayList<String> failures, final String message){
		if(!condition){failures.add(message);}
	}

	/** Command-line settings for the optional timing mode. */
	private static final class TimingConfig{
		/** Parses {@code timing=t iterations=N warmup=N old=A,B incoming=C,D} in any order. */
		static TimingConfig parse(final String[] args){
			assert(args.length>0) : "TimingConfig is only used for non-empty arguments.";
			boolean timing=false;
			int iterations=-1, warmup=-1;
			int[] oldDepths=null, incomingDepths=null;
			for(int i=0; i<args.length; i++){
				final String arg=args[i];
				final int equals=arg.indexOf('=');
				final String name=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase(java.util.Locale.ROOT);
				final String value=(equals<0 ? null : arg.substring(equals+1));
				if("timing".equals(name)){
					timing=(value==null || value.length()<1 || parse.Parse.parseBoolean(value));
				}else if("iterations".equals(name)){
					iterations=parse.Parse.parseIntKMG(value);
				}else if("warmup".equals(name)){
					warmup=parse.Parse.parseIntKMG(value);
				}else if("old".equals(name)){
					oldDepths=parseDepths(value, arg);
				}else if("incoming".equals(name)){
					incomingDepths=parseDepths(value, arg);
				}else{
					throw new IllegalArgumentException("Unknown timing argument: "+arg);
				}
			}
			if(!timing){throw new IllegalArgumentException("timing flag is required.");}
			if(iterations<1){throw new IllegalArgumentException("iterations must be positive: "+iterations);}
			if(warmup<1){throw new IllegalArgumentException("warmup must be positive: "+warmup);}
			if(oldDepths==null){throw new IllegalArgumentException("old= depth array is required.");}
			if(incomingDepths==null){throw new IllegalArgumentException("incoming= depth array is required.");}
			if(oldDepths.length!=incomingDepths.length){
				throw new IllegalArgumentException("old and incoming depth arrays differ: "+
						oldDepths.length+" != "+incomingDepths.length);
			}
			long cycleChecksum=0;
			for(int i=0; i<oldDepths.length; i++){
				final int oldDepth=oldDepths[i], incomingDepth=incomingDepths[i];
				final long sum=(long)oldDepth+incomingDepth;
				if(oldDepth<1 || incomingDepth<oldDepth || sum>=Integer.MAX_VALUE){
					throw new IllegalArgumentException("Depths must be positive, merge-eligible and finite: old="+
							oldDepth+", incoming="+incomingDepth+", sum="+sum);
				}
				cycleChecksum+=2L*sum;
			}
			return new TimingConfig(iterations, warmup, oldDepths, incomingDepths, cycleChecksum);
		}

		/** Parses a non-empty comma-delimited int list. */
		private static int[] parseDepths(final String text, final String arg){
			if(text==null || text.length()<1){throw new IllegalArgumentException("Empty depth array: "+arg);}
			int count=1;
			for(int i=0; i<text.length(); i++){
				if(text.charAt(i)==','){count++;}
			}
			final int[] depths=new int[count];
			for(int i=0, cell=0; cell<count; cell++){
				int comma=text.indexOf(',', i);
				if(comma<0){comma=text.length();}
				if(comma==i){throw new IllegalArgumentException("Empty depth in "+arg);}
				depths[cell]=parse.Parse.parseIntKMG(text.substring(i, comma));
				i=comma+1;
			}
			return depths;
		}

		private TimingConfig(final int iterations_, final int warmup_, final int[] oldDepths_,
				final int[] incomingDepths_, final long cycleChecksum_){
			iterations=iterations_;
			warmup=warmup_;
			oldDepths=oldDepths_;
			incomingDepths=incomingDepths_;
			cycleChecksum=cycleChecksum_;
		}

		final int iterations;
		final int warmup;
		final int[] oldDepths;
		final int[] incomingDepths;
		final long cycleChecksum;
	}

	/** Preallocated state for retained-edge merge timing loops with no fixture-added per-iteration allocation. */
	private static final class TimingState{
		TimingState(final TimingConfig config){
			oldDepths=config.oldDepths;
			incomingDepths=config.incomingDepths;
			cycleChecksum=config.cycleChecksum;
			left.leftEdges=new ArrayList<Edge>(1);
			right.rightEdges=new ArrayList<Edge>(1);
			left.addLeftEdge(leftOld);
			right.addRightEdge(rightOld);
		}

		/** Runs alternating left and right ordinary-depth merges and returns observed-depth checksum. */
		long run(final int iterations){
			long checksum=0;
			for(int i=0, depthIndex=0; i<iterations; i++){
				final int oldDepth=oldDepths[depthIndex], incomingDepth=incomingDepths[depthIndex];
				checksum+=runLeft(oldDepth, incomingDepth);
				checksum+=runRight(oldDepth, incomingDepth);
				if(++depthIndex==oldDepths.length){depthIndex=0;}
			}
			return checksum;
		}

		/** Computes the expected checksum without an O(iterations) second loop. */
		long expectedChecksum(final int iterations){
			final int cycles=iterations/oldDepths.length, remainder=iterations-(cycles*oldDepths.length);
			long checksum=cycleChecksum*cycles;
			for(int i=0; i<remainder; i++){
				checksum+=2L*(oldDepths[i]+incomingDepths[i]);
			}
			return checksum;
		}

		/** Resets reusable objects, runs one left-edge merge and returns the observed depth. */
		private int runLeft(final int oldDepth, final int incomingDepth){
			resetRetained(leftOld, oldDepth, LEFT_OLD);
			leftIncoming.depth=incomingDepth;
			left.addLeftEdge(leftIncoming);
			return leftOld.depth;
		}

		/** Resets reusable objects, runs one right-edge merge and returns the observed depth. */
		private int runRight(final int oldDepth, final int incomingDepth){
			resetRetained(rightOld, oldDepth, RIGHT_OLD);
			rightIncoming.depth=incomingDepth;
			right.addRightEdge(rightIncoming);
			return rightOld.depth;
		}

		/** Verifies the last retained edges after the timed loop has finished. */
		void checkEndpointState(final int iterations){
			final int depthIndex=(iterations-1)%oldDepths.length;
			final int expectedDepth=oldDepths[depthIndex]+incomingDepths[depthIndex];
			final Edge leftRetained=left.getLeftEdge(1, 0), rightRetained=right.getRightEdge(1, 1);
			if(left.leftEdgeCount()!=1 || leftRetained!=leftOld || leftOld.bases!=LEFT_INCOMING ||
					leftOld.depth!=expectedDepth){
				throw new AssertionError("Bad left timing endpoint: count="+left.leftEdgeCount()+
						", sameEdge="+(leftRetained==leftOld)+", payload="+
						(leftOld.bases==LEFT_INCOMING)+", depth="+leftOld.depth+
						", expected="+expectedDepth);
			}
			if(right.rightEdgeCount()!=1 || rightRetained!=rightOld || rightOld.bases!=RIGHT_INCOMING ||
					rightOld.depth!=expectedDepth){
				throw new AssertionError("Bad right timing endpoint: count="+right.rightEdgeCount()+
						", sameEdge="+(rightRetained==rightOld)+", payload="+
						(rightOld.bases==RIGHT_INCOMING)+", depth="+rightOld.depth+
						", expected="+expectedDepth);
			}
		}

		/** Restores retained fields that the retained-edge merge path mutates. */
		private static void resetRetained(final Edge edge, final int depth, final byte[] bases){
			edge.bases=bases;
			edge.depth=depth;
		}

		private final Contig left=new Contig(new byte[]{'A'}, 0);
		private final Contig right=new Contig(new byte[]{'A'}, 0);
		private final Edge leftOld=edge(1, 0, LEFT_OLD);
		private final Edge leftIncoming=edge(1, 0, LEFT_INCOMING);
		private final Edge rightOld=edge(1, 1, RIGHT_OLD);
		private final Edge rightIncoming=edge(1, 1, RIGHT_INCOMING);
		private final int[] oldDepths;
		private final int[] incomingDepths;
		private final long cycleChecksum;
	}

	/** Literal test-case depths with the mathematical sum and the expected stored int. */
	private static final class Case{
		Case(final int oldDepth_, final int incomingDepth_, final long expectedLong_, final int expectedDepth_){
			oldDepth=oldDepth_;
			incomingDepth=incomingDepth_;
			expectedLong=expectedLong_;
			expectedDepth=expectedDepth_;
			assert(incomingDepth>=oldDepth) :
					"Case construction must keep the incoming edge merge-eligible: old="+oldDepth+", incoming="+incomingDepth;
		}

		final int oldDepth;
		final int incomingDepth;
		final long expectedLong;
		final int expectedDepth;
	}

	private static final byte[] LEFT_OLD=payload('A');
	private static final byte[] LEFT_INCOMING=payload('C');
	private static final byte[] RIGHT_OLD=payload('G');
	private static final byte[] RIGHT_INCOMING=payload('T');

}
