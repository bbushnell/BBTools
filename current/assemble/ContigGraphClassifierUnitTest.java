package assemble;

import java.util.ArrayList;
import java.util.Arrays;

/** Deterministic unit tests for graph-derived contig QC classification. */
public class ContigGraphClassifierUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("twoSidedConnector", ContigGraphClassifierUnitTest::twoSidedConnector);
		failures+=run("oneSidedHair", ContigGraphClassifierUnitTest::oneSidedHair);
		failures+=run("hopDistance", ContigGraphClassifierUnitTest::hopDistance);
		failures+=run("sameSideOrientations", ContigGraphClassifierUnitTest::sameSideOrientations);
		failures+=run("sameAnchorOnBothEnds", ContigGraphClassifierUnitTest::sameAnchorOnBothEnds);
		failures+=run("highDepthShortIsNotNormal", ContigGraphClassifierUnitTest::highDepthShortIsNotNormal);
		failures+=run("topologyIgnoresDepth", ContigGraphClassifierUnitTest::topologyIgnoresDepth);
		failures+=run("branchedAnchorClasses", ContigGraphClassifierUnitTest::branchedAnchorClasses);
		failures+=run("selfLoop", ContigGraphClassifierUnitTest::selfLoop);
		failures+=run("nonreciprocalEdge", ContigGraphClassifierUnitTest::nonreciprocalEdge);
		failures+=run("mismatchedReciprocalPath", ContigGraphClassifierUnitTest::mismatchedReciprocalPath);
		failures+=run("redundantShortPayloadIgnored", ContigGraphClassifierUnitTest::redundantShortPayloadIgnored);
		failures+=run("depthBoundaries", ContigGraphClassifierUnitTest::depthBoundaries);
		failures+=run("threeAnchorCap", ContigGraphClassifierUnitTest::threeAnchorCap);
		failures+=run("cyclicMultiAnchor", ContigGraphClassifierUnitTest::cyclicMultiAnchor);
		failures+=run("namedHeaderIncludesClass", ContigGraphClassifierUnitTest::namedHeaderIncludesClass);
		failures+=run("cleanupEligibility", ContigGraphClassifierUnitTest::cleanupEligibility);
		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	private static int run(final String name, final Test test){
		try{
			test.run();
			System.out.println("PASS: "+name);
			return 0;
		}catch(Throwable t){
			System.err.println("FAIL: "+name);
			t.printStackTrace();
			return 1;
		}
	}

	private static void twoSidedConnector(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(1), false);
		connect(list.get(1), true, list.get(2), false);
		classify(list);
		check(list.get(0).graphLengthClass==Contig.LENGTH_LONG, "Left anchor was not long");
		check(list.get(1).graphClass==Contig.GRAPH_CONNECTED && list.get(1).graphClassHop==1,
				"Direct connector was not connected-1");
		check(list.get(2).graphLengthClass==Contig.LENGTH_LONG, "Right anchor was not long");
	}

	private static void oneSidedHair(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2));
		connect(list.get(0), true, list.get(1), false);
		classify(list);
		check(list.get(1).graphClass==Contig.GRAPH_TERMINAL && list.get(1).graphClassHop==1,
				"One-sided hair was not terminal-1");
	}

	private static void hopDistance(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(20, 2),
				contig(20, 2), contig(100, 8));
		for(int i=0; i<list.size()-1; i++){connect(list.get(i), true, list.get(i+1), false);}
		classify(list);
		check(list.get(1).graphClass==Contig.GRAPH_CONNECTED && list.get(1).graphClassHop==3,
				"First connector had wrong hop");
		check(list.get(2).graphClass==Contig.GRAPH_CONNECTED && list.get(2).graphClassHop==2,
				"Middle connector had wrong hop");
		check(list.get(3).graphClass==Contig.GRAPH_CONNECTED && list.get(3).graphClassHop==3,
				"Last connector had wrong hop");
	}

	private static void sameSideOrientations(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(1), true);
		connect(list.get(1), false, list.get(2), false);
		classify(list);
		check(list.get(1).graphClass==Contig.GRAPH_CONNECTED && list.get(1).graphClassHop==1,
				"Orientation-0/3 connector was not connected-1");
	}

	private static void sameAnchorOnBothEnds(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2));
		connect(list.get(0), false, list.get(1), false);
		connect(list.get(0), true, list.get(1), true);
		classify(list);
		check(list.get(1).graphClass==Contig.GRAPH_LOOPBACK && list.get(1).graphClassHop==1,
				"One normal anchor on both sides was not classified as loopback-1");
	}

	private static void highDepthShortIsNotNormal(){
		final ArrayList<Contig> list=list(contig(20, 100));
		classify(list);
		check(list.get(0).graphLengthClass==Contig.LENGTH_SHORT &&
				list.get(0).graphDepthClass==Contig.DEPTH_MEDIUM &&
				list.get(0).graphClass==Contig.GRAPH_UNANCHORED,
				"High-depth short contig did not retain independent short/medium/unanchored labels");
	}

	private static void topologyIgnoresDepth(){
		final ArrayList<Contig> list=list(contig(100, 20), contig(20, 2), contig(100, 20));
		connect(list.get(0), true, list.get(1), false);
		connect(list.get(1), true, list.get(2), false);
		classify(list);
		check(list.get(1).graphClass==Contig.GRAPH_CONNECTED,
				"Topology classification improperly depended on depth");
	}

	private static void branchedAnchorClasses(){
		ArrayList<Contig> list=list(contig(100, 8), contig(100, 8), contig(20, 2));
		connect(list.get(0), true, list.get(2), false);
		connect(list.get(1), true, list.get(2), false);
		classify(list);
		check(list.get(2).graphClass==Contig.GRAPH_BRANCHED_TERMINAL,
				"A 0/2+ topology was not branched-terminal");

		list=list(contig(100, 8), contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(2), false);
		connect(list.get(1), true, list.get(2), false);
		connect(list.get(2), true, list.get(3), false);
		classify(list);
		check(list.get(2).graphClass==Contig.GRAPH_BRANCHED_CONNECTED,
				"A 1/2+ topology was not branched-connected");
	}

	private static void selfLoop(){
		final ArrayList<Contig> list=list(contig(20, 2));
		list.get(0).leftCode=Tadpole.LOOP;
		classify(list);
		check(list.get(0).graphClass==Contig.GRAPH_SELF_LOOP, "Loop code was not classified as self-loop");
	}

	private static void nonreciprocalEdge(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(1), false);
		final Contig mid=list.get(1), right=list.get(2);
		mid.addRightEdge(new Edge(mid.id, right.id, 1, 1, 2, new byte[]{0}));
		classify(list);
		check(mid.graphClass==Contig.GRAPH_TERMINAL,
				"Nonreciprocal second edge did not leave the contig terminal");
	}

	private static void mismatchedReciprocalPath(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(1), false);
		final Contig mid=list.get(1), right=list.get(2);
		mid.addRightEdge(new Edge(mid.id, right.id, 4, 1, 2, new byte[]{'C','A','A','A'}));
		right.addLeftEdge(new Edge(right.id, mid.id, 4, 2, 2, new byte[]{'A','A','A','A'}));
		classify(list);
		check(mid.graphClass==Contig.GRAPH_TERMINAL,
				"Sequence-incompatible reciprocal edges created a false second anchor");
	}

	private static void redundantShortPayloadIgnored(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(20, 2), contig(100, 8));
		connect(list.get(0), true, list.get(1), false);
		final Contig mid=list.get(1), right=list.get(2);
		mid.addRightEdge(new Edge(mid.id, right.id, 1, 1, 2, new byte[]{'C'}));
		right.addLeftEdge(new Edge(right.id, mid.id, 1, 2, 2, new byte[]{'G'}));
		classify(list);
		check(mid.graphClass==Contig.GRAPH_CONNECTED,
				"Redundant represented payload bytes invalidated an exact contig join");
	}

	private static void depthBoundaries(){
		final ArrayList<Contig> list=list(contig(20, 4), contig(20, 40), contig(20, 250), contig(20, 251));
		classify(list);
		check(list.get(0).graphDepthClass==Contig.DEPTH_LOW, "Low boundary was not inclusive");
		check(list.get(1).graphDepthClass==Contig.DEPTH_MEDIUM, "Medium boundary was not inclusive");
		check(list.get(2).graphDepthClass==Contig.DEPTH_MEDIUM, "High cutoff was not included in medium");
		check(list.get(3).graphDepthClass==Contig.DEPTH_HIGH, "High class did not start above its cutoff");
	}

	private static void threeAnchorCap(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(100, 8), contig(100, 8), contig(20, 2));
		for(int i=0; i<3; i++){connect(list.get(i), true, list.get(3), false);}
		classify(list);
		check(list.get(3).graphClass==Contig.GRAPH_BRANCHED_TERMINAL,
				"Three anchors did not retain capped 2+ semantics");
	}

	private static void cyclicMultiAnchor(){
		final ArrayList<Contig> list=list(contig(100, 8), contig(100, 8), contig(100, 8),
				contig(20, 2), contig(20, 2));
		for(int i=0; i<3; i++){connect(list.get(i), true, list.get(3), false);}
		connect(list.get(3), true, list.get(4), false);
		connect(list.get(4), true, list.get(3), false);
		classify(list);
		check(list.get(3).graphClass>=0 && list.get(4).graphClass>=0,
				"A cyclic multi-anchor component was not classified");
	}

	private static void namedHeaderIncludesClass(){
		final ArrayList<Contig> list=list(contig(20, 2));
		list.get(0).name="named";
		classify(list);
		check(list.get(0).name().equals("named,class=short/low/unanchored"),
				"Named contig omitted its graph classification: "+list.get(0).name());
	}

	private static void cleanupEligibility(){
		final Contig c=contig(20, 2);
		c.graphLengthClass=Contig.LENGTH_SHORT;
		c.graphDepthClass=Contig.DEPTH_LOW;
		c.graphClass=Contig.GRAPH_TERMINAL;
		check(Tadpole.graphEvictionEligible(c, 50, 1<<Contig.DEPTH_LOW), "Eligible short/low contig was retained");
		c.graphDepthClass=Contig.DEPTH_MEDIUM;
		check(!Tadpole.graphEvictionEligible(c, 50, 1<<Contig.DEPTH_LOW), "Depth mask was ignored");
		check(Tadpole.graphEvictionEligible(c, 50, 1<<Contig.DEPTH_MEDIUM), "Selected medium depth was not eligible");
		c.graphLengthClass=Contig.LENGTH_LONG;
		check(!Tadpole.graphEvictionEligible(c, 50, 15), "Long contig was eligible for graph cleanup");
		c.graphLengthClass=Contig.LENGTH_SHORT;
		c.graphClass=Contig.GRAPH_SELF_LOOP;
		check(Tadpole.graphEvictionEligible(c, 50, 15), "Explicitly selectable self-loop was ineligible for graph cleanup");
		c.graphClass=Contig.GRAPH_UNANCHORED;
		c.graphDepthClass=Contig.DEPTH_MEDIUM;
		check(Tadpole.defaultGraphSweepEligible(c, 50, Contig.DEPTH_MEDIUM),
				"Default sweep did not include medium-depth unanchored contigs");
		c.graphDepthClass=Contig.DEPTH_HIGH;
		check(!Tadpole.defaultGraphSweepEligible(c, 50, Contig.DEPTH_MEDIUM),
				"Default sweep included high-depth unanchored contigs");
		c.graphClass=Contig.GRAPH_TERMINAL;
		c.graphDepthClass=Contig.DEPTH_SEMILOW;
		check(Tadpole.defaultGraphSweepEligible(c, 50, Contig.DEPTH_SEMILOW),
				"Default sweep did not include semilow terminal contigs");
		c.graphDepthClass=Contig.DEPTH_MEDIUM;
		check(!Tadpole.defaultGraphSweepEligible(c, 50, Contig.DEPTH_SEMILOW),
				"Default sweep included medium-depth terminal contigs");
		check(Tadpole.parseGraphTopologyMask("terminal")==
				((1<<Contig.GRAPH_TERMINAL)|(1<<Contig.GRAPH_BRANCHED_TERMINAL)),
				"Logical terminal class did not include its branched subtype");
		check(Tadpole.parseGraphTopologyMask("connected")==
				((1<<Contig.GRAPH_CONNECTED)|(1<<Contig.GRAPH_BRANCHED_CONNECTED)|
						(1<<Contig.GRAPH_MULTI_CONNECTED)),
				"Logical connected class did not include its branched subtypes");
		check(Tadpole.parseGraphTopologyMask("self-loop")==1<<Contig.GRAPH_SELF_LOOP,
				"Self-loop graph class was not selectable");
	}

	private static void classify(final ArrayList<Contig> list){
		new ContigGraphClassifier(list, 50, 3, 4, 40, 250).classify();
	}

	private static Contig contig(final int length, final float coverage){
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'N');
		final Contig c=new Contig(bases);
		c.coverage=coverage;
		return c;
	}

	private static ArrayList<Contig> list(final Contig... array){
		final ArrayList<Contig> list=new ArrayList<Contig>(array.length);
		for(int i=0; i<array.length; i++){array[i].id=i; list.add(array[i]);}
		return list;
	}

	private static void connect(final Contig a, final boolean aRight,
			final Contig b, final boolean bRight){
		final int forward=(aRight ? 1 : 0)|(bRight ? 2 : 0);
		final int reverse=(bRight ? 1 : 0)|(aRight ? 2 : 0);
		final Edge ab=new Edge(a.id, b.id, 1, forward, 2, new byte[]{'N'});
		final Edge ba=new Edge(b.id, a.id, 1, reverse, 2, new byte[]{'N'});
		if(aRight){a.addRightEdge(ab);}else{a.addLeftEdge(ab);}
		if(bRight){b.addRightEdge(ba);}else{b.addLeftEdge(ba);}
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private interface Test {void run();}
}
