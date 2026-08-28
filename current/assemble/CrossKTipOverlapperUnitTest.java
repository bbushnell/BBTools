package assemble;

import java.util.ArrayList;
import java.util.HashMap;

/** Assertion-based tests for exact low-depth tip-overlap joining. */
public class CrossKTipOverlapperUnitTest {

	public static void main(String[] args){
		BubblePopper.verbose=false;
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=false;
		BubblePopper.crossKMerge=true;
		BubblePopper.crossKMaxDepthRatio=3;
		BubblePopper.validateGraph=true;

		int failures=0;
		failures+=run("uniqueReciprocalOverlap", CrossKTipOverlapperUnitTest::uniqueReciprocalOverlap);
		failures+=run("branchTipIgnored", CrossKTipOverlapperUnitTest::branchTipIgnored);
		failures+=run("equalBestIsAmbiguous", CrossKTipOverlapperUnitTest::equalBestIsAmbiguous);
		failures+=run("reverseOrientedOverlap", CrossKTipOverlapperUnitTest::reverseOrientedOverlap);
		failures+=run("cyclicComponentRejected", CrossKTipOverlapperUnitTest::cyclicComponentRejected);
		failures+=run("graphKUnbranchedOverlap", CrossKTipOverlapperUnitTest::graphKUnbranchedOverlap);
		failures+=run("graphKBranchIgnored", CrossKTipOverlapperUnitTest::graphKBranchIgnored);
		BubblePopper.crossKMerge=false;
		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	private static int run(String name, Test test){
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

	private static void uniqueReciprocalOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		ArrayList<Contig> contigs=list(a, b);
		int pairs=new CrossKTipOverlapper(contigs, 5, 9).addEdges();
		check(pairs==1, "Expected one overlap pair, got "+pairs);
		check(a.rightEdgeCount()==1 && b.leftEdgeCount()==1, "Missing reciprocal overlap edges");
		check(a.rightEdges.get(0).overlap==7, "Wrong overlap length: "+a.rightEdges.get(0).overlap);
		int merged=popper(contigs).expand(a);
		check(merged==1, "Expected one merge, got "+merged);
		check(new String(a.bases).equals("AAAACCCCGGGGTTTT"), "Wrong merged sequence: "+new String(a.bases));
	}

	private static void branchTipIgnored(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9).addEdges()==0, "Ineligible branch tip was joined");
	}

	private static void equalBestIsAmbiguous(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		Contig c=contig(2, "CCCGGGGAAAA", true, false);
		ArrayList<Contig> contigs=list(a, b, c);
		check(new CrossKTipOverlapper(contigs, 5, 9).addEdges()==0, "Equal best overlaps were not rejected");
	}

	private static void reverseOrientedOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "AAAACCCCGGG", false, true);
		ArrayList<Contig> contigs=list(a, b);
		int pairs=new CrossKTipOverlapper(contigs, 5, 9).addEdges();
		check(pairs==1, "Expected one reverse-oriented pair, got "+pairs);
		check(a.rightEdges.get(0).destRight(), "Destination orientation was not right-facing");
		int merged=popper(contigs).expand(a);
		check(merged==1, "Expected one reverse-oriented merge, got "+merged);
		check(new String(a.bases).equals("AAAACCCCGGGGTTTT"), "Wrong reverse-oriented product: "+new String(a.bases));
	}

	private static void cyclicComponentRejected(){
		Contig a=contig(0, "GGATTAAC", true, true);
		Contig b=contig(1, "AACGGCCT", true, true);
		Contig c=contig(2, "CCTCCGGA", true, true);
		ArrayList<Contig> contigs=list(a, b, c);
		check(new CrossKTipOverlapper(contigs, 3, 3).addEdges()==0, "Cyclic overlap component was not rejected");
		for(Contig x : contigs){check(x.leftEdgeCount()==0 && x.rightEdgeCount()==0, "Cycle left graph edges");}
	}

	private static void graphKUnbranchedOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", false, false);
		a.rightCode=Tadpole.KEEP_GOING;
		b.leftCode=Tadpole.KEEP_GOING;
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9, true).addEdges()==1,
				"Unbranched graph-k overlap was not selected");
	}

	private static void graphKBranchIgnored(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", false, false);
		a.rightCode=Tadpole.F_BRANCH;
		b.leftCode=Tadpole.KEEP_GOING;
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9, true).addEdges()==0,
				"Branched graph-k overlap was selected");
	}

	private static BubblePopper popper(ArrayList<Contig> contigs){
		HashMap<Integer, ArrayList<Edge>> map=new HashMap<Integer, ArrayList<Edge>>();
		for(Contig c : contigs){
			if(c.leftEdges!=null){for(Edge e : c.leftEdges){add(map, e);}}
			if(c.rightEdges!=null){for(Edge e : c.rightEdges){add(map, e);}}
		}
		return new BubblePopper(contigs, map, 5);
	}

	private static void add(HashMap<Integer, ArrayList<Edge>> map, Edge e){
		ArrayList<Edge> list=map.get(e.destination);
		if(list==null){list=new ArrayList<Edge>(); map.put(e.destination, list);}
		list.add(e);
	}

	private static Contig contig(int id, String bases, boolean left, boolean right){
		Contig c=new Contig(bases.getBytes(), id);
		c.coverage=20;
		c.leftCode=c.rightCode=Tadpole.DEAD_END;
		c.leftBridgeEndpoint=left;
		c.rightBridgeEndpoint=right;
		return c;
	}

	private static ArrayList<Contig> list(Contig... contigs){
		ArrayList<Contig> list=new ArrayList<Contig>();
		for(Contig c : contigs){list.add(c);}
		return list;
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private interface Test {void run();}
}
