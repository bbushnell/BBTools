package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;

import assemble.Contig;
import assemble.ContigGraphClassifier;

/** Checks finite depth classes and rejection of deliberately supplied NaN coverage.
 * This direct-API fixture does not demonstrate an assembly input producing NaN.
 * @author Brian Bushnell, Shinobu */
public final class ContigGraphCoverageTest{

	/** Runs isolated single-contig cases; explicit checks remain active without assertions. */
	public static void main(final String[] args){
		int failures=0;
		failures+=check("zero", 0, Contig.DEPTH_LOW, false);
		failures+=check("medium", 40, Contig.DEPTH_MEDIUM, false);
		failures+=check("high", 251, Contig.DEPTH_HIGH, false);
		failures+=check("nan", Float.NaN, -1, true);
		if(failures>0){throw new AssertionError("Contig coverage validation failures="+failures);}
		System.out.println("CONTIG_GRAPH_COVERAGE_TEST_PASS cases=4");
	}

	/** Compares a fresh graph against a literal bucket or a specific invalid-coverage error. */
	private static int check(final String name, final float coverage, final int expectedDepth,
			final boolean expectedRejected){
		assert(expectedRejected || expectedDepth>=0) : "Finite controls require an explicit expected depth bucket.";
		final byte[] bases=new byte[20];
		Arrays.fill(bases, (byte)'A');
		final Contig contig=new Contig(bases);
		contig.id=0;
		contig.coverage=coverage;
		final ArrayList<Contig> graph=new ArrayList<Contig>(1);
		graph.add(contig);
		boolean rejected=false, diagnostic=false;
		try{
			new ContigGraphClassifier(graph, 50, 3, 4, 40, 250).classify();
		}catch(IllegalArgumentException e){
			rejected=true;
			diagnostic="NaN coverage for contig 0.".equals(e.getMessage());
		}
		final boolean correct=expectedRejected ? rejected && diagnostic :
				!rejected && contig.graphDepthClass==expectedDepth;
		System.out.println("case="+name+" expectedRejected="+expectedRejected+" rejected="+rejected+
				" expectedDepth="+expectedDepth+" actualDepth="+contig.graphDepthClass+
				" diagnostic="+diagnostic+" correct="+correct);
		return correct ? 0 : 1;
	}
}
