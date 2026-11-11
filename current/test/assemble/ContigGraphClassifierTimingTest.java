package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import assemble.Contig;
import assemble.ContigGraphClassifier;
import assemble.Edge;
import parse.Parse;
import parse.Parser;

/**
 * Times repeated direct graph classification over a deterministic mixed-depth graph.
 * The fixture measures public ContigGraphClassifier classify() passes and their
 * scratch allocation; it does not build a graph from reads or prove whole-assembler speed.
 *
 * @author Brian Bushnell, Shinobu, Jean
 */
public final class ContigGraphClassifierTimingTest{

	/** Builds one graph, warms it up, times fresh classifier passes and verifies literal classes. */
	public static void main(final String[] args){
		final Config config=Config.parse(args);
		final long constructStart=System.nanoTime();
		final ArrayList<Contig> graph=buildGraph(config.chains);
		final long constructNanos=System.nanoTime()-constructStart;

		final long warmupStart=System.nanoTime();
		runClassify(graph, config.warmup);
		final long warmupNanos=System.nanoTime()-warmupStart;
		final long warmupChecksum=(config.warmup>0 ? verifyGraph(graph, config.chains) : 0);

		final long classifyStart=System.nanoTime();
		runClassify(graph, config.repeats);
		final long classifyNanos=System.nanoTime()-classifyStart;
		final long checksum=verifyGraph(graph, config.chains);

		System.out.println("CONTIG_GRAPH_CLASSIFIER_TIMING_PASS chains="+config.chains+
				" contigs="+(config.chains*4)+" warmup="+config.warmup+
				" repeats="+config.repeats+" constructionNanos="+constructNanos+
				" warmupNanos="+warmupNanos+" classifyNanos="+classifyNanos+
				" warmupChecksum="+warmupChecksum+" checksum="+checksum);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs count classify passes, constructing a fresh classifier for each pass. */
	private static void runClassify(final ArrayList<Contig> graph, final int count){
		for(int i=0; i<count; i++){
			latestSummary=new ContigGraphClassifier(
					graph, MIN_LONG, K, LOW_CUTOFF, MEDIUM_CUTOFF, HIGH_CUTOFF).classify();
		}
	}

	/** Builds disjoint long--short--short--long chains with reciprocal overlap edges. */
	private static ArrayList<Contig> buildGraph(final int chains){
		final ArrayList<Contig> graph=new ArrayList<Contig>(chains*4);
		for(int chain=0; chain<chains; chain++){
			final int base=graph.size();
			graph.add(contig(base, LONG_LENGTH, LOW_COVERAGE));
			graph.add(contig(base+1, SHORT_LENGTH, SEMILOW_COVERAGE));
			graph.add(contig(base+2, SHORT_LENGTH, MEDIUM_COVERAGE));
			graph.add(contig(base+3, LONG_LENGTH, HIGH_COVERAGE));
			connect(graph.get(base), graph.get(base+1));
			connect(graph.get(base+1), graph.get(base+2));
			connect(graph.get(base+2), graph.get(base+3));
		}
		return graph;
	}

	/** Returns an all-A contig with consecutive id and fixed coverage. */
	private static Contig contig(final int id, final int length, final float coverage){
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'A');
		final Contig c=new Contig(bases, id);
		c.coverage=coverage;
		return c;
	}

	/** Connects a.right to b.left and b.left back to a.right through an exact overlap edge. */
	private static void connect(final Contig a, final Contig b){
		a.addRightEdge(new Edge(a.id, b.id, 0, 1, 1, null, OVERLAP));
		b.addLeftEdge(new Edge(b.id, a.id, 0, 2, 1, null, OVERLAP));
	}

	/** Verifies every contig's literal expected classes and returns a graph-field checksum. */
	private static long verifyGraph(final ArrayList<Contig> graph, final int chains){
		checkEquals("Graph size", graph.size(), chains*4);
		final int[] lengthCounts=new int[2], depthCounts=new int[4], topologyCounts=new int[8];
		long checksum=0xcbf29ce484222325L;
		for(int i=0; i<graph.size(); i++){
			final Contig c=graph.get(i);
			final int within=i&3;
			final int expectedLength=(within==0 || within==3 ? Contig.LENGTH_LONG : Contig.LENGTH_SHORT);
			final int expectedDepth=within;
			final int expectedTopology=(within==0 || within==3 ? Contig.GRAPH_TERMINAL : Contig.GRAPH_CONNECTED);
			final int expectedHop=(within==0 || within==3 ? 3 : 2);
			final int expectedLeftEdges=(within==0 ? 0 : 1);
			final int expectedRightEdges=(within==3 ? 0 : 1);
			final int leftEdges=c.leftEdgeCount();
			final int rightEdges=c.rightEdgeCount();
			checkContig(c, i, expectedLength, expectedDepth, expectedTopology, expectedHop,
					expectedLeftEdges, expectedRightEdges, leftEdges, rightEdges);
			lengthCounts[c.graphLengthClass]++;
			depthCounts[c.graphDepthClass]++;
			topologyCounts[c.graphClass]++;
			checksum=mix(checksum, c.id);
			checksum=mix(checksum, c.graphLengthClass);
			checksum=mix(checksum, c.graphDepthClass);
			checksum=mix(checksum, c.graphClass);
			checksum=mix(checksum, c.graphClassHop);
			checksum=mix(checksum, leftEdges);
			checksum=mix(checksum, rightEdges);
		}
		checkEquals("Long count", lengthCounts[Contig.LENGTH_LONG], chains*2);
		checkEquals("Short count", lengthCounts[Contig.LENGTH_SHORT], chains*2);
		for(int i=0; i<depthCounts.length; i++){
			checkBucketEquals("Depth bucket", i, depthCounts[i], chains);
		}
		checkEquals("Terminal count", topologyCounts[Contig.GRAPH_TERMINAL], chains*2);
		checkEquals("Connected count", topologyCounts[Contig.GRAPH_CONNECTED], chains*2);
		for(int i=0; i<topologyCounts.length; i++){
			if(i!=Contig.GRAPH_TERMINAL && i!=Contig.GRAPH_CONNECTED){
				checkBucketEquals("Unexpected topology count", i, topologyCounts[i], 0);
			}
		}
		return checksum;
	}

	/** Mixes one observed public graph field into the retained checksum. */
	private static long mix(final long checksum, final int value){
		return (checksum^value)*0x100000001b3L;
	}

	/** Checks one contig's literal expected fields with no success-path diagnostics. */
	private static void checkContig(final Contig c, final int index, final int expectedLength,
			final int expectedDepth, final int expectedTopology, final int expectedHop,
			final int expectedLeftEdges, final int expectedRightEdges, final int leftEdges,
			final int rightEdges){
		if(c.id!=index || c.graphLengthClass!=expectedLength || c.graphDepthClass!=expectedDepth ||
				c.graphClass!=expectedTopology || c.graphClassHop!=expectedHop ||
				leftEdges!=expectedLeftEdges || rightEdges!=expectedRightEdges){
			throw new AssertionError("Contig "+index+" expected id="+index+
					" length="+expectedLength+" depth="+expectedDepth+
					" topology="+expectedTopology+" hop="+expectedHop+
					" leftEdges="+expectedLeftEdges+" rightEdges="+expectedRightEdges+
					"; got id="+c.id+" length="+c.graphLengthClass+
					" depth="+c.graphDepthClass+" topology="+c.graphClass+
					" hop="+c.graphClassHop+" leftEdges="+leftEdges+
					" rightEdges="+rightEdges);
		}
	}

	/** Checks one named count with no success-path diagnostics. */
	private static void checkEquals(final String name, final int actual, final int expected){
		if(actual!=expected){throw new AssertionError(name+" "+actual+" != "+expected);}
	}

	/** Checks one bucket count with no success-path diagnostics. */
	private static void checkBucketEquals(final String name, final int bucket, final int actual,
			final int expected){
		if(actual!=expected){
			throw new AssertionError(name+" "+bucket+" count "+actual+" != "+expected);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parsed timing fixture dimensions. */
	private static final class Config{
		/** Parses flag=value dimensions. */
		static Config parse(final String[] args){
			int chains=8, repeats=3, warmup=1;
			for(String arg : args){
				if(Parser.isJavaFlag(arg)){continue;}
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String name=pair[0].toLowerCase(Locale.ROOT);
				final String value=pair.length>1 ? pair[1] : null;
				if(name.equals("chains") || name.equals("c")){
					chains=Parse.parseIntKMG(value);
				}else if(name.equals("repeats") || name.equals("repeat") || name.equals("r")){
					repeats=Parse.parseIntKMG(value);
				}else if(name.equals("warmup") || name.equals("w")){
					warmup=Parse.parseIntKMG(value);
				}else{throw new IllegalArgumentException("Unknown flag: "+arg);}
			}
			return new Config(chains, repeats, warmup);
		}

		/** Stores positive chain/repeat counts and a nonnegative warmup count. */
		Config(final int chains_, final int repeats_, final int warmup_){
			chains=chains_;
			repeats=repeats_;
			warmup=warmup_;
			if(chains<1 || repeats<1 || warmup<0){
				throw new IllegalArgumentException("Expected chains>=1, repeats>=1, warmup>=0; got chains="+
						chains+" repeats="+repeats+" warmup="+warmup+".");
			}
			Math.multiplyExact(chains, 4);
		}

		/** Independent chain count, timed classify repeats and untimed warmup passes. */
		final int chains, repeats, warmup;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Length cutoff and graph k, in bases. */
	private static final int MIN_LONG=64, K=31;
	/** Contig lengths and exact overlap, in bases. */
	private static final int LONG_LENGTH=96, SHORT_LENGTH=48, OVERLAP=16;
	/** One representative coverage value per depth bucket. */
	private static final float LOW_COVERAGE=0, SEMILOW_COVERAGE=10, MEDIUM_COVERAGE=100, HIGH_COVERAGE=1000;
	/** Fixed low, medium and high depth boundaries. */
	private static final float LOW_CUTOFF=4, MEDIUM_CUTOFF=40, HIGH_CUTOFF=250;
	/** Keeps each returned summary observable; overwritten once per whole-graph pass. */
	private static volatile Object latestSummary;
}
