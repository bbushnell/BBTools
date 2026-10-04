package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.ByteFile;
import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;
import ukmer.Kmer;
import shared.Timer;
import stream.Read;
import stream.ReadInputStream;

/** Deterministic topology, representation, and DOT tests for non-mutating graph exploration.
 * @author Fischl */
public final class TadpoleGraphTest {

	/** Runs in-memory tests and optionally writes small CLI fixtures to an existing directory. */
	public static void main(final String[] args){
		if(args.length==2 && args[0].equals("dense")){
			densePaths(args[1]);
			return;
		}
		if(args.length>=4 && args[0].equals("paths")){
			verifyPaths(args[1], args[2], Arrays.copyOfRange(args, 3, args.length));
			return;
		}
		if(args.length==3 && args[0].equals("gfa")){
			verifyAdjacentLinks(args[1], Integer.parseInt(args[2]));
			return;
		}
		for(int k : new int[]{11, 31, 62}){orientations(k); branches(k); limits(k);}
		final boolean packed=Kmer.PACKED;
		Kmer.PACKED=true;
		try{orientations(33); branches(33); orientations(127); branches(127);}
		finally{Kmer.PACKED=packed;}
		sharedAnchor();
		cycle();
		missingAndShort();
		dot();
		pathSpelling();
		if(args.length>0){fixtures(args[0]);}
		System.out.println("PASS: TadpoleGraph topology, depths, limits, thread parity, names, and DOT styling.");
	}

	/** Exercises all four endpoint orientations and exact representative depth arithmetic. */
	private static void orientations(final int k){
		final int sourceLength=Math.max(100, k+20), targetStart=sourceLength+80;
		final String genome=random(targetStart+sourceLength+40, 42);
		final String left=genome.substring(0, sourceLength), right=genome.substring(targetStart);
		for(int a=0; a<2; a++){
			for(int b=0; b<2; b++){
				final ArrayList<Contig> list=contigs(a==0 ? left : rc(left), b==0 ? right : rc(right));
				final byte[] before=list.get(0).bases.clone();
				final Evidence counts=new Evidence(k);
				counts.add(genome, 12);
				counts.set(genome.substring(sourceLength+25, sourceLength+25+k), 3);
				counts.set(genome.substring(sourceLength+30, sourceLength+30+k), 24);
				new TadpoleGraph(list, counts, k, 2, 400, 1000).explore(1);
				final Edge e=find(list.get(0), a==0 ? 1 : 0, 1, b);
				check(e!=null, "Missing orientation "+a+","+b+" at k="+k);
				check(e.length==80+k && e.pathMinDepth==3, "Wrong length or bottleneck depth");
				check(e.pathMaxDepth==24 && e.depth==12, "Path maximum confused with first-extension depth");
				check(Math.abs(e.pathMeanDepth-(12+3.0/(e.length+1)))<0.00001, "Wrong path mean");
				final ByteBuilder bb=new ByteBuilder();
				TadpoleDot.edge(e, k, true, bb);
				check(bb.toString().contains("tailport="+(a==0 ? 'e' : 'w')) &&
						bb.toString().contains("headport="+(b==0 ? 'w' : 'e')), "DOT ports disagree with sequence orientation");
				check(Arrays.equals(before, list.get(0).bases), "Explorer modified source bases");
				check(find(list.get(1), b, 0, a==0 ? 1 : 0)!=null, "Missing reverse connection");
			}
		}
	}

	/** Equal branches and duplicate endpoint keys must all survive, deterministically. */
	private static void branches(final int k){
		final String source=random(Math.max(100, k+20), 11);
		final String a=random(Math.max(180, k+100), 12), b=random(Math.max(180, k+100), 13);
		final String[] seq={source, a.substring(80), b.substring(80), a.substring(80), "NNNN", "ACG"};
		final Evidence counts=new Evidence(k);
		counts.add(source+a, 10);
		counts.add(source+b, 10);
		final ArrayList<Contig> one=contigs(seq), many=contigs(seq);
		new TadpoleGraph(one, counts, k, 2, 400, 10000).explore(1);
		new TadpoleGraph(many, counts, k, 2, 400, 10000).explore(3);
		for(int id=1; id<=3; id++){
			check(find(one.get(0), 1, id, 0)!=null, "Lost a branch/shared target "+id+" at k="+k);
		}
		check(describe(one).equals(describe(many)), "Thread count changed graph output");
		check(one.size()==seq.length && one.get(4).leftEdges.isEmpty(), "Isolated short contigs lost");
	}

	/** A supported continuation beyond either bound must be reported as incomplete. */
	private static void limits(final int k){
		final String genome=random(320, 4);
		final Evidence counts=new Evidence(k);
		counts.add(genome, 10);
		final String left=genome.substring(0, 100), right=genome.substring(180);
		final ArrayList<Contig> limited=contigs(left, right);
		new TadpoleGraph(limited, counts, k, 2, 5, 1000).explore(1);
		check(limited.get(0).graphRightLimited && limited.get(0).rightEdges.isEmpty(), "Distance bound silent");
		final ArrayList<Contig> capped=contigs(left, right);
		new TadpoleGraph(capped, counts, k, 2, 400, 1).explore(1);
		check(capped.get(0).graphRightLimited, "State bound silent");
		final ArrayList<Contig> exact=contigs(left, right);
		new TadpoleGraph(exact, counts, k, 2, 80+k, 1000).explore(1);
		check(find(exact.get(0), 1, 1, 0)!=null && !exact.get(0).graphRightLimited, "Exact-bound target falsely limited");
	}

	/** Exact K-base overlaps and circular same-contig overlaps are zero-step links. */
	private static void sharedAnchor(){
		final int k=11;
		final String genome=random(200, 8);
		final Evidence counts=new Evidence(k);
		counts.add(genome, 8);
		final ArrayList<Contig> list=contigs(genome.substring(0, 100), genome.substring(100-k));
		new TadpoleGraph(list, counts, k, 2, 100, 1000).explore(1);
		check(find(list.get(0), 1, 1, 0).length==0, "Shared anchor was skipped");
		check(find(list.get(0), 1, 1, 0).pathMaxDepth==8, "Zero-step maximum must be the anchor depth");
		final String circle=genome.substring(0, 100)+genome.substring(0, k);
		counts.add(circle, 8);
		final ArrayList<Contig> self=contigs(circle);
		new TadpoleGraph(self, counts, k, 2, 100, 1000).explore(1);
		check(find(self.get(0), 1, 0, 0).length==0, "Circular overlap was skipped");
	}

	/** A positive return to a palindromic source end is a self-loop, not an infinite walk. */
	private static void cycle(){
		final Evidence counts=new Evidence(4);
		counts.add("ATATATAT", 7);
		counts.set("TATA", 19);
		final ArrayList<Contig> list=contigs("CCCCATAT");
		new TadpoleGraph(list, counts, 4, 2, 30, 50).explore(1);
		final Edge e=find(list.get(0), 1, 0, 1);
		check(e!=null && e.length==2, "Palindromic positive return missing");
		check(e.pathMinDepth==7 && e.pathMeanDepth==11 && e.pathMaxDepth==19, "Cycle depth range is not measured over its full path");
		check(!list.get(0).graphRightLimited, "Closed cycle incorrectly exhausted search budget");
	}

	/** Missing counts are zero; undefined/short endpoints do not acquire invented links. */
	private static void missingAndShort(){
		final ArrayList<Contig> list=contigs("ACGTNACGT", "ACG", random(80, 9));
		new TadpoleGraph(list, new Evidence(11), 11, 2, 30, 50).explore(1);
		for(Contig c : list){
			check(c.coverage==0 && c.leftEdges.isEmpty() && c.rightEdges.isEmpty(), "Absent evidence became support");
		}
	}

	/** Plain generated-contig DOT is legacy-compatible; pretty names and attributes are explicit. */
	private static void dot(){
		final Contig c=new Contig(new byte[1000], 2);
		c.coverage=12;
		c.minCov=12;
		c.maxCov=24;
		final ByteBuilder bb=new ByteBuilder();
		TadpoleDot.node(c, 7, 31, false, false, 1000, bb);
		final String expected="\t2 [label=\"id=2\\nlen=1000\\ncov=12\\nleft="+
				Tadpole.codeStrings[c.leftCode]+"\\nright="+Tadpole.codeStrings[c.rightCode]+"\"]\n";
		check(bb.toString().equals(expected), "Plain DOT node changed: "+bb);
		bb.clear();
		final ArrayList<Contig> nodes=new ArrayList<Contig>();
		nodes.add(c);
		final double scale=TadpoleDot.logLengthPerInch(nodes);
		TadpoleDot.node(c, 7, 31, true, false, scale, bb);
		check(bb.toString().contains("label=\"contig_9") && bb.toString().contains("width=4,"), "Pretty internal name/scale missing");
		check(!bb.toString().contains("xlabel="), "Pretty label moved outside node");
		check(Math.abs(TadpoleDot.nodeWidth(100, scale)-4*Math.log(16)/Math.log(872))<0.00001,
				"Node width does not follow log(max(16,length-128))");
		check(TadpoleDot.nodeWidth(0, scale)==TadpoleDot.nodeWidth(144, scale), "Short-node log floor changed");
		check(TadpoleDot.nodeWidth(145, scale)>TadpoleDot.nodeWidth(144, scale), "Lengths above the floor must remain distinguishable");
		c.name="sample \"quoted\" \\ literal\tname";
		c.coverage=12.21f;
		bb.clear();
		TadpoleDot.node(c, 0, 31, true, true, scale, bb);
		check(bb.toString().contains("sample \\\"quoted\\\" \\\\ literal\\tname"), "DOT header escaping wrong");
		check(bb.toString().contains("depth=12,12.21,24"), "Node depth summary wrong");
		check(!bb.toString().contains("k=31"), "K redundantly repeated on node");
		final Edge e=new Edge(0, 1, 50, 3, 17, null);
		bb.clear();
		final ByteBuilder legacy=new ByteBuilder();
		e.toDot(legacy);
		TadpoleDot.edge(e, 31, false, bb);
		check(legacy.toString().equals(bb.toString()), "Plain DOT edge changed");
		e.pathMinDepth=12;
		e.pathMeanDepth=12.21;
		e.pathMaxDepth=24;
		bb.clear();
		TadpoleDot.edge(e, 31, true, bb);
		check(bb.toString().contains("depth=12,12.21,24"), "Edge depth summary wrong");
		check(!bb.toString().contains("k=31") && !bb.toString().contains("startdepth="), "Redundant pretty edge labels");
		check(TadpoleDot.penWidth(20)>TadpoleDot.penWidth(10), "Depth did not increase visual edge width");
		check(TadpoleDot.penWidth(0)==0.65 && TadpoleDot.penWidth(Integer.MAX_VALUE)<=1.3, "Edge thickness outside readable bounds");
	}

	/** Writes repeatable small external-assembly and read fixtures for launcher tests. */
	private static void fixtures(final String dir){
		final String source=random(120, 11), a=random(220, 12), b=random(220, 13);
		final ByteBuilder reads=new ByteBuilder();
		for(int i=0; i<12; i++){
			reads.append(">a_").append(i).nl().append(source+a).nl();
			reads.append(">b_").append(i).nl().append(source+b).nl();
		}
		write(dir+"/reads.fa", reads);
		final ByteBuilder assembly=new ByteBuilder();
		assembly.append(">source \"quoted\" \\ literal\n").append(source).nl();
		assembly.append(">branch A\n").append(a.substring(80)).nl();
		assembly.append(">branch B reverse\n").append(rc(b.substring(80))).nl();
		assembly.append(">short\nacg\n>ambiguous\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n");
		write(dir+"/contigs.fa", assembly);
		singleKFixtures(dir);
	}

	/** Writes two reads ending in a shared one-kmer dead-end unitig with two incoming branches. */
	private static void singleKFixtures(final String dir){
		assert(dir!=null) : "The launcher needs an output directory for the single-K assembly fixtures.";
		for(int k : new int[]{31, 33, 62, 127}){
			final String seed=random(k, 93), reverseSeed=rc(seed);
			// Table seeds use the larger strand. Its right end must be dead so
			// makeContig accepts the seed before finding the left-side branch.
			final String center=seed.compareTo(reverseSeed)>0 ? seed : reverseSeed;
			final String first=incomingBranch(center, 'A', false, 101);
			final String second=incomingBranch(center, 'C', true, 103);
			for(int reverse=0; reverse<2; reverse++){
				final ByteBuilder reads=new ByteBuilder();
				reads.append(">first\n").append(reverse==0 ? first : rc(first)).nl();
				reads.append(">second\n").append(reverse==0 ? second : rc(second)).nl();
				write(dir+"/unitig_"+k+"_"+reverse+".fa", reads);
			}
		}
	}

	/** Selects a branch whose assembled arm has the requested canonical strand. */
	private static String incomingBranch(final String center, final char previous,
			final boolean reverse, final long seed){
		assert(center!=null && center.length()>1) : "The shared terminal kmer must leave a K-1 arm overlap.";
		for(int attempt=0; attempt<256; attempt++){
			// Nearby java.util.Random seeds share leading output bits; spread them
			// out so the search can vary the arm's first base as well as its interior.
			final String read=random(center.length()+40, seed+1000003L*attempt)+previous+center;
			// The arm ends one base before the shared kmer; that kmer becomes
			// the separate length-K contig. Force one incoming edge from each end.
			final String arm=read.substring(0, read.length()-1);
			if((arm.compareTo(rc(arm))>0)==reverse){return read;}
		}
		throw new AssertionError("Could not generate both incoming source orientations for K="+center.length());
	}

	/** Checks actual sequence ends independently of Edge orientation bits on a perfect DBG. */
	private static void verifyAdjacentLinks(final String path, final int k){
		assert(k>1) : "A one-step DBG connection must overlap by K-1 bases.";
		final ArrayList<byte[]> lines=ByteFile.toLines(path);
		final HashMap<String, String> sequences=new HashMap<String, String>();
		final LineParser1 lp=new LineParser1('\t');
		for(byte[] line : lines){
			lp.set(line);
			if(lp.termEquals('S', 0)){sequences.put(lp.parseString(1), lp.parseString(2));}
		}
		int checked=0, failed=0, skipped=0;
		for(byte[] line : lines){
			lp.set(line);
			if(!lp.termEquals('L', 0)){continue;}
			int distance=-1;
			for(int i=6; i<lp.terms(); i++){
				if(lp.termStartsWith("EL:i:", i)){distance=lp.parseInt(i, 5);}
			}
			if(distance!=1){skipped++; continue;}
			final String source=sequences.get(lp.parseString(1)), dest=sequences.get(lp.parseString(3));
			check(source!=null && dest!=null, "GFA link refers to a missing segment");
			final String a=lp.termEquals('+', 2) ? source : rc(source);
			final String b=lp.termEquals('+', 4) ? dest : rc(dest);
			check(a.length()>=k && b.length()>=k, "DBG segment shorter than K");
			checked++;
			if(!a.regionMatches(true, a.length()-k+1, b, 0, k-1)){
				failed++;
				System.err.println("Wrong sequence-end overlap: "+lp.parseString(1)+lp.parseString(2)+" -> "+
						lp.parseString(3)+lp.parseString(4)+"; opposite destination fits="+
						a.regionMatches(true, a.length()-k+1, rc(b), 0, k-1));
			}
		}
		System.out.println("Adjacent GFA links checked="+checked+", failed="+failed+", nonadjacent skipped="+skipped);
		check(checked>0 && failed==0, "Actual sequence overlaps disagree with the exported endpoint orientations");
	}

	/** Tests real dense graph workers on a known linear sequence, independently of contig seeding. */
	private static void densePaths(final String dir){
		final String genome=random(1000, 901);
		write(dir+"/dense.fa", new ByteBuilder().append(">genome\n").append(genome).nl());
		int failed=0, checked=0;
		for(int k : new int[]{31, 33, 62, 127}){
			final Tadpole assembly=Tadpole.makeTadpole(new String[]{"in="+dir+"/dense.fa", "k="+k,
					"packed=t", "t=1", "mce=1", "mcs=1", "prefilter=f", "prealloc=f", "initialsize=1000"}, true);
			assembly.loadKmers(new Timer());
			for(int steps : new int[]{1, 2, 5, k+20}){
				final String source=genome.substring(0, 200), target=genome.substring(200-k+steps, 400-k+steps);
				for(int orientation=0; orientation<4; orientation++){
					final boolean sourceRight=(orientation&1)!=0, destRight=(orientation&2)!=0;
					final ArrayList<Contig> list=contigs(sourceRight ? source : rc(source), destRight ? rc(target) : target);
					assembly.initializeContigs(list);
					final AbstractProcessContigThread worker=assembly.makeProcessContigThread(list, new java.util.concurrent.atomic.AtomicInteger());
					worker.processContigs(list);
					final Edge e=find(list.get(0), sourceRight ? 1 : 0, 1, destRight ? 1 : 0);
					final boolean pass=e!=null && e.length==steps && pathOverlap(source, target, e, k) &&
							spellPath(source, target, e, k).equals(genome.substring(200-k, 200+steps));
					System.out.println("Dense path k="+k+", steps="+steps+", orientation="+orientation+", pass="+pass+
							", actual="+e);
					checked++;
					if(!pass){failed++;}
					if(pass){verifyDenseMerge(list, e, k, genome.substring(200, 200+steps), genome.substring(0, 400-k+steps));}
					// This alternative target ends on the walk's forward strand. Reaching
					//its suffix does not mean that the walk can enter its RIGHT end.
					final String outward=genome.substring(50, 200+steps);
					final ArrayList<Contig> wrong=contigs(sourceRight ? source : rc(source), destRight ? outward : rc(outward));
					assembly.initializeContigs(wrong);
					assembly.makeProcessContigThread(wrong, new java.util.concurrent.atomic.AtomicInteger()).processContigs(wrong);
					final Edge invalid=find(wrong.get(0), sourceRight ? 1 : 0, 1, destRight ? 1 : 0);
					final Edge other=find(wrong.get(0), sourceRight ? 1 : 0, 1, destRight ? 0 : 1);
					System.out.println("Outward hit k="+k+", steps="+steps+", orientation="+orientation+", rejected="+(invalid==null));
					checked++;
					if(invalid!=null || other!=null){failed++;}
				}
			}
			assembly.clearData();
		}
		System.out.println("Dense paths checked="+checked+", failed="+failed);
		check(failed==0, "Dense graph workers do not reproduce the known linear path");
	}

	/** Checks raw bridge bytes after source flipping and the complete BubblePopper splice. */
	private static void verifyDenseMerge(final ArrayList<Contig> graph, final Edge e, final int k,
			final String extension, final String expected){
		final Edge copy=new Edge(e.origin, e.destination, e.length, e.orientation, e.depth, e.bases.clone());
		if(!copy.sourceRight()){copy.flipSource();}
		check(new String(copy.bases, StandardCharsets.US_ASCII).equals(extension), "Source flipping corrupted the outward bridge payload");
		final HashMap<Integer, ArrayList<Edge>> inbound=new HashMap<Integer, ArrayList<Edge>>();
		for(Contig c : graph){
			for(int side=0; side<2; side++){
				final ArrayList<Edge> edges=side==0 ? c.leftEdges : c.rightEdges;
				if(edges==null){continue;}
				for(Edge edge : edges){
					ArrayList<Edge> list=inbound.get(edge.destination);
					if(list==null){list=new ArrayList<Edge>(); inbound.put(edge.destination, list);}
					list.add(edge);
				}
			}
		}
		check(new BubblePopper(graph, inbound, k).expand(graph.get(0))==1, "Expected one direct merge of the known path");
		check(new String(graph.get(0).bases, StandardCharsets.US_ASCII).equals(expected), "BubblePopper emitted an incorrect multibase join");
	}

	/** Runs an ordinary long-K assembly and records each edge's independently checked sequence support.
	 * Full-reference absence alone is not failure: repeats can yield valid recombinant DBG walks. */
	private static void verifyPaths(final String reference, final String report, final String[] args){
		final CapturedAssembly assembly=new CapturedAssembly(args);
		// Keep count tables until reciprocal diagnostics finish; process() would
		//release them after this same build/output stage.
		assembly.process2(Tadpole.contigMode);
		check(assembly.graph!=null, "No graph captured; request DOT/GFA with simplification disabled");
		final ArrayList<String> references=new ArrayList<String>();
		for(Read r : ReadInputStream.toReads(reference, FileFormat.FASTA, -1)){
			references.add(new String(r.bases, StandardCharsets.US_ASCII).toUpperCase(java.util.Locale.ROOT));
		}
		check(!references.isEmpty(), "Reference reader returned no sequences");
		final int k=assembly.k();
		final ByteBuilder out=new ByteBuilder();
		out.append("source\tdestination\torientation\tsteps\toverlap_ok\tfull_reference_span\tunsupported_kmers\treciprocal\tspelled_path\n");
		int total=0, longer=0, badOverlap=0, unsupported=0, missingReverse=0;
		for(Contig c : assembly.graph){
			for(int side=0; side<2; side++){
				final ArrayList<Edge> edges=side==0 ? c.leftEdges : c.rightEdges;
				if(edges==null){continue;}
				for(Edge e : edges){
					check(e.origin==c.id && e.sourceRight()==(side==1), "Edge stored on wrong source end: "+e);
					final Contig target=assembly.graph.get(e.destination);
					final String a=oriented(c, !e.sourceRight()), b=oriented(target, e.destRight());
					final boolean overlap=pathOverlap(a, b, e, k);
					final String path=spellPath(a, b, e, k);
					final boolean full=referenceContains(references, path);
					int absent=0;
					if(!full){
						for(int i=0; i+k<=path.length(); i++){
							if(!referenceContains(references, path.substring(i, i+k))){absent++;}
						}
					}
					boolean reciprocal=false;
					final ArrayList<Edge> back=e.destRight() ? target.rightEdges : target.leftEdges;
					if(back!=null){
						for(Edge reverse : back){
							if(reverse.destination==e.origin && reverse.destRight()==e.sourceRight() && reverse.length==e.length){
								final String x=oriented(target, !reverse.sourceRight()), y=oriented(c, reverse.destRight());
								if(pathOverlap(x, y, reverse, k) && spellPath(x, y, reverse, k).equals(rc(path))){reciprocal=true;}
							}
						}
					}
					out.append(e.origin).tab().append(e.destination).tab().append(e.orientation).tab().append(e.length).tab();
					out.append(overlap).tab().append(full).tab().append(absent).tab().append(reciprocal).tab().append(path).nl();
					total++;
					if(e.length>1){longer++;}
					if(!overlap){badOverlap++;}
					if(absent>0){unsupported++;}
					if(!reciprocal){
						missingReverse++;
						final Kmer end=new Kmer(k);
						final int[] counts=new int[4];
						if(e.destRight()){target.rightKmer(end); assembly.tables.fillRightCounts(end, counts);}
						else{target.leftKmer(end); assembly.tables.fillLeftCounts(end, counts);}
						int base=AminoAcid.baseToNumber[rc(path).charAt(k)], max=0;
						if(!e.destRight()){base=3-base;}
						for(int count : counts){max=Math.max(max, count);}
						System.out.println("Reverse seed "+e.origin+"->"+e.destination+": counts="+Arrays.toString(counts)+
								", chosen="+counts[base]+", max="+max+", admitted="+assembly.isJunction(max, counts[base]));
					}
				}
			}
		}
		write(report, out);
		assembly.clearData();
		System.out.println("Path audit: edges="+total+", longer="+longer+", bad_overlap="+badOverlap+
				", unsupported_paths="+unsupported+", missing_reciprocal="+missingReverse);
		// Findings are data, not an excuse to discard the report before isolating a reproducer.
		check(total>0, "Empty assembly graph cannot substantiate a path audit");
	}

	/** Checks the shared source/destination bases implied by a dense graph extension distance. */
	private static boolean pathOverlap(final String source, final String target, final Edge e, final int k){
		assert(source.length()>=k && target.length()>=k) : "Dense graph endpoints must contain their full Kmer.";
		final int shared=Math.max(0, k-e.length);
		return source.regionMatches(source.length()-shared, target, 0, shared);
	}

	/** Independently reconstructs the consumer spelling from source-strand payloads.
	 * Left exits reverse-complement the stored bytes, matching Edge.flipSource. */
	private static String spellPath(final String source, final String target, final Edge e, final int k){
		check(e.length>0 && e.overlap==0 && e.sourceTrim==0 && e.destTrim==0,
				"This diagnostic is for ordinary untrimmed dense edges: "+e);
		check(e.bases!=null && e.bases.length==e.length, "Dense edge payload length disagrees with extension steps: "+e);
		final StringBuilder path=new StringBuilder(source.substring(source.length()-k));
		final int represented=Math.min(k, e.length), novel=e.length-represented;
		for(int i=0; i<novel; i++){
			final byte base=e.bases[e.sourceRight() ? i : e.length-1-i];
			path.append((char)(e.sourceRight() ? base : AminoAcid.baseToComplementExtended[base]));
		}
		path.append(target, k-represented, k);
		return path.toString();
	}

	/** Returns the requested contig strand without changing the graph. */
	private static String oriented(final Contig c, final boolean reverse){
		final String sequence=new String(c.bases, StandardCharsets.US_ASCII);
		return reverse ? rc(sequence) : sequence;
	}

	/** Exact strand-independent reference lookup, suitable for this small diagnostic graph. */
	private static boolean referenceContains(final ArrayList<String> references, final String sequence){
		assert(!references.isEmpty()) : "Reference support requires actual input sequences, not an empty truth set.";
		final String reverse=rc(sequence);
		for(String ref : references){if(ref.contains(sequence) || ref.contains(reverse)){return true;}}
		return false;
	}

	/** Tests short overlaps, long compound payloads, all endpoint orientations, and corrupt paths. */
	private static void pathSpelling(){
		final int k=11;
		final String genome=random(150, 500);
		for(int steps : new int[]{1, 5, 11, 35}){
			final String source=genome.substring(0, 50), target=genome.substring(50-k+steps);
			for(int orientation=0; orientation<4; orientation++){
				final byte[] payload=new byte[steps];
				Arrays.fill(payload, (byte)'N');// Redundant destination payload is not emitted.
				final int novel=Math.max(0, steps-k);
				for(int i=0; i<novel; i++){
					final byte base=(byte)genome.charAt(50+i);
					payload[(orientation&1)==1 ? i : steps-1-i]=(orientation&1)==1 ? base : AminoAcid.baseToComplementExtended[base];
				}
				final Edge e=new Edge(0, 1, steps, orientation, 1, payload);
				final Contig a=contigs((orientation&1)==1 ? source : rc(source)).get(0);
				final Contig b=contigs((orientation&2)==0 ? target : rc(target)).get(0);
				final String x=oriented(a, !e.sourceRight()), y=oriented(b, e.destRight());
				check(pathOverlap(x, y, e, k), "Valid dense overlap rejected");
				check(spellPath(x, y, e, k).equals(genome.substring(50-k, 50+steps)), "Compound payload spelling incorrect");
				if(steps<k){check(!pathOverlap(x, rc(y), e, k), "Wrong destination orientation escaped overlap check");}
				if(novel>0){
					payload[e.sourceRight() ? 0 : e.length-1]='N';
					check(!genome.contains(spellPath(x, y, e, k)), "Corrupt novel base escaped sequence check");
				}
			}
		}
	}

	/** Keeps the graph for diagnostics after process() releases its tables and private list reference. */
	private static final class CapturedAssembly extends Tadpole2 {
		/** Uses ordinary Tadpole2 construction without changing traversal or output settings. */
		CapturedAssembly(final String[] args){super(args, true);}
		@Override
		public void initializeContigs(final ArrayList<Contig> contigs){
			assert(graph==null) : "The diagnostic expects one ordinary graph pass, not multiple K phases.";
			graph=contigs;
			super.initializeContigs(contigs);
		}
		private ArrayList<Contig> graph;
	}

	/** Writes fixture bytes through the BBTools output path. */
	private static void write(final String path, final ByteBuilder data){
		final ByteStreamWriter writer=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.FASTA,
				null, true, false, false, false));
		writer.start();
		writer.print(data);
		check(!writer.poisonAndWait(), "Fixture write failed: "+path);
	}

	/** Builds named contigs in stable input order. */
	private static ArrayList<Contig> contigs(final String... sequences){
		final ArrayList<Contig> list=new ArrayList<Contig>();
		for(String s : sequences){list.add(new Contig(s.getBytes(StandardCharsets.US_ASCII), "input "+list.size(), list.size()));}
		return list;
	}

	/** Finds the specified directed end-to-end connection. */
	private static Edge find(final Contig c, final int side, final int target, final int targetSide){
		final ArrayList<Edge> edges=side==0 ? c.leftEdges : c.rightEdges;
		if(edges==null){return null;}
		for(Edge e : edges){
			if(e.destination==target && e.destRight()==(targetSide==1)){return e;}
		}
		return null;
	}

	/** Produces deterministic text for complete graph parity comparisons. */
	private static String describe(final ArrayList<Contig> list){
		final ByteBuilder bb=new ByteBuilder();
		for(Contig c : list){
			TadpoleDot.node(c, 0, 11, true, true, 1000, bb);
			for(Edge e : c.leftEdges){TadpoleDot.edge(e, 11, true, bb);}
			for(Edge e : c.rightEdges){TadpoleDot.edge(e, 11, true, bb);}
		}
		return bb.toString();
	}

	/** Returns an exact reverse complement for orientation fixtures. */
	private static String rc(final String bases){return AminoAcid.reverseComplementBases(bases);}

	/** Produces deterministic DNA without fixture files in the source package. */
	private static String random(final int length, final long seed){
		final Random r=new Random(seed);
		final StringBuilder sb=new StringBuilder(length);
		for(int i=0; i<length; i++){sb.append("ACGT".charAt(r.nextInt(4)));}
		return sb.toString();
	}

	/** A failed test must terminate even when assertions have been disabled accidentally. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Small exact canonical count map for tests, independent of the production tables. */
	private static final class Evidence implements TadpoleGraph.Counts {
		/** Keeps the requested exact word length. */
		Evidence(final int k_){k=k_;}
		/** Inserts every word of a fixture at a known depth. */
		void add(final String sequence, final int depth){
			for(int i=0; i+k<=sequence.length(); i++){set(sequence.substring(i, i+k), depth);}
		}
		/** Changes a canonical word's known count. */
		void set(final String word, final int depth){map.put(canonical(word), depth);}
		@Override
		public int count(final Kmer key){return map.get(canonical(key.toString()));}
		/** Independent canonicalization for fixture count lookup. */
		private String canonical(final String s){
			final String reverse=rc(s);
			return s.compareTo(reverse)>0 ? s : reverse;
		}
		final int k;
		final ObjectIntMap<String> map=new ObjectIntMap<String>(String.class);
	}
}
