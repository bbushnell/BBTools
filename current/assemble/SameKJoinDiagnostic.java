package assemble;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Development census of same-table contig boundaries. Collection uses raw
 * eligible single-K unitigs and changes no join decisions. Reference labeling
 * is a separate invocation after the read table has been released.
 * @author Fischl
 */
public final class SameKJoinDiagnostic {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Supports config= through the native preparser, without a shell command string. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKJoinDiagnostic.class, false);
		try{new SameKJoinDiagnostic(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns only diagnostic flags; collect forwards the remaining flags to Tadpole. */
	private SameKJoinDiagnostic(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=eq<0 ? null : arg.substring(eq+1);
			if(name.equals("task")){task=value;}
			else if(name.equals("out")){output=value;}
			else if(name.equals("rawcontigs")){rawContigs=value;}
			else if(name.equals("features")){features=value;}
			else if(name.equals("controls")){controls=value;}
			else if(name.equals("controlcount")){controlCount=Parse.parseIntKMG(value);}
			else if(name.equals("controlseed")){controlSeed=Long.parseLong(value);}
			else if(name.equals("ref")){reference=value;}
			else if(name.equals("circular")){circular=Parse.parseBoolean(value);}
			else if(name.equals("context")){context=Parse.parseIntKMG(value);}
			else{nativeArgs.add(arg);}
		}
		if(output==null || output.equalsIgnoreCase("null") || context<1 ||
				(!task.equals("collect") && !task.equals("label"))){
			throw new IllegalArgumentException("Required out=; task=collect|label; context must be positive.");
		}
		if(controlCount< -1 || (controls!=null && features==null) || (controlCount!=-1 && controls==null)){
			throw new IllegalArgumentException("controls= requires features=; controlcount=-1 means auto, otherwise requires controls=.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps read-depth collection and reference truth in separate process lifetimes. */
	private void run(){
		assert(output!=null) : "Argument validation must precede opening any diagnostic output.";
		if(task.equals("label")){label();}else{collect();}
	}

	/** Runs the existing contig walk once, then borrows that same live table. */
	private void collect(){
		if(reference!=null || circular){throw new IllegalArgumentException("Truth reference/topology belongs only to task=label.");}
		final ArrayList<String> args=new ArrayList<String>(nativeArgs);
		for(String arg : args){
			final String key=arg.split("=", 2)[0].replaceFirst("^-+", "").toLowerCase(Locale.ROOT);
			if(key.startsWith("out") || key.equals("dot") || key.equals("gfa") || key.equals("khist") || key.equals("hist")
					|| key.equals("dump") || key.equals("ihist") || key.equals("inserthistogram")){
				throw new IllegalArgumentException("Use only out= and rawcontigs= for diagnostic output, not "+key);
			}
		}
		if(TadpoleMulti.hasMultipleK(args.toArray(new String[0]))){
			throw new IllegalArgumentException("This diagnostic requires a single K, not a multi-K schedule.");
		}
		// These are a documented raw-unitig diagnostic mode, not shipping defaults.
		args.addAll(Arrays.asList("mode=contig", "out=null", "showstats=f", "processcontigs=f", "pop=f",
				"sweeplen=0", "retainshortcontigs=t", "trimends=0", "trimcircular=f", "hashkmers=f"));
		final Tadpole tad=Tadpole.makeTadpole(args.toArray(new String[0]), true);
		if(tad.processContigs || tad.evictLowDepthContigs || !tad.tables().rcomp() ||
				kmer.AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN){
			throw new IllegalArgumentException("Census requires raw, canonical DNA unitigs without graph cleanup/exports.");
		}
		final ArrayList<String> inputs=new ArrayList<String>(tad.tables().in1);
		inputs.addAll(tad.tables().in2); inputs.addAll(tad.tables().extra);
		checkPaths(inputs, output, rawContigs, features, controls);
		try{
			tad.process2(Tadpole.contigMode);
			if(tad.errorState || tad.tables().errorState){throw new IllegalStateException("Assembly failed before candidate collection.");}
			final ArrayList<Contig> contigs=tad.detachContigs();
			if(rawContigs!=null){writeContigs(contigs);}
			final int k=tad.k();
			final TadpoleGraph.Counts counts=new TadpoleGraph.Counts(){
				@Override public int count(final Kmer word){return tad.bridgeCount(word);}
			};
			final long momentStart=System.nanoTime();
			final SameKJoinFeatures.Moments moments=features==null ? null : SameKJoinFeatures.moments(tad);
			if(moments!=null){
				outstream.println("SAME_K_MOMENTS unique="+moments.unique+" observations="+moments.observations+
						" squares="+moments.squares+" seconds="+((System.nanoTime()-momentStart)*1e-9));
			}
			try(final SameKJoinExport export=features==null ? null : new SameKJoinExport(features, controls, k, context, moments)){
				final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
				writer.start();
				long events;
				try{
					writer.print(HEADER+"\n");
					final ByteBuilder row=new ByteBuilder();
					final SameKJoinCollector collector=new SameKJoinCollector(contigs, counts, k);
					events=collector.collect(new SameKJoinCollector.Sink(){
						@Override
						public void accept(final long event, final int se, final int de, final Contig a, final boolean ar,
								final Contig b, final boolean br, final int[] ac, final int[] bc){
							final int overlap=k-1;
							final byte[] left=SameKJoinCollector.context(a, ar, (int)Math.max(0L, (long)a.length()-overlap-context), a.length());
							final byte[] right=SameKJoinCollector.context(b, br, 0, (int)Math.min(b.length(), (long)overlap+context));
							row.clear().append(event).tab().append(k).tab().append(se/2).tab().append(se&1);
							row.tab().append(de/2).tab().append(de&1).tab().append(a.length()).tab().append(b.length());
							row.tab().append(a.coverage, 6).tab().append(b.coverage, 6);
							row.tab().append(ar ? a.leftCode : a.rightCode).tab().append(br ? b.rightCode : b.leftCode);
							row.tab().append(SameKJoinCollector.context(a, ar, a.length()-k, a.length()));
							row.tab().append(SameKJoinCollector.context(b, br, 0, k));
							row.tab().append(tad.isJunction(max(ac, 1), second(ac, 1), max(ac, 5), second(ac, 5)) ? 1 : 0);
							row.tab().append(tad.isJunction(max(bc, 1), second(bc, 1), max(bc, 5), second(bc, 5)) ? 1 : 0);
							for(int value : ac){row.tab().append(value);}
							for(int value : bc){row.tab().append(value);}
							writer.print(row.tab().append(left).tab().append(right).nl());
							if(export!=null){export.candidate(event, a, ar, b, br, ac, bc);}
						}
					});
				}finally{
					if(writer.poisonAndWait()){throw new RuntimeException("Candidate census output failed: "+output);}
				}
				if(controls!=null){
					final int requested=controlCount<0 ? Math.toIntExact(events/9) : controlCount;
					final SameKAcceptedControls.Result result=export.controls(contigs, counts, requested, controlSeed);
					outstream.println("SAME_K_CONTROLS_PASS population="+result.population+" sampled="+result.sampled+
							" probability="+result.probability()+" seed="+controlSeed);
				}
				outstream.println("SAME_K_CENSUS_PASS unitigs="+contigs.size()+" candidates="+events+" k="+k);
			}
		}finally{tad.clearData();}
	}

	/** Outputs all retained raw unitigs, irrespective of the final mincontig policy. */
	private void writeContigs(final ArrayList<Contig> contigs){
		assert(rawContigs!=null) : "Raw sequence output is explicit, not a side effect of collecting vectors.";
		final ByteStreamWriter writer=new ByteStreamWriter(rawContigs, false, false, false);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		try{
			for(int i=0; i<contigs.size(); i++){
				row.clear().append(">unitig_").append(i).append(" cov=").append(contigs.get(i).coverage, 6);
				writer.print(row.nl().append(contigs.get(i).bases).nl());
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Raw contig output failed.");}}
	}

	/** Labels bounded pair context and, separately, censuses every occurrence of the endpoint Kmers. */
	private void label(){
		final Parser parser=new Parser();
		for(String arg : nativeArgs){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length<2 ? null : split[1];
			if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown label argument: "+arg);}
		}
		if(reference==null || parser.in1==null || rawContigs!=null || features!=null || controls!=null){
			throw new IllegalArgumentException("Labeling requires in=candidates.tsv ref=genome.fa out=labels.tsv.");
		}
		checkPaths(new ArrayList<String>(Arrays.asList(parser.in1, reference)), output);
		final ArrayList<FusionJoinDiagnostic.Join> joins=new ArrayList<FusionJoinDiagnostic.Join>();
		final ArrayList<byte[]> words=new ArrayList<byte[]>(), patterns=new ArrayList<byte[]>();
		final ByteFile input=ByteFile.makeByteFile(parser.in1, true);
		final LineParser1 lp=new LineParser1('\t');
		int k=0;
		boolean controlInput=false;
		try{
			final byte[] header=input.nextLine();
			final String schema=header==null ? "" : new String(header, StandardCharsets.US_ASCII);
			controlInput=SameKJoinExport.CONTROL_HEADER.equals(schema);
			if(!controlInput && !HEADER.equals(schema)){
				throw new IllegalArgumentException("Unrecognized same-K candidate schema.");
			}
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=(controlInput ? 10 : 36)){throw new IllegalArgumentException("Wrong same-K context row width.");}
				final int current=lp.parseInt(1);
				if(current<2 || (k!=0 && k!=current)){throw new IllegalArgumentException("Label one K per input file.");}
				k=current;
				final int event=lp.parseInt(0);
				if(event!=joins.size()+1){throw new IllegalArgumentException("Candidate IDs must be consecutive and unique.");}
				final FusionJoinDiagnostic.Join join=new FusionJoinDiagnostic.Join(event, k,
						lp.parseInt(2), controlInput ? lp.parseInt(2) : lp.parseInt(4), k-1, 0, 0,
						lp.parseByteArray(controlInput ? 8 : 34), lp.parseByteArray(controlInput ? 9 : 35));
				joins.add(join);
				words.add(lp.parseByteArray(controlInput ? 6 : 12)); words.add(lp.parseByteArray(controlInput ? 7 : 13));
				patterns.add(join.product); patterns.add(join.left); patterns.add(join.right);
			}
		}finally{if(input.close()){throw new RuntimeException("Error reading candidate census.");}}
		final ArrayList<String> refs=FusionJoinDiagnostic.sequences(reference);
		final ArrayList<byte[]> physicalRefs=new ArrayList<byte[]>();
		for(int i=0; i<refs.size(); i+=2){physicalRefs.add(refs.get(i).getBytes(StandardCharsets.US_ASCII));}
		final boolean[] topology=new boolean[physicalRefs.size()];
		Arrays.fill(topology, circular);
		final KmerBranchOracle.Result[] local=KmerBranchOracle.scan(words, physicalRefs, topology, k==0 ? 2 : k);
		final FusionJoinDiagnostic.Hit[] contextHits=FusionJoinLabeler.findAll(refs, patterns);
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		int supported=0, contradicted=0, unresolved=0;
		try{
			writer.print("#pair labels use bounded exact context, not whole-assembly truth; ambiguous placements remain unresolved\n");
			writer.print((controlInput ? "control_event" : "event")+"\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask\n");
			for(int i=0; i<joins.size(); i++){
				final FusionJoinDiagnostic.Join join=joins.get(i);
				final FusionJoinDiagnostic.Hit whole=contextHits[3*i], left=contextHits[3*i+1], right=contextHits[3*i+2];
				final String label=pairLabel(join, refs, whole, left, right, circular);
				final String target;
				if(label.equals("supported")){target="1"; supported++;}
				else if(label.equals("contradicted")){target="0"; contradicted++;}
				else{target="NA"; unresolved++;}
				row.clear().append(join.event).tab().append(k).tab().append(label).tab().append(target);
				row.tab().append(whole.count).tab().append(left.count).tab().append(right.count);
				appendTruth(row, local[2*i]); appendTruth(row, local[2*i+1]);
				writer.print(row.nl());
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Reference label output failed.");}}
		outstream.println((controlInput ? "SAME_K_CONTROL_LABEL_PASS" : "SAME_K_LABEL_PASS")+" candidates="+joins.size()+" supported="+supported+
				" contradicted="+contradicted+" unresolved="+unresolved);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Unique compatible context supports a pair; absent context alone never establishes a false join. */
	static String pairLabel(final FusionJoinDiagnostic.Join join, final ArrayList<String> refs,
			final FusionJoinDiagnostic.Hit whole, final FusionJoinDiagnostic.Hit left,
			final FusionJoinDiagnostic.Hit right, final boolean circular){
		assert(join!=null && whole!=null && left!=null && right!=null) : "Pair labels require whole and both arm hit classes.";
		if(whole.count==1 && left.count==1 && right.count==1){return "supported";}
		// Circular-origin context is not padded in this initial pair oracle. Do not
		// manufacture negatives from missing wraparound matches; local masks do wrap.
		if(circular){return "unresolved";}
		if(whole.count==0){return FusionJoinDiagnostic.label(join, refs, whole, left, right);}
		return "unresolved";
	}

	/** Raw truth is isolated from the candidate feature file. */
	private static void appendTruth(final ByteBuilder row, final KmerBranchOracle.Result truth){
		assert(truth!=null) : "Each recorded endpoint must have one reference-census result, even with zero hits.";
		row.tab().append(truth.occurrences).tab().append(truth.unknownFlanks);
		row.tab().append(truth.leftMask).tab().append(truth.rightMask).tab().append(truth.pairMask);
	}

	/** Validates aliases and refuses existing outputs, regardless of native overwrite flags. */
	private static void checkPaths(final ArrayList<String> inputs, final String... outputs){
		if(inputs.isEmpty()){throw new IllegalArgumentException("At least one input must be specified.");}
		TadpoleGraph.checkPaths(inputs.get(0), inputs, outputs);
		if(!Tools.testInputFiles(false, true, inputs.toArray(new String[0])) ||
				!Tools.testOutputFiles(false, false, false, outputs)){
			throw new IllegalArgumentException("Inputs must exist and diagnostic outputs must be new writable files.");
		}
	}

	/** Four-base maximum; no sorting or allocation in the candidate writer. */
	private static int max(final int[] depths, final int start){
		assert(start==1 || start==5) : "Nine-depth endpoint schema has two four-base blocks.";
		return Math.max(Math.max(depths[start], depths[start+1]), Math.max(depths[start+2], depths[start+3]));
	}

	/** Counts ties as a real second-highest branch. */
	private static int second(final int[] depths, final int start){
		assert(start==1 || start==5) : "Four extension counts must come from one endpoint side.";
		int first=0, second=0;
		for(int i=start; i<start+4; i++){
			if(depths[i]>=first){second=first; first=depths[i];}
			else if(depths[i]>second){second=depths[i];}
		}
		return second;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final ArrayList<String> nativeArgs=new ArrayList<String>();
	private String task="collect", output, rawContigs, reference, features, controls;
	private boolean circular;
	private int context=200;
	private int controlCount=-1;
	private long controlSeed=1;
	static final String HEADER="event\tk\tsource\tsource_end\tdest\tdest_end\tsource_len\tdest_len\tsource_cov\tdest_cov\tsource_code\tdest_code\tsource_word\tdest_word\tsource_junction\tdest_junction"+
			"\tsource_depth\tsource_rA\tsource_rC\tsource_rG\tsource_rT\tsource_lA\tsource_lC\tsource_lG\tsource_lT"+
			"\tdest_depth\tdest_rA\tdest_rC\tdest_rG\tdest_rT\tdest_lA\tdest_lC\tdest_lG\tdest_lT\tleft_context\tright_context";
}
