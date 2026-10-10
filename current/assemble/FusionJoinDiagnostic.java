package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Development-only census of selected fusion windows and local read-depth features.
 * Reference sequence supplies conservative adjacency labels, never depth features.
 * A single exact canonical table serves one query K; no assembly is modified.
 * @author Fischl
 */
public final class FusionJoinDiagnostic {

	/** Loads one count table, measures all logged pairs, then releases that table. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, FusionJoinDiagnostic.class, false);
		try{
			if(pp.args.length==2 && pp.args[0].equals("neuraltest")){
				FusionNeuralGateTest.run(pp.args[1]);
				return;
			}
			if(pp.args.length>0 && pp.args[0].equals("entropy")){
				FusionEntropyCorpus.run(Arrays.copyOfRange(pp.args, 1, pp.args.length));
				return;
			}
			if(pp.args.length>0 && pp.args[0].equals("entropytest")){
				FusionTipEntropy.selfTest();
				return;
			}
			if(pp.args.length>0 && pp.args[0].equals("plan")){
				FusionCorpusPlan.run(Arrays.copyOfRange(pp.args, 1, pp.args.length));
				return;
			}
			if(pp.args.length>0 && pp.args[0].equals("label")){
				FusionJoinLabeler.run(Arrays.copyOfRange(pp.args, 1, pp.args.length));
				return;
			}
			if(pp.args.length>=1 && pp.args[0].equals("test")){
				selfTest();
				if(pp.args.length==2){writeFixture(pp.args[1]);}
				return;
			}
			final FusionJoinDiagnostic x=new FusionJoinDiagnostic(pp.args);
			x.process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Separates diagnostic paths from the native table loader's standard options. */
	private FusionJoinDiagnostic(final String[] args){
		final ArrayList<String> tableArgs=new ArrayList<String>();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(a.equals("trace")){trace=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("assembly")){assembly=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("outvectors")){outVectors=b;}
			else if(a.equals("outhist")){outHist=b;}
			else{tableArgs.add(arg);}
		}
		if(trace==null || reference==null || assembly==null || out==null){
			throw new IllegalArgumentException("Required: in=reads trace=log ref=reference assembly=fasta out=table k=N");
		}
		if(!Tools.testInputFiles(false, true, trace, reference, assembly) ||
				!Tools.testOutputFiles(false, false, false, out, outVectors, outHist) ||
				!Tools.testForDuplicateFiles(true, trace, reference, assembly, out, outVectors, outHist)){
			throw new IllegalArgumentException("Diagnostic paths are missing, duplicate, or would overwrite output.");
		}
		// The diagnostic deliberately uses unfiltered exact counts, not compact hashes.
		tableArgs.add("prefilter=f");
		tableArgs.add("prealloc=f");
		tableArgs.add("hashonly=f");
		tableArgs.add("rcomp=t");
		tableArgs.add("minprob=0");
		tableArgs.add("minprobmain=f");
		tableArgs.add("packed=t");
		tables=new KmerTableSetU(tableArgs.toArray(new String[0]), 0);
		counts=new TadpoleGraph.Counts(){
			@Override
			public int count(final Kmer word){return tables.getCount(word);}
		};
		final ArrayList<String> inputs=new ArrayList<String>(tables.in1);
		inputs.addAll(tables.in2);
		inputs.addAll(tables.extra);
		inputs.add(trace);
		inputs.add(reference);
		// Reuse graph I/O protection, including physical aliases and stdout variants.
		TadpoleGraph.checkPaths(assembly, inputs, out, outVectors, outHist);
		k=tables.kbig;
		key=new Kmer(k);
		neighbor=new Kmer(k);
	}

	/** Borrows immutable phase counts; this instance never owns, loads or clears a table. */
	FusionJoinDiagnostic(final int k_, final TadpoleGraph.Counts counts_){
		assert(k_>0 && counts_!=null) : "Live fusion features require the current phase's count provider.";
		tables=null;
		counts=counts_;
		k=k_;
		key=new Kmer(k);
		neighbor=new Kmer(k);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes one row per selected pair, retaining unresolved labels and empty regions. */
	private void process(){
		final ArrayList<Join> joins=readJoins(trace);
		final ArrayList<String> refs=sequences(reference), contigs=sequences(assembly);
		final ByteStreamWriter[] writers=new ByteStreamWriter[3];
		boolean completed=false, writeError=false;
		try{
			tables.process(new Timer());
			final String[] paths={out, outVectors, outHist};
			for(int i=0; i<paths.length; i++){
				if(paths[i]==null){continue;}
				writers[i]=new ByteStreamWriter(paths[i], false, false, false);
				writers[i].start();
			}
			writeRows(joins, refs, contigs, writers);
			completed=true;
		}finally{
			// A malformed selected window must not leave a non-daemon writer waiting.
			for(ByteStreamWriter writer : writers){
				if(writer!=null){writeError=writer.poisonAndWait() | writeError;}
			}
			tables.clear();
			if(writeError && !completed){System.err.println("Join output also failed while handling the original error.");}
		}
		if(writeError){throw new RuntimeException("Join diagnostic output failed: "+out);}
		System.err.println("JOIN_DIAGNOSTIC_PASS pairs="+joins.size()+" k="+k);
	}

	/** Serializes the census while the caller owns unconditional writer/table cleanup. */
	private void writeRows(final ArrayList<Join> joins, final ArrayList<String> refs,
			final ArrayList<String> contigs, final ByteStreamWriter[] writers){
		final ByteStreamWriter writer=writers[0];
		final ByteBuilder row=new ByteBuilder();
		if(writers[1]!=null){writeVectorHeader(writers[1]);}
		if(writers[2]!=null){writers[2].print("event\tphase_k\tquery_k\tregion\tdepth\tpositions\n");}
		writeDiagnosticHeader(writer);
		for(Join join : joins){
			row.clear();
			final Hit whole=find(refs, join.product), left=find(refs, join.left), right=find(refs, join.right);
			final String label=label(join, refs, whole, left, right);
			row.append(join.event).tab().append(join.phase).tab().append(k).tab().append(join.source).tab();
			row.append(join.dest).tab().append(join.overlap).tab().append(join.sourceTrim).tab().append(join.destTrim);
			row.tab().append(join.left.length-join.overlap).tab().append(join.right.length-join.overlap);
			row.tab().append(label).tab().append(whole.count).tab().append(left.count).tab().append(right.count);
			row.tab().append(find(contigs, join.product).count);
			features(join, row, writers[1], writers[2]);
			row.nl();
			writer.print(row);
		}
	}

	/** Shares the exact numeric schema between offline and live-table extraction. */
	static void writeVectorHeader(final ByteStreamWriter writer){
		assert(writer!=null) : "A started writer is required for the feature schema.";
		final ByteBuilder row=new ByteBuilder();
		row.append("#schema=").append(FusionJoinFeatures.VERSION).nl().append("event\tphase_k\tquery_k");
		for(String name : FusionJoinFeatures.names()){row.tab().append(name);}
		writer.print(row.nl());
	}

	/** Retains the legacy diagnostic column order for feature-level parity tests. */
	static void writeDiagnosticHeader(final ByteStreamWriter writer){
		assert(writer!=null) : "A started writer is required for the diagnostic schema.";
		final ByteBuilder row=new ByteBuilder();
		row.append("event\tphase_k\tquery_k\tsource_id\tdest_id\toverlap\tsource_trim\tdest_trim");
		row.append("\tleft_bases\tright_bases\tlabel\twhole_hits\tleft_hits\tright_hits\tfinal_hits");
		for(String name : REGION_NAMES){
			row.append('\t').append(name).append("_n\t").append(name).append("_min\t").append(name);
			row.append("_p10\t").append(name).append("_median\t").append(name).append("_p90\t");
			row.append(name).append("_max");
		}
		row.append("\toverlap_left_ratio\toverlap_right_ratio\toverlap_highflank_ratio\tflank_ratio");
		row.append("\twindow_words\tmissing_words\tbranch_tests\tmax_alt_ratio\talt_ge_half\talt_ge_expected\n");
		writer.print(row);
	}

	/**
	 * Separates pure flank/anchor words from boundary-crossing words. The spanning
	 * subset crosses both anchor boundaries. Neighbor ratios use the overlap plus
	 * K bases per side, avoiding unrelated distant branches in the 200-base context.
	 */
	void features(final Join join, final ByteBuilder row,
			final ByteStreamWriter vectors, final ByteStreamWriter histograms){
		final float[] valuesOut=measure(join.product, join.product.length, join.left.length-join.overlap,
				join.overlap, join.sourceTrim, join.destTrim, row);
		if(vectors!=null){
			sidecar.clear().append(join.event).tab().append(join.phase).tab().append(k);
			for(float value : valuesOut){sidecar.tab().append(value, 9);}
			vectors.print(sidecar.nl());
		}
		if(histograms!=null){writeHistograms(join, values, sizes, histograms);}
	}

	/**
	 * Measures a retained window into reusable buffers for both collection and inference.
	 * The returned vector is borrowed until the next call; a null row skips text output.
	 */
	float[] measure(final byte[] bases, final int length, final int lo, final int overlap,
			final int sourceTrim, final int destTrim, final ByteBuilder row){
		final int n=length-k+1, hi=lo+overlap;
		assert(n>0 && lo>=0 && hi<=length && length<=bases.length) :
				"Retained join geometry must contain a complete query word inside the populated buffer.";
		if(values[0].length<n){values=new int[5][n];}
		Arrays.fill(sizes, 0);
		int missing=0, tested=0, half=0, dominated=0;
		double maxRatio=0;
		key.clear();
		for(int end=0; end<length; end++){
			if(AminoAcid.baseToNumber[bases[end]]<0){
				throw new IllegalArgumentException("Undefined base in retained fusion window at "+end);
			}
			key.addRight(bases[end]);
			if(end<k-1){continue;}
			final int start=end-k+1, depth=Math.max(0, counts.count(key));
			if(depth==0){missing++;}
			final int region=region(start, end+1, lo, hi);
			values[region][sizes[region]++]=depth;
			if(start<lo && end+1>hi){values[4][sizes[4]++]=depth;}
			final int first=Math.max(0, lo-k), last=Math.min(length, hi+k);
			if(start<first || end>=last){continue;}
			for(int direction=0; direction<2; direction++){
				final boolean forward=direction==1;
				final int next=forward ? end+1 : start-1;
				if(next<first || next>=last){continue;}
				final int expected=AminoAcid.baseToNumber[bases[next]];
				int intended=0, alternate=0;
				for(int b=0; b<4; b++){
					neighbor.setFrom(key);
					if(forward){neighbor.addRightNumeric(b);}else{neighbor.addLeftNumeric(b);}
					final int count=Math.max(0, counts.count(neighbor));
					if(b==expected){intended=count;}else{alternate=Math.max(alternate, count);}
				}
				tested++;
				maxRatio=Math.max(maxRatio, alternate/(double)Math.max(1, intended));
				if(alternate>0 && 2L*alternate>=intended){half++;}
				if(alternate>0 && alternate>=intended){dominated++;}
			}
		}
		for(int i=0; i<5; i++){stats(values[i], sizes[i], row);}
		if(row!=null){
			final double l=quantile(values[0], sizes[0], .5), o=quantile(values[1], sizes[1], .5);
			final double r=quantile(values[2], sizes[2], .5);
			ratio(o, l, row);
			ratio(o, r, row);
			ratio(o, Math.max(l, r), row);
			ratio(Math.max(l, r), Math.min(l, r), row);
			row.tab().append(n).tab().append(missing).tab().append(tested).tab().append(maxRatio, 6);
			row.tab().append(half).tab().append(dominated);
		}
		return encoder.fill(values, sizes, k, overlap, lo, length-hi,
				sourceTrim, destTrim, n, missing, tested, half, dominated, maxRatio);
	}

	/** Saves all observed depth bins; an empty region has NA depth and zero positions. */
	private void writeHistograms(final Join join, final int[][] values, final int[] sizes,
			final ByteStreamWriter writer){
		assert(values.length==REGION_NAMES.length) : "Histogram regions must match the diagnostic schema.";
		for(int r=0; r<values.length; r++){
			int start=0;
			do{
				int end=start;
				while(end<sizes[r] && values[r][end]==values[r][start]){end++;}
				sidecar.clear().append(join.event).tab().append(join.phase).tab().append(k);
				sidecar.tab().append(REGION_NAMES[r]).tab();
				if(sizes[r]==0){sidecar.append("NA");}else{sidecar.append(values[r][start]);}
				sidecar.tab().append(end-start).nl();
				writer.print(sidecar);
				start=end;
			}while(start<sizes[r]);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Regions are disjoint; the both-flank-spanning class is an explicit extra subset. */
	private static int region(final int start, final int end, final int lo, final int hi){
		assert(start<end && lo<=hi) : "Word and overlap coordinates use half-open intervals.";
		if(end<=lo){return 0;}
		if(start>=hi){return 2;}
		if(start>=lo && end<=hi){return 1;}
		return 3;
	}

	/** Emits sample size plus nearest-rank-on-index quantiles; absence is NA, not zero. */
	private static void stats(final int[] values, final int n, final ByteBuilder row){
		Arrays.sort(values, 0, n);
		if(row==null){return;}
		row.tab().append(n);
		for(double q : QUANTILES){
			row.tab();
			if(n==0){row.append("NA");}else{row.append((int)quantile(values, n, q));}
		}
	}

	/** Uses floor(q*(n-1)); callers sort only populated array entries first. */
	private static double quantile(final int[] values, final int n, final double q){
		assert(n>=0 && n<=values.length && q>=0 && q<=1) : "Depth quantiles need valid populated ranges.";
		return n==0 ? Double.NaN : values[(int)(q*(n-1))];
	}

	/** Undefined/zero denominators stay NA rather than manufacturing an enrichment. */
	private static void ratio(final double a, final double b, final ByteBuilder row){
		row.tab();
		if(Double.isNaN(a) || Double.isNaN(b) || b<=0){row.append("NA");}
		else{row.append(a/b, 6);}
	}

	/** Parses native byte lines; only the few retained selected windows allocate arrays. */
	static ArrayList<Join> readJoins(final String path){
		return readJoins(path, false);
	}

	/** Empty traces are legitimate for corpus assemblies without reciprocal candidates. */
	static ArrayList<Join> readJoins(final String path, final boolean allowEmpty){
		final ArrayList<Join> joins=new ArrayList<Join>();
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 parser=new LineParser1('\t');
		int phase=-1;
		for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
			parser.set(line);
			if(parser.termEquals("FusionPass", 0)){phase=parser.parseInt(1);}
			else if(parser.termEquals("OverlapPair", 0)){
				if(parser.terms()!=11 || phase<1){throw new IllegalArgumentException("Unrecognized fusion trace record.");}
				joins.add(new Join(joins.size()+1, phase, parser.parseInt(1), parser.parseInt(2),
						parser.parseInt(5), parser.parseInt(6), parser.parseInt(7),
						parser.parseByteArray(9), parser.parseByteArray(10)));
			}
		}
		if(file.close() || (!allowEmpty && joins.isEmpty())){throw new IllegalArgumentException("Unreadable/empty fusion trace: "+path);}
		return joins;
	}

	/** Retains both orientations of small truth/assembly sequences, not a second kmer table. */
	static ArrayList<String> sequences(final String path){
		final ArrayList<String> answer=new ArrayList<String>();
		final boolean force=FASTQ.FORCE_INTERLEAVED, test=FASTQ.TEST_INTERLEAVED;
		final boolean split=FastaReadInputStream.SPLIT_READS;
		final ArrayList<Read> reads;
		try{
			// Reference contigs are not the paired reads configured for the count table.
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			FastaReadInputStream.SPLIT_READS=false;
			reads=ReadInputStream.toReads(path, FileFormat.FASTA, -1);
		}finally{
			FASTQ.FORCE_INTERLEAVED=force;
			FASTQ.TEST_INTERLEAVED=test;
			FastaReadInputStream.SPLIT_READS=split;
		}
		for(Read read : reads){
			assert(read.mate==null) : "Diagnostic reference loading must retain every unpaired contig.";
			answer.add(new String(read.bases, StandardCharsets.US_ASCII).toUpperCase(java.util.Locale.ROOT));
			final byte[] reverse=read.bases.clone();
			AminoAcid.reverseComplementBasesInPlace(reverse);
			answer.add(new String(reverse, StandardCharsets.US_ASCII).toUpperCase(java.util.Locale.ROOT));
		}
		if(answer.isEmpty()){throw new IllegalArgumentException("No reference/assembly sequences: "+path);}
		return answer;
	}

	/** Counts all overlapping matches on both strands; palindromes conservatively count twice. */
	static Hit find(final ArrayList<String> sequences, final byte[] query){
		assert(query.length>0) : "An empty join arm would match every reference position.";
		final String word=new String(query, StandardCharsets.US_ASCII);
		final Hit hit=new Hit();
		for(int sequence=0; sequence<sequences.size(); sequence++){
			final String text=sequences.get(sequence);
			for(int pos=text.indexOf(word); pos>=0; pos=text.indexOf(word, pos+1)){
				hit.count++;
				hit.sequence=sequence;
				hit.pos=pos;
			}
		}
		return hit;
	}

	/** Only uniquely placed, incompatible arms on one reference establish a contradiction. */
	static String label(final Join join, final ArrayList<String> refs,
			final Hit whole, final Hit left, final Hit right){
		//TODO: Global-truth limitation: these bounded windows cannot establish repeat-copy identity.
		//Shouchella g0491 has two supported 495 bp joins through one 1185 bp repeat contig,
		//but their combined 1585 bp window is absent. Do not interpret "supported" as global chain truth.
		if(whole.count>0){return "supported";}
		if(left.count!=1 || right.count!=1 || left.sequence/2!=right.sequence/2){return "unresolved";}
		if(left.sequence!=right.sequence){return "unresolved";}
		final int delta=right.pos-left.pos-(join.left.length-join.overlap);
		// Without reference topology metadata, a possible origin wrap is not a false join.
		if(delta%refs.get(left.sequence).length()==0){return "unresolved";}
		return "contradicted";
	}

	/** Tests trimming, exact overlap, orientation-aware labels, empty regions and quantiles. */
	private static void selfTest(){
		final Join join=new Join(1, 3, 0, 1, 3, 2, 2, bytes("AACGTAANN"), bytes("NNTAACTGC"));
		check(Arrays.equals(join.product, bytes("AACGTAACTGC")), "Discarded tails entered the joined sequence");
		check(region(0, 4, 4, 7)==0 && region(4, 7, 4, 7)==1 && region(7, 10, 4, 7)==2 &&
				region(3, 8, 4, 7)==3, "Half-open region boundaries changed");
		check(Double.isNaN(quantile(new int[1], 0, .5)), "Empty anchor depth must be NA");
		check(quantile(new int[]{1, 2, 3, 4}, 4, .5)==2, "Quantile definition changed");
		final ArrayList<String> refs=new ArrayList<String>();
		refs.add("AACGTAACTGC");
		refs.add("GCAGTTACGTT");
		check(label(join, refs, find(refs, join.product), find(refs, join.left), find(refs, join.right)).equals("supported"),
				"Exact whole adjacency was not supported");
		refs.set(0, "AACGTAAGGGTAACTGC");
		refs.set(1, "GCAGTTACCCTTACGTT");
		check(label(join, refs, find(refs, join.product), find(refs, join.left), find(refs, join.right)).equals("contradicted"),
				"Unique displaced arms were not contradicted");
		refs.add("TAACTGC");
		refs.add("GCAGTTA");
		check(label(join, refs, find(refs, join.product), find(refs, join.left), find(refs, join.right)).equals("unresolved"),
				"Repeated arms must not become certain bad labels");
		boolean rejected=false;
		try{new Join(1, 3, 0, 1, 3, 0, 0, bytes("AACGTAA"), bytes("TATCTGC"));}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Nonexact anchor accepted");
		FusionJoinFeaturesTest.test();
		FusionJoinLabeler.selfTest();
		System.err.println("JOIN_DIAGNOSTIC_TEST_PASS");
	}

	/** Converts literal fixture sequences only. */
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}

	/** Writes a tiny exact-count CLI fixture, including an overlap shorter than all query Ks. */
	private static void writeFixture(final String directory){
		final java.util.Random random=new java.util.Random(191);
		final ByteBuilder genome=new ByteBuilder();
		for(int i=0; i<600; i++){genome.append("ACGT".charAt(random.nextInt(4)));}
		final String sequence=genome.toString();
		writeText(directory+"/reference.fa", ">reference\n"+sequence+"\n");
		final StringBuilder reads=new StringBuilder();
		for(int i=0; i<5; i++){reads.append('>').append(i).append('\n').append(sequence).append('\n');}
		writeText(directory+"/reads.fa", reads.toString());
		final StringBuilder trace=new StringBuilder("FusionPass\t32\t94\tfalse\n");
		final StringBuilder reverse=new StringBuilder("FusionPass\t32\t94\tfalse\n");
		for(int overlap : new int[]{64, 20}){
			trace.append("OverlapPair\t0\t1\t").append(300+overlap).append("\t300\t");
			trace.append(overlap).append("\t0\t0\t100\t").append(sequence.substring(100, 300+overlap));
			trace.append('\t').append(sequence.substring(300, 500+overlap)).append('\n');
			final byte[] source=bytes(sequence.substring(100, 300+overlap));
			final byte[] dest=bytes(sequence.substring(300, 500+overlap));
			AminoAcid.reverseComplementBasesInPlace(source);
			AminoAcid.reverseComplementBasesInPlace(dest);
			reverse.append("OverlapPair\t0\t1\t300\t").append(300+overlap).append('\t').append(overlap);
			reverse.append("\t0\t0\t100\t").append(new String(dest, StandardCharsets.US_ASCII));
			reverse.append('\t').append(new String(source, StandardCharsets.US_ASCII)).append('\n');
		}
		writeText(directory+"/trace.tsv", trace.toString());
		writeText(directory+"/trace_reverse.tsv", reverse.toString());
		FusionJoinCollector.selfTest(directory);
	}

	/** Writes small fixture text through native checked I/O. */
	private static void writeText(final String path, final String text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		writer.print(text);
		if(writer.poisonAndWait()){throw new RuntimeException("Could not write fixture: "+path);}
	}

	/** Fixture failures remain loud even if a caller mistakenly disables assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One selected edge with original trace identity and retained oriented sequences. */
	static final class Join {
		/** Reconstructs precisely the retained source plus the nonoverlapping destination. */
		Join(final int event_, final int phase_, final int source_, final int dest_, final int overlap_,
				final int st, final int dt, final byte[] sourceBases, final byte[] destBases){
			event=event_;
			phase=phase_;
			source=source_;
			dest=dest_;
			overlap=overlap_;
			sourceTrim=st;
			destTrim=dt;
			if(st<0 || dt<0 || overlap<1 || sourceBases.length-st<=overlap || destBases.length-dt<=overlap){
				throw new IllegalArgumentException("Invalid retained geometry for event "+event);
			}
			left=Arrays.copyOf(sourceBases, sourceBases.length-st);
			right=Arrays.copyOfRange(destBases, dt, destBases.length);
			for(int i=0; i<overlap; i++){
				if(left[left.length-overlap+i]!=right[i]){throw new IllegalArgumentException("Nonexact anchor at event "+event);}
			}
			product=Arrays.copyOf(left, left.length+right.length-overlap);
			System.arraycopy(right, overlap, product, left.length, right.length-overlap);
		}
		final int event, phase, source, dest, overlap, sourceTrim, destTrim;
		final byte[] left, right, product;
	}

	/** A match census; coordinates are used only when exactly one hit exists. */
	static final class Hit {
		int count, sequence=-1, pos=-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String trace, reference, assembly, out, outVectors, outHist;
	private final FusionJoinFeatures encoder=new FusionJoinFeatures();
	private final ByteBuilder sidecar=new ByteBuilder();
	private int[][] values=new int[5][0];
	private final int[] sizes=new int[5];
	private final KmerTableSetU tables;
	private final TadpoleGraph.Counts counts;
	private final Kmer key, neighbor;
	private final int k;
	private static final String[] REGION_NAMES={"left", "overlap", "right", "boundary", "span"};
	private static final double[] QUANTILES={0, .1, .5, .9, 1};
}
