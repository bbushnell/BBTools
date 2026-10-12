package cardinality;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import ukmer.Kmer;

/**
 * Estimates depth-frequency histograms at several K values in one input pass.
 * Each worker owns one DynamicDemiLog per K; only completed sketches are merged.
 * Histogram counts represent distinct canonical kmers, not reads or genomic bases.
 * Based on template.A_SampleStreamerMT.
 * @author Brian Bushnell, Fischl
 */
public class KmerHistMulti implements Accumulator<KmerHistMulti.ProcessThread>{

	/** Runs the command-line tool, closing redirected diagnostics on completion. */
	public static void main(String[] args){
		final Timer timer=new Timer();
		final KmerHistMulti tool=new KmerHistMulti(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses standard BBTools input options and histogram-specific settings. */
	public KmerHistMulti(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		parser=parse(args);
		Parser.processQuality();
		validateParams();
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		if(parser.in1==null){throw new IllegalArgumentException("Specify in=<reads>.");}
		if(parser.in2==null && parser.in1.indexOf('#')>=0 && !new File(parser.in1).exists()){
			parser.in2=parser.in1.replace("#", "2");
			parser.in1=parser.in1.replace("#", "1");
		}
		parser.in1=Tools.fixExtension(parser.in1);
		parser.in2=Tools.fixExtension(parser.in2);
		if(parser.in2!=null){FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;}
		if(parser.out1==null){parser.out1="stdout.txt";}
		if(parser.out2!=null || parser.append){
			throw new IllegalArgumentException("One histogram output is supported; out2 and append are not supported.");
		}
		Tools.testForDuplicateFiles(true, parser.in1, parser.in2, parser.qfin1, parser.qfin2, parser.out1);
		if(!Tools.testOutputFiles(parser.overwrite, false, false, parser.out1)){
			throw new IllegalArgumentException("Cannot write histogram: "+parser.out1);
		}
		ffin1=FileFormat.testInput(parser.in1, FileFormat.FASTQ, parser.extin, true, true);
		ffin2=FileFormat.testInput(parser.in2, FileFormat.FASTQ, parser.extin, true, true);
		ffout=FileFormat.testOutput(parser.out1, FileFormat.TEXT, parser.extout,
			true, parser.overwrite, false, false);
	}

	/** Applies options in order, including config-expanded command-line overrides. */
	private Parser parse(final String[] args){
		final Parser p=new Parser();
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){
				kValues=parseK(b);
			}else if(a.equals("buckets")){
				buckets=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("minprob")){
				minProb=Float.parseFloat(b);
			}else if(a.equals("histmax")){
				histMax=Parse.parseIntKMG(b);
			}else if(a.equals("raw")){
				raw=Parse.parseBoolean(b);
			}else if(a.equals("zeros")){
				zeros=Parse.parseBoolean(b);
			}else if(a.equals("khist") || a.equals("hist")){
				p.out1=b;
			}else if(p.parse(arg, a, b)){
				//Standard input, quality, threads, memory and overwrite handling.
			}else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){
				p.in1=arg;
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		assert(kValues.length>0) : "At least one sketch is required to produce histogram columns.";
		return p;
	}

	/** Validates allocation sizes and prevents silently relabeling normalized K values. */
	private void validateParams(){
		if(buckets<1 || buckets>0x40000000){
			throw new IllegalArgumentException("buckets must be in [1,1073741824]: "+buckets);
		}
		buckets=CardinalityTracker.powerOf2AtLeast(buckets);
		if(histMax<1 || histMax>Character.MAX_VALUE){
			throw new IllegalArgumentException("histmax must be in [1,65535]; DDL counters saturate at65535.");
		}
		if(!(minProb>=0 && minProb<=1)){
			throw new IllegalArgumentException("minprob must be finite and in [0,1]: "+minProb);
		}
		if(parser.samplerate!=1){
			throw new IllegalArgumentException("Subsample reads before khistmulti; samplerate would change the measured depth distribution.");
		}
		if(seed<0){seed=new java.util.Random().nextLong()&Long.MAX_VALUE;}
		for(final int k : kValues){
			if(Kmer.getKbig(k)!=k){
				throw new IllegalArgumentException("K="+k+" would be rounded by the selected packing mode; use packed=t.");
			}
		}
		assert(Integer.bitCount(buckets)==1) : "DDL indexes buckets with a bit mask: "+buckets;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes each input once; no histogram is opened until input completes successfully. */
	private void process(final Timer timer){
		final int workers=parser.workers();
		assert(workers>0) : "Streamer consumption requires at least one worker.";
		final Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, parser.qfin1, parser.qfin2,
			true, parser.maxReads, false, true, parser.threadsIn);
		final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){
			threads.add(new ProcessThread(st, createSketches(kValues, buckets, seed, minProb)));
		}
		outstream.println("K values: "+Arrays.toString(kValues)+"; buckets per sketch: "+buckets+"; seed: "+seed);
		outstream.println("Sketch array bytes: "+(5L*buckets*kValues.length*workers)+
			" (excludes reader, object and histogram overhead)");
		try{
			st.start();
			errorState|=!ThreadWaiter.startAndWait(threads, this);
		}finally{
			errorState|=ReadWrite.closeStream(st);
		}
		if(errorState){throw new RuntimeException("Input processing failed; no histogram was written.");}
		writeHistogram();
		timer.stop();
		outstream.println(Tools.timeReadsBasesProcessed(timer, readsProcessed, basesProcessed, 8));
	}

	/** Merges worker-local sketches only after all workers have terminated. */
	@Override
	public void accumulate(final ProcessThread pt){
		assert(pt.getState()==Thread.State.TERMINATED) : "Merging a live DDL would race with count updates.";
		readsProcessed+=pt.reads;
		basesProcessed+=pt.bases;
		errorState|=!pt.success;
		if(merged==null){merged=pt.sketches;}
		else{mergeSketches(merged, pt.sketches);}
	}

	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	/** Writes aligned depth bins; each last bin contains all depths at least histMax. */
	private void writeHistogram(){
		assert(merged!=null && merged.length==kValues.length) : "Output columns must match the merged K array.";
		final long[][] hist=new long[kValues.length][];
		final long[] cardinalities=new long[kValues.length];
		final double[] scales=new double[kValues.length];
		for(int i=0; i<merged.length; i++){
			hist[i]=histogram(merged[i], histMax);
			cardinalities[i]=merged[i].cardinality();
			final int filled=merged[i].filledBuckets();
			scales[i]=filled==0 ? 0 : cardinalities[i]/(double)filled;
		}
		final ByteStreamWriter writer=new ByteStreamWriter(ffout);
		writer.start();
		final ByteBuilder line=new ByteBuilder();
		writer.print("#Estimator\tDynamicDemiLog\n#Count_units\testimated_distinct_canonical_kmers\n");
		writer.print("#Scaling\tcardinality/occupied_buckets; raw counts are sampled bucket winners\n");
		writer.print("#Depth_units\tobservations_per_kmer; not_read_or_genomic_depth\n");
		writer.print("#Hash_seed\t").print(seed).nl();
		writer.print("#Minprob\t").print(minProb, 6).nl();
		writer.print("#Buckets_per_K\t").print(buckets).nl();
		writer.print("#Reads\t").print(readsProcessed).print("\n#Bases\t").print(basesProcessed).nl();
		writer.print("#Last_bin_ge\t").print(histMax).nl();
		writer.print("#Counters_saturate_at\t65535\n");
		for(int i=0; i<merged.length; i++){
			writer.print("#K\t").print(kValues[i]).print("\tCardinality\t").print(cardinalities[i]);
			writer.print("\tOccupied_buckets\t").print(merged[i].filledBuckets()).nl();
		}
		line.append("#Depth");
		for(final int k : kValues){
			line.append("\tCount_K").append(k);
			if(raw){line.append("\tRaw_K").append(k);}
		}
		writer.print(line.nl());
		for(int depth=1; depth<=histMax; depth++){
			boolean present=zeros;
			for(final long[] h : hist){present|=h[depth]>0;}
			if(!present){continue;}
			line.clear().append(depth);
			for(int i=0; i<hist.length; i++){
				line.tab().append(hist[i][depth]*scales[i], 3);
				if(raw){line.tab().append(hist[i][depth]);}
			}
			writer.print(line.nl());
		}
		if(writer.poisonAndWait()){throw new RuntimeException("Histogram output failed: "+parser.out1);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses positive, exact K values, preserving requested column order. */
	static int[] parseK(final String text){
		if(text==null || text.isEmpty()){throw new IllegalArgumentException("k requires a comma-separated list.");}
		final String[] terms=text.split(",", -1);
		final int[] result=new int[terms.length];
		for(int i=0; i<terms.length; i++){
			result[i]=Integer.parseInt(terms[i]);
			if(result[i]<1){throw new IllegalArgumentException("K must be positive: "+result[i]);}
			for(int j=0; j<i; j++){
				if(result[i]==result[j]){throw new IllegalArgumentException("Duplicate K: "+result[i]);}
			}
		}
		assert(result.length>0) : "A histogram needs at least one K column.";
		return result;
	}

	/** Allocates worker-private sketches with identical hash settings in each worker. */
	static DynamicDemiLog[] createSketches(final int[] ks, final int buckets, final long seed, final float minProb){
		assert(ks.length>0 && seed>=0) : "Resolve a common seed before constructing worker sketches.";
		final DynamicDemiLog[] sketches=new DynamicDemiLog[ks.length];
		for(int i=0; i<ks.length; i++){
			sketches[i]=DynamicDemiLog.create(buckets, ks[i], seed, minProb, true, false);
		}
		return sketches;
	}

	/** Combines corresponding sketches without retaining a second aggregate copy. */
	static void mergeSketches(final DynamicDemiLog[] target, final DynamicDemiLog[] source){
		assert(target.length==source.length) : "A merge must preserve the number of histogram columns.";
		for(int i=0; i<target.length; i++){
			assert(target[i]!=source[i] && target[i].k==source[i].k && target[i].buckets==source[i].buckets
				&& target[i].hashXor==source[i].hashXor && target[i].minProb==source[i].minProb) :
				"DDL merge requires distinct sketches with identical K, buckets, hash seed and filtering; column="+i;
			target[i].add(source[i]);
		}
	}

	/** Counts sampled bucket winners; empty buckets never contribute depth-zero kmers. */
	static long[] histogram(final DynamicDemiLog sketch, final int max){
		assert(max>0 && max<=Character.MAX_VALUE) : "Histogram bins must fit DDL's saturating 16-bit counts: "+max;
		final long[] hist=new long[max+1];
		long sum=0;
		for(final char count : sketch.counts16()){
			if(count>0){hist[Math.min(count, max)]++; sum++;}
		}
		assert(sum==sketch.filledBuckets()) : "Every occupied DDL bucket must contribute once; counts="+
			sum+", occupied="+sketch.filledBuckets();
		return hist;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes shared read batches, keeping all sketch writes local. */
	static class ProcessThread extends Thread{
		/** Retains the stream and the worker's exclusively owned sketches. */
		ProcessThread(final Streamer st_, final DynamicDemiLog[] sketches_){
			st=st_;
			sketches=sketches_;
		}

		/** A failed consumer must terminate loudly rather than leave producer queues blocked. */
		@Override
		public void run(){
			try{
				processInner();
				success=true;
			}catch(Throwable failure){
				failure.printStackTrace();
				KillSwitch.kill("khistmulti worker failed; no histogram is valid.");
			}
		}

		/** Empty batches are not EOF; hash(Read) already includes the mate. */
		private void processInner(){
			assert(sketches.length>0) : "Each worker needs at least one histogram sketch.";
			for(ListNum<Read> batch=st.nextList(); batch!=null; batch=st.nextList()){
				for(final Read r : batch){
					if(!r.validated()){r.validate(true);}
					if(r.mate!=null && !r.mate.validated()){r.mate.validate(true);}
					reads+=r.pairCount();
					bases+=r.pairLength();
					for(final DynamicDemiLog sketch : sketches){sketch.hash(r);}
				}
			}
		}

		private final Streamer st;
		final DynamicDemiLog[] sketches;
		long reads=0, bases=0;
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Parser parser;
	private final FileFormat ffin1, ffin2, ffout;
	private final PrintStream outstream;
	private int[] kValues={32, 64, 96, 124};
	private int buckets=262144;
	private int histMax=65535;
	private long seed=0;
	private float minProb=0;
	private boolean raw=false, zeros=false;
	private boolean errorState=false;
	private DynamicDemiLog[] merged;
	private long readsProcessed=0, basesProcessed=0;
}
