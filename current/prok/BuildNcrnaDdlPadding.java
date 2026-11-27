package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Learns the measured padding blacklist from current training-side windows.
 * Lifecycle derives from template.A_SampleBasic2; native whole-input iteration
 * follows RankNcrnaWindows. Existing window analysis/oracle own all verification.
 * Counts only flanks outside ALL accepted model spans, never annotation flanks.
 * Frequency scratch is leased to active inputs and survives between files.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaDdlPadding implements CallGenesIteration{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BuildNcrnaDdlPadding.class, false);
		try{new BuildNcrnaDdlPadding(pp.args, pp.outstream).process(new Timer());}
		finally{Shared.closeStream(pp.outstream);}
	}

	public BuildNcrnaDdlPadding(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){callerConfig=b;}
			else if(a.equals("sample")){sample=b;}
			else if(a.equals("binding")){binding=b;}
			else if(a.equals("reuse")){reuse=b;}
			else if(a.equals("first")){first=Parse.parseIntKMG(b);}
			else if(a.equals("count")){count=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmodels")){models=Parse.parseIntKMG(b);}
			else if(a.equals("fit")){fit=b;}
			else if(a.equals("makelist")){makeList=Parse.parseBoolean(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown padding parameter: "+arg);}
		}
		prefix=parser.out1;
		require(prefix!=null && callerConfig!=null && sample!=null && binding!=null && first>=0 && (count==-1 || count>0)
			&& threads>0 && models>0 && fit!=null && !fit.isEmpty() && fit.indexOf('\t')<0
			&& parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1,
			"Use caller=flat_config sample=training.tsv first=0 count=N out=fresh_prefix with complete contigs");
		require(Tools.testInputFiles(false, true, callerConfig, sample, binding), "Caller, training sample and worker-owned execution binding must be readable");
		for(final String suffix:new String[]{".frequencies.tsv.gz", ".inputs.tsv", ".summary.tsv", ".models.fa", ".sample.tsv", ".binding.tsv"}){
			require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Padding outputs must be fresh");
		}
		require(!makeList || Tools.testOutputFiles(false, false, false, prefix+".padding.tsv"), "The learned padding list must be fresh");
		NcrnaDdlWindowAnalysis.checkEnvironment();
		readSample();
		if(reuse!=null){readReuse();}
		require(!makeList || prior.size()+current.size()==all.size(),
			"A final padding list requires every row of the complete training sample exactly once");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Inputs and Run      ----------------*/
	/*--------------------------------------------------------------*/

	private void readSample(){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(sample)){
			RrnaResourceIO.header(in.reader, NcrnaDdlHoldout.PANEL_HEADER, sample);
			int ordinal=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==5 && p.parseInt(0)==ordinal, "The prepared sample must retain its complete ordered source schema");
				final String id=p.parseString(1), input=p.parseString(2), fasta=p.parseString(3);
				require(!sampleIndex.contains(id) && !samplePaths.contains(fasta) && fasta.startsWith("/"),
					"Each training assembly and literal FASTA may occur once in the prepared sample");
				sampleIndex.put(id, all.size());
				samplePaths.put(fasta, all.size());
				final Sample entry=new Sample(id, input, fasta);
				all.add(entry);
				if(ordinal>=first && (count<0 || (long)ordinal<first+(long)count)){
					currentByPath.put(fasta, current.size());
					current.add(entry);
				}
				ordinal++;
			}
		}
		require(first<all.size() && current.size()==(count<0 ? all.size()-first : count),
			"The explicit input slice must exist completely; never truncate a requested preparation cell");
	}

	private void readReuse(){
		// The worker supplies source/build, effective caller options and every loaded
		// resource binding. Compare before importing counts or calling any new input.
		requireSameText(sample, reuse+".sample.tsv", "Retained pilot must bind the COMPLETE ordered training sample");
		requireSameText(binding, reuse+".binding.tsv", "Retained pilot must bind the same worker-verified caller resources, options and source/build");
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(reuse+".inputs.tsv")){
			RrnaResourceIO.header(in.reader, INPUT_HEADER, reuse+".inputs.tsv");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==10, "Retained preparation inputs require the complete count/CPU schema");
				final int index=sampleIndex.get(p.parseString(0));
				require(index>=0, "A retained pilot cannot introduce an assembly outside the same prepared sample");
				final Sample entry=all.get(index);
				require(!entry.finished && !currentByPath.contains(entry.fasta) && p.termEquals(entry.input, 1) && p.termEquals(entry.fasta, 2),
					"Retained pilot and new preparation inputs must be disjoint and preserve exact assembly/path identity");
				entry.contigs=p.parseLong(3);
				entry.bases=p.parseLong(4);
				entry.windows=p.parseLong(5);
				entry.accepted=p.parseLong(6);
				entry.shortSegments=p.parseLong(7);
				entry.longSegments=p.parseLong(8);
				entry.oracleCpu=p.parseLong(9);
				require(entry.contigs>=0 && entry.bases>=0 && entry.windows>=entry.accepted && entry.accepted>=0
					&& entry.shortSegments>=entry.longSegments && entry.longSegments>=0 && entry.oracleCpu>=0,
					"Retained pilot counts and CPU must be nonnegative and conserve accepted windows");
				entry.finished=true;
				prior.add(entry);
			}
		}
		require(!prior.isEmpty(), "Retained pilot counts require a nonempty completed input audit");
		try(Input in=new Input(reuse+".frequencies.tsv.gz")){
			RrnaResourceIO.header(in.reader, FREQUENCY_HEADER, reuse+".frequencies.tsv.gz");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==3 && (p.termEquals("SHORT8", 0) || p.termEquals("MASKED10", 0)), "Unknown retained word-frequency schema/form");
				final int[] counts=p.termEquals("SHORT8", 0) ? shortCounts : longCounts;
				final int raw=p.parseInt(1), frequency=p.parseInt(2);
				require(raw>=0 && raw<counts.length && counts[raw]==0 && frequency>0
					&& (counts==shortCounts || (raw&LONG_MASK)==raw),
					"Retained word counts must be unique, positive and within their actual masked domain");
				counts[raw]=frequency;
			}
		}
	}

	private void process(final Timer timer){
		final String[] config=Parser.parseConfig(new String[]{"config="+NcrnaCallerProcess.configurationPath(new java.io.File(callerConfig))});
		final ArrayList<String> arguments=new ArrayList<String>();
		for(final String arg:config){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			require(!a.equals("config") && !a.equals("cobalaminddl") && !a.equals("cobalamingc")
				&& !a.equals("reads") && !a.equals("maxreads") && !a.equals("samplerate"),
				"Padding uses a flat, complete-contig ordinary caller config without DDL or hidden sampling");
			if(!a.equals("in") && !a.equals("in1")){arguments.add(arg);}
		}
		for(final Sample entry:current){arguments.add("in="+entry.fasta);}
		arguments.add("out=null");
		arguments.add("outstats=null");
		arguments.add("ow=f");
		arguments.add("taxonomy=f");
		arguments.add("t="+threads);
		final CallGenes caller=new CallGenes(arguments.toArray(new String[arguments.size()]));
		caller.run(timer, ()->{
			bind(caller.inputFiles());
			caller.setIteration(this);
		});
		for(final Sample entry:current){require(entry.finished, "Every selected training input must finish before learning a list");}
		for(final Counts counts:owners){
			for(int raw=0; raw<shortCounts.length; raw++){shortCounts[raw]=Math.addExact(shortCounts[raw], counts.shortCounts[raw]);}
			for(int raw=0; raw<longCounts.length; raw++){longCounts[raw]=Math.addExact(longCounts[raw], counts.longCounts[raw]);}
		}
		writeResults();
		outstream.println("NCRNA_DDL_PADDING_COMPLETE");
		outstream.println("Padding uses current-model accepted-window envelopes; retained pilot inputs were not recalled.");
	}

	private void bind(final String[] inputs){
		require(inputs!=null && inputs.length==current.size(), "The native caller must consume exactly the selected training slice");
		for(int i=0; i<inputs.length; i++){require(inputs[i].equals(current.get(i).fasta), "Native input order differs from the explicit training slice");}
		for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
			final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(i);
			if(!f.name.equals("cobalamin_rf00174")){continue;}
			require(family==null && f.library.length==models, "Bind exactly the declared current main-model roster");
			family=f;
			targetIndex=i;
		}
		require(family!=null, "The ordinary caller configuration must enable the current RF00174 family");
		if(reuse!=null){
			final NcrnaSupport.Library saved=NcrnaSupport.loadLibrary(reuse+".models.fa");
			require(saved.names.length==models, "Retained pilot model count differs from the current library");
			for(int model=0; model<models; model++){
				require(RrnaResourceIO.first(saved.names[model]).equals(RrnaResourceIO.first(family.modelNames[model]))
					&& Arrays.equals(saved.bases[model], family.library[model]),
					"Retained pilot frequencies must bind the same ordered consensus identities and bytes");
			}
		}
	}

	@Override
	public Worker create(final String input, final GeneCaller caller, final GeneModel model){
		final int index=currentByPath.get(input);
		require(index>=0 && !current.get(index).finished, "One native worker owns one declared unfinished training input");
		return new InputWorker(current.get(index), caller, model);
	}

	private synchronized Counts acquire(){
		if(!available.isEmpty()){return available.remove(available.size()-1);}
		final Counts result=new Counts();
		owners.add(result);
		return result;
	}

	private synchronized void release(final Counts counts){available.add(counts);}

	/*--------------------------------------------------------------*/
	/*----------------        Native Input Worker   ----------------*/
	/*--------------------------------------------------------------*/

	private final class InputWorker implements Worker, NcrnaDdlWindowAnalysis.Sink{
		InputWorker(final Sample sample_, final GeneCaller caller_, final GeneModel model_){
			entry=sample_;
			caller=caller_;
			model=model_;
			counts=acquire();
			analysis=new NcrnaDdlWindowAnalysis(entry.id, family, caller.ncrnaScavenger(targetIndex), this);
		}
		@Override
		public ArrayList<Orf> call(final Read read){
			analysis.begin(read);
			try{return caller.callGenes(read, model, true);}
			finally{
				analysis.end();
				entry.contigs++;
				entry.bases+=read.length();
			}
		}
		@Override
		public void window(final NcrnaDdlWindowAnalysis view, final byte[] window){
			entry.windows++;
			final NcrnaDdlOracle oracle=view.oracle();
			int left=Integer.MAX_VALUE, right=-1;
			for(int m=0; m<models; m++){
				if(oracle.start(m)<0){continue;}
				left=Math.min(left, oracle.start(m));
				right=Math.max(right, oracle.stop(m));
			}
			if(right<0){return;}
			entry.accepted++;
			final int start=view.strand()==0 ? view.start() : view.contigLength()-1-view.stop();
			final int stop=view.strand()==0 ? view.stop() : view.contigLength()-1-view.start();
			final long shortBefore=counts.shortSegments, longBefore=counts.longSegments;
			addPadding(counts, window, start, stop, view.strand(), left, right, bounds);
			entry.shortSegments+=counts.shortSegments-shortBefore;
			entry.longSegments+=counts.longSegments-longBefore;
		}
		@Override
		public void attempt(final NcrnaDdlWindowAnalysis view, final int rank, final int model,
				final int shared, final float identity, final int span, final boolean passed){}
		@Override
		public void endWindow(final NcrnaDdlWindowAnalysis view){}
		@Override
		public void finish(){
			analysis.finish();
			entry.oracleCpu=analysis.oracleCpu();
			entry.finished=true;
			release(counts);
		}
		private final Sample entry;
		private final GeneCaller caller;
		private final GeneModel model;
		private final NcrnaDdlWindowAnalysis analysis;
		private final Counts counts;
		private final int[] bounds=new int[2];
	}

	/** One pair of primitive count arrays per concurrent input lease. Native file
	 * workers are recreated between inputs; the pool retains counts across leases. */
	public static final class Counts implements CardinalityTracker.KmerConsumer{
		public void flank(final byte[] window, final int from, final int to){
			assert(window!=null && from>=0 && to>=from && to<=window.length)
				: KillSwitch.assertDie("Padding flanks must be separate valid slices outside the accepted window envelope");
			if(to-from<8){return;}
			shortSegments++;
			if(to-from>=10){longSegments++;}
			CardinalityTracker.forEachKmer(Arrays.copyOfRange(window, from, to), 8, 10, false, true, this);
		}
		@Override
		public void accept(final long word, final int k){
			assert(k==8 || k==10) : KillSwitch.assertDie("Padding receives only the declared strict8/10word forms");
			final int raw=(int)(k==8 ? word : word&LONG_MASK);
			final int[] counts=k==8 ? shortCounts : longCounts;
			assert(raw>=0 && raw<counts.length && counts[raw]<Integer.MAX_VALUE)
				: KillSwitch.assertDie("Padding occurrence counts must fit their raw-word domain and nonoverflowing integer counter");
			counts[raw]++;
		}
		public final int[] shortCounts=new int[1<<16], longCounts=new int[1<<20];
		public long shortSegments, longSegments;
	}

	/** Shared measured envelope conversion, with the two outside flanks counted
	 * separately. A literal finite fixture checks strand conversion and no bridge. */
	public static void addPadding(final Counts counts, final byte[] window, final int start, final int stop,
			final int strand, final int acceptedStart, final int acceptedStop, final int[] bounds){
		require(counts!=null && window!=null && window.length==stop-start+1,
			"Padding words must come from the exact captured strand-oriented window");
		BuildNcrnaDdlLayoutCounts.paddingBounds(start, stop, strand, acceptedStart, acceptedStop, bounds);
		counts.flank(window, 0, bounds[0]);
		counts.flank(window, bounds[1], window.length);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void writeResults(){
		copyText(sample, prefix+".sample.tsv");
		copyText(binding, prefix+".binding.tsv");
		final long[] shortWords=BuildNcrnaDdlLayoutCounts.orderedFrequencies(shortCounts);
		final long[] longWords=BuildNcrnaDdlLayoutCounts.orderedFrequencies(longCounts);
		try(Output out=new Output(prefix+".frequencies.tsv.gz")){
			out.writer.println(FREQUENCY_HEADER);
			frequencies(out, "SHORT8", shortWords);
			frequencies(out, "MASKED10", longWords);
		}
		if(makeList){
			try(Output out=new Output(prefix+".padding.tsv")){
				out.writer.println(NcrnaDdlSketchBuilder.BLACKLIST_HEADER);
				blacklist(out, "SHORT8", 8, shortWords, 256);
				blacklist(out, "MASKED10", 10, longWords, 2048);
			}
		}
		long windows=0, accepted=0, oracle=0;
		try(Output out=new Output(prefix+".inputs.tsv")){
			out.writer.println(INPUT_HEADER);
			final ByteBuilder row=new ByteBuilder();
			for(final Sample entry:all){
				if(!entry.finished){continue;}
				out.writer.print(row.clear().append(entry.id).tab().append(entry.input).tab().append(entry.fasta)
					.tab().append(entry.contigs).tab().append(entry.bases).tab().append(entry.windows).tab().append(entry.accepted)
					.tab().append(entry.shortSegments).tab().append(entry.longSegments).tab().append(entry.oracleCpu).nl());
				windows+=entry.windows;
				accepted+=entry.accepted;
				oracle+=entry.oracleCpu;
			}
		}
		try(Output out=new Output(prefix+".models.fa")){
			final ByteBuilder row=new ByteBuilder();
			for(int model=0; model<models; model++){
				out.writer.print(row.clear().append('>').append(RrnaResourceIO.first(family.modelNames[model]))
					.nl().append(family.library[model]).nl());
			}
		}
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.println("new_inputs\tretained_inputs\tmodels\twindows\taccepted_windows\tshort_words\tlong_words\toracle_cpu_ns\tblacklist_emitted\tfit");
			out.writer.print(new ByteBuilder().append(current.size()).tab().append(prior.size()).tab().append(models)
				.tab().append(windows).tab().append(accepted).tab().append(shortWords.length).tab().append(longWords.length)
				.tab().append(oracle).tab().append(makeList).tab().append(fit).nl());
		}
	}

	private static void frequencies(final Output out, final String form, final long[] words){
		final ByteBuilder row=new ByteBuilder();
		for(int rank=0; rank<words.length; rank++){
			out.writer.print(row.clear().append(form).tab().append(BuildNcrnaDdlLayoutCounts.rawWord(words, rank))
				.tab().append(words[words.length-1-rank]>>>32).nl());
		}
	}

	/** Literal binding equality; no digest or implicit claim that paths prove bytes.
	 * The experiment worker owns creation and verification of the execution binding. */
	private static void requireSameText(final String current, final String previous, final String why){
		try(Input a=new Input(current); Input b=new Input(previous)){
			for(byte[] row=a.reader.nextLine(); ; row=a.reader.nextLine()){
				final byte[] old=b.reader.nextLine();
				require(Arrays.equals(row, old), why);
				if(row==null){break;}
			}
		}
	}

	private static void copyText(final String input, final String output){
		try(Input in=new Input(input); Output out=new Output(output)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.writer.println(row);}
		}
	}

	private void blacklist(final Output out, final String form, final int k, final long[] words, final int size){
		require(words.length>=size, "The complete training fit needs enough observed words for its declared prefix: "+form);
		final ByteBuilder row=new ByteBuilder();
		for(int rank=0; rank<size; rank++){
			final long raw=BuildNcrnaDdlLayoutCounts.rawWord(words, rank);
			out.writer.print(row.clear().append("TRUE_WINDOW_PADDING").tab().append(form).tab().append(rank+1).tab()
				.append(raw).tab().append(raw|BuildNcrnaDdlLayoutCounts.wordTag(k, k==8, true)).tab()
				.append(words[words.length-1-rank]>>>32).tab().append(DynamicDemiLog.unpackKmer(raw, k)).tab().append(fit).nl());
		}
	}

	private static final class Sample{
		Sample(final String id_, final String input_, final String fasta_){
			id=id_;
			input=input_;
			fasta=fasta_;
		}
		final String id, input, fasta;
		long contigs, bases, windows, accepted, shortSegments, longSegments, oracleCpu;
		boolean finished;
	}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private String callerConfig, sample, binding, reuse, prefix, fit="TRAINING";
	private int first, count=-1, models=354, threads=1, targetIndex=-1;
	private boolean makeList=true;
	private NcrnaFamily family;
	private final PrintStream outstream;
	private final ArrayList<Sample> all=new ArrayList<Sample>();
	private final ArrayList<Sample> current=new ArrayList<Sample>();
	private final ArrayList<Sample> prior=new ArrayList<Sample>();
	private final ObjectIntMap<String> sampleIndex=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> samplePaths=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> currentByPath=new ObjectIntMap<String>(String.class);
	private final int[] shortCounts=new int[1<<16], longCounts=new int[1<<20];
	private final ArrayList<Counts> owners=new ArrayList<Counts>();
	private final ArrayList<Counts> available=new ArrayList<Counts>();
	private static final long LONG_MASK=AbstractKmerTable.makeMiddleMask(10, false, 1);
	public static final String FREQUENCY_HEADER="word_form\traw_word\tfrequency";
	public static final String INPUT_HEADER="assembly_id\tinput_id\tfasta\tcontigs\tbases\twindows\taccepted_windows\tshort_padding_segments\tlong_padding_segments\toracle_cpu_ns";
}
