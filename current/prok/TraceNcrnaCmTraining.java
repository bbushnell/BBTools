package prok;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import aligner.CovarianceModel;
import aligner.CovarianceModelConsensusMap;
import aligner.CovarianceModelConsensusMaps;
import aligner.CovarianceModelCyk;
import aligner.CovarianceModelParser;
import aligner.CovarianceModelPlacementBands;
import aligner.CovarianceModelSearch;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Saves unrestricted exact CM paths on source-qualified ordinary training calls.
 * Deterministic half-open call-index ranges let a bounded pilot be retained and
 * the remainder processed exactly once. Third-member calls stay on the input
 * axis but are not searched. No caller, truth fitting or method choice runs here.
 * Exact means unbanded; root selection retains the saved candidate's usual
 * half-shorter overlap association rather than learning a neighboring locus.
 * Worker/lifecycle follow template.A_SampleByteFileMT; native Streamer supplies
 * genomic batches, and every worker owns its aligner and reusable output buffers.
 * @author Brian Bushnell, Raiden
 */
public final class TraceNcrnaCmTraining implements Accumulator<TraceNcrnaCmTraining.ProcessThread>{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TraceNcrnaCmTraining.class, false);
		try{new TraceNcrnaCmTraining(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private TraceNcrnaCmTraining(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("bindings")){bindings=b;}
			else if(a.equals("refs")){refs=b.split(",", -1);}
			else if(a.equals("maps")){mapPaths=b.split(",", -1);}
			else if(a.equals("models")){cmPaths=b.split(",", -1);}
			else if(a.equals("expectedcalls")){expectedCalls=Parse.parseIntKMG(b);}
			else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
			else if(a.equals("first")){first=Parse.parseIntKMG(b);}
			else if(a.equals("limit")){limit=Parse.parseIntKMG(b);}
			else if(a.equals("padding")){padding=Parse.parseIntKMG(b);}
			else if(a.equals("maxcells")){maxCells=Parse.parseKMG(b);}
			else if(a.equals("plan")){plan=Parse.parseBoolean(b);}
			else{require(parser.parse(arg, a, b), "Unknown CM training-path option: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		require(in!=null && out!=null && bindings!=null && refs!=null && refs.length==3
			&& mapPaths!=null && mapPaths.length==2 && cmPaths!=null && cmPaths.length==2
			&& expectedCalls>0 && expectedInputs>0 && first>=0 && first<expectedCalls && limit>0
			&& (long)first+limit<=expectedCalls && padding>=0 && maxCells>0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && Shared.threads()>0,
			"Use complete training calls/bindings, three frozen refs, first/second-member maps/models, explicit first/limit and fresh out directory");
		require(!new File(out).exists(), "Path outputs must be fresh; a completed pilot is never overwritten");
		require(!Read.U_TO_T && !Read.IUPAC_TO_N && !Read.DOT_DASH_X_TO_N && !Read.LOWER_CASE_TO_N
			&& Read.JUNK_MODE==Read.CRASH_JUNK,
			"Genomic training inputs must preserve original ambiguity symbols and reject invalid DNA instead of rewriting it");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Preparation         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		readSources();
		readModels();
		readCalls();
		require(new File(out).mkdir(), "Could not create the fresh path-output directory");
		writePlan();
		if(plan){
			System.err.println("NCRNA_CM_TRAINING_PATH_PLAN_ONLY_NO_SEARCH");
			return;
		}
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Training costs require the per-thread CPU clock including exact traceback");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		try(Output summaries=new Output(out+"/calls.tsv.gz"); Output paths=new Output(out+"/paths.tsv.gz")){
			summary=summaries;
			pathOutput=paths;
			summary.writer.println(CALL_HEADER);
			pathOutput.writer.println(PATH_HEADER);
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(Shared.threads());
			for(int i=0; i<Shared.threads(); i++){workers.add(new ProcessThread(this, clock));}
			require(ThreadWaiter.startAndWait(workers, this) && success(), "A training-path worker failed; partial outputs are not complete");
			for(final Call call:calls){require(call.done, "Selected call missing from its exact source contig: "+call.index+" "+call.contig);}
		}
		require(processed==calls.size() && processed==finite+impossible+unsupported,
			"Every selected first/second-member training call needs one explicit exact-search outcome");
		try(Output output=new Output(out+"/summary.tsv")){
			output.writer.println("first_call_ordinal0\tinput_rows_in_range\texact_search_calls\tfinite\timpossible\tunsupported\tusable_anchors\ttrace_steps\texact_cpu_ns\tanchor_cpu_ns\twall_ns\tqualification");
			output.writer.print(row.clear().append(first).tab().append(limit).tab().append(processed).tab().append(finite).tab().append(impossible)
				.tab().append(unsupported).tab().append(anchored).tab().append(steps).tab().append(exactCpu).tab().append(anchorCpu)
				.tab().append(callWall).append("\tTRAINING_ONLY_FULL_TRACE_NOT_CALLER_ENDPOINT_ACCURACY\n"));
		}
		System.err.println("NCRNA_CM_TRAINING_PATHS_PASS");
	}

	private void readSources(){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(bindings)){
			RrnaResourceIO.header(input.reader, NcrnaDdlHoldout.PANEL_HEADER, bindings);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5 && p.parseInt(0)==sources.size(), "Training bindings retain dense selected-source order and five fields");
				final Source source=new Source(p.parseString(1), p.parseString(2), p.parseString(3));
				require(!byInput.contains(source.input), "One training input ID cannot denote two assemblies");
				byInput.put(source.input, source);
				sources.add(source);
			}
		}
		require(sources.size()==expectedInputs, "Training binding count differs from the excluded full-roster census");
	}

	private void readModels(){
		for(int member=0; member<3; member++){
			if(member<2){
				cms[member]=CovarianceModelParser.read(cmPaths[member]);
				require(MEMBERS[member].equals(cms[member].accession), "CM accession/order must match the frozen consensus member axis");
			}
			final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(refs[member]);
			for(final byte[] sequence:library.bases){Tools.toUpperCase(sequence);}
			final CovarianceModelConsensusMap[] maps=member<2 ?
				CovarianceModelConsensusMaps.read(cms[member], mapPaths[member], library.headers, library.bases) :
				new CovarianceModelConsensusMap[library.names.length];
			for(int m=0; m<library.names.length; m++){
				final Model model=new Model(models.size(), member, library.names[m], library.bases[m], maps[m]);
				require(!byModel.contains(model.name) && (member==2 || model.map!=null), "Each refined consensus needs its own unique saved CM map");
				byModel.put(model.name, model);
				models.add(model);
			}
		}
	}

	private void readCalls(){
		final LineParser1 p=new LineParser1('\t');
		final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
		int total=0;
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, PrepareNcrnaCmTraining.CALL_HEADER, in);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==13 && p.parseInt(0)==total, "Training call ordinals must be the intact complete census axis");
				final int index=total++;
				if(index<first || index>=first+limit){continue;}
				final int modelIndex=p.parseInt(2);
				require(modelIndex>=0 && modelIndex<models.size(), "Training model ordinal must index the frozen three-member axis");
				final Model model=models.get(modelIndex);
				require(p.termEquals(MEMBERS[model.member], 8) && NcrnaSupport.recordName(p.parseString(9)).equals(model.name),
					"Training winner name/member must agree with its exact consensus ordinal");
				if(model.member==2){continue;}
				final Source source=byInput.get(p.parseString(3));
				require(source!=null && p.termEquals(source.assembly, 1), "Call source must be in the already excluded training roster");
				final int start=p.parseInt(5), stop=p.parseInt(6);
				require(start>0 && stop>=start && (p.termEquals("+", 7) || p.termEquals("-", 7)), "Training calls require actual forward genomic intervals and explicit strands");
				// Saved grades frame seqids with BuildNcrnaReference.field. Decode
				// exactly once before joining the original native FASTA identifier.
				p.setBounds(4);
				final String contig=decoder.decode(line, p.a(), p.b());
				final Call call=new Call(index, source, model, contig, start, stop, p.termEquals("-", 7), p.parseString(10), p.parseString(12));
				calls.add(call);
				source.calls.add(call);
			}
		}
		require(total==expectedCalls && !calls.isEmpty(), "Complete census and a nonempty first/second-member range are required");
	}

	private void writePlan(){
		try(Output output=new Output(out+"/plan.tsv")){
			output.writer.println("member\tselected_calls\tstates\tmax_padded_length_upper_bound\tmax_cells_upper_bound\tmatrix_bytes_per_worker_upper_bound");
			for(int member=0; member<2; member++){
				long count=0, length=0;
				for(final Call call:calls){
					if(call.model.member==member){
						count++;
						length=Math.max(length, (long)call.stop-call.start+1+2L*padding);
					}
				}
				final long cells=Math.multiplyExact(cms[member].states(), Math.multiplyExact(length+1, length+1));
				output.writer.print(row.clear().append(MEMBERS[member]).tab().append(count).tab().append(cms[member].states())
					.tab().append(length).tab().append(cells).tab().append(Math.multiplyExact(cells, 8L)).nl());
				require(count==0 || cells<=maxCells, "The complete selected range exceeds maxcells before any exact search; member="+MEMBERS[member]+", cells="+cells);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Threads       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void accumulate(final ProcessThread worker){
		error|=!worker.success;
		processed+=worker.processed;
		finite+=worker.finite;
		impossible+=worker.impossible;
		unsupported+=worker.unsupported;
		anchored+=worker.anchored;
		steps+=worker.steps;
		exactCpu+=worker.exactCpu;
		anchorCpu+=worker.anchorCpu;
		callWall+=worker.callWall;
	}
	@Override
	public boolean success(){return !error;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	static final class ProcessThread extends Thread{
		ProcessThread(final TraceNcrnaCmTraining owner_, final ThreadMXBean clock_){
			owner=owner_;
			clock=clock_;
			for(int m=0; m<2; m++){scorers[m]=new CovarianceModelCyk(owner.cms[m], owner.maxCells);}
		}
		@Override
		public void run(){
			for(int index; (index=owner.next.getAndIncrement())<owner.sources.size();){
				final Source source=owner.sources.get(index);
				if(!source.calls.isEmpty()){processSource(source);}
			}
			success=true;
		}
		private void processSource(final Source source){
			final ObjectMap<String, Group> groups=new ObjectMap<String, Group>(String.class, Group.class);
			for(final Call call:source.calls){
				Group group=groups.get(call.contig);
				if(group==null){
					group=new Group();
					groups.put(call.contig, group);
				}
				group.calls.add(call);
			}
			final Streamer stream=StreamerFactory.makeStreamer(FileFormat.testInput(source.fasta, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
			stream.start();
			try{
				for(ListNum<Read> batch; (batch=stream.nextList())!=null;){
					for(final Read read:batch.list){
						final Group group=groups.get(RrnaResourceIO.first(read.id));
						if(group!=null){
							for(final Call call:group.calls){processCall(call, read.bases);}
						}
					}
					stream.returnList(batch);
				}
				require(!stream.errorState(), "A training genome failed streaming; exact path outputs are incomplete");
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
			// EOF was fully drained above. On a failed consumer, exceptionKill
			// prevents native producer-close from waiting forever on an unread queue.
			require(!ReadWrite.closeStream(stream), "Fully drained training genome must close without an input error");
		}

		private void processCall(final Call call, final byte[] genome){
			require(!call.done && genome!=null && call.stop<=genome.length, "Each selected call must bind exactly one real contig occurrence");
			final int windowStart=(int)Math.max(1L, (long)call.start-owner.padding);
			final int windowStop=(int)Math.min(genome.length, (long)call.stop+owner.padding);
			byte[] window=Arrays.copyOfRange(genome, windowStart-1, windowStop);
			if(call.minus){window=CovarianceModelSearch.reverseComplement(window);}
			Tools.toUpperCase(window);
			final int from=call.minus ? windowStop-call.stop+1 : call.start-windowStart+1;
			final int to=call.minus ? windowStop-call.start+1 : call.stop-windowStart+1;
			final long begin=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
			CovarianceModelCyk.Result result=null;
			String status="OK";
			try{result=scorers[call.model.member].searchOverlapping(window, window.length, from, to);}
			catch(CovarianceModelCyk.UnsupportedSequenceException expected){status="UNSUPPORTED_SEQUENCE";}
			final long cpu=clock.getCurrentThreadCpuTime()-begin;
			final long elapsed=System.nanoTime()-wall;
			require(cpu>=0 && elapsed>=0, "Exact score plus traceback needs monotonic thread CPU/wall measurements");
			processed++;
			exactCpu+=cpu;
			callWall+=elapsed;
			final long anchorBegin=clock.getCurrentThreadCpuTime();
			CovarianceModelPlacementBands placement=null;
			if(result==null){unsupported++;}
			else if(!Float.isFinite(result.score)){
				require(result.score==Float.NEGATIVE_INFINITY, "Invalid exact arithmetic must not become a biological impossible status");
				status="IMPOSSIBLE";
				impossible++;
			}else{
				finite++;
				if(call.model.sequence.length<=window.length){
					stats.clear();
					final float identity=QuantumAligner.alignAndTraceStatic(call.model.sequence, window, stats);
					require(Float.isFinite(identity), "Glocal anchor arithmetic must be finite independently of the exact answer");
					if(stats.matchString!=null && stats.rStop>=stats.rStart){placement=call.model.map.composeWindow(stats, 0);}
				}
				if(placement!=null){anchored++;}
			}
			final long glocalCpu=clock.getCurrentThreadCpuTime()-anchorBegin;
			anchorCpu+=glocalCpu;
			int departure=-1;
			pathRows.clear();
			if(result!=null && Float.isFinite(result.score)){
				departure=placement==null ? -1 : 0;
				for(int step=0; step<result.states.size; step++){
					final int state=result.states.get(step), left=result.from.get(step)-1, right=result.to.get(step);
					require(left>=0 && right>=left && right<=window.length, "Exact trace intervals must be consumed-prefix boundaries on this window");
					pathRows.append(call.index).tab().append(call.model.ordinal).tab().append(step).tab().append(state).tab().append(left).tab().append(right)
						.tab().append(result.choice.get(step)).tab().append(result.left.get(step)).tab().append(result.right.get(step));
					if(placement!=null){
						placement.stateBounds(state, bounds);
						for(final int bound:bounds){pathRows.tab().append(bound);}
						departure=Math.max(departure, Math.max(Math.max(bounds[0]-left, left-bounds[1]), Math.max(bounds[2]-right, right-bounds[3])));
					}else{pathRows.append("\t\t\t\t");}
					pathRows.nl();
				}
				steps+=result.states.size;
			}
			callRows.clear().append(call.index).tab().append(call.source.assembly).tab().append(call.model.ordinal).tab().append(call.source.input)
				.tab();
			BuildNcrnaReference.field(callRows, call.contig).tab().append(call.start).tab().append(call.stop).tab().append(call.minus ? '-' : '+')
				.tab().append(windowStart).tab().append(windowStop).tab().append(window.length).tab().append(status).tab();
			if(result!=null){callRows.appendSlow(result.score);}
			callRows.tab();
			if(result!=null && Float.isFinite(result.score)){callRows.append(result.targetFrom).tab().append(result.targetTo);}else{callRows.tab();}
			callRows.tab().append(placement!=null).tab().append(departure).tab().append(cpu).tab().append(glocalCpu).tab().append(elapsed)
				.tab().append(result==null ? 0 : result.cells).tab().append(call.outcome).tab().append(call.reference).nl();
			owner.write(callRows, pathRows);
			call.done=true;
		}
		final TraceNcrnaCmTraining owner;
		final ThreadMXBean clock;
		final CovarianceModelCyk[] scorers=new CovarianceModelCyk[2];
		final AlignmentStats stats=new AlignmentStats(true);
		final int[] bounds=new int[4];
		final ByteBuilder callRows=new ByteBuilder(), pathRows=new ByteBuilder();
		boolean success;
		long processed, finite, impossible, unsupported, anchored, steps, exactCpu, anchorCpu, callWall;
	}

	private synchronized void write(final ByteBuilder callRows, final ByteBuilder pathRows){
		assert(summary!=null && pathOutput!=null) : "The process owner must open both aggregate output streams before workers start";
		summary.writer.print(callRows);
		if(pathRows.length()>0){pathOutput.writer.print(pathRows);}
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Data and Fields       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final String assembly_, final String input_, final String fasta_){
			assembly=assembly_;
			input=input_;
			fasta=fasta_;
		}
		final String assembly, input, fasta;
		final ArrayList<Call> calls=new ArrayList<Call>();
	}
	private static final class Group{final ArrayList<Call> calls=new ArrayList<Call>();}
	private static final class Model{
		Model(final int ordinal_, final int member_, final String name_, final byte[] sequence_, final CovarianceModelConsensusMap map_){
			ordinal=ordinal_;
			member=member_;
			name=name_;
			sequence=sequence_;
			map=map_;
		}
		final int ordinal, member;
		final String name;
		final byte[] sequence;
		final CovarianceModelConsensusMap map;
	}
	private static final class Call{
		Call(final int index_, final Source source_, final Model model_, final String contig_, final int start_, final int stop_,
				final boolean minus_, final String outcome_, final String reference_){
			index=index_;
			source=source_;
			model=model_;
			contig=contig_;
			start=start_;
			stop=stop_;
			minus=minus_;
			outcome=outcome_;
			reference=reference_;
		}
		final int index, start, stop;
		final Source source;
		final Model model;
		final String contig, outcome, reference;
		final boolean minus;
		boolean done;
	}
	private String in, out, bindings;
	private String[] refs, mapPaths, cmPaths;
	private int expectedCalls, expectedInputs, first, limit, padding=20;
	private long maxCells=500000000;
	private boolean plan, error;
	private long processed, finite, impossible, unsupported, anchored, steps, exactCpu, anchorCpu, callWall;
	private Output summary, pathOutput;
	private final ByteBuilder row=new ByteBuilder();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private final ObjectMap<String, Source> byInput=new ObjectMap<String, Source>(String.class, Source.class);
	private final ObjectMap<String, Model> byModel=new ObjectMap<String, Model>(String.class, Model.class);
	private final CovarianceModel[] cms=new CovarianceModel[3];
	private final AtomicInteger next=new AtomicInteger();
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	public static final String CALL_HEADER="training_call_ordinal0\tassembly_id\tmodel_ordinal0\tinput_id\tseqid\tstart1\tstop1\tstrand\twindow_start1\twindow_stop1\twindow_length\tstatus\texact_bits\texact_from1\texact_to1\tusable_glocal_anchor\tmax_whole_path_departure\texact_trace_cpu_ns\tglocal_cpu_ns\tsearch_wall_ns\tmatrix_cells\tsaved_query_outcome\tmatched_reference_id";
	public static final String PATH_HEADER="training_call_ordinal0\tmodel_ordinal0\tstep\tstate\tstart_boundary0\tend_boundary0\tchoice\tleft_child_step\tright_child_step\tanchor_start_low0\tanchor_start_high0\tanchor_end_low0\tanchor_end_high0";
}
