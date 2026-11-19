package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import map.LongObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Exact selected all-model verification on every saved RF00174 window.
 * Adapted from A_SampleBasic2/A_SampleByteFileMT: byte metadata prepasses,
 * native Streamer sequence owners, coarse original-file work claims and an
 * owner-only output pass. The actual parsed caller supplies resources/controls;
 * no scheduling, claimed-region updates, competing-family calling or final path
 * runs. Scalar identity pairs are60/62/64%; every64% proposal and work count must
 * reproduce its saved window baseline before any completion report is written.
 * @author Brian Bushnell, Raiden
 */
public final class ReplayNcrnaWindows{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReplayNcrnaWindows.class, false);
		final Timer timer=new Timer();
		try{
			final ReplayNcrnaWindows tool=new ReplayNcrnaWindows(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println((tool.shuffleWindows>0 ? "NCRNA_SHUFFLED_WINDOW_CONTROL_COMPLETE" : "NCRNA_SAVED_WINDOW_ACCEPTANCE_REPLAY_COMPLETE")+" time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private ReplayNcrnaWindows(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){callerConfig=b;}
			else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
			else if(a.equals("expectedwindows")){expectedWindows=Parse.parseIntKMG(b);}
			else if(a.equals("expectedscheduled")){expectedScheduled=Parse.parseKMG(b);}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(a.equals("shufflewindows")){shuffleWindows=Parse.parseIntKMG(b);}
			else if(a.equals("shuffleseed")){shuffleSeed=Long.parseLong(b);}
			else if(a.equals("shufflecache")){shuffleCache=Parse.parseBoolean(b);}
			else{require(parser.parse(arg, a, b), "Unknown saved-acceptance replay parameter: "+arg);}
		}
		capture=parser.in1;
		prefix=parser.out1;
		threads=Shared.threads();
		require(capture!=null && prefix!=null && callerConfig!=null && expectedInputs>0
			&& expectedModels>0 && expectedWindows>=0 && expectedScheduled>=expectedWindows
			&& threads>0 && threads<=32 && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=capture_prefix caller=exact_flat_baseline_config out=fresh_prefix and exact inputs/windows/scheduled/models counts; no sampling");
		require(shuffleWindows>=0 && shuffleWindows<=expectedWindows, "The declared shuffled sample must fit the complete captured-window universe");
		if(shuffleWindows>0){RrnaResourceIO.fresh(prefix, ".shuffle.windows.tsv", ".shuffle.summary.tsv");}
		else{for(final String cut:CUT_NAMES){RrnaResourceIO.fresh(prefix+"."+cut, ".results.tsv", ".inputs.tsv", ".replay.tsv");}}
	}

	private void process(){
		final String[] baseline=Parser.parseConfig(new String[]{"config="+
			NcrnaCallerProcess.configurationPath(new File(callerConfig))});
		final ArrayList<String> flags=new ArrayList<String>();
		boolean taxonomyOff=false;
		for(final String flag:baseline){
			final String[] split=Parse.splitOnFirst(flag, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("taxonomy") || a.equals("tax")){taxonomyOff=!Parse.parseBoolean(b);}
			flags.add(flag);
		}
		require(taxonomyOff, "Resource-only replay initialization requires the original final taxonomy=f setting");
		flags.add("out=null");
		flags.add("outstats=null");
		flags.add("t="+threads);
		try(CallGenes configured=new CallGenes(flags.toArray(new String[flags.size()]))){
			int targetModels=-1;
			for(final NcrnaFamily family:GeneCaller.ncrnaFamilies){
				if(!"cobalamin_rf00174".equals(family.name)){continue;}
				require(targetModels<0 && family.idPass==.64f && family.idBorderline==.64f
					&& family.minLen==109 && family.seedMinHits==2 && family.seedDistinct,
					"The replay must start from the actual ordinary64%/109-base/two-distinct-key baseline");
				targetModels=family.library.length;
			}
			require(targetModels==expectedModels, "The exact configured reference must retain every expected original library model");
			contextFlank=configured.savedCobalaminWindowWorker(.64f).savedWindowFlank();
			loadInputs();
			final String[] configuredInputs=configured.inputFiles();
			require(configuredInputs.length==inputs.size(), "The exact caller configuration must retain every original captured input and no extras");
			for(int i=0; i<configuredInputs.length; i++){
				require(configuredInputs[i].equals(inputs.get(i).fasta), "The exact caller configuration must preserve original captured FASTA paths and order");
			}
			loadWindows();
			loadSequences(contextFlank);
			loadOptionalContexts();
			for(final OriginalInput input:inputs){input.contexts.prepare();}
			if(shuffleWindows>0){
				NcrnaShuffledWindowControl.run(configured, inputs, prefix, shuffleWindows, shuffleSeed, shuffleCache, threads);
				return;
			}
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
			for(int i=0; i<Math.min(threads, inputs.size()); i++){workers.add(new ProcessThread(this, configured));}
			ThreadWaiter.startAndWait(workers);
			if(failure.get()!=null){throw new IllegalStateException("Saved acceptance replay failed; no completion is valid", failure.get());}
			for(final ProcessThread worker:workers){require(worker.success, "All original-file workers must finish before publishing replay observations");}
			writeResults();
			writeInputsAndSummaries();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Original Metadata       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadInputs(){
		try(Input in=new Input(capture+".inputs.tsv")){
			RrnaResourceIO.header(in.reader, CAPTURE_INPUT_HEADER, capture);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==10 && p.parseInt(0)==inputs.size(), "Original capture slots must remain contiguous and complete");
				final OriginalInput input=new OriginalInput(inputs.size(), decoded(1), decoded(2),
					p.parseByteArray(1), p.parseByteArray(2), p.parseLong(3), p.parseLong(4), p.parseInt(5), p.parseLong(6));
				require(input.contigs>0 && input.bases>0 && input.expectedWindows>=0
					&& input.expectedScheduled>=input.expectedWindows && new File(input.fasta).isAbsolute(),
					"Keep every original input, including zero-window files, with exact path and recorded dimensions");
				inputs.add(input);
			}
		}
		require(inputs.size()==expectedInputs, "A replay cannot omit zero-call or background source inputs");
	}

	private void loadWindows(){
		try(Input in=new Input(capture+".windows.tsv")){
			RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.WINDOW_HEADER, capture);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==15, "Every saved candidate needs its original physical contig and oriented bounds");
				final OriginalInput input=input(1, 2, 3);
				final long key=eventKey(p, 0);
				require((int)(key>>>32)==input.slot && !windows.contains(key), "Captured event identities must remain unique in their original input namespace");
				final Window window=new Window(key, p.parseString(0), input, p.parseLong(4),
					NcrnaDdlOracle.seqid(decoded(5)), p.parseInt(6), p.parseInt(7), p.parseInt(8),
					p.parseInt(9), p.parseInt(10), p.parseInt(13));
				require(window.ordinal>=0 && window.ordinal<input.contigs && (window.strand==0 || window.strand==1)
					&& window.length>0 && (window.pass==1 || window.pass==2) && window.start>=0
					&& window.stop>=window.start && window.stop<window.length && p.parseInt(14)==window.stop-window.start+1
					&& p.parseInt(11)==window.forwardStart() && p.parseInt(12)==window.forwardStop(),
					"Saved candidate orientation, size and original-forward coordinates must agree exactly");
				windows.put(key, window);
				input.windows.add(window);
			}
		}
		require(windows.size()==expectedWindows, "Replay all captured windows, including those with no annotated gene");
		for(final OriginalInput input:inputs){require(input.windows.size()==input.expectedWindows,
			"Every original input must retain all of its saved windows");}
	}

	private void loadSequences(final int flank){
		Shared.TRIM_READ_DESCRIPTION=false;
		NcrnaDdlSavedData.sequences(capture+".windows.fa.gz", read->{
			final Window window=windows.get(eventKey(read.id));
			require(window!=null && window.event.equals(read.id) && window.context==null && read.length()==window.stop-window.start+1,
				"One saved FASTA sequence must supply each declared candidate exactly once");
			window.context=window.input.contexts.add(window.event, window.ordinal, window.name, window.length,
				window.strand, window.start, window.stop, flank, read.bases);
			for(final byte base:read.bases){
				final byte upper=base>='a' && base<='z' ? (byte)(base-32) : base;
				if(upper=='G' || upper=='C'){
					window.gc++;
					window.acgt++;
				}
				else if(upper=='A' || upper=='T'){window.acgt++;}
			}
		});
		for(final OriginalInput input:inputs){for(final Window window:input.windows){require(window.context!=null, "No captured event may lose its saved candidate bytes");}}
	}

	/** New captures may already contain complete refinement flanks. Metadata and
	 * FASTA must both exist, bind every original event, and agree with candidates. */
	private void loadOptionalContexts(){
		final boolean metadata=new File(capture+".context.tsv").exists(), sequences=new File(capture+".context.fa.gz").exists();
		require(metadata==sequences, "A context sidecar needs both its exact metadata and real FASTA bytes");
		if(!metadata){return;}
		try(Input in=new Input(capture+".context.tsv")){
			RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.CONTEXT_HEADER, capture);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==4, "Context metadata must retain event, original-oriented interval and required flank");
				final Window window=windows.get(eventKey(p, 0));
				require(window!=null && !window.sidecarMetadata && p.parseInt(1)==window.context.origin
					&& p.parseInt(2)==window.context.origin+window.context.bases.length-1 && p.parseInt(3)==contextFlank,
					"New capture context must bind the replay's exact original event and complete stored interval");
				window.sidecarMetadata=true;
			}
		}
		NcrnaDdlSavedData.sequences(capture+".context.fa.gz", read->{
			final Window window=windows.get(eventKey(read.id));
			require(window!=null && window.event.equals(read.id) && window.sidecarMetadata && !window.sidecarSequence,
				"Each real context FASTA must bind one declared event exactly once");
			window.input.contexts.supply(window.context, window.context.origin, read.bases);
			window.sidecarSequence=true;
		});
		for(final OriginalInput input:inputs){for(final Window window:input.windows){require(window.sidecarMetadata && window.sidecarSequence,
			"An optional context sidecar cannot silently omit any original candidate");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        File Workers          ----------------*/
	/*--------------------------------------------------------------*/

	private static final class ProcessThread extends Thread{
		ProcessThread(final ReplayNcrnaWindows owner_, final CallGenes configured){
			owner=owner_;
			for(int i=0; i<CUTS.length; i++){workers[i]=configured.savedCobalaminWindowWorker(CUTS[i]);}
		}
		@Override
		public void run(){
			try{
				for(int i=owner.next.getAndIncrement(); i<owner.inputs.size(); i=owner.next.getAndIncrement()){
					if(owner.failure.get()!=null){return;}
					measure(owner.inputs.get(i));
				}
				success=true;
			}catch(RuntimeException | Error invalid){owner.failure.compareAndSet(null, invalid);}
		}
		private void measure(final OriginalInput input){
			assert(input.contexts.windows()==input.windows.size()) : "One worker must own every retained window in its original source file";
			for(final Window window:input.windows){
				for(int cut=0; cut<CUTS.length; cut++){
					final NcrnaScavenger worker=workers[cut];
					final long before=worker.alignmentCount();
					window.accepted[cut]=worker.verifySavedWindow(window.name, window.context.bases, window.context.origin,
						window.length, window.strand, window.pass, window.start, window.stop, window.source,
						input.contexts.sourceFor(window.context));
					window.work[cut]=worker.alignmentCount()-before;
					require(window.work[cut]>=0, "Actual verifier work cannot run backward between saved windows");
				}
			}
		}
		private final ReplayNcrnaWindows owner;
		private final NcrnaScavenger[] workers=new NcrnaScavenger[CUTS.length];
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Original Row Order     ----------------*/
	/*--------------------------------------------------------------*/

	private void writeResults(){
		try(Input in=new Input(capture+".results.tsv"); Output low=new Output(prefix+".060.results.tsv");
			Output middle=new Output(prefix+".062.results.tsv"); Output control=new Output(prefix+".064.results.tsv")){
			final Output[] outputs={low, middle, control};
			for(final Output out:outputs){out.writer.println(RESULT_HEADER);}
			RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.RESULT_HEADER, capture);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==26, "Use the original complete result table, retaining uncaptured zero-work rejections");
				final OriginalInput input=input(0, 1, 2);
				input.scheduled++;
				final Window window=p.termEquals("NA", 12) ? null : windows.get(eventKey(p, 12));
				if(window!=null){
					require(!window.written && window.input==input && p.parseLong(3)==window.ordinal
						&& NcrnaDdlOracle.seqid(decoded(4)).equals(window.name) && p.parseInt(5)==window.strand
						&& p.parseInt(6)==window.length && p.parseInt(7)==window.pass
						&& p.parseInt(8)==window.start && p.parseInt(9)==window.stop
						&& p.parseInt(10)==window.forwardStart() && p.parseInt(11)==window.forwardStop(),
						"Every original captured result must bind its one saved candidate and both coordinate systems");
					checkControl(window);
					window.written=true;
				}else{require(p.termEquals("NA", 12) && p.termEquals("UNCAPTURED_GATE_REJECTION", 13) && p.parseLong(17)==0,
					"A missing event is allowed only for a genuinely uncaptured original zero-work gate rejection");}
				p.setBounds(13);
				final int beforeStatus=p.a();
				for(int cut=0; cut<CUTS.length; cut++){
					final NcrnaDdlOrf accepted=window==null ? null : window.accepted[cut];
					final long work=window==null ? 0 : window.work[cut];
					input.alignments[cut]+=work;
					if(accepted!=null){input.accepts[cut]++;}
					row.clear().append(line, 0, beforeStatus)
						.append(window==null ? "UNCAPTURED_GATE_REJECTION" : accepted==null ? "REJECTED_VERIFICATION" : "ACCEPTED_PROPOSAL")
						.tab().append(work).tab();
					if(accepted==null){row.append("NA\tNA\tNA\tNA\tNA\tNA");}
					else{
						BuildNcrnaReference.field(row, accepted.trnaModel).tab().append(accepted.acceptedIdentity(), 9)
							.tab().append(accepted.proposalRank()).tab().append(accepted.proposalScore(), 3)
							.tab().append(window.acceptedStart(accepted)).tab().append(window.acceptedStop(accepted));
					}
					outputs[cut].writer.print(row.tab().append(window==null ? -1 : window.gc)
						.tab().append(window==null ? -1 : window.acgt).nl());
				}
			}
		}
		long scheduled=0;
		for(final OriginalInput input:inputs){
			require(input.scheduled==input.expectedScheduled, "All original scheduled rows, including uncaptured background, must survive replay");
			scheduled+=input.scheduled;
			for(final Window window:input.windows){require(window.written, "A saved candidate without its completed original result cannot be omitted");}
		}
		require(scheduled==expectedScheduled && controlWindows==expectedWindows, "The control must verify every captured window and the entire scheduled universe");
	}

	/** Stronger than aggregate gene parity: the unchanged control reproduces
	 * each proposal's winner, refined span, identity, order/score and total work. */
	private void checkControl(final Window window){
		final NcrnaDdlOrf actual=window.accepted[2];
		final boolean expected=p.termEquals("ACCEPTED_PROPOSAL", 13);
		require(expected || p.termEquals("REJECTED_VERIFICATION", 13), "Captured baseline rows must declare an actual acceptance or verifier rejection");
		if(expected!=(actual!=null) || window.work[2]!=p.parseLong(17)){
			throw new IllegalStateException("The64% control changed acceptance/work at input "+window.input.name+", event "+window.event);
		}
		if(actual!=null){
			BuildNcrnaReference.field(text.clear(), actual.trnaModel);
			if(!same(p, 18, text.array, text.length()) || p.parseFloat(19)!=actual.acceptedIdentity()
				|| p.parseInt(20)!=actual.proposalRank()
				|| p.parseInt(22)!=window.acceptedStart(actual) || p.parseInt(23)!=window.acceptedStop(actual)){
				throw new IllegalStateException("The64% control changed the actual accepted proposal at input "+window.input.name+", event "+window.event);
			}
			// CaptureNcrnaWindows retains three decimal places for DDL score.
			// Compare that exact observation, not a rounded text float against
			// the unrounded in-memory score or an arbitrary tolerance.
			text.clear().append(actual.proposalScore(), 3);
			require(same(p, 21, text.array, text.length()), "The64% control changed the recorded three-decimal DDL score");
		}
		require(p.parseInt(24)==window.gc && p.parseInt(25)==window.acgt, "Loaded candidate bytes must reproduce the original known-base GC observation");
		controlWindows++;
	}

	private void writeInputsAndSummaries(){
		long sourceLoads=0, sourceBases=0;
		for(final OriginalInput input:inputs){
			sourceLoads+=input.contexts.sourceLoads();
			sourceBases+=input.contexts.sourceBases();
		}
		for(int cut=0; cut<CUTS.length; cut++){
			long alignments=0, accepts=0;
			try(Output out=new Output(prefix+"."+CUT_NAMES[cut]+".inputs.tsv")){
				out.writer.println(INPUT_HEADER);
				for(final OriginalInput input:inputs){
					out.writer.print(row.clear().append(input.slot).tab().append(input.nameField).tab().append(input.pathField)
						.tab().append(input.contigs).tab().append(input.bases).tab().append(input.windows.size())
						.tab().append(input.scheduled).tab().append(input.alignments[cut]).tab().append(input.accepts[cut]).nl());
					alignments+=input.alignments[cut];
					accepts+=input.accepts[cut];
				}
			}
			try(Output out=new Output(prefix+"."+CUT_NAMES[cut]+".replay.tsv")){
				out.writer.print(row.clear().append("metric\tvalue\nidentity\t").append(Float.toString(CUTS[cut]))
					.append("\ninputs\t").append(inputs.size()).append("\ncaptured_windows\t").append(expectedWindows)
					.append("\nscheduled_windows\t").append(expectedScheduled).append("\nactual_alignments\t").append(alignments)
					.append("\naccepted_proposals\t").append(accepts)
					.append("\nshared_source_file_loads_NONADDITIVE_ACROSS_CUTS\t").append(sourceLoads)
					.append("\nshared_source_bases_read_NONADDITIVE_ACROSS_CUTS\t").append(sourceBases)
					.append("\ncontrol_windows_checked\t").append(cut==2 ? controlWindows : 0)
					.append("\nmode\tREAL_WINDOW_VERIFICATION_NO_FINAL_PATH\n"));
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	static final class OriginalInput{
		OriginalInput(final int slot_, final String name_, final String fasta_, final byte[] nameField_, final byte[] pathField_,
				final long contigs_, final long bases_, final int windows_, final long scheduled_){
			slot=slot_;
			name=name_;
			fasta=fasta_;
			nameField=nameField_;
			pathField=pathField_;
			contigs=contigs_;
			bases=bases_;
			expectedWindows=windows_;
			expectedScheduled=scheduled_;
			contexts=new NcrnaSavedWindowContexts(name, fasta);
		}
		final int slot, expectedWindows;
		final long contigs, bases, expectedScheduled;
		final String name, fasta;
		final byte[] nameField, pathField;
		final NcrnaSavedWindowContexts contexts;
		final ArrayList<Window> windows=new ArrayList<Window>();
		final long[] alignments=new long[3], accepts=new long[3];
		long scheduled;
	}
	static final class Window{
		Window(final long key_, final String event_, final OriginalInput input_, final long ordinal_, final String name_,
				final int strand_, final int length_, final int pass_, final int start_, final int stop_, final int source_){
			key=key_;
			event=event_;
			input=input_;
			ordinal=ordinal_;
			name=name_;
			strand=strand_;
			length=length_;
			pass=pass_;
			start=start_;
			stop=stop_;
			source=source_;
		}
		int forwardStart(){return strand==0 ? start : length-1-stop;}
		int forwardStop(){return strand==0 ? stop : length-1-start;}
		int acceptedStart(final NcrnaDdlOrf result){return strand==0 ? context.origin+result.start : length-1-context.origin-result.stop;}
		int acceptedStop(final NcrnaDdlOrf result){return strand==0 ? context.origin+result.stop : length-1-context.origin-result.start;}
		final long key, ordinal;
		final String event, name;
		final OriginalInput input;
		final int strand, length, pass, start, stop, source;
		final NcrnaDdlOrf[] accepted=new NcrnaDdlOrf[3];
		final long[] work=new long[3];
		NcrnaSavedWindowContexts.Window context;
		int gc, acgt;
		boolean sidecarMetadata, sidecarSequence, written;
	}

	private OriginalInput input(final int slot, final int name, final int path){
		final int index=p.parseInt(slot);
		require(index>=0 && index<inputs.size(), "Every saved row must bind a declared original input slot");
		final OriginalInput input=inputs.get(index);
		require(same(p, name, input.nameField, input.nameField.length) && same(p, path, input.pathField, input.pathField.length),
			"A captured or completed event cannot substitute another original input identity/path");
		return input;
	}
	private String decoded(final int column){p.setBounds(column);return decoder.decode(p.line(), p.a(), p.b());}
	private static boolean same(final LineParser1 p, final int column, final byte[] bytes, final int length){
		p.setBounds(column);
		if(p.b()-p.a()!=length){return false;}
		for(int i=0; i<length; i++){if(p.line()[p.a()+i]!=bytes[i]){return false;}}
		return true;
	}
	private static long eventKey(final LineParser1 p, final int column){
		p.setBounds(column);
		final byte[] line=p.line();
		final int from=p.a(), to=p.b();
		int separator=from+1;
		while(separator<to && line[separator]>='0' && line[separator]<='9'){separator++;}
		require(to>from+3 && line[from]=='i' && separator>from+1 && separator<to-1 && line[separator]=='e',
			"Saved event IDs must retain native iINPUTeEVENT numeric identities");
		final int slot=Parse.parseInt(line, from+1, separator), event=Parse.parseInt(line, separator+1, to);
		require(slot>=0 && event>=0, "Saved input/event ordinals cannot be negative or overflow their original integer namespace");
		return ((long)slot<<32)|event;
	}
	/** Native FASTA IDs are already retained Strings; parse their integer pieces
	 * directly rather than allocating field substrings or byte-array copies. */
	private static long eventKey(final String name){
		require(name!=null && name.length()>=4 && name.charAt(0)=='i', "Saved FASTA IDs must retain native iINPUTeEVENT framing");
		int index=1;
		long slot=0, event=0;
		for(; index<name.length() && name.charAt(index)>='0' && name.charAt(index)<='9'; index++){
			slot=slot*10+name.charAt(index)-'0';
			require(slot<=Integer.MAX_VALUE, "Original saved input slots must fit their integer namespace");
		}
		require(index>1 && index<name.length()-1 && name.charAt(index)=='e', "A saved FASTA ID must retain both original integer ordinals");
		for(index++; index<name.length(); index++){
			final char digit=name.charAt(index);
			require(digit>='0' && digit<='9', "Saved FASTA event ordinals contain decimal digits only");
			event=event*10+digit-'0';
			require(event<=Integer.MAX_VALUE, "Original saved event ordinals must fit their integer namespace");
		}
		return (slot<<32)|event;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private String capture, prefix, callerConfig;
	private int expectedInputs, expectedWindows=-1, expectedModels, threads, contextFlank;
	private long expectedScheduled=-1, controlWindows;
	private int shuffleWindows;
	private long shuffleSeed=1701;
	private boolean shuffleCache=true;
	private final ArrayList<OriginalInput> inputs=new ArrayList<OriginalInput>();
	private final LongObjectMap<Window> windows=new LongObjectMap<Window>(Window.class);
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder(), text=new ByteBuilder();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final float[] CUTS={.60f,.62f,.64f};
	private static final String[] CUT_NAMES={"060","062","064"};
	public static final String INPUT_HEADER="input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tscheduled_windows\talignments\taccepted_proposals";
	public static final String RESULT_HEADER="input_slot\tinput_id\tfasta\tcontig_ordinal\tcontig\tstrand\tcontig_length\tpass\toriented_start0\toriented_stop0\tforward_start0\tforward_stop0\tcaptured_event_id\tstatus\ttotal_alignments\tmodel\taccepted_identity\tddl_rank\tddl_score\taccepted_forward_start0\taccepted_forward_stop0\tgc_bases\tacgt_bases";
	private static final String CAPTURE_INPUT_HEADER=INPUT_HEADER+"\tfinal_target_call_rows";
}
