package prok;

import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** A_SampleByteFileMT lifecycle with retained candidate events and one FASTA
 * stream per source. Additional CMs use their own library/map and a fresh
 * score-only worker; no endpoint, caller path or reference assignment changes.
 * The native own-CM measurement is reused, not recomputed for reproducibility.
 * @author Brian Bushnell, Raiden
 */
public final class ScoreCobalaminClanCandidates implements Accumulator<ScoreCobalaminClanCandidates.ProcessThread>{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ScoreCobalaminClanCandidates.class, false);
		try{new ScoreCobalaminClanCandidates(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private ScoreCobalaminClanCandidates(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("familyroot")){familyRoot=b;}
			else if(a.equals("scoreroot")){scoreRoot=b;}
			else if(a.equals("profiles")){profiles=b.split(",", -1);}
			else if(a.equals("expectedinputs")){expectedInputs=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("lowidentity")){low=(float)NcrnaNumbers.doubleValue(a, b);}
			else if(a.equals("highidentity")){high=(float)NcrnaNumbers.doubleValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown saved cross-CM scoring argument: "+arg);}
		}
		in=parser.in1;
		inputs=parser.in2;
		out=parser.out1;
		require(in!=null && inputs!=null && out!=null && parser.out2==null && familyRoot!=null && scoreRoot!=null
			&& profiles!=null && profiles.length==3 && expectedInputs>0 && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply in=capture.acceptance.tsv in2=capture.inputs.tsv familyroot=selected_group scoreroot=score_group profiles=main,second,third expectedinputs=N out=fresh_prefix");
		require(Float.isFinite(low) && Float.isFinite(high) && low>=0 && high>low && high<=1,
			"The explicit identity band must be finite and increasing on the caller's float axis");
		require(NcrnaCpuProfile.enabled() && NcrnaCpuProfile.cpuNow()>=0,
			"Cross-CM cost attribution requires the native worker CPU clock before any file is opened");
		require(Tools.testInputFiles(false, true, in, inputs) && Tools.testInputFiles(false, true, profiles),
			"Saved candidates, their exact ordinal/path roster and all three profile resources must be readable");
		RrnaResourceIO.fresh(out, ".scores.tsv", ".costs.tsv");
	}

	private void load(){
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		for(int member=0; member<3; member++){
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			options.resources.root=familyRoot;
			options.endpoints=false;
			options.boundary=false;
			options.boundaryOverride=false;
			final NcrnaFamily family=NcrnaFamilyLoader.load(FAMILIES[member], options, resolver);
			families[member]=family;
			final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
			overrides.root=member==0 ? familyRoot : scoreRoot;
			factories[member]=CobalaminCmEndpointRefiner.load(family, resolver, overrides,
				member==0 ? CobalaminCmEndpointRefiner.LEARNED : CobalaminCmEndpointRefiner.RETRY_64, false);
			final ObjectIntMap<String> names=new ObjectIntMap<String>(String.class);
			for(int model=0; model<family.modelNames.length; model++){
				final String name=first(family.modelNames[model]);
				require(!names.contains(name), "Saved model names must resolve to one target-library ordinal");
				names.put(name, model);
			}
			modelIds.add(names);
		}
		final String[][] names=new String[3][];
		final byte[][][] sequences=new byte[3][][];
		for(int member=0; member<3; member++){
			names[member]=families[member].modelNames;
			sequences[member]=families[member].library;
		}
		crossScorer=new CobalaminCrossMemberCmScorer(families, factories,
			NcrnaDdlProfileResource.loadClan(profiles, names, sequences));
		for(int ordinal=0; ordinal<expectedInputs; ordinal++){sources.add(null);}
		int observed=0;
		try(Input input=new Input(inputs)){
			final LineParser1 p=new LineParser1('\t');
			final byte[] header=input.reader.nextLine();
			require(header!=null, "The capture must retain its actual input-path roster");
			p.set(header);
			require(p.terms()==9 && p.termEquals("input_ordinal", 0) && p.termEquals("input_path", 1),
				"Use the native timing driver's exact ordinal/path roster, never a reordered assembly list");
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==9, "Capture input rows must retain their complete native timing schema");
				final int ordinal=p.parseInt(0);
				require(ordinal>=0 && ordinal<expectedInputs && sources.get(ordinal)==null,
					"Native timing rows arrive in completion order; each declared input ordinal must occur exactly once");
				sources.set(ordinal, new Source(p.parseString(1)));
				observed++;
			}
		}
		require(observed==expectedInputs, "The full declared input population is bound before reading events");
		for(final Source source:sources){require(source!=null, "A missing input ordinal cannot be replaced by its neighbor's completion row");}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	private void readEvents(){
		assert(sources.size()==expectedInputs && eventsRead==0) : "Candidate ordinals are meaningful only after loading the complete source roster once";
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, CobalaminCmAcceptance.HEADER, in);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==17 && p.termEquals("true", 12), "The input is an unfiltered capture with complete accepted event rows");
				final long index=eventsRead++;
				final int source=p.parseInt(0), member=member(p.parseString(3));
				final float identity=Float.intBitsToFloat(p.parseInt(8));
				require(source>=0 && source<sources.size() && Float.isFinite(identity) && identity>=low && identity<=1,
					"Every capture event must bind a real input and an identity at or above its declared capture floor");
				if(identity>=high){continue;}
				final String model=p.parseString(4);
				final int modelId=modelIds.get(member).get(model);
				final int length=p.parseInt(5), start=p.parseInt(6), stop=p.parseInt(7);
				require(modelId>=0 && length>0 && start>=0 && stop>=start && stop<length
					&& (p.termEquals("+", 2) || p.termEquals("-", 2)),
					"A band event must preserve its actual model, strand and zero-based inclusive contig coordinates");
				final float score=p.termEquals("", 9) ? Float.NaN : Float.intBitsToFloat(p.parseInt(9));
				final float ga=Float.intBitsToFloat(p.parseInt(10));
				require(Float.isFinite(ga) && validScore(score), "Saved own-CM scores retain finite, impossible and absent states separately");
				final Event event=new Event(index, line, p.parseString(1), member, modelId, length, start, stop,
					p.termEquals("-", 2), score, ga);
				sources.get(source).events.add(event);
				bandEvents++;
			}
		}
		require(bandEvents>0, "An empty identity band cannot establish cross-CM discrimination or cost");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		load();
		readEvents();
		try(Output scores=new Output(out+".scores.tsv")){
			output=scores;
			scores.writer.println(SCORE_HEADER);
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
			for(int i=0; i<Math.min(Shared.threads(), sources.size()); i++){workers.add(new ProcessThread(this));}
			final boolean joined=ThreadWaiter.startAndWait(workers, this);
			if(error!=null){throw new IllegalStateException("Saved candidate scoring failed on a worker; partial output is invalid", error);}
			require(joined && success, "Every saved-source scoring worker must finish; partial score output is not a completed measurement");
		}
		require(completed==bandEvents, "Each retained band event must be encountered exactly once on its declared contig");
		try(Output costs=new Output(out+".costs.tsv")){
			costs.writer.println("emitted_member\tscoring_member\tevents\treused\tunavailable\textra_dp_fills\textra_ranking_cpu_ns\textra_anchor_cpu_ns\textra_cm_cpu_ns\textra_total_cpu_ns");
			final ByteBuilder row=new ByteBuilder();
			for(int emitted=0; emitted<3; emitted++){
				for(int target=0; target<3; target++){
					row.clear().append(MEMBERS[emitted]).tab().append(MEMBERS[target]);
					for(int field=0; field<COST_FIELDS; field++){row.tab().append(totals[(emitted*3+target)*COST_FIELDS+field]);}
					costs.writer.print(row.nl());
				}
			}
		}
		System.err.println("COBALAMIN_CLAN_CANDIDATES_SCORED events="+completed+" capture_events="+eventsRead);
	}

	@Override
	public void accumulate(final ProcessThread worker){
		require(worker!=null, "Joined results belong to a real scoring worker");
		success&=worker.success;
		if(error==null){error=worker.error;}
		completed+=worker.completed;
		for(int i=0; i<totals.length; i++){totals[i]+=worker.costs[i];}
	}
	@Override
	public ReadWriteLock rwlock(){return null;}
	@Override
	public boolean success(){return success;}

	/*--------------------------------------------------------------*/
	/*----------------        Worker Threads        ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread{
		ProcessThread(final ScoreCobalaminClanCandidates owner_){
			owner=owner_;
			cross=owner.crossScorer.newWorker();
		}
		@Override
		public void run(){
			try{
				for(int index; !owner.stop && (index=owner.next.getAndIncrement())<owner.sources.size();){
					final Source source=owner.sources.get(index);
					if(!source.events.isEmpty()){processSource(source);}
				}
				success=true;
			}catch(Throwable failure){
				error=failure;
				owner.stop=true;
			}
		}
		private void processSource(final Source source){
			assert(source!=null && !source.events.isEmpty()) : "Only sources with saved band events need a genome stream";
			final ObjectMap<String, Group> groups=new ObjectMap<String, Group>(String.class, Group.class);
			for(final Event event:source.events){
				Group group=groups.get(event.contig);
				if(group==null){
					group=new Group();
					groups.put(event.contig, group);
				}
				group.events.add(event);
			}
			try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(
				FileFormat.testInput(source.fasta, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1))){
				input.stream.start();
				for(ListNum<Read> batch; (batch=input.stream.nextList())!=null;){
					for(final Read read:batch.list){
						final Group group=groups.get(first(read.id));
						if(group==null){continue;}
						final Orf candidate=new Orf(read.id, 0, 0, Shared.PLUS, 0, read.bases, false, Orf.RNA);
						for(int strand=0; strand<2; strand++){
							if(strand==1){AminoAcid.reverseComplementBasesInPlace(read.bases);}
							for(final Event event:group.events){if(event.minus==(strand==1)){processEvent(event, read.bases, candidate);}}
						}
						// The streamed read owns these bases; no other consumer sees the reversed array.
					}
					input.stream.returnList(batch);
				}
			}
			for(final Event event:source.events){require(event.done, "A band event's contig was absent from its exact capture input: "+event.contig);}
		}
		private void processEvent(final Event event, final byte[] bases, final Orf candidate){
			require(!event.done && bases.length==event.length, "Each candidate must bind exactly one unchanged source contig");
			final int start=event.minus ? bases.length-1-event.stop : event.start;
			final int stop=event.minus ? bases.length-1-event.start : event.stop;
			candidate.start=start;
			candidate.stop=stop;
			cross.prepare(candidate, bases);
			try{
				for(int target=0; target<3; target++){
					final int slot=(event.member*3+target)*COST_FIELDS;
					final boolean own=target==event.member;
					int model=event.model;
					long ranking=0, anchor=0, cm=0, fills=0;
					float score=event.score, ga=event.ga;
					String status=Float.isNaN(score) ? "SAVED_UNAVAILABLE" : "REUSED_OWN_CM";
					if(!own){
						score=cross.scoreChosen(target, candidate, bases);
						model=cross.model();
						ga=cross.ga();
						ranking=cross.rankingCpu();
						anchor=cross.anchorCpu();
						cm=cross.cmCpu();
						fills=cross.fills();
						status=cross.status();
					}
					final long total=own ? 0 : cross.totalCpu();
					require(candidate.start==start && candidate.stop==stop && validScore(score) && Float.isFinite(ga)
						&& ranking>=0 && anchor>=0 && cm>=0 && fills>=0 && total>=0,
						"Score-only work must preserve ends and finite monotonic CPU accounting; absent scores stay explicit");
					costs[slot]++;
					if(own){costs[slot+1]++;}
					if(Float.isNaN(score)){costs[slot+2]++;}
					costs[slot+3]+=fills;
					costs[slot+4]+=ranking;
					costs[slot+5]+=anchor;
					costs[slot+6]+=cm;
					costs[slot+7]+=total;
					row.clear().append(event.index).tab().append(event.line).tab().append(MEMBERS[target]).tab()
						.append(first(owner.families[target].modelNames[model])).tab();
					if(!Float.isNaN(score)){row.append(Float.floatToIntBits(score));}
					row.tab().append(Float.floatToIntBits(ga)).tab().append(status).tab().append(fills)
						.tab().append(ranking).tab().append(anchor).tab().append(cm).tab().append(total).nl();
					synchronized(owner.output){owner.output.writer.print(row);}
				}
			}finally{cross.clearPrepared();}
			event.done=true;
			completed++;
		}
		private final ScoreCobalaminClanCandidates owner;
		private final CobalaminCrossMemberCmScorer.Worker cross;
		private final ByteBuilder row=new ByteBuilder(1024);
		private final long[] costs=new long[3*3*COST_FIELDS];
		private long completed;
		private boolean success;
		private Throwable error;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** First library ordinal wins ties, including an all-zero profile vector. */
	public static int bestModel(final float[] scores){
		require(scores!=null && scores.length>0, "Target selection needs the complete nonempty native profile-score axis");
		return CobalaminCrossMemberCmScorer.bestModel(scores, 0, scores.length);
	}
	private static int member(final String name){
		for(int i=0; i<3; i++){if(MEMBERS[i].equals(name)){return i;}}
		throw new IllegalArgumentException("Unknown cobalamin member: "+name);
	}
	private static boolean validScore(final float score){return Float.isNaN(score) || Float.isFinite(score) || score==Float.NEGATIVE_INFINITY;}
	private static String first(final String id){return RrnaResourceIO.first(id);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final class Source{
		Source(final String fasta_){fasta=fasta_;}
		final String fasta;
		final ArrayList<Event> events=new ArrayList<Event>();
	}
	private static final class Group{final ArrayList<Event> events=new ArrayList<Event>();}
	private static final class Event{
		Event(final long index_, final byte[] line_, final String contig_, final int member_, final int model_,
				final int length_, final int start_, final int stop_, final boolean minus_, final float score_, final float ga_){
			index=index_;
			line=line_;
			contig=contig_;
			member=member_;
			model=model_;
			length=length_;
			start=start_;
			stop=stop_;
			minus=minus_;
			score=score_;
			ga=ga_;
		}
		final long index;
		final byte[] line;
		final String contig;
		final int member, model, length, start, stop;
		final boolean minus;
		final float score, ga;
		boolean done;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, inputs, out, familyRoot, scoreRoot;
	private String[] profiles;
	private int expectedInputs;
	private float low=.60f, high=.778f;
	private final NcrnaFamily[] families=new NcrnaFamily[3];
	private final CobalaminCmEndpointRefiner.Resources[] factories=new CobalaminCmEndpointRefiner.Resources[3];
	private CobalaminCrossMemberCmScorer crossScorer;
	private final ArrayList<ObjectIntMap<String>> modelIds=new ArrayList<ObjectIntMap<String>>();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final AtomicInteger next=new AtomicInteger();
	private Output output;
	private long eventsRead, bandEvents, completed;
	private final long[] totals=new long[3*3*COST_FIELDS];
	private boolean success=true;
	private Throwable error;
	private volatile boolean stop;
	private static final int COST_FIELDS=8;
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	public static final String SCORE_HEADER="capture_event_index\t"+CobalaminCmAcceptance.HEADER
		+"\tscoring_member\tanchor_model\ttarget_cm_score_bits\ttarget_ga_bits\tmeasurement_status\textra_dp_fills"
		+"\textra_ranking_cpu_ns\textra_anchor_cpu_ns\textra_cm_cpu_ns\textra_total_cpu_ns";
}
