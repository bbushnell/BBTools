package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Captures ordinary RF00174 windows before shortlisting, without another caller.
 * Adapted from template.A_SampleBasic2 for CLI/output ownership and the existing
 * SweepNcrnaShortlist whole-input integration. CallGenes owns sequence readers,
 * PGMs and workers. Each read is called once with its configured competition,
 * endpoints and path selection. Capture observes no additional alignments.
 * Metadata fields use BuildNcrnaReference's lossless escaped UTF-8 TSV codec.
 * Contig titles are the actual strings observed by the native caller.
 * Windows retain zero-based inclusive oriented and original-forward coordinates.
 * Strand 1's oriented coordinates describe the reverse-complemented contig.
 * Physical-edge retries reuse the primary event's captured window bytes.
 * Optional selected-profile results retain completed work and final calls with
 * their original input identities. GC counts use ACGT bases only; zero known
 * bases is unavailable. This collector never filters or grades the gene roster.
 * @author Brian Bushnell, Raiden
 */
public final class CaptureNcrnaWindows implements CallGenesIteration{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final CaptureNcrnaWindows tool=new CaptureNcrnaWindows(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	public CaptureNcrnaWindows(final String[] args){
		final PreParser pp=new PreParser(args, CaptureNcrnaWindows.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){
				callerConfig=b;
			}else if(a.equals("reference")){
				reference=b;
			}else if(a.equals("expectedmodels")){
				expectedModels=Parse.parseIntKMG(b);
			}else if(a.equals("results")){
				verification=Parse.parseBoolean(b);
			}else if(a.equals("replaycontext")){
				replayContext=Parse.parseBoolean(b);
			}else if(a.equals("t") || a.equals("threads")){
				threads=Parse.parseIntKMG(b);
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		prefix=parser.out1;
		require(callerConfig!=null && reference!=null && prefix!=null,
			"Supply caller=config reference=original_input_manifest out=prefix");
		require(expectedModels>0 && (threads==-1 || threads>0), "Model and explicit worker counts must be positive");
		require(!parser.overwrite && !parser.append && parser.in1==null && parser.in2==null
			&& parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Capture owns fresh reports; input selection and read limits belong in the recorded caller configuration");
		require(Tools.testInputFiles(false, true, callerConfig, reference), "Caller configuration and original-input manifest must be readable");
		require(Tools.testOutputFiles(false, false, false, prefix+".windows.tsv", prefix+".windows.fa.gz", prefix+".inputs.tsv", prefix+".summary.tsv"),
			"Capture reports must be fresh and writable");
		if(replayContext){require(Tools.testOutputFiles(false, false, false, prefix+".context.tsv", prefix+".context.fa.gz"),
			"Replay context is a separate fresh sidecar; candidate-window bytes and coordinates remain unchanged");}
		if(verification){require(Tools.testOutputFiles(false, false, false, prefix+".results.tsv", prefix+".final_calls.tsv", prefix+".calls.gff", prefix+".attempts.tsv"),
			"Selected verification results and final calls must be fresh and writable");}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final String[] config=Parser.parseConfig(new String[]{"config="+
			NcrnaCallerProcess.configurationPath(new java.io.File(callerConfig))});
		final ArrayList<String> arguments=new ArrayList<String>();
		boolean taxonomyOff=false;
		for(final String arg:config){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("taxonomy") || a.equals("tax")){taxonomyOff=!Parse.parseBoolean(b);}
			arguments.add(arg);
		}
		require(taxonomyOff, "Capture requires final taxonomy=f or tax=f in the recorded caller configuration");
		arguments.add(verification ? "out="+prefix+".calls.gff" : "out=null");
		arguments.add("outstats=null");
		if(threads>0){arguments.add("t="+threads);}
		try(Output metadata=new Output(prefix+".windows.tsv");
			Output sequences=new Output(prefix+".windows.fa.gz");
			Output contextMetadata=replayContext ? new Output(prefix+".context.tsv") : null;
			Output contextSequences=replayContext ? new Output(prefix+".context.fa.gz") : null;
			Output completed=verification ? new Output(prefix+".results.tsv") : null;
			Output attempted=verification ? new Output(prefix+".attempts.tsv") : null;
			Output retained=verification ? new Output(prefix+".final_calls.tsv") : null){
			rows=metadata;
			fasta=sequences;
			contextRows=contextMetadata;
			contextFasta=contextSequences;
			results=completed;
			attempts=attempted;
			finalCalls=retained;
			rows.writer.println(WINDOW_HEADER);
			if(replayContext){contextRows.writer.println(CONTEXT_HEADER);}
			if(verification){
				results.writer.println(RESULT_HEADER);
				attempts.writer.println(ATTEMPT_HEADER);
				finalCalls.writer.println(FINAL_HEADER);
			}
			final CallGenes nativeCaller=new CallGenes(arguments.toArray(new String[arguments.size()]));
			nativeCaller.run(timer, ()->{
				if(verification){require(nativeCaller.hasSelectedCobalaminDdl(), "results=t requires selected-format DDL profiles before starting input workers");}
				bind(nativeCaller.inputFiles());
				nativeCaller.setIteration(this);
			});
		}finally{
			rows=null;
			fasta=null;
			contextRows=null;
			contextFasta=null;
			results=null;
			attempts=null;
			finalCalls=null;
		}
		writeInputs();
		writeSummary();
		outstream.println("NCRNA_WINDOW_CAPTURE_COMPLETE inputs="+bindings.length+" windows="+totalWindows);
		outstream.println("Coordinates are zero-based inclusive in both caller-strand and original-forward orientation; physical-edge retries reuse their primary captured window.");
		outstream.println("One configured native RF00174 caller pass; no separate oracle or verification replay.");
	}

	/** Binds every configured path once; zero-window inputs keep their original ID. */
	private void bind(final String[] inputs){
		require(inputs.length>0, "The native caller must have an explicit input population");
		final ObjectMap<String, Binding> population=new ObjectMap<String, Binding>(String.class, Binding.class);
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){
			population.put(binding.fasta, binding);
		}
		bindings=new Binding[inputs.length];
		counts=new long[inputs.length][7];
		finished=new boolean[inputs.length];
		for(int i=0; i<inputs.length; i++){
			final Binding binding=population.get(inputs[i]);
			require(binding!=null && inputIndex.get(inputs[i])<0,
				"Every caller input must bind one unique original manifest FASTA: "+inputs[i]);
			bindings[i]=binding;
			inputIndex.put(inputs[i], i);
		}
		for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
			final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(i);
			if(!family.name.equals(TARGET)){continue;}
			require(targetIndex<0 && family.library.length==expectedModels,
				"Capture requires one RF00174 family with its declared consensus-model count");
			targetIndex=i;
		}
		require(targetIndex>=0, "The recorded caller configuration must enable RF00174");
	}

	@Override
	public Worker create(final String input, final GeneCaller caller, final GeneModel model){
		final int index=inputIndex.get(input);
		require(index>=0 && !finished[index], "A capture worker must own one declared original input: "+input);
		return new CaptureWorker(index, caller, model);
	}

	/** Shared native writers own mutable buffers, so publish each complete pair under one lock. */
	private synchronized void publish(final ByteBuilder metadata, final ByteBuilder sequence,
			final ByteBuilder contextMetadata, final ByteBuilder contextSequence){
		assert(rows!=null && fasta!=null && (!replayContext || contextRows!=null && contextFasta!=null))
			: "Candidate and context owners must outlive every native caller worker";
		rows.writer.print(metadata);
		fasta.writer.print(sequence);
		if(replayContext){
			contextRows.writer.print(contextMetadata);
			contextFasta.writer.print(contextSequence);
		}
		totalWindows++;
	}

	private synchronized void publishResult(final ByteBuilder row){
		assert(verification && results!=null) : "Completed selected windows must be written before the result owner closes";
		results.writer.print(row);
		totalResults++;
	}

	private synchronized void publishAttempt(final ByteBuilder row){
		assert(verification && attempts!=null) : "Attempt summaries must finish inside their native capture output scope";
		attempts.writer.print(row);
	}

	private synchronized void publishFinal(final ByteBuilder row){
		assert(verification && finalCalls!=null) : "Retained calls must keep their original input binding while the output owner is live";
		finalCalls.writer.print(row);
		totalFinalCalls++;
	}

	/** Published only after native workers and both event outputs close successfully. */
	private void writeInputs(){
		try(Output output=new Output(prefix+".inputs.tsv")){
			output.writer.println(verification ? "input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tscheduled_windows\talignments\taccepted_proposals\tfinal_target_call_rows"
				: "input_slot\tinput_id\tfasta\tcontigs\tbases\twindows");
			final ByteBuilder buffer=new ByteBuilder();
			long windows=0;
			for(int i=0; i<bindings.length; i++){
				require(finished[i], "An input did not finish; do not publish a complete capture manifest: "+bindings[i].id);
				buffer.clear().append(i).tab();
				BuildNcrnaReference.field(buffer, bindings[i].id).tab();
				BuildNcrnaReference.field(buffer, bindings[i].fasta).tab();
				buffer.append(counts[i][0]).tab().append(counts[i][1]).tab().append(counts[i][2]);
				if(verification){for(int metric=3; metric<7; metric++){buffer.tab().append(counts[i][metric]);}}
				buffer.nl();
				output.writer.print(buffer);
				windows+=counts[i][2];
			}
			require(windows==totalWindows, "Per-input capture counts must conserve every published event");
		}
	}

	/** Final completion table states the coordinate and unobserved-retry contract. */
	private void writeSummary(){
		assert(bindings!=null && counts!=null && bindings.length==counts.length)
			: "Completion accounting requires every declared input's joined worker counts";
		long contigs=0, bases=0;
		for(final long[] input:counts){
			contigs+=input[0];
			bases+=input[1];
		}
		try(Output output=new Output(prefix+".summary.tsv")){
			final ByteBuilder buffer=new ByteBuilder("item\tvalue\nfamily\t").append(TARGET).nl();
			buffer.append("models\t").append(expectedModels).append("\ninputs\t").append(bindings.length)
				.append("\ncontigs\t").append(contigs).append("\nbases\t").append(bases)
				.append("\ncaptured_primary_windows\t").append(totalWindows).nl();
			buffer.append("coordinates\tZERO_BASED_INCLUSIVE\n")
				.append("oriented_coordinates\tCALLER_STRAND_0_FORWARD_1_REVERSE_COMPLEMENT\n")
				.append("forward_coordinates\tORIGINAL_FORWARD_CONTIG\n")
				.append("capture_stage\tAFTER_SEED_GATE_BEFORE_SHORTLIST\n")
				.append("edge_retries_separately_observed\tfalse\n")
				.append("edge_retry_window_policy\tPRIMARY_CAPTURED_BYTES_REUSED\n");
			if(replayContext){buffer.append("replay_context\tREAL_ORIENTED_CANDIDATE_PLUS_MAXIMUM_REFINEMENT_FLANK\n")
				.append("replay_context_coordinates\tSEPARATE_FROM_CANDIDATE_AND_ORIGINAL_PHYSICAL_BOUNDARIES\n");}
			if(verification){
				long scheduled=0, alignments=0, accepted=0, calls=0;
				for(final long[] input:counts){
					scheduled+=input[3];
					alignments+=input[4];
					accepted+=input[5];
					calls+=input[6];
				}
				require(scheduled==totalResults, "Every scheduled selected window must retain its completed result, including uncaptured gate rejections");
				require(calls==totalFinalCalls, "Every final target call must retain its original input binding, even when different inputs reuse contig names");
				buffer.append("completed_scheduled_windows\t").append(scheduled).nl()
					.append("executed_alignments\t").append(alignments).nl()
					.append("accepted_proposal_events\t").append(accepted).nl()
					.append("final_target_call_rows\t").append(calls).nl()
					.append("verification_stage\tAFTER_ALL_WINDOW_ALIGNMENT_WORK_BEFORE_CANDIDATE_RESOLUTION_AND_PATH_SELECTION\n")
					.append("alignment_cache\tDISABLED\n")
					.append("gc_denominator\tACGT_BASES_ONLY_ZERO_IS_UNAVAILABLE\n");
			}
			output.writer.print(buffer);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Input Worker          ----------------*/
	/*--------------------------------------------------------------*/

	private final class CaptureWorker implements Worker, NcrnaWindowCaptureSink, NcrnaModelAttemptInstrumentSink{
		CaptureWorker(final int index_, final GeneCaller caller_, final GeneModel model_){
			index=index_;
			caller=caller_;
			model=model_;
			scavenger=caller.ncrnaScavenger(targetIndex);
			contextFlank=replayContext ? scavenger.savedWindowFlank() : 0;
			if(verification){
				require(scavenger.hasSelectedDdlProfiles(), "results=t requires the selected-format caller with measured accepted identity");
				scavenger.setAlignmentCache(null);
				scavenger.setModelAttemptSink(this);
			}
			scavenger.setWindowCaptureSink(this);
		}

		@Override
		public ArrayList<Orf> call(final Read read){
			require(current==null && read!=null && read.bases!=null && read.mate==null,
				"Capture requires one unchanged, unpaired native contig per call");
			current=read;
			try{
				final ArrayList<Orf> calls=caller.callGenes(read, model, true);
				if(verification){
					require(activeEvent<0 && primary==0 && rescue==0, "Every observed contig must finish its last complete window before final-call accounting");
					for(final Orf call:calls){if(TARGET.equals(call.ncrnaFamily)){retain(call);}}
				}
				return calls;
			}
			finally{
				counts[index][0]++;
				counts[index][1]+=read.length();
				current=null;
			}
		}

		/** Final path coordinates are already original-forward; never join combined
		 * GFF rows by contig name alone when different input files reuse that name. */
		private void retain(final Orf call){
			require(call instanceof NcrnaDdlOrf && current!=null, "A final selected call must keep its measured winner and current original input");
			final NcrnaDdlOrf proposal=(NcrnaDdlOrf)call;
			completedMetadata.clear().append(index).tab();
			BuildNcrnaReference.field(completedMetadata, bindings[index].id).tab();
			BuildNcrnaReference.field(completedMetadata, bindings[index].fasta).tab();
			completedMetadata.append(counts[index][0]).tab();
			BuildNcrnaReference.field(completedMetadata, current.id).tab();
			completedMetadata.append(call.start).tab().append(call.stop).tab().append(call.strand).tab();
			BuildNcrnaReference.field(completedMetadata, call.trnaModel).tab();
			completedMetadata.append(proposal.acceptedIdentity(), 9).tab().append(proposal.proposalRank())
				.tab().append(proposal.proposalScore(), 3).nl();
			publishFinal(completedMetadata);
			counts[index][6]++;
		}

		@Override
		public void capture(final String contig, final int strand, final int pass,
				final int start, final int stop, final int source, final byte[] bases, final byte[] window){
			assert(current!=null && bases==current.bases && bases.length==current.length()
				&& (strand==0 || strand==1) && start>=0 && stop>=start
				&& stop<bases.length && window.length==stop-start+1)
				: "A captured window must belong to the current original contig and retain its exact inclusive span";
			final int forwardStart=strand==0 ? start : bases.length-1-stop;
			final int forwardStop=strand==0 ? stop : bases.length-1-start;
			final long event=counts[index][2]++;
			if(verification){
				require(activeEvent<0, "The preceding window must complete before another captured event can reuse observer scratch");
				activeEvent=event;
				gc=0;
				acgt=0;
				for(final byte base:window){
					if(base=='G' || base=='C'){gc++;acgt++;}
					else if(base=='A' || base=='T'){acgt++;}
				}
			}
			metadata.clear().append('i').append(index).append('e').append(event).tab().append(index).tab();
			BuildNcrnaReference.field(metadata, bindings[index].id).tab();
			BuildNcrnaReference.field(metadata, bindings[index].fasta).tab();
			metadata.append(counts[index][0]).tab();
			BuildNcrnaReference.field(metadata, contig).tab();
			metadata.append(strand).tab().append(bases.length).tab().append(pass).tab()
				.append(start).tab().append(stop).tab().append(forwardStart).tab().append(forwardStop)
				.tab().append(source).tab().append(window.length).nl();
			sequence.clear().append('>').append('i').append(index).append('e').append(event).nl().append(window).nl();
			if(replayContext){
				final int contextStart=Math.max(0, start-contextFlank);
				final int contextStop=(int)Math.min((long)bases.length-1, (long)stop+contextFlank);
				contextMetadata.clear().append('i').append(index).append('e').append(event).tab()
					.append(contextStart).tab().append(contextStop).tab().append(contextFlank).nl();
				contextSequence.clear().append('>').append('i').append(index).append('e').append(event).nl()
					.append(bases, contextStart, contextStop-contextStart+1).nl();
			}
			publish(metadata, sequence, contextMetadata, contextSequence);
		}

		@Override
		public void modelAttempt(String name, int strand, int pass, int start, int stop,
				int rank, int model, int shared, float identity, int length, boolean accepted){
			assert(verification && current!=null && activeEvent>=0)
				: "A primary or edge attempt must belong to the currently captured native window";
			primary++;
		}

		@Override
		public void rescueAttempt(String name, int strand, int pass, int start, int stop,
				int candidateStart, int candidateStop, int rank, int model, float identity,
				int length, float hbm, boolean hbmMeasured, boolean accepted){
			assert(verification && current!=null && activeEvent>=0)
				: "A rescue attempt must retain its primary captured-window identity";
			assert(!hbmMeasured || Float.isFinite(hbm))
				: "A measured rescue HBM score is finite; this observer counts attempts without serializing that score";
			rescue++;
		}

		@Override
		public void alignmentSpan(String name, int strand, int pass, int start, int stop,
				int modelIndex, float identity, int alignedStart, int alignedStop, int stage,
				boolean genuineClip, boolean edgeFloor){
			require(verification && current!=null && activeEvent>=0 && stage>=0 && stage<3,
				"Raw attempt spans must belong to their captured native window and explicit verifier stage");
			observed[stage]++;
			if(!Float.isFinite(identity) || alignedStart<0 || alignedStop<alignedStart){return;}
			require(alignedStart>=start && alignedStop<=stop, "Raw genomic attempt endpoints must stay within the actual verifier window");
			if(bestModel[stage]<0 || identity>bestIdentity[stage]){
				bestModel[stage]=modelIndex;
				bestIdentity[stage]=identity;
				bestStart[stage]=alignedStart;
				bestStop[stage]=alignedStop;
				bestClip[stage]=genuineClip;
				bestEdgeFloor[stage]=edgeFloor;
			}
		}

		@Override
		public void completed(String contig, int strand, int pass, int start, int stop, Orf accepted, long alignments){
			if(!verification){return;}
			require(current!=null && alignments>=primary+rescue && (activeEvent>=0 || alignments==0 && accepted==null),
				"Whole-window work must include every primary/rescue attempt; only uncaptured gate rejections can have no event");
			final int length=current.length();
			final int forwardStart=strand==0 ? start : length-1-stop, forwardStop=strand==0 ? stop : length-1-start;
			completedMetadata.clear().append(index).tab();
			BuildNcrnaReference.field(completedMetadata, bindings[index].id).tab();
			BuildNcrnaReference.field(completedMetadata, bindings[index].fasta).tab();
			completedMetadata.append(counts[index][0]).tab();
			BuildNcrnaReference.field(completedMetadata, contig).tab();
			completedMetadata.append(strand).tab().append(length).tab().append(pass).tab()
				.append(start).tab().append(stop).tab().append(forwardStart).tab().append(forwardStop).tab();
			if(activeEvent>=0){completedMetadata.append('i').append(index).append('e').append(activeEvent);}
			else{completedMetadata.append("NA");}
			completedMetadata.tab().append(accepted!=null ? "ACCEPTED_PROPOSAL" : activeEvent>=0 ? "REJECTED_VERIFICATION" : "UNCAPTURED_GATE_REJECTION")
				.tab().append(primary).tab().append(rescue).tab().append(alignments-primary-rescue).tab().append(alignments).tab();
			if(accepted==null){completedMetadata.append("NA\tNA\tNA\tNA\tNA\tNA");}
			else{
				require(accepted instanceof NcrnaDdlOrf && TARGET.equals(accepted.ncrnaFamily),
					"Selected results must retain the actual target proposal subtype");
				final NcrnaDdlOrf proposal=(NcrnaDdlOrf)accepted;
				require(Float.isFinite(proposal.acceptedIdentity()), "The accepting model must retain its actual measured identity");
				BuildNcrnaReference.field(completedMetadata, accepted.trnaModel).tab();
				completedMetadata.append(proposal.acceptedIdentity(), 9).tab().append(proposal.proposalRank()).tab()
					.append(proposal.proposalScore(), 3).tab()
					.append(strand==0 ? accepted.start : length-1-accepted.stop).tab()
					.append(strand==0 ? accepted.stop : length-1-accepted.start);
				counts[index][5]++;
			}
			completedMetadata.tab().append(activeEvent<0 ? -1 : gc).tab().append(activeEvent<0 ? -1 : acgt).nl();
			publishResult(completedMetadata);
			if(activeEvent>=0){
				for(int stage=0; stage<3; stage++){
					attemptMetadata.clear().append(index).tab();
					BuildNcrnaReference.field(attemptMetadata, bindings[index].id).tab();
					attemptMetadata.append('i').append(index).append('e').append(activeEvent).tab()
						.append(stage==PRIMARY ? "PRIMARY_RAW" : stage==PHYSICAL_EDGE ? "PHYSICAL_EDGE_RAW" : "RESCUE_RAW").tab()
						.append(observed[stage]).tab().append(accepted!=null).tab();
					if(bestModel[stage]<0){attemptMetadata.append("NA\tNA\tNA\tNA\tNA\tNA\tNA");}
					else{
						BuildNcrnaReference.field(attemptMetadata, GeneCaller.ncrnaFamilies.get(targetIndex).modelNames[bestModel[stage]]).tab();
						attemptMetadata.append(bestIdentity[stage], 9).tab().append(bestStop[stage]-bestStart[stage]+1).tab()
							.append(strand==0 ? bestStart[stage] : length-1-bestStop[stage]).tab()
							.append(strand==0 ? bestStop[stage] : length-1-bestStart[stage]).tab()
							.append(bestClip[stage]).tab().append(bestEdgeFloor[stage]);
					}
					publishAttempt(attemptMetadata.nl());
				}
			}
			counts[index][3]++;
			counts[index][4]+=alignments;
			primary=0;
			rescue=0;
			activeEvent=-1;
			java.util.Arrays.fill(bestModel, -1);
			java.util.Arrays.fill(observed, 0);
		}

		@Override
		public void finish(){
			assert(current==null && !finished[index]) : "Each whole-input worker must finish once after releasing its borrowed contig";
			scavenger.setWindowCaptureSink(null);
			if(verification){
				scavenger.setModelAttemptSink(null);
				require(counts[index][4]==scavenger.alignmentCount(), "Completed-window rows must conserve every actual target alignment on this input");
			}
			finished[index]=true;
		}

		private final int index;
		private final GeneCaller caller;
		private final GeneModel model;
		private final NcrnaScavenger scavenger;
		private final int contextFlank;
		private final ByteBuilder metadata=new ByteBuilder(), sequence=new ByteBuilder();
		private final ByteBuilder contextMetadata=new ByteBuilder(), contextSequence=new ByteBuilder();
		private final ByteBuilder completedMetadata=new ByteBuilder();
		private final ByteBuilder attemptMetadata=new ByteBuilder();
		private final int[] bestModel={-1,-1,-1}, bestStart=new int[3], bestStop=new int[3];
		private final long[] observed=new long[3];
		private final float[] bestIdentity=new float[3];
		private final boolean[] bestClip=new boolean[3], bestEdgeFloor=new boolean[3];
		private long activeEvent=-1, primary, rescue;
		private int gc, acgt;
		private Read current;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	private String callerConfig, reference, prefix;
	private int threads=-1, expectedModels=508, targetIndex=-1;
	private boolean verification, replayContext;
	private final PrintStream outstream;
	private final ObjectIntMap<String> inputIndex=new ObjectIntMap<String>(String.class);
	private Binding[] bindings;
	private long[][] counts;
	private boolean[] finished;
	private Output rows, fasta, results, finalCalls, attempts, contextRows, contextFasta;
	private long totalWindows, totalResults, totalFinalCalls;
	private static final String TARGET="cobalamin_rf00174";
	public static final String WINDOW_HEADER="event_id\tinput_slot\tinput_id\tfasta\tcontig_ordinal\tcontig\tstrand\tcontig_length\tpass\toriented_start0\toriented_stop0\tforward_start0\tforward_stop0\tsource\twindow_length";
	public static final String CONTEXT_HEADER="event_id\toriented_start0\toriented_stop0\tmaximum_refinement_flank";
	public static final String RESULT_HEADER="input_slot\tinput_id\tfasta\tcontig_ordinal\tcontig\tstrand\tcontig_length\tpass\toriented_start0\toriented_stop0\tforward_start0\tforward_stop0\tcaptured_event_id\tstatus\tprimary_and_edge_alignments\trescue_alignments\tother_alignments\ttotal_alignments\tmodel\taccepted_identity\tddl_rank\tddl_score\taccepted_forward_start0\taccepted_forward_stop0\tgc_bases\tacgt_bases";
	public static final String FINAL_HEADER="input_slot\tinput_id\tfasta\tcontig_ordinal\tcontig\tforward_start0\tforward_stop0\tstrand\tmodel\taccepted_identity\tddl_rank\tddl_score";
	public static final String ATTEMPT_HEADER="input_slot\tinput_id\tcaptured_event_id\tstage\tobserved_attempts\twindow_returned_acceptance\tbest_raw_model\tbest_raw_identity\tbest_raw_span\tbest_raw_forward_start0\tbest_raw_forward_stop0\tgenuine_model_clip\tedge_identity_floor_pass";
}
