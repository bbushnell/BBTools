package prok;

import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlTruth.Group;
import prok.NcrnaDdlTruth.Parent;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaTextIO.Outputs;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Compares member-voted DDL rankings on the ordinary caller's RF00174 windows.
 * CLI/output lifecycle follows template.A_SampleBasic2 and CaptureNcrnaWindows;
 * CallGenesIteration owns the native readers, PGMs, workers and final gene paths.
 * One normal call per original Read lends windows synchronously to the accepted
 * NcrnaDdlWindowAnalysis worker. No second genome reader or caller is created.
 * Each model is verified once/window, then reused across every profile geometry.
 * All scored cobalamin parents participate in conservative attribution; only
 * qualified RF00174 parents enter rank denominators. Ambiguous spans remain
 * explicit links and never become unique-recovery credit.
 * This is an in-corpus diagnostic, not a caller-policy change or speed claim.
 * @author Brian Bushnell, Raiden
 */
public final class RankNcrnaWindows implements CallGenesIteration{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final RankNcrnaWindows tool=new RankNcrnaWindows(args);
		try{tool.process(new Timer());}
		finally{Shared.closeStream(tool.outstream);}
	}

	public RankNcrnaWindows(final String[] args){
		final PreParser pp=new PreParser(args, RankNcrnaWindows.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){callerConfig=b;}
			else if(a.equals("reference")){reference=b;}
			else if(a.equals("familymap")){familymap=b;}
			else if(a.equals("subtypes")){subtypes=b;}
			else if(a.equals("witnesses")){witnesses=b;}
			else if(a.equals("events")){events=b;}
			else if(a.equals("members")){members=b;}
			else if(a.equals("ledger")){ledger=b;}
			else if(a.equals("assignments")){assignments=b;}
			else if(a.equals("census")){census=b;}
			else if(a.equals("membercounts")){membercounts=b;}
			else if(a.equals("expectedmodels")){models=Parse.parseIntKMG(b);}
			else if(a.equals("ks")){ks=parseGrid(b, true);}
			else if(a.equals("widths")){widths=parseGrid(b, false);}
			else if(a.equals("t") || a.equals("threads")){threads=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown ranking parameter: "+arg);}
		}
		prefix=parser.out1;
		require(prefix!=null && models>0 && threads>0 && !parser.overwrite && !parser.append
			&& parser.in1==null && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Use out=fresh_prefix t=positive_count; complete inputs belong in caller=config without report sampling");
		require(Tools.testInputFiles(false, true, callerConfig, reference, familymap, subtypes, witnesses,
			events, members, ledger, assignments, census, membercounts),
			"Supply all caller, reference-support and frozen construction-member bindings");
		for(final String path:new String[]{callerConfig, reference, familymap, subtypes, witnesses,
			events, members, ledger, assignments, census, membercounts}){
			require(path!=null, "Every declared provenance input is mandatory");
		}
		NcrnaDdlWindowAnalysis.checkEnvironment();
		final String[] paths=new String[SUFFIX.length+1];
		for(int i=0; i<SUFFIX.length; i++){paths[i]=prefix+SUFFIX[i];}
		paths[SUFFIX.length]=prefix+".calls.gff";
		require(Tools.testOutputFiles(false, false, false, paths), "All ranking and ordinary-call outputs must be fresh");
	}

	/** Strictly increasing grids avoid silently counting a configuration twice. */
	private static int[] parseGrid(final String text, final boolean k){
		require(text!=null && !text.isEmpty(), "Declare a nonempty k or width grid");
		final String[] fields=text.split(",");
		final int[] result=new int[fields.length];
		for(int i=0; i<fields.length; i++){
			result[i]=Parse.parseIntKMG(fields[i]);
			require((k ? result[i]>=1 && result[i]<=31 :
				result[i]>=32 && result[i]<=2048 && Integer.bitCount(result[i])==1)
				&& (i==0 || result[i]>result[i-1]), "Grid values must be valid, distinct and increasing");
		}
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Native Run          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final String[] config=Parser.parseConfig(new String[]{"config="+
			NcrnaCallerProcess.configurationPath(new java.io.File(callerConfig))});
		final ArrayList<String> arguments=new ArrayList<String>();
		final Parser selection=new Parser();
		boolean taxonomyOff=false;
		for(final String arg:config){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			//TODO: Probable bug - Parser.parseConfig expands one level, then CallGenes'
			// PreParser expands again. A nested config can therefore hide read limits
			// from this selection guard and invalidate complete-input denominators.
			// Reject nested config here or flatten it before checking; project recipes
			// currently supply flat caller files. Parser.parseConfig:679-720.
			if(a.equals("tax") || a.equals("taxonomy")){taxonomyOff=!Parse.parseBoolean(b);}
			if(a.equals("reads") || a.equals("maxreads") || a.equals("samplerate")){
				selection.parseCommon(arg, a, b);
			}
			arguments.add(arg);
		}
		require(taxonomyOff, "The final recorded tax/taxonomy value must be false; native taxonomy changes sketch encoding");
		require(selection.maxReads<0 && selection.samplerate==1f,
			"Caller configuration must not limit or sample contigs; every declared reference target needs a complete input");
		arguments.add("out="+prefix+".calls.gff");
		arguments.add("outstats=null");
		arguments.add("ow=f");
		arguments.add("t="+threads);
		final int previous=DynamicDemiLog.exponentBits();
		try(Outputs outputs=new Outputs(SUFFIX.length)){
			for(int i=0; i<SUFFIX.length; i++){
				writers[i]=outputs.open(fileIO.FileFormat.testOutput(prefix+SUFFIX[i],
					fileIO.FileFormat.TXT, null, true, false, false, false));
			}
			headers();
			final CallGenes caller=new CallGenes(arguments.toArray(new String[arguments.size()]));
			caller.run(timer, ()->{
				DynamicDemiLog.setExponent(4);
				bind(caller.inputFiles());
				caller.setIteration(this);
			});
			final long beforeReports=cpu();
			report();
			final ByteBuilder reportCost=new ByteBuilder("final_report_owner_cpu_ns\t").append(cpu()-beforeReports).nl();
			publish(WORK, reportCost);
		}finally{
			Arrays.fill(writers, null);
			DynamicDemiLog.setExponent(previous);
		}
		outstream.println("NCRNA_DDL_WINDOW_RANKING_COMPLETE_NOT_CALLER_ADOPTED");
		outstream.println("Whole-command time includes member profiles, exhaustive verification, rankings and reports.");
	}

	private void bind(final String[] inputs){
		require(inputs!=null && inputs.length>0, "Rank the complete explicit native input population");
		for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
			final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(i);
			if(!TARGET.equals(f.name)){continue;}
			require(family==null && f.library.length==models, "Exactly one RF00174 family must retain the declared model order");
			family=f;
			targetIndex=i;
		}
		require(family!=null, "The recorded ordinary configuration must enable RF00174");
		final NcrnaDdlTruth.Input[] truth=NcrnaDdlTruth.load(reference, inputs, subtypes, witnesses,
			events, GradeNcrnaReference.loadRules(familymap));
		data=new Input[inputs.length];
		for(int i=0; i<inputs.length; i++){
			require(inputIndex.get(inputs[i])<0, "Each native input may be called exactly once");
			inputIndex.put(inputs[i], i);
			data[i]=new Input(truth[i], ks.length*widths.length*3+1, models);
		}
		final long beforeMembers=cpu();
		final NcrnaDdlMembers bound=NcrnaDdlMembers.load(members, ledger, assignments, census,
			membercounts, family.modelNames, family.library);
		memberCpu=cpu()-beforeMembers;
		profiles=new NcrnaDdlProfile[ks.length*widths.length][];
		geometryK=new int[profiles.length];
		profileGeometryCpu=new long[profiles.length];
		int geometry=0;
		for(final int k:ks){
			for(final int width:widths){
				final long before=cpu();
				geometryK[geometry]=k;
				profiles[geometry]=NcrnaDdlModelProfiles.build(bound, k, width);
				profileGeometryCpu[geometry]=cpu()-before;
				profileCpu+=profileGeometryCpu[geometry++];
				profileKeyBytes+=(long)width*models*Character.BYTES;
				profileSupportBytes+=(long)width*models*Integer.BYTES;
			}
		}
		inputMembers=bound.inputMembers();
		keptMembers=bound.keptMembers();
		orphanMembers=bound.orphanMembers();
	}

	@Override
	public Worker create(final String input, final GeneCaller caller, final GeneModel model){
		final int index=inputIndex.get(input);
		require(index>=0 && !data[index].finished, "A worker needs one exact declared input");
		return new InputWorker(index, caller, model);
	}

	/** One coarse handoff per table/window; no per-model writer lock or task. */
	private synchronized void publish(final int output, final ByteBuilder buffer){
		assert(writers[output]!=null) : "Native workers must finish before the report-output owners close";
		if(buffer.length()>0){writers[output].print(buffer);}
		buffer.clear();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Input Worker        ----------------*/
	/*--------------------------------------------------------------*/

	private final class InputWorker implements Worker, NcrnaDdlWindowAnalysis.Sink{
		InputWorker(final int index_, final GeneCaller caller_, final GeneModel model_){
			index=index_;
			result=data[index];
			caller=caller_;
			model=model_;
			analysis=new NcrnaDdlWindowAnalysis(result.truth.binding.id, family,
				caller.ncrnaScavenger(targetIndex), profiles, geometryK, this);
			for(int i=0; i<buffers.length; i++){buffers[i]=new ByteBuilder();}
			current=new LocusRanks[result.loci.length];
			flags=new boolean[result.loci.length][models];
		}

		@Override
		public ArrayList<Orf> call(final Read read){
			analysis.begin(read);
			try{return caller.callGenes(read, model, true);}
			finally{
				analysis.end();
				result.contigs++;
				result.bases+=read.length();
			}
		}

		@Override
		public void window(final NcrnaDdlWindowAnalysis view, final byte[] window){
			final long before=cpu(), rankedBefore=view.rankingCpu();
			final ByteBuilder windows=buffers[WINDOWS], sequence=buffers[SEQUENCES];
			event(windows, view).tab();
			field(windows, view.header()).tab().append(view.strand()).tab().append(view.pass()).tab()
				.append(view.start()).tab().append(view.stop()).tab()
				.append(view.strand()==0 ? view.start() : view.contigLength()-1-view.stop()).tab()
				.append(view.strand()==0 ? view.stop() : view.contigLength()-1-view.start()).tab()
				.append(view.contigLength()).tab().append(view.source()).tab().append(window.length).nl();
			sequence.append('>').append(index).append('_').append(view.event()).nl().append(window).nl();
			final Group group=result.truth.group(view.seqid());
			final NcrnaDdlOracle oracle=view.oracle();
			oracle.attribute(result.truth.binding.id, group==null ? NO_PARENTS : group.scored);
			for(int m=0; m<models; m++){
				final ByteBuilder out=buffers[ORACLE];
				event(out, view).tab().append(m).tab();
				field(out, family.modelNames[m]).tab().append(view.nativeSharedCount(m)).tab()
					.append(view.nativeShortlistRank(m)).tab().append(oracle.start(m)>=0).tab()
					.append(oracle.start(m)).tab().append(oracle.stop(m)).tab().append(oracle.partial(m)).tab();
				if(oracle.start(m)>=0){out.append(oracle.score(m), 6);}else{out.append('.');}
				out.tab().append(oracle.parentCount(m)).nl();
			}
			active=0;
			if(group!=null){
				for(final Parent parent:group.parents){
					if(parent.row.mask){continue;}
					final boolean[] compatible=scratch;
					final int count=oracle.compatible(result.truth.binding.id, parent.row.locus, compatible);
					for(int m=0; m<models; m++){
						if(!compatible[m]){continue;}
						final ByteBuilder out=buffers[LINKS];
						event(out, view).tab().append(m).tab();
						field(out, parent.row.id).tab().append(parent.target).tab()
							.append(oracle.parentCount(m)).tab().append(oracle.uniqueParent(m)==parent.scoredIndex).nl();
					}
					if(!parent.target){continue;}
					final Locus locus=parent.row.locus;
					final int start=view.strand()==0 ? view.start()+1 : view.contigLength()-view.stop();
					final int stop=view.strand()==0 ? view.stop()+1 : view.contigLength()-view.start();
					final boolean reached=locus.strand==(view.strand()==0 ? '+' : '-')
						&& start<=locus.stop && stop>=locus.start;
					if(!reached && count==0){continue;}
					final LocusRanks state=result.byParent.get(parent);
					state.windows++;
					if(count>0){state.compatibleWindows++;}
					final boolean[] unique=flags[active];
					int uniqueCount=0;
					for(int m=0; m<models; m++){
						unique[m]=compatible[m] && oracle.uniqueParent(m)==parent.scoredIndex;
						if(unique[m]){uniqueCount++;}
					}
					if(count>0 && uniqueCount==0){state.ambiguousWindows++;}
					if(uniqueCount>0){state.uniqueWindows++;}
					current[active++]=state;
					view.nativeBest(unique, rank);
					rankRow(view, parent, state, 0, rank, count, uniqueCount, window.length);
					boolean eligible=false;
					for(int m=0; m<models; m++){eligible|=unique[m] && view.nativeShortlistRank(m)>0;}
					state.eligible|=eligible;
					for(int g=0; g<profiles.length; g++){
						for(int method=0; method<3; method++){
							view.best(g, method, unique, rank);
							rankRow(view, parent, state, 1+g*3+method, rank, count, uniqueCount, window.length);
						}
					}
				}
			}
			for(int i=0; i<BUFFERS; i++){publish(i, buffers[i]);}
			result.reportCpu+=cpu()-before-(view.rankingCpu()-rankedBefore);
		}

		private void rankRow(final NcrnaDdlWindowAnalysis view, final Parent parent, final LocusRanks state,
				final int method, final int[] tuple, final int compatible, final int unique, final int windowLength){
			state.add(method, tuple[0]);
			final ByteBuilder out=buffers[RANKS];
			event(out, view).tab();
			field(out, parent.row.id).tab().append(method).tab().append(windowLength).tab()
				.append(parent.row.locus.stop-parent.row.locus.start+1).tab().append(compatible).tab()
				.append(unique).tab().append(tuple[0]).tab().append(tuple[1]).tab().append(tuple[2])
				.tab().append(tuple[3]).tab().append(compatible==0 ? "NO_COMPATIBLE_MODEL" : unique==0 ?
					"AMBIGUOUS_SPANS_ONLY" : tuple[0]<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED").nl();
		}

		@Override
		public void attempt(final NcrnaDdlWindowAnalysis view, final int rank_, final int model_,
				final int shared, final float identity, final int alignedLength, final boolean identityPass){
			final long before=cpu();
			final ByteBuilder out=buffers[ATTEMPTS];
			event(out, view).tab().append(rank_).tab().append(model_).tab().append(shared).tab()
				.append(Float.floatToRawIntBits(identity)).tab().append(alignedLength).tab().append(identityPass).nl();
			for(int p=0; p<active; p++){current[p].attempted|=flags[p][model_];}
			result.reportCpu+=cpu()-before;
		}

		@Override
		public void endWindow(final NcrnaDdlWindowAnalysis view){
			final long before=cpu();
			final ByteBuilder out=buffers[ENDS];
			event(out, view).tab().append(view.eventAttempts()).nl();
			publish(ATTEMPTS, buffers[ATTEMPTS]);
			publish(ENDS, out);
			Arrays.fill(current, 0, active, null);
			active=0;
			result.reportCpu+=cpu()-before;
		}

		private ByteBuilder event(final ByteBuilder out, final NcrnaDdlWindowAnalysis view){
			out.append(index).tab();
			field(out, result.truth.binding.id).tab().append(view.event()).tab()
				.append(view.contigOrdinal()).tab();
			return field(out, view.seqid());
		}

		@Override
		public void finish(){
			analysis.finish();
			result.windows=analysis.events();
			result.attempts=analysis.attempts();
			result.seedCpu=analysis.seedCpu();
			result.nativeCpu=analysis.nativeCpu();
			result.oracleCpu=analysis.oracleCpu();
			result.rankingCpu=analysis.rankingCpu();
			result.nativeRankCpu=analysis.nativeRankCpu();
			for(int geometry=0; geometry<profiles.length; geometry++){
				result.geometryScoreCpu[geometry]=analysis.geometryScoreCpu(geometry);
				result.geometryRankCpu[geometry]=analysis.geometryRankCpu(geometry);
			}
			result.finished=true;
		}

		private final int index;
		private final Input result;
		private final GeneCaller caller;
		private final GeneModel model;
		private final NcrnaDdlWindowAnalysis analysis;
		private final ByteBuilder[] buffers=new ByteBuilder[BUFFERS];
		private final LocusRanks[] current;
		private final boolean[][] flags;
		private final boolean[] scratch=new boolean[models];
		private final int[] rank=new int[4];
		private int active;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Final Reports       ----------------*/
	/*--------------------------------------------------------------*/

	private void headers(){
		final String event="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid";
		writers[WINDOWS].println(event+"\theader\tstrand0\tpass\toriented_start0\toriented_stop0\tforward_start0\tforward_stop0\tcontig_length\tsource\twindow_length");
		writers[ORACLE].println(event+"\tmodel_index\tmodel\tnative_shared_keys\tnative_shortlist_rank\tfully_accepted\taccepted_start1\taccepted_stop1\tpartial\tscore\tscored_parent_links");
		writers[LINKS].println(event+"\tmodel_index\tparent\ttarget\tspan_parent_links\tuniquely_attributed");
		writers[RANKS].println(event+"\tparent\tmethod\twindow_length\tlocus_length\tcompatible_models\tunique_models\trank\ttie_first\ttie_last\tmodel_index\treason");
		writers[ATTEMPTS].println(event+"\tnative_rank\tmodel_index\tshared_keys\tinitial_identity_float_bits\taligned_length\tinitial_identity_pass");
		writers[ENDS].println(event+"\tactual_attempts");
		writers[PARENTS].println("input_slot\tinput_id\tfasta\tparent\tseqid\tstart1\tstop1\tstrand\tmask\ttarget\tprior\tmulti_supported\tscored_models\tmarginal_models");
		writers[LOCI].println("input_slot\tinput_id\tparent\tmethod\twindows\tcompatible_windows\tambiguous_only_windows\tunique_windows\tbest_rank\treason\tany_unique_model_native_shortlist_eligible\tany_uniquely_attributed_oracle_model_actually_attempted");
		writers[SUMMARY].println("method\tk\tbuckets\tscore\tall_loci\tunique_oracle_rankable\tcommon_native_rankable\tno_window\tno_compatible_model\tambiguous_only\tnative_zero_support\timproved\ttied\tworsened\tdistribution\tdistribution_loci\tmedian\tp90\tp95\tp99\tmaximum\trank_limit\tcovered_all\tcovered_common");
		writers[INPUTS].println("input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tactual_attempts\ttarget_loci\tseed_cpu_ns\tnative_index_cpu_ns\toracle_cpu_ns\tranking_cpu_ns\tworker_report_cpu_ns");
		writers[LOSSES].println("input_slot\tinput_id\tmethod\trank_limit\tbaseline_reachable_loci\tcandidate_reachable_loci\tinterpretation");
		writers[WORST].println("method\tinput_slot\tinput_id\tparent\tnative_rank\tcandidate_rank\trank_increase\tprior\tmulti_supported\tlineage");
		writers[WORK].println("item\tvalue");
		writers[GEOMETRY].println("geometry\tk\tbuckets\tprofile_construction_owner_cpu_ns\tquery_and_three_scores_owner_cpu_ns\tlocus_rank_queries_owner_cpu_ns\tprofile_key_payload_bytes\tprofile_support_payload_bytes");
	}

	private void report(){
		final ByteBuilder out=new ByteBuilder();
		for(int i=0; i<data.length; i++){
			final Input input=data[i];
			require(input.finished, "Every declared input must finish before complete rank summaries are published");
			for(final Parent p:input.truth.parents){
				out.clear().append(i).tab();
				field(out, input.truth.binding.id).tab();
				field(out, input.truth.binding.fasta).tab();
				field(out, p.row.id).tab();
				field(out, p.row.locus.seqid).tab().append(p.row.locus.start).tab().append(p.row.locus.stop)
					.tab().append(p.row.locus.strand).tab().append(p.row.mask).tab().append(p.target)
					.tab().append(p.prior).tab().append(p.multi).tab().append(p.scoredModels)
					.tab().append(p.marginalModels).nl();
				publish(PARENTS, out);
			}
			out.clear().append(i).tab();
			field(out, input.truth.binding.id).tab();
			field(out, input.truth.binding.fasta).tab().append(input.contigs).tab().append(input.bases)
				.tab().append(input.windows).tab().append(input.attempts).tab().append(input.loci.length)
				.tab().append(input.seedCpu).tab().append(input.nativeCpu).tab().append(input.oracleCpu)
				.tab().append(input.rankingCpu).tab().append(input.reportCpu).nl();
			publish(INPUTS, out);
			for(final LocusRanks locus:input.loci){
				for(int method=0; method<locus.best.length; method++){
					out.append(i).tab();
					field(out, input.truth.binding.id).tab();
					field(out, locus.parent.row.id).tab().append(method).tab().append(locus.windows).tab()
						.append(locus.compatibleWindows).tab().append(locus.ambiguousWindows).tab()
						.append(locus.uniqueWindows).tab().append(locus.best[method]).tab()
						.append(reason(locus, method, locus.compatibleWindows, locus.uniqueWindows)).tab()
						.append(locus.eligible).tab().append(locus.attempted).nl();
				}
				publish(LOCI, out);
			}
		}
		summaries(out);
		long nativeRankCpu=0;
		for(final Input input:data){nativeRankCpu+=input.nativeRankCpu;}
		for(int geometry=0; geometry<profiles.length; geometry++){
			long scoreCpu=0, rankCpu=0;
			for(final Input input:data){
				scoreCpu+=input.geometryScoreCpu[geometry];
				rankCpu+=input.geometryRankCpu[geometry];
			}
			final int width=widths[geometry%widths.length];
			out.append(geometry).tab().append(geometryK[geometry]).tab().append(width).tab()
				.append(profileGeometryCpu[geometry]).tab().append(scoreCpu).tab().append(rankCpu)
				.tab().append((long)width*models*Character.BYTES).tab().append((long)width*models*Integer.BYTES).nl();
			publish(GEOMETRY, out);
		}
		out.append("member_ingestion_owner_cpu_ns\t").append(memberCpu).nl()
			.append("profile_cpu_ns\t").append(profileCpu).nl()
			.append("native_positive_rank_owner_cpu_ns\t").append(nativeRankCpu).nl()
			.append("profile_key_payload_bytes\t").append(profileKeyBytes).nl()
			.append("profile_support_payload_bytes\t").append(profileSupportBytes).nl()
			.append("profile_heap_requirement\tUNMEASURED_NOT_PAYLOAD_BYTES\n")
			.append("native_index_payload_bytes\tUNMEASURED\n")
			.append("input_members\t").append(inputMembers).nl().append("kept_members\t").append(keptMembers)
			.nl().append("orphan_members\t").append(orphanMembers).nl()
			.append("coordinates_windows\tZERO_BASED_INCLUSIVE\n")
			.append("coordinates_parents_and_oracle\tONE_BASED_INCLUSIVE_ORIGINAL_FORWARD\n")
			.append("window_reaches_parent\tSAME_STRAND_ANY_PHYSICAL_OVERLAP\n")
			.append("compatible_span\tNATIVE_RECIPROCAL_HALF_OVERLAP\n")
			.append("locus_rank\tMINIMUM_ACROSS_CAPTURED_EVENTS_NOT_ONE_CALLER_WALK\n")
			.append("ordinary_attempt_reachability\tATTEMPTED_UNIQUE_ORACLE_ACCEPTABLE_MODEL_NOT_FINAL_NATIVE_ACCEPTANCE\n")
			.append("analysis_cpu_scope\tOWNER_THREADS_ONLY_NOT_NATIVE_WRITER_OR_COMPRESSION_CPU\n")
			.append("timing\tWHOLE_COMMAND_INCLUDES_ANALYSIS_NOT_BASELINE_CALLER_TIME\n")
			.append("status\tRANKING_DIAGNOSTIC_NOT_CALLER_ADOPTED\n");
		publish(WORK, out);
	}

	/** Histograms are over true loci; common comparisons require a native rank too. */
	private void summaries(final ByteBuilder out){
		final int methods=1+profiles.length*3;
		for(int method=0; method<methods; method++){
			final long[] histogram=new long[models+1], commonHistogram=new long[models+1];
			long total=0, rankable=0, common=0, noWindow=0, noModel=0, ambiguous=0, nativeZero=0;
			long improved=0, tied=0, worsened=0;
			for(int i=0; i<data.length; i++){
				final Input input=data[i];
				final long[] originalCoverage=new long[LIMITS.length], coverage=new long[LIMITS.length];
				for(final LocusRanks locus:input.loci){
					total++;
					if(locus.windows==0){noWindow++;}
					else if(locus.compatibleWindows==0){noModel++;}
					else if(locus.uniqueWindows==0){ambiguous++;}
					if(locus.uniqueWindows>0){rankable++;}
					final int baseline=locus.best[0], candidate=locus.best[method];
					if(locus.uniqueWindows>0 && baseline<0){nativeZero++;}
					if(candidate>0){histogram[candidate]++;}
					if(baseline>0 && candidate>0){
						common++;
						commonHistogram[candidate]++;
						if(candidate<baseline){improved++;}
						else if(candidate==baseline){tied++;}
						else{worsened++;}
					}
					for(int limit=0; limit<LIMITS.length; limit++){
						if(baseline>0 && baseline<=LIMITS[limit]){originalCoverage[limit]++;}
						if(candidate>0 && candidate<=LIMITS[limit]){coverage[limit]++;}
					}
				}
				for(int limit=0; limit<LIMITS.length; limit++){
					if(method==0 || originalCoverage[limit]==0 || coverage[limit]>0){continue;}
					out.append(i).tab();
					field(out, input.truth.binding.id).tab().append(method).tab().append(LIMITS[limit]).tab()
						.append(originalCoverage[limit]).append("\t0\tLOST_LAST_RANK_REACHABLE_LOCUS_NOT_DETECTION_LOSS\n");
					publish(LOSSES, out);
				}
			}
			final long[][] histograms={histogram, commonHistogram};
			final int[][] quantiles=new int[2][FRACTIONS.length];
			final long[] distributionCounts=new long[2];
			for(int distribution=0; distribution<2; distribution++){
				for(final long count:histograms[distribution]){distributionCounts[distribution]+=count;}
				for(int q=0; q<FRACTIONS.length; q++){
					quantiles[distribution][q]=quantile(histograms[distribution], FRACTIONS[q]);
				}
			}
			long covered=0, coveredCommon=0;
			int previous=0;
			for(final int limit:LIMITS){
				for(int rank=previous+1; rank<=Math.min(limit, models); rank++){
					covered+=histogram[rank];
					coveredCommon+=commonHistogram[rank];
				}
				previous=limit;
				for(int distribution=0; distribution<2; distribution++){
					out.append(method).tab();
					geometry(out, method).tab().append(total).tab().append(rankable).tab().append(common)
						.tab().append(noWindow).tab().append(noModel).tab().append(ambiguous).tab().append(nativeZero)
						.tab().append(improved).tab().append(tied).tab().append(worsened).tab()
						.append(distribution==0 ? "ALL_POSITIVELY_RANKED_LOCI" : "COMMON_NATIVE_RANKABLE_LOCI")
						.tab().append(distributionCounts[distribution]);
					for(final int q:quantiles[distribution]){out.tab().append(q);}
					out.tab().append(limit).tab().append(covered).tab().append(coveredCommon).nl();
					publish(SUMMARY, out);
				}
			}
			worst(method, out);
		}
	}

	/** Retains the ten largest regressions per method; ties use input/parent order. */
	private void worst(final int method, final ByteBuilder out){
		if(method==0){return;}
		final LocusRanks[] worst=new LocusRanks[10];
		final int[] slots=new int[10];
		for(int input=0; input<data.length; input++){
			for(final LocusRanks locus:data[input].loci){
				if(locus.best[0]<1 || locus.best[method]<=locus.best[0]){continue;}
				for(int place=0; place<worst.length; place++){
					if(worst[place]!=null && delta(worst[place], method)>=delta(locus, method)){continue;}
					System.arraycopy(worst, place, worst, place+1, worst.length-place-1);
					System.arraycopy(slots, place, slots, place+1, slots.length-place-1);
					worst[place]=locus;
					slots[place]=input;
					break;
				}
			}
		}
		for(int place=0; place<worst.length && worst[place]!=null; place++){
			final LocusRanks locus=worst[place];
			out.append(method).tab().append(slots[place]).tab();
			field(out, data[slots[place]].truth.binding.id).tab();
			field(out, locus.parent.row.id).tab().append(locus.best[0]).tab().append(locus.best[method])
				.tab().append(delta(locus, method)).tab().append(locus.parent.prior).tab()
				.append(locus.parent.multi).append("\tNOT_SUPPLIED_BY_REGISTERED_REFERENCE\n");
			publish(WORST, out);
		}
	}

	private static int delta(final LocusRanks locus, final int method){return locus.best[method]-locus.best[0];}

	/** Nearest-rank quantile over observed positive ranks; empty means -1. */
	public static int quantile(final long[] histogram, final double fraction){
		require(histogram!=null && histogram.length>1 && histogram[0]==0 && fraction>0 && fraction<=1,
			"Rank quantiles require a positive-rank histogram and fraction in (0,1]");
		long count=0;
		for(int rank=1; rank<histogram.length; rank++){
			require(histogram[rank]>=0, "A rank-frequency count cannot be negative");
			count+=histogram[rank];
		}
		if(count==0){return -1;}
		final long target=Math.max(1, (long)Math.ceil(fraction*count));
		long seen=0;
		for(int rank=1; rank<histogram.length; rank++){
			seen+=histogram[rank];
			if(seen>=target){return rank;}
		}
		throw new IllegalStateException("A positive rank histogram must contain its own nearest-rank quantile");
	}

	private ByteBuilder geometry(final ByteBuilder out, final int method){
		if(method==0){return out.append(family.indexK).append("\t0\tNATIVE_POSITIVE_SHARED_KEYS");}
		final int g=(method-1)/3;
		return out.append(ks[g/widths.length]).tab().append(widths[g%widths.length])
			.tab().append(SCORES[(method-1)%3]);
	}

	private static String reason(final LocusRanks locus, final int method, final long compatible, final long unique){
		if(locus.windows==0){return "NO_CAPTURED_WINDOW";}
		if(compatible==0){return "NO_COMPATIBLE_MODEL";}
		if(unique==0){return "AMBIGUOUS_SPANS_ONLY";}
		return method==0 && locus.best[0]<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED";
	}

	private static ByteBuilder field(final ByteBuilder out, final String text){return BuildNcrnaReference.field(out, text);}
	private static long cpu(){
		final long value=ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime();
		require(value>=0, "The accepted worker CPU clock must remain available during profiles/reports");
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained State       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Input{
		Input(final NcrnaDdlTruth.Input truth_, final int methods, final int models){
			truth=truth_;
			loci=new LocusRanks[truth.targets];
			geometryScoreCpu=new long[(methods-1)/3];
			geometryRankCpu=new long[(methods-1)/3];
			int i=0;
			for(final Parent parent:truth.parents){
				if(!parent.target){continue;}
				final LocusRanks locus=new LocusRanks(parent, methods, models);
				loci[i++]=locus;
				byParent.put(parent, locus);
			}
			assert(i==loci.length) : "Every qualified target parent must occupy one full-denominator row";
		}
		final NcrnaDdlTruth.Input truth;
		final LocusRanks[] loci;
		final long[] geometryScoreCpu, geometryRankCpu;
		final ObjectMap<Parent, LocusRanks> byParent=new ObjectMap<Parent, LocusRanks>(Parent.class, LocusRanks.class);
		long contigs, bases, windows, attempts, seedCpu, nativeCpu, oracleCpu, rankingCpu, reportCpu, nativeRankCpu;
		boolean finished;
	}

	/** One biological target's minima across explicit events, never across parents. */
	private static final class LocusRanks{
		LocusRanks(final Parent parent_, final int methods, final int models_){
			parent=parent_;
			models=models_;
			best=new int[methods];
			Arrays.fill(best, -1);
		}
		void add(final int method, final int rank){
			require(method>=0 && method<best.length && (rank==-1 || rank>0 && rank<=models),
				"A locus rank is missing (-1) or an actual one-based position in this model library");
			if(rank>0 && (best[method]<0 || rank<best[method])){best[method]=rank;}
		}
		final Parent parent;
		final int models;
		final int[] best;
		long windows, compatibleWindows, ambiguousWindows, uniqueWindows;
		boolean eligible, attempted;
	}

	private String callerConfig, reference, familymap, subtypes, witnesses, events;
	private String members, ledger, assignments, census, membercounts, prefix;
	private int models=508, threads=1, targetIndex=-1;
	private int[] ks={5, 7, 9, 11, 13, 15, 19}, widths={32, 64, 128, 256, 512, 1024, 2048}, geometryK;
	private NcrnaFamily family;
	private NcrnaDdlProfile[][] profiles;
	private Input[] data;
	private long memberCpu, profileCpu, profileKeyBytes, profileSupportBytes;
	private long[] profileGeometryCpu;
	private int inputMembers, keptMembers, orphanMembers;
	private final PrintStream outstream;
	private final ObjectIntMap<String> inputIndex=new ObjectIntMap<String>(String.class);
	private final fileIO.ByteStreamWriter[] writers=new fileIO.ByteStreamWriter[SUFFIX.length];
	private static final String TARGET="cobalamin_rf00174";
	private static final int WINDOWS=0, SEQUENCES=1, ORACLE=2, LINKS=3, RANKS=4, ATTEMPTS=5, ENDS=6, BUFFERS=7;
	private static final int PARENTS=7, LOCI=8, SUMMARY=9, INPUTS=10, LOSSES=11, WORST=12, WORK=13, GEOMETRY=14;
	private static final String[] SUFFIX={".windows.tsv.gz", ".windows.fa.gz", ".oracle.tsv.gz",
		".links.tsv.gz", ".ranks.tsv.gz", ".attempts.tsv.gz", ".events.tsv.gz", ".parents.tsv.gz",
		".loci.tsv.gz", ".summary.tsv", ".inputs.tsv", ".last_losses.tsv", ".worst_regressions.tsv", ".work.tsv", ".geometry_cpu.tsv"};
	private static final String[] SCORES={"Q1_SSU_FLOAT", "Q1_EQUAL_NONEMPTY_BUCKETS", "Q2_ALL_RAW_KEYS"};
	private static final int[] LIMITS={1, 2, 3, 5, 10, 20, 50, 100, 200, 384, 508};
	private static final double[] FRACTIONS={.5, .9, .95, .99, 1};
	private static final Locus[] NO_PARENTS=new Locus[0];
}
