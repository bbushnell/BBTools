package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Selected;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Two forward lengths in one member-voted profile, on frozen saved windows.
 * Adapted from A_SampleByteFileMT's owner prepass, native streams, parser and
 * worker lifecycle. Each worker owns complete geometries; immutable saved data
 * is shared. Primitive scores survive once for weighting without re-sketching.
 * Identity files close before any new rank is calculated. No genome, caller,
 * exhaustive verifier, resource replacement or automatic adoption.
	 * Optional identities=prefix restores accepted values/cohorts without alignment.
 * @author Brian, Raiden
 */
public final class RankDualNcrnaWindows{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, RankDualNcrnaWindows.class, false);
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("panels")){c.panels=b;}
				else if(a.equals("selection")){c.selection=b;}
				else if(a.equals("baseline")){c.baseline=b;}
				else if(a.equals("identities")){c.identities=b;}
				else if(a.equals("members")){c.members=b;}
				else if(a.equals("ledger")){c.ledger=b;}
				else if(a.equals("assignments")){c.assignments=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("summary")){c.summary=b;}
				else if(a.equals("ref")){c.references=b;}
				else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown saved dual-length parameter: "+arg);}
			}
			require(parser.out1!=null && parser.out2==null && parser.in1==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Declare complete saved bindings and out=fresh_prefix without genomes, sampling, overwrite or append");
			c.output=parser.out1;
			c.threads=Shared.threads();
			DynamicDemiLog.setExponent(4);
			run(c);
			pp.outstream.println("NCRNA_DUAL_SAVED_RANKS_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(exponent);
			Shared.closeStream(pp.outstream);
		}
	}

	public static void run(final Config c){
		require(c!=null && c.panels!=null && c.selection!=null && c.baseline!=null && c.members!=null
			&& c.ledger!=null && c.assignments!=null && c.census!=null && c.summary!=null
			&& c.references!=null && c.output!=null && c.models>0 && c.groups>0 && c.threads>0,
			"Saved dual ranks need all original cache/reference/member bindings and positive sizes");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Saved dual profiles use assertions, fixed native exponent4 and no blacklist");
		RrnaResourceIO.fresh(c.output, SUFFIXES);
		final int zip=ReadWrite.MAX_ZIP_THREADS();
		final boolean zipSet=ReadWrite.SET_ZIP_THREADS;
		try{
			ReadWrite.setZipThreads(1);
			Shared.TRIM_READ_DESCRIPTION=false;
			new RankDualNcrnaWindows(c).process();
		}finally{
			ReadWrite.setZipThreads(zip);
			ReadWrite.SET_ZIP_THREADS=zipSet;
		}
	}

	private RankDualNcrnaWindows(final Config c_){c=c_;}

	/*--------------------------------------------------------------*/
	/*----------------         Owner Prepass        ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> references=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			references.add(read.bases);
		});
		require(names.size()==c.models, "Original reference order must contain every cached model");
		labels=new byte[c.models][];
		for(int m=0; m<c.models; m++){labels[m]=names.get(m).getBytes(StandardCharsets.UTF_8);}
		final String[] headers=names.toArray(new String[names.size()]);
		final byte[][] bases=references.toArray(new byte[references.size()][]);
		members=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments, c.census, c.summary, headers, bases);
		data.readNativeCounts(headers);
		baseline();
		alphabet();
		cohorts=c.identities==null ? new NcrnaDdlCohorts(data, NcrnaDdlIdentity.measure(data, bases))
			: NcrnaDdlCohorts.read(c.identities, data, labels, bases);
		cohorts.write(c.output, data, labels);
		controls();
		for(int geometry=0; geometry<GEOMETRIES; geometry++){
			results[geometry]=new Geometry(this, geometry);
			results[geometry].prepare();
		}
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
		for(int i=0; i<Math.min(c.threads, GEOMETRIES); i++){workers.add(new ProcessThread(this));}
		ThreadWaiter.startAndWait(workers);
		if(failure.get()!=null){throw new IllegalStateException("Dual saved ranking failed before publication", failure.get());}
		for(final ProcessThread worker:workers){require(worker.success, "Every geometry worker must finish before reporting ranks");}
		report();
	}

	/** The mirrored DNA encoder supports ACGTN, not arbitrary IUPAC symbols.
	 * Validate once in the owner, never add a branch to each geometry's hot loop. */
	private void alphabet(){
		for(final Window window:data.windows){
			final int bad=invalidBase(window.bases);
			if(bad>=0){throw new IllegalArgumentException("Unsupported saved-window base: panel="+window.input.panel.name+", event="+window.key+", offset="+bad+", byte="+window.bases[bad]);}
		}
		for(int model=0; model<c.models; model++){
			for(final byte[] member:members.members(model)){
				final int bad=invalidBase(member);
				if(bad>=0){throw new IllegalArgumentException("Unsupported assigned-member base: model="+model+", offset="+bad+", byte="+member[bad]);}
			}
		}
	}

	/** Accept the corrected replay's full roster, not legacy canonical rows
	 * silently relabelled as fixed. Only three headline controls are recomputed. */
	private void baseline(){
		referenceRanks=readBaseline(c.baseline, data, c.models);
	}

	/** Binds the accepted corrected controls in original parent/method order.
	 * Shared by saved reporters so new designs cannot silently swap baselines. */
	public static int[][] readBaseline(final String path, final NcrnaDdlSavedData data, final int models){
		require(path!=null && data!=null && models>0, "Corrected controls need original saved data and model count");
		final int[][] result=new int[3][data.parents.size()];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			NcrnaDdlSavedData.header(in, BASELINE_HEADER);
			// CompareNcrnaHashBoundary.report writes dense parent/method blocks.
			// Bind their retained UTF-8 bytes directly, without per-field strings.
			for(final Parent parent:data.parents){
				for(int method=0; method<19; method++){
					final byte[] row=in.reader.nextLine();
					require(row!=null, "Corrected baseline cannot omit any parent or original method");
					p.set(row);
					require(p.terms()==10 && fieldEquals(p, 0, parent.input.panel.nameBytes) && p.parseInt(1)==parent.input.slot
						&& fieldEquals(p, 2, parent.input.idBytes) && fieldEquals(p, 3, parent.idBytes)
						&& p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 4) && p.termEquals(parent.reason, 9),
						"Corrected baseline rows must bind original panel/input/parent/stratum/disposition in producer order");
					require(p.termEquals(BASELINE_METHODS[method], 5), "Each corrected parent must retain all nineteen original methods exactly once in order");
					for(int state=6; state<9; state++){
						final int rank=p.parseInt(state);
						require(rank==-1 || rank>0 && rank<=models, "Corrected saved ranks must be positive or explicitly unrankable");
					}
					if(method<3){result[method][parent.index]=p.parseInt(8);}
					if(method==0){require(p.parseInt(8)==parent.nativeRank, "Native minima cannot change during a DDL ambiguity correction");}
				}
			}
			require(in.reader.nextLine()==null, "Corrected baseline cannot append duplicate parents or methods");
		}
		return result;
	}

	private void controls(){
		final DDLIndexSmall index=index(NcrnaDdlModelProfiles.build(members, 10, 2048, false, true));
		final NcrnaDdlQueryKeys query=new NcrnaDdlQueryKeys(2048, 10);
		final NcrnaDdlInterleave order=new NcrnaDdlInterleave(c.models);
		final float[] scores=new float[c.models];
		final int[] counts=new int[c.models];
		final boolean[] acceptable=new boolean[c.models];
		final int[] tuple=new int[4];
		final int[][] observed=new int[3][data.parents.size()];
		for(final int[] row:observed){Arrays.fill(row, -1);}
		for(final Window window:data.windows){
			query.build(window.bases, null, false, true);
			query.count(index, counts);
			for(int model=0; model<c.models; model++){scores[model]=(float)counts[model];}
			order.build(window.nativeCounts, scores);
			for(final Edge edge:window.edges){
				for(int model=0; model<c.models; model++){acceptable[model]=accepted(edge.uniqueModels, model);}
				order.nativeBest(acceptable, tuple);
				require(Arrays.equals(tuple, edge.nativeTuple), "Current native event tuple must reproduce its exact frozen rank/tie/model");
				observed[0][edge.parent.index]=NcrnaDdlSavedData.minimum(observed[0][edge.parent.index], tuple[0]);
				order.sketchBest(acceptable, tuple);
				observed[1][edge.parent.index]=NcrnaDdlSavedData.minimum(observed[1][edge.parent.index], tuple[0]);
				order.hybridBest(acceptable, tuple);
				observed[2][edge.parent.index]=NcrnaDdlSavedData.minimum(observed[2][edge.parent.index], tuple[0]);
			}
		}
		for(int method=0; method<3; method++){
			require(Arrays.equals(observed[method], referenceRanks[method]),
				"Native/full-width untagged masked10/hybrid control minima must exactly match the accepted corrected replay");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Geometry Work        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class ProcessThread extends Thread{
		ProcessThread(final RankDualNcrnaWindows owner_){owner=owner_;}
		@Override
		public void run(){
			try{
				for(int geometry=owner.next.getAndIncrement(); geometry<GEOMETRIES; geometry=owner.next.getAndIncrement()){
					if(owner.failure.get()!=null){return;}
					owner.results[geometry].process();
				}
				success=true;
			}catch(RuntimeException | Error invalid){owner.failure.compareAndSet(null, invalid);}
		}
		private final RankDualNcrnaWindows owner;
		boolean success;
	}

	private static final class Geometry{
		Geometry(final RankDualNcrnaWindows owner_, final int id_){
			owner=owner_;
			id=id_;
			partitioned=(id&1)==1;
			width=((id>>>1)&1)==0 ? 2048 : 4096;
			longK=((id>>>2)&1)==0 ? 10 : 11;
			shortK=((id>>>3)&1)==0 ? 7 : 9;
			final long cells=(long)owner.data.windows.size()*owner.c.models;
			require(cells<=Integer.MAX_VALUE-8, "Saved half-count arrays must fit one primitive Java array; split an oversized population explicitly");
			longs=new char[(int)cells];
			shorts=partitioned ? new char[(int)cells] : null;
			minima=new int[partitioned ? 7 : 1][owner.data.parents.size()];
			for(final int[] row:minima){Arrays.fill(row, -1);}
			weights=partitioned ? new double[]{0,.25,.5,.75,1,Double.NaN,Double.NaN} : new double[]{0};
			if(partitioned){
				locusMoments=new double[owner.data.parents.size()*8];
				paired=new int[owner.data.parents.size()];
				missingSignal=new int[paired.length];
				missingBackground=new int[paired.length];
				observations=new long[paired.length*2];
				for(int direction=0; direction<2; direction++){
					for(int cell=0; cell<4; cell++){fit[direction][cell]=new NcrnaDdlMoments();}
				}
			}
		}

		/** Native member construction stays in the owner's prepass, as required
		 * by NcrnaDdlModelProfiles; workers only query immutable postings. */
		void prepare(){
			final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.buildDual(owner.members, shortK, longK, width, partitioned);
			profileIndex=index(profiles);
			if(partitioned){
				final NcrnaDdlProfile[] half=NcrnaDdlModelProfiles.build(owner.members, longK, width/2, false, true);
				for(int model=0; model<profiles.length; model++){
					for(int bucket=0; bucket<width/2; bucket++){
						require(profiles[model].keys()[bucket+width/2]==half[model].keys()[bucket],
							"Partition long-half member voting must match independent untagged half-width long-only construction");
					}
				}
				controlIndex=index(half);
			}
		}

		void process(){
			assert(profileIndex!=null) : "The owner must complete immutable native member postings before workers start";
			final NcrnaDdlDualKeys query=new NcrnaDdlDualKeys(shortK, longK, width, partitioned);
			final NcrnaDdlQueryKeys controlQuery=partitioned ? new NcrnaDdlQueryKeys(width/2, longK) : null;
			final int[] shortScores=new int[owner.c.models], longScores=new int[owner.c.models], check=new int[owner.c.models];
			final int[][] halves={shortScores, longScores};
			for(int w=0; w<owner.data.windows.size(); w++){
				final Window window=owner.data.windows.get(w);
				boolean useful=false;
				for(final Edge edge:window.edges){useful|=edge.unique>0;}
				if(useful){
					query.build(window.bases);
					if(partitioned){
						query.countHalves(profileIndex, halves);
						controlQuery.build(window.bases, null, false, true);
						controlQuery.count(controlIndex, check);
						require(Arrays.equals(longScores, check), "Partition weight0 counts must match the independent half-width long-only control on every scored window");
					}else{query.count(profileIndex, longScores);}
					for(int model=0, offset=w*owner.c.models; model<owner.c.models; model++){
						assert(longScores[model]>=0 && longScores[model]<=width
							&& (!partitioned || shortScores[model]>=0 && shortScores[model]<=width/2))
							: "DDLIndexSmall counts one match per stored bucket; saved counts must fit their actual bin areas: model="+model;
						longs[offset+model]=(char)longScores[model];
						if(partitioned){shorts[offset+model]=(char)shortScores[model];}
					}
					scored++;
				}
				if(partitioned){events(window, w*owner.c.models);}
			}
			if(partitioned){fit();}
			for(int w=0; w<owner.data.windows.size(); w++){
				for(final Edge edge:owner.data.windows.get(w).edges){
					if(edge.unique==0){continue;}
					for(int method=0; method<weights.length; method++){
						if(!finite(weights[method])){continue;}
						final int rank=best(shorts, longs, w*owner.c.models, owner.c.models, weights[method], edge.uniqueModels);
						minima[method][edge.parent.index]=NcrnaDdlSavedData.minimum(minima[method][edge.parent.index], rank);
					}
				}
			}
		}

		private void events(final Window window, final int offset){
			for(final Edge edge:window.edges){
				final int parent=edge.parent.index;
				if(edge.unique==0){
					missingSignal[parent]++;
					continue;
				}
				if(edge.unique==owner.c.models){
					missingBackground[parent]++;
					continue;
				}
				final int signal=normalizedMoments(shorts, longs, offset, owner.c.models, edge.uniqueModels, eventMoments);
				require(signal==edge.unique, "Model-level signal weights must count the exact original unique-acceptance bits");
				for(int value=0; value<8; value++){locusMoments[parent*8+value]+=eventMoments[value];}
				paired[parent]++;
				observations[parent*2]+=signal;
				observations[parent*2+1]+=owner.c.models-signal;
			}
		}

		private void fit(){
			for(final Parent parent:owner.data.parents){
				if(paired[parent.index]==0){continue;}
				for(int direction=0; direction<2; direction++){
					if(!owner.cohorts.fit(parent, direction)){continue;}
					for(int cell=0; cell<4; cell++){add(fit[direction][cell], parent.index, cell);}
				}
			}
			for(int direction=0; direction<2; direction++){
				final NcrnaDdlMoments[] f=fit[direction];
				weights[5+direction]=empiricalWeight(f[0].mean(), f[1].mean(), f[1].standardDeviation(), f[2].mean(), f[3].mean(), f[3].standardDeviation());
				reasons[direction]=f[0].loci()==0 ? "NO_PAIRED_LOW_IDENTITY_FIT_LOCUS"
					: !finite(f[1].standardDeviation()) || f[1].standardDeviation()==0 ? "SHORT_BACKGROUND_SD_UNAVAILABLE"
					: !finite(f[3].standardDeviation()) || f[3].standardDeviation()==0 ? "LONG_BACKGROUND_SD_UNAVAILABLE"
					: (f[2].mean()-f[3].mean())/f[3].standardDeviation()<=0 ? "LONG_SEPARATION_NONPOSITIVE"
					: !finite(weights[5+direction]) ? "NONFINITE_WEIGHT" : weights[5+direction]==0 ? "SHORT_SEPARATION_NONPOSITIVE" : "MEASURED";
			}
		}

		void add(final NcrnaDdlMoments out, final int parent, final int cell){
			out.addLocus(locusMoments[parent*8+cell*2], locusMoments[parent*8+cell*2+1], paired[parent], observations[parent*2+(cell&1)]);
		}
		final RankDualNcrnaWindows owner;
		final int id, shortK, longK, width;
		final boolean partitioned;
		final char[] shorts, longs;
		final int[][] minima;
		final double[] weights;
		double[] locusMoments;
		int[] paired, missingSignal, missingBackground;
		long[] observations;
		final NcrnaDdlMoments[][] fit=new NcrnaDdlMoments[2][4];
		final double[] eventMoments=new double[8];
		final String[] reasons=new String[2];
		private DDLIndexSmall profileIndex, controlIndex;
		long scored;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reporting          ----------------*/
	/*--------------------------------------------------------------*/

	private void report(){
		final ByteBuilder b=new ByteBuilder();
		final int[][] reachable=new int[3][data.selected.size()];
		final int[] after=new int[data.selected.size()];
		try(Output models=new Output(c.output+".models.tsv"); Output inputs=new Output(c.output+".inputs.tsv");
				Output methods=new Output(c.output+".methods.tsv"); Output loci=new Output(c.output+".loci.tsv.gz");
				Output summary=new Output(c.output+".summary.tsv"); Output losses=new Output(c.output+".last_losses.tsv");
				Output halves=new Output(c.output+".halves.tsv.gz"); Output moments=new Output(c.output+".moments.tsv");
				Output statistics=new Output(c.output+".statistics_loci.tsv.gz"); Output weights=new Output(c.output+".weights.tsv");
				Output work=new Output(c.output+".work.tsv")){
			models.writer.println("model_index\tmodel");
			for(int model=0; model<c.models; model++){models.writer.print(b.clear().append(model).tab().append(labels[model]).nl());}
			inputs.writer.println("panel\tinput_slot\tinput_id\tfasta\tassembly_group\ttargets\tnative_common\tnative_zero\tno_model\tno_window\tambiguous\tcaptured_windows");
			for(final Selected input:data.selected){
				inputs.writer.print(b.clear().append(input.panel.nameBytes).tab().append(input.slot).tab().append(input.idBytes).tab().append(input.pathBytes)
					.tab().append(input.group).tab().append(input.parents).tab().append(input.common).tab().append(input.zero)
					.tab().append(input.gaps).tab().append(input.noWindow).tab().append(input.ambiguous).tab().append(input.windows).nl());
			}
			methods.writer.println("method\tgeometry\tmode\tshort_k\tmasked_long_k\ttotal_buckets\tweight\tweight_double_bits\tfit_direction\tstatus");
			loci.writer.println("method\tpanel\tinput_slot\tinput_id\tparent\tseen_stratum\tidentity_bin\tgroup_side\tevaluation_role\tnative_common\tmethod_available\tbest_rank\tsource_disposition");
			summary.writer.println("method\tseen_stratum\tidentity_bin\tevaluation_role\ttarget_loci\tcommon_loci\tavailable_common\trank_sum\twithin10\twithin50\twithin100\tnative_rank_sum\tnative_within10\tnative_within50\tnative_within100\tmasked10_rank_sum\tmasked10_within10\tmasked10_within50\tmasked10_within100\thybrid_rank_sum\thybrid_within10\thybrid_within50\thybrid_within100\tnative_reachable_inputs\tnative_last_losses\tmasked10_reachable_inputs\tmasked10_last_losses\thybrid_reachable_inputs\thybrid_last_losses\tinterpretation");
			losses.writer.println("method\tseen_stratum\tidentity_bin\tevaluation_role\treference\tpanel\tinput_slot\tinput_id\tfasta\treference_reachable_loci\tinterpretation");
			halves.writer.println("geometry\tpanel\tinput_slot\tinput_id\tevent\tshort_counts_by_original_model\tlong_counts_by_original_model");
			moments.writer.println("geometry\tpopulation\tseen_stratum\tidentity_bin\tfit_direction\thalf\tclass\tlocus_weights\tevents\tmodel_observations\tmean\tpopulation_sd");
			statistics.writer.println("geometry\tpanel\tinput_slot\tinput_id\tparent\tseen_stratum\tidentity_bin\tpaired_events\tempty_signal_events\tempty_background_events\tsignal_observations\tbackground_observations\tfit_A\tfit_B\tshort_signal_mean\tshort_signal_second\tshort_background_mean\tshort_background_second\tlong_signal_mean\tlong_signal_second\tlong_background_mean\tlong_background_second");
			weights.writer.println("geometry\tfit_direction\tfit_loci\tshort_signal_mean\tshort_background_mean\tshort_background_sd\tlong_signal_mean\tlong_background_mean\tlong_background_sd\tD_short\tD_long\tweight\tweight_double_bits\tstatus");
			work.writer.println("geometry\tmode\tselected_inputs\ttarget_loci\tcaptured_windows\tscored_windows\tsaved_count_cells\tlong_half_control_windows\tinterpretation");
			for(int reference=0; reference<3; reference++){
				writeMethod(new Method(REFERENCES[reference], null, referenceRanks[reference], 0, -1, "ACCEPTED_CORRECTED_BASELINE"), methods, loci, summary, losses, reachable, after, b);
			}
			for(final Geometry g:results){
				require(g!=null, "Every one of the sixteen declared geometries must be complete before report publication");
				if(g.partitioned){writeStatistics(g, halves, moments, statistics, weights, b);}
				for(int method=0; method<g.weights.length; method++){
					final int direction=method>=5 ? method-5 : -1;
					final String name="g"+g.id+(g.partitioned ? direction<0 ? "_weight_"+g.weights[method] : "_"+DIRECTIONS[direction] : "_COMPETE");
					writeMethod(new Method(name, g, g.minima[method], g.weights[method], direction,
						direction<0 ? "FIXED" : g.reasons[direction]), methods, loci, summary, losses, reachable, after, b);
				}
				work.writer.print(b.clear().append(g.id).tab().append(g.partitioned ? "PARTITION" : "COMPETE").tab().append(data.selected.size())
					.tab().append(data.parents.size()).tab().append(data.windows.size()).tab().append(g.scored)
					.tab().append((long)g.longs.length+(g.shorts==null ? 0L : g.shorts.length)).tab().append(g.partitioned ? g.scored : 0)
					.append("\tSAVED_ORDER_AND_DIAGNOSTIC_CONTROLS_NOT_CALLER_COST\n"));
			}
		}
	}

	private void writeMethod(final Method method, final Output methods, final Output loci, final Output summary,
			final Output losses, final int[][] reachable, final int[] after, final ByteBuilder b){
		final Geometry g=method.geometry;
		b.clear().append(method.name).tab().append(g==null ? -1 : g.id).tab().append(g==null ? "REFERENCE" : g.partitioned ? "PARTITION" : "COMPETE")
			.tab().append(g==null ? -1 : g.shortK).tab().append(g==null ? -1 : g.longK).tab().append(g==null ? -1 : g.width).tab();
		if(g==null || !g.partitioned){b.append("NA");}else{number(b, method.weight);}
		methods.writer.print(b.tab().append(Double.doubleToRawLongBits(method.weight)).tab()
			.append(method.direction<0 ? "NONE" : DIRECTIONS[method.direction]).tab().append(method.status).nl());
		for(final Parent parent:data.parents){
			final int p=parent.index;
			loci.writer.print(b.clear().append(method.name).tab().append(parent.input.panel.nameBytes).tab().append(parent.input.slot)
				.tab().append(parent.input.idBytes).tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum])
				.tab().append(NcrnaDdlCohorts.BINS[cohorts.bins[p]]).tab().append(cohorts.sides[p])
				.tab().append(NcrnaDdlCohorts.ROLES[cohorts.role(parent, method.direction)]).tab().append(parent.nativeRank>0)
				.tab().append(method.available).tab().append(method.ranks[p]).tab().append(parent.reason).nl());
		}
		for(int role=method.direction<0 ? 0 : 1; role<(method.direction<0 ? 1 : NcrnaDdlCohorts.ROLES.length); role++){
			for(int state=0; state<3; state++){
				for(int bin=0; bin<6; bin++){summarize(method, state, bin, role, summary, losses, reachable, after, b);}
			}
		}
	}

	private void summarize(final Method method, final int state, final int bin, final int role,
			final Output summary, final Output losses, final int[][] reachable, final int[] after, final ByteBuilder b){
		for(final int[] counts:reachable){Arrays.fill(counts, 0);}
		Arrays.fill(baselineSums, 0);
		for(final int[] counts:baselineCovered){Arrays.fill(counts, 0);}
		Arrays.fill(after, 0);
		int targets=0, common=0, within10=0, within50=0, within100=0;
		long sum=0;
		for(final Parent parent:data.parents){
			final int p=parent.index;
			if(parent.stratum!=state || bin<5 && cohorts.bins[p]!=bin || cohorts.role(parent, method.direction)!=role){continue;}
			targets++;
			if(parent.nativeRank<=0){continue;}
			common++;
			for(int reference=0; reference<3; reference++){
				final int rank=referenceRanks[reference][p];
				require(rank>0, "Matched reference cohorts retain every original native-common locus");
				baselineSums[reference]+=rank;
				if(rank<=10){baselineCovered[reference][0]++;}
				if(rank<=50){
					baselineCovered[reference][1]++;
					reachable[reference][parent.input.index]++;
				}
				if(rank<=100){baselineCovered[reference][2]++;}
			}
			if(!method.available){continue;}
			final int rank=method.ranks[p];
			require(rank>0 && rank<=c.models, "Every available method must completely rank every original native-common locus");
			sum+=rank;
			if(rank<=10){within10++;}
			if(rank<=50){
				within50++;
				after[parent.input.index]++;
			}
			if(rank<=100){within100++;}
		}
		b.clear().append(method.name).tab().append(NcrnaDdlSavedData.STRATA[state]).tab().append(NcrnaDdlCohorts.BINS[bin])
			.tab().append(NcrnaDdlCohorts.ROLES[role]).tab().append(targets).tab().append(common).tab().append(method.available ? common : 0).tab();
		if(method.available){b.append(sum).tab().append(within10).tab().append(within50).tab().append(within100);}
		else{b.append("NA\tNA\tNA\tNA");}
		for(int reference=0; reference<3; reference++){
			b.tab().append(baselineSums[reference]);
			for(final int count:baselineCovered[reference]){b.tab().append(count);}
		}
		// Complete each summary row before reusing b for detailed loss rows.
		for(int reference=0; reference<3; reference++){
			int eligible=0, lost=0;
			for(final Selected input:data.selected){
				if(reachable[reference][input.index]>0){
					eligible++;
					if(after[input.index]==0){lost++;}
				}
			}
			b.tab().append(eligible).tab();
			if(method.available){b.append(lost);}else{b.append("NA");}
		}
		summary.writer.print(b.append("\tSAVED_NATIVE_COMMON_LOCI_NOT_CALLER_DETECTION\n"));
		if(!method.available){return;}
		for(int reference=0; reference<3; reference++){
			for(final Selected input:data.selected){
				if(reachable[reference][input.index]==0 || after[input.index]>0){continue;}
				losses.writer.print(b.clear().append(method.name).tab().append(NcrnaDdlSavedData.STRATA[state]).tab().append(NcrnaDdlCohorts.BINS[bin])
					.tab().append(NcrnaDdlCohorts.ROLES[role]).tab().append(REFERENCES[reference]).tab().append(input.panel.nameBytes)
					.tab().append(input.slot).tab().append(input.idBytes).tab().append(input.pathBytes).tab().append(reachable[reference][input.index])
					.append("\tLAST_SAVED_REACHABLE_LOCUS_LOST_AT50_NOT_CALLER_DETECTION\n"));
			}
		}
	}

	private void writeStatistics(final Geometry g, final Output halves, final Output moments, final Output statistics,
			final Output weights, final ByteBuilder b){
		for(int w=0; w<data.windows.size(); w++){
			final Window window=data.windows.get(w);
			boolean useful=false;
			for(final Edge edge:window.edges){useful|=edge.unique>0;}
			if(!useful){continue;}
			b.clear().append(g.id).tab().append(window.input.panel.nameBytes).tab().append(window.input.slot).tab().append(window.input.idBytes).tab().append(window.event);
			for(int half=0; half<2; half++){
				final char[] scores=half==0 ? g.shorts : g.longs;
				b.tab();
				for(int model=0; model<c.models; model++){
					if(model>0){b.append(',');}
					b.append((int)scores[w*c.models+model]);
				}
			}
			halves.writer.print(b.nl());
		}
		final NcrnaDdlMoments[][] cells=new NcrnaDdlMoments[15][4];
		for(int cell=0; cell<15; cell++){
			for(int kind=0; kind<4; kind++){cells[cell][kind]=new NcrnaDdlMoments();}
		}
		for(final Parent parent:data.parents){
			final int p=parent.index;
			b.clear().append(g.id).tab().append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab().append(parent.input.idBytes)
				.tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum]).tab().append(NcrnaDdlCohorts.BINS[cohorts.bins[p]])
				.tab().append(g.paired[p]).tab().append(g.missingSignal[p]).tab().append(g.missingBackground[p]).tab().append(g.observations[p*2])
				.tab().append(g.observations[p*2+1]).tab().append(g.paired[p]>0 && cohorts.fit(parent, 0)).tab().append(g.paired[p]>0 && cohorts.fit(parent, 1));
			for(int value=0; value<8; value++){
				b.tab();
				number(b, g.paired[p]>0 ? g.locusMoments[p*8+value]/g.paired[p] : Double.NaN);
			}
			statistics.writer.print(b.nl());
			if(g.paired[p]>0){
				for(int kind=0; kind<4; kind++){g.add(cells[parent.stratum*5+cohorts.bins[p]][kind], p, kind);}
			}
		}
		for(int cell=0; cell<15; cell++){
			for(int kind=0; kind<4; kind++){momentRow(g.id, "DESCRIPTIVE", NcrnaDdlSavedData.STRATA[cell/5], NcrnaDdlCohorts.BINS[cell%5], "NONE", kind, cells[cell][kind], moments, b);}
		}
		for(int direction=0; direction<2; direction++){
			final NcrnaDdlMoments[] f=g.fit[direction];
			for(int kind=0; kind<4; kind++){momentRow(g.id, "FITTING_ONLY", "SEEN_PLUS_NOT_SEEN_FIT_SIDE", "Q1_LOW", DIRECTIONS[direction], kind, f[kind], moments, b);}
			b.clear().append(g.id).tab().append(DIRECTIONS[direction]).tab().append(f[0].loci());
			for(final double value:new double[]{f[0].mean(),f[1].mean(),f[1].standardDeviation(),f[2].mean(),f[3].mean(),f[3].standardDeviation(),
					(f[0].mean()-f[1].mean())/f[1].standardDeviation(),(f[2].mean()-f[3].mean())/f[3].standardDeviation(),g.weights[5+direction]}){
				b.tab();
				number(b, value);
			}
			weights.writer.print(b.tab().append(Double.doubleToRawLongBits(g.weights[5+direction])).tab().append(g.reasons[direction]).nl());
		}
	}

	private static void momentRow(final int geometry, final String population, final String state, final String bin,
			final String direction, final int kind, final NcrnaDdlMoments value, final Output out, final ByteBuilder b){
		b.clear().append(geometry).tab().append(population).tab().append(state).tab().append(bin).tab().append(direction)
			.tab().append(kind<2 ? "SHORT" : "LONG").tab().append((kind&1)==0 ? "SIGNAL" : "BACKGROUND")
			.tab().append(value.loci()).tab().append(value.events()).tab().append(value.models()).tab();
		number(b, value.mean());
		b.tab();
		number(b, value.standardDeviation());
		out.writer.print(b.nl());
	}

	private static final class Method{
		Method(final String name_, final Geometry geometry_, final int[] ranks_, final double weight_, final int direction_, final String status_){
			name=name_;
			geometry=geometry_;
			ranks=ranks_;
			weight=weight_;
			direction=direction_;
			status=status_;
			available=finite(weight);
		}
		final String name, status;
		final Geometry geometry;
		final int[] ranks;
		final double weight;
		final int direction;
		final boolean available;
	}

	private static void number(final ByteBuilder b, final double value){
		if(finite(value)){b.appendSci(value, 16);}else{b.append("NA");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Deterministic model-order ties, preserving double empirical weights. */
	public static int best(final char[] shorts, final char[] longs, final int offset,
			final int models, final double weight, final long[] acceptable){
		require(longs!=null && offset>=0 && models>0 && (long)offset+models<=longs.length
			&& (shorts==null || shorts.length==longs.length) && finite(weight) && weight>=0
			&& acceptable!=null && acceptable.length==(models+63)/64, "Weighted ranks need aligned primitive scores, original model bits and a finite nonnegative weight");
		int selected=-1;
		double target=-1;
		for(int model=0; model<models; model++){
			final double value=longs[offset+model]+(shorts==null ? 0 : weight*shorts[offset+model]);
			require(finite(value), "A finite empirical weight must produce finite ranking scores; no silent clipping");
			if(accepted(acceptable, model) && value>target){
				selected=model;
				target=value;
			}
		}
		if(selected<0){return -1;}
		int rank=1;
		for(int model=0; model<models; model++){
			final double value=longs[offset+model]+(shorts==null ? 0 : weight*shorts[offset+model]);
			if(value>target || value==target && model<selected){rank++;}
		}
		return rank;
	}

	public static double empiricalWeight(final double signalShort, final double backgroundShort, final double sdShort,
			final double signalLong, final double backgroundLong, final double sdLong){
		if(!finite(signalShort) || !finite(backgroundShort) || !finite(sdShort) || sdShort<=0
			|| !finite(signalLong) || !finite(backgroundLong) || !finite(sdLong) || sdLong<=0){return Double.NaN;}
		final double ds=(signalShort-backgroundShort)/sdShort, dl=(signalLong-backgroundLong)/sdLong;
		if(!finite(ds) || !finite(dl) || dl<=0){return Double.NaN;}
		final double weight=ds<=0 ? 0 : ds/dl;
		return finite(weight) ? weight : Double.NaN;
	}

	public static int invalidBase(final byte[] bases){
		require(bases!=null, "DNA alphabet validation requires actual saved sequence bytes");
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b!='A' && b!='C' && b!='G' && b!='T' && b!='N' && b!='a' && b!='c' && b!='g' && b!='t' && b!='n'){return i;}
		}
		return -1;
	}

	/** One model pass computes the four class/half first and second moments.
	 * Every class has total weight one within this event. Caller records empty
	 * classes separately rather than passing them as observations of zero. */
	public static int normalizedMoments(final char[] shorts, final char[] longs, final int offset,
			final int models, final long[] acceptable, final double[] out){
		require(shorts!=null && longs!=null && shorts.length==longs.length && offset>=0 && models>1
			&& (long)offset+models<=shorts.length && acceptable!=null && acceptable.length==(models+63)/64
			&& out!=null && out.length==8, "Event moments need both half-count vectors, original model bits and eight reused output values");
		Arrays.fill(out, 0);
		int signal=0;
		for(int model=0; model<models; model++){
			final int type=accepted(acceptable, model) ? 0 : 1;
			if(type==0){signal++;}
			final double s=shorts[offset+model], l=longs[offset+model];
			out[type*2]+=s;
			out[type*2+1]+=s*s;
			out[(type+2)*2]+=l;
			out[(type+2)*2+1]+=l*l;
		}
		require(signal>0 && signal<models, "Paired event moments cannot replace an empty signal or background class with zero");
		for(int cell=0; cell<4; cell++){
			final int count=(cell&1)==0 ? signal : models-signal;
			out[cell*2]/=count;
			out[cell*2+1]/=count;
		}
		return signal;
	}

	private static DDLIndexSmall index(final NcrnaDdlProfile[] profiles){
		final char[][] keys=new char[profiles.length][];
		for(int model=0; model<keys.length; model++){keys[model]=profiles[model].keys();}
		return new DDLIndexSmall(keys);
	}
	private static boolean accepted(final long[] bits, final int model){return (bits[model>>>6]&(1L<<(model&63)))!=0;}
	private static boolean finite(final double value){return !Double.isNaN(value) && !Double.isInfinite(value);}
	private static boolean fieldEquals(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String panels, selection, baseline, identities, members, ledger, assignments, census, summary, references, output;
		public int models=508, groups=3153, threads=16;
	}
	private final Config c;
	private NcrnaDdlSavedData data;
	private NcrnaDdlMembers members;
	private NcrnaDdlCohorts cohorts;
	private byte[][] labels;
	private int[][] referenceRanks;
	private final long[] baselineSums=new long[3];
	private final int[][] baselineCovered=new int[3][3];
	private final Geometry[] results=new Geometry[GEOMETRIES];
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final int GEOMETRIES=16;
	private static final String[] REFERENCES={"native", "masked10_raw", "hybrid"};
	private static final String[] BASELINE_METHODS=new String[19];
	static{
		System.arraycopy(REFERENCES, 0, BASELINE_METHODS, 0, REFERENCES.length);
		for(int i=3; i<BASELINE_METHODS.length; i++){BASELINE_METHODS[i]="canonical_control_"+NcrnaDdlSavedData.CONTROL_METHODS[i-3];}
	}
	private static final String[] DIRECTIONS={"FIT_A_JUDGE_B", "FIT_B_JUDGE_A"};
	private static final String BASELINE_HEADER="panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\told_rank\tquery_fixed_rank\tall_fixed_rank\tsource_disposition";
	private static final String[] SUFFIXES={".identity.tsv", ".identity_bins.tsv", ".groups.tsv", ".models.tsv", ".methods.tsv", ".loci.tsv.gz", ".summary.tsv", ".last_losses.tsv", ".halves.tsv.gz", ".moments.tsv", ".statistics_loci.tsv.gz", ".weights.tsv", ".inputs.tsv", ".work.tsv"};
}
