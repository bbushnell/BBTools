package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectSet;
import map.IntLongHashMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Exact finite magnitude-threshold replay on cached, matched real calls.
 * Learns model medians without the held assembly fold; fold4 alone selects
 * thresholds. Outer rows never teach offsets, candidate breakpoints or selection.
 * No alignment, fitting, matching, caller execution or path reconstruction.
 * @author Brian, Ganyu
 */
public final class NcrnaOffsetThresholdReplay{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaOffsetThresholdReplay.class, false);
		try{
			final NcrnaOffsetThresholdReplay tool=new NcrnaOffsetThresholdReplay(pp.args);
			tool.process();
			pp.outstream.println(tool.residuals ? "NCRNA_FOLD_RESIDUALS_COMPLETE_TRAINING_ONLY_NO_OUTER_READ"
				: "NCRNA_OFFSET_THRESHOLD_REPLAY_COMPLETE_CONDITIONAL_NO_CALLER_NO_FIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaOffsetThresholdReplay(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){
				require("threshold".equals(b) || "residuals".equals(b), "mode=threshold or residuals is required");
				residuals="residuals".equals(b);
			}else if(a.equals("offsets")){offsetDirectory=b;}
			else if(a.equals("folds")){folds=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("member")){member=b;}
			else{require(parser.parse(arg, a, b), "Unknown offset-threshold argument: "+arg);}
		}
		manifest=parser.in1;
		output=parser.out1;
		require(manifest!=null && folds!=null && reference!=null && output!=null
			&& ("RF00174".equals(member) || "RF01482".equals(member)),
			"Require in=geometry_manifest folds=complete_assemblies.tsv ref=frozen_library member=RF00174|RF01482 out=fresh_directory");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Replay complete fixed cohorts without sampling or replacement");
		require(Tools.testInputFiles(false, true, manifest, folds, reference) && !new File(output).exists(),
			"Inputs must exist and output directory must be fresh");
		require(residuals ? offsetDirectory!=null && new File(offsetDirectory).isDirectory() : offsetDirectory==null,
			"mode=residuals requires offsets=accepted_threshold_result_directory; ordinary threshold mode never imports offsets");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		readFolds();
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(reference);
		for(final String name:library.names){
			require(!modelIndex.contains(name), "Frozen model identifiers must be unique");
			modelIndex.put(name, models.size());
			models.add(new Model(name));
		}
		require(!models.isEmpty(), "Offset learning needs the complete nonempty model axis");
		readDatasets();
		if(residuals){
			readOffsets();
			require(new File(output).mkdir(), "Residual histogram output must be fresh");
			writeResiduals();
			return;
		}
		for(final Model model:models){model.learn();}
		final IntList values=new IntList();
		values.add(0);
		values.add(Integer.MAX_VALUE);// Explicit never-shift control; all valid coordinate errors are bounded by this.
		for(final Model model:models){
			for(int fold=0; fold<5; fold++){
				if(model.support[fold]<SUPPORT){continue;}
				values.add(Math.abs(model.median[fold][0]));
				values.add(Math.abs(model.median[fold][1]));
			}
		}
		values.sort();
		int unique=0;
		for(int i=0; i<values.size; i++){
			if(unique==0 || values.array[i]!=values.array[unique-1]){values.array[unique++]=values.array[i];}
		}
		final int[] thresholds=Arrays.copyOf(values.array, unique);
		final boolean[] selected=new boolean[unique];
		Stats[] best=null;
		for(int t=0; t<unique; t++){
			final Stats[] trial=score(training, 4, 4, thresholds[t], 2);
			final int comparison=best==null ? 1 : compare(trial, best);
			if(comparison>0){
				Arrays.fill(selected, false);
				best=trial;
			}
			if(comparison>=0){selected[t]=true;}
		}
		require(best!=null && best[0].n>0, "Fold4 must have matched calls before selecting a threshold");
		require(new File(output).mkdir(), "Cannot create fresh threshold-result directory");
		writeOffsets();
		final ByteBuilder row=new ByteBuilder();
		try(Output curve=new Output(path("internal.tsv")); Output outer=new Output(path("outer.tsv"));
				Output decisions=new Output(path("outer_decisions.tsv"));
				Output choices=new Output(path("thresholds.tsv")); Output population=new Output(path("population.tsv"))){
			curve.writer.println(RESULT_HEADER);
			outer.writer.println(RESULT_HEADER);
			decisions.writer.println("member\tpanel\tthreshold_bp\tselected_on_fold4\tscope\tinput_id\tassembly\treference_id\tmodel\tseqid\tstrand\tcontig_length\tfixed_mask\tsupport\tmedian5\tmedian3\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\tgeometry_status\tfinal_start0\tfinal_stop0");
			choices.writer.println("member\tthreshold_bp\tstrict_comparator\tmonitor_fold\tselected\tcontrol\tselection_rule");
			for(int t=0; t<unique; t++){
				choices.writer.print(row.clear().append(member).tab().append(thresholdName(thresholds[t]))
					.append("\tabs(median)>T\t4\t").append(selected[t] ? 1 : 0).tab()
					.append(thresholds[t]==0 ? "ALL_SUPPORTED" : thresholds[t]==Integer.MAX_VALUE ? "NEVER_SHIFT" : "BREAKPOINT")
					.append("\tMAX_END_BETTER_MINUS_WORSE_THEN_MIN_WORSE_THEN_MIN_MOVED_ALL_TIES\n"));
				for(int fold=0; fold<5; fold++){
					writePolicy(curve, "fold"+fold, training, fold, fold, thresholds[t], selected[t]);
				}
				if(selected[t] || thresholds[t]==0 || thresholds[t]==Integer.MAX_VALUE){
					for(final Dataset dataset:outerSets){
						writePolicy(outer, dataset.name, dataset.calls, -2, 5, thresholds[t], selected[t]);
						writeDecisions(decisions, dataset, thresholds[t], selected[t]);
					}
				}
			}
			population.writer.println("member\ttraining_rows\touter_training_rows_excluded\tinternal_retained\tmodels\tsupport_floor\tfolds\tmonitor\tqualification");
			population.writer.print(row.clear().append(member).tab().append(trainingRows).tab().append(trainingExcluded)
				.tab().append(training.size()).tab().append(models.size()).tab().append(SUPPORT)
				.append("\t5\t4\tCONDITIONAL_REAL_CALL_REPLAY_NOT_NESTED_MONITOR_OR_CALLER_REPRODUCTION\n"));
		}
	}

	/** A missing assembly is an error, never an implicit training assignment. */
	private void readFolds(){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(folds)){
			RrnaResourceIO.header(in.reader, NcrnaEndpointFoldInputs.ASSEMBLY_HEADER, folds);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==4 && p.length(0)>0 && p.termStartsWith("/", 1), "Complete fold mapping needs assembly, absolute FASTA, role and fold");
				final String assembly=p.parseString(0);
				final int fold=p.parseInt(3);
				require(!assemblyFolds.contains(assembly) && (p.termEquals("OUTER_TEST", 2) ? fold==-1
					: p.termEquals("INTERNAL_FOLD", 2) && fold>=0 && fold<5), "One exact assembly has exactly one declared fold role");
				assemblyFolds.put(assembly, fold+1);// Zero denotes explicit outer; missing returns -1.
			}
		}
		require(assemblyFolds.size()>5, "Complete corpus mapping needs outer and internal assemblies");
	}

	private void readDatasets(){
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final ObjectSet<String> trainingCalls=new ObjectSet<String>(String.class);
		final ObjectSet<String> trainingTruth=new ObjectSet<String>(String.class);
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, DATASET_HEADER, manifest);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==4 && p.length(0)>0 && names.add(p.parseString(0)) && p.termStartsWith("/", 2)
					&& paths.add(p.parseString(2)) && p.parseLong(3)>0, "Every geometry dataset needs a unique name/path and positive accepted count");
				final boolean train=p.termEquals("TRAIN", 1);
				require(train || p.termEquals("OUTER", 1), "Dataset role must be TRAIN or OUTER");
				final Dataset dataset=new Dataset(p.parseString(0));
				// Residual teaching never opens OUTER files, even to inspect their
				// truth or geometry. Their manifest labels remain an explicit audit.
				if(train || !residuals){
					readGeometry(p.parseString(2), p.parseLong(3), train, dataset,
						train ? trainingCalls : new ObjectSet<String>(String.class), train ? trainingTruth : new ObjectSet<String>(String.class));
				}
				if(!train){outerSets.add(dataset);}
			}
		}
		require(!training.isEmpty() && outerSets.size()==2 && trainingRows==training.size()+trainingExcluded,
			"The complete support pool and both existing outer panels must be explicitly accounted");
		final int[] n=new int[5];
		for(final Call call:training){n[call.fold]++;}
		for(final int count:n){require(count>0, "Every internal fold needs real matched calls");}
	}

	private void readGeometry(final String file, final long expected, final boolean train, final Dataset dataset,
			final ObjectSet<String> seenCalls, final ObjectSet<String> seenTruth){
		final LineParser1 p=new LineParser1('\t');
		long rows=0;
		boolean complete=false;
		try(Input in=new Input(file)){
			final boolean legacyMaximum=EvaluateNcrnaRealEndpoints.readGeometryHeader(in.reader, false, file);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(row.length>0 && row[0]=='#'){
					require(p.terms()==2 && p.termEquals("#COMPLETE", 0) && p.parseLong(1)==rows
						&& rows==expected && in.reader.nextLine()==null, "Geometry cache must close exactly at the bound producer count");
					complete=true;
					break;
				}
				require(p.terms()==16 && p.parseLong(0)==++rows && p.termEquals(member, 4), "Geometry rows must be dense and belong to the selected member");
				for(int i=0; i<15; i++){require(p.length(i)>0, "Required geometry fields cannot be empty");}
				final String assembly=p.parseString(2);
				require(assemblyFolds.contains(assembly), "Unknown canonical assembly in real-call cache: "+assembly);
				final int fold=assemblyFolds.get(assembly)-1;
				require(train || fold==-1, "An outer comparison must not contain an internal training assembly");
				final String model=NcrnaSupport.recordName(p.parseString(5));
				require(modelIndex.contains(model), "Actual winning model is absent from the frozen catalog: "+model);
				final Call call=new Call();
				call.input=p.parseString(1);
				call.assembly=assembly;
				call.reference=p.parseString(3);
				call.seqid=p.parseString(6);
				call.strand=p.parseString(7);
				call.model=modelIndex.get(model);
				call.fold=fold;
				call.length=p.parseInt(8);
				call.fixed=p.parseInt(9);
				call.start=p.parseInt(10);
				call.stop=p.parseInt(11);
				call.truthStart=p.parseInt(12);
				call.truthStop=p.parseInt(13);
				call.minimum=p.parseInt(14);
				call.maximum=NcrnaTextIO.maximumLength(p, 15, legacyMaximum);
				require((p.termEquals("+", 7) || p.termEquals("-", 7)) && call.fixed>=0 && call.fixed<=3
					&& call.length>0 && call.start>=0 && call.stop>=call.start && call.stop<call.length
					&& call.truthStart>=0 && call.truthStop>=call.truthStart && call.truthStop<call.length
					&& call.minimum>0 && call.maximum>=call.minimum
					&& ((call.fixed&1)==0 || call.start==0) && ((call.fixed&2)==0 || call.stop==call.length-1),
					"Geometry must retain real physical spans, fixed ends and caller length guards");
				final String key=assembly+'\t'+p.parseString(6)+'\t'+p.parseString(7)+'\t'+call.start+'\t'+call.stop;
				final String truth=assembly+'\t'+p.parseString(6)+'\t'+p.parseString(7)+'\t'+call.truthStart+'\t'+call.truthStop;
				require(seenCalls.add(key) && seenTruth.add(truth), "A source-qualified real call and reference may be consumed once per population");
				if(train){
					trainingRows++;
					if(fold<0){
						trainingExcluded++;
						continue;
					}
					training.add(call);
					final Model m=models.get(call.model);
					m.errors[fold][0].add(call.start-call.truthStart);
					m.errors[fold][1].add(call.stop-call.truthStop);
				}else{dataset.calls.add(call);}
			}
		}
		require(complete, "Truncated geometry cache lacks a terminal completion record: "+file);
	}

	/** Each policy uses one threshold for both ends and validates the joint pair.
	 * scope0/1 holds the other raw end; scope2 evaluates both decisions together. */
	private Stats[] score(final ArrayList<Call> calls, final int filterFold, final int learnedFold, final int threshold, final int scope){
		final Stats[] stats={new Stats(), new Stats()};
		final int[] ends=new int[2];
		for(final Call call:calls){
			if(filterFold>=0 && call.fold!=filterFold){continue;}
			final Model model=models.get(call.model);
			final boolean supported=model.support[learnedFold]>=SUPPORT;
			final int status=apply(call.start, call.stop, call.length, call.fixed, call.minimum, call.maximum,
				supported ? model.median[learnedFold][0] : 0, supported ? model.median[learnedFold][1] : 0, threshold, scope, ends);
			for(int end=0; end<2; end++){
				if(scope<2 && end!=scope){continue;}
				stats[end].add(call, end, ends, status, supported);
			}
		}
		return stats;
	}

	public static int apply(final int start, final int stop, final int length, final int fixed, final int minimum,
			final int maximum, final int median5, final int median3, final int threshold, final int scope, final int[] ends){
		require(threshold>=0 && scope>=0 && scope<=2 && ends!=null && ends.length==2 && fixed>=0 && fixed<=3,
			"Threshold policy requires a nonnegative magnitude, explicit scope and two output ends");
		require(median5!=Integer.MIN_VALUE && median3!=Integer.MIN_VALUE,
			"Differences between physical nonnegative int coordinates cannot equal Integer.MIN_VALUE");
		final long five=(fixed&1)!=0 || scope==1 || Math.abs((long)median5)<=threshold ? start : (long)start-median5;
		final long three=(fixed&2)!=0 || scope==0 || Math.abs((long)median3)<=threshold ? stop : (long)stop-median3;
		final int status=EvaluateNcrnaRealEndpoints.geometry(five, three, length, minimum, maximum);
		ends[0]=status==0 ? (int)five : start;
		ends[1]=status==0 ? (int)three : stop;
		return status;
	}

	/** Selection uses joint-policy per-end counts on fold4, never outer metrics. */
	private static int compare(final Stats[] a, final Stats[] b){
		final long aw=a[0].worse+a[1].worse;
		final long bw=b[0].worse+b[1].worse;
		final long an=a[0].better+a[1].better-aw;
		final long bn=b[0].better+b[1].better-bw;
		if(an!=bn){return Long.compare(an, bn);}
		if(aw!=bw){return Long.compare(bw, aw);}
		return Long.compare(b[0].moved+b[1].moved, a[0].moved+a[1].moved);
	}

	private void writePolicy(final Output out, final String panel, final ArrayList<Call> calls,
			final int filterFold, final int learnedFold, final int threshold, final boolean selected){
		assert(learnedFold>=0 && learnedFold<=5) : "Use exactly one held-fold or full-training offset resource";
		final ByteBuilder row=new ByteBuilder();
		for(int scope=0; scope<3; scope++){
			final Stats[] stats=score(calls, filterFold, learnedFold, threshold, scope);
			for(int end=0; end<2; end++){
				if(scope<2 && end!=scope){continue;}
				stats[end].write(row.clear().append(member).tab().append(panel).tab().append(thresholdName(threshold))
					.tab().append(selected ? 1 : 0).tab().append(scope==2 ? "PAIRED_JOINT" : "SINGLE_END_OTHER_RAW")
					.tab().append(end==0 ? "5prime" : "3prime"), scope);
				out.writer.print(row);
			}
		}
	}

	private void writeOffsets(){
		for(int fold=0; fold<6; fold++){
			try(Output out=new Output(path(fold==5 ? "full.offsets.tsv" : "fold"+fold+".offsets.tsv"))){
				out.writer.println(NcrnaFixedEndpointTrials.OFFSET_HEADER);
				final ByteBuilder row=new ByteBuilder();
				for(final Model model:models){
					out.writer.print(row.clear().append(member).tab().append(model.name).tab().append(model.support[fold]).tab().append(SUPPORT)
						.tab().append(model.support[fold]>=SUPPORT ? model.median[fold][0] : 0)
						.tab().append(model.support[fold]>=SUPPORT ? model.median[fold][1] : 0).nl());
				}
			}
		}
	}

	/** Keep chosen/control per-call dispositions beside aggregate metrics so a
	 * reviewer can trace fixed, unsupported and joint-hold cases without replay. */
	private void writeDecisions(final Output out, final Dataset dataset, final int threshold, final boolean selected){
		assert(dataset!=null && !dataset.calls.isEmpty()) : "Outer audit retains every bound matched call";
		final ByteBuilder row=new ByteBuilder();
		final int[] ends=new int[2];
		for(final Call call:dataset.calls){
			final Model model=models.get(call.model);
			final boolean supported=model.support[5]>=SUPPORT;
			final int five=supported ? model.median[5][0] : 0;
			final int three=supported ? model.median[5][1] : 0;
			for(int scope=0; scope<3; scope++){
				final int status=apply(call.start, call.stop, call.length, call.fixed, call.minimum, call.maximum,
					five, three, threshold, scope, ends);
				out.writer.print(row.clear().append(member).tab().append(dataset.name).tab().append(thresholdName(threshold))
					.tab().append(selected ? 1 : 0).tab().append(scope==2 ? "PAIRED_JOINT" : scope==0 ? "SINGLE_5PRIME" : "SINGLE_3PRIME")
					.tab().append(call.input).tab().append(call.assembly).tab().append(call.reference).tab().append(model.name)
					.tab().append(call.seqid).tab().append(call.strand).tab().append(call.length).tab().append(call.fixed)
					.tab().append(model.support[5]).tab().append(five).tab().append(three)
					.tab().append(call.start).tab().append(call.stop).tab().append(call.truthStart).tab().append(call.truthStop)
					.tab().append(GEOMETRY[status]).tab().append(ends[0]).tab().append(ends[1]).nl());
			}
		}
	}

	/** Reuse the accepted fold/full offsets rather than refitting medians or
	 * repeating threshold selection. Support counts bind them to these same
	 * training witnesses and leave-fold-out policy. Unsupported shifts stay zero. */
	private void readOffsets(){
		assert(residuals && offsetDirectory!=null) : "Only residual export imports accepted offset resources";
		final LineParser1 p=new LineParser1('\t');
		for(int held=0; held<6; held++){
			final String file=new File(offsetDirectory, held==5 ? "full.offsets.tsv" : "fold"+held+".offsets.tsv").getPath();
			int count=0;
			try(Input in=new Input(file)){
				RrnaResourceIO.header(in.reader, NcrnaFixedEndpointTrials.OFFSET_HEADER, file);
				for(byte[] row; (row=in.reader.nextLine())!=null;){
					p.set(row);
					require(p.terms()==6 && count<models.size() && p.termEquals(member, 0) && p.parseInt(3)==SUPPORT,
						"Imported offsets need the exact member, model axis and support floor");
					final Model model=models.get(count++);
					require(p.termEquals(model.name, 1), "Imported offsets must preserve every frozen model in order");
					int expected=0;
					for(int fold=0; fold<5; fold++){if(fold!=held){expected=Math.addExact(expected, model.errors[fold][0].size);}}
					require(p.parseInt(2)==expected, "Imported median support differs from this exact training-fold population");
					model.support[held]=expected;
					model.median[held][0]=p.parseInt(4);
					model.median[held][1]=p.parseInt(5);
					require(expected>=SUPPORT || model.median[held][0]==0 && model.median[held][1]==0,
						"Unsupported models must retain zero applied offsets");
				}
			}
			require(count==models.size(), "Imported offsets must include unsupported models, not only observed winners");
		}
	}

	/** Residual weights come from each training witness's actual winning model
	 * and the same full physical/joint guards used by replay. Raw/held errors
	 * remain represented; each end contributes exactly one count per retained call. */
	private void writeResiduals(){
		assert(residuals) : "Residual export must not rerun threshold selection or outer evaluation";
		final ByteBuilder row=new ByteBuilder();
		final int[] ends=new int[2];
		try(Output accounting=new Output(path("training.tsv"))){
			accounting.writer.println("member\theld_fold\ttraining_retained\tsupported\tfixed5\tfixed3\tgeometry_held\tcrossed\toutside_contig\tbelow_minimum\tabove_maximum\traw_below_minimum\tqualification");
			for(int held=0; held<6; held++){
				final IntLongHashMap[] hist={new IntLongHashMap(), new IntLongHashMap()};
				final long[] bad=new long[5];
				long n=0, supported=0, fixed5=0, fixed3=0, rawShort=0;
				for(final Call call:training){
					if(call.fold==held){continue;}
					final Model model=models.get(call.model);
					final boolean support=model.support[held]>=SUPPORT;
					final int status=apply(call.start, call.stop, call.length, call.fixed, call.minimum, call.maximum,
						support ? model.median[held][0] : 0, support ? model.median[held][1] : 0, 0, 2, ends);
					hist[0].add(ends[0]-call.truthStart);
					hist[1].add(ends[1]-call.truthStop);
					n++;
					if(support){supported++;}
					if((call.fixed&1)!=0){fixed5++;}
					if((call.fixed&2)!=0){fixed3++;}
					if((long)call.stop-call.start+1<call.minimum){rawShort++;}
					bad[status]++;
				}
				require(n>0 && hist[0].sum()==n && hist[1].sum()==n,
					"Every retained training call contributes one residual per end, including zero, unsupported and held cases");
				try(Output out=new Output(path(held==5 ? "full.histogram.tsv" : "fold"+held+".histogram.tsv"))){
					out.writer.println("family\tend\tmeasure\toffset_bp\tloci");
					for(int end=0; end<2; end++){
						final int[] offsets=hist[end].toArray();
						Arrays.sort(offsets);
						for(final int offset:offsets){
							out.writer.print(row.clear().append(member).tab().append(end==0 ? "5prime" : "3prime")
								.append("\tsigned\t").append(offset).tab().append(hist[end].get(offset)).nl());
						}
					}
				}
				accounting.writer.print(row.clear().append(member).tab().append(held==5 ? "full" : Integer.toString(held))
					.tab().append(n).tab().append(supported).tab().append(fixed5).tab().append(fixed3)
					.tab().append(n-bad[0]).tab().append(bad[1]).tab().append(bad[2]).tab().append(bad[3]).tab().append(bad[4])
					.tab().append(rawShort).append("\tREAL_WINNER_FOLD_TRAINING_RESIDUALS_PHYSICAL_FIXED_AND_LENGTH_GUARDS\n"));
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Model{
		Model(final String name_){
			name=name_;
			for(int f=0; f<5; f++){
				for(int e=0; e<2; e++){errors[f][e]=new IntList();}
			}
		}
		void learn(){
			final IntList merged=new IntList();
			for(int held=0; held<6; held++){
				for(int end=0; end<2; end++){
					merged.clear();
					for(int fold=0; fold<5; fold++){
						if(fold==held){continue;}
						final IntList values=errors[fold][end];
						for(int i=0; i<values.size; i++){merged.add(values.array[i]);}
					}
					merged.sort();
					if(end==0){support[held]=merged.size;}
					else{require(support[held]==merged.size, "Both endpoint medians use the same real witness cohort");}
					median[held][end]=merged.size==0 ? 0 : merged.array[(merged.size-1)/2];
				}
			}
		}
		final String name;
		final IntList[][] errors=new IntList[5][2];
		final int[][] median=new int[6][2];
		final int[] support=new int[6];
	}
	private static final class Call{
		String input, assembly, reference, seqid, strand;
		int fold, model, length, fixed, start, stop, truthStart, truthStop, minimum, maximum;
	}
	private static final class Dataset{
		Dataset(final String name_){name=name_;}
		final String name;
		final ArrayList<Call> calls=new ArrayList<Call>();
	}
	private static final class Stats{
		void add(final Call call, final int end, final int[] ends, final int status, final boolean supported){
			final long raw=(end==0 ? call.start-call.truthStart : call.stop-call.truthStop);
			final long error=(long)ends[end]-(end==0 ? call.truthStart : call.truthStop);
			final long other=(long)ends[1-end]-(end==0 ? call.truthStop : call.truthStart);
			final long a=Math.abs(error), b=Math.abs(raw);
			n++;
			if(a==0){exact++;}
			if(a<=3){within3++;}
			if(a<=10){within10++;}
			if(a==0 && other==0){bothExact++;}
			if(a<=3 && Math.abs(other)<=3){both3++;}
			if(a<=10 && Math.abs(other)<=10){both10++;}
			if(a<b){better++;}else if(a>b){worse++;}else{tied++;}
			if(error!=raw){moved++;}
			if(b==0 && a>0){exactHarmed++;}
			if(b<=10){
				nearN++;
				nearRawSq+=(double)raw*raw;
				nearSq+=(double)error*error;
				if(a>b){nearHarmed++;}
			}
			if(b>10 && a<b){distantImproved++;}
			if(b>10 && a>b){distantHarmed++;}
			if((call.fixed&(1<<end))!=0){fixed++;}
			if(!supported){unsupported++;}
			bad[status]++;
			rawSq+=(double)raw*raw;
			sq+=(double)error*error;
			assert(better+worse+tied==n && exact<=within3 && within3<=within10 && within10<=n)
				: "Disjoint end outcomes and nested tolerances must conserve every real matched call";
		}
		void write(final ByteBuilder row, final int scope){
			require(n>0, "Do not report empty matched-call cohorts as measurements");
			row.tab().append(n).tab().append(exact).tab().append(within3).tab().append(within10);
			if(scope==2){row.tab().append(bothExact).tab().append(both3).tab().append(both10);}
			else{row.append("\tNA\tNA\tNA");}
			row.tab().append(Math.sqrt(rawSq/n), 8).tab().append(Math.sqrt(sq/n), 8)
				.tab().append(Math.sqrt(sq/n)-Math.sqrt(rawSq/n), 8)
				.tab().append(better).tab().append(worse).tab().append(tied).tab().append(moved)
				.tab().append(exactHarmed).tab().append(nearN).tab();
			if(nearN==0){row.append("NA\tNA");}
			else{row.append(Math.sqrt(nearRawSq/nearN), 8).tab().append(Math.sqrt(nearSq/nearN), 8);}
			row.tab().append(nearHarmed).tab().append(distantImproved).tab().append(distantHarmed)
				.tab().append(fixed).tab().append(unsupported).tab().append(n-bad[0])
				.tab().append(bad[1]).tab().append(bad[2]).tab().append(bad[3]).tab().append(bad[4]).tab();
			if(better+worse==0){row.append("NA");}
			else{row.append((double)better/(better+worse), 8);}
			row.append("\tCONDITIONAL_MATCHED_REAL_CALLS_NOT_CALLER_REPRODUCTION\n");
		}
		long n, exact, within3, within10, bothExact, both3, both10, better, worse, tied, moved;
		long exactHarmed, nearN, nearHarmed, distantImproved, distantHarmed, fixed, unsupported;
		double rawSq, sq, nearRawSq, nearSq;
		final long[] bad=new long[5];
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	private String path(final String name){return new File(output, name).getAbsolutePath();}
	private static String thresholdName(final int value){return value==Integer.MAX_VALUE ? "never" : Integer.toString(value);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	/*--------------------------------------------------------------*/
	/*----------------       Fields / Constants     ----------------*/
	/*--------------------------------------------------------------*/

	public static final String DATASET_HEADER="dataset\trole\tgeometry\texpected_rows";
	private static final String RESULT_HEADER="member\tpanel\tthreshold_bp\tselected_on_fold4\tscope\tend\tmatched\texact\twithin3\twithin10\tboth_exact\tboth_within3\tboth_within10\traw_rms_bp\trms_bp\trms_minus_raw_bp\tbetter\tworse\ttied\tmoved\texact_harmed\traw_near10_n\traw_near10_rms_bp\tnear10_rms_bp\tnear10_harmed\tdistant_improved\tdistant_harmed\tfixed\tunsupported\tinvalid_joint\tcrossed\toutside_contig\tbelow_minimum\tabove_maximum\tbetter_fraction\tqualification";
	private static final String[] GEOMETRY={"VALID", "CROSSED", "OUTSIDE_CONTIG", "BELOW_MINIMUM", "ABOVE_MAXIMUM"};
	private static final int SUPPORT=10;
	private String manifest, folds, reference, member, output, offsetDirectory;
	private boolean residuals;
	private long trainingRows, trainingExcluded;
	private final ObjectIntMap<String> assemblyFolds=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> modelIndex=new ObjectIntMap<String>(String.class);
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ArrayList<Call> training=new ArrayList<Call>();
	private final ArrayList<Dataset> outerSets=new ArrayList<Dataset>();
}
