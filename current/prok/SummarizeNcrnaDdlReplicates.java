package prok;

import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Pools disjoint folds before crossed hash-seed/assignment uncertainty.
 * Reads completed small result tables only; never reloads biological data.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaDdlReplicates{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaDdlReplicates.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			long[] seeds=null, assignments=null;
			String followup=null, sizegrid=null, sizeextension=null, definitions=null, gcControl=null, gcAllControl=null, gcFollowup=null, gcPower=null;
			int shortSize=-1, longSize=-1;
			boolean gc=false, gcAll=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("seeds")){seeds=numbers(b);}
				else if(a.equals("foldseeds")){assignments=numbers(b);}
				else if(a.equals("followup")){followup=b;}
				else if(a.equals("sizegrid")){sizegrid=b;}
				else if(a.equals("sizeextension")){sizeextension=b;}
				else if(a.equals("definitions")){definitions=b;}
				else if(a.equals("shortsize")){shortSize=Parse.parseIntKMG(b);}
				else if(a.equals("longsize")){longSize=Parse.parseIntKMG(b);}
				else if(a.equals("gc")){gc=Parse.parseBoolean(b);}
				else if(a.equals("gcall")){gcAll=Parse.parseBoolean(b);}
				else if(a.equals("gcallfollowup")){gcFollowup=b;}
				else if(a.equals("gcpower")){gcPower=b;}
				else if(a.equals("gccontrol")){gcControl=b;}
				else if(a.equals("gcallcontrol")){gcAllControl=b;}
				else{require(parser.parse(arg, a, b), "Unknown replicate-report flag: "+arg);}
			}
			require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use in=complete_worker_roster out=fresh_prefix;no overwrite,append,or sampling");
			if(gcPower!=null){
				writeGcCosts(parser.in1, gcPower, parser.out1);
				pp.outstream.println("NCRNA_DDL_GC_SAVED_RANK_COST_INFERENCE_COMPLETE_NOT_COUNTED_ALIGNMENTS");
				return;
			}
			if(gcAllControl!=null){
				checkGcAllControl(parser.in1, gcAllControl, parser.out1);
				pp.outstream.println("NCRNA_DDL_ALL_DATA_GC_W0_EXACT_CONTROL_PARITY_PASS");
				return;
			}
			if(gcControl!=null){
				checkGcControl(parser.in1, gcControl, parser.out1);
				pp.outstream.println("NCRNA_DDL_GC_W0_EXACT_SAVED_GRID_PARITY_PASS");
				return;
			}
			new SummarizeNcrnaDdlReplicates(seeds, assignments, followup, sizegrid, sizeextension, definitions, shortSize, longSize, gc || gcAll, gcAll, gcFollowup).process(parser.in1, parser.out1);
			pp.outstream.println("NCRNA_DDL_CROSSED_REPLICATE_REPORT_COMPLETE_NOT_CALLER");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private SummarizeNcrnaDdlReplicates(final long[] seeds_, final long[] assignments_,
			final String followup_, final String sizegrid_, final String sizeextension_, final String definitions_, final int shortSize, final int longSize, final boolean gc_, final boolean gcAll_, final String gcFollowup_){
		require(seeds_!=null && assignments_!=null && seeds_.length>=2 && assignments_.length>=2,
			"Crossed reporting requires explicit distinct hash seeds and at least two fold assignments");
		seeds=seeds_;
		assignments=assignments_;
		followup=followup_;
		sizegrid=sizegrid_;
		sizeextension=sizeextension_;
		definitions=definitions_;
		gc=gc_;
		gcAll=gcAll_;
		gcFollowup=gcFollowup_;
		require(gcFollowup==null || gcAll, "gcallfollowup imports only new settings alongside accepted all-data GC workers");
		gcVariants=gcAll ? BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS.length : BuildNcrnaDdlLayoutCounts.GC_SETTINGS.length;
		require(!gc || followup==null && sizegrid==null && sizeextension==null && definitions==null,
			"GC reporting uses only its newly measured matched controls, not an older blacklist/bonus observation");
		require((sizegrid==null || followup!=null) && (sizeextension==null || sizegrid!=null) && (definitions==null || sizegrid!=null)
			&& (definitions==null ? shortSize==-1 && longSize==-1 : shortSize>=0 && longSize>=0),
			"Size grids extend the completed follow-up;definition controls require that grid and its explicitly selected sizes");
		definitionStart=sizeextension!=null ? 160 : 138;
		methodCount=gc ? 4*gcVariants+(gcFollowup==null ? 0 : 2*BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length)
			: (sizegrid!=null ? definitionStart : followup!=null ? 68 : 50)+(definitions!=null ? 3 : 0);
		if(gc){
			for(final int bins:new int[]{512,1024}){
				for(int setting=0; setting<2*gcVariants; setting++){
					register(gcAll ? BuildNcrnaDdlLayoutCounts.gcAllMethod(bins, setting) : BuildNcrnaDdlLayoutCounts.gcMethod(bins, setting),
						new int[]{2,bins,2,1,bins==512 ? 60 : 100,256,2048});
				}
			}
			if(gcFollowup!=null){
				for(final int bins:new int[]{512,1024}){
					for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length; setting++){
						register(BuildNcrnaDdlLayoutCounts.gcExtraMethod(bins, setting), new int[]{2,bins,2,1,bins==512 ? 60 : 100,256,2048});
					}
				}
			}
		}else{
			for(int i=0; i<3; i++){
				for(final int bins:new int[]{512,1024}){
					for(int layout=0; layout<3; layout++){
						for(int filtered=0; filtered<2; filtered++){
							final int next=methods.size();
							final String name="s"+(i==2 ? 8 : 7)+"_m"+(i==0 ? 9 : 10)+"_"+LAYOUTS[layout]+"_"+bins+"_"+(filtered==0 ? "plain" : "padding");
							methods.put(name, next);
							names[next]=name;
							spec[next]=new int[]{i,bins,layout,filtered,-1};
						}
					}
				}
			}
			for(final int bonus:new int[]{0,5,1000000}){method("bonus_fixed_"+bonus, bonus);}
			for(int tenth=1; tenth<=10; tenth++){method("bonus_tenths_"+tenth, -tenth);}
			method("CONSENSUS_INDEX", Integer.MAX_VALUE);
			require(methods.size()==50, "The declared36methods plus13additional bonuses and one reference must remain complete");
			if(followup!=null){
				for(final int bins:new int[]{512,1024}){
					for(final int bonus:new int[]{60,80,120}){
						register("s8_m10_quarter_"+bins+"_padding_bonus"+bonus, new int[]{2,bins,2,1,bonus});
					}
				}
				for(final int bins:new int[]{512,1024}){
					for(int layout=0; layout<3; layout++){
						for(int filter=0; filter<2; filter++){
							register("s8_m11_"+LAYOUTS[layout]+"_"+bins+"_"+(filter==0 ? "plain" : "padding"),
								new int[]{3,bins,layout,filter,100});
						}
					}
				}
			}
			if(sizegrid!=null){
				for(final int shortPrefix:BuildNcrnaDdlLayoutCounts.GRID_SHORT){
					for(final int longPrefix:BuildNcrnaDdlLayoutCounts.GRID_LONG){
						if(shortPrefix==0 && longPrefix==0 || shortPrefix==256 && longPrefix==4096){continue;}
						register(BuildNcrnaDdlLayoutCounts.sizeMethod(shortPrefix, longPrefix),
							new int[]{2,1024,2,1,100,shortPrefix,longPrefix});
					}
				}
			}
			if(sizeextension!=null){
				for(final int[] size:BuildNcrnaDdlLayoutCounts.SIZE_EXTENSION){
					register(BuildNcrnaDdlLayoutCounts.sizeMethod(size[0], size[1]), new int[]{2,1024,2,1,100,size[0],size[1]});
				}
			}
			if(definitions!=null){
				for(final int definition:new int[]{1,3,0}){
					register("blacklist_"+BuildNcrnaDdlLayoutCounts.DEFINITION_NAMES[definition]+"_s"+shortSize+"_l"+longSize,
						new int[]{2,1024,2,1,100,shortSize,longSize});
				}
			}
			require(methods.size()==methodCount, "Each released grid adds its exact new settings and imports overlapping controls once");
		}
		values=new double[3*methodCount][assignments.length+1][seeds.length][gc ? METRICS.length : 13];
		masks=new int[3*methodCount][assignments.length+1][seeds.length];
		pooledHistograms=gc ? new long[3*methodCount][2][509] : followup==null ? null : new long[3*methodCount][2][508];
		for(final double[][][] method:values){for(final double[][] assignment:method){for(final double[] row:assignment){Arrays.fill(row, Double.NaN);}}}
	}

	private void method(final String name, final int bonus){
		final int index=methods.size();
		methods.put(name, index);
		names[index]=name;
		spec[index]=new int[]{0,512,0,0,bonus};
	}
	private void register(final String name, final int[] shape){
		assert(methods.get(name)<0 && (shape.length==5 || shape.length==7)) : "Each declared follow-up setting owns one new method ordinal and complete shape";
		final int next=methods.size();
		methods.put(name, next);
		names[next]=name;
		spec[next]=shape;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Completed Tables       ----------------*/
	/*--------------------------------------------------------------*/

	/** A deterministic lookup/arithmetic join,not another ranking or caller run.
	 * Surviving loci use their saved best-window rank among retained models.
	 * A failed locus exhausts its least-populated retained candidate list.
	 * This is an optimistic one-window proxy with no cap/patience simulation. */
	public static void writeGcCosts(final String summary, final String power, final String prefix){
		require(Tools.testInputFiles(false, true, summary, power)
			&& Tools.testOutputFiles(false, false, false, prefix+".gc_cost.tsv")
			&& Tools.testForDuplicateFiles(true, summary, power, prefix+".gc_cost.tsv"),
			"Rank inference joins completed summaries to retained-candidate counts and cannot overwrite either source");
		final ObjectMap<String,byte[]> gates=new ObjectMap<String,byte[]>(String.class, byte[].class);
		final java.util.ArrayList<String> order=new java.util.ArrayList<String>();
		final LineParser1 p=new LineParser1('\t');
		String cohort=null;
		int models=0, parents=0;
		try(Input in=new Input(power)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.GC_POWER_HEADER, power);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==20, "Cost inference requires retained-model counts and actual failed-locus candidate sums");
				if(!p.termEquals("REBUILT_FIVE_FOLD", 0) && !p.termEquals("FIXED_NOT_SEEN", 0)){continue;}
				if(cohort==null){
					cohort=p.parseString(0);
					models=p.parseInt(3);
					parents=p.parseInt(11);
				}
				require(p.termEquals(cohort, 0) && p.parseInt(3)==models && p.parseInt(11)==parents,
					"Every cost row retains one complete original population and its full candidate-model roster");
				final String kind=p.parseString(1);
				final String key=gcGateKey(kind, p.parseDouble(2));
				require(!gates.contains(key), "Each native gate contributes once to the cost lookup: "+key);
				gates.put(key, row);
				order.add(key);
			}
		}
		require(cohort!=null && gates.contains(gcGateKey("NONE", 0)), "An inferred saving needs an explicit no-filter power baseline");
		final ObjectIntMap<String> methods=new ObjectIntMap<String>(String.class);
		for(final int bins:new int[]{512,1024}){
			for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS.length; setting++){
				methods.put(BuildNcrnaDdlLayoutCounts.gcAllMethod(bins, setting), methods.size());
			}
		}
		final double[][] ranks=new double[methods.size()][2], losses=new double[methods.size()][2];
		final int[] seen=new int[methods.size()];
		try(Input in=new Input(summary)){
			RrnaResourceIO.header(in.reader, "cohort\tmethod\tmetric\tassignments\thash_seeds\tmean\tstandard_deviation\thash_variance\tsplit_variance\tinteraction_variance\tmean_standard_error\tparents\tlow_parents\tassembly_contribution", summary);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals(cohort, 0) || !p.termEquals("ALL", 13)){continue;}
				final int method=methods.get(p.parseString(1));
				if(method<0){continue;}
				final int bit=p.termEquals("rank_sum", 2) ? 1 : p.termEquals("parents_all_acceptable_filtered", 2) ? 2 : 0;
				if(bit==0){continue;}
				require((seen[method]&bit)==0 && p.parseInt(11)==parents,
					"Each imported rank/loss estimate belongs once to the same original locus population");
				seen[method]|=bit;
				final double[] values=bit==1 ? ranks[method] : losses[method];
				values[0]=Double.parseDouble(p.parseString(5));
				values[1]=Double.parseDouble(p.parseString(6));
			}
		}
		for(final int mask:seen){require(mask==3, "The accepted all-data summary retains every original window-statistic rank and loss estimate");}
		final ByteBuilder b=new ByteBuilder(4096);
		try(Output out=new Output(prefix+".gc_cost.tsv")){
			out.writer.println("cohort\tbins\tfilter\tparameter\twindows\tbaseline_models_per_window\tmean_models_retained\tmedian_models_retained\tmean_models_removed\tpercent_candidates_removed\tloci\tloci_lost\trank_available\tinferred_alignments_mean\tinferred_alignments_sd\tunfiltered_inferred_alignments_mean\tinferred_alignment_count_saved\tinferred_alignment_percent_saved\tfailed_locus_surviving_candidate_sum\tinterpretation");
			for(final int bins:new int[]{512,1024}){
				final int base=methods.get(BuildNcrnaDdlLayoutCounts.gcAllMethod(bins, 0));
				require(seen[base]==3 && losses[base][0]==0 && ranks[base][0]>0,
					"The no-GC rank baseline must be measured,positive and free of filtered failures");
				for(final String key:order){
					p.set(gates.get(key));
					final String kind=p.parseString(1);
					final double parameter=p.parseDouble(2);
					int setting=kind.equals("NONE") ? 0 : kind.equals("LOSSLESS_SD") ? 11 : -1;
					for(int s=6; s<BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS.length; s++){
						final String name=BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS[s];
						if(kind.equals("SD") && name.startsWith("filter_z") && Double.parseDouble(name.substring(8))==parameter
							|| kind.equals("MARGIN") && name.startsWith("filter_m") && Double.parseDouble(name.substring(8))==parameter){setting=s;}
					}
					final int method=setting<0 ? -1 : methods.get(BuildNcrnaDdlLayoutCounts.gcAllMethod(bins, setting));
					final boolean available=method>=0 && seen[method]==3;
					final long lost=p.parseLong(12), failedCandidates=p.parseLong(18);
					b.clear().append(cohort).tab().append(bins).tab().append(kind).tab().append(Double.toString(parameter))
						.tab().append(p.parseLong(4)).tab().append(models).tab().append(p.parseString(15)).tab().append(p.parseString(16))
						.tab().append(p.parseString(5)).tab().append(p.parseString(17)).tab().append(parents).tab().append(lost).tab().append(available);
					if(available){
						require(losses[method][0]==lost && losses[method][1]==0,
							"Fixed all-data gates lose the same loci in every hash/assignment observation;otherwise a constant failure-cost correction is invalid");
						final double inferred=ranks[method][0]-(models+1.0)*lost+failedCandidates;
						require(inferred>=0 && inferred<=((long)parents)*models,
							"One best-window inferred walk costs a surviving rank or a retained candidate list,never the failure sentinel");
						final double saved=ranks[base][0]-inferred;
						b.tab().append(inferred, 9).tab().append(ranks[method][1], 9).tab().append(ranks[base][0], 9)
							.tab().append(saved, 9).tab().append(100*saved/ranks[base][0], 9);
					}else{b.append("\tNA\tNA\t").append(ranks[base][0], 9).append("\tNA\tNA");}
					out.writer.print(b.tab().append(failedCandidates).append("\tOPTIMISTIC_ONE_BEST_WINDOW_PER_LOCUS_EXHAUST_RETAINED_CANDIDATES_ON_FAILURE_NO_CAP_OR_PATIENCE_NOT_COUNTED_ALIGNMENTS\n"));
				}
			}
		}
	}
	private static String gcGateKey(final String kind, final double parameter){
		assert(kind!=null && Double.isFinite(parameter)) : "A gate lookup binds its original filter kind and exact finite parameter";
		return kind.equals("LOSSLESS_SD") ? kind : kind+"_"+Double.toString(parameter);
	}

	/** Exact pilot control at the already measured1024-bin256/2048setting.
	 * Compares populations and all six count/rank metrics, per fold and stratum. */
	public static void checkGcControl(final String actual, final String previous, final String prefix){
		require(Tools.testInputFiles(false, true, actual, previous)
			&& Tools.testOutputFiles(false, false, false, prefix+".parity.tsv")
			&& Tools.testForDuplicateFiles(true, actual, previous, prefix+".parity.tsv"),
			"Control parity reads distinct completed tables and writes a fresh result");
		final long[][][] expected=new long[5][3][];
		final LineParser1 p=new LineParser1('\t');
		long assignment=Long.MIN_VALUE, seed=-1;
		int count=0;
		try(Input in=new Input(previous)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER, previous);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals(BuildNcrnaDdlLayoutCounts.sizeMethod(256, 2048), 3)){continue;}
				final int fold=p.parseInt(2), stratum=contribution(p);
				if(count==0){assignment=p.parseLong(0); seed=p.parseLong(1);}
				require(p.terms()==30 && p.parseLong(0)==assignment && p.parseLong(1)==seed && fold>=0 && fold<(assignment<0 ? 2 : 5)
					&& stratum<(assignment<0 ? 1 : 3) && expected[fold][stratum]==null,
					"The saved size control must contain each complete matched assignment/seed/fold/population once");
				expected[fold][stratum]=new long[8];
				for(int metric=0; metric<8; metric++){expected[fold][stratum][metric]=p.parseLong(12+metric);}
				count++;
			}
		}
		require(count==(assignment<0 ? 2 : 15), "The saved grid control must retain all judging folds and contribution populations");
		final int[][] seen=new int[5][3];
		int checked=0;
		try(Input in=new Input(actual)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER, actual);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int arm=p.termEquals(BuildNcrnaDdlLayoutCounts.gcMethod(1024, 0), 3) ? 0
					: p.termEquals(BuildNcrnaDdlLayoutCounts.gcMethod(1024, 12), 3) ? 1 : -1;
				if(arm<0){continue;}
				final int fold=p.parseInt(2), stratum=contribution(p);
				require(p.terms()==39 && p.parseLong(0)==assignment && p.parseLong(1)==seed
					&& fold>=0 && fold<(assignment<0 ? 2 : 5) && stratum<(assignment<0 ? 1 : 3)
					&& expected[fold][stratum]!=null && (seen[fold][stratum]&(1<<arm))==0,
					"New w0 controls must use the exact saved seed,assignment,fold and population");
				seen[fold][stratum]|=1<<arm;
				for(int metric=0; metric<8; metric++){
					require(p.parseLong(12+metric)==expected[fold][stratum][metric],
						"GC w0 changed a saved population/count/rank metric: fold="+fold+", population="+stratum+", field="+(12+metric));
				}
				checked++;
			}
		}
		require(checked==2*count, "Both GC-source controls must exactly reproduce every saved control row");
		try(Output out=new Output(prefix+".parity.tsv")){
			out.writer.println("assignment_seed\tddl_seed\tmatched_control_rows\tverdict");
			out.writer.println(assignment+"\t"+seed+"\t"+checked+"\tEXACT_POPULATIONS_AND_SIX_METRICS_PASS");
		}
	}

	/** Changing only GC statistics cannot change any no-GC rank or population.
	 * Both widths,sources,all folds/populations and full rank histograms agree. */
	public static void checkGcAllControl(final String actual, final String previous, final String prefix){
		require(Tools.testInputFiles(false, true, actual, previous)
			&& Tools.testOutputFiles(false, false, false, prefix+".parity.tsv")
			&& Tools.testForDuplicateFiles(true, actual, previous, prefix+".parity.tsv"),
			"All-data control parity requires distinct saved tables and a fresh receipt");
		final byte[][][][] expected=new byte[4][5][3][];
		final boolean[][][] seen=new boolean[4][5][3];
		final LineParser1 p=new LineParser1('\t'), q=new LineParser1('\t');
		long assignment=Long.MIN_VALUE, seed=-1;
		int rows=0, checked=0;
		try(Input in=new Input(previous)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER, previous);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int control=gcControlOrdinal(p, false);
				if(control<0){continue;}
				if(rows==0){assignment=p.parseLong(0); seed=p.parseLong(1);}
				final int fold=p.parseInt(2), stratum=contribution(p);
				require(p.terms()==39 && p.parseLong(0)==assignment && p.parseLong(1)==seed
					&& fold>=0 && fold<(assignment<0 ? 2 : 5) && stratum<(assignment<0 ? 1 : 3)
					&& expected[control][fold][stratum]==null, "Each saved no-GC control must retain its exact worker,fold and population once");
				expected[control][fold][stratum]=row;
				rows++;
			}
		}
		require(rows==(assignment<0 ? 8 : 60), "Both bin widths and sources must retain every saved no-GC fold/population");
		try(Input in=new Input(actual)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER, actual);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int control=gcControlOrdinal(p, true);
				if(control<0){continue;}
				final int fold=p.parseInt(2), stratum=contribution(p);
				require(p.terms()==39 && p.parseLong(0)==assignment && p.parseLong(1)==seed
					&& fold>=0 && fold<(assignment<0 ? 2 : 5) && stratum<(assignment<0 ? 1 : 3)
					&& expected[control][fold][stratum]!=null && !seen[control][fold][stratum],
					"All-data statistics cannot omit or duplicate a saved no-GC comparison cell");
				q.set(expected[control][fold][stratum]);
				for(int field=12; field<20; field++){
					require(p.parseLong(field)==q.parseLong(field), "All-data GC changed a no-GC population/count/rank metric: field="+field);
				}
				for(int field=28; field<30; field++){
					require(p.termEquals(q.parseString(field), field), "All-data GC changed a no-GC full rank histogram");
				}
				seen[control][fold][stratum]=true;
				checked++;
			}
		}
		require(checked==rows, "All-data GC must preserve every saved no-GC control cell");
		try(Output out=new Output(prefix+".parity.tsv")){
			out.writer.println("assignment_seed\tddl_seed\tmatched_control_rows\tverdict");
			out.writer.println(assignment+"\t"+seed+"\t"+checked+"\tEXACT_POPULATIONS_METRICS_AND_HISTOGRAMS_PASS");
		}
	}
	private static int gcControlOrdinal(final LineParser1 p, final boolean allData){
		for(int bins=0; bins<2; bins++){
			for(int arm=0; arm<2; arm++){
				final String name=allData ? BuildNcrnaDdlLayoutCounts.gcAllMethod(bins==0 ? 512 : 1024, arm*BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS.length)
					: BuildNcrnaDdlLayoutCounts.gcMethod(bins==0 ? 512 : 1024, arm*BuildNcrnaDdlLayoutCounts.GC_SETTINGS.length);
				if(p.termEquals(name, 3)){return bins*2+arm;}
			}
		}
		return -1;
	}

	/** GC workers contain their own matched w=0 controls at both decided widths.
	 * Rank509 is an explicit filter failure, pooled before the percentile. */
	private void processGc(final String roster, final String prefix){
		readGcRoster(roster, false);
		if(gcFollowup!=null){readGcRoster(gcFollowup, true);}
		for(int a=0; a<=assignments.length; a++){
			for(int s=0; s<seeds.length; s++){
				for(int stratum=0; stratum<(a==assignments.length ? 1 : 3); stratum++){
					for(int method=0; method<methodCount; method++){
						final int cell=method+stratum*methodCount;
						require(masks[cell][a][s]==(a==assignments.length ? 3 : 31), "GC settings require every distinct judging fold/population");
						for(int metric=9; metric<=10; metric++){
							require(values[cell][a][s][metric]==values[stratum*methodCount][a][s][metric]
								&& values[cell][a][s][metric]==values[stratum*methodCount][a][0][metric]
								&& (a==assignments.length || values[cell][a][s][metric]==values[stratum*methodCount][0][0][metric]),
								"GC changes cannot remove difficult loci from the original population or low-identity denominator");
						}
						if(a<assignments.length && stratum==0){
							for(final int metric:new int[]{0,1,2,3,4,5,9,10,18,19}){
								require(values[method][a][s][metric]==values[method+methodCount][a][s][metric]+values[method+2*methodCount][a][s][metric],
									"GC contributor/noncontributor counts,failures and ranks partition the complete judging population");
							}
						}
					}
				}
			}
		}
		write(prefix);
	}
	private void readGcRoster(final String roster, final boolean extension){
		final boolean[][] workers=new boolean[assignments.length+1][seeds.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(roster)){
			RrnaResourceIO.header(in.reader, "assignment_seed\tddl_seed\ttable", roster);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==3, "GC roster binds one assignment,hash seed and completed worker table");
				final long assignment=p.parseLong(0), seed=p.parseLong(1);
				final int a=assignment==-1 ? assignments.length : index(assignments, assignment), s=index(seeds, seed);
				require(assignment>=-1 && !workers[a][s], "GC roster cannot omit or duplicate matched seed/assignment workers");
				workers[a][s]=true;
				readGc(p.parseString(2), a, s, assignment, seed, extension);
			}
		}
		for(int a=0; a<workers.length; a++){
			for(int s=0; s<seeds.length; s++){
				require(workers[a][s], "Missing GC worker from the predeclared complete crossed roster");
			}
		}
	}
	private void readGc(final String path, final int a, final int s, final long assignment, final long seed, final boolean extension){
		final int originalMethods=4*gcVariants, start=extension ? originalMethods : 0;
		final int end=extension ? methodCount : originalMethods;
		for(final long[][] subset:pooledHistograms){for(final long[] histogram:subset){Arrays.fill(histogram, 0);}}
		final LineParser1 p=new LineParser1('\t');
		final long[] checks=new long[4];
		final int folds=assignment<0 ? 2 : 5, strata=assignment<0 ? 1 : 3;
		int rows=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==39 && p.parseLong(0)==assignment && p.parseLong(1)==seed,
					"GC rows retain exact assignment/hashseed plus ranks,failure bins and acceptable-model losses");
				final int method=methods.get(p.parseString(3)), fold=p.parseInt(2), stratum=contribution(p);
				require(method>=start && method<end && fold>=0 && fold<folds && stratum<strata,
					"Original and extension GC workers own disjoint method ranges and cannot duplicate saved controls");
				binding(p, method);
				final int index=method+methodCount*stratum;
				require((masks[index][a][s]&(1<<fold))==0, "One GC fold/population can only contribute once to its cell");
				masks[index][a][s]|=1<<fold;
				final double[] cell=values[index][a][s];
				if(Double.isNaN(cell[0])){Arrays.fill(cell, 0);}
				for(int metric=0; metric<6; metric++){
					final long value=p.parseLong(14+metric);
					require(value>=0, "GC rank/count totals must be nonnegative integers");
					cell[metric]+=value;
				}
				cell[9]+=p.parseLong(12);
				cell[10]+=p.parseLong(13);
				for(int metric=0; metric<9; metric++){
					final long value=p.parseLong(30+metric);
					require(value>=0 && (stratum==0 || value==0), "Window/model-pair counts belong to ALL only and cannot be duplicated across parent strata");
					cell[metric<5 ? 13+metric : 15+metric]+=value;
				}
				require(p.parseLong(31)<=p.parseLong(30) && p.parseLong(33)<=p.parseLong(32), "Dropped models/windows cannot exceed their original acceptable-model denominators");
				for(int low=0; low<2; low++){
					p.setBounds(28+low);
					addHistogram(row, p.a(), p.b(), pooledHistograms[index][low], checks);
					require(checks[0]==p.parseLong(12+low) && checks[1]==p.parseLong(14+low)
						&& checks[2]==p.parseLong(16+2*low) && checks[3]==p.parseLong(17+2*low),
						"GC histograms including the failure sentinel must conserve all original populations and rank/count totals");
				}
				rows++;
			}
		}
		require(rows==(end-start)*folds*strata, "Every GC worker contains all of its settings/folds/populations exactly once");
		for(int stratum=0; stratum<strata; stratum++){
			for(int method=start; method<end; method++){
				final int index=method+stratum*methodCount;
				final long[][] histograms=pooledHistograms[index];
				final double[] cell=values[index][a][s];
				cell[18]=histograms[0][508];
				cell[19]=histograms[1][508];
				if(stratum==0){
					cell[11]=percentile99(histograms[0]);
					cell[12]=percentile99(histograms[1]);
					final int control=gcControl(method);
					require(cell[13]==values[control][a][s][13] && cell[15]==values[control][a][s][15]
						&& cell[17]==values[control][a][s][17], "GC alternatives preserve the same judged-window and acceptable-model-pair census");
				}
				if(method<originalMethods && method%gcVariants<6){require(cell[14]==0 && cell[16]==0 && cell[18]==0 && cell[19]==0,
					"GC sorting never filters an acceptable model, even in the zero-SD bottom tier");}
				if(strata==3 && stratum==0){
					for(int low=0; low<2; low++){
						for(int rank=0; rank<509; rank++){
							require(histograms[low][rank]==pooledHistograms[method+methodCount][low][rank]+pooledHistograms[method+2*methodCount][low][rank],
								"GC contribution histograms partition every rank, including all-filtered failures");
						}
					}
				}
			}
		}
	}
	private int gcControl(final int method){
		assert(gc && method>=0 && method<methodCount) : "Each extension setting reuses the accepted same-width window w0 control";
		return method<4*gcVariants ? gcVariants*(method/gcVariants) : spec[method][1]==512 ? 0 : 2*gcVariants;
	}

	private void process(final String roster, final String prefix){
		final String[] outputs=sizegrid==null ? new String[]{prefix+".summary.tsv",prefix+".paired.tsv"}
			: new String[]{prefix+".summary.tsv",prefix+".paired.tsv",prefix+".sizes.tsv"};
		require(Tools.testInputFiles(false, true, roster)
			&& Tools.testOutputFiles(false, false, false, outputs)
			&& Tools.testForDuplicateFiles(true, roster, prefix+".summary.tsv", prefix+".paired.tsv", prefix+".sizes.tsv"),
			"The complete worker roster must be readable and report files fresh/distinct");
		if(gc){processGc(roster, prefix); return;}
		final boolean[][] workers=new boolean[assignments.length+1][seeds.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(roster)){
			RrnaResourceIO.header(in.reader, "assignment_seed\tddl_seed\ttable", roster);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==3, "Each completed worker needs its declared assignment,hashseed,and native result table");
				final long assignment=p.parseLong(0), seed=p.parseLong(1);
				final int a=assignment<0 ? assignments.length : index(assignments, assignment), s=index(seeds, seed);
				require(assignment>=-1 && !workers[a][s], "Worker roster cannot duplicate a seed/assignment cell");
				workers[a][s]=true;
				read(p.parseString(2), a, s, assignment, seed);
			}
		}
		if(followup!=null){readFollowupRoster(followup, 0);}
		if(sizegrid!=null){readFollowupRoster(sizegrid, 1);}
		if(sizeextension!=null){readFollowupRoster(sizeextension, 3);}
		if(definitions!=null){readFollowupRoster(definitions, 2);}
		for(int a=0; a<workers.length; a++){
			for(int s=0; s<seeds.length; s++){require(workers[a][s], "Missing worker from the full crossed seed/assignment roster");}
			for(int m=0; m<(a==assignments.length ? methodCount : 3*methodCount); m++){
				final int method=m%methodCount, stratum=m/methodCount;
				for(int s=0; s<(method==49 ? 1 : seeds.length); s++){
					require(masks[m][a][s]==(a==assignments.length ? 3 : 31), "Every method/population/cell must contain all distinct judging folds: "+names[method]);
					for(int metric=6; metric<9; metric++){values[m][a][s][metric]/=a==assignments.length ? 2 : 5;}
					require(values[m][a][s][9]==values[stratum*methodCount][a][0][9] && values[m][a][s][10]==values[stratum*methodCount][a][0][10],
						"Every method/seed must judge the identical complete parent and low-identity populations");
					if(a<assignments.length){
						require(values[m][a][s][9]==values[stratum*methodCount][0][0][9] && values[m][a][s][10]==values[stratum*methodCount][0][0][10],
							"Pooling disjoint folds must conserve all comparable loci across assignment replicates");
						for(int metric=0; metric<11; metric++){
							if(metric>=6 && metric<9){continue;}
							require(values[method][a][s][metric]==values[method+methodCount][a][s][metric]+values[method+2*methodCount][a][s][metric],
								"Original contributor/noncontributor populations must partition every pooled locus and rank/count metric");
							if(method==49){require(values[m][a][0][metric]==values[m][0][0][metric],
								"The unseeded index reference must pool to the same loci and ranks across assignments");}
						}
					}
				}
			}
		}
		write(prefix);
	}

	private void read(final String path, final int a, final int s, final long assignment, final long seed){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			final byte[] header=in.reader.nextLine();
			require(header!=null, "A completed replicate table needs its explicit schema");
			final String schema=new String(header, java.nio.charset.StandardCharsets.UTF_8);
			final boolean legacy=schema.equals(BuildNcrnaDdlLayoutCounts.LEGACY_REPLICATE_HEADER);
			require(schema.equals(BuildNcrnaDdlLayoutCounts.REPLICATE_HEADER) || legacy && assignment==-1,
				"Only accepted fixed-panel observations may use the legacy schema without contribution strata");
			int rows=0;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==(legacy ? 27 : 28) && p.parseLong(0)==assignment, "Every replicate row must bind its complete declared assignment");
				final String name=p.parseString(3);
				final int m=methods.get(name), fold=p.parseInt(2);
				require(m>=0 && fold>=0 && fold<(assignment<0 ? 2 : 5), "A worker cannot introduce an undeclared method or judging fold");
				binding(p, m);
				int stratum=0;
				if(!legacy){
					stratum=-1;
					for(int group=0; group<3; group++){
						if(p.termEquals(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[group], 27)){stratum=group;}
					}
					require(stratum>=0 && (assignment>=0 || stratum==0), "Fixed observations retain ALL only; rebuilt folds require named contribution strata");
				}
				require(m<50, "Original completed observations cannot introduce follow-up methods");
				final int cellIndex=m+methodCount*stratum;
				require(stratum>0 || p.parseInt(12)>0, "The complete judging fold must be nonempty; only contribution subsets may be empty");
				final int slot=m==49 ? 0 : s;
				require(m==49 ? p.parseLong(1)==-1 && s==0 : p.parseLong(1)==seed,
					"Only the first seed worker emits the unseeded consensus reference;all other rows retain their actual hashseed");
				require((masks[cellIndex][a][slot]&(1<<fold))==0, "A judging fold/population cannot be counted twice as another replicate");
				masks[cellIndex][a][slot]|=1<<fold;
				final double[] cell=values[cellIndex][a][slot];
				if(Double.isNaN(cell[0])){Arrays.fill(cell, 0);}
				for(int metric=0; metric<6; metric++){
					final long value=p.parseLong(14+metric);
					require(value>=0, "Rank/count observations must retain exact nonnegative integers");
					cell[metric]+=value;
				}
				for(int third=0; third<3; third++){
					final double value=p.termEquals("NA", 21+third) ? Double.NaN : p.parseDouble(21+third);
					require(m==49 || stratum>0 ? Double.isNaN(value) : Double.isFinite(value) && value>=0 && value<=1,
						"Consensus retention is a measured fraction for each size third,not an empty-bin match");
					cell[6+third]+=value;
				}
				cell[9]+=p.parseInt(12);
				cell[10]+=p.parseInt(13);
				rows++;
			}
			require(rows==(49+(s==0 ? 1 : 0))*(assignment<0 ? 2 : 15), "A completed worker must contain every released method, fold and contribution stratum exactly once");
		}
	}

	/** The follow-up roster also names each matching original observation.
	 * Controls must reproduce exact per-fold counts; only histograms are new. */
	private void readFollowupRoster(final String roster, final int kind){
		final boolean[][] workers=new boolean[assignments.length+1][seeds.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(roster)){
			RrnaResourceIO.header(in.reader, "assignment_seed\tddl_seed\toriginal_table\ttable", roster);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==4, "Follow-up workers must name their exact original paired observation");
				final long assignment=p.parseLong(0), seed=p.parseLong(1);
				final int a=assignment==-1 ? assignments.length : index(assignments, assignment), s=index(seeds, seed);
				require(assignment>=-1 && !workers[a][s], "Follow-up cannot duplicate a seed/assignment observation");
				workers[a][s]=true;
				readFollowup(p.parseString(3), p.parseString(2), a, s, assignment, seed, kind);
			}
		}
		for(final boolean[] row:workers){for(final boolean present:row){require(present, "The follow-up needs every original seed/assignment cell");}}
	}

	private void readFollowup(final String path, final String original, final int a, final int s,
			final long assignment, final long seed, final int kind){
		final int folds=assignment<0 ? 2 : 5, strata=assignment<0 ? 1 : 3;
		final long[][][][] anchors=new long[3][strata][folds][];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(original)){
			final byte[] header=in.reader.nextLine();
			require(header!=null, "Matched controls require an actual original table");
			final String schema=new String(header, java.nio.charset.StandardCharsets.UTF_8);
			final boolean legacy=schema.equals(BuildNcrnaDdlLayoutCounts.LEGACY_REPLICATE_HEADER);
			require(schema.equals(BuildNcrnaDdlLayoutCounts.REPLICATE_HEADER) || legacy && assignment==-1,
				"Only the original fixed check may lack contribution strata");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int method=methods.get(p.parseString(3));
				if(method!=29 && method!=35 && method!=34){continue;}
				binding(p, method);
				final int fold=p.parseInt(2), stratum=legacy ? 0 : contribution(p);
				require(p.parseLong(0)==assignment && p.parseLong(1)==seed && fold>=0 && fold<folds
					&& stratum<strata && anchors[controlIndex(method)][stratum][fold]==null,
					"An original bonus100 control must match the exact seed,assignment,fold and population once");
				final long[] counts=new long[8];
				for(int i=0; i<counts.length; i++){counts[i]=p.parseLong(12+i);}
				anchors[controlIndex(method)][stratum][fold]=counts;
			}
		}
		for(final long[][][] control:anchors){for(final long[][] subset:control){for(final long[] fold:subset){
			require(fold!=null, "Both original bonus100 controls must retain every judging fold/population");
		}}}
		final long[][][] histograms=pooledHistograms;
		for(final long[][] subset:histograms){for(final long[] histogram:subset){Arrays.fill(histogram, 0);}}
		final int[] seen=new int[3*methodCount];
		final long[] histogramMetrics=new long[4];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==30 && p.parseLong(0)==assignment && p.parseLong(1)==seed,
					"Every follow-up row binds its seed/assignment and exact rank distributions");
				final int method=methods.get(p.parseString(3)), fold=p.parseInt(2), stratum=contribution(p);
				require(inRound(method, kind) && fold>=0 && fold<folds && stratum<strata,
					"Each job family contains only its declared settings;overlapping controls are imported once");
				binding(p, method);
				final int cellIndex=method+methodCount*stratum;
				require((seen[cellIndex]&(1<<fold))==0, "A follow-up judging fold/population cannot occur twice");
				seen[cellIndex]|=1<<fold;
				final double[] cell=values[cellIndex][a][s];
				if(method<50){
					final long[] anchor=anchors[controlIndex(method)][stratum][fold];
					for(int i=0; i<anchor.length; i++){
						require(anchor[i]==p.parseLong(12+i), "Repeated bonus100 control changed a saved population/count/rank metric: "+names[method]+", field="+(12+i));
					}
				}else{
					if(Double.isNaN(cell[0])){Arrays.fill(cell, 0);}
					for(int metric=0; metric<6; metric++){cell[metric]+=p.parseLong(14+metric);}
					for(int third=0; third<3; third++){
						final double value=p.termEquals("NA", 21+third) ? Double.NaN : p.parseDouble(21+third);
						require(stratum>0 ? Double.isNaN(value) : Double.isFinite(value) && value>=0 && value<=1,
							"New consensus retention must be measured for ALL and unavailable for contribution subsets");
						cell[6+third]+=value;
					}
					cell[9]+=p.parseInt(12);
					cell[10]+=p.parseInt(13);
					masks[cellIndex][a][s]|=1<<fold;
				}
				for(int low=0; low<2; low++){
					p.setBounds(28+low);
					addHistogram(row, p.a(), p.b(), histograms[cellIndex][low], histogramMetrics);
					require(histogramMetrics[0]==p.parseLong(12+low) && histogramMetrics[1]==p.parseLong(14+low)
						&& histogramMetrics[2]==p.parseLong(16+2*low) && histogramMetrics[3]==p.parseLong(17+2*low),
						"Exact rank distributions must reproduce each fold's population,rank sum and top10/50 counts");
				}
			}
		}
		for(int stratum=0; stratum<strata; stratum++){
			for(int method=0; method<methodCount; method++){
				if(!inRound(method, kind)){continue;}
				final int cellIndex=method+methodCount*stratum;
				require(seen[cellIndex]==(1<<folds)-1, "Every follow-up setting needs all disjoint judging folds");
				if(stratum==0){
					for(int low=0; low<2; low++){values[cellIndex][a][s][11+low]=percentile99(histograms[cellIndex][low]);}
				}
				if(strata==3){
					for(int low=0; low<2; low++){
						for(int i=0; i<508; i++){
							require(histograms[method][low][i]==histograms[method+methodCount][low][i]+histograms[method+2*methodCount][low][i],
								"Contributor histograms must partition every pooled rank bin,not just aggregate sums");
						}
					}
				}
			}
		}
	}
	private static int controlIndex(final int method){
		assert(method==29 || method==35 || method==34) : "Only the two bonus100 controls and unfiltered winner are explicitly repeated";
		return method==29 ? 0 : method==35 ? 1 : 2;
	}
	private boolean inRound(final int method, final int kind){
		assert(kind>=0 && kind<=3) : "Each roster is nearby bonuses,original grid,selected definitions or the22-cell extension";
		return kind==0 ? method==29 || method==35 || method>=50 && method<68
			: kind==1 ? method==34 || method>=68 && method<138
			: kind==3 ? method>=138 && method<160 : method>=definitionStart && method<definitionStart+3;
	}

	private static int contribution(final LineParser1 p){
		for(int group=0; group<3; group++){if(p.termEquals(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[group], 27)){return group;}}
		throw new IllegalArgumentException("Follow-up requires the declared original assembly-contribution population");
	}

	/** Parses a comma-separated byte range into caller-owned pooled bins/scratch. */
	public static void addHistogram(final byte[] row, final int start, final int stop,
			final long[] pooled, final long[] metrics){
		require(row!=null && start>=0 && stop>start && stop<=row.length && pooled!=null
			&& pooled.length>0 && metrics!=null && metrics.length==4,
			"Rank histogram pooling needs one complete byte range,rank bins and count/sum/top10/top50 scratch");
		Arrays.fill(metrics, 0);
		int from=start, rank=0;
		for(int end=start; end<=stop; end++){
			if(end<stop && row[end]!=','){continue;}
			require(end>from && rank<pooled.length, "A rank histogram cannot contain empty or excess bins");
			final long count=Parse.parseLong(row, from, end);
			require(count>=0, "Rank histogram counts cannot be negative");
			pooled[rank]=Math.addExact(pooled[rank], count);
			metrics[0]=Math.addExact(metrics[0], count);
			metrics[1]=Math.addExact(metrics[1], Math.multiplyExact(rank+1L, count));
			if(rank<10){metrics[2]=Math.addExact(metrics[2], count);}
			if(rank<50){metrics[3]=Math.addExact(metrics[3], count);}
			rank++;
			from=end+1;
		}
		require(rank==pooled.length, "Rank histograms require exactly one bin for every original model rank");
	}

	/** Histogram index zero represents rank one. Pool folds before calling. */
	public static int percentile99(final long[] counts){
		require(counts!=null && counts.length>0, "Nearest-rank percentile requires an exact rank histogram");
		long total=0;
		for(final long count:counts){require(count>=0, "Rank histogram counts cannot be negative"); total=Math.addExact(total, count);}
		require(total>0, "Nearest-rank percentile needs a nonempty population");
		final long target=total-total/100;
		long cumulative=0;
		for(int i=0; i<counts.length; i++){
			cumulative+=counts[i];
			if(cumulative>=target){return i+1;}
		}
		throw new IllegalStateException("Rank histogram failed to conserve its own population");
	}

	private boolean hasMetric(final int method, final int stratum, final int metric){
		if(gc){return metric<6 || metric>=11 && (stratum==0 || metric==18 || metric==19);}
		if(metric<6){return true;}
		if(metric<9){return method!=49 && stratum==0;}
		return metric>=11 && followup!=null && stratum==0 && (method==29 || method==35 || method>=50 || sizegrid!=null && method==34);
	}

	private void binding(final LineParser1 p, final int method){
		final int[] shape=spec[method];
		if(gc){
			require(p.parseInt(4)==8 && p.parseInt(5)==10 && p.termEquals("quarter", 6)
				&& p.parseInt(7)==shape[1] && p.parseInt(8)==256 && p.parseInt(9)==2048
				&& p.parseInt(10)==shape[4] && p.parseDouble(11)==0,
				"GC rows must preserve the decided words, quarter layout, blacklist, bins and their corresponding bonus");
			return;
		}
		final int pair=shape[0], shortK=pair>=2 ? 8 : 7, longK=pair==0 ? 9 : pair==3 ? 11 : 10;
		final int bonus=method<36 || method==49 ? 100 : shape[4]<0 ? 0 : shape[4];
		final double scale=method>=39 && method<49 ? -shape[4]/10.0 : 0;
		require(p.parseInt(4)==shortK && p.parseInt(5)==longK && p.termEquals(LAYOUTS[shape[2]], 6)
			&& p.parseInt(7)==shape[1] && p.parseInt(8)==(shape.length>5 ? shape[5] : shape[3]==0 ? 0 : pair>=2 ? 256 : 64)
			&& p.parseInt(9)==(shape.length>5 ? shape[6] : shape[3]==0 ? 0 : 4096) && p.parseInt(10)==bonus
			// Parse.parseDouble multiplies by decimalInvMult;allow the writer's nine-place rounding.
			&& Math.abs(p.parseDouble(11)-scale)<=0.000000000500001,
			"Method identity must agree with actual word lengths,layout,width,blacklist and weighting: "+names[method]);
		require(p.parseInt(12)>=0 && p.parseInt(13)>=0 && p.parseInt(13)<=p.parseInt(12)
			&& p.parseInt(20)>=0 && Double.isFinite(p.parseDouble(24)) && p.parseDouble(24)>=0
			&& Double.isFinite(p.parseDouble(25)) && p.parseDouble(25)>=0,
			"A completed fold needs measured nonempty populations,empty-model accounting and finite nonnegative timing");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reports            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeSizes(final String prefix){
		final int[][] sizes=sizeextension==null ? BuildNcrnaDdlLayoutCounts.SIZE_GRID : BuildNcrnaDdlLayoutCounts.EXPANDED_SIZE_GRID;
		final int count=sizes.length;
		final int[] indexes=new int[count], shortSizes=new int[count], longSizes=new int[count];
		final double[][][] ranks=new double[count][][];
		int next=0;
		for(final int[] size:sizes){
			final int shortSize=size[0], longSize=size[1];
			final int method=shortSize==0 && longSize==0 ? 34 : shortSize==256 && longSize==4096 ? 35
				: methods.get(BuildNcrnaDdlLayoutCounts.sizeMethod(shortSize, longSize));
			require(method>=0, "Each declared size cell must have one measured or imported method");
			indexes[next]=method;
			shortSizes[next]=shortSize;
			longSizes[next]=longSize;
			ranks[next]=cells(method, -1, 0, false);
			next++;
		}
		final int[] choice=chooseSize(ranks, shortSizes, longSizes);
		chosenSizeMethod=indexes[choice[1]];
		require(definitions==null || spec[definitionStart][5]==shortSizes[choice[1]] && spec[definitionStart][6]==longSizes[choice[1]],
			"Alternative definitions must use the size selected from the full-population grid by the predeclared one-paired-SE rule");
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+".sizes.tsv")){
			out.writer.println("cohort\tshort_blacklist\tlong_blacklist\trank_mean\trank_sd\tlow_rank_mean\tlow_rank_sd\tp99_mean\tp99_sd\tlow_p99_mean\tlow_p99_sd\trank_difference_from_full_best\tpaired_standard_error\twithin_one_se_full\tbest_full\tchosen_full");
			for(int fixed=0; fixed<2; fixed++){
				for(int i=0; i<count; i++){
					final double[] fullDifference=moments(difference(ranks[i], ranks[choice[0]]));
					b.clear().append(fixed==0 ? "REBUILT_FIVE_FOLD" : "FIXED_NOT_SEEN").tab().append(shortSizes[i]).tab().append(longSizes[i]);
					for(final int metric:new int[]{0,1,11,12}){
						final double[] stats=moments(cells(indexes[i], -1, metric, fixed!=0));
						b.tab().append(stats[0], 9).tab().append(stats[1], 9);
					}
					final double[] paired=moments(cells(indexes[i], indexes[choice[0]], 0, fixed!=0));
					out.writer.print(b.tab().append(paired[0], 9).tab().append(paired[5], 9)
						.tab().append(fullDifference[0]<=fullDifference[5]).tab().append(i==choice[0]).tab().append(i==choice[1]).nl());
				}
			}
		}
	}

	/** Best full mean defines the reference. Choose the fewest entries within
	 * one paired SE;equal total sizes prefer fewer short-word entries. */
	public static int[] chooseSize(final double[][][] ranks, final int[] shortSizes, final int[] longSizes){
		require(ranks!=null && ranks.length>0 && shortSizes!=null && longSizes!=null
			&& shortSizes.length==ranks.length && longSizes.length==ranks.length,
			"Size selection requires every candidate's matched cells and literal prefix lengths");
		final double[] means=new double[ranks.length];
		int best=0;
		for(int i=0; i<ranks.length; i++){
			require(shortSizes[i]>=0 && longSizes[i]>=0, "Blacklist entry counts cannot be negative");
			means[i]=moments(ranks[i])[0];
			if(means[i]<means[best] || means[i]==means[best] && smaller(i, best, shortSizes, longSizes)){best=i;}
		}
		int chosen=best;
		for(int i=0; i<ranks.length; i++){
			final double[] paired=moments(difference(ranks[i], ranks[best]));
			if(paired[0]<=paired[5] && smaller(i, chosen, shortSizes, longSizes)){chosen=i;}
		}
		return new int[]{best,chosen};
	}
	private static boolean smaller(final int left, final int right, final int[] shortSizes, final int[] longSizes){
		final long a=(long)shortSizes[left]+longSizes[left], b=(long)shortSizes[right]+longSizes[right];
		return a<b || a==b && shortSizes[left]<shortSizes[right];
	}
	private static double[][] difference(final double[][] left, final double[][] right){
		require(left!=null && right!=null && left.length==right.length, "Paired cells must use identical assignment rosters");
		final double[][] result=new double[left.length][];
		for(int a=0; a<left.length; a++){
			require(left[a]!=null && right[a]!=null && left[a].length==right[a].length,
				"Paired cells must retain every matched hash seed within each assignment");
			result[a]=new double[left[a].length];
			for(int s=0; s<result[a].length; s++){result[a][s]=left[a][s]-right[a][s];}
		}
		return result;
	}

	private void write(final String prefix){
		if(sizegrid!=null){writeSizes(prefix);}
		final ByteBuilder b=new ByteBuilder();
		try(Output summary=new Output(prefix+".summary.tsv"); Output paired=new Output(prefix+".paired.tsv")){
			summary.writer.println("cohort\tmethod\tmetric\tassignments\thash_seeds\tmean\tstandard_deviation\thash_variance\tsplit_variance\tinteraction_variance\tmean_standard_error\tparents\tlow_parents\tassembly_contribution");
			paired.writer.println("cohort\tleft\tright\tmetric\tassignments\thash_seeds\tmean_paired_difference\tstandard_error\thash_variance\tsplit_variance\tinteraction_variance\tassembly_contribution");
			for(int fixed=0; fixed<2; fixed++){
				final String cohort=fixed==0 ? "REBUILT_FIVE_FOLD" : "FIXED_NOT_SEEN";
				for(int stratum=0; stratum<(fixed==0 ? 3 : 1); stratum++){
					for(int m=0; m<methodCount; m++){
						final int cell=m+methodCount*stratum;
						for(int metric=0; metric<(gc ? METRICS.length : 13); metric++){
							if(!hasMetric(m, stratum, metric)){continue;}
							final double[][] cells=cells(cell, -1, metric, fixed!=0);
							final double[] stats=moments(cells);
							b.clear().append(cohort).tab().append(names[m]).tab().append(METRICS[metric]).tab().append(cells.length).tab().append(cells[0].length);
							for(final double value:stats){b.tab().append(value, 9);}
							final double[] population=values[cell][fixed==0 ? 0 : assignments.length][0];
							summary.writer.print(b.tab().append((long)population[9]).tab().append((long)population[10])
								.tab().append(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[stratum]).nl());
						}
						for(int other=m+1; other<methodCount; other++){
							if(!paired(m, other)){continue;}
							for(int metric=0; metric<(gc ? METRICS.length : 13); metric++){
								if(metric>=6 && metric<11 || !hasMetric(m, stratum, metric) || !hasMetric(other, stratum, metric)){continue;}
								final double[][] cells=cells(cell, other+methodCount*stratum, metric, fixed!=0);
								final double[] stats=moments(cells);
								b.clear().append(cohort).tab().append(names[m]).tab().append(names[other]).tab().append(METRICS[metric])
									.tab().append(cells.length).tab().append(cells[0].length).tab().append(stats[0], 9).tab().append(stats[5], 9);
								for(int i=2; i<=4; i++){b.tab().append(stats[i], 9);}
								paired.writer.print(b.tab().append(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[stratum]).nl());
							}
						}
					}
				}
			}
		}
	}

	private double[][] cells(final int m, final int other, final int metric, final boolean fixed){
		final boolean unseeded=!gc && m%methodCount==49, otherUnseeded=!gc && other>=0 && other%methodCount==49;
		final double[][] result=new double[fixed ? 1 : assignments.length][unseeded && other<0 ? 1 : seeds.length];
		for(int a=0; a<result.length; a++){
			for(int s=0; s<result[a].length; s++){
				final int assignment=fixed ? assignments.length : a;
				result[a][s]=values[m][assignment][unseeded ? 0 : s][metric]
					-(other<0 ? 0 : values[other][assignment][otherUnseeded ? 0 : s][metric]);
			}
		}
		return result;
	}

	private boolean paired(final int left, final int right){
		if(gc){
			if(right<4*gcVariants){return left/gcVariants==right/gcVariants && left%gcVariants==0;}
			final int variant=(right-4*gcVariants)%BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length;
			final int control=gcControl(right);
			return left==control || left==control+12+variant/3;
		}
		if(right>=definitionStart){return left==35 || left==chosenSizeMethod || left>=definitionStart;}
		if(right>=68){return left==35 || left==34 || left>=definitionStart;}
		if(right>=50){
			if(left==49){return true;}
			if(left>=36 && left<50){return false;}
			int differences=0;
			for(int axis=0; axis<5; axis++){
				final int l=axis==4 && left<36 ? 100 : spec[left][axis];
				if(l!=spec[right][axis]){differences++;}
			}
			return differences==1;
		}
		if(right==49){return true;}
		if(right>=36){return left==0;}
		int differences=0;
		for(int axis=0; axis<4; axis++){if(spec[left][axis]!=spec[right][axis]){differences++;}}
		return differences==1 && !(spec[left][2]!=spec[right][2] && spec[left][2]>0 && spec[right][2]>0);
	}

	/** Balanced crossed random-effects ANOVA,one observation per seed/assignment.
	 * Negative estimated components are clipped at zero;interaction is retained.
	 * Returns mean,sample SD,hash variance,split variance,interaction variance,SE. */
	public static double[] moments(final double[][] cells){
		require(cells!=null && cells.length>0 && cells[0]!=null && cells[0].length>0,
			"Crossed moments require a nonempty rectangular seed/assignment grid");
		final int a=cells.length, s=cells[0].length;
		final double[] row=new double[a], column=new double[s];
		double mean=0;
		for(int i=0; i<a; i++){
			require(cells[i]!=null && cells[i].length==s, "All assignment replicates must use the same hashseed set");
			for(int j=0; j<s; j++){
				require(Double.isFinite(cells[i][j]), "No missing or nonfinite cell may silently reduce the crossed denominator");
				mean+=cells[i][j];
				row[i]+=cells[i][j]/s;
				column[j]+=cells[i][j]/a;
			}
		}
		mean/=a*(double)s;
		double total=0, split=0, hash=0, interaction=0;
		for(int i=0; i<a; i++){
			split+=s*square(row[i]-mean);
			for(int j=0; j<s; j++){
				total+=square(cells[i][j]-mean);
				interaction+=square(cells[i][j]-row[i]-column[j]+mean);
			}
		}
		for(int j=0; j<s; j++){hash+=a*square(column[j]-mean);}
		final double residual=a>1 && s>1 ? interaction/((a-1.0)*(s-1.0)) : 0;
		final double hashVariance=s>1 ? Math.max(0, (hash/(s-1)-residual)/a) : 0;
		final double splitVariance=a>1 ? Math.max(0, (split/(a-1)-residual)/s) : 0;
		final double se=Math.sqrt(hashVariance/s+splitVariance/a+residual/(a*(double)s));
		return new double[]{mean,a*s>1 ? Math.sqrt(total/(a*(double)s-1)) : 0,hashVariance,splitVariance,residual,se};
	}
	private static double square(final double value){return value*value;}
	private static long[] numbers(final String text){
		require(text!=null && !text.isEmpty(), "Explicit replicate seed lists cannot be empty");
		final String[] terms=text.split(",", -1);
		final long[] values=new long[terms.length];
		for(int i=0; i<terms.length; i++){
			values[i]=Long.parseLong(terms[i]);
			require(values[i]>=0, "Replicate seeds cannot use a random-seed sentinel");
			for(int j=0; j<i; j++){require(values[j]!=values[i], "A replicate seed cannot occur twice");}
		}
		return values;
	}
	private static int index(final long[] values, final long value){
		for(int i=0; i<values.length; i++){if(values[i]==value){return i;}}
		throw new IllegalArgumentException("Undeclared replicate seed: "+value);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final long[] seeds, assignments;
	private final String followup, sizegrid, sizeextension, definitions, gcFollowup;
	private final int methodCount, definitionStart, gcVariants;
	private final boolean gc, gcAll;
	private final String[] names=new String[163];
	private final int[][] spec=new int[163][];
	private final ObjectIntMap<String> methods=new ObjectIntMap<String>(String.class);
	private final double[][][][] values;
	private final int[][][] masks;
	private final long[][][] pooledHistograms;
	private int chosenSizeMethod=-1;
	private static final String[] LAYOUTS={"mixed","half","quarter"};
	private static final String[] METRICS={"rank_sum","low_rank_sum","top10","top50","low_top10","low_top50",
		"consensus_retained_small","consensus_retained_medium","consensus_retained_large","parents","low_parents",
		"rank_p99","low_rank_p99", "acceptable_model_window_pairs","dropped_acceptable_model_window_pairs",
		"windows_with_acceptable_models","windows_all_acceptable_filtered","missing_query_gc_windows",
		"parents_all_acceptable_filtered","low_parents_all_acceptable_filtered", "window_stat_fallback_model_fold_cases",
		"family_member_fallback_model_fold_cases","zero_sd_model_fold_cases","missing_effective_gc_model_fold_cases"};
}
