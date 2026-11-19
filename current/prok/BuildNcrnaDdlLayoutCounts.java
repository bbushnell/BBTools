package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import fileIO.ReadWrite;
import kmer.AbstractKmerTable;
import map.IntHashMap;
import map.LongHashSet;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Builds standard DDL reference objects by member voting and native consensus.
 * Query windows are native DDL objects; all scores use the shared DDL API.
 * Only transient reference-bin votes are selected here. Normal counts build
 * no index; optional parity/timing uses the existing small index with native
 * objects. No sketch codec, raw-probe path or custom container is used.
 * @author Brian Bushnell, Raiden, Ganyu
 */
public final class BuildNcrnaDdlLayoutCounts implements CardinalityTracker.KmerConsumer{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BuildNcrnaDdlLayoutCounts.class, false);
		try{
			final Config c=parse(pp.args);
			run(c);
			timer.stop();
			pp.outstream.println("Time: "+timer);
			pp.outstream.println(c.mergeSeeds!=null ? "NCRNA_DDL_CROSS_SEED_SUMMARY_COMPLETE_NOT_CALLER_ADOPTED"
				: c.mergeStudy!=null ? "NCRNA_DDL_STUDY_MERGE_COMPLETE_NOT_CALLER_ADOPTED"
				: c.study!=null ? "NCRNA_DDL_STUDY_COMPLETE_NOT_CALLER_ADOPTED"
				: c.checkIndex==null ? "NCRNA_DDL_LAYOUT_COUNTS_COMPLETE_NOT_FITTED_OR_ADOPTED"
				: "NCRNA_DDL_INDEX_CHECK_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private static Config parse(final String[] args){
		assert(args!=null) : "PreParser supplies the expanded native flag=value arguments";
		final Config c=new Config();
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("recipes")){c.recipes=b;}
			else if(a.equals("panels")){c.panels=b;}
			else if(a.equals("selection")){c.selection=b;}
			else if(a.equals("ref")){c.references=b;}
			else if(a.equals("members")){c.members=b;}
			else if(a.equals("ledger")){c.ledger=b;}
			else if(a.equals("assignments")){c.assignments=b;}
			else if(a.equals("census")){c.census=b;}
			else if(a.equals("summary")){c.summary=b;}
			else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
			else if(a.equals("checkindex")){c.checkIndex=b;}
			else if(a.equals("indexseconds")){c.indexSeconds=Parse.parseIntKMG(b);}
			else if(a.equals("study")){c.study=b;}
			else if(a.equals("gcfits")){c.gcFits=b;}
			else if(a.equals("gcallfit")){c.gcAllFit=b;}
			else if(a.equals("gcalllimits")){c.gcAllLimits=b;}
			else if(a.equals("gcadditions")){c.gcAdditions=Parse.parseBoolean(b);}
			else if(a.equals("identities")){c.identities=b;}
			else if(a.equals("mergestudy")){c.mergeStudy=b;}
			else if(a.equals("mergeseeds")){c.mergeSeeds=b;}
			else if(a.equals("seedbaseline")){c.seedBaseline=b;}
			else if(a.equals("sourcegroups")){c.sourceGroups=b;}
			else if(a.equals("genomes")){c.genomes=b;}
			else if(a.equals("foldseed")){c.foldSeed=Long.parseLong(b);}
			else if(a.equals("fixedcheck")){c.fixedCheck=Parse.parseBoolean(b);}
			else if(a.equals("indexreference")){c.indexReference=Parse.parseBoolean(b);}
			else if(a.equals("folds")){c.folds=Parse.parseIntKMG(b);}
			else if(a.equals("shortprefixes")){
				final String[] sizes=b.split(",", -1);
				require(sizes.length==3, "shortprefixes declares one padding prefix per released word pair");
				for(int i=0; i<3; i++){c.shortPrefixes[i]=Parse.parseIntKMG(sizes[i]);}
			}
			else if(a.equals("longprefix")){c.longPrefix=Parse.parseIntKMG(b);}
			else if(a.equals("ddlseed")){c.seed=Long.parseLong(b);}
			else if(a.equals("wordtag")){
				require("xor".equals(b) || "leading".equals(b), "wordtag=xor|leading names the pre-hash length namespace");
				c.oldTags=b.equals("xor");
			}
			else{require(parser.parse(arg, a, b), "Unknown native DDL count parameter: "+arg);}
		}
		require(parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1 && Shared.threads()==1,
			"Use the complete named member/window inputs, fresh output and t=1; no raw-profile input or sampling");
		c.output=parser.out1;
		return c;
	}

	public static void run(final Config c){
		if(c!=null && c.mergeSeeds!=null){
			mergeSeeds(c.mergeSeeds, c.seedBaseline, c.output);
			return;
		}
		if(c!=null && c.mergeStudy!=null){
			mergeStudy(c.mergeStudy, c.output);
			return;
		}
		require(c!=null && c.models>0 && c.models<=65536 && c.groups>0 && c.recipes!=null
			&& c.panels!=null && c.selection!=null && c.references!=null && c.members!=null
			&& c.ledger!=null && c.assignments!=null && c.census!=null && c.summary!=null && c.output!=null,
			"Native DDL counting needs the complete original model/member/window roster and declared recipes");
		require(c.study==null || (c.study.equals("baseline") || c.study.equals("blacklist") || c.study.equals("bestblacklist") || gcDescription(c) || gcAllLosses(c) || gcAllPower(c) || replicated(c) || seededStudy(c) || partitionStudy(c))
			&& c.identities!=null && c.checkIndex==null, "Studies need frozen identities and exclude index timing");
		require(c.seed>=0, "The declared DDL seed must be nonnegative,not a random-seed sentinel");
		final String[] ordinary={c.recipes, c.panels, c.selection, c.references, c.members, c.ledger, c.assignments, c.census, c.summary};
		final String[] inputs=c.study==null ? ordinary : Arrays.copyOf(ordinary, ordinary.length+3);
		if(c.study!=null){
			inputs[ordinary.length]=c.identities+".identity.tsv";
			inputs[ordinary.length+1]=c.identities+".groups.tsv";
			inputs[ordinary.length+2]=c.identities+".identity_bins.tsv";
		}
		require(c.indexSeconds==0 || c.indexSeconds>=10, "Query timing needs at least10CPU seconds per path; zero selects finite parity only");
		require(c.checkIndex!=null || c.indexSeconds==10, "indexseconds applies only with an explicit checkindex recipe");
		require(!replicated(c) || c.oldTags && c.foldSeed>=0 && c.sourceGroups!=null && c.genomes!=null
			&& c.folds>=2 && c.folds<=5 && c.longPrefix>=0 && c.longPrefix<=32768
			&& c.shortPrefixes!=null && c.shortPrefixes.length==3
			&& c.shortPrefixes[0]>=0 && c.shortPrefixes[0]<=1024 && c.shortPrefixes[1]>=0 && c.shortPrefixes[1]<=1024
			&& c.shortPrefixes[2]>=0 && c.shortPrefixes[2]<=16384
			&& Tools.testInputFiles(false, true, c.sourceGroups, c.genomes),
			"Replicated folds require resolved source groups/genomes,nonnegative foldseed,and the shared XOR word namespace");
		require(!gcLosses(c) || c.gcFits!=null && Tools.testInputFiles(false, true, c.gcFits),
			"GC loss diagnostics require the accepted worker roster and reuse its recorded assignments");
		require(!gcAllData(c) || gcAllLosses(c) || c.gcAllFit!=null && c.gcAllLimits!=null
			&& Tools.testInputFiles(false, true, c.gcAllFit, c.gcAllLimits),
			"All-data ranking requires the one shared full-precision GC fit and measured finite lossless limits");
		require(!c.gcAdditions || "replicategcall".equals(c.study), "gcadditions extends only the accepted all-data window-GC ranking");
		final String[] outputs=gcAllPower(c) ? gcPowerOutputs(c)
			: gcAllLosses(c) ? new String[]{c.output+".gc_all_fit.tsv", c.output+".gc_all_limits.tsv", c.output+".gc_all_losses.tsv"}
			: gcLosses(c) ? new String[]{c.output+".gc_curve.tsv", c.output+".gc_choices.tsv"}
			: gcDescription(c) ? new String[]{c.output+".gc_models.tsv", c.output+".gc_family.tsv"}
			: indexStudy(c) ? new String[]{c.output+".index.tsv"}
			: replicated(c) ? gcStudy(c) && !gcAllData(c) ? new String[]{c.output+".replicates.tsv", c.output+".folds.tsv", c.output+".cohorts.tsv", c.output+".models.tsv", c.output+".gc_fit.tsv"}
				: new String[]{c.output+".replicates.tsv", c.output+".folds.tsv", c.output+".cohorts.tsv", c.output+".models.tsv"}
			: c.study==null ? new String[]{c.output+".layouts.tsv", c.output+(c.checkIndex==null ? ".counts.tsv.gz" : ".index.tsv")}
			: c.study.equals("baseline") ? new String[]{c.output+".study.tsv"}
			: partitionStudy(c) ? new String[]{c.output+".study.tsv", c.output+(c.study.equals("partitionseeds") ? ".seeds.summary.tsv" : ".background_groups.tsv")}
			: seededStudy(c) ? new String[]{c.output+".study.tsv", c.output+".seeds.summary.tsv",c.output+".loci.tsv",c.output+".locus_ranges.tsv",c.output+".tails.tsv"}
			: new String[]{c.output+".study.tsv", c.output+".blacklist.tsv", c.output+".background_groups.tsv"};
		require(Tools.testInputFiles(false, true, inputs) && Tools.testOutputFiles(false, false, false, outputs),
			"All bound native-study inputs must be readable and both outputs fresh");
		final String[] files=Arrays.copyOf(inputs, inputs.length+outputs.length);
		System.arraycopy(outputs, 0, files, inputs.length, outputs.length);
		require(Tools.testForDuplicateFiles(true, files), "Native count inputs and outputs cannot alias");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && !DynamicDemiLog.blacklistExists(), "Native parity requires assertions and no unrequested global blacklist");
		final int exponent=DynamicDemiLog.exponentBits(), zip=ReadWrite.MAX_ZIP_THREADS();
		final long seed=CardinalityTracker.defaultSeed;
		final boolean zipSet=ReadWrite.SET_ZIP_THREADS, trim=Shared.TRIM_READ_DESCRIPTION;
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Native reference construction needs an actual CPU clock");
		final boolean enabled=clock.isThreadCpuTimeEnabled();
		try{
			DynamicDemiLog.setExponent(4);
			CardinalityTracker.defaultSeed=c.seed;
			ReadWrite.setZipThreads(1);
			Shared.TRIM_READ_DESCRIPTION=false;
			if(!enabled){clock.setThreadCpuTimeEnabled(true);}
			new BuildNcrnaDdlLayoutCounts(c, clock).process();
		}finally{
			DynamicDemiLog.setExponent(exponent);
			CardinalityTracker.defaultSeed=seed;
			ReadWrite.setZipThreads(zip);
			ReadWrite.SET_ZIP_THREADS=zipSet;
			Shared.TRIM_READ_DESCRIPTION=trim;
			if(!enabled){clock.setThreadCpuTimeEnabled(false);}
		}
	}

	private BuildNcrnaDdlLayoutCounts(final Config config, final ThreadMXBean clock_){
		c=config;
		clock=clock_;
		sketchSeed=c.seed;
		oldTags=c.oldTags;
		rankHistograms=histogramStudy(c) ? new int[3][2][c.models+(gcStudy(c) ? 2 : 1)] : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Process            ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final ArrayList<Row> recipes=readRecipes();
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{names.add(read.id); sequences.add(read.bases);});
		require(names.size()==c.models, "The reference FASTA must retain all original model ordinals");
		final String[] headers=names.toArray(new String[names.size()]);
		final byte[][] consensuses=sequences.toArray(new byte[sequences.size()][]);
		final NcrnaDdlFolds folds=replicated(c) ? new NcrnaDdlFolds(data, c.models, c.sourceGroups, c.genomes) : null;
		final NcrnaDdlMembers members=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments,
			c.census, c.summary, headers, consensuses, folds);
		if(gcDescription(c)){
			describeGc(members, headers, data);
			return;
		}
		if(gcAllLosses(c)){
			allDataGcLosses(members, headers, data);
			return;
		}
		if(gcAllPower(c)){
			allDataGcPower(members, headers, data);
			return;
		}
		if(replicated(c)){
			if(gcLosses(c)){
				gcLossCurves(members, headers, consensuses, data, folds);
				return;
			}
			replicates(recipes, members, headers, consensuses, data, folds);
			return;
		}
		if(c.study!=null){
			final String selected=c.study.equals("bestblacklist") ? BEST_PAIRED : WINNER;
			for(final Row row:recipes){
				if(row.id.equals(selected)){
					study(c.study.equals("variance1024") ? new Row(row, 1024) : row, members, headers, consensuses, data);
					return;
				}
			}
			throw new IllegalArgumentException("Study roster must contain its exact selected512-bin recipe: "+selected);
		}
		try(Output metadata=new Output(c.output+".layouts.tsv");
				Output counts=c.checkIndex==null ? new Output(c.output+".counts.tsv.gz") : null){
			metadata.writer.println(NcrnaDdlLayoutCounts.LAYOUT_HEADER_V2);
			if(counts!=null){counts.writer.println(NcrnaDdlLayoutCounts.COUNT_HEADER_V2);}
			for(final Row row:recipes){
				if(c.checkIndex!=null && !row.id.equals(c.checkIndex)){continue;}
				System.err.println("NATIVE_DDL_RECIPE_START\t"+row.id);
				consensusCpu=consensusWall=0;
				final DynamicDemiLog[][] references=buildReferences(row, members, consensuses);
				writeMetadata(row, metadata);
				if(c.checkIndex==null){writeCounts(row, references, data, counts);}
				else{checkIndex(row, references[0], data);}
				System.err.println("NATIVE_DDL_RECIPE_COMPLETE\t"+row.id);
			}
		}
	}

	private ArrayList<Row> readRecipes(){
		final ArrayList<Row> rows=new ArrayList<Row>();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		int controls=0;
		try(Input input=new Input(c.recipes)){
			RrnaResourceIO.header(input.reader, RECIPE_HEADER, c.recipes);
			for(byte[] text; (text=input.reader.nextLine())!=null;){
				p.set(text);
				require(p.terms()==15, "A native recipe needs all declared word/part/consensus fields");
				final Row row=new Row(p);
				require(ids.add(row.id), "Duplicate native recipe identifier: "+row.id);
				if(row.control){controls++;}
				rows.add(row);
			}
		}
		require(!rows.isEmpty() && rows.size()<=100 && controls==1,
			"Use one explicit legacy standard control and at most100 native DDL recipes");
		if(c.checkIndex!=null){
			Row selected=null;
			for(final Row row:rows){if(row.id.equals(c.checkIndex)){selected=row;}}
			require(selected!=null && selected.terms==1, "checkindex must name one declared single-DDL recipe; no partition/index extension is inferred");
		}
		return rows;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Native Construction     ----------------*/
	/*--------------------------------------------------------------*/

	/** Describes all retained members and accepted saved windows without folds.
	 * A model/window contributes once, irrespective of its number of parent links. */
	private void describeGc(final NcrnaDdlMembers members, final String[] headers, final NcrnaDdlSavedData data){
		final double[][][] stats={NcrnaDdlSavedData.gcStatistics(c.models),
			NcrnaDdlSavedData.gcStatistics(c.models), NcrnaDdlSavedData.gcStatistics(c.models)};
		for(int model=0; model<c.models; model++){
			for(final byte[] sequence:members.members(model)){NcrnaDdlSavedData.addGc(stats[0][model], Tools.calcGC(sequence));}
		}
		data.readNativeCounts(headers, stats, false);
		try(Output out=new Output(c.output+".gc_models.tsv"); Output family=new Output(c.output+".gc_family.tsv")){
			out.writer.println("model_ordinal\tmodel\tsequence_source\tcount\tmissing_gc\tmean_gc\tsd_gc\tmin_gc\tmax_gc");
			family.writer.println("sequence_source\tmodels_with_gc\tmodels_without_gc\tmean_model_gc\tsd_model_gc\tmin_model_gc\tmax_model_gc\tmedian_within_model_sd\tzero_sd_models\tinterpretation");
			for(int source=0; source<stats.length; source++){
				final double[] means=NcrnaDdlSavedData.gcStatistics(1)[0], deviations=new double[c.models];
				int found=0, zero=0;
				for(int model=0; model<c.models; model++){
					final double[] row=stats[source][model];
					require(headers[model].indexOf('\t')<0 && headers[model].indexOf('\n')<0,
						"Full model headers must remain intact in one descriptive TSV field");
					line.clear().append(model).tab().append(headers[model]).tab().append(GC_SOURCES[source])
						.tab().append((long)row[0]).tab().append((long)row[5]);
					appendGcStatistics(line, row);
					out.writer.print(line.nl());
					if(row[0]>0){
						NcrnaDdlSavedData.addGc(means, row[1]);
						final double sd=NcrnaDdlSavedData.gcSd(row);
						deviations[found++]=sd;
						if(sd==0){zero++;}
					}
					if(source==2){
						require(row[0]+row[5]==stats[1][model][0]+stats[1][model][5],
							"Padded and model-specific trimmed descriptions must conserve every accepted model/window pair");
					}
				}
				Arrays.sort(deviations, 0, found);
				line.clear().append(GC_SOURCES[source]).tab().append(found).tab().append(c.models-found);
				appendGcStatistics(line, means);
				line.tab();
				if(found==0){line.append("NA");}
				else{line.append((deviations[(found-1)/2]+deviations[found/2])/2, 9);}
				family.writer.print(line.tab().append(zero).append("\tALL_SAVED_DATA_NO_HELDOUTS_MODEL_MEANS_EQUALLY_WEIGHTED_SAMPLE_SD\n"));
			}
		}
		System.err.println("NCRNA_DDL_ALL_MODEL_GC_DESCRIPTION_COMPLETE\tmodels="+c.models+"\tsaved_windows="+data.windows.size());
	}
	private static void appendGcStatistics(final ByteBuilder out, final double[] stats){
		assert(out!=null && stats!=null && stats.length==6) : "GC output preserves its actual Welford count and missingness";
		if(stats[0]==0){out.append("\tNA\tNA\tNA\tNA");}
		else{
			out.tab().append(stats[1], 9).tab().append(NcrnaDdlSavedData.gcSd(stats), 9)
				.tab().append(stats[3], 9).tab().append(stats[4], 9);
		}
	}

	/** Sums disjoint judging halves; never averages rank sums or their means.
	 * All rows remain in the input. This is reporting only,no biological reload. */
	public static void mergeStudy(final String input, final String prefix){
		require(input!=null && prefix!=null && Tools.testInputFiles(false, true, input)
			&& Tools.testOutputFiles(false, false, false, prefix+".merged.tsv")
			&& Tools.testForDuplicateFiles(true, input, prefix+".merged.tsv"),
			"Study merge requires readable input and a distinct fresh output");
		final ObjectMap<String, long[]> groups=new ObjectMap<String, long[]>(String.class, long[].class);
		final ObjectMap<String, double[]> gc=new ObjectMap<String, double[]>(String.class, double[].class);
		final ArrayList<String> order=new ArrayList<String>();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, STUDY_HEADER, input);
			for(byte[] text; (text=in.reader.nextLine())!=null;){
				p.set(text);
				require(p.terms()==25 && p.parseInt(3)>0 && p.parseInt(4)>0 && p.parseInt(5)>0,
					"A study merge needs complete nonempty measured populations");
				final String id=p.parseString(0)+"\t"+p.parseInt(1)+"\t"+p.parseInt(2);
				long[] values=groups.get(id);
				if(values==null){
					values=new long[13];
					values[0]=p.parseInt(3);
					groups.put(id, values);
					gc.put(id, new double[]{Double.NaN,Double.NaN,Double.NaN,Double.NaN,p.parseDouble(24)});
					order.add(id);
				}
				final int fit=p.parseInt(20), judge=p.parseInt(21), bit=fit<0 ? 4 : 1<<fit;
				require(fit>=-1 && fit<=1 && judge==(fit<0 ? -1 : 1-fit)
					&& (values[12]&bit)==0 && values[0]==p.parseInt(3),
					"Each study setting has one library row or two distinct opposite judging halves,with unchanged model roster");
				values[12]|=bit;
				for(int i=1; i<=11; i++){values[i]=Math.addExact(values[i], p.parseLong(i+3));}
				final double[] composition=gc.get(id);
				require(composition[4]==p.parseDouble(24), "Complementary halves must bind identical family composition");
				final int offset=fit<=0 ? 0 : 2;
				composition[offset]=p.parseInt(1)==0 ? Double.NaN : p.parseDouble(22);
				composition[offset+1]=p.parseInt(2)==0 ? Double.NaN : p.parseDouble(23);
			}
		}
		require(!order.isEmpty(), "Study merge cannot publish a result without measured rows");
		for(final String id:order){require(groups.get(id)[12]==3 || groups.get(id)[12]==4, "Missing or mixed judging half for "+id);}
		try(Output out=new Output(prefix+".merged.tsv")){
			out.writer.println(MERGED_STUDY_HEADER);
			final ByteBuilder b=new ByteBuilder();
			for(final String id:order){
				final long[] v=groups.get(id);
				b.clear().append(id);
				for(int i=0; i<=11; i++){b.tab().append(v[i]);}
				b.tab().append(v[11]/(double)v[1], 9);
				for(final double value:gc.get(id)){b.tab(); appendGc(b, value);}
				out.writer.print(b.append("\tDISJOINT_JUDGING_COUNTS_SUMMED_NOT_CALLER\n"));
			}
		}
	}

	/** Summarizes matched seed observations against their own unfiltered XOR
	 * anchors. The roster binds seed to saved merged table;no counts are rerun. */
	public static void mergeSeeds(final String roster, final String anchors, final String prefix){
		require(roster!=null && anchors!=null && prefix!=null
			&& Tools.testInputFiles(false, true, roster, anchors)
			&& Tools.testOutputFiles(false, false, false, prefix+".seed_summary.tsv"),
			"Cross-seed summary needs a seed/table roster,paired-seed baseline and fresh output");
		final LongList seedList=new LongList();
		final ArrayList<String> paths=new ArrayList<String>();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(roster)){
			RrnaResourceIO.header(in.reader, "ddl_seed\tmerged_table", roster);
			for(byte[] text; (text=in.reader.nextLine())!=null;){
				p.set(text);
				require(p.terms()==2 && p.parseLong(0)>=0, "Seed roster rows bind one explicit nonnegative seed and table");
				final long seed=p.parseLong(0);
				for(int i=0; i<seedList.size; i++){require(seedList.get(i)!=seed, "A seed cannot be counted twice");}
				seedList.add(seed);
				final String path=p.parseString(1);
				require(path.startsWith("/"), "Seed rosters bind absolute native result paths");
				paths.add(path);
			}
		}
		require(paths.size()>=2 && paths.size()<=8, "A matched cross-seed report needs two through eight distinct seeded runs");
		require(Tools.testInputFiles(false, true, paths.toArray(new String[paths.size()])), "Every listed seeded table must be readable");
		final long[] seeds=seedList.toArray();
		final double[][] baseline=new double[seeds.length][];
		final int[][] population=new int[seeds.length][];
		final boolean comparisonAnchors;
		try(Input in=new Input(anchors)){
			final byte[] header=in.reader.nextLine();
			comparisonAnchors=Arrays.equals(header, "ddl_seed\tcomparison_table".getBytes(java.nio.charset.StandardCharsets.UTF_8));
			if(comparisonAnchors){fitAnchors(in, seeds, baseline, population);}
			else{
				require(Arrays.equals(header, (STUDY_HEADER+"\tddl_seed").getBytes(java.nio.charset.StandardCharsets.UTF_8)),
					"Seed baseline is a measured study or an explicit saved-comparison roster");
			}
			for(byte[] text; !comparisonAnchors && (text=in.reader.nextLine())!=null;){
				p.set(text);
				require(p.terms()==26, "Paired seed anchors retain all scientific fields and the explicit seed");
				if(!p.termEquals("OLD_XOR_HIGH_BITS", 0)){continue;}
				for(int s=0; s<seeds.length; s++){
					if(p.parseLong(25)!=seeds[s]){continue;}
					require(baseline[s]==null, "Duplicate unfiltered XOR anchor for a selected seed");
					baseline[s]=studyMetrics(p);
					population[s]=new int[]{p.parseInt(3),p.parseInt(4),p.parseInt(5),p.parseInt(6),p.parseInt(10)};
				}
			}
		}
		final int metrics=comparisonAnchors ? 6 : SUMMARY_METRICS.length;
		final ObjectMap<String, double[]> values=new ObjectMap<String, double[]>(String.class, double[].class);
		final ArrayList<String> order=new ArrayList<String>();
		for(int s=0; s<seeds.length; s++){
			require(baseline[s]!=null, "Every selected seed needs its matching saved unfiltered baseline");
			final String sidecar=new java.io.File(paths.get(s)).getParent()+"/ddl_seed.txt";
			try(Input seedFile=new Input(sidecar)){
				final byte[] text=seedFile.reader.nextLine();
				require(text!=null, "A seeded merged table needs its producer's seed sidecar");
				p.set(text);
				require(p.terms()==1 && p.parseLong(0)==seeds[s] && seedFile.reader.nextLine()==null,
					"Roster seed must match the original producer's seed sidecar");
			}
			try(Input in=new Input(paths.get(s))){
				RrnaResourceIO.header(in.reader, MERGED_STUDY_HEADER, paths.get(s));
				for(byte[] text; (text=in.reader.nextLine())!=null;){
					p.set(text);
					require(p.terms()==22 && p.parseInt(3)==population[s][0] && p.parseInt(4)>0
						&& (population[s][1]<0 || p.parseInt(4)==population[s][1])
						&& p.parseInt(5)==population[s][2] && p.parseInt(6)==population[s][3] && p.parseInt(10)==population[s][4],
						"Seeded settings and their unfiltered anchors must retain identical model/window/parent populations");
					population[s][1]=p.parseInt(4);
					final String setting=p.parseString(0)+"\t"+p.parseInt(1)+"\t"+p.parseInt(2);
					final double[] row=studyMetrics(p);
					for(int metric=0; metric<metrics; metric++){
						final String key=setting+"\t"+SUMMARY_METRICS[metric];
						double[] observations=values.get(key);
						if(observations==null){
							observations=new double[seeds.length];
							Arrays.fill(observations, Double.NaN);
							values.put(key, observations);
							order.add(key);
						}
						require(Double.isNaN(observations[s]), "A setting/metric/seed observation cannot be duplicated");
						observations[s]=row[metric];
					}
				}
			}
			require(Arrays.equals(population[s], population[0]), "Hash seeds cannot change the fixed model/window/parent populations");
		}
		require(!order.isEmpty(), "Cross-seed summary cannot publish an empty experiment");
		for(final String key:order){
			for(final double value:values.get(key)){require(!Double.isNaN(value), "Missing matched seed observation for "+key);}
		}
		try(Output out=new Output(prefix+".seed_summary.tsv")){
			out.writer.println("definition\tshort_blacklist\tlong_blacklist\tmetric\tseeds\tmean\tminimum\tmaximum\tunfiltered_mean\timproved_seeds\ttied_seeds\tworsened_seeds\tpaired_improvement_mean\tpaired_improvement_minimum\tpaired_improvement_maximum");
			final ByteBuilder b=new ByteBuilder();
			for(int row=0; row<order.size(); row++){
				final String key=order.get(row);
				final int metric=row%metrics;
				final boolean lowerBetter=metric<2 || metric==6;
				double sum=0, base=0, min=Double.POSITIVE_INFINITY, max=Double.NEGATIVE_INFINITY;
				double gainSum=0, gainMin=Double.POSITIVE_INFINITY, gainMax=Double.NEGATIVE_INFINITY;
				int improved=0, tied=0;
				final double[] observations=values.get(key);
				for(int s=0; s<seeds.length; s++){
					final double value=observations[s], original=baseline[s][metric];
					sum+=value;
					base+=original;
					min=Math.min(min, value);
					max=Math.max(max, value);
					final double gain=lowerBetter ? original-value : value-original;
					gainSum+=gain;
					gainMin=Math.min(gainMin, gain);
					gainMax=Math.max(gainMax, gain);
					if(value==original){tied++;}
					else if(lowerBetter ? value<original : value>original){improved++;}
				}
				out.writer.print(b.clear().append(key).tab().append(seeds.length).tab().append(sum/seeds.length, 9)
					.tab().append(min, 9).tab().append(max, 9).tab().append(base/seeds.length, 9)
					.tab().append(improved).tab().append(tied).tab().append(seeds.length-improved-tied)
					.tab().append(gainSum/seeds.length, 9).tab().append(gainMin, 9).tab().append(gainMax, 9).nl());
			}
		}
	}

	/** Reuses completed8/10 rows;their seed sidecars and frozen protocol remain binding.
	 * Fit rows have no window/nonzero-pair count. Windows are checked across every
	 * filtered setting/seed;only the six measured ranking metrics are compared. */
	private static void fitAnchors(final Input roster, final long[] seeds, final double[][] baseline, final int[][] population){
		assert(roster!=null && seeds.length==baseline.length && seeds.length==population.length)
			: "Saved comparison anchors pair one exact recipe with every selected DDL seed";
		final LineParser1 p=new LineParser1('\t');
		String cohorts=null;
		for(byte[] text; (text=roster.reader.nextLine())!=null;){
			p.set(text);
			require(p.terms()==2, "Comparison baseline roster contains only DDL seed and absolute completed table");
			final long seed=p.parseLong(0);
			int slot=-1;
			for(int s=0; s<seeds.length; s++){if(seeds[s]==seed){slot=s;}}
			require(slot>=0 && baseline[slot]==null && p.parseString(1).startsWith("/"), "Comparison anchors need each selected seed exactly once");
			final String path=p.parseString(1), sidecar=new java.io.File(path).getParent()+"/ddl_seed.txt";
			try(Input in=new Input(sidecar)){
				final byte[] value=in.reader.nextLine();
				require(value!=null, "Saved comparison anchor needs its actual producer's DDL seed sidecar");
				p.set(value);
				require(p.terms()==1 && p.parseLong(0)==seed && in.reader.nextLine()==null, "Comparison roster and original DDL seed must agree");
			}
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, SummarizeNcrnaLayoutFits.HEADER, path);
				for(byte[] value; (value=in.reader.nextLine())!=null;){
					p.set(value);
					if(!p.termEquals(BEST_PAIRED, 0)){continue;}
					require(baseline[slot]==null && p.terms()==44 && p.termEquals("RF00174", 1) && p.parseInt(2)==8 && p.parseInt(3)==10
						&& p.termEquals("MIDDLE_RIGHT", 4) && p.termEquals("MIXED", 5) && p.parseInt(6)==512
						&& p.parseInt(7)==0 && p.parseInt(8)==0 && p.termEquals("POOL_BOTH", 9) && p.parseInt(10)==0
						&& p.parseInt(11)==100 && p.termEquals("member_total", 12) && p.parseInt(13)==2
						&& p.termEquals("NATIVE_DDL", 14) && p.termEquals("WRAP", 15)
						&& p.parseInt(34)==5000 && p.parseLong(35)==20261113 && p.termEquals("QUERY_DDL", 37),
						"Filter-only anchors must be the completed8/masked10mixed512pooled100 recipe and unchanged weight-search protocol");
					if(cohorts==null){cohorts=p.parseString(36);}
					require(!cohorts.isEmpty() && p.termEquals(cohorts, 36), "All seed anchors must share their frozen cohort source");
					population[slot]=new int[]{p.parseInt(31),-1,p.parseInt(32),p.parseInt(32)-p.parseInt(42),p.parseInt(33)};
					final int[] fields={43,26,38,39,40,41};
					baseline[slot]=new double[fields.length];
					for(int i=0; i<fields.length; i++){
						final double metric=p.parseDouble(fields[i]);
						require(Double.isFinite(metric) && metric>=0, "Completed comparison anchor metrics must be finite and nonnegative");
						baseline[slot][i]=metric;
					}
				}
			}
			require(baseline[slot]!=null, "Every saved comparison table must contain the exact unfiltered8/10recipe");
		}
	}
	private static double[] studyMetrics(final LineParser1 p){
		final double[] result=new double[SUMMARY_FIELDS.length];
		for(int i=0; i<result.length; i++){
			result[i]=p.parseDouble(SUMMARY_FIELDS[i]);
			require(Double.isFinite(result[i]) && result[i]>=0, "Scientific seed metrics must be finite and nonnegative");
		}
		return result;
	}

	/** Each actual part is an ordinary native DDL object from the existing factory. */
	private DynamicDemiLog[] createParts(final Row row, final int count){
		final DynamicDemiLog[] parts=new DynamicDemiLog[count];
		for(int part=0; part<parts.length; part++){
			assert(powerOfTwo(row.width(part))) : "Native DDL insertion uses width-1 as its low-bit mask; each actual member/consensus part must be a power of two: "+row.id+", part="+part+", width="+row.width(part);
			parts[part]=DynamicDemiLog.create(row.width(part), row.wordLength(part), sketchSeed, 0, false, false);
		}
		return parts;
	}

	/** The only custom reference decision is the voted key retained in each bin. */
	private DynamicDemiLog[][] buildReferences(final Row row, final NcrnaDdlMembers members, final byte[][] consensuses){
		final DynamicDemiLog[][] references=new DynamicDemiLog[row.terms][c.models];
		final DynamicDemiLog[] voter=createParts(row, row.memberTerms);
		for(int model=0; model<c.models; model++){
			final byte[][] sequences=members.members(model);
			require(sequences.length>0 && sequences.length<=Integer.MAX_VALUE-row.bonus,
				"Each real member casts one native winner vote per bin; the declared consensus bonus must not overflow");
			final IntHashMap[] votes=new IntHashMap[row.memberTerms];
			for(int part=0; part<votes.length; part++){votes[part]=new IntHashMap(row.width(part));}
			int voters=0;
			for(int member=0; member<sequences.length; member++){
				if(excludedMembers!=null && excludedMembers[model][member]){continue;}
				voters++;
				final byte[] sequence=sequences[member];
				for(int part=0; part<row.memberTerms; part++){
					fill(voter[part], sequence, row, part, true);
					final char[] keys=voter[part].maxArray();
					for(int bin=0; bin<keys.length; bin++){
						if(keys[bin]!=0){votes[part].increment((bin<<16)|keys[bin]);}
					}
				}
			}
			final int bonus=replicateBonus<0 ? row.bonus : consensusBonus(replicateBonus, replicateScale, voters);
			assert(bonus>=0 && bonus<=Integer.MAX_VALUE-voters) : "Consensus plus actual member votes must fit the native integer vote map";
			if(row.bonus>0){
				//TODO: Probable bug - this historical study bonuses only consensus
				//DDL maxima. Brian's explicit consensus-weight order permits several
				//consensus bin/key candidates, then member support/greatest-key ties
				//(DECISIONS NCRNA_DDL_EXPLICIT_CONSENSUS_WEIGHTS_20261113). The new
				//builder uses all distinct eligible consensus pairs. Preserve these
				//old measured artifacts with this limitation; no historical rerun here.
				final long cpu=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
				for(int part=0; part<row.memberTerms; part++){
					fill(voter[part], consensuses[model], row, part, row.consensusMode.equals("POOL_BOTH"));
					final char[] keys=voter[part].maxArray();
					for(int bin=0; bin<keys.length; bin++){
						if(keys[bin]!=0 && bonus>0){votes[part].increment((bin<<16)|keys[bin], bonus);}
					}
				}
				consensusCpu+=clock.getCurrentThreadCpuTime()-cpu;
				consensusWall+=System.nanoTime()-wall;
			}
			for(int part=0; part<row.memberTerms; part++){
				references[part][model]=chooseReference(votes[part], row.width(part), model, row.wordLength(part), !row.control);
			}
			if(retention!=null){
				final char[] consensusKeys=voter[0].maxArray(), keptKeys=references[0][model].maxArray();
				int matches=0;
				for(int bin=0; bin<keptKeys.length; bin++){if(consensusKeys[bin]!=0 && consensusKeys[bin]==keptKeys[bin]){matches++;}}
				retention[sizeThirds[model]]+=matches/(double)keptKeys.length;
			}
			if(row.consensusBins>0){
				final int part=row.memberTerms;
				assert(powerOfTwo(row.consensusBins)) : "The separate consensus must retain native low-bit binning: "+row.id+", width="+row.consensusBins;
				final long cpu=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
				final DynamicDemiLog consensus=DynamicDemiLog.create(row.consensusBins, row.longK, sketchSeed, 0, false, false);
				fill(consensus, consensuses[model], row, part, row.consensusMode.equals("BOTH"));
				references[part][model]=consensus;
				consensusCpu+=clock.getCurrentThreadCpuTime()-cpu;
				consensusWall+=System.nanoTime()-wall;
			}
		}
		return references;
	}

	/** Existing primitive vote map, followed by an existing DDL-object wrapper.
	 * Neither support counts nor the construction arrays survive in the reference. */
	public static DynamicDemiLog chooseReference(final IntHashMap votes, final int bins,
			final int model, final int k, final boolean greatestTie){
		final char[] selected=new char[bins];
		final int[] support=new int[bins], keys=votes.keys(), values=votes.values();
		final int invalid=votes.invalid();
		for(int i=0; i<keys.length; i++){
			if(keys[i]==invalid){continue;}
			final int bin=keys[i]>>>16, key=keys[i]&65535, count=values[i];
			assert(bin<bins && key>0 && count>0) : "Reference selection sees only nonempty native winner votes in its actual part";
			if(count>support[bin] || count==support[bin] && (greatestTie ? key>selected[bin] : key<selected[bin])){
				support[bin]=count;
				selected[bin]=(char)key;
			}
		}
		return DynamicDemiLog.fromArray(selected, model, null, k, -1);
	}

	/** Fixed zero is meaningful;only scaled recipes round upward with minimum one. */
	public static int consensusBonus(final int fixed, final double scale, final int voters){
		require(fixed>=0 && Double.isFinite(scale) && scale>=0 && scale<=1 && voters>=0,
			"Consensus weights require a nonnegative fixed bonus or declared fraction of actual remaining voters");
		return scale>0 ? Math.max(1, (int)Math.ceil(scale*voters)) : fixed;
	}

	/** Feeds the existing strict word iterator into native hash-and-store.
	 * Masking uses the existing BBTools middle-mask helper after span validation. */
	private void fill(final DynamicDemiLog ddl, final byte[] bases, final Row row, final int part, final boolean includeShort){
		ddl.clear();
		active=ddl;
		activePartitionBits=row.partitionBits;
		activeShortK=row.shortK;
		try{
			if(part<row.memberTerms && row.shortBins>0){
				if(part==0){
					if(includeShort){words(bases, row.shortK, false, row.control || oldTags ? 0 : 1L<<(2*row.shortK));}
				}else{words(bases, row.longK, row.maskedLong, row.control || oldTags ? 0 : 1L<<(2*row.longK));}
			}else{
				if(includeShort){words(bases, row.shortK, false, wordTag(row.shortK, true, row.control || oldTags));}
				words(bases, row.longK, row.maskedLong, wordTag(row.longK, false, row.control || oldTags));
			}
		}finally{active=null; activePartitionBits=0; activeShortK=0;}
	}
	private void words(final byte[] bases, final int k, final boolean masked, final long tag){
		activeMask=AbstractKmerTable.makeMiddleMask(k, false, masked ? 1 : 0);
		activeTag=tag;
		CardinalityTracker.forEachKmer(bases, k, k, false, true, this);
	}
	@Override
	public void accept(final long word, final int k){
		final int raw=(int)(word&activeMask);
		if(wordFrequencies!=null){
			assert(raw>=0 && raw<wordFrequencies.length) : "Study frequencies use the declared forward or middle-masked raw-word domain";
			require(wordFrequencies[raw]<Integer.MAX_VALUE, "Raw-word occurrence count must fit the primitive frequency table");
			wordFrequencies[raw]++;
			return;
		}
		assert(active!=null) : "Native insertion requires the currently constructed DDL object";
		final long tagged=(word&activeMask)|activeTag;
		if(activePartitionBits==0){storeWord(active, tagged, blacklist);}
		else if(blacklist==null || !blacklist.contains(tagged)){
			assert((active.buckets==512 || active.buckets==1024) && DynamicDemiLog.exponentBits()==4)
				: "Routing excludes actual bin selectors and4+12mantissa windows;the wrapped exponent is not unconditionally independent";
			final long hash=Tools.hash64shift(tagged^sketchSeed);
			final int partition=partitionFor(hash, active.buckets, activePartitionBits, k==activeShortK);
			active.storeHashed(tagged, hash, partition, activePartitionBits);
		}
	}

	/** Routing excludes bucket selectors and all possible4+12mantissa windows.
	 * The wrapped exponent can depend on this region when hash>>>36==0.
	 * Public for the routing-bit fixture; the original native hash is unchanged. */
	public static int partitionFor(final long hash, final int buckets, final int partitionBits, final boolean shortWord){
		assert((buckets==512 || buckets==1024) && (partitionBits==1 || partitionBits==2))
			: "The released single-DDL partition experiment reserves one or two selector bits at512/1024bins";
		if(shortWord){return 0;}
		if(partitionBits==1){return 1;}
		final int binBits=Integer.numberOfTrailingZeros(buckets), routingBits=36-binBits;
		final long routing=(hash>>>binBits)&((1L<<routingBits)-1);
		return 1+(int)((routing*3)>>>routingBits);
	}

	/** Shared namespace for sketch insertion,blacklist membership and reporting. */
	public static long wordTag(final int k, final boolean shortWord, final boolean xor){
		assert(k>0 && k<=31) : "Packed nucleotide words need an explicit span before a leading length bit is added";
		return xor ? (shortWord ? SHORT_TAG : LONG_TAG) : 1L<<(2*k);
	}
	/** Rejects literal tagged words before hashing;false means no sketch change. */
	public static boolean storeWord(final DynamicDemiLog ddl, final long tagged, final LongHashSet blacklist){
		assert(ddl!=null) : "The primitive sink inserts into its caller-owned native sketch";
		if(blacklist!=null && blacklist.contains(tagged)){return false;}
		ddl.hashAndStore(tagged);
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native Scoring        ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean replicated(final Config c){return "replicates".equals(c.study) || indexStudy(c) || histogramStudy(c) || gcLosses(c);}
	private static boolean gcLosses(final Config c){return "gclosses".equals(c.study);}
	private static boolean gcAllLosses(final Config c){return "gcalllosses".equals(c.study);}
	private static boolean gcAllPower(final Config c){return "gcallpower".equals(c.study);}
	private static boolean gcAllData(final Config c){return gcAllLosses(c) || gcAllPower(c) || "replicategcall".equals(c.study);}
	private static boolean indexStudy(final Config c){return "replicateindex".equals(c.study);}
	private static boolean gcDescription(final Config c){return "gcdescribe".equals(c.study);}
	private static boolean gcStudy(final Config c){return "replicategc".equals(c.study) || "replicategcall".equals(c.study);}
	private static boolean followup(final Config c){return "replicatefollowup".equals(c.study);}
	private static boolean sizeStudy(final Config c){return "replicatesizes".equals(c.study);}
	private static boolean extendedSizes(final Config c){return "replicatesizesextended".equals(c.study);}
	private static boolean definitionStudy(final Config c){return "replicatedefinitions".equals(c.study);}
	private static boolean histogramStudy(final Config c){return followup(c) || sizeStudy(c) || extendedSizes(c) || definitionStudy(c) || gcStudy(c);}

	/** One hash seed and one whole five-fold assignment per independent worker.
	 * Rows retain disjoint folds; the report pools them before estimating variance. */
	private void replicates(final ArrayList<Row> recipes, final NcrnaDdlMembers members,
			final String[] headers, final byte[][] consensuses, final NcrnaDdlSavedData data,
			final NcrnaDdlFolds folds){
		final byte[][] names=new byte[headers.length][];
		for(int i=0; i<names.length; i++){names[i]=headers[i].getBytes(java.nio.charset.StandardCharsets.UTF_8);}
		final NcrnaDdlCohorts cohorts=NcrnaDdlCohorts.read(c.identities, data, names, consensuses);
		if(!c.fixedCheck){folds.prepare(members, cohorts);}
		if(indexStudy(c)){
			indexPercentiles(data, cohorts, folds);
			return;
		}
		data.readNativeCounts(headers, null, gcStudy(c));
		if(gcAllData(c)){readAllGcStatistics(members);}
		final Row[] selected=new Row[followup(c) ? 4 : 3];
		for(final Row row:recipes){
			if(row.control || row.terms!=1 || row.memberBins!=512 || row.shortBins!=0 || !row.maskedLong
				|| !row.consensusMode.equals("POOL_BOTH") || row.bonus!=100){continue;}
			final int pair=row.shortK==7 && row.longK==9 ? 0 : row.shortK==7 && row.longK==10 ? 1
				: row.shortK==8 && row.longK==10 ? 2 : -1;
			if(pair>=0){
				require(selected[pair]==null, "The replicated roster needs each released word pair exactly once");
				selected[pair]=row;
			}
		}
		for(int i=0; i<3; i++){require(selected[i]!=null, "The replicated roster must contain all three pooled100,mixed512 recipes");}
		if(followup(c)){selected[3]=new Row(selected[2], 512, 0, 11);}
		try(Output output=new Output(c.output+".replicates.tsv"); Output assignment=new Output(c.output+".folds.tsv");
				Output metadata=new Output(c.output+".cohorts.tsv"); Output models=new Output(c.output+".models.tsv");
				Output gcFit=gcStudy(c) && !gcAllData(c) ? new Output(c.output+".gc_fit.tsv") : null){
			output.writer.println(gcStudy(c) ? GC_REPLICATE_HEADER : histogramStudy(c) ? FOLLOWUP_REPLICATE_HEADER : REPLICATE_HEADER);
			if(gcFit!=null){gcFit.writer.println(GC_FIT_HEADER);}
			assignment.writer.println("input_index\tassembly_group\tjudge_fold\tassignment_seed\tcomparable_group");
			metadata.writer.println("cohort\tlow_identity_cut\tcut_float_bits\tfinite_identities\tinterpretation");
			final float cut=c.fixedCheck ? cohorts.cuts[0] : folds.lowCut;
			metadata.writer.print(line.clear().append(c.fixedCheck ? "FIXED_NOT_SEEN" : "ALL_COMPARABLE").tab()
				.append(cut, 9).tab().append(Float.floatToRawIntBits(cut)).tab().append(c.fixedCheck ? cohorts.finiteNotSeen : folds.finiteIdentities)
				.append("\tNEAREST_RANK_QUARTILE_FROM_FROZEN_IDENTITIES_BEFORE_NEW_SCORES\n"));
			models.writer.println("fold\tmodel_ordinal\tvoting_members\texcluded_members\tsize_third");
			final int[] original=c.fixedCheck ? folds.fixed(cohorts, assignment) : folds.assign(c.foldSeed, c.folds, assignment);
			for(int judge=0; judge<(c.fixedCheck ? 2 : c.folds); judge++){
				inputSides=new int[original.length];
				for(int i=0; i<original.length; i++){inputSides[i]=original[i]==judge ? 1 : 0;}
				fitSide=0;
				excludedMembers=c.fixedCheck ? null : folds.exclusions(original, judge, members);
				final int[] voters=votingCounts(members);
				final long[] ordered=new long[c.models];
				for(int model=0; model<c.models; model++){ordered[model]=((long)voters[model]<<32)|model;}
				Arrays.sort(ordered);
				sizeThirds=new int[c.models];
				for(int i=0; i<ordered.length; i++){sizeThirds[(int)ordered[i]]=Math.min(2, (int)(3L*i/c.models));}
				for(int model=0; model<c.models; model++){
					models.writer.print(line.clear().append(judge).tab().append(model).tab().append(voters[model])
						.tab().append(members.members(model).length-voters[model]).tab().append(sizeThirds[model]).nl());
				}
				int empty=0;
				for(final int count:voters){if(count==0){empty++;}}
				if(gcStudy(c)){
					gcReplicates(selected[2], judge, members, consensuses, data, cohorts, folds, empty, output, gcFit);
					continue;
				}
				if(sizeStudy(c) || extendedSizes(c) || definitionStudy(c)){
					blacklistReplicates(selected[2], judge, members, consensuses, data, cohorts, folds, empty, output);
					continue;
				}
				if(c.indexReference && !followup(c)){
					final int[] ranks=new int[data.parents.size()];
					for(final Parent parent:data.parents){ranks[parent.index]=parent.nativeRank;}
					emitReplicate("CONSENSUS_INDEX", selected[0], 0, 0, 100, 0, judge, ranks,
						data, cohorts, folds, empty, null, 0, 0, output);
				}
				for(int pair=0; pair<selected.length; pair++){
					if(followup(c) && pair<2){continue;}
					final Row base=selected[pair];
					final int shortSize=c.shortPrefixes[Math.min(pair, 2)], longSize=c.longPrefix;
					final long[] shortWords=shortSize>0 ? frequencies(base.shortK, false, 2, members, data) : new long[0];
					final long[] longWords=longSize>0 ? frequencies(base.longK, true, 2, members, data) : new long[0];
					for(final int bins:new int[]{512,1024}){
						for(int partition=0; partition<=2; partition++){
							final Row row=new Row(base, bins, partition);
							for(int filtered=0; filtered<2; filtered++){
								if(followup(c) && pair==2 && (partition!=2 || filtered==0)){continue;}
								blacklist=null;
								if(filtered!=0){
									blacklist=new LongHashSet(Math.max(2, shortSize+longSize));
									addBlacklist(shortWords, shortSize, row.shortK, true);
									addBlacklist(longWords, longSize, row.longK, false);
								}
								final String method="s"+row.shortK+"_m"+row.longK+"_"+PARTITION_NAMES[partition]+"_"+bins+"_"+(filtered==0 ? "plain" : "padding");
								final int[] bonuses=followup(c) && pair==2 ? FOLLOWUP_BONUSES : new int[]{100};
								for(final int bonus:bonuses){
									runReplicate(bonus==100 ? method : method+"_bonus"+bonus, row,
										filtered==0 ? 0 : shortSize, filtered==0 ? 0 : longSize,
										bonus, 0, judge, members, consensuses, data, cohorts, folds, empty, output);
								}
							}
						}
					}
				}
				blacklist=null;
				if(followup(c)){continue;}
				for(final int bonus:new int[]{0,5,1000000}){
					runReplicate("bonus_fixed_"+bonus, selected[0], 0, 0, bonus, 0, judge,
						members, consensuses, data, cohorts, folds, empty, output);
				}
				for(int tenth=1; tenth<=10; tenth++){
					runReplicate("bonus_tenths_"+tenth, selected[0], 0, 0, 0, tenth/10.0, judge,
						members, consensuses, data, cohorts, folds, empty, output);
				}
			}
		}finally{
			excludedMembers=null;
			retention=null;
			replicateBonus=-1;
			replicateScale=0;
			blacklist=null;
			fitSide=-1;
		}
	}

	/** Saved index ranks and each population's global low-identity label are fixed.
	 * Pooling disjoint folds restores this same histogram for EVERY assignment;
	 * no reference construction,query scoring or repeated index measurement occurs. */
	private void indexPercentiles(final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts, final NcrnaDdlFolds folds){
		final long[][] histogram=new long[2][c.models];
		final long[] counts=new long[2], sums=new long[2];
		for(final Parent parent:data.parents){
			if(!replicateParent(parent)){continue;}
			final int rank=parent.nativeRank;
			require(rank>0 && rank<=c.models, "The reference percentile uses the original positive saved rank of each comparable locus");
			final boolean low=c.fixedCheck ? cohorts.bins[parent.index]==0 : folds.low[parent.index];
			histogram[0][rank-1]++;
			counts[0]++;
			sums[0]+=rank;
			if(low){
				histogram[1][rank-1]++;
				counts[1]++;
				sums[1]+=rank;
			}
		}
		try(Output out=new Output(c.output+".index.tsv")){
			out.writer.println("cohort\tparents\tlow_parents\trank_sum\tlow_rank_sum\trank_p99\tlow_rank_p99\tassignment_sd\tinterpretation");
			out.writer.print(line.clear().append(c.fixedCheck ? "FIXED_NOT_SEEN" : "REBUILT_FIVE_FOLD")
				.tab().append(counts[0]).tab().append(counts[1]).tab().append(sums[0]).tab().append(sums[1])
				.tab().append(SummarizeNcrnaDdlReplicates.percentile99(histogram[0]))
				.tab().append(SummarizeNcrnaDdlReplicates.percentile99(histogram[1]))
				.append("\t0\tEXACTLY_ASSIGNMENT_INVARIANT_SAVED_RANKS_NOT_RESCORED\n"));
		}
	}

	/** Size jobs omit the winner already measured by the nearby-bonus worker.
	 * Alternative definitions run later at the explicitly selected prefix sizes. */
	private void blacklistReplicates(final Row base, final int judge, final NcrnaDdlMembers members,
			final byte[][] consensuses, final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts,
			final NcrnaDdlFolds folds, final int empty, final Output output){
		assert(base.shortK==8 && base.longK==10) : "The size/definition study holds the released8/masked10 word pair fixed";
		final Row row=new Row(base, 1024, 2);
		for(final int definition:definitionStudy(c) ? new int[]{1,3,0} : new int[]{2}){
			final int shortMax=definitionStudy(c) ? c.shortPrefixes[2] : extendedSizes(c) ? 16384 : 1024;
			final int longMax=definitionStudy(c) ? c.longPrefix : extendedSizes(c) ? 1024 : 32768;
			final long[] shortWords=shortMax>0 ? frequencies(8, false, definition, members, data) : new long[0];
			final long[] longWords=longMax>0 ? frequencies(10, true, definition, members, data) : new long[0];
			final int[][] sizes=definitionStudy(c) ? new int[][]{{shortMax,longMax}} : extendedSizes(c) ? SIZE_EXTENSION : SIZE_GRID;
			for(final int[] size:sizes){
				final int shortSize=size[0], longSize=size[1];
				if(sizeStudy(c) && shortSize==256 && longSize==4096){continue;}
				blacklist=new LongHashSet(Math.max(2, shortSize+longSize));
				addBlacklist(shortWords, shortSize, 8, true);
				addBlacklist(longWords, longSize, 10, false);
				final String name=definitionStudy(c) ? "blacklist_"+DEFINITION_NAMES[definition]+"_s"+shortSize+"_l"+longSize
					: shortSize==0 && longSize==0 ? "s8_m10_quarter_1024_plain" : sizeMethod(shortSize, longSize);
				runReplicate(name, row, shortSize, longSize, 100, 0, judge,
					members, consensuses, data, cohorts, folds, empty, output);
			}
		}
		blacklist=null;
	}
	public static String sizeMethod(final int shortSize, final int longSize){
		assert(shortSize>=0 && longSize>=0) : "A blacklist size method names nonnegative prefix lengths";
		return "s8_m10_quarter_1024_padding_s"+shortSize+"_l"+longSize;
	}
	private static int[][] gridSizes(final boolean extended){
		final int base=GRID_SHORT.length*GRID_LONG.length;
		final int[][] result=new int[base+(extended ? SIZE_EXTENSION.length : 0)][];
		int next=0;
		for(final int s:GRID_SHORT){
			for(final int l:GRID_LONG){result[next++]=new int[]{s,l};}
		}
		if(extended){for(final int[] size:SIZE_EXTENSION){result[next++]=size;}}
		assert(next==result.length) : "The combined grid preserves all72original settings and appends exactly22new settings";
		return result;
	}

	private int[] votingCounts(final NcrnaDdlMembers members){
		final int[] result=new int[c.models];
		for(int model=0; model<c.models; model++){
			for(int i=0; i<members.members(model).length; i++){
				if(excludedMembers==null || !excludedMembers[model][i]){result[model]++;}
			}
		}
		return result;
	}

	/** Fit both GC sources from the opposite fold in the original accumulation order. */
	private double[][][] fitGc(final int judge, final NcrnaDdlMembers members, final NcrnaDdlSavedData data,
			final Output fitOut, final long[][] fallbackCounts, final long[][] sourceCounts){
		assert(gcAllData(c) || fitSide==0 && inputSides!=null && (excludedMembers==null)==c.fixedCheck)
			: "GC fitting uses the declared fit side and the matching fixed-check or rebuilt-member exclusions";
		final double[][] stats=NcrnaDdlSavedData.gcStatistics(c.models);
		for(final Window window:data.windows){
			if(!gcAllData(c) && inputSides[window.input.index]!=fitSide || window.acceptedModelBits==null){continue;}
			for(int word=0; word<window.acceptedModelBits.length; word++){
				for(long bits=window.acceptedModelBits[word]; bits!=0; bits&=bits-1){
					final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
					NcrnaDdlSavedData.addGc(stats[model], window.gc);
				}
			}
		}
		final double[][] memberStats=NcrnaDdlSavedData.gcStatistics(c.models);
		final double[] family=NcrnaDdlSavedData.gcStatistics(1)[0], deviations=new double[c.models];
		double minimum=Double.POSITIVE_INFINITY, maximum=Double.NEGATIVE_INFINITY, offset=0;
		int observed=0, paired=0;
		for(int model=0; model<c.models; model++){
			for(int member=0; member<members.members(model).length; member++){
				if(!gcAllData(c) && excludedMembers!=null && excludedMembers[model][member]){continue;}
				NcrnaDdlSavedData.addGc(memberStats[model], Tools.calcGC(members.members(model)[member]));
			}
			if(memberStats[model][0]>0){
				NcrnaDdlSavedData.addGc(family, memberStats[model][1]);
				deviations[observed++]=NcrnaDdlSavedData.gcSd(memberStats[model]);
				minimum=Math.min(minimum, memberStats[model][3]);
				maximum=Math.max(maximum, memberStats[model][4]);
				if(stats[model][0]>0){offset+=stats[model][1]-memberStats[model][1]; paired++;}
			}
		}
		offset=paired>0 ? offset/paired : 0;
		Arrays.sort(deviations, 0, observed);
		final double familySd=observed==0 ? Double.NaN : (deviations[(observed-1)/2]+deviations[observed/2])/2;
		final double familyMean=observed==0 ? Double.NaN : family[1];
		final double[][][] effective=new double[2][c.models][4];
		for(int model=0; model<c.models; model++){
			final double[] member=effective[1][model], raw=memberStats[model];
			sourceCounts[0][model]=(long)stats[model][0];
			sourceCounts[1][model]=(long)raw[0];
			member[0]=raw[0]>0 ? raw[1] : familyMean;
			member[1]=raw[0]>0 ? NcrnaDdlSavedData.gcSd(raw) : familySd;
			member[2]=raw[0]>0 ? raw[3] : observed==0 ? Double.NaN : minimum;
			member[3]=raw[0]>0 ? raw[4] : observed==0 ? Double.NaN : maximum;
			final double[] window=effective[0][model];
			window[0]=stats[model][0]>0 ? stats[model][1] : member[0]+offset;
			window[1]=stats[model][0]>0 ? NcrnaDdlSavedData.gcSd(stats[model]) : member[1];
			window[2]=stats[model][0]>0 ? stats[model][3] : member[2]+offset;
			window[3]=stats[model][0]>0 ? stats[model][4] : member[3]+offset;
			for(int arm=0; arm<2; arm++){
				if(arm==0 && stats[model][0]==0){fallbackCounts[arm][0]++;}
				if(raw[0]==0 && (arm==1 || stats[model][0]==0)){fallbackCounts[arm][1]++;}
				if(effective[arm][model][1]==0){fallbackCounts[arm][2]++;}
				if(Double.isNaN(effective[arm][model][0])){fallbackCounts[arm][3]++;}
				if(fitOut==null){continue;}
				line.clear().append(judge).tab().append(model).tab().append(GC_ARMS[arm])
					.tab().append((long)stats[model][0]).tab().append((long)raw[0])
					.tab().append((long)stats[model][5]).tab().append((long)raw[5]);
				for(final double value:effective[arm][model]){
					line.tab();
					if(Double.isNaN(value)){line.append("NA");}
					else if(gcAllData(c)){line.append(Double.toString(value));}
					else{line.append(value, 9);}
				}
				fitOut.writer.print(line.tab().append(arm==0 && stats[model][0]==0)
					.tab().append(raw[0]==0 && (arm==1 || stats[model][0]==0)).tab().append(arm==0 ? offset : 0, 9)
					.tab().append(paired).nl());
			}
		}
		return effective;
	}

	/** Load the same full-precision all-data table for full and fixed ranking.
	 * No query-population GC or per-fold statistic is computed in this path. */
	private void readAllGcStatistics(final NcrnaDdlMembers members){
		require(gcAllData(c) && c.gcAllFit!=null && c.gcAllLimits!=null,
			"Both full and fixed ranking populations must load the same complete all-data GC statistics and limits");
		allGcStatistics=new double[2][c.models][4];
		allGcFallback=new long[2][4];
		allGcLosslessZ=new double[2];
		allGcWindowCounts=new long[c.models];
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try(Input in=new Input(c.gcAllFit)){
			RrnaResourceIO.header(in.reader, GC_FIT_HEADER, c.gcAllFit);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int model=rows/2, arm=rows&1;
				require(p.terms()==15 && model<c.models && p.parseInt(0)==-1 && p.parseInt(1)==model
					&& p.termEquals(GC_ARMS[arm], 2) && p.parseLong(3)>=0 && p.parseLong(4)>=0
					&& p.parseLong(4)+p.parseLong(6)==members.members(model).length,
					"Shared GC statistics contain every original model/source once,all retained members,and no judging-fold fit");
				final double[] stats=allGcStatistics[arm][model];
				if(arm==0){allGcWindowCounts[model]=p.parseLong(3);}
				for(int metric=0; metric<4; metric++){
					stats[metric]=p.termEquals("NA", metric+7) ? Double.NaN : Double.parseDouble(p.parseString(metric+7));
				}
				require(Double.isNaN(stats[0]) || Double.isFinite(stats[0]) && stats[1]>=0 && Double.isFinite(stats[1])
					&& stats[2]<=stats[3], "Shared GC statistics preserve exact finite means,SDs and ranges,including shifted all-data fallback ranges");
				if(p.termEquals("true", 11)){allGcFallback[arm][0]++;}
				if(p.termEquals("true", 12)){allGcFallback[arm][1]++;}
				if(stats[1]==0){allGcFallback[arm][2]++;}
				if(Double.isNaN(stats[0])){allGcFallback[arm][3]++;}
				rows++;
			}
		}
		require(rows==2*c.models, "A truncated GC table cannot silently cause a per-population or per-fold refit");
		try(Input in=new Input(c.gcAllLimits)){
			RrnaResourceIO.header(in.reader, GC_ALL_LIMIT_HEADER, c.gcAllLimits);
			for(int arm=0; arm<2; arm++){
				final byte[] row=in.reader.nextLine();
				require(row!=null, "Each GC statistics source needs its measured all-data lossless limit");
				p.set(row);
				require(p.terms()==10 && p.termEquals(GC_ARMS[arm], 0) && p.parseLong(1)>0,
					"All-data lossless limits preserve both source identities and complete accepted-pair populations");
				allGcLosslessZ[arm]=p.termEquals("INF", 7) ? Double.POSITIVE_INFINITY
					: p.termEquals("NA", 7) ? Double.NaN : Double.parseDouble(p.parseString(7));
				require(gcAllPower(c) || c.gcAdditions || Double.isFinite(allGcLosslessZ[arm]) && allGcLosslessZ[arm]>=0,
					"This ranking roster requires a measured finite all-data lossless SD width for each source;unavailable limits must not be silently substituted");
			}
			require(in.reader.nextLine()==null, "Shared GC limits cannot append another population or fit fold");
		}
	}

	/** All-data GC is fitted once,independently of ranking folds and DDL hashes. */
	private void allDataGcLosses(final NcrnaDdlMembers members, final String[] headers, final NcrnaDdlSavedData data){
		assert(gcAllData(c) && excludedMembers==null && inputSides==null)
			: "Brian's all-data GC instruction excludes every ranking holdout from the statistics path";
		data.readNativeCounts(headers, null, true);
		final long[][] counts=new long[2][c.models], fallback=new long[2][4];
		try(Output fits=new Output(c.output+".gc_all_fit.tsv"); Output limits=new Output(c.output+".gc_all_limits.tsv");
				Output losses=new Output(c.output+".gc_all_losses.tsv")){
			fits.writer.println(GC_FIT_HEADER);
			final double[][][] stats=fitGc(-1, members, data, fits, fallback, counts);
			limits.writer.println(GC_ALL_LIMIT_HEADER);
			losses.writer.println("gc_source\tfilter\tparameter\tavailable\taccepted_pairs\tpairs_removed\twindows_with_accepted_models\twindows_lost\tcomparable_loci\tloci_lost\tinterpretation");
			for(int arm=0; arm<2; arm++){
				long pairs=0, finite=0, missing=0, infinite=0;
				double maximum=0, losslessZ=0, margin=0;
				for(final Window window:data.windows){
					if(window.acceptedModelBits==null){continue;}
					for(int word=0; word<window.acceptedModelBits.length; word++){
						for(long bits=window.acceptedModelBits[word]; bits!=0; bits&=bits-1){
							final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
							final double[] modelStats=stats[arm][model];
							pairs++;
							if(Double.isNaN(window.gc) || Double.isNaN(modelStats[0])){missing++; continue;}
							final double distance=gcDistance(window.gc, modelStats[0], modelStats[1]);
							if(Double.isInfinite(distance)){infinite++; losslessZ=Double.POSITIVE_INFINITY;}
							else{
								finite++;
								maximum=Math.max(maximum, distance);
								double z=distance;
								while(!gcWithin(window.gc, modelStats, z, false)){z=Math.nextUp(z);}
								losslessZ=Math.max(losslessZ, z);
							}
							double m=Math.max(0, Math.max(modelStats[2]-window.gc, window.gc-modelStats[3]));
							while(!gcWithin(window.gc, modelStats, m, true)){m=Math.nextUp(m);}
							margin=Math.max(margin, m);
						}
					}
				}
				require(pairs==finite+missing+infinite, "Every accepted all-data model/window pair has one finite,missing or infinite GC distance");
				line.clear().append(GC_ARMS[arm]).tab().append(pairs).tab().append(finite).tab().append(missing).tab().append(infinite).tab();
				line.append(finite==0 ? "NA" : Double.toString(maximum)).tab()
					.append(infinite>0 ? "INF" : finite==0 ? "NA" : Double.toString(maximum)).tab()
					.append(infinite>0 ? "INF" : finite==0 ? "NA" : Double.toString(losslessZ)).tab()
					.append(finite+infinite==0 ? "NA" : Double.toString(margin));
				limits.writer.print(line.append("\tALL_DATA_NO_HELDOUTS_EXACT_FILTER_BOUNDS\n"));
				for(final double z:GC_WIDTHS){allDataGcLossRow(losses, data, stats[arm], arm, "SD", z, false, true);}
				allDataGcLossRow(losses, data, stats[arm], arm, "LOSSLESS_SD", losslessZ, false, infinite==0 && finite>0);
				for(final double m:new double[]{0,0.01,0.02,0.05}){allDataGcLossRow(losses, data, stats[arm], arm, "MARGIN", m, true, true);}
				allDataGcLossRow(losses, data, stats[arm], arm, "LOSSLESS_MARGIN", margin, true, finite+infinite>0);
			}
		}
		System.err.println("NCRNA_DDL_ALL_DATA_GC_LOSSES_COMPLETE_NO_HOLDOUTS_OR_DDL_SCORING");
	}
	private static String[] gcPowerOutputs(final Config c){
		assert(c!=null && gcAllPower(c)) : "GC power has its own fresh output family and never replaces the accepted fit table";
		return c.fixedCheck ? new String[]{c.output+".gc_power.tsv", c.output+".gc_power_histogram.tsv", c.output+".gc_survival.tsv"}
			: new String[]{c.output+".gc_power.tsv", c.output+".gc_power_histogram.tsv", c.output+".gc_survival.tsv",
				c.output+".gc_model_limits.tsv", c.output+".gc_model_quantiles.tsv", c.output+".gc_model_histogram.tsv",
				c.output+".gc_distances.tsv", c.output+".gc_distance_links.tsv", c.output+".gc_outliers.tsv"};
	}

	/** Gate work is independent of DDL bins,seeds and ranking assignments.
	 * The fixed invocation reads its original selection,not NOT_SEEN labels from
	 * the expanded full population. Both invocations load the same accepted fit. */
	private void allDataGcPower(final NcrnaDdlMembers members, final String[] headers, final NcrnaDdlSavedData data){
		assert(excludedMembers==null && inputSides==null) : "GC power loads the frozen population without constructing a ranking fold";
		data.readNativeCounts(headers, null, true);
		readAllGcStatistics(members);
		if(!c.fixedCheck){gcModelLimits(headers, data);}
		require(Double.isFinite(allGcLosslessZ[0]), "The requested power table includes the accepted finite exact window lossless SD bound");
		final int settings=GC_POWER_Z.length+GC_ALL_MARGINS.length+2;
		final double[] parameters=new double[settings];
		final String[] kinds=new String[settings];
		kinds[0]="NONE";
		for(int i=0; i<GC_POWER_Z.length; i++){
			parameters[i+1]=GC_POWER_Z[i];
			kinds[i+1]="SD";
		}
		for(int i=0; i<GC_ALL_MARGINS.length; i++){
			parameters[GC_POWER_Z.length+i+1]=GC_ALL_MARGINS[i];
			kinds[GC_POWER_Z.length+i+1]="MARGIN";
		}
		parameters[settings-1]=allGcLosslessZ[0];
		kinds[settings-1]="LOSSLESS_SD";
		final long[][][] histogram=new long[2][settings][c.models+1], totals=new long[2][settings][6];
		final boolean[][] survived=new boolean[settings][data.parents.size()];
		final int[][] minimumKept=new int[settings][data.parents.size()];
		for(final int[] row:minimumKept){Arrays.fill(row, c.models+1);}
		final boolean[] seen=new boolean[data.parents.size()], kept=new boolean[c.models];
		for(final Window window:data.windows){
			boolean judged=false;
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){judged|=replicateParent(edge.parent);}
			for(int setting=0; setting<settings; setting++){
				final boolean margin=kinds[setting].equals("MARGIN");
				final double parameter=parameters[setting];
				int removed=0, acceptable=0, retained=0;
				for(int model=0; model<c.models; model++){
					kept[model]=setting==0 || gcWithin(window.gc, allGcStatistics[0][model], parameter, margin);
					if(!kept[model]){removed++;}
					if(window.acceptedModelBits!=null && (window.acceptedModelBits[model>>>6]&(1L<<(model&63)))!=0){
						acceptable++;
						if(kept[model]){retained++;}
					}
				}
				for(int population=0; population<(judged ? 2 : 1); population++){
					histogram[population][setting][removed]++;
					final long[] row=totals[population][setting];
					row[0]++;
					row[1]+=removed;
					row[2]+=acceptable;
					row[3]+=acceptable-retained;
					if(acceptable>0){
						row[4]++;
						if(retained==0){row[5]++;}
					}
				}
				for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
					if(!replicateParent(edge.parent)){continue;}
					final int parent=edge.parent.index;
					seen[parent]=true;
					minimumKept[setting][parent]=Math.min(minimumKept[setting][parent], c.models-removed);
					for(int word=0; word<edge.uniqueModels.length; word++){
						for(long bits=edge.uniqueModels[word]; bits!=0; bits&=bits-1){
							if(kept[(word<<6)+Long.numberOfTrailingZeros(bits)]){survived[setting][parent]=true;}
						}
					}
				}
			}
		}
		try(Output power=new Output(c.output+".gc_power.tsv"); Output bins=new Output(c.output+".gc_power_histogram.tsv");
				Output loci=new Output(c.output+".gc_survival.tsv")){
			power.writer.println(GC_POWER_HEADER);
			bins.writer.println("cohort\tfilter\tparameter\tmodels_removed\twindows");
			loci.writer.println("cohort\tfilter\tparameter\tpanel\tinput_slot\tinput_id\tparent_id\tacceptable_model_survives\tminimum_models_retained_across_saved_windows");
			final String cohort=c.fixedCheck ? "FIXED_NOT_SEEN" : "REBUILT_FIVE_FOLD";
			for(int setting=0; setting<settings; setting++){
				final String kind=kinds[setting];
				final double parameter=parameters[setting];
				long parents=0, lost=0, failedCandidates=0;
				for(final Parent parent:data.parents){
					if(!replicateParent(parent)){continue;}
					require(seen[parent.index], "Power retains every original comparable locus across all its saved windows");
					parents++;
					if(!survived[setting][parent.index]){
						lost++;
						failedCandidates+=minimumKept[setting][parent.index];
					}
					loci.writer.print(line.clear().append(cohort).tab().append(kind).tab().append(Double.toString(parameter)).tab()
						.append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab().append(parent.input.idBytes).tab()
						.append(parent.idBytes).tab().append(survived[setting][parent.index]).tab().append(minimumKept[setting][parent.index]).nl());
				}
				require(setting!=0 || lost==0, "The no-filter baseline retains every original comparable locus and cannot charge a failure sentinel");
				for(int population=0; population<2; population++){
					final String name=population==1 ? cohort : c.fixedCheck ? "FIXED_ALL_SAVED_WINDOWS" : "FULL_ALL_SAVED_WINDOWS";
					final long[] row=totals[population][setting];
					long windows=0, removed=0;
					for(int count=0; count<=c.models; count++){
						final long frequency=histogram[population][setting][count];
						windows+=frequency;
						removed+=frequency*count;
						bins.writer.print(line.clear().append(name).tab().append(kind).tab().append(Double.toString(parameter))
							.tab().append(count).tab().append(frequency).nl());
					}
					require(windows==row[0] && removed==row[1] && windows>0,
						"The complete0..models histogram conserves the explicit window denominator and removed-model sum");
					final long baseline=windows*c.models;
					final double medianRemoved=gcRemovedMedian(histogram[population][setting]);
					power.writer.print(line.clear().append(name).tab().append(kind).tab().append(Double.toString(parameter))
						.tab().append(c.models).tab().append(windows).tab().append(Double.toString((double)removed/windows))
						.tab().append(Double.toString(medianRemoved)).tab().append(row[2])
						.tab().append(row[3]).tab().append(row[4]).tab().append(row[5]).tab().append(parents).tab().append(lost)
						.tab().append(baseline).tab().append(baseline-removed).tab().append(Double.toString((double)(baseline-removed)/windows))
						.tab().append(Double.toString(c.models-medianRemoved)).tab().append(Double.toString(100.0*removed/baseline))
						.tab().append(failedCandidates)
						.append("\tALL_DATA_WINDOW_GC_FIXED_STATISTICS_NOT_CALLER_ACCURACY\n"));
				}
			}
		}
		if(!c.fixedCheck){writeGcOutliers(headers, data, survived, parameters, kinds);}
		System.err.println("NCRNA_DDL_GC_POWER_COMPLETE_SHARED_ALL_DATA_STATISTICS_NO_DDL_SCORING");
	}

	/** Ordinary midpoint median over a complete count histogram,including zero. */
	public static double gcRemovedMedian(final long[] histogram){
		require(histogram!=null && histogram.length>1, "Removed-model median needs the explicit0..models count domain");
		long n=0;
		for(final long count:histogram){
			require(count>=0, "A window histogram cannot contain a negative frequency");
			n+=count;
		}
		require(n>0, "An empty window population has no removed-model median");
		final long first=(n+1)/2, second=n/2+1;
		long cumulative=0;
		int left=-1;
		for(int removed=0; removed<histogram.length; removed++){
			cumulative+=histogram[removed];
			if(left<0 && cumulative>=first){left=removed;}
			if(cumulative>=second){return (left+removed)*0.5;}
		}
		throw new IllegalArgumentException("The removed-model histogram must contain both central windows");
	}

	private void gcModelLimits(final String[] headers, final NcrnaDdlSavedData data){
		assert(!c.fixedCheck && allGcStatistics!=null) : "Per-model limits use every full-data accepted window and the unchanged accepted fit";
		final long[] windows=new long[c.models], finite=new long[c.models];
		final double[][] maximum=new double[2][c.models];
		gcPairDistances(headers, data, windows, finite, maximum);
		try(Output models=new Output(c.output+".gc_model_limits.tsv"); Output quantiles=new Output(c.output+".gc_model_quantiles.tsv");
				Output histogram=new Output(c.output+".gc_model_histogram.tsv")){
			models.writer.println("model_ordinal\tmodel_id\taccepted_windows\tfinite_gc_windows\tmissing_gc_windows\tmaximum_standardized_distance\tinclusive_lossless_z\twindow_count_stratum\tstatus");
			for(int model=0; model<c.models; model++){
				require(finite[model]==allGcWindowCounts[model], "Own-window limits must retain the original full-data GC observation census for every model");
				models.writer.print(line.clear().append(model).tab().append(headers[model].getBytes(java.nio.charset.StandardCharsets.UTF_8))
					.tab().append(windows[model]).tab().append(finite[model]).tab().append(windows[model]-finite[model]).tab()
					.append(finite[model]==0 ? "NA" : Double.isInfinite(maximum[0][model]) ? "INF" : Double.toString(maximum[0][model])).tab()
					.append(finite[model]==0 ? "NA" : Double.isInfinite(maximum[1][model]) ? "INF" : Double.toString(maximum[1][model]))
					.tab().append(GC_WINDOW_STRATA[gcWindowStratum(windows[model])])
					.tab().append(windows[model]==0 ? "NO_WINDOWS" : finite[model]==0 ? "NO_FINITE_GC" : "OWN_ACCEPTED_WINDOWS").nl());
			}
			quantiles.writer.println("statistic\twindow_count_stratum\tmodels\tfinite_models\tinfinite_models\tmissing_models\tmin\tp10\tp25\tp50\tp75\tp90\tp99\tmax\tinterpretation");
			histogram.writer.println("statistic\twindow_count_stratum\tlower_inclusive\tupper_exclusive\tmodels\tinterpretation");
			for(int metric=0; metric<2; metric++){
				for(int stratum=-1; stratum<4; stratum++){
					final double[] values=new double[c.models];
					int n=0, population=0, infinite=0, missing=0;
					for(int model=0; model<c.models; model++){
						if(stratum<0 ? windows[model]==0 : gcWindowStratum(windows[model])!=stratum){continue;}
						population++;
						if(finite[model]==0){missing++;}
						else if(Double.isInfinite(maximum[metric][model])){infinite++;}
						else{values[n++]=maximum[metric][model];}
					}
					Arrays.sort(values, 0, n);
					final String name=metric==0 ? "MAXIMUM_STANDARDIZED_DISTANCE" : "INCLUSIVE_LOSSLESS_Z";
					final String group=stratum<0 ? "ALL_WITH_WINDOWS" : GC_WINDOW_STRATA[stratum];
					line.clear().append(name).tab().append(group).tab().append(population).tab().append(n).tab().append(infinite).tab().append(missing);
					for(final int percentile:new int[]{0,10,25,50,75,90,99,100}){
						final int ordinal=(int)Math.max(0, ((long)(n+infinite)*percentile+99)/100-1);
						line.tab().append(n+infinite==0 ? "NA" : ordinal>=n ? "INF" : Double.toString(values[ordinal]));
					}
					quantiles.writer.print(line.append("\tNEAREST_RANK_INCLUDING_INF_MISSING_SEPARATE\n"));
					final int bins=n==0 ? 0 : (int)Math.floor(values[n-1]*2)+1;
					int start=0;
					for(int bin=0; bin<bins; bin++){
						final double lower=bin*0.5;
						int end=start;
						while(end<n && values[end]<lower+0.5){end++;}
						histogram.writer.print(line.clear().append(name).tab().append(group).tab().append(Double.toString(lower)).tab()
							.append(Double.toString(lower+0.5)).tab().append(end-start).append("\tHALF_SD_BINS_ZERO_COUNTS_INCLUDED\n"));
						start=end;
					}
					histogram.writer.print(line.clear().append(name).tab().append(group).append("\tINF\tINF\t").append(infinite).append("\tINFINITE_DISTANCES\n"));
					histogram.writer.print(line.clear().append(name).tab().append(group).append("\tNA\tNA\t").append(missing).append("\tNO_FINITE_GC_OBSERVATIONS\n"));
				}
			}
		}
	}
	/** Each accepted pair is measured once;parent links are a separate relation. */
	private void gcPairDistances(final String[] headers, final NcrnaDdlSavedData data,
			final long[] windows, final long[] finite, final double[][] maximum){
		assert(!c.fixedCheck && allGcStatistics!=null) : "The pair ledger describes the original full accepted-window population and frozen all-data statistics";
		int pairs=0;
		for(final Window window:data.windows){pairs=Math.addExact(pairs, window.acceptedModels);}
		gcPairOffsets=new int[data.windows.size()+1];
		gcPairModels=new int[pairs];
		gcPairDistances=new double[pairs];
		final byte[][] names=new byte[c.models][];
		final byte[][][] text=new byte[c.models][4][];
		for(int model=0; model<c.models; model++){
			names[model]=headers[model].getBytes(java.nio.charset.StandardCharsets.UTF_8);
			for(int field=0; field<4; field++){
				text[model][field]=Double.toString(allGcStatistics[0][model][field]).getBytes(java.nio.charset.StandardCharsets.US_ASCII);
			}
		}
		int pair=0;
		try(Output ledger=new Output(c.output+".gc_distances.tsv"); Output links=new Output(c.output+".gc_distance_links.tsv")){
			ledger.writer.println("pair_ordinal\tpanel\tinput_slot\tinput_id\tgenome_path\twindow_event\tcontig_ordinal\tcontig_id\twindow_start1\twindow_end1\tstrand\tmodel_ordinal\tmodel_id\twindow_gc\tmodel_mean_gc\tmodel_sd_gc\tmodel_min_gc\tmodel_max_gc\tmodel_window_count\tsd_distance\tpasses_4_sd\tpasses_4_25_sd\tpasses_5_sd\tpasses_8_58_sd");
			links.writer.println("pair_ordinal\tparent_ordinal\tlocus_id\tlocus_start1\tlocus_end1\tcomparable_locus\tmodel_uniquely_acceptable_for_locus");
			for(int w=0; w<data.windows.size(); w++){
				final Window window=data.windows.get(w);
				gcPairOffsets[w]=pair;
				if(window.acceptedModelBits!=null){
					final String windowGc=Double.isNaN(window.gc) ? "NA" : Double.toString(window.gc);
					for(int word=0; word<window.acceptedModelBits.length; word++){
						for(long bits=window.acceptedModelBits[word]; bits!=0; bits&=bits-1){
							final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
							final double[] stats=allGcStatistics[0][model];
							windows[model]++;
							double z=Double.NaN;
							if(!Double.isNaN(window.gc)){
								finite[model]++;
								z=gcDistance(window.gc, stats[0], stats[1]);
								maximum[0][model]=Math.max(maximum[0][model], z);
								double inclusive=z;
								if(Double.isFinite(inclusive)){while(!gcWithin(window.gc, stats, inclusive, false)){inclusive=Math.nextUp(inclusive);}}
								maximum[1][model]=Math.max(maximum[1][model], inclusive);
							}
							require(pair<pairs, "Accepted-pair bits cannot exceed the native per-window accepted-model census");
							gcPairModels[pair]=model;
							gcPairDistances[pair]=z;
							line.clear().append(pair).tab();
							gcWindowIdentity(line, window);
							line.tab().append(model).tab().append(names[model]).tab().append(windowGc);
							for(final byte[] field:text[model]){line.tab().append(field);}
							line.tab().append(allGcWindowCounts[model]).tab().append(Double.isNaN(z) ? "NA" : Double.isInfinite(z) ? "INF" : Double.toString(z));
							for(final double limit:GC_LOOKUP_Z){line.tab().append(gcWithin(window.gc, stats, limit, false));}
							ledger.writer.print(line.nl());
							for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
								final Parent parent=edge.parent;
								final boolean unique=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;
								links.writer.print(line.clear().append(pair).tab().append(parent.index).tab().append(parent.idBytes)
									.tab().append(parent.start).tab().append(parent.stop).tab().append(NcrnaDdlFolds.comparable(parent)).tab().append(unique).nl());
							}
							pair++;
						}
					}
				}
				gcPairOffsets[w+1]=pair;
			}
		}
		require(pair==pairs, "One ledger row per accepted window/model pair must conserve the native count independently of locus-link multiplicity");
	}
	private static void gcWindowIdentity(final ByteBuilder out, final Window window){
		assert(out!=null && window!=null) : "Ledger and case rows share stable original input,event,contig and one-based coordinate identity";
		out.append(window.input.panel.nameBytes).tab().append(window.input.slot).tab().append(window.input.idBytes).tab()
			.append(window.input.pathBytes).tab().append(window.event).tab().append(window.ordinal).tab().append(window.seqidBytes)
			.tab().append(window.start+1).tab().append(window.stop+1).tab().append(window.strand==0 ? "+" : "-");
	}

	/** Lookup only:all distances below were already computed for the pair ledger. */
	private void writeGcOutliers(final String[] headers, final NcrnaDdlSavedData data,
			final boolean[][] survived, final double[] parameters, final String[] kinds){
		assert(gcPairOffsets!=null && !c.fixedCheck) : "Outlier lookup consumes the one full-data pair-distance ledger,not a new scoring run";
		final byte[][][] taxonomy=gcRecordedTaxonomy(data);
		try(Output out=new Output(c.output+".gc_outliers.tsv")){
			out.writer.println("case_type\tsd_limit\tpair_ordinal\tpanel\tinput_slot\tinput_id\tgenome_path\twindow_event\tcontig_ordinal\tcontig_id\twindow_start1\twindow_end1\tstrand\tlocus_id\tmodel_id\twindow_gc\tmodel_mean_gc\tmodel_sd_gc\tmodel_min_gc\tmodel_max_gc\tmodel_window_count\tsd_distance\tlocus_status\trecorded_taxid\trecorded_domain\trecorded_phylum\trecorded_genome_label");
			for(final double limit:GC_LOOKUP_Z){
				int setting=-1;
				for(int i=0; i<parameters.length; i++){if(kinds[i].equals("SD") && parameters[i]==limit){setting=i;}}
				require(setting>=0, "Every requested outlier threshold must have the same native power/survival gate");
				for(int w=0; w<data.windows.size(); w++){
					final Window window=data.windows.get(w);
					int closest=-1;
					boolean retained=false;
					for(int pair=gcPairOffsets[w]; pair<gcPairOffsets[w+1]; pair++){
						final int model=gcPairModels[pair];
						if(closest<0 || gcPairDistances[pair]<gcPairDistances[closest]){closest=pair;}
						final boolean kept=gcWithin(window.gc, allGcStatistics[0][model], limit, false);
						retained|=kept;
						if(limit==8.58 && !kept){gcOutlierRow(out, "PAIR", limit, pair, window, null, headers, survived[setting], taxonomy);}
					}
					if(limit!=8.58 && closest>=0 && !retained){
						gcOutlierRow(out, "WINDOW", limit, closest, window, null, headers, survived[setting], taxonomy);
					}
				}
				if(limit!=4){continue;}
				for(final Parent parent:data.parents){
					if(!NcrnaDdlFolds.comparable(parent) || survived[setting][parent.index]){continue;}
					int closest=-1, selectedWindow=-1;
					for(int w=0; w<data.windows.size(); w++){
						final Window window=data.windows.get(w);
						if(window.input!=parent.input){continue;}
						for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
							if(edge.parent!=parent){continue;}
							for(int pair=gcPairOffsets[w]; pair<gcPairOffsets[w+1]; pair++){
								final int model=gcPairModels[pair];
								if((edge.uniqueModels[model>>>6]&(1L<<(model&63)))==0){continue;}
								if(closest<0 || gcPairDistances[pair]<gcPairDistances[closest]){
									closest=pair;
									selectedWindow=w;
								}
							}
						}
					}
					require(closest>=0, "A lost comparable locus still has its original accepted model/window pair in the ledger");
					gcOutlierRow(out, "LOCUS", limit, closest, data.windows.get(selectedWindow), parent, headers, survived[setting], taxonomy);
				}
			}
		}
	}
	private void gcOutlierRow(final Output out, final String kind, final double limit, final int pair,
			final Window window, final Parent selected, final String[] headers, final boolean[] survived, final byte[][][] taxonomy){
		assert(pair>=0 && pair<gcPairModels.length) : "Each exceptional case cites an actual retained pair-ledger ordinal";
		final int model=gcPairModels[pair];
		final double[] stats=allGcStatistics[0][model];
		final ByteBuilder status=new ByteBuilder();
		for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
			final Parent parent=edge.parent;
			if(!NcrnaDdlFolds.comparable(parent) || selected!=null && parent!=selected){continue;}
			if(status.length()>0){status.append("; ");}
			status.append(parent.idBytes).append(survived[parent.index] ? " survives" : " lost");
			if(survived[parent.index]){
				boolean sameWindow=false;
				for(int word=0; word<edge.uniqueModels.length; word++){
					for(long bits=edge.uniqueModels[word]; bits!=0; bits&=bits-1){
						sameWindow|=gcWithin(window.gc, allGcStatistics[0][(word<<6)+Long.numberOfTrailingZeros(bits)], limit, false);
					}
				}
				status.append(sameWindow ? " through another acceptable model in this window" : " through another saved window");
			}
		}
		if(status.length()==0){status.append("No comparable locus linked to this saved window");}
		line.clear().append(kind).tab().append(Double.toString(limit)).tab().append(pair).tab();
		gcWindowIdentity(line, window);
		line.tab().append(selected==null ? new byte[]{'.'} : selected.idBytes).tab()
			.append(headers[model].getBytes(java.nio.charset.StandardCharsets.UTF_8)).tab().append(Double.toString(window.gc));
		for(final double value:stats){line.tab().append(Double.toString(value));}
		line.tab().append(allGcWindowCounts[model]).tab().append(Double.toString(gcPairDistances[pair])).tab().append(status);
		final byte[][] tax=taxonomy[window.input.index];
		for(int field=0; field<4; field++){line.tab().append(tax==null ? "not recorded in supplied genome manifest".getBytes(java.nio.charset.StandardCharsets.US_ASCII) : tax[field]);}
		out.writer.print(line.nl());
	}
	private byte[][][] gcRecordedTaxonomy(final NcrnaDdlSavedData data){
		assert(gcPairOffsets!=null) : "Recorded taxonomy is looked up only for genomes with an exceptional accepted pair";
		final byte[][][] result=new byte[data.selected.size()][][];
		if(c.genomes==null){return result;}
		final boolean[] needed=new boolean[result.length];
		final IntList inputs=new IntList();
		for(int w=0; w<data.windows.size(); w++){
			final Window window=data.windows.get(w);
			for(int pair=gcPairOffsets[w]; pair<gcPairOffsets[w+1]; pair++){
				if(!gcWithin(window.gc, allGcStatistics[0][gcPairModels[pair]], 4, false) && !needed[window.input.index]){
					needed[window.input.index]=true;
					inputs.add(window.input.index);
				}
			}
		}
		if(inputs.size==0){return result;}
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(c.genomes)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, c.genomes);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==11, "Taxonomy lookup preserves the recorded genome/taxid/domain/phylum/fasta schema");
				for(int i=0; i<inputs.size; i++){
					final int input=inputs.get(i);
					if(!p.termEquals(data.selected.get(input).path, 7)){continue;}
					final byte[][] tax={p.parseByteArray(3),p.parseByteArray(4),p.parseByteArray(5),p.parseByteArray(1)};
					require(result[input]==null || Arrays.deepEquals(result[input], tax), "A literal genome path cannot silently acquire conflicting recorded taxonomy");
					result[input]=tax;
				}
			}
		}
		return result;
	}

	private static int gcWindowStratum(final long windows){
		assert(windows>=0) : "Accepted window counts determine the requested model-count strata before inspecting distances";
		return windows==0 ? 0 : windows<10 ? 1 : windows<30 ? 2 : 3;
	}

	/** Finite z uses the same inclusive arithmetic as the ranking filter.
	 * Lossless limits are increased by a representable step only when rounding
	 * would otherwise reject an actual accepted boundary observation. */
	public static boolean gcWithin(final double value, final double[] stats, final double parameter, final boolean margin){
		assert(stats!=null && stats.length==4 && Double.isFinite(parameter) && parameter>=0)
			: "A GC bound uses an actual finite nonnegative SD width or GC-fraction margin;infinite lossless widths are unavailable";
		if(Double.isNaN(value) || Double.isNaN(stats[0])){return true;}
		final double lower=margin ? stats[2]-parameter : stats[0]-parameter*stats[1];
		final double upper=margin ? stats[3]+parameter : stats[0]+parameter*stats[1];
		return value>=lower && value<=upper;
	}
	private void allDataGcLossRow(final Output out, final NcrnaDdlSavedData data, final double[][] stats,
			final int arm, final String kind, final double parameter, final boolean margin, final boolean available){
		assert(stats.length==c.models && arm>=0 && arm<2) : "All-data loss counts retain both complete model rosters and original locus denominators";
		final boolean[] kept=new boolean[c.models], seen=new boolean[data.parents.size()], survived=new boolean[seen.length];
		long pairs=0, removed=0, windows=0, lostWindows=0, loci=0, lostLoci=0;
		for(final Window window:data.windows){
			int accepted=0, retained=0;
			for(int model=0; model<c.models; model++){
				kept[model]=!available || gcWithin(window.gc, stats[model], parameter, margin);
				if(window.acceptedModelBits==null || (window.acceptedModelBits[model>>>6]&(1L<<(model&63)))==0){continue;}
				accepted++;
				if(kept[model]){retained++;}
			}
			pairs+=accepted;
			removed+=accepted-retained;
			if(accepted>0){windows++; if(retained==0){lostWindows++;}}
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
				if(!NcrnaDdlFolds.comparable(edge.parent)){continue;}
				final int parent=edge.parent.index;
				seen[parent]=true;
				for(int word=0; word<edge.uniqueModels.length; word++){
					for(long bits=edge.uniqueModels[word]; bits!=0; bits&=bits-1){
						if(kept[(word<<6)+Long.numberOfTrailingZeros(bits)]){survived[parent]=true;}
					}
				}
			}
		}
		for(final Parent parent:data.parents){
			if(!NcrnaDdlFolds.comparable(parent)){continue;}
			require(seen[parent.index], "Every original comparable locus must retain a saved window in all-data GC losses");
			loci++;
			if(!survived[parent.index]){lostLoci++;}
		}
		if(kind.startsWith("LOSSLESS") && available){require(removed==0 && lostWindows==0 && lostLoci==0,
			"A reported all-data lossless bound must preserve every accepted model/window pair and its loci");}
		line.clear().append(GC_ARMS[arm]).tab().append(kind).tab().append(Double.isInfinite(parameter) ? "INF" : Double.toString(parameter))
			.tab().append(available).tab().append(pairs).tab().append(available ? Long.toString(removed) : "NA")
			.tab().append(windows).tab().append(available ? Long.toString(lostWindows) : "NA")
			.tab().append(loci).tab().append(available ? Long.toString(lostLoci) : "NA");
		out.writer.print(line.append("\tALL_DATA_NO_HELDOUTS_NOT_CALLER_ACCURACY\n"));
	}

	/** Loss curves depend on GC and saved assignments, never on DDL hashes or bins. */
	private void gcLossCurves(final NcrnaDdlMembers members, final String[] headers,
			final byte[][] consensuses, final NcrnaDdlSavedData data, final NcrnaDdlFolds folds){
		final byte[][] names=new byte[headers.length][];
		for(int i=0; i<names.length; i++){names[i]=headers[i].getBytes(java.nio.charset.StandardCharsets.UTF_8);}
		final NcrnaDdlCohorts cohorts=NcrnaDdlCohorts.read(c.identities, data, names, consensuses);
		if(!c.fixedCheck){folds.prepare(members, cohorts);}
		data.readNativeCounts(headers, null, true);
		final LineParser1 p=new LineParser1('\t');
		final LongHashSet assignments=new LongHashSet(16);
		int observations=0;
		try(Input roster=new Input(c.gcFits); Output curve=new Output(c.output+".gc_curve.tsv");
				Output choices=new Output(c.output+".gc_choices.tsv")){
			RrnaResourceIO.header(roster.reader, "assignment_seed\tddl_seed\ttable", c.gcFits);
			curve.writer.println("assignment_seed\tgc_source\tminimum_observations\tz\tfit_pairs_lost\tjudging_pairs_lost\tjudging_windows_lost\tjudging_loci_lost\tjudging_low_loci_lost\tinterpretation");
			choices.writer.println("assignment_seed\tgc_source\tminimum_observations\tfold\tselection\tz\tfit_pairs_lost\tjudging_pairs_lost\tjudging_windows_lost\tjudging_loci_lost\tjudging_low_loci_lost\tgate_models\twindow_fallback_below_minimum\tinterpretation");
			for(byte[] row; (row=roster.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==3, "The accepted roster binds assignment,hash seed and native table path");
				final long assignment=p.parseLong(0);
				if((assignment<0)!=c.fixedCheck || p.parseLong(1)!=12345){continue;}
				require(assignments.add(assignment), "GC-only losses use each recorded assignment once, without hash-seed duplicates");
				final String table=p.parseString(2), suffix=".replicates.tsv";
				require(table.endsWith(suffix), "The accepted worker must retain its original table prefix");
				final int[] original=readGcAssignment(table.substring(0, table.length()-suffix.length())+".folds.tsv", assignment, data);
				final int count=c.fixedCheck ? 2 : c.folds;
				if(!c.fixedCheck){folds.useSavedAssignment(original, count);}
				final long[][][][] pooled=new long[2][2][5][GC_WIDTHS.length+1];
				for(int judge=0; judge<count; judge++){
					inputSides=new int[original.length];
					for(int i=0; i<original.length; i++){inputSides[i]=original[i]==judge ? 1 : 0;}
					fitSide=0;
					excludedMembers=c.fixedCheck ? null : folds.exclusions(original, judge, members);
					final long[][] rawCounts=new long[2][c.models], fallback=new long[2][4];
					final double[][][] effective=fitGc(judge, members, data, null, fallback, rawCounts);
					final long[][][][] losses=new long[2][2][][];
					for(int arm=0; arm<2; arm++){
						for(int policy=0; policy<2; policy++){
							final int minimum=policy==0 ? 0 : 30;
							losses[arm][policy]=gcLossHistogram(data, cohorts, folds, effective[arm], rawCounts[arm], minimum);
							int gated=0, fallbackBelow=0;
							for(int model=0; model<c.models; model++){
								if(rawCounts[arm][model]>=minimum && !Double.isNaN(effective[arm][model][0])){gated++;}
								if(arm==0 && rawCounts[0][model]==0 && rawCounts[0][model]<minimum){fallbackBelow++;}
							}
							int selected=0;
							while(selected<GC_WIDTHS.length && losses[arm][policy][0][selected]>0){selected++;}
							writeGcChoice(choices, assignment, arm, minimum, judge, "FIT_ZERO_PAIR_LOSS", selected,
								losses[arm][policy], gated, fallbackBelow);
							for(int metric=0; metric<5; metric++){
								for(int z=0; z<GC_WIDTHS.length; z++){pooled[arm][policy][metric][z]+=losses[arm][policy][metric][z];}
							}
						}
					}
					checkGcLossParity(table, judge, losses);
					System.err.println("NCRNA_DDL_GC_LOSS_FOLD_COMPLETE\tassignment="+assignment+"\tfold="+judge);
				}
				for(int arm=0; arm<2; arm++){
					for(int policy=0; policy<2; policy++){
						final int minimum=policy==0 ? 0 : 30;
						for(int z=0; z<GC_WIDTHS.length; z++){
							line.clear().append(assignment).tab().append(GC_ARMS[arm]).tab().append(minimum).tab().append(GC_WIDTHS[z], 2);
							for(int metric=0; metric<5; metric++){line.tab().append(pooled[arm][policy][metric][z]);}
							curve.writer.print(line.append("\tDISJOINT_JUDGING_FOLDS_POOLED_FIT_PAIRS_REPEATED_TRAINING_EXPOSURES\n"));
						}
						for(final int allowed:new int[]{1,5,10}){
							int selected=0;
							while(selected<GC_WIDTHS.length && pooled[arm][policy][3][selected]>allowed){selected++;}
							writeGcChoice(choices, assignment, arm, minimum, -1, "DESCRIPTIVE_JUDGING_LOSS_AT_MOST_"+allowed,
								selected, pooled[arm][policy], -1, -1);
						}
					}
				}
				observations++;
			}
		}
		require(observations==(c.models==508 ? c.fixedCheck ? 1 : 10 : 1),
			"GC loss diagnostics must cover the complete predeclared assignment roster once");
		System.err.println("NCRNA_DDL_GC_LOSS_CURVES_AND_EXACT_PARITY_PASS\tassignments="+observations);
	}

	private int[] readGcAssignment(final String path, final long assignment, final NcrnaDdlSavedData data){
		final int[] result=new int[data.selected.size()];
		Arrays.fill(result, -1);
		final IntHashMap byGroup=new IntHashMap();
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "input_index\tassembly_group\tjudge_fold\tassignment_seed\tcomparable_group", path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5 && p.parseInt(0)==rows && rows<result.length,
					"Saved assignments retain every selected input exactly once in its original order");
				final int group=p.parseInt(1), fold=p.parseInt(2);
				require(group==data.selected.get(rows).group && fold>=0 && fold<(c.fixedCheck ? 2 : c.folds)
					&& p.parseLong(3)==(c.fixedCheck ? 20261113 : assignment),
					"Saved fold reuse cannot change assembly bindings,declared folds or the assignment seed");
				require(!byGroup.contains(group) || byGroup.get(group)==fold, "One assembly group must remain wholly in its recorded fold");
				byGroup.put(group, fold);
				result[rows++]=fold;
			}
		}
		require(rows==result.length, "A truncated saved assignment would silently change the fit/judging populations");
		return result;
	}

	/** First-passing-width histograms count each accepted pair once and each locus
	 * once across all its uniquely acceptable models and saved windows. */
	private long[][] gcLossHistogram(final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts,
			final NcrnaDdlFolds folds, final double[][] stats, final long[] observations, final int minimum){
		assert(stats.length==c.models && observations.length==c.models) : "Loss curves preserve the same complete model roster as the accepted scoring study";
		final int end=GC_WIDTHS.length;
		final long[][] hist=new long[5][end+1];
		final int[] parents=new int[data.parents.size()], first=new int[c.models];
		final boolean[] seen=new boolean[parents.length];
		Arrays.fill(parents, end);
		for(final Window window:data.windows){
			final boolean judging=inputSides[window.input.index]==1;
			boolean comparable=false;
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){comparable|=replicateParent(edge.parent);}
			if(judging && !comparable){continue;}
			Arrays.fill(first, -1);
			int best=end, accepted=0;
			if(window.acceptedModelBits!=null){
				for(int word=0; word<window.acceptedModelBits.length; word++){
					for(long bits=window.acceptedModelBits[word]; bits!=0; bits&=bits-1){
						final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
						final int z=gcFirstWidth(window.gc, stats[model], observations[model], minimum);
						first[model]=z;
						best=Math.min(best, z);
						hist[judging ? 1 : 0][z]++;
						accepted++;
					}
				}
			}
			if(!judging){continue;}
			if(accepted>0){hist[2][best]++;}
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
				if(!replicateParent(edge.parent)){continue;}
				final int parent=edge.parent.index;
				seen[parent]=true;
				for(int word=0; word<edge.uniqueModels.length; word++){
					for(long bits=edge.uniqueModels[word]; bits!=0; bits&=bits-1){
						final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
						if(first[model]<0){first[model]=gcFirstWidth(window.gc, stats[model], observations[model], minimum);}
						parents[parent]=Math.min(parents[parent], first[model]);
					}
				}
			}
		}
		for(final Parent parent:data.parents){
			if(inputSides[parent.input.index]!=1 || !replicateParent(parent)){continue;}
			require(seen[parent.index], "Every comparable judging locus must retain a saved window in the loss denominator");
			hist[3][parents[parent.index]]++;
			if(c.fixedCheck ? cohorts.bins[parent.index]==0 : folds.low[parent.index]){hist[4][parents[parent.index]]++;}
		}
		for(final long[] values:hist){
			// A pair with first-pass index j is lost at every tested index below j.
			long lost=values[end];
			for(int z=end-1; z>=0; z--){
				final long passing=values[z];
				values[z]=lost;
				lost+=passing;
			}
		}
		return hist;
	}

	/** Uses the same inclusive mean±z*SD comparisons as gcKeep; no division,
	 * epsilon or rounded fit-table values can change an exact boundary decision. */
	public static int gcFirstWidth(final double value, final double[] stats, final long observations, final int minimum){
		assert(stats!=null && stats.length==4 && observations>=0 && minimum>=0) : "A gate requires its own source count,not a fallback source count";
		if(observations<minimum || Double.isNaN(value) || Double.isNaN(stats[0])){return 0;}
		int lo=0, hi=GC_WIDTHS.length;
		while(lo<hi){
			final int middle=(lo+hi)>>>1;
			final double width=GC_WIDTHS[middle]*stats[1];
			if(value>=stats[0]-width && value<=stats[0]+width){hi=middle;}
			else{lo=middle+1;}
		}
		return lo;
	}

	private void writeGcChoice(final Output out, final long assignment, final int arm, final int minimum,
			final int fold, final String selection, final int selected, final long[][] losses,
			final int gated, final int fallbackBelow){
		assert(selected>=0 && selected<=GC_WIDTHS.length) : "No lossless width through8 is explicitly unavailable,not silently capped";
		line.clear().append(assignment).tab().append(GC_ARMS[arm]).tab().append(minimum).tab().append(fold).tab().append(selection).tab();
		if(selected==GC_WIDTHS.length){line.append("NA\tNA\tNA\tNA\tNA\tNA");}
		else{
			line.append(GC_WIDTHS[selected], 2);
			for(int metric=0; metric<5; metric++){line.tab().append(losses[metric][selected]);}
		}
		out.writer.print(line.tab().append(gated).tab().append(fallbackBelow).tab()
			.append(fold<0 ? "JUDGING_TUNING_ONLY" : "FIT_SELECTED_OPPOSITE_FOLD_JUDGED").nl());
	}

	private void checkGcLossParity(final String table, final int judge, final long[][][][] losses){
		final LineParser1 p=new LineParser1('\t');
		final boolean[][][] seen=new boolean[2][2][3];
		int rows=0;
		try(Input in=new Input(table)){
			RrnaResourceIO.header(in.reader, GC_REPLICATE_HEADER, table);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==39, "Exact loss parity requires the accepted GC table including failure histograms");
				if(p.parseInt(2)!=judge || !p.termEquals("ALL", 27)){continue;}
				for(int bins=0; bins<2; bins++){
					for(int arm=0; arm<2; arm++){
						for(int setting=6; setting<9; setting++){
							if(!p.termEquals(gcMethod(bins==0 ? 512 : 1024, arm*12+setting), 3)){continue;}
							final int z=(setting-5)*4;
							require(!seen[bins][arm][setting-6] && p.parseLong(31)==losses[arm][0][1][z]
								&& p.parseLong(33)==losses[arm][0][2][z]
								&& lastGcHistogram(p.parseString(28))==losses[arm][0][3][z]
								&& lastGcHistogram(p.parseString(29))==losses[arm][0][4][z],
								"GC-only loss curves must exactly match accepted z2/3/4 pair,window,locus and low-locus losses: "+table+", fold="+judge+", arm="+arm+", setting="+setting);
							seen[bins][arm][setting-6]=true;
							rows++;
						}
					}
				}
			}
		}
		require(rows==12, "Exact loss parity needs both bin counts,both GC sources and all three previously measured SD filters");
	}
	private static long lastGcHistogram(final String histogram){
		assert(histogram!=null && histogram.lastIndexOf(',')>=0) : "The final rank-histogram bin is the all-acceptable-filtered sentinel";
		return Long.parseLong(histogram.substring(histogram.lastIndexOf(',')+1));
	}

	private void gcReplicates(final Row base, final int judge, final NcrnaDdlMembers members,
			final byte[][] consensuses, final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts,
			final NcrnaDdlFolds folds, final int empty, final Output out, final Output fitOut){
		assert(base.shortK==8 && base.longK==10 && fitSide==0) : "GC judging preserves the decided forward8/masked10 recipe and opposite fit side";
		final boolean allData=gcAllData(c);
		final int variants=c.gcAdditions ? GC_EXTRA_SETTINGS.length : allData ? GC_ALL_SETTINGS.length : GC_SETTINGS.length;
		final long[][] fallbackCounts=allData ? allGcFallback : new long[2][4], sourceCounts=new long[2][c.models];
		final double[][][] effective=allData ? allGcStatistics : fitGc(judge, members, data, fitOut, fallbackCounts, sourceCounts);
		final long[] shortWords=frequencies(8, false, 2, members, data), longWords=frequencies(10, true, 2, members, data);
		blacklist=new LongHashSet(Math.max(2, c.shortPrefixes[2]+c.longPrefix));
		addBlacklist(shortWords, c.shortPrefixes[2], 8, true);
		addBlacklist(longWords, c.longPrefix, 10, false);
		for(final int bins:new int[]{512,1024}){
			final long cpu=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
			final int bonus=bins==512 ? 60 : 100;
			replicateBonus=bonus;
			replicateScale=0;
			retention=new double[3];
			final Row row=new Row(base, bins, 2);
			final DynamicDemiLog[] references=buildReferences(row, members, consensuses)[0];
			final DynamicDemiLog query=createParts(row, 1)[0];
			final int methods=(c.gcAdditions ? 1 : 2)*variants;
			final int[][] ranks=new int[methods][data.parents.size()];
			for(final int[] values:ranks){Arrays.fill(values, -1);}
			final long[][] counts=new long[methods][9];
			for(int method=0; method<methods; method++){
				System.arraycopy(fallbackCounts[method/variants], 0, counts[method], 5, 4);
			}
			final double[] scores=new double[c.models], original=new double[c.models];
			final double[][] distance=new double[2][c.models];
			final boolean[] kept=new boolean[c.models];
			for(final Window window:data.windows){
				if(inputSides[window.input.index]!=1){continue;}
				boolean judged=false;
				for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){judged|=replicateParent(edge.parent);}
				if(!judged){continue;}
				fill(query, window.bases, row, 0, true);
				for(int model=0; model<c.models; model++){
					if(Shared.SIMD){simd.Vector.compareDDL(query.maxArray(), references[model].maxArray(), comparison);}
					else{DynamicDemiLog.compareDetailed(query.maxArray(), references[model].maxArray(), comparison);}
					original[model]=comparison[1];
					for(int arm=0; arm<2; arm++){
						distance[arm][model]=gcDistance(window.gc, effective[arm][model][0], effective[arm][model][1]);
					}
				}
				for(int setting=0; setting<methods; setting++){
					final int arm=setting/variants, variant=setting%variants;
					int acceptable=0, retained=0;
					for(int model=0; model<c.models; model++){
						kept[model]=c.gcAdditions ? gcExtraKeep(window.gc, effective[0][model], variant)
							: variant<6 || (allData ? gcAllKeep(window.gc, effective[arm][model], variant, allGcLosslessZ[arm])
							: gcKeep(window.gc, effective[arm][model], variant));
						scores[model]=c.gcAdditions ? gcExtraScore(original[model], distance[0][model], variant)
							: variant<6 ? gcScore(original[model], distance[arm][model], GC_WEIGHTS[variant]) : original[model];
						if(window.acceptedModelBits!=null && (window.acceptedModelBits[model>>>6]&(1L<<(model&63)))!=0){
							acceptable++;
							if(kept[model]){retained++;}
						}
					}
					counts[setting][0]+=acceptable;
					counts[setting][1]+=acceptable-retained;
					if(acceptable>0){counts[setting][2]++;}
					if(acceptable>0 && retained==0){counts[setting][3]++;}
					if(Double.isNaN(window.gc)){counts[setting][4]++;}
					for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
						if(!replicateParent(edge.parent)){continue;}
						final int parent=edge.parent.index, rank=gcRank(scores, edge.uniqueModels, kept);
						if(ranks[setting][parent]<0 || rank<ranks[setting][parent]){ranks[setting][parent]=rank;}
					}
				}
			}
			final double seconds=(clock.getCurrentThreadCpuTime()-cpu)/1e9, elapsed=(System.nanoTime()-wall)/1e9;
			require(c.gcAdditions || Arrays.equals(ranks[0], ranks[variants]), "The two GC source arms must have exactly identical no-GC controls");
			for(int setting=0; setting<methods; setting++){
				gcWindowMetrics=counts[setting];
				emitReplicate(c.gcAdditions ? gcExtraMethod(bins, setting) : allData ? gcAllMethod(bins, setting) : gcMethod(bins, setting), row, c.shortPrefixes[2], c.longPrefix, bonus, 0, judge,
					ranks[setting], data, cohorts, folds, empty, retention, setting==0 ? seconds : 0, setting==0 ? elapsed : 0, out);
			}
			System.err.println("NCRNA_DDL_GC_FOLD_COMPLETE\tfold="+judge+"\tbins="+bins);
		}
		gcWindowMetrics=null;
		blacklist=null;
	}
	public static String gcMethod(final int bins, final int setting){
		require((bins==512 || bins==1024) && setting>=0 && setting<2*GC_SETTINGS.length, "GC methods name one decided bin count,statistics source and declared setting");
		return "gc_"+bins+"_"+GC_ARMS[setting/12]+"_"+GC_SETTINGS[setting%12];
	}
	public static String gcAllMethod(final int bins, final int setting){
		require((bins==512 || bins==1024) && setting>=0 && setting<2*GC_ALL_SETTINGS.length,
			"All-data GC methods declare their source,bin count and measured filter setting");
		return "gc_all_"+bins+"_"+GC_ARMS[setting/GC_ALL_SETTINGS.length]+"_"+GC_ALL_SETTINGS[setting%GC_ALL_SETTINGS.length];
	}
	public static boolean gcAllKeep(final double value, final double[] stats, final int setting, final double losslessZ){
		assert(setting>=6 && setting<GC_ALL_SETTINGS.length) : "All-data filtering uses fixed2/3/4/6/8SD,the measured lossless width,or one declared GC margin";
		return gcWithin(value, stats, setting<11 ? GC_ALL_Z[setting-6] : setting==11 ? losslessZ : GC_ALL_MARGINS[setting-12], setting>=12);
	}
	public static String gcExtraMethod(final int bins, final int setting){
		require((bins==512 || bins==1024) && setting>=0 && setting<GC_EXTRA_SETTINGS.length,
			"GC additions name only new window-statistic filters or range-plus-penalty combinations");
		return "gc_all_"+bins+"_window_"+GC_EXTRA_SETTINGS[setting];
	}
	public static boolean gcExtraKeep(final double value, final double[] stats, final int setting){
		assert(setting>=0 && setting<GC_EXTRA_SETTINGS.length) : "The six combinations preserve the original inclusive range gate";
		return gcWithin(value, stats, GC_ALL_MARGINS[setting/3], true);
	}
	public static double gcExtraScore(final double count, final double distance, final int setting){
		assert(setting>=0 && setting<GC_EXTRA_SETTINGS.length) : "Only combined settings change equal-bin scores after filtering";
		return gcScore(count, distance, GC_WEIGHTS[1+setting%3]);
	}
	/** A missing fit/query has no GC evidence. Zero SD is exact, without epsilon. */
	public static double gcDistance(final double value, final double mean, final double sd){
		assert(Double.isNaN(value) || value>=0 && value<=1) : "GC distance uses a fraction and its actual fit-side mean/SD";
		if(Double.isNaN(mean) || Double.isNaN(value)){return 0;}
		return sd==0 ? value==mean ? 0 : Double.POSITIVE_INFINITY : Math.abs(value-mean)/sd;
	}
	public static double gcScore(final double count, final double distance, final double weight){
		assert(Double.isFinite(count) && count>=0 && distance>=0 && Double.isFinite(weight) && weight>=0)
			: "GC ordering subtracts a declared nonnegative penalty from an actual equal-bin count";
		return weight==0 ? count : count-weight*distance;
	}
	public static boolean gcKeep(final double value, final double[] stats, final int setting){
		assert(stats!=null && stats.length==4 && setting>=6 && setting<12) : "GC filters use effective fit-side mean,SD,min,max";
		if(Double.isNaN(stats[0]) || Double.isNaN(value)){return true;}
		final double lower=setting<9 ? stats[0]-(setting-4)*stats[1] : stats[2]-GC_MARGINS[setting-9];
		final double upper=setting<9 ? stats[0]+(setting-4)*stats[1] : stats[3]+GC_MARGINS[setting-9];
		return value>=lower && value<=upper;
	}
	/** Complete ordering permits negative scores and the explicit zero-SD bottom
	 * tier. Filters remove models; no surviving acceptable model returns models+1. */
	public static int gcRank(final double[] scores, final long[] allowed, final boolean[] kept){
		require(scores!=null && scores.length>0 && allowed!=null && allowed.length==(scores.length+63)/64
			&& kept!=null && kept.length==scores.length, "GC ranking preserves the full model roster and acceptable bits");
		int selected=-1;
		for(int model=0; model<scores.length; model++){
			assert(Double.isFinite(scores[model]) || scores[model]==Double.NEGATIVE_INFINITY) : "Only zero-SD mismatches may occupy the deterministic bottom tier; NaN is forbidden";
			if(kept[model] && (allowed[model>>>6]&(1L<<(model&63)))!=0
				&& (selected<0 || scores[model]>scores[selected])){selected=model;}
		}
		if(selected<0){return scores.length+1;}
		int rank=1;
		for(int model=0; model<scores.length; model++){
			if(kept[model] && (scores[model]>scores[selected] || scores[model]==scores[selected] && model<selected)){rank++;}
		}
		return rank;
	}

	private void runReplicate(final String method, final Row row, final int shortSize, final int longSize,
			final int bonus, final double scale, final int judge, final NcrnaDdlMembers members,
			final byte[][] consensuses, final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts,
			final NcrnaDdlFolds folds, final int empty, final Output output){
		final long cpu=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
		replicateBonus=bonus;
		replicateScale=scale;
		retention=new double[3];
		final DynamicDemiLog[] references=buildReferences(row, members, consensuses)[0];
		final DynamicDemiLog query=createParts(row, 1)[0];
		final int[] ranks=new int[data.parents.size()];
		Arrays.fill(ranks, -1);
		final double[] scores=new double[c.models];
		for(final Window window:data.windows){
			if(inputSides[window.input.index]!=1){continue;}
			boolean judged=false;
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){judged|=replicateParent(edge.parent);}
			if(!judged){continue;}
			fill(query, window.bases, row, 0, true);
			for(int model=0; model<c.models; model++){
				if(Shared.SIMD){simd.Vector.compareDDL(query.maxArray(), references[model].maxArray(), comparison);}
				else{DynamicDemiLog.compareDetailed(query.maxArray(), references[model].maxArray(), comparison);}
				scores[model]=comparison[1];
			}
			for(final prok.NcrnaDdlSavedData.Edge edge:window.edges){
				final int parent=edge.parent.index;
				if(!replicateParent(edge.parent)){continue;}
				final int rank=NcrnaPartitionWeightSearch.firstRank(scores, edge.uniqueModels);
				if(rank>0 && (ranks[parent]<0 || rank<ranks[parent])){ranks[parent]=rank;}
			}
		}
		emitReplicate(method, row, shortSize, longSize, bonus, scale, judge, ranks, data, cohorts,
			folds, empty, retention, (clock.getCurrentThreadCpuTime()-cpu)/1e9, (System.nanoTime()-wall)/1e9, output);
		System.err.println("NATIVE_DDL_REPLICATE_COMPLETE\t"+method+"\tfold="+judge+"\tseed="+sketchSeed);
	}

	private boolean replicateParent(final Parent parent){
		return NcrnaDdlFolds.comparable(parent) && (!c.fixedCheck || parent.stratum==0);
	}

	private void emitReplicate(final String method, final Row row, final int shortSize, final int longSize,
			final int bonus, final double scale, final int judge, final int[] ranks,
			final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts, final NcrnaDdlFolds folds,
			final int empty, final double[] retained, final double cpu, final double wall, final Output output){
		final long[][] totals=new long[3][8];
		final int[][][] histograms=rankHistograms;
		if(histograms!=null){for(final int[][] subset:histograms){for(final int[] histogram:subset){Arrays.fill(histogram, 0);}}}
		for(final Parent parent:data.parents){
			if(inputSides[parent.input.index]!=1 || !replicateParent(parent)){continue;}
			final int rank=ranks[parent.index];
			require(rank>0 && rank<=c.models+(gcStudy(c) ? 1 : 0),
				"Every comparable locus retains its rank; GC filtering uses models+1 only for loss of all acceptable models");
			final boolean low=c.fixedCheck ? cohorts.bins[parent.index]==0 : folds.low[parent.index];
			addRank(totals[0], rank, low);
			addRank(totals[folds.originalContributor(parent) ? 1 : 2], rank, low);
			if(histograms!=null){
				final int subset=folds.originalContributor(parent) ? 1 : 2;
				histograms[0][0][rank]++;
				histograms[subset][0][rank]++;
				if(low){
					histograms[0][1][rank]++;
					histograms[subset][1][rank]++;
				}
			}
		}
		require(totals[0][0]>0, "A replicate cannot silently report an empty judging fold");
		for(int metric=0; metric<8; metric++){
			assert(totals[0][metric]==totals[1][metric]+totals[2][metric])
				: "Original contributor and noncontributor loci partition the complete judging population and every rank/count metric";
		}
		for(int stratum=0; stratum<(c.fixedCheck ? 1 : 3); stratum++){
			line.clear().append(c.fixedCheck ? -1 : c.foldSeed).tab().append(method.equals("CONSENSUS_INDEX") ? -1 : sketchSeed)
				.tab().append(judge).tab().append(method).tab().append(row.shortK).tab().append(row.longK)
				.tab().append(PARTITION_NAMES[row.partitionBits]).tab().append(row.memberBins)
				.tab().append(shortSize).tab().append(longSize).tab().append(bonus).tab().append(scale, 9);
			for(final long value:totals[stratum]){line.tab().append(value);}
			line.tab().append(empty);
			for(int third=0; third<3; third++){
				int count=0;
				for(final int size:sizeThirds){if(size==third){count++;}}
				line.tab();
				if(stratum>0 || retained==null || count==0){line.append("NA");}
				else{line.append(retained[third]/count, 9);}
			}
			line.tab().append(stratum==0 ? cpu : 0, 9).tab().append(stratum==0 ? wall : 0, 9)
				.append("\tSAVED_ACCEPTANCE_ORACLE_AND_FULL_CONSENSUS_RETAINED_NOT_CALLER\t").append(CONTRIBUTION[stratum]);
			if(histograms!=null){
				for(final int[] histogram:histograms[stratum]){
					line.tab();
					for(int rank=1; rank<histogram.length; rank++){
						if(rank>1){line.append(',');}
						line.append(histogram[rank]);
					}
				}
			}
			if(gcStudy(c)){
				for(int metric=0; metric<9; metric++){line.tab().append(stratum==0 ? gcWindowMetrics[metric] : 0);}
			}
			line.nl();
			output.writer.print(line);
		}
	}

	private static void addRank(final long[] totals, final int rank, final boolean low){
		assert(totals.length==8 && rank>0) : "Judging summaries contain eight conserved metrics from positive acceptable-model ranks";
		totals[0]++;
		totals[2]+=rank;
		if(rank<=10){totals[4]++;}
		if(rank<=50){totals[5]++;}
		if(low){
			totals[1]++;
			totals[3]+=rank;
			if(rank<=10){totals[6]++;}
			if(rank<=50){totals[7]++;}
		}
	}

	/** One fixed512-bin recipe; existing frozen cohorts and unit-score ranker.
	 * Background means no fully accepted model in the COMPLETE saved oracle,
	 * not merely no target-parent link. Blacklists are exploratory on saved data. */
	private void study(final Row row, final NcrnaDdlMembers members, final String[] headers,
			final byte[][] consensuses, final NcrnaDdlSavedData data){
		final boolean best=c.study.equals("bestblacklist");
		require(row.shortK==(best ? 8 : 7) && row.longK==(best ? 10 : 9) && row.maskedLong && row.memberBins==(c.study.equals("variance1024") ? 1024 : 512)
			&& row.shortBins==0 && row.consensusBins==0 && row.bonus==100
			&& row.consensusMode.equals("POOL_BOTH") && row.terms==1 && !row.control,
			"Study preserves its explicit7/masked9 or8/masked10 pooled100 recipe;only the variance control changes512bins to1024");
		require(!best || oldTags, "The paired8/masked10 blacklist study preserves the accepted XOR namespace");
		final byte[][] names=new byte[headers.length][];
		for(int i=0; i<names.length; i++){names[i]=headers[i].getBytes(java.nio.charset.StandardCharsets.UTF_8);}
		final NcrnaDdlCohorts cohorts=NcrnaDdlCohorts.read(c.identities, data, names, consensuses);
		familyGc=familyGc(members);
		if(partitionStudy(c)){
			partitionStudy(row, members, headers, consensuses, data, cohorts);
			return;
		}
		try(Output output=new Output(c.output+".study.tsv")){
			output.writer.println(seededStudy(c) ? STUDY_HEADER+"\tddl_seed" : STUDY_HEADER);
			if(seededStudy(c)){
				seedStudy(row, members, consensuses, data, cohorts, output);
				return;
			}
			if(c.study.equals("baseline")){
				studyRow(row, oldTags ? "OLD_XOR_BASELINE" : "LEADING_BIT_BASELINE", 0, 0, members, consensuses, data, cohorts, output);
				return;
			}
			data.readNativeCounts(headers);
			backgroundSides(data, cohorts);
			try(Output lists=new Output(c.output+".blacklist.tsv")){
				lists.writer.println("definition\tword_form\trank\traw_word\ttagged_word\tfrequency\tforward_sequence\tfit_side");
				final int[] shortSizes=best ? new int[]{64,256} : SHORT_SIZES;
				final int[] longSizes=best ? new int[]{1024,4096} : LONG_SIZES;
				for(int definition=best ? 2 : 0; definition<(best ? 3 : 4); definition++){
					final String name=DEFINITIONS[definition];
					for(fitSide=definition==0 ? -1 : 0; fitSide<(definition==0 ? 0 : 2); fitSide++){
					final long[] shortWords=frequencies(row.shortK, false, definition, members, data);
					final long[] longWords=frequencies(row.longK, true, definition, members, data);
					writeBlacklist(name, "SHORT"+row.shortK, row.shortK, shortWords, shortSizes[shortSizes.length-1], true, lists);
					writeBlacklist(name, "MASKED"+row.longK, row.longK, longWords, longSizes[longSizes.length-1], false, lists);
					for(int a=best ? 0 : -1; a<shortSizes.length; a++){
						for(int b=best ? 0 : -1; b<longSizes.length; b++){
							if(a<0 && b<0){continue;}
							final int shortSize=a<0 ? 0 : shortSizes[a], longSize=b<0 ? 0 : longSizes[b];
							blacklist=new LongHashSet(shortSize+longSize+1);
							addBlacklist(shortWords, shortSize, row.shortK, true);
							addBlacklist(longWords, longSize, row.longK, false);
							shortGc=prefixGc(shortWords, shortSize, row.shortK, false);
							longGc=prefixGc(longWords, longSize, row.longK, true);
							require(blacklist.size()==shortSize+longSize, "The selected length tags keep the two raw blacklist forms disjoint");
							studyRow(row, name, shortSize, longSize, members, consensuses, data, cohorts, output);
						}
					}
					}
				}
			}finally{blacklist=null;}
		}
	}

	/** Reuses native DDL's literal-word LongHashSet membership mechanism.
	 * This study rejects BEFORE hashing; native hashAndStore checks AFTER hashing.
	 * Its FASTA loader canonicalizes a single untagged k, so it cannot load these
	 * forward, mixed-length, middle-masked words without changing their meaning. */
	private void addBlacklist(final long[] words, final int size, final int k, final boolean shortWord){
		require(size>=0 && size<=words.length, "Requested blacklist prefix exceeds observed positive-frequency words: k="+k+", requested="+size+", available="+words.length);
		for(int rank=0; rank<size; rank++){blacklist.add(rawWord(words, rank)|wordTag(k, shortWord, oldTags));}
	}

	/** Frequency descending, then forward raw-word ascending; primitive sort only.
	 * Presence counts at most once per model across ALL that model's members. */
	private long[] frequencies(final int k, final boolean masked, final int definition,
			final NcrnaDdlMembers members, final NcrnaDdlSavedData data){
		assert((masked ? k==9 || k==10 || k==11 : k==7 || k==8) && definition>=0 && definition<4)
			: "Only released short7/8 or masked9/10/11 forms and declared library/background corpora are counted";
		if(definition==0){
			final byte[][][] sequences=new byte[c.models][][];
			for(int model=0; model<c.models; model++){sequences[model]=members.members(model);}
			return memberPresenceFrequencies(k, masked, sequences, excludedMembers);
		}
		final int[] counts=new int[1<<(2*k)];
		wordFrequencies=counts;
		try{
			int backgrounds=0;
			for(final Window window:data.windows){
				require(window.nativeModels==c.models, "Background membership requires the complete original cached oracle");
				if(inputSides[window.input.index]!=fitSide){continue;}
				if(window.acceptedModels==0){
					if(definition==1 || definition==3){words(window.bases, k, masked, 0); backgrounds++;}
				}else if(definition==2 || definition==3){
					// Exclude the envelope of all accepted spans,not just target links.
					paddingBounds(window.start, window.stop, window.strand, window.acceptedStart, window.acceptedStop, padding);
					final int from=padding[0], to=padding[1];
					if(from>=k){words(Arrays.copyOfRange(window.bases, 0, from), k, masked, 0); backgrounds++;}
					if(window.length-to>=k){words(Arrays.copyOfRange(window.bases, to, window.length), k, masked, 0); backgrounds++;}
				}
			}
			require(backgrounds>0, "Background blacklist cannot be learned from an empty no-acceptable-model population");
			System.err.println("NATIVE_DDL_BACKGROUND_SEGMENTS\t"+DEFINITIONS[definition]+"\tfit="+fitSide+"\tn="+backgrounds+"\tk="+k);
		}finally{wordFrequencies=null;}
		return orderedFrequencies(counts);
	}

	/** Actual positive-member definition: native forward words,one presence per
	 * model,with all held-out/exact-duplicate members removed before counting.
	 * Returns frequency/raw-word packed values in the same prefix order as backgrounds. */
	public static long[] memberPresenceFrequencies(final int k, final boolean masked,
			final byte[][][] models, final boolean[][] excluded){
		require((masked ? k>=9 && k<=11 : k==7 || k==8) && models!=null && models.length>0
			&& (excluded==null || excluded.length==models.length),
			"Member-presence lists need a released word form and matching original model/exclusion rosters");
		final MemberPresence counter=new MemberPresence(k, masked);
		for(int model=0; model<models.length; model++){
			require(models[model]!=null && (excluded==null || excluded[model]!=null && excluded[model].length==models[model].length),
				"Member exclusions must bind every original model/member slot before any word is counted");
			counter.ordinal=model+1;
			for(int member=0; member<models[model].length; member++){
				if(excluded!=null && excluded[model][member]){continue;}
				require(models[model][member]!=null, "A retained voting member must have its actual sequence bytes");
				CardinalityTracker.forEachKmer(models[model][member], k, k, false, true, counter);
			}
		}
		return orderedFrequencies(counter.counts);
	}
	private static long[] orderedFrequencies(final int[] counts){
		assert(counts!=null && counts.length>0) : "Prefix sorting requires the actual finite raw-word domain";
		int positive=0;
		for(final int count:counts){if(count>0){positive++;}}
		final long[] ordered=new long[positive];
		int next=0;
		for(int raw=0; raw<counts.length; raw++){
			if(counts[raw]>0){ordered[next++]=((long)counts[raw]<<32)|(0x7fffffffL-raw);}
		}
		Arrays.sort(ordered);
		return ordered;
	}
	private static long rawWord(final long[] words, final int rank){
		assert(rank>=0 && rank<words.length) : "Blacklist prefixes retain only observed positive-frequency words";
		return 0x7fffffffL-(words[words.length-1-rank]&0xffffffffL);
	}
	private void writeBlacklist(final String name, final String form, final int k,
			final long[] words, final int size, final boolean shortWord, final Output out){
		require(size<=words.length, "Requested maximum blacklist must have enough observed words: "+form);
		for(int rank=0; rank<size; rank++){
			final long raw=rawWord(words, rank);
			line.clear().append(name).tab().append(form).tab().append(rank+1).tab().append(raw)
				.tab().append(raw|wordTag(k, shortWord, oldTags)).tab().append(words[words.length-1-rank]>>>32)
				.tab().append(DynamicDemiLog.unpackKmer(raw, k)).tab().append(fitSide).nl();
			out.writer.print(line);
		}
	}

	/** Shared DDL construction/comparison and shared saved-locus ranking.
	 * Scores include zero; ties keep ascending original model ordinals. */
	private NcrnaPartitionWeightSearch.Metrics studyRow(final Row row, final String definition, final int shortSize, final int longSize,
			final NcrnaDdlMembers members, final byte[][] consensuses, final NcrnaDdlSavedData data,
			final NcrnaDdlCohorts cohorts, final Output output){
		//TODO: Probable bug - older standalone study/index diagnostics call Vector directly and cannot run on Java8.
		// The replicated path uses Shared.SIMD and the allocation-free native scalar fallback.
		final long cpu=clock.getCurrentThreadCpuTime(), wall=System.nanoTime();
		final DynamicDemiLog[] references=buildReferences(row, members, consensuses)[0];
		final DynamicDemiLog query=createParts(row, 1)[0];
		final char[][] counts=new char[data.windows.size()][c.models];
		final int[][] parents=new int[counts.length][];
		final long[][][] allowed=new long[counts.length][][];
		final int[] strata=new int[data.parents.size()];
		for(final Parent parent:data.parents){strata[parent.index]=parent.stratum;}
		long nonzero=0;
		int judgedWindows=0;
		for(int w=0; w<counts.length; w++){
			final Window window=data.windows.get(w);
			parents[w]=new int[window.edges.size()];
			allowed[w]=new long[window.edges.size()][];
			for(int e=0; e<window.edges.size(); e++){
				parents[w][e]=window.edges.get(e).parent.index;
				allowed[w][e]=window.edges.get(e).uniqueModels;
			}
			if(fitSide>=0 && inputSides[window.input.index]!=1-fitSide){continue;}
			judgedWindows++;
			fill(query, window.bases, row, 0, true);
			for(int model=0; model<c.models; model++){
				simd.Vector.compareDDL(query.maxArray(), references[model].maxArray(), comparison);
				final int score=comparison[1];
				assert(score>=0 && score<=row.memberBins) : "A native match count cannot exceed its declared bin count";
				counts[w][model]=(char)score;
				if(score>0){nonzero++;}
			}
		}
		final NcrnaPartitionWeightSearch engine=new NcrnaPartitionWeightSearch(c.models, 1, counts,
			parents, allowed, strata, cohorts.sides, cohorts.bins, 1);
		final NcrnaPartitionWeightSearch.Metrics m=engine.evaluate(new double[]{1}).metrics[fitSide<0 ? 2 : 1-fitSide];
		if(seededStudy(c)){lastStudyRanks=engine.ranks(new double[]{1});}
		line.clear().append(definition).tab().append(shortSize).tab().append(longSize).tab().append(c.models)
			.tab().append(judgedWindows).tab().append(m.parents).tab().append(m.unavailable)
			.tab().append(m.availableRankSum).tab().append(m.within10).tab().append(m.within50)
			.tab().append(m.lowParents).tab().append(m.lowAvailableRankSum).tab().append(m.lowWithin10)
			.tab().append(m.lowWithin50).tab().append(nonzero).tab().append(nonzero/(double)judgedWindows, 9)
			.tab().append((clock.getCurrentThreadCpuTime()-cpu)/1e9, 9)
			.tab().append((System.nanoTime()-wall)/1e9, 9).tab().append(oldTags ? "OLD_XOR_HIGH_BITS_RAW_PREHASH" : "LEADING_BIT_RAW_PREHASH")
			.append("\tEXPLORATORY_SAVED_WINDOWS_NOT_CALLER")
			.tab().append(fitSide).tab().append(fitSide<0 ? -1 : 1-fitSide).tab();
		appendGc(line, shortSize==0 ? Double.NaN : shortGc);
		line.tab();
		appendGc(line, longSize==0 ? Double.NaN : longGc);
		line.tab().append(familyGc, 9);
		if(seededStudy(c) || c.study.equals("partitionseeds")){line.tab().append(sketchSeed);}
		line.nl();
		output.writer.print(line);
		System.err.print("NATIVE_DDL_STUDY_ROW\t"+line.toString());
		return m;
	}

	/** Only new partitioned rows are measured; accepted mixed rows remain their paired control. */
	private void partitionStudy(final Row mixed, final NcrnaDdlMembers members, final String[] headers,
			final byte[][] consensuses, final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts){
		require(oldTags, "The paired partition experiment preserves the accepted XOR pre-hash namespaces");
		final Row[] rows={new Row(mixed, 512, 1), new Row(mixed, 512, 2)};
		final boolean filtered=c.study.equals("partitionblacklist");
		if(filtered){data.readNativeCounts(headers); backgroundSides(data, cohorts);}
		try(Output output=new Output(c.output+".study.tsv");
				Output summary=filtered ? null : new Output(c.output+".seeds.summary.tsv")){
			output.writer.println(filtered ? STUDY_HEADER : STUDY_HEADER+"\tddl_seed");
			if(summary!=null){summary.writer.println("tag_form\tmetric\tseeds\tmean\tminimum\tmaximum");}
			for(final Row row:rows){
				final String name=row.partitionBits==1 ? "SINGLE_DDL_256_256" : "SINGLE_DDL_128_384";
				if(filtered){
					for(fitSide=0; fitSide<2; fitSide++){
						final long[] shortWords=frequencies(7, false, 2, members, data);
						final long[] longWords=frequencies(9, true, 2, members, data);
						blacklist=new LongHashSet(64+4096+1);
						addBlacklist(shortWords, 64, 7, true);
						addBlacklist(longWords, 4096, 9, false);
						shortGc=prefixGc(shortWords, 64, 7, false);
						longGc=prefixGc(longWords, 4096, 9, true);
						studyRow(row, name, 64, 4096, members, consensuses, data, cohorts, output);
					}
				}else{
					final long[] sum=new long[6], min=new long[6], max=new long[6];
					Arrays.fill(min, Long.MAX_VALUE);
					for(final long seed:STUDY_SEEDS){
						sketchSeed=seed;
						CardinalityTracker.defaultSeed=seed;
						final NcrnaPartitionWeightSearch.Metrics m=studyRow(row, name, 0, 0, members, consensuses, data, cohorts, output);
						final long[] values={m.availableRankSum,m.lowAvailableRankSum,m.within10,m.within50,m.lowWithin10,m.lowWithin50};
						for(int i=0; i<values.length; i++){
							sum[i]+=values[i];
							min[i]=Math.min(min[i], values[i]);
							max[i]=Math.max(max[i], values[i]);
						}
					}
					for(int i=0; i<sum.length; i++){
						summary.writer.print(line.clear().append(name).tab().append(SEED_METRICS[i]).tab().append(STUDY_SEEDS.length)
							.tab().append(sum[i]/(double)STUDY_SEEDS.length, 9).tab().append(min[i]).tab().append(max[i]).nl());
					}
				}
			}
		}finally{blacklist=null; fitSide=-1;}
	}
	private static boolean partitionStudy(final Config c){return "partitionseeds".equals(c.study) || "partitionblacklist".equals(c.study);}

	/** Paired seeds isolate hash-choice variability from the requested tag change.
	 * Every row uses the identical fixed recipe,data,unit ranker and seed family. */
	private void seedStudy(final Row row, final NcrnaDdlMembers members, final byte[][] consensuses,
			final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts, final Output output){
		final boolean variance=c.study.equals("variance1024");
		final long[] seeds=variance ? Arrays.copyOf(STUDY_SEEDS, 4) : STUDY_SEEDS;
		try(Output summary=new Output(c.output+".seeds.summary.tsv"); Output loci=new Output(c.output+".loci.tsv");
				Output ranges=new Output(c.output+".locus_ranges.tsv"); Output tails=new Output(c.output+".tails.tsv")){
			summary.writer.println("tag_form\tmetric\tseeds\tmean\tminimum\tmaximum");
			loci.writer.println("tag_form\tbins\tddl_seed\tpanel\tinput_slot\tinput_id\tparent\tassembly_group\tlow_identity\trank");
			line.clear().append("tag_form\tbins\tpanel\tinput_id\tparent\tlow_identity\tminimum_rank\tmaximum_rank\trank_range");
			for(final long seed:seeds){line.append("\tseed_").append(seed);}
			ranges.writer.print(line.nl());
			tails.writer.println("tag_form\tbins\tminimum_sum_seed\tmaximum_sum_seed\tminimum_rank_sum\tmaximum_rank_sum\tnet_rank_spread\tpositive_delta_loci\tnegative_delta_loci\tpositive_rank_units\tnegative_rank_units\ttop5_positive_units\tloci_for_80pct_positive_units");
			for(int tag=0; tag<(variance ? 1 : 2); tag++){
				oldTags=tag==0;
				final String name=oldTags ? "OLD_XOR_HIGH_BITS" : "LEADING_BIT";
				final long[] sum=new long[6], min=new long[6], max=new long[6];
				final int[][] ranks=new int[seeds.length][];
				final long[] totals=new long[seeds.length];
				Arrays.fill(min, Long.MAX_VALUE);
				for(int s=0; s<seeds.length; s++){
					final long seed=seeds[s];
					sketchSeed=seed;
					CardinalityTracker.defaultSeed=seed;
					final NcrnaPartitionWeightSearch.Metrics m=studyRow(row, name, 0, 0, members, consensuses, data, cohorts, output);
					ranks[s]=lastStudyRanks;
					totals[s]=m.availableRankSum;
					long locusSum=0;
					for(final Parent parent:data.parents){
						if(parent.stratum!=0){continue;}
						final int rank=ranks[s][parent.index];
						if(rank>0){locusSum+=rank;}
						loci.writer.print(line.clear().append(name).tab().append(row.memberBins).tab().append(seed)
							.tab().append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab().append(parent.input.idBytes)
							.tab().append(parent.idBytes).tab().append(parent.input.group).tab().append(cohorts.bins[parent.index]==0)
							.tab().append(rank).nl());
					}
					require(locusSum==m.availableRankSum, "Per-locus output must recombine to the exact reported not-seen rank sum");
					final long[] v={m.availableRankSum,m.lowAvailableRankSum,m.within10,m.within50,m.lowWithin10,m.lowWithin50};
						for(int i=0; i<v.length; i++){
							sum[i]+=v[i];
							min[i]=Math.min(min[i], v[i]);
							max[i]=Math.max(max[i], v[i]);
						}
				}
				for(int i=0; i<sum.length; i++){
					summary.writer.print(line.clear().append(name).tab().append(SEED_METRICS[i]).tab().append(seeds.length)
						.tab().append(sum[i]/(double)seeds.length, 9).tab().append(min[i]).tab().append(max[i]).nl());
				}
				seedTails(name, row.memberBins, seeds, ranks, totals, data, cohorts, ranges, tails);
			}
		}
	}

	private void seedTails(final String tag, final int bins, final long[] seeds, final int[][] ranks,
			final long[] totals, final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts,
			final Output ranges, final Output tails){
		int minimum=0, maximum=0;
		for(int s=1; s<seeds.length; s++){
			if(totals[s]<totals[minimum]){minimum=s;}
			if(totals[s]>totals[maximum]){maximum=s;}
		}
		final LongList ordered=new LongList(), deltas=new LongList();
		long positive=0, negative=0;
		int negativeLoci=0;
		for(final Parent parent:data.parents){
			if(parent.stratum!=0){continue;}
			int low=ranks[0][parent.index], high=low;
			for(final int[] values:ranks){
				final int rank=values[parent.index];
				require((rank<0)==(low<0), "Hash choice cannot change the frozen acceptable-model population");
				low=Math.min(low, rank);
				high=Math.max(high, rank);
			}
			ordered.add(((long)(high-low)<<32)|(0x7fffffffL-parent.index));
			final int delta=ranks[maximum][parent.index]-ranks[minimum][parent.index];
			if(delta>0){positive+=delta; deltas.add(delta);}
			else if(delta<0){negative-=delta; negativeLoci++;}
		}
		final long[] sorted=ordered.toArray(), differences=deltas.toArray();
		Arrays.sort(sorted);
		Arrays.sort(differences);
		for(int i=sorted.length-1; i>=0; i--){
			final int index=(int)(0x7fffffffL-(sorted[i]&0xffffffffL));
			final Parent parent=data.parents.get(index);
			int low=ranks[0][index], high=low;
			for(final int[] values:ranks){low=Math.min(low, values[index]); high=Math.max(high, values[index]);}
			line.clear().append(tag).tab().append(bins).tab().append(parent.input.panel.nameBytes).tab().append(parent.input.idBytes)
				.tab().append(parent.idBytes).tab().append(cohorts.bins[index]==0).tab().append(low).tab().append(high)
				.tab().append(high-low);
			for(final int[] values:ranks){line.tab().append(values[index]);}
			ranges.writer.print(line.nl());
		}
		long top5=0, cumulative=0;
		int for80=0;
		for(int i=differences.length-1; i>=0; i--){
			if(i>=differences.length-5){top5+=differences[i];}
			if(cumulative<(4*positive+4)/5){cumulative+=differences[i]; for80++;}
		}
		require(positive-negative==totals[maximum]-totals[minimum], "Tail attribution must exactly conserve the aggregate seed rank-sum spread");
		tails.writer.print(line.clear().append(tag).tab().append(bins).tab().append(seeds[minimum]).tab().append(seeds[maximum])
			.tab().append(totals[minimum]).tab().append(totals[maximum]).tab().append(positive-negative)
			.tab().append(differences.length).tab().append(negativeLoci).tab().append(positive).tab().append(negative)
			.tab().append(top5).tab().append(for80).nl());
	}
	private static boolean seededStudy(final Config c){return "seeds".equals(c.study) || "variance1024".equals(c.study);}

	/** Preserve the frozen target-group halves. Previously unassigned background
	 * groups get their own deterministic split; none can enter both directions. */
	private void backgroundSides(final NcrnaDdlSavedData data, final NcrnaDdlCohorts cohorts){
		final IntHashMap sides=new IntHashMap(c.groups);
		for(final Parent parent:data.parents){
			if(parent.stratum==0){
				final int side=cohorts.sides[parent.index];
				require(sides.get(parent.input.group)<0 || sides.get(parent.input.group)==side,
					"All parents of one assembly group retain the same frozen half");
				sides.put(parent.input.group, side);
			}
		}
		final IntList extra=new IntList();
		for(final NcrnaDdlSavedData.Selected input:data.selected){if(sides.get(input.group)<0){extra.add(input.group);}}
		extra.sort();
		final int[] groups=extra.toArray(), split=NcrnaDdlCohorts.split(groups);
		for(int i=0; i<groups.length; i++){sides.put(groups[i], split[i]);}
		inputSides=new int[data.selected.size()];
		try(Output output=new Output(c.output+".background_groups.tsv")){
			output.writer.println("assembly_group\tside\tassignment");
			for(final NcrnaDdlSavedData.Selected input:data.selected){
				final int side=sides.get(input.group);
				require(side==0 || side==1, "Every selected background group needs exactly one deterministic half");
				inputSides[input.index]=side;
				output.writer.print(line.clear().append(input.group).tab().append(side).tab()
					.append(Arrays.binarySearch(groups, input.group)<0 ? "FROZEN_TARGET_HALF" : "ADDITIONAL_BACKGROUND_SPLIT").nl());
			}
		}
	}
	private static double familyGc(final NcrnaDdlMembers members){
		long gc=0, defined=0;
		for(int m=0; m<members.models(); m++){
			for(final byte[] bases:members.members(m)){
				for(final byte raw:bases){
					final int b=raw&~32;
					if(b=='A' || b=='C' || b=='T' || b=='G'){defined++; if(b=='C' || b=='G'){gc++;}}
				}
			}
		}
		require(defined>0, "Family composition needs defined bases from the actual construction members");
		return gc/(double)defined;
	}
	public static double prefixGc(final long[] words, final int size, final int k, final boolean masked){
		assert(words!=null && size>=0 && size<=words.length && k>(masked ? 1 : 0) && k<=15)
			: "GC reads an observed raw-word prefix and divides by only its contributing nucleotide positions";
		if(size==0){return Double.NaN;}
		long gc=0;
		for(int rank=0; rank<size; rank++){gc+=Long.bitCount(rawWord(words, rank)&0x5555555555555555L);}
		return gc/(double)((long)size*(k-(masked ? 1 : 0)));
	}
	private static void appendGc(final ByteBuilder b, final double value){
		if(Double.isNaN(value)){b.append("NA");}
		else{
			assert(value>=0 && value<=1) : "GC fractions count defined contributing bases only";
			b.append(value, 9);
		}
	}

	/** Converts an original-forward one-based accepted envelope to the half-open
	 * occupied interval in a zero-based strand-oriented saved window. Only the
	 * prefix/suffix outside these bounds may enter the padding background. */
	public static void paddingBounds(final int start, final int stop, final int strand,
			final int acceptedStart, final int acceptedStop, final int[] bounds){
		require(start>=0 && stop>=start && (strand==0 || strand==1) && acceptedStart>0
			&& acceptedStop>=acceptedStart && acceptedStop>start && acceptedStart<=stop+1
			&& bounds!=null && bounds.length>=2,
			"Padding requires an overlapping original-forward accepted envelope and two caller-owned bounds");
		final int length=stop-start+1;
		final int left=strand==0 ? acceptedStart-1-start : stop-acceptedStop+1;
		final int right=strand==0 ? acceptedStop-start : stop-acceptedStart+2;
		bounds[0]=Tools.mid(0, left, length);
		bounds[1]=Tools.mid(0, right, length);
		assert(bounds[0]<=bounds[1]) : "Strand conversion must preserve the accepted envelope before extracting either padding flank";
	}

	private void writeCounts(final Row row, final DynamicDemiLog[][] references, final NcrnaDdlSavedData data, final Output output){
		final DynamicDemiLog[] query=createParts(row, row.terms);
		for(int w=0; w<data.windows.size(); w++){
			final Window window=data.windows.get(w);
			for(int part=0; part<row.terms; part++){
				fill(query[part], window.bases, row, part, part<row.memberTerms || row.consensusMode.equals("BOTH"));
			}
			line.clear().append(row.id).tab().append(w).tab().append(window.input.panel.nameBytes)
				.tab().append(window.input.slot).tab().append(window.event);
			for(int model=0; model<c.models; model++){
				for(int part=0; part<4; part++){
					int count=0;
					if(part<row.terms){
						require(query[part].buckets==references[part][model].buckets, "Scoring compares matching native DDL parts");
						simd.Vector.compareDDL(query[part].maxArray(), references[part][model].maxArray(), comparison);
						count=comparison[1];
						assert(count>=0 && count<=row.width(part)) : "Shared nonempty equal-bin counts are bounded by this DDL's width";
					}
					line.tab().append(count);
				}
			}
			output.writer.print(line.nl());
		}
	}

	private void writeMetadata(final Row row, final Output output){
		line.clear().append(2).tab().append(row.id).append("\tRF00174\t").append(row.shortK).tab().append(row.longK)
			.tab().append(row.maskedLong ? "MIDDLE_RIGHT" : "NONE").tab().append(row.shortBins>0 ? "PARTITIONED" : "MIXED")
			.tab().append(row.memberBins).tab().append(row.shortBins).tab().append(row.shortBins>0 ? row.memberBins-row.shortBins : 0)
			.tab().append(row.consensusMode).tab().append(row.consensusBins).tab().append(row.bonus)
			.append("\t4\t12\tWRAP\t").append(row.termNames).tab().append(row.shortBins>0 ? 1 : 0).tab().append(c.models)
			.tab().append(2L*c.models*row.memberBins).tab().append(2L*c.models*row.consensusBins)
			.append("\t0\t0\t0\t0\t0\t0")
			.tab().append(consensusCpu/1e9, 9).tab().append(consensusWall/1e9, 9)
			.tab().append(row.control ? "LEGACY_PARTITIONS" : "NATIVE_DDL").append("\tQUERY_DDL\n");
		output.writer.print(line);
	}

	/** Checks every saved window/model, including zero matches, before measuring.
	 * Native inputs are built once; only query work enters the timed sweeps. */
	private void checkIndex(final Row row, final DynamicDemiLog[] references, final NcrnaDdlSavedData data){
		require(row.terms==1 && references.length==c.models && !data.windows.isEmpty(),
			"Saved parity requires one standard DDL per original reference and the complete nonempty window roster");
		final DDLIndexSmall index=new DDLIndexSmall(references);
		final DynamicDemiLog[] queries=new DynamicDemiLog[data.windows.size()];
		final int[] counts=new int[c.models];
		long expected=0;
		for(int w=0; w<queries.length; w++){
			final DynamicDemiLog query=DynamicDemiLog.create(row.memberBins, row.longK, sketchSeed, 0, false, false);
			fill(query, data.windows.get(w).bases, row, 0, true);
			queries[w]=query;
			index.query(query, counts);
			for(int model=0; model<references.length; model++){
				simd.Vector.compareDDL(query.maxArray(), references[model].maxArray(), comparison);
				require(counts[model]==comparison[1], "Index/direct score mismatch: window="+w+", model="+model
					+", index="+counts[model]+", direct="+comparison[1]);
				expected+=counts[model];
			}
		}
		final long pairs=Math.multiplyExact((long)queries.length, references.length);
		System.err.println("NATIVE_DDL_INDEX_ALL_PAIRS_PASS\t"+queries.length+"\t"+references.length+"\t"+pairs);
		try(Output output=new Output(c.output+".index.tsv")){
			output.writer.println(INDEX_HEADER);
			if(c.indexSeconds==0){
				indexRow(output, row.id, "PARITY_ONLY", queries.length, pairs, 0, 0, expected, index);
				return;
			}
			require(Shared.SIMD && simd.Vector.simd256 && c.indexSeconds>=10,
				"Measured direct comparison must dispatch Vector.compareDDL to SIMDLogLog, with at least ten CPU seconds per observation");
			System.err.println("DIRECT_TIMING_METHOD\tsimd.Vector.compareDDL(char[],char[],int[]) -> SIMDLogLog.compareDetailed");
			final long minimum=1_000_000_000L*c.indexSeconds;
			for(int round=0; round<3; round++){
				for(int order=0; order<2; order++){
					final boolean indexed=(order==(round&1));
					long calls=0, cpu=0, checksum=0;
					do{
						final long start=clock.getCurrentThreadCpuTime();
						for(final DynamicDemiLog query:queries){
							if(indexed){
								index.query(query, counts);
								for(final int count:counts){checksum+=count;}
							}else{
								final char[] keys=query.maxArray();
								for(final DynamicDemiLog reference:references){
									simd.Vector.compareDDL(keys, reference.maxArray(), comparison);
									checksum+=comparison[1];
								}
							}
						}
						cpu+=clock.getCurrentThreadCpuTime()-start;
						calls+=queries.length;
					}while(cpu<minimum);
					require(checksum==Math.multiplyExact(expected, calls/queries.length),
						"Every timed path must consume all checked model scores on each complete window sweep");
					indexRow(output, row.id, (indexed ? "SMALL_INDEX_OBJECT_QUERY_R" : "DIRECT_SIMD_VECTOR_COMPARE_DDL_R")+(round+1),
						queries.length, pairs, calls, cpu, checksum, index);
				}
			}
		}
		System.err.println("NATIVE_DDL_INDEX_PARITY_AND_QUERY_TIMING_COMPLETE");
	}

	/** One row per measured query path; checksum consumes every returned score.
	 * Finite parity-only fixtures report no timing instead of a startup estimate. */
	private void indexRow(final Output output, final String id, final String method, final int windows,
			final long pairs, final long calls, final long cpu, final long checksum, final DDLIndexSmall index){
		assert(windows>0 && pairs==(long)windows*c.models && calls>=0 && cpu>=0)
			: "A query means one complete original window against every original model";
		line.clear().append(id).tab().append(c.models).tab().append(windows).tab().append(pairs).tab()
			.append(method).tab().append(calls).tab();
		if(calls==0){line.append("NA\tNA");}
		else{line.append(cpu/1e9, 9).tab().append(cpu/(1000.0*calls), 9);}
		line.tab().append(checksum).tab().append(index.postingPayloadBytes()).tab().append(index.pointerSlots())
			.tab().append(Shared.SIMD).tab().append(calls==0 ? "FINITE_PARITY_NOT_TIMING"
				: "QUERY_AND_RESULT_CHECKSUM_ONLY_CONSTRUCTION_HASHING_IO_EXCLUDED").nl();
		output.writer.print(line);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Metadata           ----------------*/
	/*--------------------------------------------------------------*/

	/** One reusable primitive callback;no per-word or per-member objects. */
	private static final class MemberPresence implements CardinalityTracker.KmerConsumer{
		MemberPresence(final int k, final boolean masked){
			counts=new int[1<<(2*k)];
			seen=new int[counts.length];
			mask=AbstractKmerTable.makeMiddleMask(k, false, masked ? 1 : 0);
		}
		@Override
		public void accept(final long word, final int k){
			final int raw=(int)(word&mask);
			assert(ordinal>0 && raw>=0 && raw<counts.length) : "Native member words must bind one actual model and its packed word domain";
			if(seen[raw]==ordinal){return;}
			seen[raw]=ordinal;
			counts[raw]++;
		}
		final int[] counts, seen;
		final long mask;
		int ordinal;
	}

	public static final class Config{
		public String recipes, panels, selection, references, members, ledger, assignments, census, summary, output, checkIndex;
		public String study, identities, mergeStudy, mergeSeeds, seedBaseline;
		public String sourceGroups, genomes, gcFits, gcAllFit, gcAllLimits;
		public long foldSeed=20261114;
		public boolean fixedCheck, indexReference, gcAdditions;
		public int folds=5, longPrefix=4096;
		public int[] shortPrefixes={64,64,256};
		public long seed=12345;
		public boolean oldTags=true;
		public int models=508, groups=3153;
		public int indexSeconds=10;
	}

	private static final class Row{
		Row(final Row source, final int width){
			this(source, width, 0);
		}
		Row(final Row source, final int width, final int partitionBits_){
			this(source, width, partitionBits_, source.longK);
		}
		Row(final Row source, final int width, final int partitionBits_, final int longK_){
			require(source.memberTerms==1 && source.consensusBins==0 && source.terms==1
				&& powerOfTwo(width) && (partitionBits_==0 || (width==512 || width==1024) && (partitionBits_==1 || partitionBits_==2)),
				"Selector partition controls use one native512/1024-bin DDL with one or two partition bits");
			id=source.id+"_bins"+width;
			shortK=source.shortK;
			require(longK_==source.longK || source.shortK==8 && source.longK==10 && longK_==11,
				"The follow-up changes only the declared8/masked10 pair to8/masked11");
			longK=longK_;
			maskedLong=source.maskedLong;
			memberBins=width;
			shortBins=0;
			partitionBits=partitionBits_;
			consensusBins=0;
			consensusMode=source.consensusMode;
			bonus=source.bonus;
			memberTerms=terms=1;
			termNames=source.termNames;
			control=false;
		}
		Row(final LineParser1 p){
			partitionBits=0;
			id=p.parseString(0);
			shortK=p.parseInt(1);
			longK=p.parseInt(2);
			maskedLong=p.termEquals("MIDDLE_RIGHT", 3);
			final boolean partitioned=p.termEquals("PARTITIONED", 4);
			memberBins=p.parseInt(5);
			shortBins=p.parseInt(6);
			consensusMode=p.parseString(7);
			consensusBins=p.parseInt(8);
			bonus=p.parseInt(9);
			control=p.termEquals("true", 13);
			require(!id.isEmpty() && shortK>0 && shortK<longK && longK<=15
				&& (maskedLong || p.termEquals("NONE", 3)) && (partitioned || p.termEquals("MIXED", 4)),
				"Each native recipe needs an ID, actual word lengths, mask and member mode");
			require(memberBins>0 && memberBins<=32768 && (partitioned
				? shortBins>0 && shortBins<memberBins && powerOfTwo(shortBins) && powerOfTwo(memberBins-shortBins)
				: shortBins==0 && powerOfTwo(memberBins)), "Every actual member part must be a standard power-of-two DDL");
			final boolean separate=consensusMode.equals("BOTH") || consensusMode.equals("LONG_ONLY");
			final boolean pooled=consensusMode.equals("POOL_BOTH") || consensusMode.equals("POOL_LONG");
			require(separate ? powerOfTwo(consensusBins) && consensusBins<=32768 && bonus==0
				: pooled ? consensusBins==0 && bonus>0 : consensusMode.equals("NONE") && consensusBins==0 && bonus==0,
				"Separate consensus uses a standard DDL; pooled bonus and member-only modes own no separate part");
			require(p.parseInt(10)==4 && p.parseInt(11)==12 && (control || p.termEquals("false", 13))
				&& p.termEquals(control ? "LEGACY_PARTITIONS" : "NATIVE_DDL", 12) && p.termEquals("QUERY_DDL", 14),
				"Native objects use declared4+12 keys, an explicit control and query-DDL only; custom/raw paths are withdrawn");
			if(control){
				require(shortK==7 && longK==10 && maskedLong && partitioned && memberBins==384 && shortBins==128
					&& consensusMode.equals("BOTH") && consensusBins==128 && bonus==0,
					"The legacy standard control retains exact128+256+128 geometry, word hashes and smaller-key vote ties");
			}
			memberTerms=partitioned ? 2 : 1;
			terms=memberTerms+(separate ? 1 : 0);
			termNames=partitioned ? separate ? "member_short,member_long,consensus_total" : "member_short,member_long"
				: separate ? "member_total,consensus_total" : "member_total";
		}
		int width(final int part){return part>=memberTerms ? consensusBins : shortBins==0 ? memberBins : part==0 ? shortBins : memberBins-shortBins;}
		int wordLength(final int part){return shortBins>0 && part==0 ? shortK : longK;}
		final String id, consensusMode, termNames;
		final int shortK, longK, memberBins, shortBins, consensusBins, bonus, memberTerms, terms, partitionBits;
		final boolean maskedLong, control;
	}

	private static boolean powerOfTwo(final int value){return value>0 && (value&(value-1))==0;}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Config c;
	private final ThreadMXBean clock;
	private final int[][][] rankHistograms;
	private DynamicDemiLog active;
	private LongHashSet blacklist;
	private int[] wordFrequencies;
	private int activePartitionBits;
	private int activeShortK;
	private boolean[][] excludedMembers;
	private int replicateBonus=-1;
	private double replicateScale;
	private double[] retention;
	private long[] gcWindowMetrics;
	private double[][][] allGcStatistics;
	private long[][] allGcFallback;
	private double[] allGcLosslessZ;
	private long[] allGcWindowCounts;
	private int[] gcPairOffsets, gcPairModels;
	private double[] gcPairDistances;
	private int[] sizeThirds;
	private int fitSide=-1;
	private int[] inputSides;
	private int[] lastStudyRanks;
	private double shortGc=Double.NaN, longGc=Double.NaN, familyGc;
	private long activeMask, activeTag, consensusCpu, consensusWall;
	private long sketchSeed;
	private boolean oldTags;
	private final int[] comparison=new int[4];
	private final int[] padding=new int[2];
	private final ByteBuilder line=new ByteBuilder(8192);
	private static final long SHORT_TAG=1L<<62, LONG_TAG=1L<<63;
	private static final String WINNER="native_s7_l9_mixed512_pooled_both_w100";
	private static final String BEST_PAIRED="native_s8_l10_mixed512_pooled_both_w100";
	private static final int[] SHORT_SIZES={16,64,256,1024}, LONG_SIZES={64,256,1024,4096};
	private static final String[] DEFINITIONS={"MODEL_PRESENCE", "FALSE_SEEDED_WINDOWS", "TRUE_WINDOW_PADDING", "POOLED_BACKGROUND"};
	private static final String[] PARTITION_NAMES={"mixed","half","quarter"};
	private static final String[] GC_SOURCES={"MEMBERS", "ACCEPTED_PADDED_WINDOWS", "MODEL_ACCEPTED_SPANS"};
	public static final String LEGACY_REPLICATE_HEADER="assignment_seed\tddl_seed\tfold\tmethod\tshort_k\tlong_k\tlayout\tbins\tshort_blacklist\tlong_blacklist\tfixed_bonus\tscaled_bonus\tparents\tlow_parents\trank_sum\tlow_rank_sum\ttop10\ttop50\tlow_top10\tlow_top50\tempty_fit_models\tconsensus_retained_small\tconsensus_retained_medium\tconsensus_retained_large\tcpu_seconds\telapsed_seconds\tinterpretation";
	public static final String REPLICATE_HEADER=LEGACY_REPLICATE_HEADER+"\tassembly_contribution";
	public static final String FOLLOWUP_REPLICATE_HEADER=REPLICATE_HEADER+"\trank_histogram\tlow_rank_histogram";
	public static final String GC_REPLICATE_HEADER=FOLLOWUP_REPLICATE_HEADER
		+"\tacceptable_model_window_pairs\tdropped_acceptable_model_window_pairs\twindows_with_acceptable_models\twindows_all_acceptable_filtered\tmissing_query_gc_windows"
		+"\twindow_stat_fallback_model_fold_cases\tfamily_member_fallback_model_fold_cases\tzero_sd_model_fold_cases\tmissing_effective_gc_model_fold_cases";
	public static final String[] GC_SETTINGS={"sort_w0","sort_w0.5","sort_w1","sort_w2","sort_w4","sort_w8", "filter_z2","filter_z3","filter_z4", "filter_m0","filter_m0.02","filter_m0.05"};
	public static final String[] GC_ALL_SETTINGS={"sort_w0","sort_w0.5","sort_w1","sort_w2","sort_w4","sort_w8",
		"filter_z2","filter_z3","filter_z4","filter_z6","filter_z8","filter_lossless_z",
		"filter_m0","filter_m0.01","filter_m0.02","filter_m0.05"};
	public static final String[] GC_EXTRA_SETTINGS={"range_m0_sort_w0.5","range_m0_sort_w1","range_m0_sort_w2",
		"range_m0.01_sort_w0.5","range_m0.01_sort_w1","range_m0.01_sort_w2"};
	private static final String[] GC_WINDOW_STRATA={"NO_WINDOWS","LT10","10_TO29","GE30"};
	private static final double[] GC_LOOKUP_Z={4,4.25,5,8.58};
	public static final double[] GC_POWER_Z={1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,
		4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.58};
	private static final double[] GC_ALL_Z={2,3,4,6,8}, GC_ALL_MARGINS={0,0.01,0.02,0.05};
	public static final String[] GC_ARMS={"window","member"};
	public static final String GC_FIT_HEADER="fold\tmodel_ordinal\tgc_source\twindow_count\tmember_count\tmissing_window_gc\tmissing_member_gc\tmean_gc\tsd_gc\tmin_gc\tmax_gc\twindow_fallback\tfamily_member_fallback\twindow_member_offset\toffset_paired_models";
	public static final String GC_ALL_LIMIT_HEADER="gc_source\taccepted_model_window_pairs\tfinite_distances\tmissing_gc_pairs\tinfinite_zero_sd_pairs\tmaximum_finite_distance\tmaximum_distance\tlossless_z\tlossless_margin\tinterpretation";
	public static final String GC_POWER_HEADER="cohort\tfilter\tparameter\tmodels\twindows\tmean_models_removed\tmedian_models_removed\taccepted_pairs\tpairs_removed\taccepted_windows\twindows_lost\tloci\tloci_lost"
		+"\tbaseline_candidate_pairs\tcandidate_pairs_retained\tmean_models_retained\tmedian_models_retained\tpercent_candidates_removed\tfailed_locus_surviving_candidate_sum\tinterpretation";
	public static final double[] GC_WEIGHTS={0,0.5,1,2,4,8}, GC_MARGINS={0,0.02,0.05};
	public static final double[] GC_WIDTHS={1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.5,7,7.5,8};
	private static final int[] FOLLOWUP_BONUSES={60,80,100,120};
	public static final int[] GRID_SHORT={0,16,32,64,128,256,512,1024};
	public static final int[] GRID_LONG={0,256,512,1024,2048,4096,8192,16384,32768};
	public static final int[][] SIZE_EXTENSION={
		{2048,0},{2048,128},{2048,256},{2048,512},{2048,1024},
		{4096,0},{4096,128},{4096,256},{4096,512},{4096,1024},
		{8192,0},{8192,128},{8192,256},{8192,512},{8192,1024},
		{16384,0},{16384,128},{16384,256},{16384,512},{16384,1024},
		{512,128},{1024,128}};
	public static final int[][] SIZE_GRID=gridSizes(false), EXPANDED_SIZE_GRID=gridSizes(true);
	public static final String[] DEFINITION_NAMES={"members","false","padding","pooled"};
	public static final String[] CONTRIBUTION={"ALL", "CONTRIBUTOR", "NONCONTRIBUTOR"};
	private static final long[] STUDY_SEEDS={12345,1,7,42,123,1337,20261114,8675309};
	private static final String[] SEED_METRICS={"available_rank_sum","low_available_rank_sum","within10","within50","low_within10","low_within50"};
	private static final String[] SUMMARY_METRICS={"available_rank_sum","low_available_rank_sum","within10","within50","low_within10","low_within50","mean_nonzero_models_per_window"};
	private static final int[] SUMMARY_FIELDS={7,11,8,9,12,13,15};
	public static final String STUDY_HEADER="definition\tshort_blacklist\tlong_blacklist\tmodels\twindows\tnot_seen_parents\tunavailable\tavailable_rank_sum\twithin10\twithin50\tlow_parents\tlow_available_rank_sum\tlow_within10\tlow_within50\tnonzero_pairs\tmean_nonzero_models_per_window\trow_cpu_seconds\trow_elapsed_seconds\tword_encoding\tinterpretation\tfit_side\tjudge_side\tshort_prefix_gc\tlong_prefix_gc\tfamily_member_gc";
	public static final String MERGED_STUDY_HEADER="definition\tshort_blacklist\tlong_blacklist\tmodels\twindows\tnot_seen_parents\tunavailable\tavailable_rank_sum\twithin10\twithin50\tlow_parents\tlow_available_rank_sum\tlow_within10\tlow_within50\tnonzero_pairs\tmean_nonzero_models_per_window\tlist0_short_gc\tlist0_long_gc\tlist1_short_gc\tlist1_long_gc\tfamily_member_gc\tinterpretation";
	public static final String RECIPE_HEADER="layout_id\tshort_k\tlong_k\tlong_mask\tmember_mode\tmember_bins\tmember_short_bins\tconsensus_mode\tconsensus_bins\tconsensus_bonus\texponent_bits\tmantissa_bits\thash_namespace\tlegacy_control\tquery_form";
	public static final String INDEX_HEADER="layout_id\tmodels\twindows\tchecked_window_model_pairs\tmethod\tquery_calls\tquery_cpu_seconds\tcpu_us_per_window\tresult_checksum\tposting_payload_bytes\tpointer_slots\tsimd\tinterpretation";
}
