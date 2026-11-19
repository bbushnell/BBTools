package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import ddl.DDLIndexSmall;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaPartitionWeightSearch.Fit;
import prok.NcrnaPartitionWeightSearch.Metrics;
import prok.NcrnaPartitionWeightSearch.Trial;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Saved-count, assembly-group cross-fit of compact partition weights.
 * Uses the native template.A_SampleByteFile parse/validate/process lifecycle.
 * Loaded resources and original oracle bits are immutable; no model construction,
 * identity measurement, alignment or genome calling is performed.
 * @author Ganyu
 */
public final class OptimizeNcrnaPartitionWeights{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, OptimizeNcrnaPartitionWeights.class, false);
		final int threads=Shared.threads();
		try{
			Shared.setThreads(1);
			final Config c=parse(pp.args);
			require(Shared.threads()==1, "Partition fitting is a single-owner saved-count pass; use t=1");
			run(c);
			timer.stop();
			pp.outstream.println("Time: "+timer);
			pp.outstream.println("NCRNA_PARTITION_WEIGHT_SEARCH_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			Shared.setThreads(threads);
			Shared.closeStream(pp.outstream);
		}
	}

	private static Config parse(final String[] args){
		assert(args!=null) : "PreParser expands complete native flag=value arguments";
		final Config c=new Config();
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("panels")){c.panels=b;}
			else if(a.equals("selection")){c.selection=b;}
			else if(a.equals("ref")){c.references=b;}
			else if(a.equals("identities")){c.identities=b;}
			else if(a.equals("ranks")){c.ranks=b;}
			else if(a.equals("counts")){c.counts=b;}
			else if(a.equals("layouts")){c.layouts=b;}
			else if(a.equals("control")){c.control=b;}
			else if(a.equals("fitlayouts")){
				require(b!=null, "fitlayouts requires an explicit comma-separated ID list");
				c.fitLayouts=b;
			}
			else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
			else if(a.equals("randomdraws")){c.randomDraws=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown partition-weight parameter: "+arg);}
		}
		require(!parser.overwrite && !parser.append && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Use complete saved populations with fresh output, no sampling, overwrite or append");
		c.left=parser.in1;
		c.right=parser.in2;
		c.output=parser.out1;
		return c;
	}

	public static void run(final Config c){
		require(c!=null && c.panels!=null && c.selection!=null
			&& c.references!=null && c.identities!=null && c.ranks!=null && c.output!=null && c.models>0 && c.groups>0 && c.randomDraws>0,
			"Bind saved panels/selection, reference, frozen identities/ranks and expected counts");
		require(c.layouts==null ? c.left!=null && c.right!=null && c.control==null && c.fitLayouts==null
			: c.counts!=null && c.control!=null && c.left==null && c.right==null,
			"Use either both built resources or layouts/counts with an explicit legacy control; do not mix input modes");
		final String[] files={c.left, c.right, c.panels, c.selection, c.references, c.ranks,
			c.identities+".identity.tsv", c.identities+".groups.tsv", c.identities+".identity_bins.tsv", c.counts,
			c.layouts, c.output+".trials.tsv.gz", c.output+".fits.tsv", c.output+".decision.tsv", c.output+".counts.tsv.gz",
			c.layouts==null ? null : c.output+".layouts.tsv", c.layouts==null ? null : c.output+".timing.tsv",
			c.layouts==null ? null : c.output+".comparison.tsv", c.layouts==null ? null : c.output+".pareto.tsv",
			c.layouts==null ? null : c.output+".choices.tsv"};
		require(Tools.testInputFiles(false, true, Arrays.copyOfRange(files, 0, 11)), "All saved input bindings must be readable");
		require(Tools.testOutputFiles(false, false, false, Arrays.copyOfRange(files, 11, files.length)), "Search outputs must be fresh");
		require(Tools.testForDuplicateFiles(true, files), "Inputs and output paths must not alias");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Saved data and geometry invariants require assertions");
		final int zip=ReadWrite.MAX_ZIP_THREADS();
		final boolean zipSet=ReadWrite.SET_ZIP_THREADS;
		final boolean trim=Shared.TRIM_READ_DESCRIPTION;
		try{
			ReadWrite.setZipThreads(1);
			Shared.TRIM_READ_DESCRIPTION=false;
			new OptimizeNcrnaPartitionWeights(c).process();
		}finally{
			ReadWrite.setZipThreads(zip);
			ReadWrite.SET_ZIP_THREADS=zipSet;
			Shared.TRIM_READ_DESCRIPTION=trim;
		}
	}

	private OptimizeNcrnaPartitionWeights(final Config config){c=config;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			sequences.add(read.bases);
		});
		require(names.size()==c.models, "Resource roster must preserve every original cached model");
		final String[] headers=names.toArray(new String[names.size()]);
		final byte[][] bases=sequences.toArray(new byte[sequences.size()][]);
		referenceBases=bases;
		final byte[][] ids=new byte[c.models][];
		for(int m=0; m<ids.length; m++){ids[m]=headers[m].getBytes(StandardCharsets.UTF_8);}
		cohorts=NcrnaDdlCohorts.read(c.identities, data, ids, bases);
		final int[] strata=new int[data.parents.size()];
		for(final Parent p:data.parents){strata[p.index]=p.stratum;}
		final int[][] edges=new int[data.windows.size()][];
		final long[][][] masks=new long[edges.length][][];
		for(int w=0; w<edges.length; w++){
			final Window window=data.windows.get(w);
			edges[w]=new int[window.edges.size()];
			masks[w]=new long[edges[w].length][];
			for(int e=0; e<edges[w].length; e++){
				edges[w][e]=window.edges.get(e).parent.index;
				masks[w][e]=window.edges.get(e).uniqueModels;
			}
		}
		if(c.layouts!=null){
			processLayouts(edges, masks, strata);
			return;
		}
		final NcrnaDdlProfileResource[] resources={NcrnaDdlProfileResource.load(c.left, headers, bases),
			NcrnaDdlProfileResource.load(c.right, headers, bases)};
		for(int v=0; v<2; v++){
			require(resources[v].compact() && resources[v].buckets()==WIDTHS[v], "Only ordered 512/1024 plain7 compact resources enter this search");
		}
		final char[][][] counts=c.counts==null ? capture(resources) : readCounts();
		final NcrnaPartitionWeightSearch[] engines=new NcrnaPartitionWeightSearch[4];
		final int[][] unitRanks=new int[2][];
		for(int v=0; v<4; v++){
			final int terms=3;
			engines[v]=new NcrnaPartitionWeightSearch(c.models, terms, searchCounts(counts[v], terms), edges, masks, strata, cohorts.sides, cohorts.bins, c.randomDraws);
			if(v<2){unitRanks[v]=engines[v].ranks(new double[]{1, 1, 1});}
		}
		checkUnitRanks(unitRanks, METHODS);
		try(Output trials=new Output(c.output+".trials.tsv.gz"); Output fits=new Output(c.output+".fits.tsv");
				Output decisions=new Output(c.output+".decision.tsv")){
			writeHeaders(trials, fits, decisions, false);
			for(int v=0; v<4; v++){
				fitAndWrite(geometry(v), WIDTHS[v%2], engines[v], trials, fits, decisions, null);
			}
		}
	}

	/** Validate the entire interchange before publishing trials. Retain only the
	 * explicit control on the first pass, then fit one layout at a time. This
	 * bounds count memory by one layout rather than the complete recipe grid. */
	private void processLayouts(final int[][] edges, final long[][][] masks, final int[] strata){
		final NcrnaDdlLayoutCounts manifest=new NcrnaDdlLayoutCounts(c.layouts, c.models);
		final boolean[] selected=selectedLayouts(manifest);
		validateLayouts(manifest, edges, masks, strata);
		final SummarizeNcrnaLayoutFits summary=new SummarizeNcrnaLayoutFits(c.layouts, c.randomDraws, c.identities);
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Layout fitting times require current-thread CPU accounting");
		final boolean enabled=clock.isThreadCpuTimeEnabled();
		try{
			if(!enabled){clock.setThreadCpuTimeEnabled(true);}
			manifest.write(c.output+".layouts.tsv");
			try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(c.counts, manifest, data);
					Output trials=new Output(c.output+".trials.tsv.gz"); Output fits=new Output(c.output+".fits.tsv");
					Output decisions=new Output(c.output+".decision.tsv"); Output timing=new Output(c.output+".timing.tsv")){
				writeHeaders(trials, fits, decisions, true);
				timing.writer.println("layout_id\tmember_bins\tterms\tscore_evaluations\tfit_thread_cpu_seconds\tfit_elapsed_seconds\tscope");
				final ByteBuilder timeRow=new ByteBuilder();
				for(int v=0; v<manifest.layouts.length; v++){
					final NcrnaDdlLayoutCounts.Layout layout=manifest.layouts[v];
					final char[][] counts=reader.next(selected[v]);
					if(!selected[v]){continue;}
					final NcrnaPartitionWeightSearch engine=new NcrnaPartitionWeightSearch(c.models, layout.terms,
							prefixCounts(counts, layout.terms), edges, masks, strata, cohorts.sides, cohorts.bins, c.randomDraws, layout.referenceTerm);
					System.err.println("NCRNA_LAYOUT_FIT_START\t"+layout.id);
					final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
					require(cpuStart>=0, "Fitting CPU clock must be enabled before timing");
					final int evaluations=fitAndWrite(layout.id, layout.memberBins, engine, trials, fits, decisions, summary);
					final long cpu=clock.getCurrentThreadCpuTime()-cpuStart, wall=System.nanoTime()-wallStart;
					require(cpu>=0 && wall>=0, "Fitting durations must be nonnegative monotonic measurements");
					timeRow.clear().append(layout.id).tab().append(layout.memberBins).tab().append(layout.terms)
						.tab().append(evaluations).tab().append(cpu/1e9, 9).tab().append(wall/1e9, 9)
						.append("\tFIT_AND_REPORT_ENQUEUE_EXCLUDES_INPUT_AND_CONTROL_VALIDATION").nl();
					timing.writer.print(timeRow);
					System.err.print("NCRNA_LAYOUT_FIT_COMPLETE\t"+timeRow.toString());
				}
				reader.finish();
			}
			summary.write(c.output);
		}finally{if(!enabled){clock.setThreadCpuTimeEnabled(false);}}
	}

	/** Select only complete predeclared searches, never subsets of their trials.
	 * Manifest order remains the output order regardless of flag order. */
	private boolean[] selectedLayouts(final NcrnaDdlLayoutCounts manifest){
		final boolean[] selected=new boolean[manifest.layouts.length];
		if(c.fitLayouts==null){
			Arrays.fill(selected, true);
			return selected;
		}
		for(final String id:c.fitLayouts.split(",", -1)){
			require(!id.isEmpty(), "fitlayouts requires nonempty exact layout IDs");
			int found=-1;
			for(int v=0; v<manifest.layouts.length; v++){
				if(manifest.layouts[v].id.equals(id)){found=v;}
			}
			require(found>=0, "Unknown selected layout: "+id);
			require(!selected[found], "Duplicate selected layout: "+id);
			selected[found]=true;
		}
		return selected;
	}

	private void validateLayouts(final NcrnaDdlLayoutCounts manifest, final int[][] edges,
			final long[][][] masks, final int[] strata){
		final int control=manifest.control(c.control);
		final boolean nativeQueries=manifest.layouts[control].queryForm.equals("QUERY_DDL");
		char[][] controlCounts=null;
		try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(c.counts, manifest, data)){
			for(int v=0; v<manifest.layouts.length; v++){
				final char[][] rows=reader.next(!nativeQueries && v==control);
				if(v==control){controlCounts=rows;}
			}
			reader.finish();
		}
		if(nativeQueries){
			for(final NcrnaDdlLayoutCounts.Layout layout:manifest.layouts){
				require(layout.queryForm.equals("QUERY_DDL") && (layout.legacy || layout.hashNamespace.equals("NATIVE_DDL")),
					"Native object fitting cannot mix a raw probe or historical custom-mapping count row into this study");
			}
			return;// BuildNcrnaDdlLayoutCountsTest checks native control scores independently; old RAW ranks are a different query form.
		}
		final NcrnaDdlLayoutCounts.Layout layout=manifest.layouts[control];
		final NcrnaPartitionWeightSearch controlEngine=new NcrnaPartitionWeightSearch(c.models, layout.terms,
			prefixCounts(controlCounts, layout.terms), edges, masks, strata, cohorts.sides, cohorts.bins, c.randomDraws, layout.referenceTerm);
		final double[] unit=new double[layout.terms];
		Arrays.fill(unit, 1);
		checkUnitRanks(new int[][]{controlEngine.ranks(unit)}, new String[]{METHODS[0]});
	}

	private static void writeHeaders(final Output trials, final Output fits, final Output decisions, final boolean layouts){
		final String prefix=layouts ? "layout_id\tmember_bins" : "geometry\tbins";
		trials.writer.println(prefix+"\tfit\ttrial\tstage\tw0\tw1\tw2\tw3\tpopulation\tparents\tunavailable\tcapped_sum\tavailable_rank_sum\twithin10\tlow_parents\tlow_capped_sum\tlow_available_rank_sum\twithin50\tlow_within10\tlow_within50");
		fits.writer.println(prefix+"\tfit\tw0\tw1\tw2\tw3\tfit_capped_sum\tjudge_capped_sum\tjudge_unit_sum\tjudge_gain\tevaluated_trials\tnear_trials\tw0_min\tw0_max\tw1_min\tw1_max\tw2_min\tw2_max\tw3_min\tw3_max\tunit_within_one_percent\trange_scope");
		decisions.writer.println(prefix+"\tboth_judge_gains_positive\ta_to_b_gain\tb_to_a_gain\tfit_pair_log_distance\trefit_w0\trefit_w1\trefit_w2\trefit_w3\tinterpretation");
	}

	/** Both input modes use exactly the same fitting, judgment and reporting. */
	private static int fitAndWrite(final String geometry, final int bins, final NcrnaPartitionWeightSearch engine,
			final Output trials, final Output fits, final Output decisions, final SummarizeNcrnaLayoutFits summary){
		final Fit a=engine.fit(0);
		final Fit b=engine.fit(1);
		writeFit(geometry, bins, a, trials, fits);
		writeFit(geometry, bins, b, trials, fits);
		final long gainA=gain(a), gainB=gain(b);
		final boolean passed=gainA>0 && gainB>0;
		final Fit refit=passed ? engine.fit(-1) : null;
		if(refit!=null){writeFit(geometry, bins, refit, trials, fits);}
		double distance=0;
		for(int term=0; term<a.winner.weights.length; term++){
			final double delta=StrictMath.log(a.winner.weights[term]/b.winner.weights[term]);
			distance+=delta*delta;
		}
		final ByteBuilder line=new ByteBuilder();
		line.append(geometry).tab().append(bins).tab().append(passed).tab().append(gainA).tab().append(gainB)
			.tab().append(Double.toString(StrictMath.sqrt(distance)));
		appendWeights(line, passed ? refit.winner : a.trials.get(0));
		line.tab().append(!passed ? "UNIT_WEIGHTS_NO_TWO_WAY_GAIN" : a.unitNear && b.unitNear
			? "REFIT_NOT_INDEPENDENTLY_JUDGED;UNIT_WITHIN_ONE_PERCENT_BOTH_FITS" : "REFIT_NOT_INDEPENDENTLY_JUDGED").nl();
		decisions.writer.print(line);
		if(summary!=null){summary.add(geometry, a, b);}
		int initial=0;
		while(initial<a.trials.size() && initial<b.trials.size() && a.trials.get(initial)==b.trials.get(initial)){initial++;}
		assert(initial>0 && initial<=engine.randomDraws+1)
			: "Both fit histories start with the same shared unit/random Trial objects; refinements are independent";
		// Unit/random evaluations are shared across fits; refinements are local.
		return a.trials.size()+b.trials.size()-initial+(refit==null ? 0 : refit.trials.size()-initial);
	}

	/** Historical resources have one mixed consensus term. Old saved origins
	 * collapse into its total; no winner-origin metadata enters the scorer. */
	private char[][] searchCounts(final char[][] source, final int terms){
		if(terms==4){return source;}
		final char[][] result=new char[source.length][];
		for(int w=0; w<source.length; w++){
			result[w]=new char[3*c.models];
			for(int m=0; m<c.models; m++){
				result[w][3*m]=source[w][4*m];
				result[w][3*m+1]=source[w][4*m+1];
				final int consensus=source[w][4*m+2]+source[w][4*m+3];
				require(consensus<=65535, "Historical consensus subcounts must fit one plain-DDL total");
				result[w][3*m+2]=(char)consensus;
			}
		}
		return result;
	}

	/** Interchange reader normalizes each row to its declared active prefix. */
	private char[][] prefixCounts(final char[][] source, final int terms){
		require(terms>=1 && terms<=4, "Only actual active score terms enter a weight search");
		if(terms==4){return source;}
		final char[][] result=new char[source.length][];
		for(int w=0; w<source.length; w++){
			result[w]=new char[terms*c.models];
			for(int m=0; m<c.models; m++){
				for(int term=0; term<terms; term++){result[w][terms*m+term]=source[w][4*m+term];}
				for(int term=terms; term<4; term++){require(source[w][4*m+term]==0, "Inactive count slots must stay zero");}
			}
		}
		return result;
	}

	/** Count once, and save primitive components for subsequent independent re-ranking. */
	private char[][][] capture(final NcrnaDdlProfileResource[] resources){
		assert(resources.length==2) : "Both fixed widths share the same immutable saved-window roster";
		final char[][][] counts=new char[4][data.windows.size()][];
		final int[][] parts={new int[c.models], new int[c.models], new int[c.models]};
		final int[] total=new int[c.models];
		final ByteBuilder line=new ByteBuilder();
		try(Output out=new Output(c.output+".counts.tsv.gz")){
			out.writer.println(COUNT_HEADER);
			for(int v=0; v<4; v++){
				final int width=WIDTHS[v%2];
				final char[][] keys=new char[c.models][width];
				final NcrnaDdlPartitionKeys query=new NcrnaDdlPartitionKeys(7, width, v<2);
				for(int m=0; m<c.models; m++){
					for(int k=0; k<width; k++){keys[m][k]=resources[v%2].key(m, k);}
					final char[] consensus=query.consensusKeys(referenceBases[m]);
					for(int q=0; q<width/4; q++){
						final int bucket=3*width/4+q;
						if(v<2){require(consensus[bucket]==keys[m][bucket], "Consensus reconstruction must preserve each accepted stored key");}
						else{keys[m][bucket]=consensus[bucket];}
					}
				}
				final DDLIndexSmall index=new DDLIndexSmall(keys);
				for(int w=0; w<data.windows.size(); w++){
					final Window window=data.windows.get(w);
					query.build(window.bases);
					query.countParts(index, parts);
					query.count(index, total);
					final char[] row=new char[4*c.models];
					counts[v][w]=row;
					line.clear().append(geometry(v)).tab().append(width).tab().append(w).tab().append(window.input.panel.nameBytes)
						.tab().append(window.input.slot).tab().append(window.event);
					for(int m=0; m<c.models; m++){
						require(parts[0][m]+parts[1][m]+parts[2][m]==total[m], "Partition decomposition must exactly conserve the unit score");
						for(int part=0; part<4; part++){
							final int n=part<2 ? parts[part][m] : part==2 ? 0 : parts[2][m];
							require(n>=0 && n<=(part==1 ? width/2 : width/4), "Integer counts must remain within their own partition");
							row[4*m+part]=(char)n;
							line.tab().append(n);
						}
					}
					out.writer.print(line.nl());
				}
			}
		}
		return counts;
	}

	private char[][][] readCounts(){
		final char[][][] counts=new char[4][data.windows.size()][];
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(c.counts)){
			RrnaResourceIO.header(input.reader, COUNT_HEADER, c.counts);
			for(int v=0; v<4; v++){
				for(int w=0; w<data.windows.size(); w++){
					final byte[] text=input.reader.nextLine();
					require(text!=null, "Saved partition counts must retain every original window");
					p.set(text);
					final Window window=data.windows.get(w);
					require(p.terms()==6+4*c.models && p.termEquals(geometry(v), 0) && p.parseInt(1)==WIDTHS[v%2] && p.parseInt(2)==w
						&& same(p, 3, window.input.panel.nameBytes) && p.parseInt(4)==window.input.slot
						&& p.parseLong(5)==window.event, "Count rows must bind original geometry/width/window/panel/input/event order and complete model components");
					final char[] row=new char[4*c.models];
					counts[v][w]=row;
					for(int i=0; i<row.length; i++){
						final int n=p.parseInt(6+i);
						require(n>=0 && n<=(i%4==1 ? WIDTHS[v%2]/2 : WIDTHS[v%2]/4)
							&& (v<2 || i%4!=2 || n==0), "Saved count exceeds its partition or introduces a forbidden consensus7 term");
						row[i]=(char)n;
					}
					for(int m=0; m<c.models; m++){
						require(row[4*m+2]+row[4*m+3]<=WIDTHS[v%2]/4, "Consensus origin terms together cannot count a stored bin twice");
					}
				}
			}
			require(input.reader.nextLine()==null, "Saved counts cannot append duplicate or extra windows");
		}
		return counts;
	}

	private void checkUnitRanks(final int[][] ranks, final String[] methods){
		assert(ranks.length==methods.length) : "Each independently reconstructed control has one saved method identity";
		final int[] next=new int[methods.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(c.ranks)){
			RrnaResourceIO.header(in.reader, RANK_HEADER, c.ranks);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==10, "Accepted partition ranks retain their ten-column schema");
				int v=-1;
				for(int method=0; method<methods.length; method++){
					if(p.termEquals(methods[method], 0)){v=method;}
				}
				if(v<0){continue;}
				final int i=next[v]++;
				require(i<data.parents.size(), "Accepted ranks contain extra or duplicate selected parents");
				final Parent parent=data.parents.get(i);
				require(same(p, 1, parent.input.panel.nameBytes) && p.parseInt(2)==parent.input.slot
					&& same(p, 3, parent.input.idBytes) && same(p, 4, parent.idBytes)
					&& p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 5)
					&& p.termEquals(parent.nativeRank>0 ? "true" : "false", 7)
					&& p.termEquals(parent.reason, 9) && p.parseInt(8)==ranks[v][i],
					"Unit weights must reproduce every accepted parent rank and original identity before fitting: parent="+i);
			}
		}
		for(final int count:next){require(count==data.parents.size(), "Every unit-weight control must retain every original parent");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean same(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static long gain(final Fit f){
		assert(f.side==0 || f.side==1) : "Only an independently fitted half has an opposite judging half";
		return f.trials.get(0).metrics[1-f.side].cappedSum-f.winner.metrics[1-f.side].cappedSum;
	}
	private static void writeFit(final String geometry, final int bins, final Fit f, final Output trials, final Output fits){
		final String role=f.side<0 ? "REFIT" : f.side==0 ? "A" : "B";
		final ByteBuilder line=new ByteBuilder();
		for(int i=0; i<f.trials.size(); i++){
			final Trial t=f.trials.get(i);
			for(int side=0; side<3; side++){
				final Metrics m=t.metrics[side];
				line.clear().append(geometry).tab().append(bins).tab().append(role).tab().append(i).tab().append(t.stage);
				appendWeights(line, t);
				line.tab().append(side==2 ? "ALL_NOT_SEEN" : side==0 ? "A" : "B").tab().append(m.parents)
					.tab().append(m.unavailable).tab().append(m.cappedSum).tab().append(m.availableRankSum)
					.tab().append(m.within10).tab().append(m.lowParents).tab().append(m.lowCappedSum)
					.tab().append(m.lowAvailableRankSum).tab().append(m.within50)
					.tab().append(m.lowWithin10).tab().append(m.lowWithin50).nl();
				trials.writer.print(line);
			}
		}
		final Trial t=f.winner;
		line.clear().append(geometry).tab().append(bins).tab().append(role);
		appendWeights(line, t);
		line.tab().append(t.metrics[f.side<0 ? 2 : f.side].cappedSum).tab();
		if(f.side<0){line.append("NA\tNA\tNA");}
		else{
			line.append(t.metrics[1-f.side].cappedSum).tab().append(f.trials.get(0).metrics[1-f.side].cappedSum).tab().append(gain(f));
		}
		line.tab().append(f.trials.size()).tab().append(f.nearTrials);
		for(int term=0; term<4; term++){
			final int source=term<t.weights.length ? term : -1;
			line.tab().append(source<0 ? "0.0" : Double.toString(f.range[2*source]));
			line.tab().append(source<0 ? "0.0" : Double.toString(f.range[2*source+1]));
		}
		fits.writer.print(line.tab().append(f.unitNear).append("\tEVALUATED_FIT_TRIALS_ONLY_NOT_CONFIDENCE_INTERVAL\n"));
	}
	private static String geometry(final int variant){return variant<2 ? "CONSENSUS_BOTH_LENGTHS" : "CONSENSUS10_ONLY";}
	private static void appendWeights(final ByteBuilder line, final Trial trial){
		final double[] weights=trial.normalized();
		for(int term=0; term<4; term++){
			final int source=term<weights.length ? term : -1;
			line.tab().append(source<0 ? "0.0" : Double.toString(weights[source]));
		}
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String left, right, panels, selection, references, identities, ranks, counts, layouts, control, fitLayouts, output;
		public int models=508, groups=3153;
		public int randomDraws=5000;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Config c;
	private NcrnaDdlSavedData data;
	private NcrnaDdlCohorts cohorts;
	private byte[][] referenceBases;
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final int[] WIDTHS={512, 1024};
	private static final String[] METHODS={"partition_plain7_masked10_512", "partition_plain7_masked10_1024"};
	private static final String COUNT_HEADER="geometry\twidth\twindow\tpanel\tinput_slot\tevent\tmodel_order_member7_member10_consensus7_consensus10";
	private static final String RANK_HEADER="method\tpanel\tinput_slot\tinput_id\tparent\tseen_stratum\tidentity_bin\tnative_common\tbest_rank\tsource_disposition";
}
