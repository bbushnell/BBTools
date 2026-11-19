package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import map.ObjectMap;
import structures.LongList;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaPartitionWeightSearch.Fit;
import prok.NcrnaPartitionWeightSearch.Metrics;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Reports cross-fitted quality against stored member payload, never refit cost.
 * Native ByteFile lifecycle follows template.A_SampleByteFile. The optimizer
 * supplies actual fitted metrics; the CLI merges disjoint completed summaries.
 * @author Ganyu, Raiden
 */
public final class SummarizeNcrnaLayoutFits{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaLayoutFits.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			int expected=0;
			String layouts=null;
			String seedRoster=null, reference=null;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("expectedlayouts")){expected=Parse.parseIntKMG(b);}
				else if(a.equals("layouts")){layouts=b;}
				else if(a.equals("seedroster")){seedRoster=b;}
				else if(a.equals("reference")){reference=b;}
				else{require(parser.parse(arg, a, b), "Unknown layout-summary parameter: "+arg);}
			}
			require((seedRoster!=null ? parser.in1==null && reference!=null : parser.in1!=null && reference==null)
				&& parser.out1!=null && layouts!=null && expected>0 && !parser.overwrite && !parser.append
				&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
				"Use in=comma-separated comparisons, layouts=canonical metadata, out=fresh prefix and expectedlayouts=complete roster size");
			if(seedRoster==null){merge(parser.in1.split(",", -1), layouts, parser.out1, expected);}
			else{mergeSeeds(seedRoster, layouts, reference, parser.out1, expected);}
			pp.outstream.println(seedRoster==null ? "NCRNA_DDL_LAYOUT_SUMMARY_COMPLETE_NOT_CALLER_ADOPTED"
				: "NCRNA_DDL_MATCHED_RECIPE_SEEDS_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Metadata is small; count streams and trial histories are never copied. */
	public SummarizeNcrnaLayoutFits(final String layouts, final int draws, final String cohorts){
		require(layouts!=null && draws>0 && cohorts!=null && !cohorts.isEmpty(), "Summary binds executed metadata, frozen cohort source and random budget");
		randomDraws=draws;
		cohortSource=cohorts;
		loadMetadata(layouts);
	}

	private SummarizeNcrnaLayoutFits(){
		randomDraws=0;
		cohortSource=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void loadMetadata(final String layouts){
		require(layouts!=null && !layouts.isEmpty(), "The exact expected layout roster must be supplied");
		try(Input input=new Input(layouts)){
			final Columns columns=new Columns(input.reader.nextLine());
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				columns.set(line);
				final Metadata item=new Metadata(columns);
				require(metadata.put(item.id, item)==null, "Summary metadata cannot repeat a layout ID");
			}
		}
		require(!metadata.isEmpty(), "Summary metadata needs at least one declared layout");
	}

	/** Only opposite judging halves enter the front; low-identity is reporting. */
	public void add(final String id, final Fit a, final Fit b){
		final Metadata m=metadata.get(id);
		require(m!=null && a.side==0 && b.side==1, "Summary requires both independently fitted halves of a declared layout");
		final Metrics ua=a.trials.get(0).metrics[1];
		final Metrics ub=b.trials.get(0).metrics[0];
		final Metrics fa=a.winner.metrics[1];
		final Metrics fb=b.winner.metrics[0];
		final long unit=Math.addExact(ua.cappedSum, ub.cappedSum);
		final long fitted=Math.addExact(fa.cappedSum, fb.cappedSum);
		final long gainA=ua.cappedSum-fa.cappedSum;
		final long gainB=ub.cappedSum-fb.cappedSum;
		final boolean pass=gainA>0 && gainB>0;
		final Metrics permittedA=pass ? fa : ua;
		final Metrics permittedB=pass ? fb : ub;
		final long unitLow=Math.addExact(ua.lowAvailableRankSum, ub.lowAvailableRankSum);
		final long fittedLow=Math.addExact(fa.lowAvailableRankSum, fb.lowAvailableRankSum);
		final long unitLowCap=Math.addExact(ua.lowCappedSum, ub.lowCappedSum);
		final long fittedLowCap=Math.addExact(fa.lowCappedSum, fb.lowCappedSum);
		final ByteBuilder out=new ByteBuilder();
		for(final String value:m.description){out.append(value).tab();}
		out.append(m.bytes).tab().append(m.bytes/m.models).tab().append(unit).tab().append(fitted)
			.tab().append(gainA).tab().append(gainB).tab().append(pass).tab().append(pass ? fitted : unit)
			.tab().append(unitLow).tab().append(fittedLow).tab().append(pass ? fittedLow : unitLow)
			.tab().append(unitLowCap).tab().append(fittedLowCap).tab().append(pass ? fittedLowCap : unitLowCap)
			.tab().append(a.unitNear && b.unitNear).tab().append(m.models)
			.tab().append(Math.addExact(ua.parents, ub.parents)).tab().append(Math.addExact(ua.lowParents, ub.lowParents))
			.tab().append(randomDraws).tab().append(NcrnaPartitionWeightSearch.SEED).tab().append(cohortSource).tab().append(m.queryForm)
			.tab().append(Math.addExact(permittedA.within10, permittedB.within10))
			.tab().append(Math.addExact(permittedA.within50, permittedB.within50))
			.tab().append(Math.addExact(permittedA.lowWithin10, permittedB.lowWithin10))
			.tab().append(Math.addExact(permittedA.lowWithin50, permittedB.lowWithin50))
			.tab().append(Math.addExact(permittedA.parents-permittedA.unavailable, permittedB.parents-permittedB.unavailable))
			.tab().append(Math.addExact(permittedA.availableRankSum, permittedB.availableRankSum));
		addRow(out.toBytes());
	}

	public static void merge(final String[] inputs, final String layouts, final String output, final int expected){
		require(inputs!=null && inputs.length>0 && layouts!=null && output!=null && expected>0,
			"Merge requires canonical layout metadata and a complete explicit result roster");
		final ArrayList<String> paths=new ArrayList<String>();
		paths.add(layouts);
		for(final String input:inputs){
			require(input!=null && !input.isEmpty(), "Comparison input paths cannot be empty");
			paths.add(input);
		}
		for(final String suffix:SUFFIXES){paths.add(output+suffix);}
		require(Tools.testForDuplicateFiles(true, paths.toArray(new String[paths.size()])), "Comparison inputs and outputs cannot alias");
		final SummarizeNcrnaLayoutFits result=new SummarizeNcrnaLayoutFits();
		result.loadMetadata(layouts);
		require(result.metadata.size()==expected, "Expected layout count must equal the canonical metadata roster");
		for(final String path:inputs){
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, HEADER, path);
				for(byte[] line; (line=input.reader.nextLine())!=null;){result.addRow(line);}
			}
		}
		require(result.rows.size()==expected, "Merged comparison must contain exactly the declared number of distinct layouts");
		result.write(output);
	}

	private void addRow(final byte[] line){
		final Row row=new Row(line);
		final Metadata expected=metadata.get(row.id);
		require(expected!=null, "Comparison contains an undeclared layout ID: "+row.id);
		expected.validate(row);
		require(ids.add(row.id), "A completed layout cannot be repeated across result batches: "+row.id);
		if(!rows.isEmpty()){
			final Row first=rows.get(0);
			require(row.family.equals(first.family) && row.schema.equals(first.schema) && row.cohorts.equals(first.cohorts)
				&& row.models==first.models && row.parents==first.parents
				&& row.availableParents==first.availableParents
				&& row.lowParents==first.lowParents && row.draws==first.draws && row.seed==first.seed,
				"One front must share representation version, family, model/judge populations and the same search budget/seed");
		}
		rows.add(row);
	}

	/** The row's existing seed field is the WEIGHT-search seed,not the DDL hash
	 * seed. A separate producer-bound roster identifies each sketch replicate. */
	public static void mergeSeeds(final String roster, final String layouts, final String reference,
			final String output, final int expected){
		require(roster!=null && layouts!=null && reference!=null && output!=null && expected>0,
			"Matched recipe comparison requires seed/table roster,metadata,reference recipe and complete expected size");
		RrnaResourceIO.fresh(output, ".seed_summary.tsv", ".seed_choices.tsv");
		final SummarizeNcrnaLayoutFits canonical=new SummarizeNcrnaLayoutFits();
		canonical.loadMetadata(layouts);
		require(canonical.metadata.size()==expected, "Matched seed metadata must retain the complete declared recipe roster");
		final ArrayList<String> paths=new ArrayList<String>(), order=new ArrayList<String>();
		final LongList seeds=new LongList();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(roster)){
			RrnaResourceIO.header(in.reader, "ddl_seed\tcomparison_table", roster);
			for(byte[] text; (text=in.reader.nextLine())!=null;){
				p.set(text);
				require(p.terms()==2 && p.parseLong(0)>=0 && p.parseString(1).startsWith("/"), "Each DDL seed needs one absolute completed comparison path");
				for(int i=0; i<seeds.size; i++){require(seeds.get(i)!=p.parseLong(0), "A DDL seed cannot be counted twice");}
				seeds.add(p.parseLong(0));
				paths.add(p.parseString(1));
			}
		}
		require(seeds.size>=2 && seeds.size<=8, "Matched recipe results need two through eight distinct seeds");
		final ObjectMap<String, Row[]> replicates=new ObjectMap<String, Row[]>(String.class, Row[].class);
		for(int s=0; s<seeds.size; s++){
			try(Input seedFile=new Input(new java.io.File(paths.get(s)).getParent()+"/ddl_seed.txt")){
				final byte[] text=seedFile.reader.nextLine();
				require(text!=null, "A comparison needs its producer's explicit DDL seed sidecar");
				p.set(text);
				require(p.terms()==1 && p.parseLong(0)==seeds.get(s) && seedFile.reader.nextLine()==null, "Comparison and roster DDL seeds must agree");
			}
			try(Input in=new Input(paths.get(s))){
				RrnaResourceIO.header(in.reader, HEADER, paths.get(s));
				for(byte[] text; (text=in.reader.nextLine())!=null;){
					final Row row=new Row(text);
					final Metadata meta=canonical.metadata.get(row.id);
					require(meta!=null, "Seeded comparison cannot introduce an undeclared recipe");
					meta.validate(row);
					Row[] values=replicates.get(row.id);
					if(values==null){values=new Row[seeds.size]; replicates.put(row.id, values); order.add(row.id);}
					require(values[s]==null, "A recipe cannot repeat within one seeded result");
					values[s]=row;
				}
			}
		}
		require(order.size()==expected && replicates.get(reference)!=null, "Every expected recipe and the explicit paired reference must be present");
		final Row[] anchors=replicates.get(reference);
		for(final String id:order){
			for(int s=0; s<seeds.size; s++){
				final Row row=replicates.get(id)[s], anchor=anchors[s];
				require(row!=null && anchor!=null && row.models==anchor.models && row.parents==anchor.parents
					&& row.availableParents==anchor.availableParents && row.lowParents==anchor.lowParents
					&& row.draws==anchor.draws && row.seed==anchor.seed && row.cohorts.equals(anchor.cohorts),
					"Matched recipes require every seed and identical judging populations,cohorts,and weight-search protocol");
				require(anchors[0]!=null && anchor.models==anchors[0].models && anchor.parents==anchors[0].parents
					&& anchor.availableParents==anchors[0].availableParents && anchor.lowParents==anchors[0].lowParents
					&& anchor.draws==anchors[0].draws && anchor.seed==anchors[0].seed && anchor.cohorts.equals(anchors[0].cohorts),
					"DDL seed replicates must use the identical frozen population and weight-search protocol");
			}
		}
		final String[] best=new String[4];
		final long[][] bestValues=new long[4][];
		try(Output all=new Output(output+".seed_summary.tsv"); Output choices=new Output(output+".seed_choices.tsv")){
			all.writer.println("layout_id\tmetric\tseeds\tmean\tminimum\tmaximum\treference_mean\tpaired_improvement_mean\tpaired_improvement_minimum\tpaired_improvement_maximum\timproved_seeds\ttied_seeds\tworsened_seeds");
			choices.writer.println("objective\tlayout_id\tseeds\tmean_available_rank_sum\tmean_low_rank_sum\tmean_within10\tmean_within50\tmean_low_within10\tmean_low_within50\tstored_member_bytes_per_model");
			final ByteBuilder b=new ByteBuilder();
			for(final String id:order){
				final long[] totals=new long[6];
				for(int metric=0; metric<6; metric++){
					long min=Long.MAX_VALUE, max=Long.MIN_VALUE, base=0, gainSum=0, gainMin=Long.MAX_VALUE, gainMax=Long.MIN_VALUE;
					int wins=0, ties=0;
					for(int s=0; s<seeds.size; s++){
						final long value=seedMetric(replicates.get(id)[s], metric), original=seedMetric(anchors[s], metric);
						final long gain=metric<2 ? original-value : value-original;
						totals[metric]+=value;
						base+=original;
						min=Math.min(min, value);
						max=Math.max(max, value);
						gainSum+=gain;
						gainMin=Math.min(gainMin, gain);
						gainMax=Math.max(gainMax, gain);
						if(gain>0){wins++;}else if(gain==0){ties++;}
					}
					all.writer.print(b.clear().append(id).tab().append(SEED_METRICS[metric]).tab().append(seeds.size)
						.tab().append(totals[metric]/(double)seeds.size, 9).tab().append(min).tab().append(max)
						.tab().append(base/(double)seeds.size, 9).tab().append(gainSum/(double)seeds.size, 9)
						.tab().append(gainMin).tab().append(gainMax).tab().append(wins).tab().append(ties).tab().append(seeds.size-wins-ties).nl());
				}
				for(int objective=0; objective<best.length; objective++){
					if(best[objective]==null || (objective<2 ? totals[objective]<bestValues[objective][objective]
						: totals[objective]>bestValues[objective][objective])
						|| totals[objective]==bestValues[objective][objective] && id.compareTo(best[objective])<0){
						best[objective]=id;
						bestValues[objective]=totals.clone();
					}
				}
			}
			for(int objective=0; objective<best.length; objective++){
				b.clear().append(SEED_METRICS[objective]).tab().append(best[objective]).tab().append(seeds.size);
				for(final long value:bestValues[objective]){b.tab().append(value/(double)seeds.size, 9);}
				final Row row=replicates.get(best[objective])[0];
				choices.writer.print(b.tab().append(row.bytes/row.models).nl());
			}
		}
	}
	private static long seedMetric(final Row row, final int metric){
		assert(metric>=0 && metric<6) : "Matched summary preserves the six declared rank/top-k outcomes";
		return metric==0 ? row.availableRankSum : metric==1 ? row.low : metric==2 ? row.within10
			: metric==3 ? row.within50 : metric==4 ? row.lowWithin10 : row.lowWithin50;
	}

	public void write(final String output){
		require(output!=null && !rows.isEmpty(), "A comparison requires completed fits before publishing results");
		RrnaResourceIO.fresh(output, SUFFIXES);
		long best=Long.MAX_VALUE;
		for(final Row row:rows){best=Math.min(best, row.cost);}
		try(Output all=new Output(output+SUFFIXES[0]); Output front=new Output(output+SUFFIXES[1]);
				Output choices=new Output(output+SUFFIXES[2])){
			all.writer.println(HEADER);
			front.writer.println(HEADER);
			for(final Row row:rows){
				all.writer.println(row.line);
				boolean dominated=false;
				for(final Row other:rows){
					if(other.bytes<=row.bytes && other.cost<=row.cost && (other.bytes<row.bytes || other.cost<row.cost)){dominated=true;}
				}
				if(!dominated){front.writer.println(row.line);}
			}
			choices.writer.println("within_percent\tbest_permitted_judged_cost\tlayout_id\tstored_bytes_per_model\tpermitted_judged_cost\tpermitted_low_rank_sum\tpermitted_low_capped_sum\tscope");
			final ByteBuilder line=new ByteBuilder();
			for(final int percent:new int[]{1, 3}){
				Row chosen=null;
				for(final Row row:rows){
					if(Math.multiplyExact(row.cost, 100)>Math.multiplyExact(best, 100+percent)){continue;}
					if(chosen==null || row.bytes<chosen.bytes || row.bytes==chosen.bytes && (row.cost<chosen.cost
						|| row.cost==chosen.cost && row.id.compareTo(chosen.id)<0)){chosen=row;}
				}
				require(chosen!=null, "The best judged layout must lie within its own relative threshold");
				choices.writer.print(line.clear().append(percent).tab().append(best).tab().append(chosen.id)
					.tab().append(chosen.bytes/chosen.models).tab().append(chosen.cost).tab().append(chosen.low)
					.tab().append(chosen.lowCap).append("\tDEVELOPMENT_CROSSFIT_NOT_CALLER_OR_UNTOUCHED_TEST\n"));
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Metadata{
		Metadata(final Columns c){
			id=c.text("layout_id");
			queryForm=c.optional("query_form", "RAW_KEYS");
			require(queryForm.equals("RAW_KEYS") || queryForm.equals("QUERY_DDL"), "The query representation must be explicit and supported");
			models=c.number("model_count");
			final String schema=c.optional("schema_version", "1");
			require(schema.equals("1") || schema.equals("2"), "Only declared historical or plain-key metadata versions may be summarized");
			bytes=Math.addExact(c.number("member_key_payload_bytes"), schema.equals("1") ? c.number("member_origin_payload_bytes") : 0);
			require(models>0 && models<=65536 && bytes>0 && bytes%models==0, "Stored member payload must divide exactly across the bound model roster");
			require(schema.equals("1") || bytes==Math.multiplyExact(2*models, c.number("member_bins")),
				"Every plain-key stored model is exactly two bytes per member bin");
			description=new String[]{id, c.text("family"), c.optional("short_k", "NA_VERSION1"), c.optional("long_k", "NA_VERSION1"),
				c.optional("long_mask", "NA_VERSION1"), c.text("member_mode"), c.text("member_bins"),
				c.alias("member_short_bins", "member7_bins"), c.alias("member_long_bins", "member10_bins"), c.text("consensus_mode"),
				c.text("consensus_bins"), c.optional("consensus_bonus", "0"), c.optional("score_terms", c.text("consensus_mode").equals("BOTH")
					? "member_short,member_long,historical_consensus_short,historical_consensus_long" : "member_short,member_long,historical_consensus_long"),
				schema, c.optional("hash_namespace", "HISTORICAL_VERSION1"), c.text("overflow_policy")};
		}
		void validate(final Row row){
			final LineParser1 p=new LineParser1('\t');
			p.set(row.line);
			for(int i=0; i<description.length; i++){
				require(p.termEquals(description[i], i), "Comparison disagrees with canonical layout metadata: "+id+", field "+i);
			}
			require(row.models==models && row.bytes==bytes && p.termEquals(queryForm, 37),
				"Comparison model count, stored payload and query form must match canonical metadata: "+id);
		}
		final String id, queryForm;
		final String[] description;
		final long models, bytes;
	}

	private static final class Columns{
		Columns(final byte[] header){
			require(header!=null, "A metadata input must contain its named-field header");
			names=new String(header, StandardCharsets.UTF_8).split("\t", -1);
			final HashSet<String> unique=new HashSet<String>();
			for(final String name:names){require(!name.isEmpty() && unique.add(name), "Metadata column names must be unique and nonempty");}
		}
		void set(final byte[] line){
			p.set(line);
			require(p.terms()==names.length, "Metadata rows must match their declared schema");
		}
		int index(final String name){
			for(int i=0; i<names.length; i++){
				if(names[i].equals(name)){return i;}
			}
			return -1;
		}
		String text(final String name){
			final int i=index(name);
			require(i>=0, "Missing metadata column: "+name);
			return p.parseString(i);
		}
		String optional(final String name, final String fallback){return index(name)<0 ? fallback : text(name);}
		String alias(final String name, final String old){return index(name)<0 ? text(old) : text(name);}
		long number(final String name){
			final int i=index(name);
			require(i>=0, "Missing metadata number: "+name);
			return p.parseLong(i);
		}
		final String[] names;
		final LineParser1 p=new LineParser1('\t');
	}

	private static final class Row{
		Row(final byte[] input){
			final LineParser1 p=new LineParser1('\t');
			p.set(input);
			require(p.terms()==44, "Comparison rows require the extended top-k and available-population header; historical missing counts cannot become zeros");
			id=p.parseString(0);
			family=p.parseString(1);
			schema=p.parseString(13);
			bytes=p.parseLong(16);
			models=p.parseLong(31);
			parents=p.parseLong(32);
			lowParents=p.parseLong(33);
			draws=p.parseLong(34);
			seed=p.parseLong(35);
			cohorts=p.parseString(36);
			within10=p.parseLong(38);
			within50=p.parseLong(39);
			lowWithin10=p.parseLong(40);
			lowWithin50=p.parseLong(41);
			availableParents=p.parseLong(42);
			availableRankSum=p.parseLong(43);
			require(p.termEquals("RAW_KEYS", 37) || p.termEquals("QUERY_DDL", 37), "Merged rows must retain their query form");
			cost=p.parseLong(23);
			low=p.parseLong(26);
			lowCap=p.parseLong(29);
			final long unit=p.parseLong(18);
			final long fitted=p.parseLong(19);
			final long gainA=p.parseLong(20);
			final long gainB=p.parseLong(21);
			final boolean pass=gainA>0 && gainB>0;
			require(!id.isEmpty() && !family.isEmpty() && !cohorts.isEmpty() && (schema.equals("1") || schema.equals("2"))
				&& models>0 && models<=65536 && bytes>0 && bytes%models==0
				&& p.parseLong(17)==bytes/models && parents>0 && lowParents>=0 && lowParents<=parents && draws>0,
				"Comparison payload and population accounting must be positive, integral and conserved");
			require(p.termEquals(pass ? "true" : "false", 22) && unit>=0 && fitted>=0 && Math.subtractExact(unit, fitted)==Math.addExact(gainA, gainB)
				&& cost==(pass ? fitted : unit) && cost>=0 && low>=0 && lowCap>=0,
				"The permitted cost follows the two-positive-judge-gains rule, never the all-data refit");
			require(unit<=Math.multiplyExact(parents, NcrnaPartitionWeightSearch.CAP)
				&& fitted<=Math.multiplyExact(parents, NcrnaPartitionWeightSearch.CAP)
				&& lowCap<=Math.multiplyExact(lowParents, NcrnaPartitionWeightSearch.CAP),
				"Capped costs cannot exceed their declared parent population times the fixed cap");
			require(low==p.parseLong(pass ? 25 : 24) && lowCap==p.parseLong(pass ? 28 : 27),
				"Low-identity summaries must use the same permitted choice as the main cost");
			require(availableParents>=0 && availableParents<=parents && availableRankSum>=availableParents
				&& availableRankSum<=Math.multiplyExact(availableParents, models)
				&& within10>=0 && within10<=within50 && within50<=availableParents
				&& lowWithin10>=0 && lowWithin10<=lowWithin50 && lowWithin50<=lowParents
				&& lowWithin10<=within10 && lowWithin50<=within50,
				"Top-k counts use the same judging parents and low-identity subset as the permitted rank cost");
			require(schema.equals("1") || bytes==Math.multiplyExact(2*models, p.parseLong(6)),
				"Merged plain-key results cannot add bytes beyond their member bins");
			line=input.clone();
		}
		final String id, family, schema, cohorts;
		final long models, parents, lowParents, draws, seed, bytes, cost, low, lowCap;
		final long within10, within50, lowWithin10, lowWithin50;
		final long availableParents, availableRankSum;
		final byte[] line;
	}

	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int randomDraws;
	private final String cohortSource;
	private final HashMap<String, Metadata> metadata=new HashMap<String, Metadata>();
	private final HashSet<String> ids=new HashSet<String>();
	private final ArrayList<Row> rows=new ArrayList<Row>();

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final String[] SUFFIXES={".comparison.tsv", ".pareto.tsv", ".choices.tsv"};
	private static final String[] SEED_METRICS={"available_rank_sum","low_available_rank_sum","within10","within50","low_within10","low_within50"};
	public static final String HEADER="layout_id\tfamily\tshort_k\tlong_k\tlong_mask\tmember_mode\tmember_bins\tmember_short_bins\tmember_long_bins\tconsensus_mode\tconsensus_bins\tconsensus_bonus\tscore_terms\tschema_version\thash_namespace\toverflow_policy\tstored_member_bytes\tstored_bytes_per_model\tunit_judged_cost\tfitted_judged_cost\ta_to_b_gain\tb_to_a_gain\tboth_gains_positive\tpermitted_judged_cost\tunit_low_rank_sum\tfitted_low_rank_sum\tpermitted_low_rank_sum\tunit_low_capped_sum\tfitted_low_capped_sum\tpermitted_low_capped_sum\tunit_within_one_percent_both_fits\tmodel_count\tjudge_parents\tlow_identity_parents\trandom_draws\tseed\tcohort_source\tquery_form"+
		"\tpermitted_within10\tpermitted_within50\tpermitted_low_within10\tpermitted_low_within50\tavailable_judge_parents\tpermitted_available_rank_sum";
}
