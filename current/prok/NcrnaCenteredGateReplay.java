package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;

/** Replays paired real-call gate caches with physical guards. Internal fold4
 * alone chooses each endpoint cutoff; outer calls use frozen choices. No
 * network inference, fitting, rematching or caller execution occurs here.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateReplay{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateReplay.class, false);
		try{
			new NcrnaCenteredGateReplay(pp.args).process();
			pp.outstream.println("NCRNA_CENTERED_GATE_REPLAY_COMPLETE_CONDITIONAL_NO_INFERENCE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaCenteredGateReplay(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){
				require("select".equals(b) || "evaluate".equals(b), "mode=select or evaluate is required");
				select=b.equals("select");
			}else if(a.equals("member")){member=b;}
			else if(a.equals("seed")){seed=Integer.parseInt(b);}
			else if(a.equals("expected")){expected=Long.parseLong(b);}
			else if(a.equals("choices")){choices=b;}
			else{require(parser.parse(arg, a, b), "Unknown real centered-gate replay option: "+arg);}
		}
		five=parser.in1;
		three=parser.in2;
		output=parser.out1;
		require(five!=null && three!=null && output!=null && expected>0 && seed>=58 && seed<=60
			&& ("RF00174".equals(member) || "RF01482".equals(member)) && (select ? choices==null : choices!=null),
			"Bind both complete end caches, member, gate seed, accepted call count and fresh output; evaluate requires frozen choices");
		require(!parser.overwrite && !parser.append && parser.out2==null && parser.maxReads<0 && parser.samplerate==1
			&& !new File(output).exists() && Tools.testInputFiles(false, true, five, three, choices), "Coordinate replay preserves complete populations without replacement");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Cache Validation       ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		readCalls();
		require(!calls.isEmpty(), "Selected replay population must contain actual matched calls");
		if(!select){
			final String panel=calls.get(0).dataset;
			for(final Call call:calls){require(call.dataset.equals(panel), "Each outer replay keeps one original panel denominator; do not pool panels");}
		}
		require(new File(output).mkdir(), "Cannot create fresh replay output directory");
		if(select){selectCutoffs();}else{readChoices();}
		writeCurves();
		writeResults();
	}

	private void readCalls(){
		final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		long rows=0;
		final long[] scored=new long[2];
		boolean complete=false;
		try(Input first=new Input(five); Input second=new Input(three)){
			gateNetworks[0]=readHeader(first, 0);
			gateNetworks[1]=readHeader(second, 1);
			for(byte[] line; (line=first.reader.nextLine())!=null;){
				final byte[] other=second.reader.nextLine();
				require(other!=null, "Paired caches must retain the identical complete real-call cohort");
				a.set(line);
				b.set(other);
				if(a.termEquals("#COMPLETE", 0)){
					require(a.terms()==3 && b.terms()==3 && b.termEquals("#COMPLETE", 0)
						&& a.parseLong(1)==rows && b.parseLong(1)==rows && rows==expected
						&& a.parseLong(2)==scored[0] && b.parseLong(2)==scored[1]
						&& first.reader.nextLine()==null && second.reader.nextLine()==null,
						"Both caches must close once at accepted matched-call and separate inference counts");
					complete=true;
					break;
				}
				require(a.terms()==53 && b.terms()==53, "The real gate score layout has53explicit fields");
				a.setBounds(21);
				b.setBounds(21);
				require(a.b()==b.b(), "Both ends must share exact dataset/call/model/raw/truth/center geometry");
				for(int i=0; i<a.b(); i++){require(line[i]==other[i], "Both ends must retain identical real-call geometry before the end-specific fields");}
				final int fold=a.parseInt(1);
				require(select ? fold>=0 && fold<5 : fold==-1, "Fold4 selection cannot consume outer calls, and outer evaluation cannot consume internal calls");
				final Call call=new Call();
				call.dataset=a.parseString(0);
				call.assembly=a.parseString(4);
				call.reference=a.parseString(5);
				require(seen.add(call.assembly+'\t'+call.reference), "Each assembly-qualified reference contributes once to this panel");
				require(a.termEquals(member, 6) && (a.termEquals("+", 9) || a.termEquals("-", 9)), "Cached member and gene strand must be explicit");
				call.length=a.parseInt(10);
				call.fixed=a.parseInt(11);
				call.minimum=a.parseInt(16);
				call.maximum=NcrnaTextIO.maximumLength(a, 17, legacyMaximum[0]);
				require(call.maximum==NcrnaTextIO.maximumLength(b, 17, legacyMaximum[1]), "Paired caches must retain the same family maximum");
				call.support=a.parseLong(20);
				call.centerBad=geometryCode(a, 21);
				require(call.length>0 && call.fixed>=0 && call.fixed<=3 && call.minimum==109 && call.maximum==Integer.MAX_VALUE
					&& call.support>=0, "Real cobalamin caches retain physical length bounds, fixed bits and model support");
				for(int end=0; end<2; end++){
					call.raw[end]=a.parseInt(12+end);
					call.truth[end]=a.parseInt(14+end);
					call.center[end]=a.parseInt(18+end);
				}
				require(validPair(call.raw, call.length) && validPair(call.truth, call.length) && validPair(call.center, call.length),
					"Raw, truth and retained center pairs must fit the actual contig; short raw spans stay in the cohort");
				require(call.centerBad==0 ? geometry(call, call.center[0], call.center[1])==0
					: Arrays.equals(call.raw, call.center), "Invalid proposed centers retain the entire raw pair");
				for(int end=0; end<2; end++){
					readEnd(end==0 ? a : b, call, end);
					if(call.hasGateScore(end)){scored[end]++;}
				}
				rows++;
				if(!select || fold==4){calls.add(call);}
			}
		}
		require(complete, "A truncated gate cache cannot be interpreted as a complete result");
	}

	private String readHeader(final Input in, final int end){
		final byte[] magic=in.reader.nextLine();
		require(magic!=null, "Gate cache requires a version header");
		final LineParser1 version=new LineParser1('\t');
		version.set(magic);
		legacyCache[end]=version.termEquals(ScoreNcrnaCenteredGate.LEGACY_MAGIC, 0);
		legacyMaximum[end]=legacyCache[end] || version.termEquals(ScoreNcrnaCenteredGate.NUMERIC_MAXIMUM_MAGIC, 0);
		require(version.terms()==1 && (legacyMaximum[end] || version.termEquals(ScoreNcrnaCenteredGate.MAGIC, 0)),
			"Unsupported gate-cache version");
		RrnaResourceIO.header(in.reader, "#schema\t"+NcrnaCenteredGateFeatures.SCHEMA, output);
		RrnaResourceIO.header(in.reader, "#member\t"+member, output);
		RrnaResourceIO.header(in.reader, "#end\t"+ENDS[end], output);
		RrnaResourceIO.header(in.reader, "#seed\t"+seed, output);
		RrnaResourceIO.header(in.reader, "#partition\t"+(select ? "internal" : "outer"), output);
		final byte[] row=in.reader.nextLine();
		require(row!=null, "Gate cache requires its actual fitted network binding");
		final LineParser1 p=new LineParser1('\t');
		p.set(row);
		require(p.terms()==2 && p.termEquals("#network", 0) && p.termStartsWith("/", 1), "Fitted gate provenance must name an absolute network path");
		final String net=p.parseString(1);
		if(!legacyCache[end]){RrnaResourceIO.header(in.reader, ScoreNcrnaCenteredGate.OPTIONAL_HEADER, output);}
		RrnaResourceIO.header(in.reader, ScoreNcrnaCenteredGate.cacheHeader(legacyMaximum[end]).toString(), output);
		return net;
	}

	private void readEnd(final LineParser1 p, final Call call, final int end){
		require(p.termEquals(ENDS[end], 22), "Each paired cache has its declared endpoint");
		final String status=p.parseString(23);
		call.proposal[end]=p.parseLong(25);
		call.missing[end]=status.equals("MISSING_FEATURE_HOLD");
		final boolean available=!absent(p, 31, legacyCache[end]);
		if(available){
			call.heads[end]=new float[22];
			for(int i=0; i<22; i++){
				call.heads[end][i]=Float.intBitsToFloat(p.parseInt(31+i));
				require(Float.isFinite(call.heads[end][i]), "Cached first-stage evidence must remain finite");
			}
		}else{for(int i=31; i<53; i++){require(absent(p, i, legacyCache[end]), "Unavailable first-stage evidence must be absent in every field");}}
		final int neural=available ? NcrnaCenteredGateFeatures.adjustment(call.heads[end]) : 0;
		require(!available || call.centerBad==0 && (call.fixed&(1<<end))==0, "Fixed or invalid-center ends must never carry network inference");
		require(p.parseInt(26)==neural && p.parseLong(27)==(long)call.center[end]-call.raw[end], "Cached centered displacements must agree with original coordinates and scores");
		if(status.equals("PROPOSAL")){
			require(available && p.parseLong(24)>0 && !absent(p, 30, legacyCache[end]) && call.centerBad==0 && (call.fixed&(1<<end))==0
				&& call.proposal[end]==(long)call.center[end]+neural && call.proposal[end]!=call.raw[end], "Only real unfixed whole-move candidates carry a gate score");
			call.score[end]=Float.intBitsToFloat(p.parseInt(30));
			require(Float.isFinite(call.score[end]) && p.parseInt(28)==(NcrnaCenteredGateFeatures.improves(call.raw[end], call.proposal[end], call.truth[end]) ? 1 : 0),
				"Gate score and strict-improvement label must remain valid");
		}else{
			require(p.parseLong(24)==0 && absent(p, 30, legacyCache[end]) && p.parseInt(28)==0 && call.proposal[end]==call.raw[end]
				&& (status.equals("FIXED") && (call.fixed&(1<<end))!=0 && !available
				|| status.equals("INVALID_CENTER_RAW_HOLD") && call.centerBad!=0 && !available
				|| call.missing[end] && call.centerBad==0 && !available
				|| status.equals("SAME_AS_RAW_HOLD") && available && (long)call.center[end]+neural==call.raw[end]),
				"Deterministic holds must preserve their actual reason and stay at RAW");
		}
		require(geometryCode(p, 29)==geometry(call, end==0 ? call.proposal[end] : call.raw[0], end==1 ? call.proposal[end] : call.raw[1]),
			"Cached single-end geometry must use the other RAW endpoint");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Cutoffs and Results     ----------------*/
	/*--------------------------------------------------------------*/

	private void selectCutoffs(){
		final ByteBuilder row=new ByteBuilder();
		{
			for(int end=0; end<2; end++){
				final FloatList list=new FloatList();
				for(final Call call:calls){if(call.hasGateScore(end)){list.add(call.score[end]);}}
				final float[] cuts=list.toArray();
				Arrays.sort(cuts);
				final ArrayList<Choice> points=new ArrayList<Choice>();
				for(int i=0; i<cuts.length; i++){
					if(i>0 && cuts[i]==cuts[i-1]){continue;}
					points.add(new Choice("CURVE", end, "POINT", cuts[i], false, single(end, 4, cuts[i], false), -1));
				}
				points.add(new Choice("CURVE", end, "POINT", 0, true, single(end, 4, 0, true), -1));
				curvePoints.addAll(points);
				Stats best=null;
				for(final Choice point:points){
					if(best==null || compare(point.stats, best)>0){best=point.stats;}
				}
				for(final Choice point:points){
					if(compare(point.stats, best)==0){selected.add(new Choice("BEST", end, "ALL_TIED_BEST", point.cutoff, point.never, point.stats, -1));}
				}
				addMatched(points, end, "RATE_MARGIN0", single(end, 2, 0, false).moved);
				addMatched(points, end, "RATE_MARGIN01", single(end, 3, 0, false).moved);
			}
		}
		try(Output out=new Output(path("choices.tsv"))){
			out.writer.print(new ByteBuilder(CHOICE_MAGIC).append("\n#member\t").append(member).append("\n#seed\t").append(seed)
				.append("\n#network5\t").append(gateNetworks[0]).append("\n#network3\t").append(gateNetworks[1]).nl());
			out.writer.println(CHOICE_HEADER);
			for(final Choice point:selected){
				row.clear().append(point.policy).tab().append(ENDS[point.end]).tab().append(point.side).tab().append(point.stats==null ? 0 : 1).tab();
				if(point.stats==null){row.append("\t\t");}
				else{
					appendCutoff(row, point.cutoff, point.never, true);
					row.tab();
					appendCutoff(row, point.cutoff, point.never, false);
					row.tab().append(point.stats.moved);
				}
				row.tab();
				if(point.target>=0){row.append(point.target);}
				row.tab();
				if(point.stats==null){row.tab();}else{row.append(point.stats.better).tab().append(point.stats.worse);}
				row.tab();
				if(point.stats!=null){row.append(point.never ? "NEVER" : "THRESHOLD");}
				out.writer.print(row.nl());
			}
		}
	}

	/** Presentation of every observed cutoff, with frozen monitor choices marked.
	 * Outer curves describe saved calls only: they never create or modify choices.
	 * Include the original monitor cutoffs even when no outer score equals them. */
	private void writeCurves(){
		if(!select){
			for(int end=0; end<2; end++){
				final FloatList values=new FloatList();
				for(final Call call:calls){if(call.hasGateScore(end)){values.add(call.score[end]);}}
				for(final Choice choice:selected){
					if(choice.end==end && choice.stats!=null && !choice.never){values.add(choice.cutoff);}
				}
				final float[] cutoffs=values.toArray();
				Arrays.sort(cutoffs);
				for(int i=0; i<cutoffs.length; i++){
					if(i>0 && cutoffs[i]==cutoffs[i-1]){continue;}
					curvePoints.add(new Choice("CURVE", end, "POINT", cutoffs[i], false, single(end, 4, cutoffs[i], false), -1));
				}
				curvePoints.add(new Choice("CURVE", end, "POINT", 0, true, single(end, 4, 0, true), -1));
			}
		}
		final ByteBuilder row=new ByteBuilder();
		try(Output out=new Output(path("curves.tsv"))){
			out.writer.println("member\tseed\tend\tcutoff_bits_if_threshold\tcutoff_if_threshold\tmatched\tbetter\tworse\ttied\tmoved\traw_rms_bp\trms_bp\trms_minus_raw_bp\texact\tmove_rate\tqualification\tmonitor_marks\tpanel\tscope\tcutoff_mode");
			for(final Choice point:curvePoints){
				final Stats stats=point.stats;
				require(stats.n==calls.size(), "Cutoff curves preserve the full frozen cohort, including deterministic holds");
				row.clear().append(member).tab().append(seed).tab().append(ENDS[point.end]).tab();
				appendCutoff(row, point.cutoff, point.never, true);
				row.tab();
				appendCutoff(row, point.cutoff, point.never, false);
				row.tab().append(stats.n).tab().append(stats.better).tab().append(stats.worse).tab().append(stats.tied).tab().append(stats.moved)
					.tab().append(stats.rawRms(), 8).tab().append(stats.rms(), 8).tab().append(stats.rms()-stats.rawRms(), 8)
					.tab().append(stats.exact).tab().append((double)stats.moved/stats.n, 8).tab()
					.append(select ? "INTERNAL_MONITOR_SELECTION" : "OUTER_DESCRIPTIVE_NO_SELECTION").tab();
				boolean marked=false;
				for(final Choice choice:selected){
					if(choice.end!=point.end || choice.stats==null || choice.never!=point.never || !choice.never && choice.cutoff!=point.cutoff){continue;}
					if(marked){row.append(';');}
					row.append(choice.policy).append(':').append(choice.side);
					marked=true;
				}
				row.tab().append(select ? "fold4" : calls.get(0).dataset).append("\tSINGLE_OTHER_RAW\t").append(point.never ? "NEVER" : "THRESHOLD");
				out.writer.print(row.nl());
			}
		}
	}

	/** Counts attainable PHYSICAL single-end moves, not score-rank positions.
	 * Invalid proposals have zero weight. A control may exceed the gate's maximum
	 * attainable count; retain the available bracket and an explicit missing side. */
	private void addMatched(final ArrayList<Choice> points, final int end, final String policy, final long target){
		long below=-1, above=Long.MAX_VALUE;
		for(final Choice point:points){
			final long moves=point.stats.moved;
			if(moves<=target && moves>below){below=moves;}
			if(moves>=target && moves<above){above=moves;}
		}
		for(final Choice point:points){
			if(point.stats.moved==below){selected.add(new Choice(policy, end, below==above ? "EXACT" : "LOW", point.cutoff, point.never, point.stats, target));}
			if(above!=below && point.stats.moved==above){selected.add(new Choice(policy, end, "HIGH", point.cutoff, point.never, point.stats, target));}
		}
		if(below<0){selected.add(new Choice(policy, end, "LOW_UNAVAILABLE", 0, false, null, target));}
		if(above==Long.MAX_VALUE){selected.add(new Choice(policy, end, "HIGH_UNAVAILABLE", 0, false, null, target));}
	}

	private void readChoices(){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(choices)){
			final byte[] magic=input.reader.nextLine();
			require(magic!=null, "Frozen choices require a version header");
			p.set(magic);
			final boolean legacy=p.termEquals(LEGACY_CHOICE_MAGIC, 0);
			require(p.terms()==1 && (legacy || p.termEquals(CHOICE_MAGIC, 0)), "Unsupported frozen-choice version");
			RrnaResourceIO.header(input.reader, "#member\t"+member, choices);
			RrnaResourceIO.header(input.reader, "#seed\t"+seed, choices);
			RrnaResourceIO.header(input.reader, "#network5\t"+gateNetworks[0], choices);
			RrnaResourceIO.header(input.reader, "#network3\t"+gateNetworks[1], choices);
			RrnaResourceIO.header(input.reader, legacy ? LEGACY_CHOICE_HEADER : CHOICE_HEADER, choices);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==(legacy ? 10 : 11) && (p.termEquals("5prime", 1) || p.termEquals("3prime", 1)), "Frozen choices require their complete schema and endpoint");
				final String policy=p.parseString(0), side=p.parseString(2);
				require(policy.equals("BEST") || policy.equals("RATE_MARGIN0") || policy.equals("RATE_MARGIN01"), "Unknown frozen cutoff policy");
				final int end=p.termEquals("5prime", 1) ? 0 : 1;
				final long target=absent(p, 7, legacy) ? -1 : p.parseLong(7);
				require(policy.equals("BEST") ? target==-1 : target>=0, "Matched movement targets must come from the internal control population");
				if(p.termEquals("0", 3)){
					require(!policy.equals("BEST") && (side.equals("LOW_UNAVAILABLE") || side.equals("HIGH_UNAVAILABLE"))
						&& absent(p, 4, legacy) && absent(p, 5, legacy) && absent(p, 6, legacy)
						&& absent(p, 8, legacy) && absent(p, 9, legacy) && (legacy || p.length(10)==0), "Unavailable brackets cannot invent a cutoff or monitoring counts");
					selected.add(new Choice(policy, end, side, 0, false, null, target));
					continue;
				}
				require(p.termEquals("1", 3), "A cutoff availability flag must be0or1");
				final boolean never=legacy ? p.parseInt(4)==0x7f800000 : p.termEquals("NEVER", 10);
				final float cut;
				if(never){
					require(legacy ? p.termEquals("Infinity", 5) : p.length(4)==0 && p.length(5)==0,
						"Never-move has no numeric cutoff; legacy infinity must be exact");
					cut=0;
				}else{
					require(legacy || p.termEquals("THRESHOLD", 10), "An available gate choice needs an explicit cutoff mode");
					cut=Float.intBitsToFloat(p.parseInt(4));
					require(Float.isFinite(cut) && cut==Float.parseFloat(p.parseString(5)), "Frozen finite cutoff bits and readable cutoff must agree");
				}
				final Stats stats=new Stats();
				stats.moved=p.parseLong(6);
				stats.better=p.parseLong(8);
				stats.worse=p.parseLong(9);
				require(stats.moved>=0 && stats.better>=0 && stats.worse>=0 && stats.better+stats.worse<=stats.moved,
					"Monitoring movement and outcome counts must be physically possible");
				require(!never || stats.moved==0 && stats.better==0 && stats.worse==0,
					"A frozen never-move choice must retain zero actual movements and changed errors");
				selected.add(new Choice(policy, end, side, cut, never, stats, target));
			}
		}
		for(int end=0; end<2; end++){
			boolean found=false;
			for(final Choice point:selected){found|=point.end==end && point.policy.equals("BEST") && point.stats!=null;}
			require(found, "Both endpoints need their frozen monitor-selected best cutoff set");
		}
	}
	private void writeResults(){
		final ByteBuilder row=new ByteBuilder();
		try(Output out=new Output(path("summary.tsv")); Output absent=new Output(path("unavailable.tsv"))){
			out.writer.println(RESULT_HEADER);
			absent.writer.println("member\tseed\tpolicy\tend\tside\tmonitor_target_moves\tqualification");
			for(int policy=0; policy<CONTROL_NAMES.length; policy++){
				final int calculation=policy==4 ? ALL_PROPOSALS : policy;
				for(int end=0; end<2; end++){
					writeStats(out, row, CONTROL_NAMES[policy], null, null, false, end, single(end, calculation, 0, false));
				}
				final Stats[] pair=joint(calculation, 0, false, 0, false);
				for(int end=0; end<2; end++){
					writeStats(out, row, CONTROL_NAMES[policy], null, null, true, end, pair[end]);
				}
			}
			for(final Choice choice:selected){
				if(choice.stats==null){
					absent.writer.print(row.clear().append(member).tab().append(seed).tab().append(choice.policy).tab().append(ENDS[choice.end])
						.tab().append(choice.side).tab().append(choice.target).append("\tMONITOR_UNATTAINABLE_NOT_RETUNED_ON_OUTER\n"));
					continue;
				}
				writeStats(out, row, choice.policy, choice.end==0 ? choice : null, choice.end==1 ? choice : null,
					false, choice.end, single(choice.end, 4, choice.cutoff, choice.never));
			}
			for(final Choice a:selected){
				if(a.end!=0 || a.stats==null){continue;}
				for(final Choice b:selected){
					if(b.end!=1 || b.stats==null || !a.policy.equals(b.policy)){continue;}
					final Stats[] pair=joint(4, a.cutoff, a.never, b.cutoff, b.never);
					for(int end=0; end<2; end++){
						writeStats(out, row, a.policy, a, b, true, end, pair[end]);
					}
				}
			}
		}
	}

	private Stats[] joint(final int policy, final float fiveCutoff, final boolean fiveNever, final float threeCutoff, final boolean threeNever){
		final Stats[] stats={new Stats(), new Stats()};
		final int[] result=new int[2];
		for(final Call call:calls){
			final long start=proposal(call, 0, policy, fiveCutoff, fiveNever);
			final long stop=proposal(call, 1, policy, threeCutoff, threeNever);
			final int bad;
			if(policy==2 || policy==3 || policy==5){
				// The existing first-stage caller tests both neural adjustments
				// from one guarded C and restores C when that pair is invalid.
				bad=start==call.center[0] && stop==call.center[1] ? 0 : geometry(call, start, stop);
				result[0]=bad==0 ? (int)start : call.center[0];
				result[1]=bad==0 ? (int)stop : call.center[1];
			}else{
				bad=apply(call.raw[0], call.raw[1], start, stop, true, true, call.length, call.minimum, call.maximum, result);
			}
			for(int end=0; end<2; end++){stats[end].add(call, result, end, bad);}
		}
		return stats;
	}

	private void writeStats(final Output out, final ByteBuilder row, final String policy,
			final Choice fiveChoice, final Choice threeChoice, final boolean joint, final int end, final Stats stats){
		require(stats.n==calls.size(), "Every replay policy must retain exactly the same matched-call cohort");
		row.clear().append(member).tab().append(seed).tab().append(select ? "fold4" : calls.get(0).dataset)
			.tab().append(policy).tab().append(fiveChoice==null ? "" : fiveChoice.side).tab().append(threeChoice==null ? "" : threeChoice.side).tab();
		if(fiveChoice!=null){appendCutoff(row, fiveChoice.cutoff, fiveChoice.never, false);}
		row.tab();
		if(threeChoice!=null){appendCutoff(row, threeChoice.cutoff, threeChoice.never, false);}
		row.tab().append(joint ? "PAIRED_JOINT" : "SINGLE_OTHER_RAW")
			.tab().append(ENDS[end]).tab().append(stats.n).tab().append(stats.exact).tab().append(stats.within3).tab().append(stats.within10);
		if(joint){row.tab().append(stats.bothExact).tab().append(stats.both3).tab().append(stats.both10);}
		else{row.append("\t\t\t");}
		row.tab().append(stats.rawRms(), 8).tab().append(stats.rms(), 8).tab().append(stats.rms()-stats.rawRms(), 8)
			.tab().append(stats.better).tab().append(stats.worse).tab().append(stats.tied).tab().append(stats.moved)
			.tab().append(stats.exactHarmed).tab().append(stats.near).tab();
		if(stats.near>0){row.append(Math.sqrt(stats.nearRawSquared/stats.near), 8).tab().append(Math.sqrt(stats.nearSquared/stats.near), 8);}
		else{row.tab();}
		row.tab().append(stats.nearHarmed).tab().append(stats.distantImproved).tab().append(stats.distantHarmed)
			.tab().append(stats.fixed).tab().append(stats.unsupported).tab().append(stats.missing).tab().append(stats.centerHeld)
			.tab().append(stats.n-stats.invalid[0]).tab().append(stats.invalid[1]).tab().append(stats.invalid[2])
			.tab().append(stats.invalid[3]).tab().append(stats.invalid[4]).tab();
		if(stats.better+stats.worse>0){row.append((double)stats.better/(stats.better+stats.worse), 8);}
		row.tab().append(Double.toString(stats.rawSquared)).tab().append(Double.toString(stats.squared))
			.tab().append(Double.toString(stats.nearRawSquared)).tab().append(Double.toString(stats.nearSquared))
			.append("\tCONDITIONAL_REAL_CALLS_NOT_CALLER_OR_NESTED_MONITOR\t")
			.append(fiveChoice==null ? "" : fiveChoice.never ? "NEVER" : "THRESHOLD").tab()
			.append(threeChoice==null ? "" : threeChoice.never ? "NEVER" : "THRESHOLD").nl();
		out.writer.print(row);
	}

	private Stats single(final int end, final int policy, final float cutoff, final boolean never){
		final Stats stats=new Stats();
		final int[] result=new int[2];
		for(final Call call:calls){
			final long proposed=proposal(call, end, policy, cutoff, never);
			final int bad=apply(call.raw[0], call.raw[1], end==0 ? proposed : call.raw[0], end==1 ? proposed : call.raw[1],
				true, true, call.length, call.minimum, call.maximum, result);
			stats.add(call, result, end, bad);
		}
		return stats;
	}

	private static long proposal(final Call call, final int end, final int policy, final float cutoff, final boolean never){
		assert(policy>=0 && policy<=ALL_PROPOSALS) : "Replay policies distinguish the numeric gate from the all-proposal control";
		if(policy==0){return call.raw[end];}
		if(policy==1){return call.center[end];}
		if(policy==ALL_PROPOSALS){return call.hasGateScore(end) ? call.proposal[end] : call.raw[end];}
		if(policy==4){return !never && call.hasGateScore(end) && call.score[end]>=cutoff ? call.proposal[end] : call.raw[end];}
		if(policy==5){return (long)call.center[end]+(call.heads[end]==null ? 0
			: NcrnaEndpointConfidence.selectFull(call.heads[end], 10, 0, 0, false)-10);}
		return (long)call.center[end]+(call.heads[end]==null ? 0
			: NcrnaEndpointConfidence.selectFull(call.heads[end], 10, policy==2 ? 0 : .1f)-10);
	}

	private static int compare(final Stats a, final Stats b){
		int comparison=Long.compare(a.better-a.worse, b.better-b.worse);
		if(comparison==0){comparison=Long.compare(b.worse, a.worse);}
		if(comparison==0){comparison=Long.compare(b.moved, a.moved);}
		return comparison;
	}

	/** All accepted gate movements are tested together before either end changes.
	 * Rejection and invalid joint geometry both restore the original RAW pair. */
	public static int apply(final int rawStart, final int rawStop, final long proposalStart, final long proposalStop,
			final boolean takeStart, final boolean takeStop, final int length, final int minimum, final int maximum, final int[] result){
		require(rawStart>=0 && rawStop>=rawStart && rawStop<length && result!=null && result.length==2,
			"Gate replay requires a real raw pair and an owned two-end result buffer");
		final long start=takeStart ? proposalStart : rawStart;
		final long stop=takeStop ? proposalStop : rawStop;
		final int bad=start==rawStart && stop==rawStop ? 0 : EvaluateNcrnaRealEndpoints.geometry(start, stop, length, minimum, maximum);
		result[0]=bad==0 ? (int)start : rawStart;
		result[1]=bad==0 ? (int)stop : rawStop;
		return bad;
	}
	/** Only the versioned legacy input accepts NA; new optional numeric cells are empty. */
	private static boolean absent(final LineParser1 p, final int field, final boolean legacy){
		return legacy ? p.termEquals("NA", field) : p.length(field)==0;
	}

	/** Never-move is a mode, with no numeric threshold to serialize. */
	private static void appendCutoff(final ByteBuilder row, final float cutoff, final boolean never, final boolean bits){
		require(Float.isFinite(cutoff), "A computed or selected cutoff must be finite before output");
		if(!never){
			if(bits){row.append(Float.floatToRawIntBits(cutoff));}else{row.append(Float.toString(cutoff));}
		}
	}

	private static boolean validPair(final int[] pair, final int length){return pair[0]>=0 && pair[1]>=pair[0] && pair[1]<length;}
	private static int geometry(final Call call, final long start, final long stop){
		return EvaluateNcrnaRealEndpoints.geometry(start, stop, call.length, call.minimum, call.maximum);
	}
	private static int geometryCode(final LineParser1 p, final int field){
		for(int i=0; i<GEOMETRY.length; i++){if(p.termEquals(GEOMETRY[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown physical-geometry status");
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private String path(final String file){return new File(output, file).getPath();}

	/*--------------------------------------------------------------*/
	/*----------------          Data Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Call{
		String dataset, assembly, reference;
		int length, fixed, minimum, maximum, centerBad;
		long support;
		final int[] raw=new int[2], truth=new int[2], center=new int[2];
		final long[] proposal=new long[2];
		/** Read only when proposal differs from raw; that existing state proves inference. */
		final float[] score=new float[2];
		boolean hasGateScore(final int end){return proposal[end]!=raw[end];}
		final float[][] heads=new float[2][];
		final boolean[] missing=new boolean[2];
	}
	private static final class Choice{
		Choice(final String policy_, final int end_, final String side_, final float cutoff_, final boolean never_, final Stats stats_, final long target_){
			policy=policy_;
			end=end_;
			side=side_;
			require(Float.isFinite(cutoff_), "Choice storage cannot contain a nonfinite cutoff");
			cutoff=cutoff_;
			never=never_;
			stats=stats_;
			target=target_;
		}
		final String policy, side;
		final int end;
		final float cutoff;
		final boolean never;
		final Stats stats;
		final long target;
	}
	private static final class Stats{
		void add(final Call call, final int[] result, final int end, final int bad){
			final long raw=(long)call.raw[end]-call.truth[end];
			final long error=(long)result[end]-call.truth[end];
			final long original=Math.abs(raw), after=Math.abs(error);
			n++;
			if(after<original){better++;}else if(after>original){worse++;}else{tied++;}
			if(error!=raw){moved++;}
			if(after==0){exact++;}
			if(after<=3){within3++;}
			if(after<=10){within10++;}
			final long other=Math.abs((long)result[1-end]-call.truth[1-end]);
			if(after==0 && other==0){bothExact++;}
			if(after<=3 && other<=3){both3++;}
			if(after<=10 && other<=10){both10++;}
			rawSquared+=(double)raw*raw;
			squared+=(double)error*error;
			if(original==0 && after>0){exactHarmed++;}
			if(original<=10){
				near++;
				nearRawSquared+=(double)raw*raw;
				nearSquared+=(double)error*error;
				if(after>original){nearHarmed++;}
			}else if(after<original){distantImproved++;}else if(after>original){distantHarmed++;}
			if((call.fixed&(1<<end))!=0){fixed++;}
			if(call.support<10){unsupported++;}
			if(call.missing[end]){missing++;}
			if(call.centerBad!=0){centerHeld++;}
			invalid[bad]++;
			assert(better+worse+tied==n) : "Every original matched call must remain in each paired outcome denominator";
		}
		double rms(){return Math.sqrt(squared/n);}
		double rawRms(){return Math.sqrt(rawSquared/n);}
		long n, better, worse, tied, moved, exact, within3, within10, bothExact, both3, both10;
		long exactHarmed, near, nearHarmed, distantImproved, distantHarmed, fixed, unsupported, missing, centerHeld;
		double squared, rawSquared, nearSquared, nearRawSquared;
		final long[] invalid=new long[5];
	}
	private String five, three, output, choices, member;
	private int seed=-1;
	private long expected=-1;
	private boolean select=true;
	private final String[] gateNetworks=new String[2];
	private final boolean[] legacyCache=new boolean[2];
	private final boolean[] legacyMaximum=new boolean[2];
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private final ArrayList<Choice> selected=new ArrayList<Choice>();
	private final ArrayList<Choice> curvePoints=new ArrayList<Choice>();
	private static final int ALL_PROPOSALS=6;
	private static final String CHOICE_MAGIC="#NCRNA_CENTERED_GATE_CHOICES_V2";
	private static final String LEGACY_CHOICE_MAGIC="#NCRNA_CENTERED_GATE_CHOICES_V1";
	private static final String LEGACY_CHOICE_HEADER="policy\tend\tside\tavailable\tcutoff_bits\tcutoff\tmonitor_moved\ttarget_moved\tmonitor_better\tmonitor_worse";
	/** Availability owns the monitoring cells; numeric cutoff cells exist only in THRESHOLD mode. */
	private static final String CHOICE_HEADER="policy\tend\tside\tavailable\tcutoff_bits_if_threshold\tcutoff_if_threshold\tmonitor_moved_if_available\ttarget_moved_if_rate\tmonitor_better_if_available\tmonitor_worse_if_available\tcutoff_mode";
	public static final String LEGACY_RESULT_HEADER="member\tseed\tpanel\tpolicy\tside5\tside3\tcutoff5\tcutoff3\tscope\tend\tmatched\texact\twithin3\twithin10\tboth_exact\tboth_within3\tboth_within10\traw_rms_bp\trms_bp\trms_minus_raw_bp\tbetter\tworse\ttied\tmoved\texact_harmed\traw_near10_n\traw_near10_rms_bp\tnear10_rms_bp\tnear10_harmed\tdistant_improved\tdistant_harmed\tfixed\tunsupported\tmissing_features\tcenter_held\tinvalid_joint\tcrossed\toutside_contig\tbelow_minimum\tabove_maximum\tbetter_fraction\traw_sum_squared_bp2\tsum_squared_bp2\tnear_raw_sum_squared_bp2\tnear_sum_squared_bp2\tqualification";
	public static final String RESULT_HEADER="member\tseed\tpanel\tpolicy\tside5\tside3\tcutoff5_if_threshold\tcutoff3_if_threshold\tscope\tend\tmatched\texact\twithin3\twithin10\tboth_exact_if_joint\tboth_within3_if_joint\tboth_within10_if_joint\traw_rms_bp\trms_bp\trms_minus_raw_bp\tbetter\tworse\ttied\tmoved\texact_harmed\traw_near10_n\traw_near10_rms_bp_if_nonempty\tnear10_rms_bp_if_nonempty\tnear10_harmed\tdistant_improved\tdistant_harmed\tfixed\tunsupported\tmissing_features\tcenter_held\tinvalid_joint\tcrossed\toutside_contig\tbelow_minimum\tabove_maximum\tbetter_fraction_if_changed\traw_sum_squared_bp2\tsum_squared_bp2\tnear_raw_sum_squared_bp2\tnear_sum_squared_bp2\tqualification\tcutoff5_mode\tcutoff3_mode";
	private static final String[] CONTROL_NAMES={"RAW", "CENTER_ONLY", "CENTER_MARGIN0", "CENTER_MARGIN01", "ALL_WHOLE_PROPOSALS", "IGNORE_HEAD_DIAGNOSTIC"};
	private static final String[] ENDS={"5prime", "3prime"};
	private static final String[] GEOMETRY={"VALID", "CROSSED", "OUTSIDE_CONTIG", "BELOW_MINIMUM", "ABOVE_MAXIMUM"};
}
