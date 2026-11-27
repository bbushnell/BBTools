package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import map.IntHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Compares saved cobalamin endpoints and sequence-bound CM input padding.
 * Uses the A_SampleByteFile argument/I/O lifecycle and native retained-key maps.
 * References are the grader's recorded assignments; nothing is rematched.
 * Alias-group parity zero selects each model/end policy; parity one judges it.
 * Explicit fitall pools the same eligible saved evidence for a final export;
 * its endpoint estimates are in-sample and have no independent judging half.
 * Spliced endpoints retain their recorded reference and path membership.
 * They are estimates, not outputs of a coupled covariance search or caller.
 * @author Brian Bushnell, Raiden
 */
public final class AnalyzeCobalaminInputPadding{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AnalyzeCobalaminInputPadding.class, false);
		try{
			new AnalyzeCobalaminInputPadding(pp.args).process();
			pp.outstream.println("COBALAMIN_INPUT_PADDING_ANALYSIS_COMPLETE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private AnalyzeCobalaminInputPadding(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("names")){names=list(b);}
			else if(a.equals("traces")){traces=list(b);}
			else if(a.equals("profiles")){profiles=list(b);}
			else if(a.equals("padding")){padding=list(b);}
			else if(a.equals("library")){libraryPath=b;}
			else if(a.equals("bindings")){bindings=b;}
			else if(a.equals("audit")){audit=b;}
			else if(a.equals("basepadding")){basePaddingPath=b;}
			else if(a.equals("minselectioncalls")){minSelectionCalls=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("minselectiongroups")){minSelectionGroups=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("exportselected")){exportSelected=Parse.parseBoolean(b);}
			else if(a.equals("fitall")){fitAll=Parse.parseBoolean(b);}
			else if(a.equals("expectedinputs")){expectedInputs=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("cmcpu")){cpuValues=list(b);}
			else if(a.equals("expectedharm5")){harmValues=list(b);}
			else{require(parser.parse(arg, a, b), "Unknown input-padding analysis argument: "+arg);}
		}
		require(parser.in1!=null && parser.in2==null && parser.out1!=null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Supply complete normalized witness prefixes and fresh out=prefix without filtering");
		grades=list(parser.in1);
		out=parser.out1;
		require(grades.length>=2 && grades.length<=8 && names!=null && names.length==grades.length
			&& traces!=null && traces.length==grades.length && profiles!=null && profiles.length==grades.length
			&& padding!=null && padding.length==grades.length-1 && expectedInputs>0,
			"One baseline and each candidate require distinct names, witnesses, traces and profiles; candidates also require padding");
		require(libraryPath!=null && bindings!=null && audit!=null && Tools.testInputFiles(false, true, libraryPath, bindings, audit),
			"The exact consensus, ordered panel and resolved-alias audit must be readable");
		for(int side=0; side<grades.length; side++){
			require(names[side].matches("[A-Za-z0-9_]+"), "Arm names must be printable identifiers");
			for(int j=0; j<side; j++){require(!names[j].equals(names[side]), "Each saved arm needs a distinct name");}
			require(Tools.testInputFiles(false, true, grades[side]+".witnesses.tsv", grades[side]+".genomes.tsv", traces[side], profiles[side]),
				"Every arm requires its accepted witnesses, version-four trace and independent CM invocation profile");
		}
		require(cpuValues==null || cpuValues.length==padding.length,
			"Optional cmcpu= has one observed CM-search CPU value per candidate, excluding anchors");
		require(harmValues==null || harmValues.length==padding.length,
			"Optional expectedharm5= has one accepted pairwise count per candidate");
		RrnaResourceIO.fresh(out, ".coverage.tsv", ".model_ends.tsv", ".harms.tsv", ".estimates.tsv", ".costs.tsv");
		require(minSelectionCalls>0 && minSelectionGroups>0, "Model/end selection needs positive observation and distinct-source-group minima");
		require(!exportSelected || basePaddingPath!=null && grades.length==2,
			"Export selected padding separately against one explicit base, never an accidental multi-arm intersection");
		require(!fitAll || exportSelected,
			"fitall is a final selected-padding export on all eligible saved pairs, with no independent judging half");
		if(basePaddingPath!=null){require(Tools.testInputFiles(false, true, basePaddingPath), "The explicit baseline padding resource must be readable");}
		if(exportSelected){RrnaResourceIO.fresh(out, ".selected."+names[1]+".tsv");}
		Shared.setThreads(1);
	}

	private void process(){
		readModels();
		if(basePaddingPath!=null){
			final NcrnaCmInputPadding.Data base=NcrnaCmInputPadding.read(basePaddingPath, library.names, library.bases);
			require(base.fixed, "A declared uniform base must be a fixed resource, not another selected or percentile policy");
			base5=base.pad5(0);
			base3=base.pad3(0);
		}
		readPanel();
		readAudit();
		arms=new Arm[grades.length];
		for(int side=0; side<arms.length; side++){
			final Arm arm=new Arm(names[side], models.length, panel.size());
			arms[side]=arm;
			if(side>0){
				arm.padding=NcrnaCmInputPadding.read(padding[side-1], library.names, library.bases);
				if(exportSelected){
					require(!arm.padding.fixed && !arm.padding.selected && "P98_PLUS10".equals(arm.padding.rule),
						"Selected-base export combines the original supported P98_PLUS10 fit, not a fixed or previously selected mask");
				}
				if(cpuValues!=null){
					arm.cpu=Double.parseDouble(cpuValues[side-1]);
					require(Double.isFinite(arm.cpu) && arm.cpu>=0, "Observed search CPU must be finite nonnegative seconds");
				}
			}
			readTrace(arm, traces[side]);
			checkProfile(arm, profiles[side]);
			readWitnesses(arm, grades[side]+".witnesses.tsv");
			checkGenomes(arm, grades[side]+".genomes.tsv");
		}
		join();
		write();
		if(exportSelected){
			final Arm arm=arms[1];
			final boolean[] selected5=new boolean[models.length], selected3=new boolean[models.length];
			for(int m=0; m<models.length; m++){
				selected5[m]=choose(arm, m, 0);
				selected3[m]=choose(arm, m, 1);
			}
			NcrnaCmInputPadding.writeSelected(out+".selected."+arm.name+".tsv", base5, base3, arm.padding, selected5, selected3);
		}
	}

	private void readModels(){
		library=NcrnaSupport.loadLibrary(libraryPath);
		models=new Model[library.names.length];
		for(int i=0; i<models.length; i++){
			final Model model=new Model(i, library.names[i]);
			nameKey(model.name);
			require(!byModel.contains(probe), "Consensus names uniquely identify the sequence-bound padding axis");
			models[i]=model;
			byModel.put(probe.copy(), model);
		}
	}
	private void readPanel(){
		try(Input input=new Input(bindings)){
			RrnaResourceIO.header(input.reader, PANEL, bindings);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5 && p.parseInt(0)==panel.size(), "Panel ordinals must be complete and in their accepted order");
				probe.set(p, 2);
				require(!byInput.contains(probe), "Panel input IDs must be unique");
				final Source source=new Source(panel.size(), p.parseString(1), p.parseString(2), p.parseString(3));
				byInput.put(probe.copy(), source);
				panel.add(source);
			}
		}
		require(panel.size()==expectedInputs, "The complete panel includes every declared negative and positive genome");
	}
	private void readAudit(){
		try(Input input=new Input(audit)){
			RrnaResourceIO.header(input.reader, AUDIT, audit);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==8 && p.parseInt(4)>=0, "Alias audit rows need the declared nonnegative resolved group");
				probe.set(p, 1);
				final Source source=byInput.get(probe);
				if(source==null){continue;}
				require(source.group<0 && p.termEquals(source.assembly, 0),
					"Each panel input binds exactly one audited assembly and resolved group before policy selection");
				source.group=p.parseInt(4);
			}
		}
		for(final Source source:panel){require(source.group>=0, "Panel input lacks an audited policy-selection group: "+source.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Saved Evidence        ----------------*/
	/*--------------------------------------------------------------*/

	private void readTrace(final Arm arm, final String path){
		final ObjectMap<Target, Call> invocations=new ObjectMap<Target, Call>(Target.class, Call.class);
		final ObjectMap<Target, Terminal> terminals=new ObjectMap<Target, Terminal>(Target.class, Terminal.class);
		final boolean[] seen=new boolean[panel.size()];
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "#ncrna_cm_proposals\t4", path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				final Source source=ordinal(p.parseInt(1));
				if(p.termEquals("I", 0)){
					require(p.terms()==3 && !seen[source.ordinal] && p.termEquals(source.fasta, 2),
						"Each trace input must bind its exact panel FASTA once");
					seen[source.ordinal]=true;
					continue;
				}
				require(seen[source.ordinal], "The input binding must precede trace evidence");
				if(p.termEquals("C", 0)){
					require(p.terms()==37 && p.termEquals("cobalamin_rf00174", 7),
						"Input-padding analysis requires the accepted main-member version-four trace scope");
					traceKey();
					require(!invocations.contains(probe), "A saved request identity occurs once");
					final Call call=new Call(source, model(8), p.parseString(5), p.parseInt(6),
						p.parseInt(9), p.parseInt(10), p.parseInt(11), bit(32));
					if(call.selected){
						require(p.parseInt(23)>=0 && p.parseInt(24)>=p.parseInt(23) && p.parseInt(24)<call.length,
							"Selected effective sense bounds must remain inside their actual source contig");
						call.finalStart=p.parseInt(25)+1;
						call.finalStop=p.parseInt(26)+1;
						require(call.finalStart==genomicStart(call, p.parseInt(23), p.parseInt(24))
							&& call.finalStop==genomicStop(call, p.parseInt(23), p.parseInt(24)),
							"A selected final interval must equal the effective sense interval converted once to genomic coordinates");
						geometryKey(source, call.seqid, call.strand, call.model, call.finalStart, call.finalStop);
						require(!arm.selected.contains(probe), "An emitted witness needs a unique exact selected trace geometry");
						arm.selected.put(probe.copy(), call);
					}
					traceKey();
					invocations.put(probe.copy(), call);
					arm.calls.add(call);
					terminalKey(false);
					Terminal terminal=terminals.get(probe);
					if(terminal==null){
						terminal=new Terminal(call.seqid, call.length);
						terminals.put(probe.copy(), terminal);
					}
					require(!terminal.closed && terminal.seqid.equals(call.seqid) && terminal.length==call.length,
						"Candidates cannot follow their terminal or change a native worker/read contig identity");
					call.terminal=terminal;
					terminal.calls++;
					if(bit(31)){terminal.returned++;}
					if(call.selected){terminal.selected++;}
					if(!bit(14)){terminal.unknown++;}
				}else if(p.termEquals("M", 0)){
					require(p.terms()==21, "The saved CM measurement schema has changed");
					traceKey();
					final Call call=invocations.get(probe);
					require(call!=null && !call.measured && p.parseInt(4)==call.start && p.parseInt(5)==call.stop
						&& p.parseInt(6)==call.stop-call.start+1,
						"A CM measurement binds one original post-trim interval, without transferring ordinals between arms");
					call.measured=true;
					call.decision=p.parseString(20);
					call.cmLength=p.parseInt(11);
					call.attempted=!p.termEquals("", 12);
					require(call.cmLength>0 && p.parseInt(10)==library.bases[call.model].length,
						"Saved CM and consensus lengths must bind the exact serving model");
					arm.measured++;
					if(call.attempted){arm.attempted++;}
				}else{
					require(p.termEquals("T", 0) && p.terms()==11, "Unknown trace record or changed terminal schema");
					terminalKey(true);
					Terminal terminal=terminals.get(probe);
					if(terminal==null){
						terminal=new Terminal(p.parseString(4), p.parseInt(5));
						terminals.put(probe.copy(), terminal);
					}
					require(!terminal.closed && terminal.calls==p.parseInt(6) && terminal.returned==p.parseInt(7)
						&& terminal.selected==p.parseInt(8) && terminal.unknown==p.parseInt(9) && p.parseInt(10)==0
						&& p.termEquals(terminal.seqid, 4) && p.parseInt(5)==terminal.length,
						"Each actual contig terminal conserves its candidates and selection states");
					terminal.closed=true;
				}
			}
		}
		for(final boolean present:seen){require(present, "The trace omitted a declared panel input");}
		// Each candidate's terminal was checked at T; an absent T is caught here.
		for(final Call call:arm.calls){
			require(call.terminal!=null && call.terminal.closed,
				"Every candidate must belong to a closed contig terminal");
		}
	}

	private void checkProfile(final Arm arm, final String path){
		try(Input input=new Input(path)){
			byte[] header;
			do{header=input.reader.nextLine();}while(header!=null && header.length>0 && header[0]=='#');
			require(header!=null, "The independent invocation profile requires a header");
			p.set(header);
			require(p.terms()>=2 && p.termEquals("member", 0) && p.termEquals("calls", 1),
				"The profile must identify actual CM invocations");
			final byte[] row=input.reader.nextLine();
			require(row!=null, "The main member requires one independent invocation count");
			p.set(row);
			require((p.termEquals("RF00174", 0) || p.termEquals("cobalamin_rf00174", 0)) && p.parseLong(1)==arm.measured
				&& input.reader.nextLine()==null,
				"Actual main-member M rows must reconcile with the independent profile; omitted whole contigs cannot pass");
		}
	}

	private void readWitnesses(final Arm arm, final String path){
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, WITNESSES, path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==12, "Normalized endpoint witnesses require their declared twelve-column schema");
				probe.set(p, 0);
				final Source source=byInput.get(probe);
				require(source!=null, "Every witness must belong to the bound panel");
				final Witness witness=new Witness(source, p.parseString(1), p.parseString(2),
					CompareCobalaminBenchmarkGrades.modelId(p.parseString(3)), p.parseString(4),
					p.parseInt(5), p.parseInt(6), p.parseInt(7), p.parseInt(8), strand(9),
					p.parseLong(10), p.parseLong(11));
				require(witness.start>0 && witness.stop>=witness.start && witness.refStart>0 && witness.refStop>=witness.refStart
					&& witness.d5==(witness.strand==0 ? (long)witness.start-witness.refStart : (long)witness.refStop-witness.stop)
					&& witness.d3==(witness.strand==0 ? (long)witness.stop-witness.refStop : (long)witness.refStart-witness.start),
					"Witness signed offsets must agree with their recorded one-based genomic bounds and strand");
				require(witness.member.equals("RF00174") || witness.member.equals("RF01482") || witness.member.equals("RF01689"),
					"Every retained witness must identify one known cobalamin member");
				if(witness.main()){
					nameKey(witness.modelName);
					final Model model=byModel.get(probe);
					require(model!=null, "A main-member witness must bind the exact sequence-bound model axis");
					witness.model=model.ordinal;
					geometryKey(source, witness.seqid, witness.strand, witness.model, witness.start, witness.stop);
					witness.call=arm.selected.get(probe);
					require(witness.call!=null && !witness.call.linked, "A recorded TP needs one unused exact selected trace link");
					witness.call.linked=true;
					arm.mainWitnesses++;
				}
				referenceKey(source, witness.reference);
				require(!arm.references.contains(probe), "A graded reference is consumed at most once in an arm");
				arm.references.put(probe.copy(), witness);
				arm.witnesses.add(witness);
				arm.matchedByInput[source.ordinal]++;
			}
		}
	}
	private void checkGenomes(final Arm arm, final String path){
		final boolean[] seen=new boolean[panel.size()];
		int n=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "input_id\tcalls\tmatched\tmasked\tunmatched", path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				probe.set(p, 0);
				final Source source=byInput.get(probe);
				require(p.terms()==5 && source!=null && !seen[source.ordinal]
					&& p.parseLong(1)>=0 && p.parseLong(3)>=0 && p.parseLong(4)>=0
					&& p.parseLong(1)==p.parseLong(2)+p.parseLong(3)+p.parseLong(4)
					&& p.parseLong(2)==arm.matchedByInput[source.ordinal],
					"Each complete genome table must conserve its accepted outcomes and every recorded TP witness, including negative inputs");
				seen[source.ordinal]=true;
				n++;
			}
		}
		require(n==panel.size(), "The accepted genome table must include the complete diagnostic panel");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixed Comparisons      ----------------*/
	/*--------------------------------------------------------------*/

	private void join(){
		for(final Witness before:arms[0].witnesses){
			final Witness[] pair=new Witness[arms.length];
			pair[0]=before;
			boolean complete=true;
			for(int side=1; side<arms.length; side++){
				referenceKey(before.source, before.reference);
				final Witness after=arms[side].references.get(probe);
				pair[side]=after;
				if(after==null){
					complete=false;
					continue;
				}
				require(before.source==after.source && before.seqid.equals(after.seqid) && before.strand==after.strand
					&& before.refStart==after.refStart && before.refStop==after.refStop,
					"A paired reference must keep its input, genomic identity and literal reference bounds across arms");
				require(before.main()==after.main(), "A main/minor member switch needs a separate policy design, not an implicit padding splice");
				if(before.main()){
					require(before.call.length==after.call.length, "The same frozen input and contig cannot change length between the saved caller arms");
				}
				arms[side].paired++;
				if(nearHarm(before.d5, after.d5)){arms[side].harm5++;}
				if(!before.main()){
					require(before.d5==after.d5 && before.d3==after.d3,
						"These main-member-only splices require the minor members' measured endpoints to remain unchanged");
				}
			}
			if(complete){common.add(pair);}
		}
		require(!common.isEmpty(), "No fixed all-arm reference intersection exists for judging endpoint policies");
		for(int side=1; side<arms.length; side++){
			final Arm arm=arms[side];
			if(harmValues!=null){require(arm.harm5==Long.parseLong(harmValues[side-1]), "Pairwise prior-close5 harms disagree with the accepted comparison");}
			for(final Witness[] pair:common){
				final Witness before=pair[0], after=pair[side];
				if(!after.main()){continue;}
				for(int end=0; end<2; end++){
					arm.ends[after.model][end].add(before, after, end, pad(arm, after.model, end), base(end), fitAll);
				}
			}
		}
	}

	private void write(){
		try(Output coverage=new Output(out+".coverage.tsv"); Output ends=new Output(out+".model_ends.tsv");
			Output harms=new Output(out+".harms.tsv"); Output estimates=new Output(out+".estimates.tsv");
			Output costs=new Output(out+".costs.tsv")){
			coverage.writer.println("arm\twitnesses\tmain_witnesses\tpaired_with_baseline\tpairwise_prior_close5_harms\tall_arm_common"
				+(fitAll ? "\tfitting_common\tindependent_judgment_common" : "\tselection_common\tjudgment_common"));
			ends.writer.println(schema("arm\tcandidate_model\tend\tconfigured_pad\tconfigured_relation20\tselection_n\tselection_baseline_sse\tselection_candidate_sse\tuse_candidate_pad\tcommon_n\tcommon_baseline_sse\tcommon_candidate_sse\tbetter\tworse\ttied\tprior_close_harms\tcrossing3\tmodel_switches\tclipped_pad_sum\tbaseline_clipped_pad_sum\tbetter_clipped_less20\tbetter_clipped_equal20\tbetter_clipped_more20\tworse_clipped_less20\tworse_clipped_equal20\tworse_clipped_more20\tbetter_clipped_delta_sum\tworse_clipped_delta_sum")+(basePaddingPath==null ? "" : "\tselection_groups"));
			harms.writer.println(schema("arm\tinput_id\treference_id\tseqid\tstrand\tend\tbaseline_model\tcandidate_model\tbaseline_error\tcandidate_error\tmovement\tcrossing3\tbaseline_configured_pad\tcandidate_configured_pad\tbaseline_clipped_pad\tcandidate_clipped_pad\tcandidate_same_interval_pad20\tresolved_group\tbaseline_cm_decision\tcandidate_cm_decision\tbaseline_original_start0\tbaseline_original_stop0\tcandidate_original_start0\tcandidate_original_stop0\tcontig_length\tqualification"));
			estimates.writer.println("arm\tpolicy\tcohort"+(fitAll ? "\tin_sample_n" : "\tjudgment_n")
				+"\tboth_exact\texact5\texact3\twithin3_5\twithin3_3\tprior_close_harms5\tprior_close_harms3\tcrossing3_5\tcrossing3_3\trms5\trms3\tqualification");
			costs.writer.println("arm\tpolicy\tattempted_requests"+(fitAll ? "\tin_sample_attempted_requests" : "\tjudgment_attempted_requests")
				+"\tfull_dense_work_proxy\tpolicy_dense_work_proxy\tproxy_ratio\tobserved_full_cm_cpu_seconds\tprojected_cm_cpu_seconds\tqualification");
			long choosing=0;
			for(final Witness[] pair:common){if(fitAll || (pair[0].source.group&1)==0){choosing++;}}
			for(int side=0; side<arms.length; side++){
				final Arm arm=arms[side];
				coverage.writer.print(buffer.clear().append(arm.name).tab().append(arm.witnesses.size()).tab().append(arm.mainWitnesses)
					.tab().append(side==0 ? arm.witnesses.size() : arm.paired).tab().append(arm.harm5)
					.tab().append(common.size()).tab().append(choosing).tab().append(common.size()-choosing).nl());
				if(side==0){continue;}
				writeEnds(arm, ends);
				for(final Witness before:arms[0].witnesses){
					referenceKey(before.source, before.reference);
					final Witness after=arm.references.get(probe);
					if(after==null || !after.main()){continue;}
					for(int end=0; end<2; end++){
						if(nearHarm(before.error(end), after.error(end))){writeHarm(arm, before, after, end, harms);}
					}
				}
				for(int policy=0; policy<4; policy++){
					final Metric[] metrics=new Metric[COHORTS.length];
					for(int cohort=0; cohort<metrics.length; cohort++){metrics[cohort]=new Metric();}
					for(final Witness[] pair:common){
						final Witness before=pair[0], after=pair[side];
						if(!fitAll && (before.source.group&1)!=1){continue;}
						final long d5=use(arm, after, policy, 0) ? after.d5 : before.d5;
						final long d3=use(arm, after, policy, 1) ? after.d3 : before.d3;
						metrics[0].add(before, d5, d3);
						metrics[before.modelName.equals(after.modelName) ? 1 : 2].add(before, d5, d3);
						final boolean configured=after.main() && (pad(arm, after.model, 0)!=base5 || pad(arm, after.model, 1)!=base3);
						final boolean actual=after.main() && (clipped(after.call, 0, pad(arm, after.model, 0))!=clipped(after.call, 0, base5)
							|| clipped(after.call, 1, pad(arm, after.model, 1))!=clipped(after.call, 1, base3));
						metrics[configured ? 4 : 3].add(before, d5, d3);
						metrics[actual ? 6 : 5].add(before, d5, d3);
					}
					for(int cohort=0; cohort<metrics.length; cohort++){writeMetric(arm, policy, cohort, metrics[cohort], estimates);}
					writeCost(arm, policy, costs);
				}
			}
		}
	}

	private void writeEnds(final Arm arm, final Output output){
		for(final Model model:models){
			for(int end=0; end<2; end++){
				final EndStat s=arm.ends[model.ordinal][end];
				final int padding=pad(arm, model.ordinal, end);
				buffer.clear().append(arm.name).tab().append(model.name).tab().append(END_NAMES[end])
					.tab().append(padding).tab().append(padding<base(end) ? "LESS" : padding==base(end) ? "EQUAL" : "MORE")
					.tab().append(s.selectionN).tab().append(s.oldSquares)
					.tab().append(s.newSquares).tab().append(choose(arm, model.ordinal, end)).tab().append(s.n)
					.tab().append(s.allOldSquares).tab().append(s.allNewSquares).tab().append(s.better)
					.tab().append(s.worse).tab().append(s.n-s.better-s.worse).tab().append(s.harms).tab().append(s.crossings)
					.tab().append(s.switches).tab().append(s.clipped).tab().append(s.baselineClipped)
					.tab().append(s.betterLess).tab().append(s.betterEqual).tab().append(s.betterMore)
					.tab().append(s.worseLess).tab().append(s.worseEqual).tab().append(s.worseMore)
					.tab().append(s.betterDelta).tab().append(s.worseDelta);
				if(basePaddingPath!=null){buffer.tab().append(s.selectionGroups.size());}
				output.writer.print(buffer.nl());
			}
		}
	}
	private void writeHarm(final Arm arm, final Witness before, final Witness after, final int end, final Output output){
		final int configured=pad(arm, after.model, end);
		output.writer.print(buffer.clear().append(arm.name).tab().append(before.source.id).tab().append(before.reference)
			.tab().append(before.seqid).tab().append(before.strand==0 ? '+' : '-').tab().append(END_NAMES[end])
			.tab().append(before.modelName).tab().append(after.modelName).tab().append(before.error(end)).tab().append(after.error(end))
			.tab().append(after.error(end)-before.error(end)).tab().append(Math.abs(after.error(end))>3).tab().append(base(end))
			.tab().append(configured).tab().append(clipped(before.call, end, base(end))).tab().append(clipped(after.call, end, configured))
			.tab().append(clipped(after.call, end, base(end))).tab().append(before.source.group)
			.tab().append(before.call.decision).tab().append(after.call.decision)
			.tab().append(before.call.start).tab().append(before.call.stop).tab().append(after.call.start).tab().append(after.call.stop)
			.tab().append(after.call.length)
			.append("\tPAIRED_RECORDED_REFERENCE_NO_REMATCHING\n"));
	}
	private void writeMetric(final Arm arm, final int policy, final int cohort, final Metric s, final Output output){
		buffer.clear().append(arm.name).tab().append(label(fitAll && policy==3 ? "ALL_ELIGIBLE_MODEL_END" : POLICIES[policy]))
			.tab().append(label(COHORTS[cohort])).tab().append(s.n)
			.tab().append(s.both).tab().append(s.exact5).tab().append(s.exact3).tab().append(s.within5).tab().append(s.within3)
			.tab().append(s.harm5).tab().append(s.harm3).tab().append(s.cross5).tab().append(s.cross3).tab();
		if(s.n>0){buffer.append(Math.sqrt(s.square5/s.n), 8);}
		buffer.tab();
		if(s.n>0){buffer.append(Math.sqrt(s.square3/s.n), 8);}
		output.writer.print(buffer.append(fitAll ? "\tALL_ELIGIBLE_FIXED_REFERENCE_ENDPOINT_SPLICE_IN_SAMPLE_NOT_COUPLED_CALLER\n"
			: "\tGROUP_PARITY1_FIXED_REFERENCE_ENDPOINT_SPLICE_NOT_COUPLED_CALLER\n"));
	}
	private void writeCost(final Arm arm, final int policy, final Output output){
		double full=0, projected=0;
		long n=0, judgment=0;
		for(final Call call:arm.calls){
			if(!call.attempted){continue;}
			n++;
			if(fitAll || (call.source.group&1)==1){judgment++;}
			full+=work(call, pad(arm, call.model, 0), pad(arm, call.model, 1));
			final int pad5=use(arm, call.model, policy, 0) ? pad(arm, call.model, 0) : base5;
			final int pad3=use(arm, call.model, policy, 1) ? pad(arm, call.model, 1) : base3;
			projected+=work(call, pad5, pad3);
		}
		require(n==arm.attempted && Double.isFinite(full) && full>0 && Double.isFinite(projected),
			"Cost ratios use the same nonempty observed attempted-request population, not fabricated searches");
		final double ratio=projected/full;
		buffer.clear().append(arm.name).tab().append(label(fitAll && policy==3 ? "ALL_ELIGIBLE_MODEL_END" : POLICIES[policy])).tab().append(n).tab().append(judgment)
			.tab().append(full, 0).tab().append(projected, 0).tab().append(ratio, 8).tab();
		if(Double.isFinite(arm.cpu)){buffer.append(arm.cpu, 6);}
		buffer.tab();
		if(Double.isFinite(arm.cpu)){buffer.append(arm.cpu*ratio, 6);}
		output.writer.print(buffer.append("\tCM_LENGTH_TIMES_TRIANGULAR_INPUT_SPAN_UNCERTAIN_SPARSE_RETRY_PROXY_ANCHOR_EXCLUDED\n"));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private boolean use(final Arm arm, final Witness witness, final int policy, final int end){
		return witness.main() && use(arm, witness.model, policy, end);
	}
	private boolean use(final Arm arm, final int model, final int policy, final int end){
		assert(policy>=0 && policy<4 && end>=0 && end<2) : "Policy and sense-end indexes address the four declared endpoint splices";
		return policy==1 || policy==2 && end==1 || policy==3 && choose(arm, model, end);
	}
	private boolean choose(final Arm arm, final int model, final int end){
		final EndStat s=arm.ends[model][end];
		if(s.selectionN<minSelectionCalls || s.selectionGroups.size()<minSelectionGroups || s.newSquares>=s.oldSquares){return false;}
		if(basePaddingPath==null && !exportSelected){return true;}
		return !arm.padding.fixed && !arm.padding.selected
			&& arm.padding.finiteCalls(model)>=arm.padding.minCalls && arm.padding.sourceGroups(model)>=arm.padding.minGroups
			&& pad(arm, model, end)!=base(end);
	}
	private int base(final int end){
		assert(end==0 || end==1) : "Only the two sense-oriented ends have declared base flanks";
		return end==0 ? base5 : base3;
	}
	private String schema(final String text){return basePaddingPath==null ? text : text.replace("20", "_base");}
	private String label(final String text){return basePaddingPath==null ? text : text.replace("20", "_BASE");}
	private static double work(final Call call, final int pad5, final int pad3){
		final long span=(long)call.stop-call.start+1+clipped(call, 0, pad5)+clipped(call, 1, pad3);
		assert(span>0 && call.cmLength>0) : "Validated original spans and CM lengths define the nominal dense work proxy";
		return call.cmLength*(double)span*(span+1)*0.5;
	}
	private static int clipped(final Call call, final int end, final int padding){
		assert(call!=null && padding>=0) : "Only a linked main-member request and nonnegative padding define a clipped flank";
		return Math.min(padding, end==0 ? call.start : call.length-1-call.stop);
	}
	private int pad(final Arm arm, final int model, final int end){return end==0 ? arm.padding.pad5(model) : arm.padding.pad3(model);}
	private static boolean nearHarm(final long before, final long after){return Math.abs(before)<=3 && Math.abs(after)>Math.abs(before);}
	private static int genomicStart(final Call call, final int start, final int stop){return (call.strand==0 ? start : call.length-1-stop)+1;}
	private static int genomicStop(final Call call, final int start, final int stop){return (call.strand==0 ? stop : call.length-1-start)+1;}
	private Source ordinal(final int ordinal){
		require(ordinal>=0 && ordinal<panel.size(), "Trace ordinals must address the exact accepted panel");
		return panel.get(ordinal);
	}
	private boolean bit(final int field){
		require(p.termEquals("0", field) || p.termEquals("1", field), "Trace flags must be explicit zero or one");
		return p.termEquals("1", field);
	}
	private int strand(final int field){
		require(p.termEquals("+", field) || p.termEquals("-", field), "Normalized witnesses require an explicit strand");
		return p.termEquals("+", field) ? 0 : 1;
	}
	private int model(final int field){
		probe.set(p, field);
		final Model model=byModel.get(probe);
		require(model!=null, "Trace model name must occur on the exact consensus/padding axis");
		return model.ordinal;
	}
	private void traceKey(){
		key.clear().append(p.parseInt(1)).tab().append(p.parseInt(2)).tab().append(p.parseLong(3));
		probe.set(key.array, 0, key.length());
	}
	private void terminalKey(final boolean terminal){
		key.clear().append(p.parseInt(1)).tab().append(p.parseInt(2)).tab().append(p.parseLong(terminal ? 3 : 4));
		probe.set(key.array, 0, key.length());
	}
	private void nameKey(final String name){
		key.clear().append(name);
		probe.set(key.array, 0, key.length());
	}
	private void referenceKey(final Source source, final String reference){
		key.clear().append(source.ordinal).tab().append(reference);
		probe.set(key.array, 0, key.length());
	}
	private void geometryKey(final Source source, final String seqid, final int strand, final int model, final int start, final int stop){
		key.clear().append(source.ordinal).tab().append(seqid).tab().append(strand).tab().append(model).tab().append(start).tab().append(stop);
		probe.set(key.array, 0, key.length());
	}
	private static String[] list(final String value){
		require(value!=null && !value.isEmpty(), "List-valued options require explicit nonempty values");
		final String[] values=value.split(",", -1);
		for(final String item:values){require(!item.isEmpty(), "Empty list elements cannot silently drop an arm");}
		return values;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained State       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final int ordinal_, final String assembly_, final String id_, final String fasta_){
			ordinal=ordinal_;
			assembly=assembly_;
			id=id_;
			fasta=fasta_;
		}
		final int ordinal;
		final String assembly, id, fasta;
		int group=-1;
	}
	private static final class Model{
		Model(final int ordinal_, final String name_){
			ordinal=ordinal_;
			name=name_;
		}
		final int ordinal;
		final String name;
	}
	private static final class Terminal{
		Terminal(final String seqid_, final int length_){
			seqid=seqid_;
			length=length_;
		}
		final String seqid;
		final int length;
		long calls, returned, selected, unknown;
		boolean closed;
	}
	private static final class Call{
		Call(final Source source_, final int model_, final String seqid_, final int strand_, final int length_,
			final int start_, final int stop_, final boolean selected_){
			require(strand_>=0 && strand_<=1 && start_>=0 && stop_>=start_ && stop_<length_,
				"Original CM request coordinates are zero-based inclusive in the source sense frame");
			source=source_;
			model=model_;
			seqid=seqid_;
			strand=strand_;
			length=length_;
			start=start_;
			stop=stop_;
			selected=selected_;
		}
		final Source source;
		final int model, strand, length, start, stop;
		final String seqid;
		final boolean selected;
		int finalStart, finalStop, cmLength;
		boolean measured, attempted, linked;
		String decision="NOT_INVOKED";
		Terminal terminal;
	}
	private static final class Witness{
		Witness(final Source source_, final String reference_, final String member_, final String modelName_, final String seqid_,
			final int start_, final int stop_, final int refStart_, final int refStop_, final int strand_, final long d5_, final long d3_){
			source=source_;
			reference=reference_;
			member=member_;
			modelName=modelName_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			refStart=refStart_;
			refStop=refStop_;
			strand=strand_;
			d5=d5_;
			d3=d3_;
		}
		boolean main(){return member.equals("RF00174");}
		long error(final int end){return end==0 ? d5 : d3;}
		final Source source;
		final String reference, member, modelName, seqid;
		final int start, stop, refStart, refStop, strand;
		final long d5, d3;
		int model=-1;
		Call call;
	}
	private static final class Arm{
		Arm(final String name_, final int models, final int inputs){
			name=name_;
			matchedByInput=new int[inputs];
			ends=new EndStat[models][2];
			for(int m=0; m<models; m++){
				for(int end=0; end<2; end++){ends[m][end]=new EndStat();}
			}
		}
		final String name;
		final EndStat[][] ends;
		final int[] matchedByInput;
		final ArrayList<Call> calls=new ArrayList<Call>();
		final ArrayList<Witness> witnesses=new ArrayList<Witness>();
		final ObjectMap<Target, Call> selected=new ObjectMap<Target, Call>(Target.class, Call.class);
		final ObjectMap<Target, Witness> references=new ObjectMap<Target, Witness>(Target.class, Witness.class);
		NcrnaCmInputPadding.Data padding;
		long measured, attempted, mainWitnesses, paired, harm5;
		double cpu=Double.NaN;
	}
	private static final class EndStat{
		void add(final Witness before, final Witness after, final int end, final int padding, final int base, final boolean fitAll){
			final long old=before.error(end), now=after.error(end);
			final int actual=clipped(after.call, end, padding), fixed=clipped(after.call, end, base);
			final int delta=actual-fixed;
			if(fitAll || (before.source.group&1)==0){
				selectionN++;
				selectionGroups.add(before.source.group);
				oldSquares=Math.addExact(oldSquares, old*old);
				newSquares=Math.addExact(newSquares, now*now);
			}
			n++;
			allOldSquares=Math.addExact(allOldSquares, old*old);
			allNewSquares=Math.addExact(allNewSquares, now*now);
			if(Math.abs(now)<Math.abs(old)){
				better++;
				betterDelta+=delta;
				if(delta<0){betterLess++;}else if(delta==0){betterEqual++;}else{betterMore++;}
			}
			if(Math.abs(now)>Math.abs(old)){
				worse++;
				worseDelta+=delta;
				if(delta<0){worseLess++;}else if(delta==0){worseEqual++;}else{worseMore++;}
			}
			if(nearHarm(old, now)){
				harms++;
				if(Math.abs(now)>3){crossings++;}
			}
			if(!before.modelName.equals(after.modelName)){switches++;}
			clipped+=actual;
			baselineClipped+=clipped(before.call, end, base);
		}
		final IntHashSet selectionGroups=new IntHashSet(4);
		long n, selectionN, better, worse, harms, crossings, switches, clipped, baselineClipped;
		long betterLess, betterEqual, betterMore, worseLess, worseEqual, worseMore, betterDelta, worseDelta;
		long oldSquares, newSquares, allOldSquares, allNewSquares;
	}
	private static final class Metric{
		void add(final Witness before, final long d5, final long d3){
			n++;
			if(d5==0){exact5++;}
			if(d3==0){exact3++;}
			if(d5==0 && d3==0){both++;}
			if(Math.abs(d5)<=3){within5++;}
			if(Math.abs(d3)<=3){within3++;}
			if(nearHarm(before.d5, d5)){
				harm5++;
				if(Math.abs(d5)>3){cross5++;}
			}
			if(nearHarm(before.d3, d3)){
				harm3++;
				if(Math.abs(d3)>3){cross3++;}
			}
			square5+=(double)d5*d5;
			square3+=(double)d3*d3;
		}
		long n, both, exact5, exact3, within5, within3, harm5, harm3, cross5, cross3;
		double square5, square3;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String[] grades, names, traces, profiles, padding, cpuValues, harmValues;
	private String libraryPath, bindings, audit, out;
	private String basePaddingPath;
	private int base5=20, base3=20, minSelectionCalls=1, minSelectionGroups=1;
	private boolean exportSelected, fitAll;
	private int expectedInputs;
	private NcrnaSupport.Library library;
	private Model[] models;
	private Arm[] arms;
	private final ArrayList<Source> panel=new ArrayList<Source>();
	private final ArrayList<Witness[]> common=new ArrayList<Witness[]>();
	private final ObjectMap<Target, Source> byInput=new ObjectMap<Target, Source>(Target.class, Source.class);
	private final ObjectMap<Target, Model> byModel=new ObjectMap<Target, Model>(Target.class, Model.class);
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder(), buffer=new ByteBuilder();
	private static final String[] END_NAMES={"5prime", "3prime"};
	private static final String[] POLICIES={"KEEP20_SPLICE", "FULL_CANDIDATE", "THREE_PRIME_ONLY", "SELECTION_HALF_MODEL_END"};
	private static final String[] COHORTS={"ALL", "SAME_MODEL", "MODEL_SWITCH",
		"CONFIGURED_BOTH20", "CONFIGURED_ANY_DIFFERENT", "CLIPPED_SAME_AS20", "CLIPPED_DIFFERENT_FROM20"};
	private static final String PANEL="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	private static final String AUDIT="assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled";
	private static final String WITNESSES="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
}
