package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;

import map.ObjectMap;
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
import structures.DoubleList;

/** Template.A_SampleByteFile argument/I/O shape, with native retained-key maps.
 * Reads CM proposal traces and saved normalized grades only. Selected queries
 * join by exact source, contig, strand, model and final coordinates; no rematching.
 * The primary endpoint oracle is a UNIQUE scored reference with reciprocal-half
 * overlap to the ORIGINAL accepted post-trim span, chosen without proposal ends.
 * Missing/ambiguous oracles remain explicit. Guard tests concern only finite
 * proposals the unchanged caller actually allowed, never existing length holds.
 * Predeclared scalar tests report discrimination and exact-end cost; this tool
 * does not install a guard or predict its downstream detection/path effects.
 * @author Brian Bushnell, Raiden
 */
public final class AnalyzeCobalaminCmProposals{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AnalyzeCobalaminCmProposals.class, false);
		try{new AnalyzeCobalaminCmProposals(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private AnalyzeCobalaminCmProposals(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("bindings")){bindings=b;}
			else if(a.equals("mode")){
				require("analysis".equals(b) || "scoreguard".equals(b) || "acceptance".equals(b), "Use mode=analysis, scoreguard or acceptance");
				scoreGuard="scoreguard".equals(b);
				acceptance="acceptance".equals(b);
			}
			else if(a.equals("audit")){audit=b;}
			else if(a.equals("rawgrades")){rawGrades=b.split(",");}
			else if(a.equals("nativecm")){nativeCm=Parse.parseBoolean(b);}
			else if(a.equals("eventassignments")){eventAssignments=Parse.parseBoolean(b);}
			else if(a.equals("members")){
				require(b!=null && !b.isEmpty(), "members= requires explicit distinct RF00174/RF01482/RF01689 IDs");
				memberMask=0;
				for(final String name:b.split(",", -1)){
					final int bit=1<<memberId(name);
					require((memberMask&bit)==0, "Each explicitly traced CM member must occur once");
					memberMask|=bit;
				}
			}
			else if(a.equals("profiles")){profiles=b.split(",");}
			else if(a.equals("acceptancescores")){acceptanceScores=b.split(",");}
			else if(a.equals("acceptancecosts")){acceptanceCosts=b.split(",");}
			else if(a.equals("names")){names=b.split(",");}
			else if(a.equals("expectedinputs")){expectedInputs=NcrnaNumbers.intValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown CM proposal analysis argument: "+arg);}
		}
		require(parser.in1!=null && parser.in2!=null && parser.out1!=null && bindings!=null && names!=null
			&& parser.out2==null && !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=trace[,trace] in2=grade_prefix[,prefix] names=label[,label] bindings=panel.sources.tsv out=fresh_prefix");
		require(profiles!=null, "profiles= requires each trace's actual CM invocation profile to detect omitted whole contigs");
		traces=parser.in1.split(",");
		grades=parser.in2.split(",");
		out=parser.out1;
		require(traces.length==grades.length && traces.length==names.length && traces.length==profiles.length && expectedInputs>0,
			"Every trace needs its own exact grade prefix and distinct label, on one declared panel");
		for(int i=0; i<names.length; i++){
			require(!names[i].isEmpty() && !names[i].equals("ALL"), "Labels must be nonempty and cannot use the pooled label ALL");
			for(int j=0; j<i; j++){require(!names[i].equals(names[j]), "Each input arm needs a distinct report label");}
			require(Tools.testInputFiles(false, true, traces[i], profiles[i], grades[i]+".queries.tsv", grades[i]+".references.tsv", grades[i]+".genomes.tsv"), "Trace, CM counters and normalized grades must exist");
		}
		require(Tools.testInputFiles(false, true, bindings), "The exact panel bindings must be readable");
		require(!acceptance || audit!=null && Tools.testInputFiles(false, true, audit),
			"Acceptance curves require the resolved-source-group audit before any selection/judging comparison");
		require(!acceptance || rawGrades!=null && rawGrades.length==traces.length && Tools.testInputFiles(false, true, rawGrades),
			"Acceptance reporting requires each exact grade.overlap.raw.tsv to retain recorded TP_INEXACT labels merged by the normalized exporter");
		require((acceptanceScores==null)==(acceptanceCosts==null), "Candidate scores and independent acceptance costs must be supplied together");
		require(acceptanceScores==null || acceptance && acceptanceScores.length==traces.length && acceptanceCosts.length==traces.length
			&& Tools.testInputFiles(false, true, acceptanceScores) && Tools.testInputFiles(false, true, acceptanceCosts),
			"Each acceptance trace needs its own exact score-event and cost tables; other analysis modes retain their original schema");
		final String[] outputs=acceptance ? new String[]{out+".acceptance_calls.tsv", out+".acceptance_curves.tsv", out+".acceptance_coverage.tsv", out+".acceptance_native.tsv"}
			: scoreGuard ? new String[]{out+".score_guard_calls.tsv", out+".score_guard_summary.tsv"}
			: new String[]{out+".calls.tsv", out+".bins.tsv", out+".tests.tsv", out+".summary.tsv"};
		require(Tools.testOutputFiles(false, false, false, outputs), "All requested proposal-analysis outputs must be fresh");
		require(!eventAssignments || acceptance && acceptanceScores!=null && traces.length==1,
			"Event assignments require one exact acceptance capture; event indices are local to its saved score file");
		if(eventAssignments){RrnaResourceIO.fresh(out, ".event_assignments.tsv");}
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reading           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		if(acceptance){
			processAcceptance();
			return;
		}
		if(scoreGuard){
			processScoreGuard();
			return;
		}
		final Stats pooled=new Stats();
		try(Output calls=new Output(out+".calls.tsv"); Output bins=new Output(out+".bins.tsv");
			Output tests=new Output(out+".tests.tsv"); Output summary=new Output(out+".summary.tsv")){
			calls.writer.println(CALL_HEADER);
			bins.writer.println("arm\tpopulation\tvariable\tbin\tn\tworse_any\tworse_5prime\tworse_3prime\toriginal_exact_ends\tproposed_exact_ends\toriginal_within3_ends\tproposed_within3_ends\tselected");
			tests.writer.println("arm\ttest\tn\treject_worse\treject_other\tkeep_worse\tkeep_other\tproposed_exact_ends_kept\tproposed_exact_ends_rejected\toriginal_exact_ends_restored\tyouden\tbest_youden\tqualification");
			summary.writer.println("arm\tmeasure\tcount\tqualification");
			for(int side=0; side<traces.length; side++){
				reset();
				readBindings();
				readTrace(traces[side]);
				checkProfile(profiles[side]);
				readReferences(grades[side]);
				readQueries(grades[side]);
				final Stats stats=new Stats();
				for(final Call call:order){analyze(call, names[side], stats, calls);}
				checkSelectedLinks();
				stats.write(names[side], bins, tests, summary);
				pooled.add(stats);
			}
			pooled.write("ALL", bins, tests, summary);
		}
		System.err.println("COBALAMIN_CM_PROPOSAL_ANALYSIS_COMPLETE");
	}
	/** Reports the actual GA holds, without changing a proposal or replaying a path. */
	private void processScoreGuard(){
		try(Output calls=new Output(out+".score_guard_calls.tsv"); Output summary=new Output(out+".score_guard_summary.tsv")){
			calls.writer.println(CALL_HEADER+"\tmasked_original_oracles");
			summary.writer.println(GuardCounts.HEADER);
			for(int side=0; side<traces.length; side++){
				reset();
				readBindings();
				readTrace(traces[side]);
				checkProfile(profiles[side]);
				readReferences(grades[side]);
				readQueries(grades[side]);
				final GuardCounts all=new GuardCounts(), selected=new GuardCounts(), unselected=new GuardCounts(), assigned=new GuardCounts();
				for(final Call call:order){
					analyze(call, names[side], null, null);
					if(!call.decision.equals("SCORE_REJECTED")){continue;}
					require(call.measured && call.proposedStart>=0 && call.scored && Float.isFinite(call.score) && call.score<call.ga
						&& call.effectiveStart==call.start && call.effectiveStop==call.stop,
						"SCORE_REJECTED must retain the original interval and the finite below-GA proposal captured by CobalaminCmEndpointRefiner.refine");
					all.add(call, call.oracle, call.oracleStatus);
					(call.selected ? selected : unselected).add(call, call.oracle, call.oracleStatus);
					if(call.query!=null && call.query.assigned!=null){assigned.add(call, call.query.assigned, "RECORDED_TP_REFERENCE");}
					writeCall(call, names[side], call.oracle, call.oracleStatus, false, calls);
				}
				checkSelectedLinks();
				require(all.held==selected.held+unselected.held && assigned.held<=selected.held,
					"Selected and unselected actual holds partition the complete request cohort; recorded TPs are a selected subset");
				all.write(names[side], "ALL_HELD_ORIGINAL_ORACLE", summary);
				selected.write(names[side], "SELECTED_HELD_ORIGINAL_ORACLE", summary);
				unselected.write(names[side], "UNSELECTED_HELD_ORIGINAL_ORACLE", summary);
				assigned.write(names[side], "SELECTED_HELD_RECORDED_TP_REFERENCE", summary);
			}
		}
		System.err.println("COBALAMIN_CM_SCORE_GUARD_ANALYSIS_COMPLETE");
	}
	private void checkSelectedLinks(){
		for(final Genome genome:genomes){
			for(final Query query:genome.queries){
				require((memberMask&(1<<query.member))==0 || query.linked,
					"Every emitted query in the explicitly traced member scope needs one exact selected invocation link");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Acceptance Curves      ----------------*/
	/*--------------------------------------------------------------*/

	/** Describes fixed emitted calls only; it never predicts the caller's new path
	 * or recalls candidates that the measured identity gate excluded. Curves have
	 * one row per selection-observed cutoff and retain selection/judging counts separately.
	 */
	private void processAcceptance(){
		try(Output calls=new Output(out+".acceptance_calls.tsv"); Output curves=new Output(out+".acceptance_curves.tsv");
			Output coverage=new Output(out+".acceptance_coverage.tsv"); Output offsets=new Output(out+".acceptance_native.tsv");
			Output events=eventAssignments ? new Output(out+".event_assignments.tsv") : null){
			if(events!=null){events.writer.println(EVENT_ASSIGNMENT_HEADER);}
			calls.writer.println("arm\tinput_id\tresolved_group\tsplit\tcontig\tstrand\tstart\tstop\tmember\tmodel\toutcome\tmatched_inexact\treference_id\taccepted_identity\tendpoint_identity\tcm_bits\tga_bits\tcm_minus_ga\tcm_over_ga\tcm_decision\tnative_score_status\tnative_cm_bits\tnative_ga_bits\tanchored_minus_native\thbm_score\thbm_status\tboundary_different");
			curves.writer.println("arm\tscope\tmodel\tvariable\tcutoff_inclusive\tsplit\tmatched_kept\tmatched_removed\tunmatched_kept\tunmatched_removed\tinexact_kept\tinexact_removed\tmasked_kept\tmasked_removed\tmissing_matched\tmissing_unmatched\tmissing_inexact\tmissing_masked\tqualification");
			coverage.writer.println("arm\tmember\tsplit\tmatched\tunmatched\tinexact\tmasked\ttraced\tidentity_available\tcm_available\tnative_pair_available\thbm_available\tqualification");
			offsets.writer.println("arm\tsplit\tsubset\tpairs\tga_mismatches\tmean_bits\tsample_sd_bits\tmin_bits\tp05_bits\tmedian_bits\tp95_bits\tmax_bits\tqualification");
			for(int side=0; side<traces.length; side++){
				reset();
				readBindings();
				readAcceptanceGroups();
				readTrace(traces[side]);
				if(acceptanceScores!=null){readAcceptanceScores(acceptanceScores[side], acceptanceCosts[side]);}
				checkProfile(profiles[side]);
				readReferences(grades[side]);
				readQueries(grades[side]);
				readRecordedInexact(rawGrades[side]);
				if(nativeCm){readNativeScores();}
				final ArrayList<Call> selected=new ArrayList<Call>();
				final long[][][] counts=new long[3][2][9];
				final NativeOffsets[][] nativeOffsets=new NativeOffsets[2][2];
				for(int split=0; split<2; split++){
					for(int subset=0; subset<2; subset++){nativeOffsets[split][subset]=new NativeOffsets();}
				}
				for(final Genome genome:genomes){
					for(final Query query:genome.queries){
						final long[] c=counts[query.member][genome.group&1];
						c[category(query)]++;
						if(inexact(query)){c[2]++;}
					}
				}
				for(final Call call:order){
					final boolean minorBeforeAcceptance=acceptanceScores!=null && member(call.member)!=0 && !call.returned && !call.selected;
					require(call.measured || call.captured || minorBeforeAcceptance || !call.invoked && !call.known,
						"A known CM invocation or fixed-end skip requires its captured M record, as in endpoint analysis");
					require(call.proposedStart<0 || call.scored && Float.isFinite(call.score),
						"Finite proposal coordinates require a finite actual search score, even when endpoints are ignored");
					if(!call.selected){continue;}
					linkSelected(call);
					require(!call.scored || call.ga>0, "Score/GA acceptance curves require the measured model's positive gathering threshold");
					bindNativeScore(call);
					selected.add(call);
					final long[] c=counts[member(call.member)][call.genome.group&1];
					c[4]++;
					if(Float.isFinite(call.acceptedIdentity)){c[5]++;}
					if(call.scored){c[6]++;}
					if(Double.isFinite(call.nativeBits) && call.scored && Float.isFinite(call.score)){
						c[7]++;
						nativeOffsets[call.genome.group&1][0].add(call);
						if(inexact(call.query)){nativeOffsets[call.genome.group&1][1].add(call);}
					}
					writeAcceptanceCall(names[side], call, calls);
					if(events!=null){writeEventAssignment(call, events);}
				}
				checkSelectedLinks();
				for(int split=0; split<2; split++){
					for(int subset=0; subset<2; subset++){
						nativeOffsets[split][subset].write(names[side], split, subset, offsets);
					}
				}
				for(int m=0; m<3; m++){
					for(int split=0; split<2; split++){
						final long[] c=counts[m][split];
						require(c[2]<=c[0] && c[4]<=c[0]+c[1]+c[3] && c[5]<=c[4] && c[6]<=c[4] && c[7]<=c[6],
							"Inexact calls are a matched subset; score availability is confined to traced emitted calls");
						require((memberMask&(1<<m))==0 || c[4]==c[0]+c[1]+c[3],
							"Every emitted query in the declared trace scope must appear once in acceptance coverage");
						row.clear().append(names[side]).tab().append(MEMBERS[m]).tab().append(splitName(split));
						for(int i=0; i<9; i++){row.tab().append(c[i]);}
						coverage.writer.print(row.append("\tRECORDED_EMITTED_CALLS_HBM_NOT_STORED_IN_CM_TRACE\n"));
					}
				}
				final Call[] sorted=selected.toArray(new Call[selected.size()]);
				for(int axis=0; axis<3; axis++){
					writeAcceptanceCurves(names[side], sorted, axis, false, curves);
					writeAcceptanceCurves(names[side], sorted, axis, true, curves);
				}
			}
		}
		System.err.println("COBALAMIN_CM_ACCEPTANCE_CURVES_COMPLETE");
	}
	/** Native input/worker/candidate identity is the join key. Coordinates alone
	 * collide across repeated windows and cannot establish the scoring source. */
	private void readAcceptanceScores(final String scores, final String costs){
		final long[][][] counts=new long[expectedInputs][3][5];
		long eventIndex=0;
		try(Input input=new Input(scores)){
			RrnaResourceIO.header(input.reader, CobalaminCmAcceptance.HEADER, scores);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==17, "Acceptance candidate schema must retain the native trace identity and score-window bounds");
				final int ordinal=p.parseInt(0), worker=p.parseInt(15), member=memberColumn(3);
				final long id=p.parseLong(16);
				key.clear().append(ordinal).tab().append(worker).tab().append(id);
				probe.set(key.array, 0, key.length());
				final Call call=byTrace.get(probe);
				require(call!=null && !call.captured && worker>=0 && id>=0 && member(call.member)==member
					&& p.termEquals(call.seqid, 1) && strand(2)==call.strand && p.termEquals(call.model, 4)
					&& p.parseInt(5)==call.length && p.parseInt(6)+1==genomicStart(call, call.effectiveStart, call.effectiveStop)
					&& p.parseInt(7)+1==genomicStop(call, call.effectiveStart, call.effectiveStop),
					"Each acceptance event must bind exactly one actual native candidate, model and effective genomic interval");
				final float identity=Float.intBitsToFloat(p.parseInt(8)), ga=Float.intBitsToFloat(p.parseInt(10));
				final boolean scored=!p.termEquals("", 9), reused=p.termEquals("REUSED_ENDPOINT", 11);
				final float score=scored ? Float.intBitsToFloat(p.parseInt(9)) : Float.NaN;
				require(Float.isFinite(identity) && identity>=0 && identity<=1 && Float.isFinite(ga) && ga>0
					&& (!scored || Float.isFinite(score) || score==Float.NEGATIVE_INFINITY)
					&& (!Float.isFinite(call.acceptedIdentity) || Float.floatToIntBits(call.acceptedIdentity)==p.parseInt(8)),
					"Acceptance identity must equal the verifier's recorded bits; missing and impossible scores remain distinct");
				require(reused || p.termEquals("SCORE_ONLY", 11), "Calibration captures must actually request or reuse scores, never identity-only filtering");
				require(member!=0 || reused==call.measured, "The main scorer reuses every recorded endpoint result, including unavailable searches");
				require(member!=0 || call.measured || !call.invoked && !call.known,
					"A main endpoint invocation cannot lose its M record and masquerade as a new score-only request");
				require(!reused || call.measured && call.scored==scored && Float.floatToIntBits(call.ga)==p.parseInt(10)
					&& (!scored || Float.floatToIntBits(call.score)==p.parseInt(9)),
					"A reused endpoint score must match that exact candidate's existing M record bit for bit");
				require(p.termEquals("true", 12), "Calibration captures retain every admitted candidate; filter results require a separate caller grade");
				final boolean window=!p.termEquals("", 13);
				require(window==!p.termEquals("", 14) && (!scored || window), "Measured scores require their actual target window; unavailable searches may have none");
				if(window){require(p.parseInt(13)>=0 && p.parseInt(14)>=p.parseInt(13) && p.parseInt(14)<call.length,
					"A captured score window must remain within its original genomic contig");}
				call.captured=true;
				call.captureEventIndex=eventIndex++;
				call.acceptedIdentity=identity;
				call.scored=scored;
				call.score=score;
				call.ga=ga;
				if(!call.measured){call.decision="SCORE_ONLY";}
				final long[] c=counts[ordinal][member];
				c[0]++;
				if(scored){c[1]++;}
				c[reused ? 2 : 3]++;
				if(!reused && member==0){acceptanceExtraMain++;}
			}
		}
		for(final Call call:order){require(!call.returned && !call.selected || call.captured,
			"Every returned or emitted candidate needs its exact acceptance event; earlier native minor boundary failures do not reach that stage");}
		try(Input input=new Input(costs)){
			RrnaResourceIO.header(input.reader, CobalaminCmAcceptance.COST_HEADER, costs);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11, "Acceptance costs distinguish requests, actual DP fills and CPU subdivisions");
				final int ordinal=p.parseInt(0), member=memberColumn(1);
				ordinal(ordinal);
				final long[] c=counts[ordinal][member];
				require(c[4]++==0 && p.parseLong(2)==c[0] && p.parseLong(3)==c[0] && p.parseLong(4)==c[1]
					&& p.parseLong(5)==c[2] && p.parseLong(6)==c[3],
					"Independent per-input/member totals must conserve every captured, kept, measured, reused and extra-request event");
				for(int field=7; field<11; field++){require(p.parseLong(field)>=0, "Actual DP and worker CPU counters cannot be negative");}
				require(p.parseLong(8)>=p.parseLong(9) && p.parseLong(8)-p.parseLong(9)>=p.parseLong(10),
					"CM and anchor CPU subdivide extra acceptance CPU; they must not be added twice to caller time");
			}
		}
		for(final long[][] genome:counts){
			for(final long[] member:genome){require(member[4]==1, "Every input, including negatives, needs all three independent member cost rows");}
		}
	}
	private int memberColumn(final int column){
		for(int member=0; member<MEMBERS.length; member++){if(p.termEquals(MEMBERS[member], column)){return member;}}
		throw new IllegalArgumentException("Acceptance rows require an explicit cobalamin member");
	}
	private void readAcceptanceGroups(){
		try(Input input=new Input(audit)){
			RrnaResourceIO.header(input.reader, "assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled", audit);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==8 && p.parseInt(4)>=0, "Resolved source groups must be present before splitting acceptance evidence");
				probe.set(p, 1);
				final Genome genome=byInput.get(probe);
				if(genome==null){continue;}
				require(genome.group<0 && p.termEquals(genome.assembly, 0),
					"Every panel input binds its exact audited assembly and group once; the audit may name the original FASTA behind a corpus symlink");
				genome.group=p.parseInt(4);
			}
		}
		for(final Genome genome:genomes){require(genome.group>=0, "A panel input lacks its resolved source group: "+genome.id);}
	}
	/** Optional canonical GFF evidence is joined by the grader's fixed reference ID.
	 * An overlay-only reference stays unavailable; geometry is never used to rematch.
	 */
	private void readNativeScores(){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(final Genome genome:genomes){
			boolean needed=false;
			for(final Query query:genome.queries){
				if(query.assigned!=null && (memberMask&(1<<query.member))!=0){needed=true;break;}
			}
			if(!needed){continue;}
			for(final NcrnaReferenceGff.Row nativeRow:codec.readFamily(genome.referencePath, genome.id, "cobalamin")){
				for(final Reference ref:genome.refs){
					if(!ref.id.equals(nativeRow.id)){continue;}
					require(!ref.nativeSeen && ref.seqid.equals(nativeRow.locus.seqid) && ref.start==nativeRow.locus.start
						&& ref.stop==nativeRow.locus.stop && (ref.strand==0 ? '+' : '-')==nativeRow.locus.strand,
						"Reference IDs encode input/family/geometry; native measurement lookup must not change that geometry or rematch the grade");
					ref.nativeSeen=true;
					ref.nativeScores=nativeRow.scores;
				}
			}
		}
	}
	/** The inclusive grader freezes strict matches, then admits residual overlaps.
	 * Its recorded TP_INEXACT class is not the same as differing endpoints.
	 */
	private void readRecordedInexact(final String path){
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, GradeNcrnaReference.RAW_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==18, "Recorded strict/inexact assignment needs the canonical raw grade schema");
				if(!p.termEquals("QUERY", 1) || !p.termEquals("cobalamin", 2)){continue;}
				final Genome genome=input();
				final int start=p.parseInt(4), stop=p.parseInt(5), strand=strand(6);
				Query found=null;
				for(final Query query:genome.queries){
					if(query.start==start && query.stop==stop && query.strand==strand
						&& p.termEquals(query.seqid, 3) && rawModelEquals(query.model)){
						require(found==null && !query.rawLinked, "A raw grader label must bind one distinct normalized query exactly");
						found=query;
					}
				}
				require(found!=null, "Raw grade query has no exact normalized input/contig/strand/model/coordinate counterpart");
				final boolean residual=p.termEquals("TP_INEXACT", 7), matched=residual || p.termEquals("TP", 7);
				require(matched==(found.assigned!=null) && (!matched || p.termEquals(found.reference, 10))
					&& (found.outcome.equals("MASKED") ? p.termEquals("MASKED", 7) : matched || p.termEquals("FP", 7)),
					"Recorded strict/inexact labels and reference IDs must preserve the exact normalized grade, without rematching");
				found.rawLinked=true;
				found.residualInexact=residual;
			}
		}
		for(final Genome genome:genomes){
			for(final Query query:genome.queries){require(query.rawLinked, "Every normalized cobalamin query needs its recorded raw strict/inexact label");}
		}
	}
	/** Compares the first whitespace-delimited model token without field Strings. */
	private boolean rawModelEquals(final String model){
		p.setBounds(15);
		final int start=p.a(), limit=p.b();
		final byte[] line=p.line();
		int end=start;
		while(end<limit && line[end]>32){end++;}
		if(end-start!=model.length()){return false;}
		for(int i=start; i<end; i++){if(line[i]!=model.charAt(i-start)){return false;}}
		return true;
	}
	private static boolean inexact(final Query query){return query.residualInexact;}
	private static boolean boundaryDifferent(final Query query){
		return query.assigned!=null && (query.start!=query.assigned.start || query.stop!=query.assigned.stop);
	}
	private static int category(final Query query){
		return query.assigned!=null ? 0 : query.outcome.equals("MASKED") ? 3 : 1;
	}
	private void bindNativeScore(final Call call){
		assert(call.query!=null) : "Native scores attach only after the exact emitted-query join";
		final Reference ref=call.query.assigned;
		if(!nativeCm){call.nativeStatus="NOT_REQUESTED";}
		else if(ref==null){call.nativeStatus="NO_MATCHED_REFERENCE";}
		else if(!ref.nativeSeen){call.nativeStatus="REFERENCE_NOT_IN_BOUND_GFF";}
		else if(ref.nativeScores==null){call.nativeStatus="NO_SCORE_ATTRIBUTES";}
		else if(ref.nativeScores.missing){call.nativeStatus="NO_CURRENT_NATIVE_HIT";}
		else{
			final NcrnaReferenceGff.Scores scores=ref.nativeScores;
			final String cm=MEMBERS[member(call.member)];
			int from=0, found=-1;
			for(int index=0; index<scores.size(); index++){
				final int comma=scores.models.indexOf(',', from), to=comma<0 ? scores.models.length() : comma;
				if(to-from==cm.length() && scores.models.regionMatches(from, cm, 0, cm.length())){
					require(found<0, "A native score list cannot ambiguously repeat the requested CM identity");
					found=index;
				}
				from=to+1;
			}
			if(found<0){call.nativeStatus="NO_SAME_CM_MEASUREMENT";}
			else{
				call.nativeBits=scores.bitsTenths(found)/10.0;
				call.nativeGa=Double.parseDouble(scores.gaBits(found));
				call.nativeStatus="SAME_CM_PRINTED_TENTH_BIT_REFERENCE_SPAN";
			}
		}
	}
	private void writeAcceptanceCall(final String arm, final Call call, final Output output){
		final Query q=call.query;
		row.clear().append(arm).tab().append(call.genome.id).tab().append(call.genome.group).tab().append(splitName(call.genome.group&1))
			.tab().append(call.seqid).tab().append(call.strand==0 ? '+' : '-').tab().append(call.finalStart).tab().append(call.finalStop)
			.tab().append(MEMBERS[q.member]).tab().append(call.model).tab().append(q.outcome).tab().append(inexact(q)).tab().append(q.reference).tab();
		appendAvailable(row, call.acceptedIdentity).tab();
		appendAvailable(row, call.endpointIdentity).tab();
		appendAvailable(row, call.scored ? call.score : Double.NaN).tab();
		appendAvailable(row, call.measured || call.captured ? call.ga : Double.NaN).tab();
		appendAvailable(row, acceptanceValue(call, 1)).tab();
		appendAvailable(row, acceptanceValue(call, 2)).tab().append(call.decision).tab().append(call.nativeStatus).tab();
		appendAvailable(row, call.nativeBits).tab();
		appendAvailable(row, call.nativeGa).tab();
		appendAvailable(row, call.scored && Float.isFinite(call.score) ? call.score-call.nativeBits : Double.NaN)
			.append("\tNA\tNOT_STORED_IN_CM_TRACE\t").append(boundaryDifferent(q)).nl();
		output.writer.print(row);
	}
	/** Identity-preserving companion for the OR fitter. This exports the already
	 * checked native event/selected-query join, never a new coordinate matcher. */
	private void writeEventAssignment(final Call call, final Output output){
		require(call.captured && call.captureEventIndex>=0 && call.query!=null && call.query.linked,
			"An event assignment requires one captured native event and its exact selected query");
		row.clear().append(call.captureEventIndex).tab().append(call.genome.id).tab().append(call.genome.group)
			.tab().append(call.ordinal).tab().append(call.worker).tab().append(call.id).tab().append(call.seqid)
			.tab().append(call.strand==0 ? '+' : '-').tab().append(call.finalStart).tab().append(call.finalStop)
			.tab().append(MEMBERS[call.query.member]).tab().append(call.model).tab().append(call.query.outcome)
			.tab().append(call.query.reference).tab().append(Float.floatToIntBits(call.acceptedIdentity)).tab();
		if(call.scored){row.append(Float.floatToIntBits(call.score));}
		output.writer.print(row.tab().append(Float.floatToIntBits(call.ga)).nl());
	}
	/** A sorted sweep computes every inclusive cutoff in O(N log N), with ties
	 * removed together. No threshold-by-call rescans, boxing or transient keys.
	 */
	private void writeAcceptanceCurves(final String arm, final Call[] calls, final int axis, final boolean perModel, final Output output){
		Arrays.sort(calls, new Comparator<Call>(){
			@Override
			public int compare(final Call a, final Call b){
				final int model=perModel ? a.model.compareTo(b.model) : 0;
				return model!=0 ? model : Double.compare(acceptanceValue(a, axis), acceptanceValue(b, axis));
			}
		});
		for(int first=0; first<calls.length; ){
			int end=first+1;
			while(end<calls.length && (!perModel || calls[end].model.equals(calls[first].model))){end++;}
			final long[][] totals=new long[2][4], kept=new long[2][4], missing=new long[2][4];
			for(int i=first; i<end; i++){
				final Call c=calls[i];
				final int split=c.genome.group&1, kind=category(c.query);
				totals[split][kind]++;
				final long[] target=Double.isNaN(acceptanceValue(c, axis)) ? missing[split] : kept[split];
				target[kind]++;
				if(inexact(c.query)){
					totals[split][2]++;
					target[2]++;
				}
			}
			for(int i=first; i<end && !Double.isNaN(acceptanceValue(calls[i], axis)); ){
				final double cutoff=acceptanceValue(calls[i], axis);
				int tiedEnd=i;
				boolean fromSelection=false;
				while(tiedEnd<end && Double.compare(acceptanceValue(calls[tiedEnd], axis), cutoff)==0){
					fromSelection|=(calls[tiedEnd].genome.group&1)==0;
					tiedEnd++;
				}
				for(int split=0; fromSelection && split<2; split++){
					row.clear().append(arm).tab().append(perModel ? "MODEL" : "FAMILY").tab().append(perModel ? calls[first].model : "ALL")
						.tab().append(ACCEPTANCE_AXES[axis]).tab();
					appendAvailable(row, cutoff).tab().append(splitName(split));
					for(int kind=0; kind<4; kind++){
						final long removed=totals[split][kind]-kept[split][kind]-missing[split][kind];
						require(removed>=0, "A cutoff sweep partitions observed values into kept and removed, with missing separate");
						row.tab().append(kept[split][kind]).tab().append(removed);
					}
					for(final long n:missing[split]){row.tab().append(n);}
					output.writer.print(row.append("\tFIXED_EMITTED_CALL_FILTER_NO_PATH_OR_RECALL_PREDICTION\n"));
				}
				while(i<end && Double.compare(acceptanceValue(calls[i], axis), cutoff)==0){
					final Call c=calls[i++];
					final long[] counts=kept[c.genome.group&1];
					counts[category(c.query)]--;
					if(inexact(c.query)){counts[2]--;}
				}
			}
			first=end;
		}
	}
	private static double acceptanceValue(final Call call, final int axis){
		assert(axis>=0 && axis<3) : "Acceptance axes are verifier identity, CM bit margin and CM/GA ratio";
		return axis==0 ? call.acceptedIdentity : !call.scored ? Double.NaN : axis==1 ? call.score-(double)call.ga : call.score/(double)call.ga;
	}
	private static ByteBuilder appendAvailable(final ByteBuilder buffer, final double value){
		return Double.isNaN(value) ? buffer.append("NA") : buffer.appendSlow(value);
	}
	private static String splitName(final int split){return split==0 ? "SELECTION" : "JUDGING";}

	/** Retains primitive differences only; Welford's update avoids subtracting two
	 * large squared sums when the offset distribution has a narrow spread.
	 */
	private static final class NativeOffsets{
		void add(final Call call){
			assert(call.query.assigned!=null && call.scored && Float.isFinite(call.score) && Double.isFinite(call.nativeBits))
				: "Native comparisons require the same recorded TP reference and two measured finite scores";
			final double value=call.score-call.nativeBits;
			values.add(value);
			final double delta=value-mean;
			mean+=delta/values.size;
			m2+=delta*(value-mean);
			if(call.ga!=(float)call.nativeGa){gaMismatch++;}
		}
		void write(final String arm, final int split, final int subset, final Output output){
			assert(split>=0 && split<2 && subset>=0 && subset<2) : "Offset tables distinguish group parity and the overlapping matched-inexact subset";
			Arrays.sort(values.array, 0, values.size);
			final ByteBuilder b=new ByteBuilder();
			b.append(arm).tab().append(splitName(split)).tab().append(subset==0 ? "MATCHED" : "MATCHED_INEXACT")
				.tab().append(values.size).tab().append(gaMismatch).tab();
			appendAvailable(b, values.size==0 ? Double.NaN : mean).tab();
			appendAvailable(b, values.size<2 ? Double.NaN : Math.sqrt(Math.max(0, m2/(values.size-1))));
			for(final double q:OFFSET_QUANTILES){
				b.tab();
				appendAvailable(b, values.size==0 ? Double.NaN : values.array[Math.max(0, (int)Math.ceil(q*values.size)-1)]);
			}
			output.writer.print(b.append("\tANCHORED_CALLER_MINUS_PRINTED_CMSEARCH_REFERENCE_SPAN_NEAREST_RANK_QUANTILES\n"));
		}
		final DoubleList values=new DoubleList(16);
		double mean, m2;
		int gaMismatch;
	}
	private void reset(){
		genomes.clear();
		byInput.clear();
		byTrace.clear();
		byContig.clear();
		contigs.clear();
		order.clear();
		acceptanceExtraMain=0;
	}
	private void readBindings(){
		try(Input in=new Input(bindings)){
			RrnaResourceIO.header(in.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", bindings);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5 && p.parseInt(0)==genomes.size(), "Panel bindings must preserve every input ordinal exactly once");
				final Genome genome=new Genome(p.parseString(2), p.parseString(3), acceptance ? p.parseString(1) : null,
					acceptance && nativeCm ? p.parseString(4) : null);
				probe.set(p, 2);
				require(!byInput.contains(probe), "Panel input IDs must be unique");
				byInput.put(probe.copy(), genome);
				genomes.add(genome);
			}
		}
		require(genomes.size()==expectedInputs, "CM diagnostic panel must contain every declared input, including negative genomes");
	}
	private void readTrace(final String path){
		int inputs=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "#ncrna_cm_proposals\t4", path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("I", 0)){
					final Genome genome=ordinal(p.parseInt(1));
					require(p.terms()==3 && !genome.seen && p.termEquals(genome.fasta, 2), "Trace input must bind its exact declared FASTA once");
					genome.seen=true;
					inputs++;
				}else if(p.termEquals("C", 0)){
					require(p.terms()==37, "CM trace candidate schema changed");
					traceKey();
					require(!byTrace.contains(probe), "Candidate invocation IDs must be unique within input/worker");
					final Call call=new Call(ordinal(p.parseInt(1)), p.parseInt(1), p.parseInt(2), p.parseLong(3), p.parseString(5),
						p.parseInt(6), p.parseString(7), p.parseString(8), p.parseInt(9), p.parseInt(10), p.parseInt(11));
					require((memberMask&(1<<member(call.member)))!=0, "A trace cannot silently contain members outside its declared scope");
					call.effectiveStart=p.parseInt(23);
					call.effectiveStop=p.parseInt(24);
					call.returned=bit(31);
					call.selected=bit(32);
					call.known=bit(14);
					call.invoked=bit(16);
					if(acceptance && call.known){
						call.acceptedIdentity=Float.intBitsToFloat(p.parseInt(18));
						call.endpointIdentity=Float.intBitsToFloat(p.parseInt(19));
						require(Float.isFinite(call.acceptedIdentity) && call.acceptedIdentity>=0 && call.acceptedIdentity<=1
							&& Float.isFinite(call.endpointIdentity) && call.endpointIdentity>=0 && call.endpointIdentity<=1,
							"Known acceptance features retain the actual verifier and endpoint-feature identity fractions");
					}
					require(call.genome.seen, "An input binding precedes all of its diagnostic candidate records");
					final Contig contig=contig(call.ordinal, call.worker, p.parseLong(4), call.seqid, call.length);
					require(!contig.closed, "A candidate cannot follow its own completed contig terminal");
					contig.calls++;
					if(call.returned){contig.returned++;}
					if(call.selected){contig.selected++;}
					if(!call.known){contig.unknown++;}
					if(call.selected){
						call.finalStart=p.parseInt(25)+1;
						call.finalStop=p.parseInt(26)+1;
						require(call.finalStart==genomicStart(call, call.effectiveStart, call.effectiveStop)
							&& call.finalStop==genomicStop(call, call.effectiveStart, call.effectiveStop), "Selected trace must preserve its effective endpoints through strand conversion");
					}
					traceKey();
					byTrace.put(probe.copy(), call);
					order.add(call);
				}else if(p.termEquals("M", 0)){
					require(p.terms()==21, "CM proposal schema changed");
					traceKey();
					final Call call=byTrace.get(probe);
					require(call!=null && !call.measured && p.parseInt(4)==call.start && p.parseInt(5)==call.stop
						&& p.parseInt(6)==call.stop-call.start+1, "A proposal binds one existing candidate and its original accepted span");
					call.measured=true;
					if(!p.termEquals("", 7)){
						call.proposedStart=p.parseInt(7);
						call.proposedStop=p.parseInt(8);
						require(call.proposedStart>=0 && call.proposedStop>=call.proposedStart && call.proposedStop<call.length
							&& p.parseInt(9)==call.proposedStop-call.proposedStart+1, "Raw proposal coordinates and inclusive length must agree");
					}else{require(p.termEquals("", 8) && p.termEquals("", 9), "Unavailable proposals have no invented endpoints or length");}
					call.consensusLength=p.parseInt(10);
					call.cmLength=p.parseInt(11);
					call.scored=!p.termEquals("", 12);
					if(call.scored){call.score=Float.intBitsToFloat(p.parseInt(12));}
					call.ga=Float.intBitsToFloat(p.parseInt(13));
					call.width=p.parseInt(14);
					call.contacts16=p.parseInt(15);
					call.contacts32=p.parseInt(16);
					call.anchor=p.parseInt(17);
					call.anchorOrigin=p.termEquals("", 18) ? -1 : p.parseInt(18);
					call.unavailable=bit(19);
					call.decision=p.parseString(20);
					require(call.consensusLength>0 && call.cmLength>0 && Float.isFinite(call.ga)
						&& (!call.scored || Float.isFinite(call.score) || call.score==Float.NEGATIVE_INFINITY), "Only finite GA and finite/impossible attempted scores are valid");
				}else{
					require(p.termEquals("T", 0) && p.terms()==11, "Unknown CM trace record or changed contig terminal schema");
					final Contig contig=contig(p.parseInt(1), p.parseInt(2), p.parseLong(3), p.parseString(4), p.parseInt(5));
					require(!contig.closed && contig.calls==p.parseInt(6) && contig.returned==p.parseInt(7)
						&& contig.selected==p.parseInt(8) && contig.unknown==p.parseInt(9) && p.parseInt(10)==0,
						"A contig terminal must conserve every actual C record and its returned/selected/unknown fate exactly once");
					contig.closed=true;
				}
			}
		}
		require(inputs==expectedInputs, "The diagnostic must account for the entire panel input roster");
		for(final Contig contig:contigs){require(contig.closed, "Every candidate must belong to a closed source/worker/read/contig terminal");}
	}
	private Contig contig(final int ordinal, final int worker, final long read, final String name, final int length){
		require(worker>=0 && read>=0 && length>=0 && ordinal(ordinal).seen, "Contig terminals require their real input, worker and read identities");
		key.clear().append(ordinal).tab().append(worker).tab().append(read);
		probe.set(key.array, 0, key.length());
		Contig contig=byContig.get(probe);
		if(contig==null){
			contig=new Contig(name, length);
			byContig.put(probe.copy(), contig);
			contigs.add(contig);
		}
		require(contig.name.equals(name) && contig.length==length, "One source/worker/read identity cannot change its contig or length");
		return contig;
	}
	private void checkProfile(final String path){
		final long[] observed=new long[3];
		final long[] held=new long[3];
		final boolean[] seen=new boolean[3];
		for(final Call call:order){
			if(call.measured){observed[member(call.member)]++;}
			if(call.decision.equals("SCORE_REJECTED")){held[member(call.member)]++;}
		}
		observed[0]+=acceptanceExtraMain;
		try(Input in=new Input(path)){
			byte[] header;
			do{header=in.reader.nextLine();}while(header!=null && header.length>0 && header[0]=='#');
			require(header!=null, "CM invocation profile must have its actual scalar schema");
			p.set(header);
			require(p.terms()>=2 && p.termEquals("member", 0) && p.termEquals("calls", 1), "CM profile must identify actual per-member invocation totals");
			int scoreColumn=-1;
			if(scoreGuard){
				for(int i=2; i<p.terms(); i++){
					if(p.termEquals("rejected_score", i)){
						require(scoreColumn<0, "The independent profile must declare its rejected-score counter once");
						scoreColumn=i;
					}
				}
				require(scoreColumn>=0, "Score-guard reporting requires the actual rejected_score profile counter");
			}
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final int m=member(p.parseString(0));
				require((memberMask&(1<<m))!=0 && !seen[m] && p.parseLong(1)==observed[m],
					"Captured M counts must equal independent CM invocations for exactly the declared members");
				if(scoreGuard){require(p.parseLong(scoreColumn)==held[m], "All SCORE_REJECTED requests must reconcile with the independent caller counter");}
				seen[m]=true;
			}
		}
		for(int m=0; m<seen.length; m++){
			require(seen[m]==((memberMask&(1<<m))!=0) || acceptanceScores!=null && !seen[m] && observed[m]==0,
				"Every declared endpoint scorer needs a profile row; explicit acceptance costs account for score-only minor scorers");
		}
	}
	private void readReferences(final String prefix){
		try(Input in=new Input(prefix+".references.tsv")){
			RrnaResourceIO.header(in.reader, "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored", prefix);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==7, "Normalized reference schema changed");
				final Genome genome=input();
				final Reference ref=new Reference(p.parseString(1), p.parseString(2), p.parseInt(3), p.parseInt(4), strand(5), bool(6));
				for(final Reference previous:genome.refs){require(!previous.id.equals(ref.id), "Normalized reference IDs must be unique within each input");}
				genome.refs.add(ref);
			}
		}
	}
	private void readQueries(final String prefix){
		try(Input in=new Input(prefix+".queries.tsv")){
			RrnaResourceIO.header(in.reader, "input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id", prefix);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==10, "Normalized query schema changed");
				final Genome genome=input();
				final Query query=new Query(p.parseString(1), p.parseInt(2), p.parseInt(3), strand(4), memberId(p.parseString(5)),
					CompareCobalaminBenchmarkGrades.modelId(p.parseString(6)), p.parseString(7), p.parseString(9));
				final int outcome=query.outcome.equals("TP") ? 0 : query.outcome.equals("MASKED") ? 1
					: query.outcome.equals("UNMATCHED_OTHER") || query.outcome.equals("UNMATCHED_DUPLICATE") ? 2 : -1;
				final boolean eligible=bool(8);
				require(outcome>=0 && query.start>0 && query.stop>=query.start && (outcome!=0 || eligible)
					&& (outcome==0)==!(query.reference.equals(".") || query.reference.equals("NA")), "Queries need explicit known outcomes, valid bounds and reference IDs only for recorded TPs");
				genome.outcomes[outcome]++;
				if(outcome==0){
					Reference assigned=null;
					for(final Reference ref:genome.refs){if(ref.id.equals(query.reference)){assigned=ref;}}
					// The inclusive exporter normalizes TP_INEXACT to TP after
					// freezing strict matches. Its residual pairs need positive
					// overlap, not the separate original-oracle half-overlap rule.
					require(assigned!=null && assigned.scored && !assigned.consumed && assigned.strand==query.strand
						&& assigned.seqid.equals(query.seqid) && Math.min(assigned.stop, query.stop)>=Math.max(assigned.start, query.start),
						"A recorded strict/inclusive TP must bind one unconsumed scored reference on the same contig/strand with positive overlap; never rematch it");
					assigned.consumed=true;
					query.assigned=assigned;
				}
				genome.queries.add(query);
			}
		}
		int n=0;
		try(Input in=new Input(prefix+".genomes.tsv")){
			RrnaResourceIO.header(in.reader, "input_id\tcalls\tmatched\tmasked\tunmatched", prefix);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final Genome genome=input();
				require(p.terms()==5 && !genome.graded && genome.queries.size()==p.parseInt(1)
					&& p.parseInt(1)==p.parseInt(2)+p.parseInt(3)+p.parseInt(4)
					&& genome.outcomes[0]==p.parseInt(2) && genome.outcomes[1]==p.parseInt(3) && genome.outcomes[2]==p.parseInt(4),
					"Complete graded query populations must conserve matched, masked and unmatched outcomes for each input");
				genome.graded=true;
				n++;
			}
		}
		require(n==expectedInputs, "Grade population must include every diagnostic input");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fixed Oracle          ----------------*/
	/*--------------------------------------------------------------*/

	private void analyze(final Call call, final String arm, final Stats stats, final Output output){
		if(stats!=null){stats.counts[0]++;}
		if(call.selected){
			if(stats!=null){stats.counts[1]++;}
			linkSelected(call);
		}
		analyzeEndpoint(call, arm, stats, output);
	}
	/** Exact saved query linkage is shared by endpoint and acceptance reports. */
	private void linkSelected(final Call call){
		assert(call.selected) : "Only an emitted trace invocation can consume a normalized query link";
		Query found=null;
		for(final Query query:call.genome.queries){
			if(query.strand==call.strand && query.start==call.finalStart && query.stop==call.finalStop
				&& query.member==member(call.member) && query.seqid.equals(call.seqid) && query.model.equals(call.model)){
				require(found==null && !query.linked, "A selected invocation must bind one distinct emitted query exactly");
				found=query;
			}
		}
		require(found!=null, "Selected CM invocation has no exact source/model/coordinate grade join");
		found.linked=true;
		call.query=found;
	}
	private void analyzeEndpoint(final Call call, final String arm, final Stats stats, final Output output){
		if(!call.measured){
			if(stats!=null){stats.counts[2]++;}
			require(!call.invoked && !call.known, "Every known CM boundary invocation or fixed-end skip must contain its actual proposal record");
			if(output!=null){writeCall(call, arm, null, "NO_CM_RECORD", false, output);}
			return;
		}
		if(stats!=null){
			stats.counts[3]++;
			if(call.unavailable){stats.counts[4]++;}
			if(call.contacts32!=0){stats.counts[5]++;}
			if(call.width==-1){stats.counts[6]++;}
			if(call.scored && call.score==Float.NEGATIVE_INFINITY){stats.counts[7]++;}
		}
		final int start=genomicStart(call, call.start, call.stop), stop=genomicStop(call, call.start, call.stop);
		Reference oracle=null;
		int eligible=0;
		for(final Reference ref:call.genome.refs){
			if(ref.strand==call.strand && ref.seqid.equals(call.seqid)){
				final int overlap=Math.max(0, Math.min(stop, ref.stop)-Math.max(start, ref.start)+1);
				if(2L*overlap>=Math.max(stop-start+1, ref.stop-ref.start+1)){
					if(ref.scored){
						oracle=ref;
						eligible++;
					}
					else{call.maskedOracles++;}
				}
			}
		}
		final String status=eligible==0 ? "NO_ORIGINAL_ORACLE" : eligible>1 ? "AMBIGUOUS_ORIGINAL_ORACLE" : "UNIQUE_ORIGINAL_ORACLE";
		if(stats!=null){stats.counts[eligible==0 ? 8 : eligible>1 ? 9 : 10]++;}
		if(eligible!=1){oracle=null;}
		call.oracle=oracle;
		call.oracleStatus=status;
		final boolean finite=call.proposedStart>=0;
		final boolean applicable=call.decision.equals("APPLIED") || call.decision.equals("UNCHANGED");
		require(!finite || call.scored && Float.isFinite(call.score), "Finite proposed coordinates require a finite actual attempted score");
		require(!applicable || finite, "An applied CM proposal cannot have missing coordinates");
		if(oracle!=null && finite){
			call.old5=call.strand==0 ? (long)start-oracle.start : (long)oracle.stop-stop;
			call.old3=call.strand==0 ? (long)stop-oracle.stop : (long)oracle.start-start;
			final int a=genomicStart(call, call.proposedStart, call.proposedStop), b=genomicStop(call, call.proposedStart, call.proposedStop);
			call.new5=call.strand==0 ? (long)a-oracle.start : (long)oracle.stop-b;
			call.new3=call.strand==0 ? (long)b-oracle.stop : (long)oracle.start-a;
			if(stats!=null){stats.observe(call, applicable);}
		}
		if(output!=null){writeCall(call, arm, oracle, status, applicable && oracle!=null && finite, output);}
	}
	private void writeCall(final Call c, final String arm, final Reference ref, final String oracle, final boolean testable, final Output out){
		row.clear().append(arm).tab().append(c.genome.id).tab().append(c.ordinal).tab().append(c.worker).tab().append(c.id)
			.tab().append(c.seqid).tab().append(c.strand==0 ? '+' : '-').tab().append(c.member).tab().append(c.model)
			.tab().append(c.selected).tab().append(c.query==null ? "NOT_EMITTED" : c.query.outcome)
			.tab().append(c.query==null ? "." : c.query.reference).tab().append(genomicStart(c, c.start, c.stop))
			.tab().append(genomicStop(c, c.start, c.stop)).tab().append(c.stop-c.start+1).tab();
		if(c.proposedStart>=0){row.append(genomicStart(c, c.proposedStart, c.proposedStop));}
		row.tab();
		if(c.proposedStart>=0){row.append(genomicStop(c, c.proposedStart, c.proposedStop));}
		row.tab();
		if(c.proposedStart>=0){row.append(c.proposedStop-c.proposedStart+1);}
		row.tab().append(genomicStart(c, c.effectiveStart, c.effectiveStop)).tab().append(genomicStop(c, c.effectiveStart, c.effectiveStop)).tab();
		if(c.measured){row.append(c.consensusLength).tab().append(c.cmLength);}else{row.tab();}
		row.tab();
		if(c.proposedStart>=0){row.append(c.ratio(0), 8);}
		row.tab();
		if(c.proposedStart>=0){row.append(c.ratio(1), 8);}
		row.tab();
		if(c.proposedStart>=0){row.append(c.ratio(2), 8);}
		row.tab();
		if(c.scored){row.appendSlow(c.score);}
		row.tab();
		if(c.measured){row.appendSlow(c.ga);}
		row.tab();
		if(c.scored && Float.isFinite(c.score)){row.append(c.score-c.ga, 8);}
		row.tab();
		if(c.measured){row.append(c.width).tab().append(c.contacts16).tab().append(c.contacts32).tab().append(c.anchor);}
		else{row.append("\t\t\t");}
		row.tab();
		if(c.anchorOrigin>=0){row.append(c.anchorOrigin);}
		row.tab();
		if(c.measured){row.append(c.unavailable);}
		row.tab().append(c.measured ? c.decision : "NO_CM_RECORD")
			.tab().append(oracle).tab().append(ref==null ? "." : ref.id).tab().append(testable);
		if(ref!=null && c.proposedStart>=0){
			row.tab().append(c.old5).tab().append(c.new5).tab().append(c.old3).tab().append(c.new3)
				.tab().append(Math.abs(c.new5)>Math.abs(c.old5)).tab().append(Math.abs(c.new3)>Math.abs(c.old3))
				.tab().append(endClass(c.old5)).tab().append(endClass(c.new5)).tab().append(endClass(c.old3)).tab().append(endClass(c.new3));
		}else{row.append("\t\t\t\t\t\t\t\t\t\t");}
		if(scoreGuard){row.tab().append(c.maskedOracles);}
		out.writer.print(row.nl());
	}
	private Genome ordinal(final int ordinal){
		require(ordinal>=0 && ordinal<genomes.size(), "Trace ordinal leaves declared panel");
		return genomes.get(ordinal);
	}
	private Genome input(){
		probe.set(p, 0);
		final Genome genome=byInput.get(probe);
		require(genome!=null, "Grade row names an input absent from the exact diagnostic panel");
		return genome;
	}
	private void traceKey(){
		key.clear().append(p.parseInt(1)).tab().append(p.parseInt(2)).tab().append(p.parseLong(3));
		probe.set(key.array, 0, key.length());
	}
	private boolean bit(final int column){
		require(p.termEquals("0", column) || p.termEquals("1", column), "Trace flags must explicitly distinguish measured false/true");
		return p.termEquals("1", column);
	}
	private boolean bool(final int column){
		require(p.termEquals("false", column) || p.termEquals("true", column), "Normalized grade flags require explicit booleans");
		return p.termEquals("true", column);
	}
	private int strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Normalized strand must be explicitly plus or minus");
		return p.termEquals("-", column) ? 1 : 0;
	}
	private static int member(final String name){
		for(int i=0; i<MEMBER_NAMES.length; i++){if(name.equals(MEMBER_NAMES[i])){return i;}}
		throw new IllegalArgumentException("Unknown CM profile member: "+name);
	}
	private static int memberId(final String name){
		for(int i=0; i<MEMBERS.length; i++){if(name.equals(MEMBERS[i])){return i;}}
		throw new IllegalArgumentException("Use distinct RF00174/RF01482/RF01689 member IDs: "+name);
	}
	private static int genomicStart(final Call c, final int start, final int stop){return c.strand==0 ? start+1 : c.length-stop;}
	private static int genomicStop(final Call c, final int start, final int stop){return c.strand==0 ? stop+1 : c.length-start;}
	private static String endClass(final long error){return error==0 ? "EXACT" : Math.abs(error)<=3 ? "WITHIN3" : "BEYOND3";}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Accumulation         ----------------*/
	/*--------------------------------------------------------------*/

	/** The first three scopes partition request selection; the recorded-TP scope
	 * is a separate fixed-assignment view, never an additional set of detections.
	 */
	private static final class GuardCounts{
		void add(final Call call, final Reference ref, final String status){
			assert(call!=null && call.proposedStart>=0 && call.decision.equals("SCORE_REJECTED")) :
				"Guard summaries consume only the actual finite score holds checked by processScoreGuard";
			held++;
			if(!call.selected){notEmitted++;}
			else if(call.query.outcome.equals("TP")){selectedTP++;}
			else if(call.query.outcome.equals("MASKED")){selectedMasked++;}
			else{selectedUnmatched++;}
			if(ref==null){
				if(status.equals("AMBIGUOUS_ORIGINAL_ORACLE")){ambiguous++;}
				else{
					noReference++;
					if(call.maskedOracles>0){maskedOnly++;}
				}
				return;
			}
			scorable++;
			final int a=genomicStart(call, call.start, call.stop), b=genomicStop(call, call.start, call.stop);
			final int x=genomicStart(call, call.proposedStart, call.proposedStop), y=genomicStop(call, call.proposedStart, call.proposedStop);
			final long old5=call.strand==0 ? (long)a-ref.start : (long)ref.stop-b;
			final long old3=call.strand==0 ? (long)b-ref.stop : (long)ref.start-a;
			final long new5=call.strand==0 ? (long)x-ref.start : (long)ref.stop-y;
			final long new3=call.strand==0 ? (long)y-ref.stop : (long)ref.start-x;
			final boolean oldExact=old5==0 && old3==0, newExact=new5==0 && new3==0;
			final boolean newNear=Math.abs(new5)<=3 && Math.abs(new3)<=3;
			if(oldExact){originalBothExact++;}else{originalInexact++;}
			if(newExact){proposalBothExact++;}
			if(newNear){proposalBothWithin3++;}
			if(oldExact && !newExact){originalExactSpoiled++;}
			if(!oldExact){
				if(newExact){inexactToExact++;}
				else if(newNear){inexactToWithin3Nonexact++;}
				else{inexactToBeyond3++;}
				if(Math.abs(new5)>Math.abs(old5) || Math.abs(new3)>Math.abs(old3)){inexactAnyEndWorse++;}
			}
			if(Math.abs(new5)<Math.abs(old5)){better5++;}
			if(Math.abs(new5)>Math.abs(old5)){worse5++;}
			if(Math.abs(new3)<Math.abs(old3)){better3++;}
			if(Math.abs(new3)>Math.abs(old3)){worse3++;}
			if(old5==0){originalExact5++;}
			if(old3==0){originalExact3++;}
			if(new5==0){proposalExact5++;}
			if(new3==0){proposalExact3++;}
			if(Math.abs(old5)<=3){originalWithin5++;}
			if(Math.abs(old3)<=3){originalWithin3++;}
			if(Math.abs(new5)<=3){proposalWithin5++;}
			if(Math.abs(new3)<=3){proposalWithin3++;}
		}
		void write(final String arm, final String scope, final Output output){
			require(held==scorable+noReference+ambiguous && held==selectedTP+selectedMasked+selectedUnmatched+notEmitted
				&& scorable==originalBothExact+originalInexact && originalInexact==inexactToExact+inexactToWithin3Nonexact+inexactToBeyond3,
				"Every held request retains one oracle and selection disposition; inexact proposal classes must conserve the fixed original-inexact population");
			final ByteBuilder b=new ByteBuilder();
			output.writer.print(b.append(arm).tab().append(scope).tab().append(held).tab().append(scorable).tab().append(noReference)
				.tab().append(ambiguous).tab().append(maskedOnly).tab().append(selectedTP).tab().append(selectedMasked).tab().append(selectedUnmatched)
				.tab().append(notEmitted).tab().append(originalBothExact).tab().append(originalInexact).tab().append(proposalBothExact)
				.tab().append(proposalBothWithin3).tab().append(originalExactSpoiled).tab().append(inexactToExact).tab().append(inexactToWithin3Nonexact)
				.tab().append(inexactToBeyond3).tab().append(inexactAnyEndWorse).tab().append(better5).tab().append(worse5).tab().append(better3).tab().append(worse3)
				.tab().append(originalExact5).tab().append(originalExact3).tab().append(proposalExact5).tab().append(proposalExact3)
				.tab().append(originalWithin5).tab().append(originalWithin3).tab().append(proposalWithin5).tab().append(proposalWithin3)
				.append("\tFIXED_REFERENCE_REQUEST_COUNTERFACTUAL_NO_PATH_REPLAY_NOT_EMITTED_GAINS\n"));
		}
		long held, scorable, noReference, ambiguous, maskedOnly, selectedTP, selectedMasked, selectedUnmatched, notEmitted;
		long originalBothExact, originalInexact, proposalBothExact, proposalBothWithin3, originalExactSpoiled;
		long inexactToExact, inexactToWithin3Nonexact, inexactToBeyond3, inexactAnyEndWorse, better5, worse5, better3, worse3;
		long originalExact5, originalExact3, proposalExact5, proposalExact3, originalWithin5, originalWithin3, proposalWithin5, proposalWithin3;
		static final String HEADER="arm\tscope\theld\tscorable\tno_scored_reference\tambiguous_scored_reference\tno_scored_reference_masked_overlap"
			+"\tselected_recorded_tp\tselected_masked\tselected_unmatched\tnot_emitted\toriginal_both_exact\toriginal_inexact"
			+"\tproposal_both_exact\tproposal_both_within3\toriginal_exact_spoiled\tinexact_to_both_exact\tinexact_to_both_within3_nonexact"
			+"\tinexact_to_any_end_beyond3\tinexact_any_end_worse\tbetter5\tworse5\tbetter3\tworse3\toriginal_exact5\toriginal_exact3"
			+"\tproposal_exact5\tproposal_exact3\toriginal_within3_5\toriginal_within3_3\tproposal_within3_5\tproposal_within3_3\tqualification";
	}

	private static final class Stats{
		Stats(){
			for(int population=0; population<2; population++){
				for(int axis=0; axis<AXES.length; axis++){
					bins[population][axis]=new Bucket[axis<3 ? RATIO_BINS.length+1 : axis==3 ? SCORE_BINS.length+1 : axis==4 ? 3 : axis==7 ? 4 : 2];
					for(int i=0; i<bins[population][axis].length; i++){bins[population][axis][i]=new Bucket();}
				}
			}
		}
		void observe(final Call c, final boolean applicable){
			assert(c.proposedStart>=0) : "Endpoint error bins use finite recorded proposals and an existing fixed oracle";
			final int population=applicable ? 0 : 1;
			for(int axis=0; axis<AXES.length; axis++){
				final int index=axis<3 ? bin(c.ratio(axis), RATIO_BINS) : axis==3 ? bin(c.score-c.ga, SCORE_BINS)
					: axis==4 ? c.member.equals("cobalamin_rf00174") ? 0 : c.member.equals("cobalamin_rf01482") ? 1 : 2
					: axis==5 ? c.contacts32==0 ? 0 : 1 : axis==6 ? c.unavailable ? 1 : 0 : c.anchor;
				bins[population][axis][index].add(c);
			}
			if(!applicable){return;}
			counts[11]++;
			final boolean bad=Math.abs(c.new5)>Math.abs(c.old5) || Math.abs(c.new3)>Math.abs(c.old3);
			if(bad){counts[12]++;}
			final int exact=(c.new5==0 ? 1 : 0)+(c.new3==0 ? 1 : 0), original=(c.old5==0 ? 1 : 0)+(c.old3==0 ? 1 : 0);
			for(int i=0; i<TESTS.length; i++){
				final boolean reject=i<CUTS.length ? c.ratio(0)<CUTS[i] : i<2*CUTS.length ? c.ratio(2)<CUTS[i-CUTS.length]
					: i==2*CUTS.length ? c.score<c.ga : i==2*CUTS.length+1 ? c.contacts32!=0 : c.unavailable;
				final long[] t=tests[i];
				t[reject ? bad ? 0 : 1 : bad ? 2 : 3]++;
				t[reject ? 5 : 4]+=exact;
				if(reject){t[6]+=original;}
			}
		}
		void add(final Stats other){
			assert(other!=null && other!=this) : "Pooled statistics merge independent arm populations once";
			for(int i=0; i<counts.length; i++){counts[i]+=other.counts[i];}
			for(int i=0; i<tests.length; i++){for(int j=0; j<7; j++){tests[i][j]+=other.tests[i][j];}}
			for(int p=0; p<2; p++){
				for(int a=0; a<AXES.length; a++){
					for(int i=0; i<bins[p][a].length; i++){bins[p][a][i].merge(other.bins[p][a][i]);}
				}
			}
		}
		void write(final String arm, final Output binOutput, final Output testOutput, final Output summary){
			final ByteBuilder b=new ByteBuilder();
			for(int i=0; i<counts.length; i++){
				summary.writer.print(b.clear().append(arm).tab().append(MEASURES[i]).tab().append(counts[i])
					.append("\tACTUAL_INVOCATIONS_SELECTED_QUERY_JOINS_FIXED_ORIGINAL_ORACLE\n"));
			}
			for(int p=0; p<2; p++){
				for(int axis=0; axis<AXES.length; axis++){
					for(int i=0; i<bins[p][axis].length; i++){
						final Bucket bucket=bins[p][axis][i];
						if(bucket.n==0){continue;}
						b.clear().append(arm).tab().append(p==0 ? "BASE_APPLICABLE_UNIQUE_ORIGINAL_ORACLE" : "EXISTING_HOLD_UNIQUE_ORIGINAL_ORACLE").tab().append(AXES[axis]).tab();
						if(axis<4){
							final double[] edges=axis<3 ? RATIO_BINS : SCORE_BINS;
							if(i==0){b.append("below_").append(edges[0], 4);}
							else if(i==edges.length){b.append("at_least_").append(edges[i-1], 4);}
							else{b.append(edges[i-1], 4).append("_to_below_").append(edges[i], 4);}
						}else if(axis==4){b.append(MEMBERS[i]);}else{b.append(i);}
						binOutput.writer.print(b.tab().append(bucket.n).tab().append(bucket.bad).tab().append(bucket.worse5).tab().append(bucket.worse3)
							.tab().append(bucket.oldExact).tab().append(bucket.newExact).tab().append(bucket.oldWithin).tab().append(bucket.newWithin).tab().append(bucket.selected).nl());
					}
				}
			}
			int best=-1;
			double bestJ=-Double.MAX_VALUE;
			for(int i=0; i<tests.length; i++){
				final long[] t=tests[i];
				if(t[0]+t[2]==0 || t[1]+t[3]==0){continue;}
				final double j=t[0]/(double)(t[0]+t[2])-t[1]/(double)(t[1]+t[3]);
				if(j>bestJ || j==bestJ && t[5]<tests[best][5]){
					best=i;
					bestJ=j;
				}
			}
			for(int i=0; i<tests.length; i++){
				final long[] t=tests[i];
				b.clear().append(arm).tab().append(TESTS[i]).tab().append(t[0]+t[1]+t[2]+t[3]);
				for(final long v:t){b.tab().append(v);}
				b.tab();
				if(t[0]+t[2]>0 && t[1]+t[3]>0){b.append(t[0]/(double)(t[0]+t[2])-t[1]/(double)(t[1]+t[3]), 8);}
				testOutput.writer.print(b.tab().append(i==best).append("\tPROPOSAL_ONLY_NO_PATH_REPLAY_NO_DETECTION_PREDICTION\n"));
			}
		}
		final long[] counts=new long[MEASURES.length];
		final long[][] tests=new long[TESTS.length][7];
		final Bucket[][][] bins=new Bucket[2][AXES.length][];
	}
	private static final class Bucket{
		void add(final Call c){
			assert(c.proposedStart>=0 && c.cmLength>0) : "Endpoint buckets consume finite recorded proposals with validated model dimensions";
			n++;
			final boolean a=Math.abs(c.new5)>Math.abs(c.old5), b=Math.abs(c.new3)>Math.abs(c.old3);
			if(a || b){bad++;}
			if(a){worse5++;}
			if(b){worse3++;}
			oldExact+=(c.old5==0 ? 1 : 0)+(c.old3==0 ? 1 : 0);
			newExact+=(c.new5==0 ? 1 : 0)+(c.new3==0 ? 1 : 0);
			oldWithin+=(Math.abs(c.old5)<=3 ? 1 : 0)+(Math.abs(c.old3)<=3 ? 1 : 0);
			newWithin+=(Math.abs(c.new5)<=3 ? 1 : 0)+(Math.abs(c.new3)<=3 ? 1 : 0);
			if(c.selected){selected++;}
		}
		void merge(final Bucket b){
			assert(b!=null && b!=this) : "Pooled bins merge distinct arm populations once";
			n+=b.n;
			bad+=b.bad;
			worse5+=b.worse5;
			worse3+=b.worse3;
			oldExact+=b.oldExact;
			newExact+=b.newExact;
			oldWithin+=b.oldWithin;
			newWithin+=b.newWithin;
			selected+=b.selected;
		}
		long n, bad, worse5, worse3, oldExact, newExact, oldWithin, newWithin, selected;
	}
	private static int bin(final double value, final double[] edges){
		assert(Double.isFinite(value)) : "Proposal bins require an observed finite ratio or score margin";
		int i=0;
		while(i<edges.length && value>=edges[i]){i++;}
		return i;
	}
	private static String[] tests(){
		final String[] names=new String[2*CUTS.length+3];
		for(int i=0; i<CUTS.length; i++){
			names[i]="proposed_original_ratio_below_"+CUTS[i];
			names[i+CUTS.length]="proposed_cm_ratio_below_"+CUTS[i];
		}
		names[2*CUTS.length]="score_below_GA";
		names[2*CUTS.length+1]="width32_contact";
		names[2*CUTS.length+2]="anchor_unavailable";
		return names;
	}

	/*--------------------------------------------------------------*/
	/*----------------             State            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		Genome(final String id_, final String fasta_, final String assembly_, final String referencePath_){
			id=id_;
			fasta=fasta_;
			assembly=assembly_;
			referencePath=referencePath_;
		}
		final String id, fasta, assembly, referencePath;
		int group=-1;
		final ArrayList<Reference> refs=new ArrayList<Reference>();
		final ArrayList<Query> queries=new ArrayList<Query>();
		boolean seen, graded;
		final int[] outcomes=new int[3];
	}
	private static final class Contig{
		Contig(final String name_, final int length_){
			name=name_;
			length=length_;
		}
		final String name;
		final int length;
		int calls, returned, selected, unknown;
		boolean closed;
	}
	private static final class Reference{
		Reference(final String id_, final String seqid_, final int start_, final int stop_, final int strand_, final boolean scored_){
			id=id_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			scored=scored_;
			require(start>0 && stop>=start, "Reference coordinates must be positive inclusive intervals");
		}
		final String id, seqid;
		final int start, stop, strand;
		final boolean scored;
		boolean consumed;
		boolean nativeSeen;
		NcrnaReferenceGff.Scores nativeScores;
	}
	private static final class Query{
		Query(final String seqid_, final int start_, final int stop_, final int strand_, final int member_, final String model_, final String outcome_, final String reference_){
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			member=member_;
			model=model_;
			outcome=outcome_;
			reference=reference_;
		}
		final String seqid, model, outcome, reference;
		final int start, stop, strand, member;
		boolean linked, rawLinked, residualInexact;
		Reference assigned;
	}
	private static final class Call{
		Call(final Genome genome_, final int ordinal_, final int worker_, final long id_, final String seqid_, final int strand_,
				final String member_, final String model_, final int length_, final int start_, final int stop_){
			genome=genome_;
			ordinal=ordinal_;
			worker=worker_;
			id=id_;
			seqid=seqid_;
			strand=strand_;
			member=member_;
			model=model_;
			length=length_;
			start=start_;
			stop=stop_;
			require(start>=0 && stop>=start && stop<length && (strand==0 || strand==1) && NcrnaFamilyRegistry.isCobalamin(member), "Candidate coordinates/member must be actual caller metadata");
		}
		double ratio(final int axis){return (proposedStop-proposedStart+1)/(double)(axis==0 ? stop-start+1 : axis==1 ? consensusLength : cmLength);}
		final Genome genome;
		final int ordinal, worker, length, start, stop, strand;
		final long id;
		long captureEventIndex=-1;
		final String seqid, member, model;
		int proposedStart=-1, proposedStop=-1, effectiveStart, effectiveStop, finalStart, finalStop, consensusLength, cmLength, width, contacts16, contacts32, anchor, anchorOrigin=-1, maskedOracles;
		boolean measured, scored, known, invoked, selected, returned, unavailable, captured;
		float score, ga;
		float acceptedIdentity=Float.NaN, endpointIdentity=Float.NaN;
		double nativeBits=Double.NaN, nativeGa=Double.NaN;
		String nativeStatus;
		String decision="NO_CM_RECORD";
		long old5, old3, new5, new3;
		Query query;
		Reference oracle;
		String oracleStatus;
	}
	private String bindings, out, audit;
	private String[] traces, grades, names, profiles, rawGrades, acceptanceScores, acceptanceCosts;
	private long acceptanceExtraMain;
	private int expectedInputs=1000;
	private int memberMask=7;
	private boolean scoreGuard, acceptance, nativeCm=true, eventAssignments;
	public static final String EVENT_ASSIGNMENT_HEADER="capture_event_index\tinput_id\tresolved_group\tinput_ordinal\tworker\ttrace_id\tcontig\tstrand\tstart\tstop\temitted_member\tmodel\toutcome\treference_id\tidentity_bits\town_score_bits\town_ga_bits";
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder(), row=new ByteBuilder();
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Call> order=new ArrayList<Call>();
	private final ArrayList<Contig> contigs=new ArrayList<Contig>();
	private final ObjectMap<Target, Genome> byInput=new ObjectMap<Target, Genome>(16, Target.class, Genome.class);
	private final ObjectMap<Target, Call> byTrace=new ObjectMap<Target, Call>(16, Target.class, Call.class);
	private final ObjectMap<Target, Contig> byContig=new ObjectMap<Target, Contig>(16, Target.class, Contig.class);
	private static final String CALL_HEADER="arm\tinput_id\tinput_ordinal0\tworker\ttrace_id\tcontig\tstrand\tmember\tmodel\tselected\tgrade_outcome\tgrade_reference"
		+"\toriginal_start1\toriginal_stop1\toriginal_length\tproposed_start1\tproposed_stop1\tproposed_length\teffective_start1\teffective_stop1\tconsensus_length\tcm_length"
		+"\tproposed_original_ratio\tproposed_consensus_ratio\tproposed_cm_ratio\tcm_bits\tGA_bits\tcm_minus_GA\twidth\tcontacts16\tcontacts32\tanchor_kind\tanchor_origin0\tanchor_unavailable"
		+"\tdecision\toracle_status\toracle_reference\tguard_testable\toriginal_d5\tproposed_d5\toriginal_d3\tproposed_d3\tworse5\tworse3"
		+"\toriginal_5prime_class\tproposed_5prime_class\toriginal_3prime_class\tproposed_3prime_class";
	private static final double[] RATIO_BINS={.5, .67, .75, .8, .9, 1, 1.1, 1.25};
	private static final double[] SCORE_BINS={-20, -10, 0, 10, 20};
	private static final double[] CUTS={.5, .6, .67, .75, .8, .9, .95};
	private static final String[] TESTS=tests();
	private static final String[] AXES={"proposed_original_ratio", "proposed_consensus_ratio", "proposed_cm_ratio", "score_minus_GA", "member", "width32_contact", "anchor_unavailable", "anchor_kind"};
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] ACCEPTANCE_AXES={"accepted_identity", "cm_minus_ga", "cm_over_ga"};
	private static final double[] OFFSET_QUANTILES={0, .05, .5, .95, 1};
	private static final String[] MEMBER_NAMES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	private static final String[] MEASURES={"candidates", "selected_exact_query_joins", "missing_cm_record", "cm_invocations", "anchor_unavailable", "width32_contact", "exact_search", "impossible",
		"no_original_oracle", "ambiguous_original_oracle", "unique_original_oracle", "applicable_unique_original_oracle", "worse_any_applicable"};
}
