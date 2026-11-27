package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
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

/** Compares normalized, already graded cobalamin benchmark views.
 * Uses template.A_SampleByteFile parsing and native retained-key maps. No caller
 * or matching algorithm runs. Same-strand physical overlap is reported beside
 * existing one-to-one outcomes; it never changes a graded TP, FP, FN or mask.
 * Genome call presence and scored-reference detection are separate measures.
 * Optional radius= adds paired endpoint errors on references matched in both
 * inputs. Full paired RMS is primary; the same OFF-selected within-radius
 * subset is reported separately. Member cohorts stay fixed by the OFF label.
 * Optional trace= and bindings= replay the four predeclared confidence margins
 * on the SAME paired reference cohort. Exact final-call identities must join
 * saved selected candidates; this is conditional diagnostic evidence, not a
 * replacement caller or matcher. Unknown eligibility makes selection unavailable.
 *
 * Optional subset= contains three selected consensus FASTAs from the SAME
 * original generation as in1. Name membership then reports whether a previous
 * winning model survived subset deletion; it cannot establish cross-generation
 * model equivalence. All retained rows are report state, not transient per-field
 * parsing objects. A single thread owns the complete comparison and streams.
 * @author Raiden
 */
public final class CompareCobalaminBenchmarkGrades{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CompareCobalaminBenchmarkGrades.class, false);
		try{
			for(final String arg:pp.args){
				if(arg.equalsIgnoreCase("mode=report")){
					new CobalaminBoundaryReport(pp.args).process();
					return;
				}
			}
			new CompareCobalaminBenchmarkGrades(pp.args).process();
		}
		finally{Shared.closeStream(pp.outstream);}
	}
	private CompareCobalaminBenchmarkGrades(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("subset")){subset=b.split(",");}
			else if(a.equals("changesonly")){changesOnly=Parse.parseBoolean(b);}
			else if(a.equals("radius")){radius=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("trace")){trace=b;}
			else if(a.equals("bindings")){bindings=b;}
			else if(a.equals("replaymargins")){
				require(b!=null && !b.isEmpty(), "replaymargins requires an explicit comma-separated numeric list");
				final String[] values=b.split(",", -1);
				replayMargins=new float[values.length];
				for(int i=0; i<values.length; i++){replayMargins[i]=NcrnaNumbers.floatValue(a, values[i]);}
			}
			else{require(parser.parse(arg, a, b), "Unknown normalized-comparison option: "+arg);}
		}
		prefixes[0]=parser.in1;
		prefixes[1]=parser.in2;
		out=parser.out1;
		require(prefixes[0]!=null && prefixes[1]!=null && out!=null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Supply complete in1/in2 normalized grade prefixes and fresh out=prefix");
		for(final String prefix:prefixes){
			require(Tools.testInputFiles(false, true, prefix+".genomes.tsv", prefix+".queries.tsv", prefix+".references.tsv", prefix+".witnesses.tsv"), "Both normalized grades must be complete and readable");
		}
		require(Tools.testOutputFiles(false, false, false, out+".summary.tsv", out+".genomes.tsv", out+".loci.tsv", out+".overlaps.tsv", out+".lengths.tsv", out+".call_lengths.tsv"), "Comparison outputs must be fresh");
		require(radius>=-1, "Optional radius must be nonnegative, or omitted for unchanged historical reports");
		if(radius>=0){
			require(Tools.testOutputFiles(false, false, false, out+".paired_endpoints.tsv", out+".paired_near_endpoints.tsv",
				out+".paired_population.tsv", out+".paired_loci.tsv", out+".paired_details.tsv"), "Paired endpoint outputs must be fresh");
			pairedStats=new PairStats[4][2];
			for(int cohort=0; cohort<4; cohort++){
				for(int end=0; end<2; end++){pairedStats[cohort][end]=new PairStats();}
			}
		}
		require((trace==null)==(bindings==null), "Replay requires both trace= and its exact source bindings=");
		require(trace!=null || replayMargins==null, "replaymargins cannot be silently ignored without trace=/bindings=");
		if(trace!=null){
			require(radius>0 && Tools.testInputFiles(false, true, trace, bindings), "Replay needs a positive radius and readable trace/source bindings");
			require(Tools.testOutputFiles(false, false, false, out+".replay_loci.tsv", out+".replay_endpoints.tsv",
				out+".replay_selection.tsv", out+".replay_accounting.tsv"), "Replay outputs must be fresh");
		}
		Shared.setThreads(1);
		require(!changesOnly || radius<0 && trace==null && subset==null,
			"changesonly uses saved outcomes alone; paired/replay/subset options would be ignored");
		if(changesOnly){CobalaminCallChangeAudit.checkOutputs(out);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reading           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		if(subset!=null){
			require(subset.length==3, "Subset membership requires all three SAME-GENERATION selected consensus files");
			for(final String path:subset){
				for(final String name:NcrnaSupport.loadLibrary(path).names){require(kept.add(name), "Subset model names must be unique across the clan");}
			}
		}
		for(int side=0; side<2; side++){
			readGenomes(side);
			readReferences(side);
			readQueries(side);
			readWitnesses(side);
		}
		validate();
		if(changesOnly){
			new CobalaminCallChangeAudit(out).write(genomes);
			return;
		}
		if(trace!=null){
			replay=new CobalaminBoundaryReplay(trace, bindings, radius, replayMargins);
			for(final Genome genome:genomes){
				replay.input(genome.id);
				for(final Query query:genome.after){
					query.trace=replay.query(genome.id, query.seqid, query.strand, query.member, query.model, query.start, query.stop);
				}
			}
			replay.linked();
		}
		write();
		if(replay!=null){replay.write(out);}
	}

	private void readGenomes(final int side){
		try(Input in=input(side, ".genomes.tsv", GENOMES)){
			int count=0;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5, "Normalized input schema changed");
				probe.set(p, 0);
				Genome genome=byInput.get(probe);
				if(side==0){
					require(genome==null, "Repeated baseline input");
					genome=new Genome(p.parseString(0));
					byInput.put(probe.copy(), genome);
					genomes.add(genome);
				}else{
					require(genome!=null && !genome.seen, "Candidate input is duplicated or absent from baseline");
					genome.seen=true;
				}
				for(int i=0; i<4; i++){
					genome.expected[side][i]=p.parseLong(i+1);
					require(genome.expected[side][i]>=0, "Genome counts cannot be negative");
				}
				require(genome.expected[side][0]==genome.expected[side][1]+genome.expected[side][2]+genome.expected[side][3], "Per-genome calls must partition into matched, masked and unmatched");
				count++;
			}
			require(count==genomes.size() && count>0, "Comparison input populations must agree, including negative genomes");
		}
	}
	private void readReferences(final int side){
		try(Input in=input(side, ".references.tsv", REFERENCES)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==7, "Normalized reference schema changed");
				final Genome genome=genome();
				referenceKey(genome, 1);
				Reference ref=byReference.get(probe);
				if(side==0){
					require(ref==null, "Repeated genome-qualified reference ID");
					ref=new Reference(genome, p.parseString(1), p.parseString(2), p.parseInt(3), p.parseInt(4), strand(5), flag(6));
					require(ref.start>0 && ref.stop>=ref.start, "References require ordered positive coordinates");
					byReference.put(probe.copy(), ref);
					genome.references.add(ref);
					references.add(ref);
				}else{
					require(ref!=null && !ref.seen && p.termEquals(ref.seqid, 2) && p.parseInt(3)==ref.start && p.parseInt(4)==ref.stop
						&& strand(5)==ref.strand && flag(6)==ref.scored, "Reference population, coordinates or disposition changed");
					ref.seen=true;
				}
			}
		}
	}
	private void readQueries(final int side){
		try(Input in=input(side, ".queries.tsv", QUERIES)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==10, "Normalized query schema changed");
				final Genome genome=genome();
				final Query query=new Query(p.parseString(1), p.parseInt(2), p.parseInt(3), strand(4),
					p.parseString(5), p.parseString(6), p.parseString(7), p.parseString(9));
				require(query.start>0 && query.stop>=query.start, "Query coordinates must be ordered positive intervals");
				member(query.member);
				require(query.outcome.equals("TP") || query.outcome.equals("MASKED") || query.outcome.equals("UNMATCHED_DUPLICATE")
					|| query.outcome.equals("UNMATCHED_OTHER"), "Unknown normalized query outcome");
				require(query.outcome.equals("TP")==!(query.reference.equals(".") || query.reference.equals("NA")), "Only a recorded TP may name an assigned reference");
				genome.queries(side).add(query);
				genome.actual[side][0]++;
				genome.actual[side][query.outcome.equals("TP") ? 1 : query.outcome.equals("MASKED") ? 2 : 3]++;
			}
		}
	}
	private void readWitnesses(final int side){
		try(Input in=input(side, ".witnesses.tsv", WITNESSES)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==12, "Normalized witness schema changed");
				final Genome genome=genome();
				referenceKey(genome, 1);
				final Reference ref=byReference.get(probe);
				require(ref!=null && ref.scored && ref.winners[side]==null && p.termEquals(ref.seqid, 4)
					&& p.parseInt(7)==ref.start && p.parseInt(8)==ref.stop && strand(9)==ref.strand, "TP witness must bind one original scored reference");
				Query winner=null;
				for(final Query query:genome.queries(side)){
					if(query.outcome.equals("TP") && query.reference.equals(ref.id)){
						require(winner==null && p.termEquals(query.member, 2) && p.termEquals(query.model, 3)
							&& p.parseInt(5)==query.start && p.parseInt(6)==query.stop && query.seqid.equals(ref.seqid)
							&& query.strand==ref.strand, "TP witness and query must identify the same unique recorded call");
						winner=query;
					}
				}
				require(winner!=null, "Recorded witness lacks its TP query");
				ref.winners[side]=winner;
			}
		}
	}
	private void validate(){
		for(final Genome genome:genomes){
			require(genome.seen, "Candidate omitted an input");
			for(int side=0; side<2; side++){
				for(int i=0; i<4; i++){require(genome.expected[side][i]==genome.actual[side][i], "Detailed query counts disagree with the complete genome table");}
				long matched=0;
				for(final Reference ref:genome.references){if(ref.winners[side]!=null){matched++;}}
				require(matched==genome.actual[side][1], "TP query and witness counts disagree");
			}
		}
		for(final Reference ref:references){require(ref.seen, "Candidate omitted a reference");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try(Output loci=new Output(out+".loci.tsv"); Output overlaps=new Output(out+".overlaps.tsv");
			Output inputs=new Output(out+".genomes.tsv"); Output lengths=new Output(out+".lengths.tsv");
			Output callLengths=new Output(out+".call_lengths.tsv");
			Output paired=radius<0 ? null : new Output(out+".paired_loci.tsv")){
			loci.writer.println("input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model");
			overlaps.writer.println("input_id\treference_id\tchange\tside\temitted_member\tmodel\tcall_start\tcall_stop\tcall_outcome\tassigned_reference\tintersection_bp\treference_bp\tcall_bp\tfraction_reference\tfraction_call\tgeometry");
			inputs.writer.println("input_id\tbaseline_calls\tcandidate_calls\tbaseline_matched\tcandidate_matched\traw_call_presence\tscored_detection_presence");
			lengths.writer.println(LENGTH_VERSION);
			lengths.writer.println("side\tinput_id\temitted_member\tmodel\toutcome\treference_id\treference_status\tcall_bp\treference_bp\tintersection_bp\tcall_to_reference_length\tqualification");
			callLengths.writer.println(CALL_LENGTH_VERSION);
			callLengths.writer.println("side\tinput_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tcall_bp\treference_overlap_pairs\tqualification");
			if(paired!=null){paired.writer.println(PAIRED_LOCI_HEADER);}
			for(final Reference ref:references){
				if(!ref.scored){continue;}
				final Query before=ref.winners[0], after=ref.winners[1];
				if(before!=null && after!=null){
					retained++;
					labels[member(before.member)][member(after.member)]++;
					if(!before.member.equals(after.member)){switches++;}
					if(paired!=null){paired(ref, before, after, paired);}
					continue;
				}
				if(before==null && after==null){continue;}
				final boolean loss=before!=null;
				if(loss){lost++;}else{gained++;}
				final String change=loss ? "LOST_GRADED_REFERENCE" : "GAINED_GRADED_REFERENCE";
				boolean present=false;
				for(final Query query:ref.genome.queries(loss ? 1 : 0)){if(intersection(ref, query)>0){present=true;}}
				if(loss && !present){absent++;}
				buffer.clear().append(ref.genome.id).tab().append(ref.id).tab().append(ref.seqid).tab().append(ref.start).tab().append(ref.stop)
					.tab().append(ref.strand).tab().append(change).tab().append(present ? "OVERLAPPING_CLAN_CALL_PRESENT" : "NO_OVERLAPPING_CLAN_CALL").tab();
				winner(buffer, before);
				buffer.tab().append(before==null || subset==null ? "NA" : kept.contains(modelId(before.model)) ? "true" : "false").tab();
				winner(buffer, after);
				loci.writer.print(buffer.nl());
				for(int side=0; side<2; side++){
					for(final Query query:ref.genome.queries(side)){
						final int bp=intersection(ref, query);
						if(bp<1){continue;}
						final int rl=ref.stop-ref.start+1, ql=query.stop-query.start+1;
						overlaps.writer.print(buffer.clear().append(ref.genome.id).tab().append(ref.id).tab().append(change).tab().append(side)
							.tab().append(query.member).tab().append(query.model).tab().append(query.start).tab().append(query.stop).tab().append(query.outcome)
							.tab().append(query.reference).tab().append(bp).tab().append(rl).tab().append(ql).tab().append(bp/(double)rl, 8).tab().append(bp/(double)ql, 8)
							.tab().append(2L*bp>=Math.max(rl, ql) ? "RECIPROCAL_HALF_GEOMETRY" : "BELOW_RECIPROCAL_HALF").nl());
					}
				}
			}
			for(final Genome genome:genomes){
				final int raw=presence(genome.actual[0][0], genome.actual[1][0]);
				final int scored=presence(genome.actual[0][1], genome.actual[1][1]);
				presence[0][raw]++;
				presence[1][scored]++;
				inputs.writer.print(buffer.clear().append(genome.id).tab().append(genome.actual[0][0]).tab().append(genome.actual[1][0])
					.tab().append(genome.actual[0][1]).tab().append(genome.actual[1][1]).tab().append(PRESENCE[raw]).tab().append(PRESENCE[scored]).nl());
				for(int side=0; side<2; side++){
					for(final Query query:genome.queries(side)){
						int overlapPairs=0;
						for(final Reference ref:genome.references){
							if(query.outcome.equals("TP") && ref.winners[side]!=query){continue;}
							final int bp=intersection(ref, query);
							if(bp<1){continue;}
							overlapPairs++;
							length(lengths, side, genome, query, ref, bp);
						}
						callLengths.writer.print(buffer.clear().append(side).tab().append(genome.id).tab().append(query.seqid)
							.tab().append(query.start).tab().append(query.stop).tab().append(query.strand)
							.tab().append(query.member).tab().append(query.model).tab().append(query.outcome)
							.tab().append(query.stop-query.start+1).tab().append(overlapPairs)
							.append("\tEVERY_CALL_ONCE_REFERENCE_PAIRS_ARE_NOT_EXTRA_TPS\n"));
					}
				}
			}
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println(SUMMARY_VERSION);
			summary.writer.println(SUMMARY_HEADER);
			count(summary, "graded_scored_references", "RETAINED", retained, "UNCHANGED_ONE_TO_ONE_MATCHES");
			count(summary, "graded_scored_references", "LOST", lost, "UNCHANGED_ONE_TO_ONE_MATCHES");
			count(summary, "graded_scored_references", "GAINED", gained, "UNCHANGED_ONE_TO_ONE_MATCHES");
			count(summary, "lost_references_without_any_candidate_overlap", "LOST", absent, "ADD_BACK_ELIGIBLE_PRESENCE_CLASS_ONLY");
			count(summary, "retained_references_with_changed_emitted_member", "RETAINED", switches, "LABEL_SWITCH_NOT_SUBTYPE_TRUTH");
			for(int i=0; i<2; i++){
				for(int state=0; state<PRESENCE.length; state++){
					count(summary, i==0 ? "genomes_with_any_clan_call" : "genomes_with_scored_reference_detection",
						PRESENCE[state], presence[i][state], "GENOME_PRESENCE_NOT_LOCUS_RECALL");
				}
			}
			for(int before=0; before<3; before++){
				for(int after=0; after<3; after++){
					count(summary, MEMBERS[before]+"_to_"+MEMBERS[after], "RETAINED", labels[before][after],
						"RETAINED_REFERENCE_EMITTED_LABELS_NOT_SUBTYPE_TRUTH");
				}
			}
		}
		if(pairedStats!=null){writePairedSummary();}
	}
	private void count(final Output output, final String measure, final String state, final long n, final String qualification){
		require(n>=0, "A comparison population count cannot be negative: "+measure);
		output.writer.print(buffer.clear().append(measure).tab().append(state).tab().append(n).tab().append(qualification).nl());
	}

	/** Uses only the same reference's two validated winners. Member cohorts stay
	 * fixed by OFF's label even when ON switches labels; lost/gained references
	 * remain in the existing separate locus ledger, never in these denominators.
	 */
	private void paired(final Reference ref, final Query before, final Query after, final Output output){
		final long a5=ref.strand=='-' ? (long)ref.stop-before.stop : (long)before.start-ref.start;
		final long z5=ref.strand=='-' ? (long)ref.stop-after.stop : (long)after.start-ref.start;
		final long a3=ref.strand=='-' ? (long)ref.start-before.start : (long)before.stop-ref.stop;
		final long z3=ref.strand=='-' ? (long)ref.start-after.start : (long)after.stop-ref.stop;
		final int cohort=member(before.member);
		final boolean changedLabel=!before.member.equals(after.member);
		if(replay!=null){replay.pair(after.trace, ref.id, cohort, ref.start, ref.stop, a5, a3);}
		for(int end=0; end<2; end++){
			final long a=end==0 ? a5 : a3, z=end==0 ? z5 : z3;
			pairedStats[cohort][end].add(a, z, radius, changedLabel);
			pairedStats[3][end].add(a, z, radius, changedLabel);
		}
		output.writer.print(buffer.clear().append(ref.genome.id).tab().append(ref.id).tab().append(ref.seqid).tab().append(ref.strand)
			.tab().append(before.member).tab().append(after.member).tab().append(a5).tab().append(z5).tab().append(a3).tab().append(z3)
			.tab().append(Math.abs(a5)<=radius ? "true" : "false").tab().append(Math.abs(a3)<=radius ? "true" : "false").nl());
	}

	/** Full paired RMS is primary; the OFF-selected small-error subset is secondary. */
	private void writePairedSummary(){
		try(Output output=new Output(out+".paired_endpoints.tsv"); Output near=new Output(out+".paired_near_endpoints.tsv");
			Output population=new Output(out+".paired_population.tsv"); Output details=new Output(out+".paired_details.tsv")){
			output.writer.println(PAIRED_VERSION);
			output.writer.println(PAIRED_HEADER);
			near.writer.println(NEAR_VERSION);
			near.writer.println(NEAR_HEADER);
			population.writer.println(POPULATION_VERSION);
			population.writer.println(POPULATION_HEADER);
			details.writer.println(DETAILS_VERSION);
			details.writer.println(DETAILS_HEADER);
			for(int cohort=0; cohort<4; cohort++){
				for(int end=0; end<2; end++){
					final PairStats stats=pairedStats[cohort][end];
					require(stats.n==stats.better+stats.worse+stats.unchanged && stats.worse==stats.worseWithin+stats.worseBeyond
						&& stats.withinN>=0 && stats.withinN<=stats.n,
						"Every paired end must have exactly one error-change disposition and every worsening one OFF-radius class");
					pairKey(cohort, end, stats.n);
					population.writer.print(buffer.tab().append(stats.withinN).nl());
					if(stats.n==0){continue;}// No distribution exists; its zero population remains explicit above.
					writeDetails(details, cohort, end, stats);
					pairKey(cohort, end, stats.n).tab();
					rms(buffer, stats.offSquares, stats.n);
					rms(buffer.tab(), stats.onSquares, stats.n);
					buffer.tab().append(stats.better).tab().append(stats.worse).tab().append(stats.unchanged)
						.tab().append(stats.worseWithin).tab().append(stats.worseBeyond);
					output.writer.print(buffer.tab().append(stats.switches)
						.append("\tSAME_SCORED_REFERENCE_IDS_OFF_EMITTED_MEMBER_COHORT_FULL_RMS_PRIMARY\n"));
					if(stats.withinN>0){
						pairKey(cohort, end, stats.n).tab().append(stats.withinN).tab();
						rms(buffer, stats.offWithinSquares, stats.withinN);
						rms(buffer.tab(), stats.onWithinSquares, stats.withinN);
						near.writer.print(buffer.tab().append(stats.worseWithin).append("\tSAME_OFF_SELECTED_SUBSET\n"));
					}
				}
			}
		}
	}
	/** All measures use this same recorded reference intersection. A sign flip
	 * moves an end even when its absolute distance is unchanged.
	 */
	private void writeDetails(final Output output, final int cohort, final int end, final PairStats stats){
		require(stats.n>0 && stats.offExact<=stats.offWithin3 && stats.offWithin3<=stats.n
			&& stats.onExact<=stats.onWithin3 && stats.onWithin3<=stats.n
			&& stats.moved>=stats.better+stats.worse && stats.moved<=stats.n
			&& Double.isFinite(stats.offAbsolute) && Double.isFinite(stats.onAbsolute),
			"Paired detail counts share N; every distance change requires movement, but equal-distance sign flips also move");
		pairKey(cohort, end, stats.n).tab().append(stats.offExact).tab().append(stats.onExact)
			.tab().append(stats.offWithin3).tab().append(stats.onWithin3).tab();
		rms(buffer, stats.offSquares, stats.n);
		rms(buffer.tab(), stats.onSquares, stats.n);
		output.writer.print(buffer.tab().append(stats.offAbsolute/stats.n, 8).tab().append(stats.onAbsolute/stats.n, 8)
			.tab().append(stats.moved).tab().append(stats.better).tab().append(stats.worse).tab().append(stats.unchanged)
			.tab().append(stats.switches).append("\tSAME_SCORED_REFERENCE_IDS_OFF_MEMBER_COHORT_NO_REMATCHING\n"));
	}
	private ByteBuilder pairKey(final int cohort, final int end, final long n){
		assert(cohort>=0 && cohort<4 && end>=0 && end<2) : "Paired statistics use the three OFF-member cohorts plus their clan and two ends";
		return buffer.clear().append(cohort<3 ? MEMBERS[cohort] : "cobalamin").tab().append(end==0 ? "5prime" : "3prime")
			.tab().append(radius).tab().append(n);
	}
	private static void rms(final ByteBuilder buffer, final double squares, final long count){
		require(count>0 && Double.isFinite(squares) && squares>=0,
			"RMS requires an existing population and finite nonnegative squared errors; omit empty rows, never hide arithmetic failures");
		buffer.append(Math.sqrt(squares/count), 8);
	}

	private void length(final Output out, final int side, final Genome genome, final Query query, final Reference ref, final int bp){
		require(ref!=null && bp>0, "Reference length statistics exist only for an actual overlapping reference");
		final int ql=query.stop-query.start+1, rl=ref.stop-ref.start+1;
		buffer.clear().append(side).tab().append(genome.id).tab().append(query.member).tab().append(query.model).tab().append(query.outcome)
			.tab().append(ref.id).tab().append(ref.scored ? "SCORED" : "MASK")
			.tab().append(ql).tab().append(rl).tab().append(bp).tab();
		buffer.append(ql/(double)rl, 8);
		out.writer.print(buffer.append("\tMATCHED_UNIQUE_OR_UNMATCHED_OVERLAP_PAIRS_NOT_EXTRA_TPS\n"));
	}
	private static void winner(final ByteBuilder buffer, final Query query){
		if(query==null){buffer.append("NA\tNA");}else{buffer.append(query.member).tab().append(query.model);}
	}
	private static int intersection(final Reference ref, final Query query){
		if(ref.strand!=query.strand || !ref.seqid.equals(query.seqid)){return 0;}
		return Math.max(0, Math.min(ref.stop, query.stop)-Math.max(ref.start, query.start)+1);
	}
	private static int presence(final long before, final long after){return before>0 ? after>0 ? 0 : 1 : after>0 ? 2 : 3;}
	private static int member(final String value){
		for(int i=0; i<MEMBERS.length; i++){if(value.equals(MEMBERS[i])){return i;}}
		throw new IllegalArgumentException("Unknown emitted member in a normalized grade: "+value);
	}
	/** Native FASTA IDs exclude descriptions. Grade fields escape literal percent
	 * once; native library headers reject control characters, so no control escape
	 * can belong to a valid source model identifier. Full report text stays intact.
	 */
	static String modelId(final String model){return NcrnaSupport.recordName(model.replace("%25", "%"));}
	private Input input(final int side, final String suffix, final String header){
		final String path=prefixes[side]+suffix;
		final Input in=new Input(path);
		try{
			RrnaResourceIO.header(in.reader, header, path);
			return in;
		}
		catch(RuntimeException | Error failure){
			try{in.close();}catch(RuntimeException | Error cleanup){failure.addSuppressed(cleanup);}
			throw failure;
		}
	}
	private Genome genome(){
		probe.set(p, 0);
		final Genome genome=byInput.get(probe);
		require(genome!=null, "Detailed normalized row names an undeclared input");
		return genome;
	}
	private void referenceKey(final Genome genome, final int column){
		key.clear().append(genome.id).tab();
		p.setBounds(column);
		key.append(p.line(), p.a(), p.b()-p.a());
		probe.set(key.array, 0, key.length());
	}
	private char strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Explicit strand required");
		return p.termEquals("-", column) ? '-' : '+';
	}
	private boolean flag(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Explicit boolean required");
		return p.termEquals("true", column);
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------             State            ----------------*/
	/*--------------------------------------------------------------*/

	static final class Genome{
		Genome(final String id_){id=id_;}
		ArrayList<Query> queries(final int side){return side==0 ? before : after;}
		final String id;
		final long[][] expected=new long[2][4], actual=new long[2][4];
		final ArrayList<Query> before=new ArrayList<Query>(), after=new ArrayList<Query>();
		final ArrayList<Reference> references=new ArrayList<Reference>();
		boolean seen;
	}
	static final class Reference{
		Reference(final Genome genome_, final String id_, final String seqid_, final int start_, final int stop_, final char strand_, final boolean scored_){
			genome=genome_;
			id=id_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			scored=scored_;
		}
		final Genome genome;
		final String id, seqid;
		final int start, stop;
		final char strand;
		final boolean scored;
		final Query[] winners=new Query[2];
		boolean seen;
	}
	static final class Query{
		Query(final String seqid_, final int start_, final int stop_, final char strand_, final String member_, final String model_, final String outcome_, final String reference_){
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			member=member_;
			model=model_;
			outcome=outcome_;
			reference=reference_;
		}
		final String seqid, member, model, outcome, reference;
		final int start, stop;
		final char strand;
		CobalaminBoundaryReplay.Call trace;
	}
	private String out;
	private boolean changesOnly;
	private String trace, bindings;
	private CobalaminBoundaryReplay replay;
	private float[] replayMargins;
	private int radius=-1;
	private PairStats[][] pairedStats;
	private String[] subset;
	private long retained, lost, gained, absent, switches;
	private final long[][] presence=new long[2][4];
	private final long[][] labels=new long[3][3];
	private final String[] prefixes=new String[2];
	private final ObjectMap<Target, Genome> byInput=new ObjectMap<Target, Genome>(16, Target.class, Genome.class);
	private final ObjectMap<Target, Reference> byReference=new ObjectMap<Target, Reference>(16, Target.class, Reference.class);
	private final ObjectSet<String> kept=new ObjectSet<String>(String.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Reference> references=new ArrayList<Reference>();
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder buffer=new ByteBuilder(), key=new ByteBuilder();
	private static final String[] PRESENCE={"RETAINED", "LOST", "GAINED", "ABSENT_BOTH"};
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String GENOMES="input_id\tcalls\tmatched\tmasked\tunmatched";
	private static final String REFERENCES="input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored";
	private static final String QUERIES="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id";
	private static final String WITNESSES="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
	public static final String SUMMARY_VERSION="#COBALAMIN_COMPARISON_COUNTS_V2";
	public static final String SUMMARY_HEADER="measure\tstate\tcount\tqualification";
	public static final String LENGTH_VERSION="#COBALAMIN_REFERENCE_LENGTHS_V2";
	public static final String CALL_LENGTH_VERSION="#COBALAMIN_CALL_LENGTHS_V2";
	public static final String PAIRED_VERSION="#COBALAMIN_PAIRED_ENDPOINTS_V2";
	public static final String PAIRED_LOCI_HEADER="input_id\treference_id\tseqid\tstrand\toff_member\ton_member\toff_d5\ton_d5\toff_d3\ton_d3\toff5_within_radius\toff3_within_radius";
	public static final String PAIRED_HEADER="cohort\tend\tradius_bp\tpaired_n\toff_rms_nt\ton_rms_nt\tbetter\tworse\tunchanged\tworse_off_within_radius\tworse_off_beyond_radius\tlabel_switches\tqualification";
	public static final String NEAR_VERSION="#COBALAMIN_PAIRED_NEAR_ENDPOINTS_V2";
	public static final String NEAR_HEADER="cohort\tend\tradius_bp\tpaired_n\toff_within_radius_n\toff_within_radius_rms_nt\ton_same_subset_rms_nt\tworse_off_within_radius\tqualification";
	public static final String POPULATION_VERSION="#COBALAMIN_PAIRED_POPULATION_V2";
	public static final String POPULATION_HEADER="cohort\tend\tradius_bp\tpaired_n\toff_within_radius_n";
	public static final String DETAILS_VERSION="#COBALAMIN_PAIRED_DETAILS_V1";
	public static final String DETAILS_HEADER="cohort\tend\tradius_bp\tpaired_n\toff_exact\ton_exact\toff_within3\ton_within3"
		+"\toff_rms_nt\ton_rms_nt\toff_mean_absolute_nt\ton_mean_absolute_nt\tmoved\tbetter\tworse\tunchanged\tlabel_switches\tqualification";

	/** Primitive accumulation on one fixed paired cohort/end; no transient rows. */
	private static final class PairStats{
		void add(final long off, final long on, final int radius, final boolean switched){
			final double a=Math.abs((double)off), z=Math.abs((double)on);
			n++;
			offSquares+=a*a;
			onSquares+=z*z;
			offAbsolute+=a;
			onAbsolute+=z;
			if(off==0){offExact++;}
			if(on==0){onExact++;}
			if(a<=3){offWithin3++;}
			if(z<=3){onWithin3++;}
			if(off!=on){moved++;}
			if(z<a){better++;}
			else if(z>a){
				worse++;
				if(a<=radius){worseWithin++;}else{worseBeyond++;}
			}else{unchanged++;}
			if(a<=radius){
				withinN++;
				offWithinSquares+=a*a;
				onWithinSquares+=z*z;
			}
			if(switched){switches++;}
			assert(Double.isFinite(offSquares) && Double.isFinite(onSquares))
				: "Paired coordinate differences are widened before squaring; finite RMS must not be corrupted by arithmetic overflow";
		}
		long n, better, worse, unchanged, worseWithin, worseBeyond, withinN, switches;
		long offExact, onExact, offWithin3, onWithin3, moved;
		double offSquares, onSquares, offWithinSquares, onWithinSquares;
		double offAbsolute, onAbsolute;
	}
}
