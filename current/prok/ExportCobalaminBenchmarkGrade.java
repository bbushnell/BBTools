package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
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

/** Exports already assigned cobalamin union matches for benchmark reporting.
 * Adapted from template.A_SampleByteFile, using native byte streams and reused
 * parsers/builders. Accepts completed canonical grades or the historical
 * CobalaminPanelGrade FULL_CM/union view. No calling, rematching, truth editing
 * or biological subtype inference occurs. Emitted members partition union calls.
 * Explicit inexact mode exports the grader's recorded overlap matches without
 * assigning them again; default canonical and historical exports stay strict.
 * Optional mode=members summarizes disjoint canonical overlap grades with
 * separate emitted-query and representative-reference member axes.
 *
 * Retains only selected query/reference rows and one counter record per input,
 * including zero-call inputs. Coordinates come from actual reference rows;
 * recorded signed errors must agree. The nine-column matches contain both the
 * clan and emitted-member views for NcrnaEndpointSummary; do not sum those views.
 * Each instance, its parsers and stream owners are confined to one thread.
 * @author Raiden
 * @author Diona
 */
public final class ExportCobalaminBenchmarkGrade{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ExportCobalaminBenchmarkGrade.class, false);
		try{new ExportCobalaminBenchmarkGrade(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	private ExportCobalaminBenchmarkGrade(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("legacy")){legacy=Parse.parseBoolean(b);}
			else if(a.equals("inexact")){inexact=Parse.parseBoolean(b);}
			else if(a.equals("mode")){
				require("export".equals(b) || "members".equals(b), "mode must be export or members");
				memberMode="members".equals(b);
			}
			//TODO: Probable bug - canonical export accepts panel= but never applies
			// it; only readLegacyReferences consumes this field. Reject or implement
			// canonical subsets before allowing this option to select a population.
			else if(a.equals("panel")){panel=b;}
			else if(a.equals("timing")){timing=b;}
			else if(a.equals("callerlog")){callerLog=b;}
			else if(a.equals("ref")){modelFiles=b.split(",", -1);}
			else if(a.equals("expectedmodels")){modelCounts=Parse.parseIntArray(b, ",");}
			else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown benchmark-export parameter: "+arg);}
		}
		prefix=parser.in1;
		out=parser.out1;
		require(prefix!=null && out!=null && (!legacy || panel!=null), "Supply in=completed_grade_prefix out=fresh_prefix; legacy=t also requires panel=original_panel.loci.tsv");
		require(!legacy || !inexact, "inexact=t requires a canonical overlap grade; legacy grades retain their recorded matching");
		require((timing==null)==(callerLog==null), "Optional timing=GNU_time_triple and callerlog=complete_caller_stderr must be supplied together");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Export complete grades to fresh outputs without sampling or append");
		Shared.setThreads(1);
		if(memberMode){
			require(inexact && !legacy && panel==null && timing==null && callerLog==null
				&& modelFiles==null && modelCounts==null && expectedInputs<0,
				"mode=members requires canonical inexact=t grades and no legacy, timing or model-roster options");
			memberPrefixes=prefix.split(",", -1);
			for(final String part:memberPrefixes){
				require(!part.isEmpty() && Tools.testInputFiles(false, true, part+".inputs.tsv",
					part+".overlap.raw.tsv", part+".overlap.summary.tsv"), "Member-summary partitions require complete overlap grades");
			}
			RrnaResourceIO.fresh(out, ".member_outcomes.tsv");
			return;
		}
		require(Tools.testInputFiles(false, true, prefix+(inexact ? ".overlap.summary.tsv" : ".summary.tsv"),
			prefix+(legacy ? ".calls.tsv" : inexact ? ".overlap.raw.tsv" : ".raw.tsv"),
			prefix+(legacy ? ".genomes.tsv" : ".inputs.tsv")), "Completed grade reports must be readable");
		require(Tools.testOutputFiles(false, false, false, out+".members.tsv", out+".genomes.tsv", out+".matches.tsv",
			out+".witnesses.tsv", out+".queries.tsv", out+".references.tsv", out+".metrics.tsv"), "Benchmark exports must be fresh");
		if(modelFiles!=null || modelCounts!=null || expectedInputs>=0){
			require(modelFiles!=null && modelFiles.length==3 && modelCounts!=null && modelCounts.length==3
				&& expectedInputs>0 && !legacy, "Model-win reporting requires three ordered ref=FASTAs, expectedmodels=MAIN,SECOND,THIRD and expectedinputs=N on a canonical grade");
			for(final int count:modelCounts){require(count>0, "Every model-win roster must have its positive complete model count");}
			require(Tools.testInputFiles(false, true, modelFiles), "Complete model-win rosters must be readable and distinct");
			RrnaResourceIO.fresh(out, ".model_wins.tsv", ".model_win_summary.tsv");
		}
	}

	/** One partition at a time reuses the ordinary reader and summary checks.
	 * The argument constructor validates all inputs and the fresh output path. */
	private ExportCobalaminBenchmarkGrade(final String partition){
		assert(partition!=null && !partition.isEmpty()) : "A member summary must retain its original grade prefix";
		prefix=partition;
		inexact=true;
		memberMode=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reading           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		if(memberMode){
			summarizeMembers();
			return;
		}
		readInputs();
		if(legacy){readLegacyReferences();}
		readRows();
		checkSummary();
		if(timing!=null){readMeasurements();}
		if(modelFiles!=null){countModelWins();}
		write();
		if(modelFiles!=null){writeModelWins();}
	}

	/** Aggregates disjoint saved partitions, never a new biological match. Local
	 * feature line numbers may repeat; original input IDs may not.
	 * Counters are [axis][member][slot]: QUERY axis0 slots0..3 are TP,
	 * TP_INEXACT, MASKED, FP; REFERENCE axis1 slots0..3 are TP, TP_INEXACT,
	 * FN, DONT_CARE. The same layout is used by countMemberOutcomes. */
	private void summarizeMembers(){
		assert(memberMode && memberPrefixes!=null && memberPrefixes.length>0) : "Member aggregation needs explicitly declared grade partitions";
		final ObjectSet<Target> inputs=new ObjectSet<Target>(Target.class);
		final long[][][] counts=new long[2][3][4];
		for(final String path:memberPrefixes){
			final ExportCobalaminBenchmarkGrade part=new ExportCobalaminBenchmarkGrade(path);
			part.readInputs();
			for(final Genome genome:part.genomes){
				final byte[] id=genome.id.getBytes(StandardCharsets.UTF_8);
				probe.set(id, 0, id.length);
				if(!inputs.add(probe.copy())){
					throw new IllegalArgumentException("Repeated input ID across member-summary partitions: "+genome.id+" in "+path);
				}
			}
			part.readRows();
			part.checkSummary();
			part.countMemberOutcomes(counts);
		}
		try(Output rows=new Output(out+".member_outcomes.tsv")){
			rows.writer.println("axis\tmember\tinputs\ttotal\ttp\ttp_inexact\tmatched\tunmatched\tmasked\tmissed\tdont_care\tqualification");
			for(int axis=0; axis<2; axis++){
				for(int member=0; member<3; member++){
					final long[] c=counts[axis][member];
					buffer.clear().append(axis==0 ? "QUERY_EMITTED_MEMBER" : "REFERENCE_REPRESENTATIVE_MEMBER")
						.tab().append(MEMBERS[member]).tab().append(inputs.size()).tab().append(c[0]+c[1]+c[2]+c[3])
						.tab().append(c[0]).tab().append(c[1]).tab().append(c[0]+c[1]);
					if(axis==0){buffer.tab().append(c[3]).tab().append(c[2]).append("\tNA\tNA");}
					else{buffer.append("\tNA\tNA\t").append(c[2]).tab().append(c[3]);}
					rows.writer.print(buffer.append("\tRECORDED_UNION_ASSIGNMENTS_NOT_SUBTYPE_ACCURACY\n"));
				}
			}
		}
	}

	/** The two axes partition the same union assignments. REFERENCE model is
	 * its recorded representative member, never inferred from the emitted model. */
	private void countMemberOutcomes(final long[][][] totals){
		assert(memberMode && inexact && !legacy) : "Member axes must use recorded canonical overlap assignments";
		final long[][][] counts=new long[2][3][4];
		for(final Query query:queries){
			p.set(query.row);
			final int outcome=query.bucket==0 ? (p.termEquals("TP_INEXACT", 7) ? 1 : 0) : query.bucket==1 ? 2 : 3;
			counts[0][query.member][outcome]++;
			final long line=p.parseLong(9);
			key.clear().append(query.genome.id).tab().append(line);
			probe.set(key.array, 0, key.length());
			require(line>0 && seen.add(probe.copy()), "Query line IDs must be unique within an input; repeated feature rows cannot count twice");
			if(query.bucket!=0){continue;}
			referenceKey(query.genome, 10);
			final Reference ref=byReference.get(probe);
			require(ref!=null && ref.scored && !ref.consumed, "Every member-summary TP must consume exactly one declared scored reference");
			require(ref.outcome==outcome && ref.queryLine==line,
				"Recorded QUERY/REFERENCE match class and reciprocal query line must agree; no rematching is allowed");
			final int start=p.parseInt(4), stop=p.parseInt(5);
			require(p.termEquals(ref.seqid, 3) && strand(6)==ref.strand && start>0 && stop>=start
				&& p.parseInt(16)==(ref.strand=='-' ? ref.stop-stop : start-ref.start)
				&& p.parseInt(17)==(ref.strand=='-' ? ref.start-start : stop-ref.stop),
				"Recorded member-summary witness must retain its contig, strand and signed endpoint offsets");
			ref.consumed=true;
		}
		long matched=0, inexactCount=0, missed=0, masked=0;
		for(final Reference ref:references){
			require(ref.consumed==(ref.outcome<2), "Every recorded matched reference must have its reciprocal QUERY witness; FN and masks cannot be consumed");
			counts[1][ref.member][ref.outcome]++;
			if(ref.outcome<2){matched++;}
			if(ref.outcome==1){inexactCount++;}
			if(ref.outcome==2){missed++;}
			if(ref.outcome==3){masked++;}
		}
		require(matched+missed==scoredReferences && inexactCount==inexactMatches
			&& missed==missedReferences && masked==maskedReferences,
			"Reference member counts must conserve the recorded union matched/inexact/missed/masked totals");
		for(int axis=0; axis<2; axis++){
			for(int member=0; member<3; member++){
				for(int outcome=0; outcome<4; outcome++){totals[axis][member][outcome]+=counts[axis][member][outcome];}
			}
		}
	}

	/** Bind all current models, including zero-call/zero-win models. Counts refer
	 * to already assigned scored QUERY matches, not construction membership. */
	private void countModelWins(){
		require(genomes.size()==expectedInputs, "Model-win input population differs from the declared complete panel");
		final ObjectMap<Target, ModelWins> byId=new ObjectMap<Target, ModelWins>(Target.class, ModelWins.class);
		for(int member=0; member<3; member++){
			final boolean trimmed=Shared.TRIM_READ_DESCRIPTION;
			final NcrnaSupport.Library library;
			try{
				Shared.TRIM_READ_DESCRIPTION=false;
				library=NcrnaSupport.loadLibrary(modelFiles[member]);
			}finally{Shared.TRIM_READ_DESCRIPTION=trimmed;}
			require(library.names.length==modelCounts[member], "Model-win roster count differs for "+MEMBERS[member]);
			for(int i=0; i<library.names.length; i++){
				final String id=library.names[i];
				final ModelWins model=new ModelWins(member, i, id, library.headers[i]);
				final byte[] key=id.getBytes(StandardCharsets.UTF_8);
				probe.set(key, 0, key.length);
				require(id.startsWith(PREFIXES[member]) && byId.put(probe.copy(), model)==null, "Model-win IDs must be unique and belong to their declared member");
				modelWins.add(model);
			}
		}
		final long[] totals=new long[3];
		for(final Query query:queries){
			p.set(query.row);
			p.setBounds(15);
			int end=p.a();
			while(end<p.b() && (p.line()[end]&255)>32){end++;}
			probe.set(p.line(), p.a(), end);
			final ModelWins model=byId.get(probe);
			if(model==null || model.member!=query.member){
				throw new IllegalArgumentException("Every emitted model must belong to the supplied current roster: "+p.parseString(15));
			}
			if(query.bucket==0){
				model.wins++;
				totals[model.member]++;
			}
		}
		for(int member=0; member<3; member++){
			require(totals[member]==members[member][0], "Model-win counts must conserve each emitted member's assigned matches");
		}
	}

	private void writeModelWins(){
		try(Output rows=new Output(out+".model_wins.tsv"); Output summary=new Output(out+".model_win_summary.tsv")){
			rows.writer.println("emitted_member\tmember_model_ordinal0\tmodel_id\tmatched_wins\tmodel_header");
			for(final ModelWins model:modelWins){
				rows.writer.print(buffer.clear().append(MEMBERS[model.member]).tab().append(model.ordinal).tab().append(model.id)
					.tab().append(model.wins).tab().append(model.header).nl());
			}
			summary.writer.println("emitted_member\tinputs\tmodels\tmatched_wins\tmin\tq1\tmedian\tq3\tmax\tzero_wins\tfewer_than_15\tfewer_than_46\tquantile_rule\tqualification");
			for(int member=-1; member<3; member++){
				final long[] values=new long[member<0 ? modelWins.size() : modelCounts[member]];
				int next=0, zero=0, under15=0, under46=0;
				long sum=0;
				for(final ModelWins model:modelWins){
					if(member>=0 && model.member!=member){continue;}
					values[next++]=model.wins;
					sum+=model.wins;
					if(model.wins==0){zero++;}
					if(model.wins<15){under15++;}
					if(model.wins<46){under46++;}
				}
				assert(next==values.length) : "Every declared model, including zeros, belongs to exactly one member distribution";
				Arrays.sort(values);
				buffer.clear().append(member<0 ? "ALL" : MEMBERS[member]).tab().append(genomes.size()).tab().append(values.length)
					.tab().append(sum).tab().append(values[0]).tab().append(quantile(values, .25), 6, true)
					.tab().append(quantile(values, .5), 6, true).tab().append(quantile(values, .75), 6, true)
					.tab().append(values[values.length-1]).tab().append(zero).tab().append(under15).tab().append(under46)
					.append("\tLINEAR_INTERPOLATION_AT_(N-1)*P\t")
					.append(inexact ? "INCLUSIVE_ASSIGNED_QUERY_MATCHES" : "STRICT_ASSIGNED_QUERY_MATCHES")
					.append("_NOT_CONSTRUCTION_SUPPORT_OR_UNIQUE_MODEL_COVERAGE\n");
				summary.writer.print(buffer);
			}
		}
	}

	private static double quantile(final long[] sorted, final double fraction){
		assert(sorted.length>0 && fraction>=0 && fraction<=1) : "Quantiles include all declared models in sorted nonempty populations";
		final double position=(sorted.length-1)*fraction;
		final int lower=(int)position, upper=Math.min(lower+1, sorted.length-1);
		return sorted[lower]+(position-lower)*(sorted[upper]-sorted[lower]);
	}

	/** Copies final measurement fields; never estimates CPU or repeats calling. */
	private void readMeasurements(){
		try(Input in=new Input(timing)){
			final byte[] row=in.reader.nextLine();
			require(row!=null && in.reader.nextLine()==null, "Native timing requires exactly one user/system/wall triple");
			p.set(row);
			require(p.terms()==3, "Native timing requires three tab-separated seconds fields");
			for(int i=0; i<3; i++){
				final double value=p.parseDouble(i);
				require(Double.isFinite(value) && value>=0, "Native timing values must be finite nonnegative seconds");
				times[i]=p.parseString(i);
			}
		}
		try(Input in=new Input(callerLog)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				for(int i=0; i<3; i++){
					if(!p.termEquals(ALIGNMENT_LABELS[i], 0)){continue;}
					require(p.terms()==2 && alignments[i]<0, "Each member needs exactly one native alignment counter");
					alignments[i]=Long.parseLong(p.parseString(1).trim());
					require(alignments[i]>=0, "Native alignment counts cannot be negative");
				}
			}
		}
		for(final long value:alignments){require(value>=0, "Complete caller log lacks a cobalamin member alignment count");}
	}

	/** Keeps every declared input, including genomes with no clan calls. */
	private void readInputs(){
		final String path=prefix+(legacy ? ".genomes.tsv" : ".inputs.tsv");
		try(Input in=new Input(path)){
			if(legacy){RrnaResourceIO.header(in.reader, LEGACY_GENOMES, path);}
			else{RrnaResourceIO.header(in.reader, GradeNcrnaReference.INPUT_HEADER, path);}
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==(legacy ? 13 : 10), "Input report schema changed");
				if(legacy && !legacyUnion()){continue;}
				final int column=legacy ? 2 : 0;
				probe.set(p, column);
				require(byInput.get(probe)==null, "A completed grade must declare each input exactly once");
				final Genome genome=new Genome(p.parseString(column));
				byInput.put(probe.copy(), genome);
				genomes.add(genome);
			}
		}
		require(!genomes.isEmpty(), "An export requires the complete nonempty input population");
	}

	/** Legacy full-core coordinates stay separate from source annotation intervals. */
	private void readLegacyReferences(){
		try(Input in=new Input(panel)){
			RrnaResourceIO.header(in.reader, CobalaminPanelPrep.HEADER, panel);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==10, "Legacy panel requires ten fields");
				final boolean scored=p.parseInt(7)>0;
				if(!scored){require(p.parseInt(7)==0 && p.parseInt(8)==0, "Legacy unresolved masks must not claim full-core endpoints");}
				addReference(1, 0, 2, scored ? 7 : 5, scored ? 8 : 6, 9, scored);
			}
		}
	}

	/** Reads the recorded union outcomes; excluded canonical families are skipped. */
	private void readRows(){
		final String path=prefix+(legacy ? ".calls.tsv" : inexact ? ".overlap.raw.tsv" : ".raw.tsv");
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, legacy ? LEGACY_CALLS : GradeNcrnaReference.RAW_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==(legacy ? 14 : 18), "Completed detailed grade schema changed");
				if(legacy){
					if(!legacyUnion()){continue;}
				}else{
					if(!p.termEquals("cobalamin", 2)){continue;}
					if(p.termEquals("REFERENCE", 1)){
						final Reference ref=addReference(0, 8, 3, 4, 5, 6, !p.termEquals("DONT_CARE", 7));
						if(memberMode){
							for(int i=0; i<3; i++){
								if(p.termEquals(MEMBERS[i], 15)){
									ref.member=i;
									break;
								}
							}
							require(ref.member>=0, "Raw reference model must name its representative RF00174/RF01482/RF01689 member");
							ref.outcome=p.termEquals("TP", 7) ? 0 : p.termEquals("TP_INEXACT", 7) ? 1
								: p.termEquals("FN", 7) ? 2 : p.termEquals("DONT_CARE", 7) ? 3 : -1;
							require(ref.outcome>=0, "Raw reference has an unknown recorded member-summary outcome");
							if(ref.outcome<2){
								ref.queryLine=p.parseLong(11);
								require(ref.queryLine>0, "Matched reference must retain its original query line");
							}else{
								require(p.termEquals(".", 11), "Missed and masked references cannot name a scored query witness");
							}
						}
						continue;
					}
					require(p.termEquals("QUERY", 1), "Selected canonical row must be QUERY or REFERENCE");
				}
				addQuery(row);
			}
		}
	}

	/** One retained reference owns its genome-qualified identity and coordinates. */
	private Reference addReference(final int input, final int id, final int seqid, final int start,
			final int stop, final int strand, final boolean scored){
		final Genome genome=genome(input);
		referenceKey(genome, id);
		require(byReference.get(probe)==null, "Reference IDs must be unique within each genome; legacy ambiguity cannot be guessed");
		final Reference ref=new Reference(genome, p.parseString(id), p.parseString(seqid), p.parseInt(start), p.parseInt(stop), strand(strand), scored);
		require(ref.start>0 && ref.stop>=ref.start, "Reference endpoints must be ordered one-based inclusive coordinates");
		byReference.put(probe.copy(), ref);
		references.add(ref);
		if(scored){scoredReferences++;}else{maskedReferences++;}
		return ref;
	}

	/** Retains one selected query row; no per-field transient Strings are created. */
	private void addQuery(final byte[] row){
		final Genome genome=genome(legacy ? 3 : 0);
		final int member=member(legacy ? 9 : 15);
		final int outcome=legacy ? 10 : 7;
		final boolean endpointInexact=inexact && p.termEquals("TP_INEXACT", outcome);
		final boolean tp=p.termEquals("TP", outcome) || endpointInexact;
		final boolean masked=p.termEquals(legacy ? "UNSCORED" : "MASKED", outcome);
		require(tp || masked || p.termEquals(legacy ? "APPARENT_FP" : "FP", outcome), "Union query has an unknown recorded outcome");
		final boolean duplicate;
		if(legacy){
			key.clear().append(genome.id).tab();
			field(key, 4).tab();
			field(key, 8).tab().append(p.parseInt(6)).tab().append(p.parseInt(7));
			probe.set(key.array, 0, key.length());
			duplicate=!seen.add(probe.copy());
		}else{
			require(p.termEquals("true", 13) || p.termEquals("false", 13), "Duplicate flag must retain its explicit canonical value");
			duplicate=p.termEquals("true", 13);
		}
		require(!tp || !duplicate, "A duplicate query cannot consume a scored reference");
		final int bucket=tp ? 0 : masked ? 1 : duplicate ? 2 : 3;
		members[member][bucket]++;
		genome.counts[bucket]++;
		if(duplicate){duplicates++;}
		if(endpointInexact){inexactMatches++;}
		queries.add(new Query(genome, member, bucket, row));
	}

	/** Existing summaries reconcile the complete emitted-member partition. */
	private void checkSummary(){
		if(inexact){checkInexactSummary();return;}
		final String path=prefix+".summary.tsv";
		boolean found=false;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, legacy ? LEGACY_SUMMARY : GradeNcrnaReference.SUMMARY_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(legacy ? !legacyUnion() : !p.termEquals("cobalamin", 0)){continue;}
				require(!found && p.terms()==(legacy ? 16 : 17), "Exactly one joint union summary must exist");
				found=true;
				final long[] sum=new long[4];
				for(final long[] member:members){for(int i=0; i<4; i++){sum[i]+=member[i];}}
				expected5=p.parseLong(legacy ? 13 : 10);
				expected3=p.parseLong(legacy ? 14 : 11);
				expectedBoth=p.parseLong(legacy ? 15 : 12);
				missedReferences=p.parseLong(legacy ? 7 : 4);
				precision=p.parseString(legacy ? 10 : 13);
				recall=p.parseString(legacy ? 11 : 14);
				require(scoredReferences==p.parseLong(legacy ? 3 : 5) && scoredReferences==sum[0]+missedReferences,
					"Normalized scored references must reconcile native TP/FN counts");
				if(!legacy){require(maskedReferences==p.parseLong(6), "Canonical mask-reference count must agree with the native summary");}
				require(p.parseLong(legacy ? 2 : 1)==genomes.size()
					&& p.parseLong(legacy ? 5 : 2)==sum[0]
					&& p.parseLong(legacy ? 8 : 7)==sum[1]
					&& p.parseLong(legacy ? 6 : 3)==sum[2]+sum[3]
					&& p.parseLong(legacy ? 9 : 8)==duplicates, "Emitted-member counts disagree with the completed union grade");
			}
		}
		require(found, "The completed grade lacks its selected union summary");
	}

	/** Reconciles existing inclusive assignments; never recomputes overlap matches.
	 * The ordinary input roster retains zero-call genomes independently of the
	 * per-family inclusive count rows. Endpoint exactness still uses actual bounds. */
	private void checkInexactSummary(){
		final String path=prefix+".overlap.summary.tsv";
		boolean found=false;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, GradeNcrnaReference.INEXACT_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(!p.termEquals("cobalamin", 1)){continue;}
				require(!found && p.terms()==19 && p.termEquals("ALL", 0), "Exactly one complete inclusive union summary must exist");
				found=true;
				final long[] sum=new long[4];
				for(final long[] member:members){for(int i=0; i<4; i++){sum[i]+=member[i];}}
				final long matched=sum[0], unmatched=sum[2]+sum[3];
				missedReferences=p.parseLong(9);
				require(p.parseLong(8)==inexactMatches && inexactMatches<=matched
					&& p.parseLong(3)==matched-inexactMatches, "Recorded inexact match counts disagree with the inclusive raw assignments");
				require(missedReferences>=0 && p.parseLong(2)==genomes.size() && p.parseLong(7)==matched
					&& p.parseLong(10)==unmatched && p.parseLong(11)==sum[1]
					&& p.parseLong(12)==scoredReferences && scoredReferences==matched+missedReferences
					&& p.parseLong(13)==maskedReferences && p.parseLong(14)==duplicates,
					"Inclusive emitted-member, mask and reference counts must reconcile the complete recorded union");
				expected5=p.parseLong(15);
				expected3=p.parseLong(16);
				expectedBoth=p.parseLong(17);
				precision=matched+unmatched==0 ? "NA" : Double.toString(matched/(double)(matched+unmatched));
				recall=scoredReferences==0 ? "NA" : Double.toString(matched/(double)scoredReferences);
			}
		}
		require(found, "The completed overlap grade lacks its selected union summary");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes source-qualified views; summary counters are checked before completion. */
	private void write(){
		try(Output matches=new Output(out+".matches.tsv"); Output witness=new Output(out+".witnesses.tsv");
			Output inputs=new Output(out+".genomes.tsv"); Output counts=new Output(out+".members.tsv");
			Output allQueries=new Output(out+".queries.tsv"); Output allReferences=new Output(out+".references.tsv")){
			matches.writer.println("#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3");
			witness.writer.println("input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3");
			allQueries.writer.println("input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id");
			allReferences.writer.println("input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored");
			inputs.writer.println("input_id\tcalls\tmatched\tmasked\tunmatched");
			counts.writer.println("emitted_member\tcalls\tmatched\tmasked\tunmatched_duplicate\tunmatched_other\tqualification");
			long exact5=0, exact3=0, both=0;
			for(final Query query:queries){
				p.set(query.row);
				buffer.clear().append(query.genome.id).tab();
				field(buffer, legacy ? 4 : 3).tab().append(p.parseInt(legacy ? 6 : 4)).tab().append(p.parseInt(legacy ? 7 : 5))
					.tab().append(strand(legacy ? 8 : 6)).tab().append(MEMBERS[query.member]).tab();
				field(buffer, legacy ? 9 : 15).tab().append(OUTCOMES[query.bucket]).tab();
				if(legacy){buffer.append("NA");}
				else{field(buffer, 14);}
				buffer.tab();
				field(buffer, legacy ? 11 : 10);
				allQueries.writer.print(buffer.nl());
				if(query.bucket!=0){continue;}
				referenceKey(query.genome, legacy ? 11 : 10);
				final Reference ref=byReference.get(probe);
				require(ref!=null && ref.scored && !ref.consumed, "Every TP must consume exactly one declared scored reference");
				ref.consumed=true;
				final int cs=p.parseInt(legacy ? 6 : 4), ce=p.parseInt(legacy ? 7 : 5);
				final char strand=strand(legacy ? 8 : 6);
				require(p.termEquals(ref.seqid, legacy ? 4 : 3) && strand==ref.strand, "TP query and recorded reference must share contig and strand");
				final int d5=p.parseInt(legacy ? 12 : 16), d3=p.parseInt(legacy ? 13 : 17);
				require(cs>0 && ce>=cs && d5==(strand=='-' ? ref.stop-ce : cs-ref.start)
					&& d3==(strand=='-' ? ref.start-cs : ce-ref.stop), "Stored gene-sense offsets disagree with actual matched coordinates");
				if(d5==0){exact5++;}
				if(d3==0){exact3++;}
				if(d5==0 && d3==0){both++;}
				match(matches, "cobalamin", ref, cs, ce, d5, d3);
				match(matches, MEMBERS[query.member], ref, cs, ce, d5, d3);
				buffer.clear().append(ref.genome.id).tab().append(ref.id).tab().append(MEMBERS[query.member]).tab();
				field(buffer, legacy ? 9 : 15).tab().append(ref.seqid).tab().append(cs).tab().append(ce)
					.tab().append(ref.start).tab().append(ref.stop).tab().append(ref.strand).tab().append(d5).tab().append(d3).nl();
				witness.writer.print(buffer);
			}
			for(final Reference ref:references){
				allReferences.writer.print(buffer.clear().append(ref.genome.id).tab().append(ref.id).tab().append(ref.seqid)
					.tab().append(ref.start).tab().append(ref.stop).tab().append(ref.strand).tab().append(ref.scored).nl());
			}
			require(exact5==expected5 && exact3==expected3 && both==expectedBoth, "Exported TP endpoints disagree with the native joint summary");
			for(final Genome genome:genomes){
				final long[] c=genome.counts;
				inputs.writer.print(buffer.clear().append(genome.id).tab().append(c[0]+c[1]+c[2]+c[3]).tab().append(c[0]).tab().append(c[1]).tab().append(c[2]+c[3]).nl());
			}
			for(int member=0; member<3; member++){
				final long[] c=members[member];
				counts.writer.print(buffer.clear().append(MEMBERS[member]).tab().append(c[0]+c[1]+c[2]+c[3])
					.tab().append(c[0]).tab().append(c[1]).tab().append(c[2]).tab().append(c[3])
					.append("\tEMITTED_MEMBER_PARTITION_OF_UNION_OUTCOMES_NOT_SUBTYPE_ACCURACY\n"));
			}
		}
		try(Output metrics=new Output(out+".metrics.tsv")){
			final long[] sum=new long[4];
			for(final long[] member:members){for(int i=0; i<4; i++){sum[i]+=member[i];}}
			metrics.writer.println("inputs\tscored_references\tmasked_references\tcalls\tmatched\tunmatched\tmissed\tmasked_calls\tduplicate_queries\texact5\texact3\tboth_exact\tprecision\trecall\tqualification\tcpu_user_seconds\tcpu_system_seconds\twall_seconds\talign_rf00174\talign_rf01482\talign_rf01689");
			buffer.clear().append(genomes.size()).tab().append(scoredReferences).tab().append(maskedReferences)
				.tab().append(queries.size()).tab().append(sum[0]).tab().append(sum[2]+sum[3]).tab().append(missedReferences).tab().append(sum[1])
				.tab().append(duplicates).tab().append(expected5).tab().append(expected3).tab().append(expectedBoth).tab().append(precision).tab().append(recall)
				.tab().append(legacy ? "LEGACY_FULL_CM_UNION" : inexact ? "CANONICAL_NEW_CM_UNION_INEXACT" : "CANONICAL_NEW_CM_UNION");
			for(final String time:times){buffer.tab().append(time);}
			for(final long value:alignments){buffer.tab();if(value<0){buffer.append("NA");}else{buffer.append(value);}}
			metrics.writer.print(buffer.nl());
		}
	}

	/** Length-prefixes the input ID so repeated contig names cannot alias genomes. */
	private void match(final Output out, final String family, final Reference ref,
			final int cs, final int ce, final int d5, final int d3){
		assert(ref.consumed) : "Only the recorded one-to-one TP witness may enter an endpoint view";
		out.writer.print(buffer.clear().append(family).tab().append(ref.genome.id.length()).append(':').append(ref.genome.id).append('/').append(ref.seqid)
			.tab().append(ref.strand).tab().append(ref.start).tab().append(ref.stop).tab().append(cs).tab().append(ce).tab().append(d5).tab().append(d3).nl());
	}

	private Genome genome(final int column){
		probe.set(p, column);
		final Genome genome=byInput.get(probe);
		require(genome!=null, "Detailed row names an undeclared genome");
		return genome;
	}
	private void referenceKey(final Genome genome, final int column){
		key.clear().append(genome.id).tab();
		field(key, column);
		probe.set(key.array, 0, key.length());
	}
	private ByteBuilder field(final ByteBuilder dest, final int column){
		p.setBounds(column);
		return dest.append(p.line(), p.a(), p.b()-p.a());
	}
	private char strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Explicit plus/minus strand is required");
		return p.termEquals("-", column) ? '-' : '+';
	}
	private int member(final int column){
		for(int i=0; i<MEMBERS.length; i++){
			if(p.termStartsWith(PREFIXES[i], column)){return i;}
		}
		throw new IllegalArgumentException("Unknown emitted cobalamin member: "+p.parseString(column));
	}
	private boolean legacyUnion(){return p.termEquals("FULL_CM", 0) && p.termEquals("union", 1);}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------             State            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		Genome(final String id_){id=id_;}
		final String id;
		final long[] counts=new long[4];
	}
	private static final class Reference{
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
		boolean consumed;
		int member=-1, outcome=-1;
		long queryLine=-1;
	}
	private static final class Query{
		Query(final Genome genome_, final int member_, final int bucket_, final byte[] row_){
			genome=genome_;
			member=member_;
			bucket=bucket_;
			row=row_;
		}
		final Genome genome;
		final int member, bucket;
		final byte[] row;
	}
	private static final class ModelWins{
		ModelWins(final int member_, final int ordinal_, final String id_, final String header_){
			member=member_;
			ordinal=ordinal_;
			id=id_;
			header=header_;
		}
		final int member, ordinal;
		final String id, header;
		long wins;
	}
	private String[] modelFiles;
	private int[] modelCounts;
	private int expectedInputs=-1;
	private final ArrayList<ModelWins> modelWins=new ArrayList<ModelWins>();
	private String prefix, out, panel;
	private boolean legacy, inexact, memberMode;
	private String[] memberPrefixes;
	private long inexactMatches;
	private long duplicates, expected5, expected3, expectedBoth;
	private long scoredReferences, maskedReferences, missedReferences;
	private String precision, recall;
	private String timing, callerLog;
	private final String[] times={"NA", "NA", "NA"};
	private final long[] alignments={-1, -1, -1};
	private final long[][] members=new long[3][4];
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Query> queries=new ArrayList<Query>();
	private final ArrayList<Reference> references=new ArrayList<Reference>();
	private final ObjectMap<Target, Genome> byInput=new ObjectMap<Target, Genome>(16, Target.class, Genome.class);
	private final ObjectMap<Target, Reference> byReference=new ObjectMap<Target, Reference>(16, Target.class, Reference.class);
	private final ObjectSet<Target> seen=new ObjectSet<Target>(Target.class);
	private final Target probe=new Target();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder key=new ByteBuilder(), buffer=new ByteBuilder();
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] PREFIXES={"cobalamin_RF00174_", "cobalamin_RF01482_", "cobalamin_RF01689_"};
	private static final String[] OUTCOMES={"TP", "MASKED", "UNMATCHED_DUPLICATE", "UNMATCHED_OTHER"};
	private static final String[] ALIGNMENT_LABELS={"ncRNA Alignments (cobalamin_rf00174):", "ncRNA Alignments (cobalamin_rf01482):", "ncRNA Alignments (cobalamin_rf01689):"};
	private static final String LEGACY_GENOMES="view\tscope\tgenome\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\texact5\texact3\tboth_exact";
	private static final String LEGACY_CALLS="view\tscope\tcall_index\tgenome\tseqid\temitted_family\tstart1\tstop1\tstrand\tmodel\toutcome\ttruth_id\td5\td3";
	private static final String LEGACY_SUMMARY="view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact";
}
