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
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaSubtypeEvidence.Support;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Describes frozen assignments, never new matches. Keeps a stored representative
 * label separate from protected-prior/confident annotation support, and separates
 * a descriptive member-band prediction from the actual global filtered caller.
 * @author Brian Bushnell, Diona
 */
public final class SummarizeCobalaminMemberAssignments{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeCobalaminMemberAssignments.class, false);
		try{
			new SummarizeCobalaminMemberAssignments(pp.args).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeCobalaminMemberAssignments(final String[] args){
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("reference")){manifest=b;}
			else if(a.equals("raw")){raw=b;}
			else if(a.equals("grades")){grades=b.split(",", -1);}
			else if(a.equals("mapping")){mapping=b;}
			else if(a.equals("evidence")){evidencePaths=b;}
			else if(a.equals("choice")){choice=b;}
			else if(a.equals("turnover")){turnover=b;}
			else if(a.equals("losses")){losses=b;}
			else if(a.equals("member")){member=b;}
			else if(a.equals("nativepass")){nativePass=(float)NcrnaNumbers.doubleValue(a, b);}
			else if(a.equals("expectedinputs")){expectedInputs=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedreferences")){expectedReferences=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedmasked")){expectedMasked=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedmemberlost")){expectedMemberLost=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectednosame")){expectedNoSame=NcrnaNumbers.intValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown assignment diagnostic argument: "+arg);}
		}
		before=parser.in1;
		low=parser.in2;
		out=parser.out1;
		memberIndex=emitted(member);
		require(before!=null && low!=null && out!=null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1 && Shared.threads()==1 && memberIndex<3
			&& Float.isFinite(nativePass) && nativePass>=0 && nativePass<=1,
			"Use complete saved in/in2, one explicit member/native fallback, and fresh out=; no caller or matcher is invoked");
		require(expectedInputs>0 && expectedReferences>=0 && expectedMasked>=0 && expectedMemberLost>=0
			&& expectedNoSame>=0 && grades!=null && grades.length==3,
			"Declare complete input/reference/loss cohorts and baseline, lower-floor, actual-filter grade prefixes");
		require(Tools.testInputFiles(false, true, before, low, manifest, raw, mapping, evidencePaths, choice, turnover, losses),
			"All saved assignment and annotation inputs must exist");
		for(String prefix:grades){require(Tools.testInputFiles(false, true, prefix+".genomes.tsv"), "Complete genome accounting is required");}
		require(Tools.testInputFiles(false, true, grades[2]+".queries.tsv"), "Actual filtered queries are required");
		RrnaResourceIO.fresh(out, ".matrix.tsv", ".details.tsv", ".annotations.tsv", ".cohorts.tsv");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		readInputs();
		readReferences();
		readAnnotations();
		readPolicy();
		readGenomeCounts();
		try(Output details=new Output(out+".details.tsv")){
			detail=details;
			detail.writer.println("cohort\tinput_id\temitted_member\tmodel\treference_id\trepresentative_raw\trepresentative_label\tannotation_label\tconfident_members\tmarginal_members\tprotected_prior\tbaseline_reference_matched\tbaseline_same_member\tqualification");
			readScores(0, before);
			readScores(1, low);
			readActual();
			readGlobalLosses();
			finishMember();
		}
		for(Genome genome:genomes){
			for(int side=0; side<3; side++){
				for(int kind=0; kind<3; kind++){
					require(genome.expected[side][kind]==genome.actual[side][kind], "Every saved input/outcome count must reconcile: "+genome.id);
				}
			}
		}
		writeMatrices();
		System.err.println("COBALAMIN_MEMBER_ASSIGNMENTS_SUMMARIZED");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Reference Evidence     ----------------*/
	/*--------------------------------------------------------------*/

	private void readInputs(){
		try(Input input=new Input(manifest)){
			RrnaResourceIO.header(input.reader, "input_id\tfasta\tgff", manifest);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==3 && !p.termEquals("", 0) && !p.termEquals("", 2), "Reference manifest must bind exact input IDs to GFF metadata");
				final Genome genome=new Genome(p.parseString(0), p.parseString(2));
				require(!byInput.contains(genome.id), "Duplicate manifest input");
				byInput.put(genome.id, genome);
				genomes.add(genome);
			}
		}
		require(genomes.size()==expectedInputs, "Complete manifest input denominator changed");
	}

	private void readReferences(){
		int scored=0;
		int masked=0;
		try(Input input=new Input(raw)){
			RrnaResourceIO.header(input.reader, GradeNcrnaReference.RAW_HEADER, raw);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==18, "Raw grade schema changed");
				if(!p.termEquals("cobalamin", 2)){continue;}
				if(p.termEquals("QUERY", 1)){
					if(p.termEquals("TP", 7) || p.termEquals("TP_INEXACT", 7)){
						final Genome genome=genome(0);
						final long queryLine=p.parseLong(9);
						refKey(genome.id, Long.toString(queryLine));
						require(queryLine>0 && originalQueryLines.add(probe.copy()), "A raw query cannot be assigned to two references");
						refKey(genome.id, p.parseString(10));
						require(!originalCalls.contains(probe), "Raw matched references must have one original query assignment");
						originalCalls.put(probe.copy(), new OriginalCall(p.parseString(3), p.parseInt(4), p.parseInt(5), strand(6), RrnaResourceIO.first(p.parseString(15)), queryLine));
					}
					continue;
				}
				require(p.termEquals("REFERENCE", 1), "Unknown raw grade role");
				final Genome genome=genome(0);
				final String outcome=p.parseString(7);
				require(outcome.equals("TP") || outcome.equals("TP_INEXACT") || outcome.equals("FN") || outcome.equals("DONT_CARE"), "Unknown reference disposition");
				final Reference ref=new Reference(genome, p.parseString(8), p.parseString(3), p.parseInt(4), p.parseInt(5), strand(6),
					outcome.equals("DONT_CARE"), p.termEquals(".", 15) ? null : p.parseString(15), outcome.startsWith("TP"), outcome.startsWith("TP") ? p.parseLong(11) : 0);
				require(ref.start>0 && ref.stop>=ref.start && !ref.id.equals(".") && !ref.id.isEmpty(), "References require complete immutable identities and geometry");
				refKey(genome.id, ref.id);
				require(!byReference.contains(probe), "Duplicate input-qualified reference");
				byReference.put(probe.copy(), ref);
				references.add(ref);
				genome.refs.add(ref);
				if(ref.mask){masked++;}else{scored++;}
			}
		}
		require(scored==expectedReferences && masked==expectedMasked, "Frozen scored/masked reference denominators changed");
		int matched=0;
		for(Reference ref:references){
			refKey(ref.genome.id, ref.id);
			final OriginalCall call=originalCalls.get(probe);
			require(ref.oldMatched==(call!=null) && (!ref.oldMatched || ref.queryLine==call.line), "Raw query/reference assignments must be reciprocal");
			if(ref.oldMatched){matched++;}
		}
		require(matched==originalCalls.size(), "An original matched query cannot lack its reference row");
	}

	private void readAnnotations(){
		final ObjectMap<String, Evidence> sources=new ObjectMap<String, Evidence>(String.class, Evidence.class);
		final ArrayList<Evidence> order=new ArrayList<Evidence>();
		try(Input input=new Input(evidencePaths)){
			RrnaResourceIO.header(input.reader, "witnesses\tevents", evidencePaths);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==2, "Evidence sources explicitly bind witness and event tables");
				final Evidence source=new Evidence(p.parseString(0), p.parseString(1));
				require(!sources.contains(source.witnesses) && Tools.testInputFiles(false, true, source.witnesses, source.events),
					"Missing, duplicate or unreadable required sidecars must fail, not become UNKNOWN");
				sources.put(source.witnesses, source);
				order.add(source);
			}
		}
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(Genome genome:genomes){
			if(genome.refs.isEmpty()){continue;}
			final ArrayList<Row> rows=codec.readFamily(genome.gff, genome.id, "cobalamin");
			require(rows.size()==genome.refs.size(), "Saved and current canonical cobalamin populations differ: "+genome.id);
			for(Row row:rows){
				final Reference ref=reference(genome.id, row.id);
				require(ref!=null && ref.annotation==null && ref.seqid.equals(row.locus.seqid) && ref.start==row.locus.start
					&& ref.stop==row.locus.stop && ref.strand==row.locus.strand && ref.mask==row.mask
					&& same(ref.representative, row.model), "Canonical annotation identity/geometry/status/representative changed");
				ref.annotation=row;
				final Evidence source=sources.get(row.evidence);
				require(source!=null, "Canonical GFF names an undeclared required evidence source: "+row.evidence);
				source.refs.add(ref);
				if(!source.seen.contains(genome.id)){
					source.seen.put(genome.id, genome);
					source.inputs.add(genome.id);
				}
			}
		}
		try(Output output=new Output(out+".annotations.tsv")){
			output.writer.println("input_id\treference_id\tscored\trepresentative_raw\trepresentative_label\tannotation_label\tconfident_members\tmarginal_members\tprotected_prior\tmeasurement_models\twitnesses");
			for(Evidence source:order){
				if(source.refs.isEmpty()){continue;}
				final NcrnaSubtypeEvidence evidence=new NcrnaSubtypeEvidence(mapping, source.witnesses, source.events, source.inputs);
				require(evidence.size()==3, "This diagnostic requires the explicit three-member annotation axis");
				for(int i=0; i<3; i++){require(MEMBERS[i].equals(evidence.name(i)) && MEMBERS[i].equals(evidence.referenceModel(i)), "Annotation axis order changed");}
				for(Reference ref:source.refs){
					final Support support=evidence.describe(ref.annotation);
					ref.supported=ref.mask ? 0 : support.scored;
					ref.marginal=support.marginal;
					ref.protectedPrior=support.priorType>=0;
					ref.representativeLabel=representativeLabel(ref.representative);
					ref.annotationLabel=label(ref.supported);
					output.writer.print(row.clear().append(ref.genome.id).tab().append(ref.id).tab().append(!ref.mask)
						.tab().append(text(ref.representative)).tab().append(LABELS[ref.representativeLabel])
						.tab().append(LABELS[ref.annotationLabel]).tab().append(members(ref.supported))
						.tab().append(members(ref.marginal)).tab().append(ref.protectedPrior)
						.tab().append(ref.annotation.scores==null ? "NA" : ref.annotation.scores.models).tab().append(source.witnesses).nl());
				}
				evidence.finish();
			}
		}
		for(Reference ref:references){require(ref.annotation!=null && ref.annotationLabel>=0, "Every frozen reference needs a bound annotation description");}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Assignments      ----------------*/
	/*--------------------------------------------------------------*/

	private void readPolicy(){
		int seen=0;
		try(Input input=new Input(choice)){
			p.set(input.reader.nextLine());
			require(p.terms()==27 && p.termEquals("low_identity", 0) && p.termEquals("cm_minus_ga", 2), "Use the completed member choice table");
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==27 && p.termEquals("true", 4), "Choice rows must describe the selected saved point");
				final int split=split(3);
				require((seen&(1<<split))==0, "Repeated choice population");
				final float lo=p.parseFloat(0);
				final float hi=p.parseFloat(1);
				final double cutoff=p.parseDouble(2);
				final String status=p.parseString(25);
				require(Float.isFinite(lo) && Float.isFinite(hi) && lo>=0 && hi>=lo && hi<=1 && Double.isFinite(cutoff)
					&& (status.equals("JOINT_IMPROVEMENT") || status.equals("NO_JOINT_IMPROVEMENT")), "Preserve a finite saved choice and its recommendation status");
				if(seen==0){
					floor=lo;
					upper=hi;
					margin=cutoff;
					policyStatus=status;
				}
				require(floor==lo && upper==hi && margin==cutoff && policyStatus.equals(status), "Choice populations disagree");
				seen|=1<<split;
			}
		}
		require(seen==7, "Selection, judging and pooled choice rows are required");
		seen=0;
		try(Input input=new Input(turnover)){
			RrnaResourceIO.header(input.reader, TURNOVER_HEADER, turnover);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11 && p.termEquals(member, 0), "Turnover must bind exactly the declared member");
				final int split=split(1);
				require((seen&(1<<split))==0 && p.parseDouble(8)==margin && p.termEquals(policyStatus, 9), "Turnover must use the exact saved policy");
				for(int i=0; i<6; i++){expectedTurnover[split][i]=p.parseLong(i+2);}
				seen|=1<<split;
			}
		}
		require(seen==7 && expectedTurnover[2][3]==expectedMemberLost, "The requested member-loss cohort differs from the saved turnover");
	}

	private void readGenomeCounts(){
		for(int side=0; side<3; side++){
			int count=0;
			try(Input input=new Input(grades[side]+".genomes.tsv")){
				RrnaResourceIO.header(input.reader, GENOME_HEADER, grades[side]);
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==5, "Complete per-genome accounting is required");
					final Genome genome=genome(0);
					require(!genome.seen[side], "Duplicate genome accounting");
					genome.seen[side]=true;
					long total=0;
					for(int k=0; k<3; k++){
						genome.expected[side][k]=p.parseLong(k+2);
						require(genome.expected[side][k]>=0, "Negative query count");
						total+=genome.expected[side][k];
					}
					require(total==p.parseLong(1), "Query classes must partition every genome's calls");
					count++;
				}
			}
			require(count==expectedInputs, "Negative and positive genomes must both remain explicit");
		}
	}

	private void readScores(final int side, final String path){
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, ACCEPTANCE_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==27, "Saved caller score schema changed");
				final Genome genome=genome(1);
				final int group=p.parseInt(2);
				require(group>=0 && (genome.group<0 || genome.group==group) && p.termEquals((group&1)==0 ? "SELECTION" : "JUDGING", 3), "An input retains its original source-group split");
				genome.group=group;
				final int kind=kind(10);
				final Call call=call(genome, 8, 9, 4, 6, 7, 5, 12, kind);
				genome.actual[side][kind]++;
				if(kind==0){
					require(!call.ref.scoreSeen[side], "Duplicate assigned reference in saved calls");
					call.ref.scoreSeen[side]=true;
					if(side==0){
						refKey(genome.id, call.ref.id);
						final OriginalCall original=originalCalls.get(probe);
						require(original!=null && original.model.equals(call.model) && p.termEquals(original.seqid, 4)
							&& original.start==p.parseInt(6) && original.stop==p.parseInt(7) && original.strand==strand(5),
							"Original score/model/geometry must equal the recorded raw query assignment");
						call.ref.before=call;
					}
				}
				if(side==0 || call.member!=memberIndex){continue;}
				final float identity=p.parseFloat(13);
				require(Float.isFinite(identity) && identity>=0 && identity<=1, "Missing or invalid saved identity");
				final float bits=p.termEquals("NA", 15) ? Float.NaN : p.termEquals("-Infinity", 15) ? Float.NEGATIVE_INFINITY : p.parseFloat(15);
				final float ga=p.parseFloat(16);
				require((Float.isNaN(bits) || Float.isFinite(bits) || bits==Float.NEGATIVE_INFINITY) && Float.isFinite(ga) && ga>0, "Invalid saved score or GA");
				final boolean kept=identity>=upper || identity>=floor && (Float.isNaN(bits) ? identity>=nativePass : (double)bits-ga>=margin);
				if(!kept){continue;}
				if(kind==0){
					call.ref.memberKept=call;
					emit(1, call);
				}else if(kind==2){
					emit(2, call);
					memberUnmatched[group&1]++;
				}
				else{memberMasked++;}
			}
		}
		if(side==0){
			for(Reference ref:references){require(ref.oldMatched==(ref.before!=null), "Saved baseline score assignments differ from the frozen raw grade");}
		}
	}

	private void readActual(){
		try(Input input=new Input(grades[2]+".queries.tsv")){
			RrnaResourceIO.header(input.reader, QUERY_HEADER, grades[2]);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==10, "Actual filtered-query schema changed");
				final Genome genome=genome(0);
				final int kind=kind(7);
				final Call call=call(genome, 5, 6, 1, 2, 3, 4, 9, kind);
				genome.actual[2][kind]++;
				if(kind==0){
					require(p.termEquals("true", 8) && !call.ref.actualMatched, "Actual matched references must be eligible and unique");
					call.ref.actualMatched=true;
					emit(4, call);
				}else if(kind==2){emit(5, call);}
				else{actualMasked++;}
			}
		}
	}

	private void readGlobalLosses(){
		try(Input input=new Input(losses)){
			RrnaResourceIO.header(input.reader, LOSS_HEADER, losses);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==14, "Saved global-loss diagnostic schema changed");
				final Genome genome=genome(0);
				final Reference ref=reference(genome.id, p.parseString(1));
				require(ref!=null && ref.before!=null && !ref.actualMatched && !ref.lossSeen
					&& p.termEquals(ref.before.rawMember, 2) && p.termEquals(ref.before.model, 3),
					"Global loss evidence must bind its original exact model/reference assignment");
				ref.lossSeen=true;
				globalLost++;
				if(p.termEquals("NO_SAME_CM_MEASUREMENT", 8)){
					ref.noSame=true;
					emit(3, ref.before);
				}
			}
		}
		for(Reference ref:references){require(ref.lossSeen==(ref.before!=null && !ref.actualMatched), "Full global losses must reconcile before taking the NO_SAME subset");}
		require(totals[3]==expectedNoSame, "The requested NO_SAME loss subset changed");
	}

	private void finishMember(){
		final long[][] actual=new long[2][6];
		for(Reference ref:references){
			final boolean old=ref.before!=null && ref.before.member==memberIndex;
			final boolean kept=ref.memberKept!=null;
			if(!old && !kept){continue;}
			final int split=ref.genome.group&1;
			if(old){actual[split][0]++;}
			if(kept){
				actual[split][1]++;
				actual[split][old ? 2 : 4]++;
			}else{
				actual[split][3]++;
				emit(0, ref.before);
				if(ref.noSame){lossOverlap++;}
			}
		}
		for(int split=0; split<2; split++){actual[split][5]=memberUnmatched[split];}
		for(int split=0; split<3; split++){
			for(int field=0; field<6; field++){
				final long value=split==2 ? actual[0][field]+actual[1][field] : actual[split][field];
				require(value==expectedTurnover[split][field], "Exact member-local assignments disagree with the saved turnover");
			}
		}
		require(totals[0]==expectedMemberLost, "Member loss denominator changed");
	}

	private Call call(final Genome genome, final int memberColumn, final int modelColumn, final int seqidColumn,
			final int startColumn, final int stopColumn, final int strandColumn, final int refColumn, final int kind){
		final String rawMember=p.parseString(memberColumn);
		final String model=RrnaResourceIO.first(p.parseString(modelColumn));
		final int emitted=emitted(rawMember);
		if(emitted<3){require(model.startsWith("cobalamin_"+MEMBERS[emitted]+"_"), "Emitted model and member disagree");}
		Reference ref=null;
		if(kind==0){
			ref=reference(genome.id, p.parseString(refColumn));
			final int start=p.parseInt(startColumn);
			final int stop=p.parseInt(stopColumn);
			require(ref!=null && !ref.mask && p.termEquals(ref.seqid, seqidColumn) && ref.strand==strand(strandColumn)
				&& start>0 && stop>=start && Math.min(stop, ref.stop)>=Math.max(start, ref.start),
				"A saved TP must retain its exact input/reference binding and positive same-strand overlap; no rematching");
		}else{require(p.termEquals(".", refColumn) || p.termEquals("NA", refColumn), "Unmatched/masked queries have no matched reference");}
		return new Call(genome, emitted, rawMember, model, ref);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void emit(final int cohort, final Call call){
		final Reference ref=call.ref;
		counts[0][cohort][call.member][ref==null ? NO_MATCH : ref.representativeLabel]++;
		counts[1][cohort][call.member][ref==null ? NO_MATCH : ref.annotationLabel]++;
		totals[cohort]++;
		row.clear().append(COHORTS[cohort]).tab().append(call.genome.id).tab().append(call.rawMember).tab().append(call.model).tab();
		if(ref==null){row.append("NA\tNA\tNO_MATCHED_REFERENCE\tNO_MATCHED_REFERENCE\tNA\tNA\tfalse\tfalse\tfalse");}
		else{
			row.append(ref.id).tab().append(text(ref.representative)).tab().append(LABELS[ref.representativeLabel])
				.tab().append(LABELS[ref.annotationLabel]).tab().append(members(ref.supported)).tab().append(members(ref.marginal))
				.tab().append(ref.protectedPrior).tab().append(ref.before!=null)
				.tab().append(ref.before!=null && ref.before.member==call.member);
		}
		detail.writer.print(row.tab().append(QUALIFICATION).nl());
	}

	private void writeMatrices(){
		try(Output output=new Output(out+".matrix.tsv")){
			output.writer.println("axis\tcohort\temitted_member\tRF00174\tRF01482\tRF01689\tAMBIGUOUS\tUNKNOWN\tNO_MATCHED_REFERENCE\ttotal\tqualification");
			for(int axis=0; axis<2; axis++){
				for(int cohort=0; cohort<6; cohort++){
					long cohortTotal=0;
					for(int emitted=0; emitted<4; emitted++){
						final long[] values=counts[axis][cohort][emitted];
						long total=0;
						for(long value:values){total+=value;}
						if(total==0){continue;}
						row.clear().append(axis==0 ? "STORED_REPRESENTATIVE" : "ANNOTATION_SUPPORT").tab().append(COHORTS[cohort])
							.tab().append(emitted==3 ? "UNKNOWN" : MEMBERS[emitted]);
						for(long value:values){row.tab().append(value);}
						output.writer.print(row.tab().append(total).tab().append(QUALIFICATION).nl());
						cohortTotal+=total;
					}
					require(cohortTotal==totals[cohort], "Every record must occur exactly once on each annotation axis");
				}
			}
		}
		try(Output output=new Output(out+".cohorts.tsv")){
			output.writer.println("cohort\trecords\tselection_status\tqualification");
			for(int i=0; i<6; i++){
				output.writer.print(row.clear().append(COHORTS[i]).tab().append(totals[i]).tab()
					.append(i<3 ? policyStatus : "ACTUAL_GLOBAL_FILTER_CALLER").tab().append(QUALIFICATION).nl());
			}
			output.writer.println("GLOBAL_ALL_LOST_REFERENCES\t"+globalLost+"\tACTUAL_GLOBAL_FILTER_CALLER\tNO_SAME_IS_A_SUBSET_NOT_THE_FULL_LOSS_SET");
			output.writer.println("LOSS_COHORT_INTERSECTION\t"+lossOverlap+"\tDESCRIPTIVE\tCOHORTS_OVERLAP_DO_NOT_SUM");
			output.writer.println("MEMBER_MASKED_EXCLUDED\t"+memberMasked+"\t"+policyStatus+"\tMASKED_IS_NOT_UNMATCHED");
			output.writer.println("GLOBAL_MASKED_EXCLUDED\t"+actualMasked+"\tACTUAL_GLOBAL_FILTER_CALLER\tMASKED_IS_NOT_UNMATCHED");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private Genome genome(final int column){
		final Genome result=byInput.get(p.parseString(column));
		require(result!=null, "Saved row names an undeclared input");
		return result;
	}
	private Reference reference(final String input, final String id){
		refKey(input, id);
		return byReference.get(probe);
	}
	private void refKey(final String input, final String id){
		key.clear().append(input).tab().append(id);
		probe.set(key.array, 0, key.length());
	}
	private int kind(final int column){
		if(p.termEquals("TP", column)){return 0;}
		if(p.termEquals("MASKED", column)){return 1;}
		require(p.termEquals("UNMATCHED_OTHER", column) || p.termEquals("UNMATCHED_DUPLICATE", column), "Unknown query outcome");
		return 2;
	}
	private int split(final int column){
		if(p.termEquals("SELECTION", column)){return 0;}
		if(p.termEquals("JUDGING", column)){return 1;}
		require(p.termEquals("ALL", column), "Unknown saved split");
		return 2;
	}
	private char strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Strand must be explicit");
		return p.termEquals("+", column) ? '+' : '-';
	}
	private static int emitted(final String value){
		for(int i=0; i<3; i++){if(MEMBERS[i].equals(value)){return i;}}
		return 3;
	}
	private static int representativeLabel(final String value){
		if(value==null || value.equals(".") || value.equals("NA")){return UNKNOWN;}
		int bits=0;
		for(String token:value.split("[,;|]", -1)){
			final int member=emitted(token.trim());
			if(member==3){return UNKNOWN;}
			bits|=1<<member;
		}
		return label(bits);
	}
	private static int label(final int bits){return bits==0 ? UNKNOWN : Integer.bitCount(bits)>1 ? AMBIGUOUS : Integer.numberOfTrailingZeros(bits);}
	/** The explicit three-member mapping limits support masks to 0..7. */
	private static String members(final int bits){return MEMBER_SETS[bits];}
	private static String text(final String value){return value==null ? "NA" : value;}
	private static boolean same(final String a, final String b){return a==null ? b==null : a.equals(b);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		Genome(String id_, String gff_){
			id=id_;
			gff=gff_;
		}
		final String id;
		final String gff;
		final ArrayList<Reference> refs=new ArrayList<Reference>();
		final long[][] expected=new long[3][3];
		final long[][] actual=new long[3][3];
		final boolean[] seen=new boolean[3];
		int group=-1;
	}
	private static final class Reference{
		Reference(Genome genome_, String id_, String seqid_, int start_, int stop_, char strand_, boolean mask_, String representative_, boolean oldMatched_, long queryLine_){
			genome=genome_;
			id=id_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			mask=mask_;
			representative=representative_;
			oldMatched=oldMatched_;
			queryLine=queryLine_;
		}
		final Genome genome;
		final String id;
		final String seqid;
		final int start;
		final int stop;
		final char strand;
		final boolean mask;
		final String representative;
		final boolean oldMatched;
		final long queryLine;
		Row annotation;
		Call before;
		Call memberKept;
		final boolean[] scoreSeen=new boolean[2];
		boolean actualMatched;
		boolean lossSeen;
		boolean noSame;
		boolean protectedPrior;
		int supported;
		int marginal;
		int representativeLabel=-1;
		int annotationLabel=-1;
	}
	private static final class Call{
		Call(Genome genome_, int member_, String rawMember_, String model_, Reference ref_){
			genome=genome_;
			member=member_;
			rawMember=rawMember_;
			model=model_;
			ref=ref_;
		}
		final Genome genome;
		final int member;
		final String rawMember;
		final String model;
		final Reference ref;
	}
	private static final class Evidence{
		Evidence(String witnesses_, String events_){
			witnesses=witnesses_;
			events=events_;
		}
		final String witnesses;
		final String events;
		final ArrayList<String> inputs=new ArrayList<String>();
		final ArrayList<Reference> refs=new ArrayList<Reference>();
		final ObjectMap<String, Genome> seen=new ObjectMap<String, Genome>(String.class, Genome.class);
	}
	private static final class OriginalCall{
		OriginalCall(String seqid_, int start_, int stop_, char strand_, String model_, long line_){
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			model=model_;
			line=line_;
		}
		final String seqid;
		final int start;
		final int stop;
		final char strand;
		final String model;
		final long line;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private String before;
	private String low;
	private String[] grades;
	private String manifest;
	private String raw;
	private String mapping;
	private String evidencePaths;
	private String choice;
	private String turnover;
	private String losses;
	private String out;
	private String member="RF01482";
	private int memberIndex;
	private float nativePass=Float.NaN;
	private float floor;
	private float upper;
	private double margin;
	private String policyStatus;
	private int expectedInputs=-1;
	private int expectedReferences=-1;
	private int expectedMasked=-1;
	private int expectedMemberLost=-1;
	private int expectedNoSame=-1;
	private long memberMasked;
	private long actualMasked;
	private long globalLost;
	private long lossOverlap;
	private Output detail;
	private final ObjectMap<String, Genome> byInput=new ObjectMap<String, Genome>(String.class, Genome.class);
	private final ObjectMap<Target, Reference> byReference=new ObjectMap<Target, Reference>(Target.class, Reference.class);
	private final ObjectMap<Target, OriginalCall> originalCalls=new ObjectMap<Target, OriginalCall>(Target.class, OriginalCall.class);
	private final ObjectSet<Target> originalQueryLines=new ObjectSet<Target>(Target.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Reference> references=new ArrayList<Reference>();
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder row=new ByteBuilder();
	private final long[][][][] counts=new long[2][6][4][6];
	private final long[] totals=new long[6];
	private final long[] memberUnmatched=new long[2];
	private final long[][] expectedTurnover=new long[3][6];

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int AMBIGUOUS=3;
	private static final int UNKNOWN=4;
	private static final int NO_MATCH=5;
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] MEMBER_SETS={"NA", "RF00174", "RF01482", "RF00174,RF01482",
		"RF01689", "RF00174,RF01689", "RF01482,RF01689", "RF00174,RF01482,RF01689"};
	private static final String[] LABELS={"RF00174", "RF01482", "RF01689", "AMBIGUOUS", "UNKNOWN", "NO_MATCHED_REFERENCE"};
	private static final String[] COHORTS={"MEMBER_PREDICTED_LOST_BASELINE", "MEMBER_PREDICTED_ALL_MATCHES", "MEMBER_PREDICTED_UNMATCHED", "GLOBAL_ACTUAL_LOST_NO_SAME_SUBSET", "GLOBAL_ACTUAL_ALL_MATCHES", "GLOBAL_ACTUAL_UNMATCHED"};
	private static final String QUALIFICATION="FROZEN_ONE_TO_ONE_ASSIGNMENTS;ALL_MATCHES_INCLUDE_GAINS;ANNOTATION_SUPPORT_NOT_SUBTYPE_TRUTH;NO_REMATCH_OR_REFIT";
	public static final String ACCEPTANCE_HEADER="arm\tinput_id\tresolved_group\tsplit\tcontig\tstrand\tstart\tstop\tmember\tmodel\toutcome\tmatched_inexact\treference_id\taccepted_identity\tendpoint_identity\tcm_bits\tga_bits\tcm_minus_ga\tcm_over_ga\tcm_decision\tnative_score_status\tnative_cm_bits\tnative_ga_bits\tanchored_minus_native\thbm_score\thbm_status\tboundary_different";
	public static final String QUERY_HEADER="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id";
	public static final String GENOME_HEADER="input_id\tcalls\tmatched\tmasked\tunmatched";
	public static final String TURNOVER_HEADER="member\tsplit\tbaseline_matched\tmatched_kept\told_references_retained\told_references_lost\tnew_references_gained\tunmatched_kept\tcm_minus_ga\tselection_status\tqualification";
	public static final String LOSS_HEADER="input_id\treference_id\tmember\tmodel\tlast_detection_genome\tidentity\tanchored_bits\tanchored_ga\tnative_status\tnative_bits\tnative_ga\tours_would_pass\tnative_would_pass\tclassification";
}
