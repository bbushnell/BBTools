package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaReferenceMatcher.Result;
import prok.NcrnaSubtypeEvidence.Support;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;

/** Optional subtype reports, separate from the unchanged union grading path.
 * Both assignments use the same loaded per-input queries and support. Exclusive
 * masks are applied BEFORE matching; filtering recovery outcomes would produce
 * false misses when a query can be reassigned to a single-supported reference.
 * Instances are thread-confined. Outputs are borrowed and closed by the main
 * grader before summary publication. No GFF or global sidecar is reread here.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaSubtypeGrade{

	NcrnaSubtypeGrade(final String mapping, final String witnesses, final String events,
		final Rule[] rules, final List<String> inputs){
		evidence=new NcrnaSubtypeEvidence(mapping, witnesses, events, rules, inputs);
		stats=new Stats[2][evidence.size()];
		for(int view=0; view<2; view++){
			for(int type=0; type<evidence.size(); type++){stats[view][type]=new Stats();}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Matching           ----------------*/
	/*--------------------------------------------------------------*/

	/** Derives support once and runs each scoring view over one input's loci. */
	void grade(final String input, final List<Row> references, final List<Call> calls,
		final float overlap, final Output raw, final Output supportOutput){
		assert(raw!=null && supportOutput!=null) : "Subtype detail writers are owned by the main grader's checked output scope";
		final ArrayList<Row> parents=new ArrayList<Row>();
		final ArrayList<Support> support=new ArrayList<Support>();
		for(final Row row:references){
			if(!NcrnaSubtypeEvidence.FAMILY.equals(row.locus.family)){continue;}
			final Support item=evidence.support(row);
			parents.add(row);
			support.add(item);
			if(row.mask){maskedParents++;}
			else{
				scoredParents++;
				if(item.multi()){multiParents++;}
				else{singleParents++;}
			}
			if(item.priorType>=0){protectedParents++;}
			writeSupport(row, item, supportOutput);
		}
		for(int type=0; type<evidence.size(); type++){
			final ArrayList<Call> selected=new ArrayList<Call>();
			final ArrayList<Locus> queryLoci=new ArrayList<Locus>();
			for(final Call call:calls){
				if(call.included && evidence.queryType(call.ruleIndex)==type){
					selected.add(call);
					queryLoci.add(call.locus);
				}
			}
			for(int view=0; view<2; view++){
				gradeView(input, parents, support, selected, queryLoci, type, view, overlap, raw);
			}
		}
	}

	private void gradeView(final String input, final ArrayList<Row> parents,
		final ArrayList<Support> support, final ArrayList<Call> calls, final ArrayList<Locus> queries,
		final int type, final int view, final float overlap, final Output out){
		assert(parents.size()==support.size() && queries.size()==calls.size())
			: "Subtype projection preserves original row indices and authoritative parent support";
		final Stats count=stats[view][type];
		final ArrayList<Locus> truth=new ArrayList<Locus>();
		final IntList original=new IntList();
		final boolean[] temporaryMasks=new boolean[parents.size()];
		for(int i=0; i<parents.size(); i++){
			final Row row=parents.get(i);
			final Support item=support.get(i);
			if(!NcrnaSubtypeEvidence.includedInView(row, item, type)){continue;}
			final boolean mask=NcrnaSubtypeEvidence.maskedInView(row, item, type, view==EXCLUSIVE);
			temporaryMasks[truth.size()]=mask;
			truth.add(row.locus);
			original.add(i);
			if(mask){count.maskReferences++;}
			else{count.scoredReferences++;}
		}
		final boolean[] masks=Arrays.copyOf(temporaryMasks, truth.size());
		final Result match=NcrnaReferenceMatcher.match(queries, ones(queries.size()), truth, ones(truth.size()), masks, overlap);
		count.queries+=calls.size();
		count.tp+=match.matched;
		count.fp+=match.falsePositives;
		count.fn+=match.falseNegatives;
		count.maskCalls+=match.maskedCalls;
		count.duplicates+=match.duplicates;
		for(int i=0; i<calls.size(); i++){
			final Call call=calls.get(i);
			final int paired=match.callToReference[i], masked=match.callMaskReference[i];
			final int witness=paired>=0 ? paired : masked;
			final Row row=witness<0 ? null : parents.get(original.get(witness));
			final Support item=witness<0 ? null : support.get(original.get(witness));
			start(view, type, input, "QUERY");
			buffer.append(".\t").append(call.line).tab();
			field(call.model).tab();
			locus(call.locus);
			buffer.append(paired>=0 ? "TP" : masked>=0 ? "MASKED" : "FP").tab();
			field(paired<0 ? "." : row.id).append("\t.\t");
			field(masked<0 ? "." : row.id).tab().append(match.duplicateCalls[i]).tab().append(match.scoredEligible[i]).tab();
			supportFields(item);
			if(paired<0){buffer.append(".\t.\n");}
			else{
				final boolean minus=call.locus.strand=='-';
				final int d5=minus ? row.locus.stop-call.locus.stop : call.locus.start-row.locus.start;
				final int d3=minus ? row.locus.start-call.locus.start : call.locus.stop-row.locus.stop;
				buffer.append(d5).tab().append(d3).nl();
				if(d5==0){count.exact5++;}
				if(d3==0){count.exact3++;}
				if(d5==0 && d3==0){count.exactBoth++;}
			}
			flush(out);
		}
		for(int i=0; i<truth.size(); i++){
			final Row row=parents.get(original.get(i));
			start(view, type, input, "REFERENCE");
			field(row.id).append("\t.\t");
			field(row.model==null ? "." : row.model).tab();
			locus(row.locus);
			buffer.append(masks[i] ? "DONT_CARE" : match.referenceToCall[i]>=0 ? "TP" : "FN").append("\t.\t");
			if(match.referenceToCall[i]<0){buffer.append('.');}
			else{buffer.append(calls.get(match.referenceToCall[i]).line);}
			buffer.append("\t.\t.\t.\t");
			supportFields(support.get(original.get(i)));
			buffer.append(".\t.\n");
			flush(out);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	/** Called after checked detail closure; evidence or count failure prevents summary. */
	void writeSummary(final String prefix, final int inputs, final float overlap){
		evidence.finish();
		require(singleParents+multiParents==scoredParents, "Distinct single/multi-support parents conserve the scored union population");
		try(Output out=new Output(prefix+SUFFIXES[3])){
			out.writer.println("inputs\tscored_parents\tsingle_supported_parents\tmulti_supported_parents\tmasked_parents\tprotected_parents\tqualification");
			buffer.append(inputs).tab().append(scoredParents).tab().append(singleParents).tab().append(multiParents)
				.tab().append(maskedParents).tab().append(protectedParents).append("\tDISTINCT_PARENT_COUNTS_NOT_SUMMED_SUBTYPE_MEMBERSHIPS\n");
			flush(out);
		}
		try(Output out=new Output(prefix+SUFFIXES[2])){
			out.writer.println(SUMMARY_HEADER);
			for(int view=0; view<2; view++){
				for(int type=0; type<evidence.size(); type++){
					final Stats count=stats[view][type];
					require(count.tp+count.fn==count.scoredReferences && count.tp+count.fp+count.maskCalls==count.queries,
						"Every subtype view must conserve scored references and original subtype queries");
					buffer.append(VIEWS[view]).tab();
					field(evidence.name(type)).tab().append(inputs).tab().append(count.tp).tab().append(count.fp).tab().append(count.fn)
						.tab().append(count.scoredReferences).tab().append(count.maskReferences).tab().append(count.maskCalls)
						.tab().append(count.duplicates).tab().append(count.exact5).tab().append(count.exact3).tab().append(count.exactBoth).tab();
					ratio(count.tp, count.tp+count.fp).tab();
					ratio(count.tp, count.tp+count.fn).tab().append(overlap, 6).tab().append(QUALIFICATIONS[view]).nl();
					flush(out);
				}
			}
		}
	}

	private void writeSupport(final Row row, final Support item, final Output out){
		field(row.input).tab();
		field(row.id).tab();
		locus(row.locus);
		buffer.append(row.mask ? "dont_care" : "scored").tab();
		field(row.model==null ? "." : row.model).tab().append(item.priorType>=0).tab();
		supportFields(item);
		buffer.setLength(buffer.length()-1);
		buffer.nl();
		flush(out);
	}
	private void start(final int view, final int type, final String input, final String role){
		assert(buffer.length()==0) : "Subtype rows flush before the next row starts";
		buffer.append(VIEWS[view]).tab();
		field(evidence.name(type)).tab();
		field(input).tab().append(role).tab();
	}
	private void locus(final Locus locus){
		field(locus.seqid).tab().append(locus.start).tab().append(locus.stop).tab().append(locus.strand).tab();
	}
	private void supportFields(final Support support){
		if(support==null){buffer.append(".\t.\t.\t");}
		else{
			bits(support.scored).tab();
			bits(support.marginal).tab().append(support.multi()).tab();
		}
	}
	private ByteBuilder bits(final int bits){
		if(bits==0){return buffer.append('.');}
		boolean first=true;
		for(int i=0; i<evidence.size(); i++){
			if((bits&(1<<i))==0){continue;}
			if(!first){buffer.append(',');}
			field(evidence.name(i));
			first=false;
		}
		return buffer;
	}
	private ByteBuilder ratio(final long numerator, final long denominator){
		return denominator==0 ? buffer.append("NA") : buffer.append(numerator/(double)denominator, 6);
	}
	private ByteBuilder field(final String text){return BuildNcrnaReference.field(buffer, text);}
	private void flush(final Output out){
		assert(buffer.length()>0) : "A subtype report write must contain a complete row";
		out.writer.print(buffer);
		buffer.clear();
	}
	private static int[] ones(final int size){
		final int[] ids=new int[size];
		Arrays.fill(ids, 1);
		return ids;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Stats{
		long queries, tp, fp, fn, scoredReferences, maskReferences, maskCalls, duplicates, exact5, exact3, exactBoth;
	}
	private final NcrnaSubtypeEvidence evidence;
	private final Stats[][] stats;
	private final ByteBuilder buffer=new ByteBuilder();
	private long scoredParents, singleParents, multiParents, maskedParents, protectedParents;
	private static final int EXCLUSIVE=1;
	static final String[] SUFFIXES={".subtype_raw.tsv", ".subtype_support.tsv", ".subtype_summary.tsv", ".subtype_counts.tsv"};
	static final String[] VIEWS={"MODEL_SUPPORT_RECOVERY", "EXCLUSIVE_SUBTYPE_ACCURACY"};
	static final String[] QUALIFICATIONS={"model-support recovery, not a missed gene, not wrong-subtype proof; NONADDITIVE", "single-subtype reference-support restriction, not unique-biological-subtype proof; NONADDITIVE"};
	static final String RAW_HEADER="view\tsubtype\tinput_id\trole\tparent_id\tquery_line\tmodel\tseqid\tstart\tstop\tstrand\toutcome\tmatched_parent_id\tmatched_query_line\tmask_parent_id\tduplicate\tscored_eligible\tscored_support\tmarginal_support\tmulti_supported\td5\td3";
	static final String SUPPORT_HEADER="input_id\tparent_id\tseqid\tstart\tstop\tstrand\tstatus\trepresentative_model\tprotected_prior\tscored_support\tmarginal_support\tmulti_supported";
	static final String SUMMARY_HEADER="view\tsubtype\tinputs\tTP\tFP\tFN\tscored_references\tmasked_references\tmasked_calls\tduplicate_queries\texact_5prime\texact_3prime\texact_both\tprecision\trecall\toverlap\tqualification";
}
