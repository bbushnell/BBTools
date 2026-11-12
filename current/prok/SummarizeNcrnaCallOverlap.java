package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Describes apparent false calls using a saved GradeNcrnaReference raw report.
 * Adapted from template.A_SampleByteFile: native parsing, byte streams and checked
 * closure. One input's selected-family rows are retained at a time. No sequences,
 * alignments, reference admission or changed grading decisions are produced.
 *
 * Coordinates are one-based inclusive. Text fields retain the original report's
 * percent escaping and UTF-8. Overlap fractions describe call coverage, not
 * reciprocal overlap. Distance is the number of intervening bases; negative
 * values mean an overlap with an originally matched call. Nearest distances and
 * model names do not establish a tandem or biological truth. Numeric verifier
 * identity is unavailable in this input format and is explicitly reported as NA.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaCallOverlap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final SummarizeNcrnaCallOverlap tool=new SummarizeNcrnaCallOverlap(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	public SummarizeNcrnaCallOverlap(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("family")){
				family=b;
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(input!=null && prefix!=null && family!=null && !family.isEmpty(),
			"Supply in=grade.raw.tsv family=exact_report_family out=fresh_prefix");
		familyBytes=family.getBytes(StandardCharsets.UTF_8);
		require(!parser.overwrite && !parser.append, "Audit outputs must be fresh; overwrite/append are unsupported");
		require(Tools.testInputFiles(false, true, input), "Saved raw report must be readable");
		require(Tools.testOutputFiles(false, false, false, prefix+".calls.tsv", prefix+".summary.tsv", prefix+".distances.tsv"),
			"Audit output paths must be fresh and writable");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		try(Input in=new Input(input); Output out=new Output(prefix+".calls.tsv")){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.RAW_HEADER, input);
			out.writer.println(CALL_HEADER);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==18, "Expected complete native raw grading rows");
				if(!matches(familyBytes, 2)){continue;}
				if(inputId==null || !matches(inputBytes, 0)){
					if(inputId!=null){report(out);}
					inputId=text(0);
					inputBytes=inputId.getBytes(StandardCharsets.UTF_8);
					require(!inputId.isEmpty() && seenInputs.get(inputId)<0,
						"Each exact input's family rows must be contiguous: "+inputId);
					seenInputs.put(inputId, seenInputs.size());
					queries.clear();
					references.clear();
				}
				final boolean query=parser.termEquals("QUERY", 1);
				require(query || parser.termEquals("REFERENCE", 1), "Unknown native report row role");
				final Row row=new Row(text(3), parser.parseInt(4), parser.parseInt(5),
					text(6), text(7), text(query ? 9 : 8), text(15), query && parser.termEquals("true", 13));
				if(query){
					require(row.outcome.equals("TP") || row.outcome.equals("FP") || row.outcome.equals("MASKED")
						|| row.outcome.equals("EXCLUDED"), "Unknown query outcome");
					if(row.outcome.equals("EXCLUDED")){continue;}
					queries.add(row);
				}else{
					require(row.outcome.equals("TP") || row.outcome.equals("FN") || row.outcome.equals("DONT_CARE"),
						"Unknown reference outcome");
					references.add(row);
				}
			}
			if(inputId!=null){report(out);}
		}
		writeSummary();
		timer.stop();
		outstream.println("Inspected "+falseCalls+" apparent false calls; original grades unchanged. Time: "+timer);
	}

	/** Quadratic comparisons are bounded to one original input and one family.
	 * Reference and query witnesses use coordinate/ID tie-breaking, not input order.
	 */
	private void report(final Output out){
		assert(inputId!=null) : "Report groups retain exact original input identity";
		final long priorTrue=trueCalls, priorRecovered=recoveredReferences;
		queries.sort(Row::compare);
		for(final Row ref:references){
			if(ref.outcome.equals("TP")){recoveredReferences++;}
			else if(ref.outcome.equals("FN")){missedReferences++;}
			else{maskedReferences++;}
		}
		for(final Row call:queries){
			if(call.outcome.equals("TP")){trueCalls++;continue;}
			if(call.outcome.equals("MASKED")){maskedCalls++;continue;}
			falseCalls++;
			if(call.duplicate){duplicateFalseCalls++;}
			Row best=null, mask=null, nearest=null;
			int overlap=0, maskOverlap=0;
			long nearestGap=Long.MAX_VALUE;
			for(final Row ref:references){
				final int n=intersection(call, ref);
				if(ref.outcome.equals("DONT_CARE")){
					if(better(ref, n, mask, maskOverlap)){mask=ref;maskOverlap=n;}
				}else if(better(ref, n, best, overlap)){best=ref;overlap=n;}
			}
			for(final Row other:queries){
				if(!other.outcome.equals("TP") || !sameTarget(call, other)){continue;}
				final long gap=(long)Math.max(call.start, other.start)-Math.min(call.stop, other.stop)-1;
				if(gap<nearestGap || gap==nearestGap && Row.compare(other, nearest)<0){
					nearest=other;
					nearestGap=gap;
				}
			}
			final String geometry;
			if(2L*overlap>=call.length()){
				overlappingHalf++;
				geometry="AT_LEAST_HALF_CALL_OVERLAP";
			}else if(overlap>0){
				partialOverlap++;
				geometry="POSITIVE_OVERLAP_BELOW_HALF";
			}else{
				zeroOverlap++;
				geometry="NO_SCORED_REFERENCE_OVERLAP";
				if(nearest!=null && nearestGap<0){overlapsMatchedCall++;}
				if(nearest!=null && nearestGap>=0){
					for(int i=0; i<DISTANCES.length; i++){
						if(nearestGap<=DISTANCES[i]){nearCounts[i]++;}
					}
				}
			}
			int eligibleCalls=0;
			if(best!=null){
				for(final Row other:queries){
					if(2L*intersection(other, best)>=other.length()){eligibleCalls++;}
				}
			}
			raw(buffer.clear(), inputId).tab();
			raw(buffer, family).tab();
			raw(buffer, call.seqid).tab();
			buffer.append(call.start).tab().append(call.stop).tab().append(call.strand).tab();
			raw(buffer, call.id).tab();
			raw(buffer, call.model).tab().append(call.duplicate).tab().append(geometry).tab();
			raw(buffer, best==null ? "." : best.id).tab().append(overlap).tab().append(overlap/(double)call.length(), 8).tab();
			buffer.append(best==null ? "." : best.outcome).tab().append(eligibleCalls).tab();
			raw(buffer, mask==null ? "." : mask.id).tab().append(maskOverlap).tab();
			raw(buffer, nearest==null ? "." : nearest.id).tab();
			if(nearest==null){buffer.append("NA\t.\t.\t.");}
			else{
				buffer.append(nearestGap).tab().append(nearest.start).tab().append(nearest.stop).tab();
				raw(buffer, nearest.model);
			}
			out.writer.print(buffer.append("\tNA\tGEOMETRY_ONLY_NOT_RECLASSIFIED\n"));
		}
		require(trueCalls-priorTrue==recoveredReferences-priorRecovered,
			"Original one-to-one query/reference TP rows must agree within each input: "+inputId);
	}

	/** Only positive physical intersection supplies an overlap witness. */
	private static boolean better(final Row row, final int overlap, final Row old, final int previous){
		return overlap>0 && (old==null || overlap>previous || overlap==previous && Row.compare(row, old)<0);
	}
	private static boolean sameTarget(final Row a, final Row b){return a.seqid.equals(b.seqid) && a.strand.equals(b.strand);}
	private static int intersection(final Row a, final Row b){
		return sameTarget(a, b) ? Math.max(0, Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1) : 0;
	}
	/** Already escaped native fields must not be escaped again. ASCII allocates
	 * nothing; unusual retained Unicode uses the same UTF-8 fallback as the writer.
	 */
	private static ByteBuilder raw(final ByteBuilder out, final String text){
		assert(text!=null) : "Absent report witnesses are explicit dots, never null";
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)>127){return out.append(text.getBytes(StandardCharsets.UTF_8));}
		}
		return out.append(text);
	}
	/** Text is retained only for the current input's result rows; no field boxing. */
	private String text(final int column){
		parser.setBounds(column);
		final String result=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(result.indexOf('\uFFFD')<0, "Malformed UTF-8 would change a report identity");
		return result;
	}
	/** Compare cached UTF-8 identities without a transient String for each row. */
	private boolean matches(final byte[] token, final int column){
		assert(token!=null) : "Group/family identity is cached before byte comparison";
		parser.setBounds(column);
		if(parser.b()-parser.a()!=token.length){return false;}
		for(int i=0; i<token.length; i++){
			if(token[i]!=parser.line()[parser.a()+i]){return false;}
		}
		return true;
	}

	private void writeSummary(){
		require(falseCalls==overlappingHalf+partialOverlap+zeroOverlap, "Overlap classes must conserve every original FP");
		try(Output out=new Output(prefix+".distances.tsv")){
			out.writer.println("family\tmaximum_intervening_bases\tzero_reference_overlap_FPs_near_original_TP\tqualification");
			for(int i=0; i<DISTANCES.length; i++){
				raw(buffer.clear(), family).tab().append(DISTANCES[i]).tab().append(nearCounts[i]);
				out.writer.print(buffer.append("\tCUMULATIVE_GEOMETRY_NOT_TANDEM_ADJUDICATION\n"));
			}
		}
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.println("family\tinputs_with_family_rows\toriginal_TP_calls\toriginal_FP_calls\toriginal_masked_calls\toriginal_recovered_references\toriginal_missed_references\tmasked_references\tFP_at_least_half_call_overlap\tFP_positive_overlap_below_half\tFP_zero_scored_overlap\tFP_exact_duplicates\tzero_reference_overlap_FPs_overlapping_original_TP\tqualification");
			raw(buffer.clear(), family).tab().append(seenInputs.size()).tab().append(trueCalls).tab().append(falseCalls).tab();
			buffer.append(maskedCalls).tab().append(recoveredReferences).tab().append(missedReferences).tab().append(maskedReferences).tab();
			buffer.append(overlappingHalf).tab().append(partialOverlap).tab().append(zeroOverlap).tab().append(duplicateFalseCalls).tab();
			out.writer.print(buffer.append(overlapsMatchedCall).append("\tORIGINAL_GRADES_UNCHANGED_IDENTITY_UNAVAILABLE\n"));
		}
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Rows             ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Row{
		Row(final String seqid_, final int start_, final int stop_, final String strand_,
			final String outcome_, final String id_, final String model_, final boolean duplicate_){
			require(start_>0 && stop_>=start_ && (strand_.equals("+") || strand_.equals("-")),
				"Overlap arithmetic needs positive inclusive coordinates and an explicit strand");
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			outcome=outcome_;
			id=id_;
			model=model_;
			duplicate=duplicate_;
		}
		int length(){return stop-start+1;}
		static int compare(final Row a, final Row b){
			if(b==null){return -1;}
			int c=a.seqid.compareTo(b.seqid);
			if(c==0){c=a.strand.compareTo(b.strand);}
			if(c==0){c=Integer.compare(a.start, b.start);}
			if(c==0){c=Integer.compare(a.stop, b.stop);}
			if(c==0){c=a.id.compareTo(b.id);}
			return c;
		}
		final String seqid, strand, outcome, id, model;
		final int start, stop;
		final boolean duplicate;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, prefix, family, inputId;
	private byte[] familyBytes, inputBytes;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectIntMap<String> seenInputs=new ObjectIntMap<String>(String.class);
	private final ArrayList<Row> queries=new ArrayList<Row>(), references=new ArrayList<Row>();
	private long trueCalls, falseCalls, maskedCalls, recoveredReferences, missedReferences, maskedReferences;
	private long overlappingHalf, partialOverlap, zeroOverlap, duplicateFalseCalls, overlapsMatchedCall;
	private final long[] nearCounts=new long[DISTANCES.length];
	private static final int[] DISTANCES={0, 10, 25, 50, 100, 200, 500};
	private static final String CALL_HEADER="input_id\tfamily\tseqid\tstart\tstop\tstrand\tquery_line\tmodel\texact_duplicate\toverlap_class\tbest_scored_reference\toverlap_bases\tcall_overlap_fraction\toriginal_reference_outcome\tcalls_half_covered_by_witness_including_masked\tbest_mask_reference\tmask_overlap_bases\tnearest_original_TP_line\tintervening_bases\tnearest_TP_start\tnearest_TP_stop\tnearest_TP_model\tverifier_identity\tqualification";
}
