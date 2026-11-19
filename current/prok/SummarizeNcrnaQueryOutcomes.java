package prok;

import java.util.ArrayList;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Attributes already matched query outcomes to emitted model prefixes.
 * Adapted from template.A_SampleByteFile, using native byte streams and a reused
 * parser/buffer. Reads the completed GradeNcrnaReference raw/summary reports and
 * their original explicit family map; no GFF, DNA or reference payload is read.
 * No matching, biological subtype assignment or score selection occurs here.
 *
 * TP, MASKED, unmatched duplicate and other unmatched are disjoint for included
 * calls. Duplicate/scored-eligibility flags remain independent counters; masked
 * duplicates remain masked and cannot become TP. Excluded rules are reconciled separately;
 * zero-query included prefixes survive. Family totals must agree with the original joint summary
 * before an output opens. Small rule objects survive; query fields never allocate.
 * Text fields retain the joint report's percent escaping. Each instance and its
 * native stream owners belong to one thread.
 * @author Raiden
 */
public final class SummarizeNcrnaQueryOutcomes{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, SummarizeNcrnaQueryOutcomes.class, false);
		try{
			new SummarizeNcrnaQueryOutcomes(pp.args).process();
			timer.stop();
			pp.outstream.println("NCRNA_QUERY_OUTCOMES_SUMMARIZED time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaQueryOutcomes(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("familymap")){familyMap=b;}
			else if(a.equals("family")){family=b;}
			else if(a.equals("expected")){
				final String[] values=b.split(",");
				expected=new long[values.length];
				for(int i=0; i<values.length; i++){
					expected[i]=Parse.parseKMG(values[i]);
					require(expected[i]>=0, "Expected included-query counts cannot be negative");
				}
			}
			else{require(parser.parse(arg, a, b), "Unknown query-summary option: "+arg);}
		}
		prefix=parser.in1;
		output=parser.out1;
		require(prefix!=null && output!=null && familyMap!=null && family!=null, "Supply in=completed_grade_prefix familymap=original_map family=scoring_family out=fresh_summary.tsv");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Summary requires complete reports and a fresh output");
		require(Tools.testInputFiles(false, true, prefix+".raw.tsv", prefix+".summary.tsv", familyMap)
			&& Tools.testOutputFiles(false, false, false, output), "Cannot read completed reports/map or create fresh output");
		require(Tools.testForDuplicateFiles(true, prefix+".raw.tsv", prefix+".summary.tsv", familyMap, output), "Summary must not overwrite an input");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads and validates all declarations and outcomes before creating the report. */
	private void process(){
		loadRules();
		readQueries();
		checkSummary();
		int selected=0;
		for(final Bucket bucket:buckets){
			if(!selected(bucket)){continue;}
			if(expected!=null){
				require(selected<expected.length && queryCount(bucket)==expected[selected], "Included-query count differs for model prefix "+bucket.rule.modelPrefix);
			}
			selected++;
		}
		require(selected>0 && (expected==null || selected==expected.length), "Selected included family and expected prefix counts must agree with the original map");
		try(Output out=new Output(output)){
			out.writer.println(HEADER);
			for(final Bucket bucket:buckets){
				if(!selected(bucket)){continue;}
				final Rule rule=bucket.rule;
				buffer.clear();
				BuildNcrnaReference.field(buffer, rule.family).tab();
				BuildNcrnaReference.field(buffer, rule.modelPrefix).tab();
				BuildNcrnaReference.field(buffer, rule.displayLabel).tab().append(rule.include).tab().append(queryCount(bucket));
				for(final long count:bucket.counts){buffer.tab().append(count);}
				out.writer.print(buffer.append("\tEMITTED_MODEL_ATTRIBUTION_AFTER_JOINT_MATCH_NOT_SUBTYPE_ACCURACY\n"));
			}
		}
	}

	/** Loads the original projection, retaining excluded and zero-query rules. */
	private void loadRules(){
		for(final Rule rule:GradeNcrnaReference.loadRules(familyMap)){
			for(final Bucket previous:buckets){
				require(!previous.rule.modelPrefix.equals(rule.modelPrefix), "Repeated model prefix: "+rule.modelPrefix);
			}
			buckets.add(new Bucket(rule));
		}
		require(!buckets.isEmpty(), "An outcome summary requires the original nonempty family map");
	}

	/** Consumes the recorded joint outcome; flags never replace the recorded match. */
	private void readQueries(){
		final String path=prefix+".raw.tsv";
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.RAW_HEADER, path);
			long line=1;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				line++;
				p.set(row);
				require(p.terms()==18, "Joint raw report requires eighteen fields");
				if(p.termEquals("REFERENCE", 1)){continue;}
				require(p.termEquals("QUERY", 1) && p.parseLong(9)>0, "Expected a query with its original positive GFF line number");
				Bucket selected=null;
				for(final Bucket bucket:buckets){
					if(!matches(bucket.prefixToken, 15, true)){continue;}
					if(selected!=null){throw new IllegalArgumentException("Query model matches multiple prefixes at report line "+line+": "+p.parseString(15));}
					selected=bucket;
				}
				if(selected==null || !matches(selected.familyToken, 2, false)){
					throw new IllegalArgumentException("Query model/family is absent from the original projection at report line "+line+": "+p.parseString(15));
				}
				final long[] counts=selected.counts;
				if(!selected.rule.include){
					require(p.termEquals("EXCLUDED", 7) && p.termEquals(".", 13) && p.termEquals(".", 14), "Excluded query must retain the grader's unscored status");
					counts[4]++;
					continue;
				}
				final boolean duplicate=flag(13), eligible=flag(14);
				if(duplicate){counts[5]++;}
				if(eligible){counts[6]++;}
				if(p.termEquals("TP", 7)){
					require(eligible && !duplicate && !p.termEquals(".", 10) && p.termEquals(".", 12),
						"Recorded TP needs scored eligibility, a matched reference, no duplicate flag and no mask reference");
					counts[0]++;
				}else if(p.termEquals("MASKED", 7)){
					require(!eligible && p.termEquals(".", 10) && !p.termEquals(".", 12), "Scored-first matching cannot mask a scored-eligible or paired query");
					counts[1]++;
				}else{
					require(p.termEquals("FP", 7) && p.termEquals(".", 10) && p.termEquals(".", 12), "Included unpaired/unmasked query requires FP outcome");
					counts[duplicate ? 2 : 3]++;
				}
			}
		}
	}

	/** Reconciles every prefix's counters with the completed joint family summary. */
	private void checkSummary(){
		final String path=prefix+".summary.tsv";
		long inputs=-1;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.SUMMARY_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==17 && p.parseLong(1)>0, "Completed grade summary requires seventeen fields and positive inputs");
				if(inputs<0){inputs=p.parseLong(1);}
				require(inputs==p.parseLong(1), "Joint summary families disagree on their input population");
				final long[] sum=new long[7];
				int found=0;
				for(final Bucket bucket:buckets){
					if(!matches(bucket.familyToken, 0, false)){continue;}
					require(!bucket.summarySeen, "Repeated joint summary family: "+bucket.rule.family);
					bucket.summarySeen=true;
					found++;
					for(int i=0; i<sum.length; i++){sum[i]+=bucket.counts[i];}
				}
				require(found>0 && sum[0]==p.parseLong(2) && sum[2]+sum[3]==p.parseLong(3)
					&& sum[1]==p.parseLong(7) && sum[5]==p.parseLong(8) && sum[4]==p.parseLong(9),
					"Emitted-prefix outcomes disagree with the original joint summary");
			}
		}
		for(final Bucket bucket:buckets){require(bucket.summarySeen, "Missing completed summary family: "+bucket.rule.family);}
	}

	/** Only the exact boolean spelling emitted by the grader is accepted. */
	private boolean flag(final int column){
		if(!p.termEquals("true", column) && !p.termEquals("false", column)){
			throw new IllegalArgumentException("Expected a literal true/false flag at column "+column);
		}
		return p.termEquals("true", column);
	}
	/** Matches retained UTF-8/escaped tokens directly, including non-ASCII map text. */
	private boolean matches(final byte[] token, final int column, final boolean prefix){
		assert(token!=null && column>=0 && column<p.terms()) : "Validated report schemas bound every projected token comparison";
		final int length=p.setBounds(column);
		if(prefix ? length<token.length : length!=token.length){return false;}
		final byte[] row=p.line();
		for(int i=0; i<token.length; i++){if(row[p.a()+i]!=token[i]){return false;}}
		if(!prefix || length==token.length){return true;}
		// NcrnaReferenceCalls.projectModel requires a WHOLE prefix. The grade
		// escapes a native tab as %09, while underscore and space stay literal.
		final int after=p.a()+token.length;
		return row[after]=='_' || row[after]==' ' || (length-token.length>=3
			&& row[after]=='%' && row[after+1]=='0' && row[after+2]=='9');
	}
	/** Excluded families are reconciled separately and never enter selected call totals. */
	private boolean selected(final Bucket bucket){return bucket.rule.include && bucket.rule.family.equals(family);}
	private static long queryCount(final Bucket bucket){
		final long[] c=bucket.counts;
		return c[0]+c[1]+c[2]+c[3];
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Bucket{
		Bucket(final Rule rule_){
			rule=rule_;
			familyToken=escaped(rule.family);
			prefixToken=escaped(rule.modelPrefix);
		}
		private static byte[] escaped(final String text){
			return BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();
		}
		final Rule rule;
		final byte[] familyToken, prefixToken;
		/** TP, masked, unmatched-duplicate, unmatched-other, excluded, duplicate, eligible. */
		final long[] counts=new long[7];
		boolean summarySeen;
	}
	private String prefix, output, familyMap, family;
	private long[] expected;
	private final ArrayList<Bucket> buckets=new ArrayList<Bucket>();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	public static final String HEADER="scoring_family\tmodel_prefix\tdisplay_label\tincluded\tincluded_queries\tTP\tMASKED\tunmatched_duplicate\tunmatched_other\texcluded\tduplicate_flag\tscored_eligible_flag\tqualification";
}
