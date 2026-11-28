package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.ReportNcrnaLossContext;
import shared.Shared;
import structures.ByteBuilder;

/** Finite, external-data-free loss-context checks. No caller or alignment.
 * @author Diona
 */
public final class ReportNcrnaLossContextTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaLossContextTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Use native out=/JVM/thread options");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Fixtures require a fresh directory");
			final String root=new File(parser.out1).getAbsolutePath();
			check(new File(root).mkdir(), "Preserve prior fixture output");
			headerIdentity();
			pipeline(root);
			pp.outstream.println("NCRNA_LOSS_CONTEXT_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void headerIdentity(){
		final ByteBuilder b=new ByteBuilder();
		ReportNcrnaLossContext.decodeAttribute("ID=mid+%257C;input_id=in%2B1", "ID", b);
		check(b.toString().equals("mid+%7C"), "Literal plus survives and escaped percent is decoded exactly once");
		b.clear();
		ReportNcrnaLossContext.decodeAttribute("ID=x;input_id=in%2b1", "input_id", b);
		check(b.toString().equals("in+1"), "Both hex cases bind the same identity");
		fails(()->ReportNcrnaLossContext.decodeAttribute("ID=x%GG;input_id=i", "ID", new ByteBuilder()), "Malformed escapes fail loudly");
		fails(()->ReportNcrnaLossContext.decodeAttribute("ID=x%;input_id=i", "ID", new ByteBuilder()), "Truncated escapes fail loudly");
		check(ReportNcrnaLossContext.gap(101, 200, 250, 349)==49
			&& ReportNcrnaLossContext.gap(250, 349, 101, 200)==49
			&& ReportNcrnaLossContext.gap(1, 100, 101, 200)==0
			&& ReportNcrnaLossContext.gap(1, 100, 80, 150)==0,
			"Inclusive gap distinguishes no intervening bases from coordinate differences");
	}

	private static void pipeline(final String root){
		final String before=root+"/before", after=root+"/after", comparison=root+"/comparison";
		final String seeds=root+"/seeds.fa", fasta=root+"/scored.fa", lengths=root+"/contigs.tsv", manifest=root+"/cases.tsv";
		write(seeds, ">seed\n"+SEED+"\n");
		write(before+".family.tsv", "#NCRNA_NORMALIZED_FAMILY_V1\nfamily\tFMN\n");
		write(after+".family.tsv", "#NCRNA_NORMALIZED_FAMILY_V1\nfamily\tFMN\n");
		final StringBuilder refs=new StringBuilder(ReportNcrnaLossContext.REFERENCE_HEADER+"\n");
		final StringBuilder losses=new StringBuilder(SAVED_LOSS_HEADER+"\n");
		final StringBuilder queries=new StringBuilder(ReportNcrnaLossContext.QUERY_HEADER+"\n");
		for(int i=0; i<IDS.length; i++){
			refs.append(reference(IDS[i], CONTIGS[i], STARTS[i], STOPS[i], STRANDS[i]));
			losses.append(loss(i));
			queries.append(query(CONTIGS[i], STARTS[i], STOPS[i], STRANDS[i], "TP", IDS[i]));
		}
		refs.append(reference("neighbor", "tandem", 250, 349, '+'));
		for(int i=7; i<101; i++){refs.append(reference("dummy"+i, "dummy"+i, 101, 200, '+'));}
		queries.append(query("tandem", 250, 349, '+', "TP", "neighbor"));
		final String afterQueries=ReportNcrnaLossContext.QUERY_HEADER+"\n"
			+query("tandem", 250, 349, '+', "TP", "neighbor")
			+query("tandem", 181, 210, '+', "UNMATCHED_OTHER", "NA");
		write(before+".references.tsv", refs.toString());
		write(before+".queries.tsv", queries.toString());
		write(after+".queries.tsv", afterQueries);
		write(comparison+".loci.tsv", losses.toString());
		write(lengths, ReportNcrnaLossContext.LENGTH_HEADER+"\n"
			+"in+1\tedge\t20\nin+1\ttandem\t1000\nin+1\treverse\t1000\nin+1\tmiddle\t1000\nin+1\tlong\t1000\n");
		final String[] sequences={repeat('N', 20), SEED+repeat('N', 83), repeat('N', 83)+SEED,
			repeat('N', 40)+SEED+repeat('N', 43), repeat('N', 40)+SEED+repeat('N', 43), SEED+repeat('N', 266)+SEED};
		final StringBuilder fa=new StringBuilder();
		for(int i=0; i<IDS.length; i++){fa.append(fastaRecord(i, sequences[i]));}
		write(fasta, fa.toString());
		write(manifest, manifest(before, after, comparison, seeds, 6));
		run(manifest, fasta, lengths, root+"/valid", 6);
		checkLoci(root+"/valid.loci.tsv");
		checkSummary(root+"/valid.summary.tsv");

		// No supplied contig lengths must never silently turn unknown edges into NONE.
		run(manifest, fasta, null, root+"/no_lengths", 6);
		final Table noLengths=read(root+"/no_lengths.loci.tsv");
		final String[] middle=noLengths.find("reference_id", "mid+%7C");
		check(noLengths.get(middle, "at_contig_edge").equals("NA")
			&& noLengths.get(middle, "none_of_these").equals("false")
			&& noLengths.get(middle, "unassessed_without_positive_class").equals("true"),
			"Unknown right-boundary evidence remains unavailable, not a negative class");

		write(root+"/duplicate.loci.tsv", losses.toString()+loss(1));
		write(root+"/duplicate.tsv", manifest(before, after, root+"/duplicate", seeds, 7));
		fails(()->run(root+"/duplicate.tsv", fasta, lengths, root+"/bad_duplicate", 6), "Duplicate lost references cannot inflate a cohort");
		write(root+"/missing.fa", fa.toString().replace(fastaRecord(2, sequences[2]), ""));
		fails(()->run(manifest, root+"/missing.fa", lengths, root+"/bad_missing", 5), "A missing target sequence is not reported as zero seed hits");
		fails(()->run(manifest, fasta, lengths, root+"/bad_count", 7), "The complete declared reference FASTA population is checked during its one read");
		write(root+"/bad_lengths.tsv", ReportNcrnaLossContext.LENGTH_HEADER+"\nin+1\tlong\t399\n");
		fails(()->run(manifest, fasta, root+"/bad_lengths.tsv", root+"/bad_span", 6), "A feature stop cannot exceed its explicitly supplied contig length");

		write(root+"/empty.loci.tsv", SAVED_LOSS_HEADER+"\n");
		write(root+"/empty.tsv", manifest(before, after, root+"/empty", seeds, 0));
		run(root+"/empty.tsv", fasta, lengths, root+"/zero", 6);
		final Table zero=read(root+"/zero.summary.tsv");
		for(final String[] r:zero.rows){check(zero.get(r, "unique_lost_references").equals("0") && zero.get(r, "fraction_all_losses").equals("NA"), "Zero losses remain an empty denominator, not a measured zero fraction");}
	}

	private static void checkLoci(final String path){
		final Table t=read(path);
		check(t.rows.size()==6, "Every unique lost reference remains in the detailed report");
		for(final String[] r:t.rows){check(t.get(r, "scored_length_population").equals("101") && t.get(r, "length_p1").equals("100") && t.get(r, "length_p99").equals("100"), "Length tails use the complete scored cohort, not the lost subset");}
		final String[] shortRow=t.find("reference_id", "short");
		check(t.get(shortRow, "short_below_p1").equals("true") && t.get(shortRow, "at_contig_edge").equals("true") && t.get(shortRow, "no_forward_seed").equals("true"), "One lost reference may have several independently observed classes");
		final String[] left=t.find("reference_id", "left");
		check(t.get(left, "tandem_like_scored_reference").equals("true") && t.get(left, "near_reference_count").equals("1") && t.get(left, "nearest_reference_gap").equals("49"), "Only the other reference supports the tandem-like context");
		check(t.get(left, "near_call_before_count").equals("1") && t.get(left, "near_call_after_count").equals("1") && t.get(left, "overlapping_other_calls_after").equals("1"), "Self and overlapping same-locus calls do not become nearby distinct copies");
		check(t.get(left, "keys_5prime_only").equals("true") && t.get(left, "seed_min_center0_sense").equals("8"), "First-terminal seed uses its full forward occurrence center");
		final String[] right=t.find("reference_id", "right");
		check(t.get(right, "strand").equals("-") && t.get(right, "keys_3prime_only").equals("true") && t.get(right, "seed_min_center0_sense").equals("91"), "Saved minus-strand gene-sense sequence must not be reverse-complemented again");
		check(t.get(t.find("reference_id", "mid+%7C"), "none_of_these").equals("true"), "Complete negative evidence is NONE and escaped identity remains exact");
		check(t.get(t.find("reference_id", "unknown"), "unassessed_without_positive_class").equals("true"), "Unavailable contig length is not NONE");
		final String[] longRow=t.find("reference_id", "long");
		check(t.get(longRow, "long_above_p99").equals("true") && t.get(longRow, "seed_occurrences").equals("2") && t.get(longRow, "distinct_forward_seeds").equals("1")
			&& t.get(longRow, "keys_5prime_only").equals("false") && t.get(longRow, "keys_3prime_only").equals("false"), "Repeated copies of one key at opposite ends are not an end-only key population");
	}
	private static void checkSummary(final String path){
		final Table t=read(path);
		check(t.rows.size()==12, "Every declared descriptive class retains a row");
		for(final String[] r:t.rows){check(t.get(r, "unique_lost_references").equals("6") && t.get(r, "lost_with_class").equals("1"), "Overlapping class counts keep the same unique-loss denominator");}
		check(t.get(t.find("class", "AT_CONTIG_EDGE"), "assessable_losses").equals("5"), "Unavailable edges are excluded from assessed counts");
	}

	private static void run(final String manifest, final String fasta, final String lengths, final String output, final int records){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+manifest, "referencefa="+fasta, "orientation=sense", "expectedrecords="+records, "out="+output, "endpercent=25", "t=1"));
		if(lengths!=null){args.add("contiglengths="+lengths);}
		ReportNcrnaLossContext.runForTesting(args.toArray(new String[args.size()]));
	}
	private static String manifest(final String before, final String after, final String comparison, final String seeds, final int expected){return ReportNcrnaLossContext.MANIFEST_HEADER+"\nFMN\tunion_p99_plus10\t"+before+'\t'+after+'\t'+comparison+'\t'+seeds+'\t'+expected+'\n';}
	private static String reference(final String id, final String contig, final int start, final int stop, final char strand){return "in+1\t"+id+'\t'+contig+'\t'+start+'\t'+stop+'\t'+strand+"\ttrue\n";}
	private static String query(final String contig, final int start, final int stop, final char strand, final String outcome, final String ref){return "in+1\t"+contig+'\t'+start+'\t'+stop+'\t'+strand+"\tFMN\tfixture_model\t"+outcome+"\ttrue\t"+ref+'\n';}
	private static String loss(final int i){return "in+1\t"+IDS[i]+'\t'+CONTIGS[i]+'\t'+STARTS[i]+'\t'+STOPS[i]+'\t'+STRANDS[i]+"\tLOST_GRADED_REFERENCE\t"+(i==1 ? "OVERLAPPING_FAMILY_CALL_PRESENT" : "NO_OVERLAPPING_FAMILY_CALL")+"\tFMN\tfixture_model\tNA\tNA\tNA\n";}
	private static String fastaRecord(final int i, final String bases){return ">ID="+(i==3 ? "mid+%257C" : IDS[i])+";input_id=in%2B1;reference_family=FMN;reference_status=scored\n"+bases+'\n';}
	private static String repeat(final char c, final int count){
		final char[] out=new char[count];
		Arrays.fill(out, c);
		return new String(out);
	}
	private static void write(final String path, final String text){try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private static void fails(final Runnable action, final String message){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, message);
	}
	private static void check(final boolean pass, final String message){if(!pass){throw new AssertionError(message);}}
	private static Table read(final String path){
		try(Input input=new Input(path)){
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			final ArrayList<String> columns=new ArrayList<String>();
			for(int i=0; i<p.terms(); i++){columns.add(p.parseString(i));}
			final Table result=new Table(columns);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				check(p.terms()==columns.size(), "Detailed report schema and actual row width must agree");
				final String[] fields=new String[p.terms()];
				for(int i=0; i<fields.length; i++){fields[i]=p.parseString(i);}
				result.rows.add(fields);
			}
			return result;
		}
	}
	private static final class Table{
		Table(final ArrayList<String> columns_){columns=columns_;}
		String get(final String[] row, final String field){
			final int i=columns.indexOf(field);
			check(i>=0, "Expected named output column: "+field);
			return row[i];
		}
		String[] find(final String field, final String value){
			for(final String[] row:rows){if(get(row, field).equals(value)){return row;}}
			throw new AssertionError("Missing fixture row "+value);
		}
		final ArrayList<String> columns;
		final ArrayList<String[]> rows=new ArrayList<String[]>();
	}
	private static final String SEED="ACGTACGTACGTACGTA";
	/** Literal generic comparator contract, independent of the diagnostic parser constant. */
	private static final String SAVED_LOSS_HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_family_union\tbaseline_model\tbaseline_model_in_subset\tcandidate_family_union\tcandidate_model";
	private static final String[] IDS={"short","left","right","mid+%7C","unknown","long"};
	private static final String[] CONTIGS={"edge","tandem","reverse","middle","unknown","long"};
	private static final int[] STARTS={1,101,501,101,101,101}, STOPS={20,200,600,200,200,400};
	private static final char[] STRANDS={'+','+','-','+','+','+'};
}
