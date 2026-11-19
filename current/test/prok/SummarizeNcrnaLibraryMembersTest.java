package test.prok;

import java.io.File;
import java.util.Locale;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaRosterBaseline;
import prok.SummarizeNcrnaLibraryMembers;
import shared.Shared;
import structures.ByteBuilder;

/** Literal purity boundary, member/model weighting and original-gene duplicates.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaLibraryMembersTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaLibraryMembersTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown member-summary fixture argument: "+arg);
			}
			check(parser.out1!=null && new File(parser.out1).mkdir(), "Use a fresh fixture output directory");
			run(parser.out1);
			pp.outstream.println("NCRNA_LIBRARY_MEMBER_SHAPES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final String members=root+"/members.fa", consensus=root+"/consensus.fa", assignments=root+"/assignments.tsv";
		final String originals=root+"/originals.tsv", prefix=root+"/report";
		final ByteBuilder bases=new ByteBuilder();
		try(Output sequences=new Output(members); Output models=new Output(consensus);
			Output table=new Output(assignments); Output genes=new Output(originals)){
			models.writer.println(">model0\n"+repeat(bases, 100)+"\n>model1\n"+repeat(bases, 200)+"\n>model2\n"+repeat(bases, 50));
			table.writer.println("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus");
			genes.writer.println(PrepareNcrnaRosterBaseline.ORIGINAL_HEADER);
			for(int id=0; id<23; id++){
				final boolean full=id==10 || id==20;
				final int model=id<11 ? 0 : id<21 ? 1 : 2;
				final int length=full ? 100 : 50;
				final String name="RF00174_X_"+id+" cmfrom="+(full ? 1 : 45)+" cmto="+(full ? 85 : 65)
					+" cmlen=85 incomplete_model_coverage=false";
				sequences.writer.println(">"+name+"\n"+repeat(bases, length));
				table.writer.println(name+"\t-1\tRNA\tgroup\tmodel"+model+"\tkept");
				original(genes, id, id, full);
			}
			original(genes, 23, 0, false);// One duplicate original must not gain a builder vote.
			original(genes, 24, -1, false);// An excluded original still remains in its own denominator.
		}
		final String[] flags={"in="+members, "originals="+originals, "consensus="+consensus,
			"assignments="+assignments, "labels=literal", "out="+prefix, "expectedmembers=23", "expectedoriginals=25"};
		SummarizeNcrnaLibraryMembers.main(flags);
		final ObjectMap<String, Double> values=new ObjectMap<String, Double>(String.class, Double.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(prefix+".literal.summary.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==2 && !values.contains(p.parseString(0)), "Each summary metric must occur once");
				values.put(p.parseString(0), p.parseDouble(1));
			}
		}
		equal(values, "models", 3);
		equal(values, "construction_members", 23);
		equal(values, "both_ends_short_members", 21);
		// A permissive saved tolerance can label these geometric short spans complete.
		equal(values, "incomplete_coverage_members", 0);
		equal(values, "exact_full_span_members", 2);
		equal(values, "models_with_short_members", 3);
		equal(values, "short_only_models", 1);
		equal(values, "models_over_90_percent_short", 2);// 10/11 and2/2 pass;9/10 must not.
		equal(values, "short_models_with_full_span_members", 2);
		equal(values, "short_members_mixed_with_full_span", 19);
		equal(values, "short_members_mixed_with_any_other_class", 19);
		equal(values, "original_both_ends_short_genes", 23);
		equal(values, "original_short_genes_outside_construction", 1);
		equal(values, "original_short_genes_in_models_with_full_span_members", 20);
		equal(values, "median_across_models_of_member_length_ratio_5th", .5);
		equal(values, "median_across_models_of_member_length_ratio_50th", .5);
		equal(values, "median_across_models_of_member_length_ratio_95th", 1);
		equal(values, "5th_percentile_across_model_median_length_ratios", .25);
		equal(values, "95th_percentile_across_model_median_length_ratios", 1);
		final String[] wrong=flags.clone();
		wrong[5]="out="+root+"/wrong";
		wrong[6]="expectedmembers=24";
		boolean rejected=false;
		try{SummarizeNcrnaLibraryMembers.main(wrong);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An incomplete assignment population must not produce a silently smaller summary");
	}

	private static String repeat(final ByteBuilder bases, final int count){
		bases.clear();
		for(int i=0; i<count; i++){bases.append('A');}
		return bases.toString();
	}

	private static void original(final Output out, final int id, final int representative, final boolean full){
		final String[] fields=new String[29];
		java.util.Arrays.fill(fields, "0");
		fields[0]=""+id;
		fields[4]="1";
		fields[5]=full ? "100" : "50";
		fields[7]="true";
		fields[13]=""+id;
		fields[14]=""+representative;
		fields[18]=full ? "1" : "45";
		fields[19]=full ? "85" : "65";
		fields[20]="85";
		fields[17]=representative<0 ? "1" : "0";
		fields[21]="false";
		fields[22]=representative<0 ? "EXCLUDED_NATIVE_TRUNCATION" : "ELIGIBLE_NATIVE_HIT";
		out.writer.println(String.join("\t", fields));
	}

	private static void equal(final ObjectMap<String, Double> values, final String key, final double expected){
		check(values.contains(key) && values.get(key)==expected, "Literal member/model/original weighting differs for "+key);
	}

	private static void check(final boolean condition, final String reason){
		if(!condition){throw new AssertionError(reason);}
	}
}
