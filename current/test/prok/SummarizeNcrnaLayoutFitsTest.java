package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.NcrnaDdlLayoutCounts;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaLayoutFits;

/** Literal storage/quality frontier and opposite-half rejection witnesses.
 * @author Ganyu
 */
final class SummarizeNcrnaLayoutFitsTest{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	static void run(final String root){
		final String first=root+"/summary_a.tsv";
		final String second=root+"/summary_b.tsv";
		final String layouts=root+"/summary_layouts.tsv";
		try(Output out=new Output(layouts)){
			out.writer.println(NcrnaDdlLayoutCounts.LAYOUT_HEADER_V2);
			out.writer.println(metadata("best", 1024));
			out.writer.println(metadata("one_percent", 512));
			out.writer.println(metadata("three_percent", 256));
			out.writer.println(metadata("dominated", 512));
			out.writer.println(metadata("smallest", 128));
			out.writer.println(metadata("failed_fit", 2048));
		}
		try(Output out=new Output(first)){
			out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(row("best", 1024, 1000, 1000, 0, 0));
			out.writer.println(row("one_percent", 512, 1010, 1010, 0, 0));
			out.writer.println(row("three_percent", 256, 1030, 1030, 0, 0));
		}
		try(Output out=new Output(second)){
			out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(row("dominated", 512, 1100, 1100, 0, 0));
			out.writer.println(row("smallest", 128, 1100, 1100, 0, 0));
			out.writer.println(row("failed_fit", 2048, 1000, 900, 120, -20));
		}
		final String merged=root+"/summary_merge";
		SummarizeNcrnaLayoutFits.main(new String[]{"in="+first+","+second, "layouts="+layouts, "out="+merged, "expectedlayouts=6"});
		final ArrayList<String> ids=new ArrayList<String>();
		try(Input in=new Input(merged+".pareto.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){ids.add(new String(line, StandardCharsets.UTF_8).split("\t", 2)[0]);}
		}
		check(ids.equals(Arrays.asList("best", "one_percent", "three_percent", "smallest")),
			"A cheap poor row may be nondominated, but a failed opposite-half fit cannot enter at its attempted lower cost");
		try(Input in=new Input(merged+".choices.tsv")){
			in.reader.nextLine();
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).startsWith("1\t1000\tone_percent\t512\t1010\t"),
				"The one-percent threshold is inclusive and selects the smallest eligible payload");
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).startsWith("3\t1000\tthree_percent\t256\t1030\t"),
				"The three-percent threshold is inclusive and cannot select the smaller out-of-budget row");
			check(in.reader.nextLine()==null, "Exactly the two requested size/quality thresholds are reported");
		}
		boolean failed=false;
		try{SummarizeNcrnaLayoutFits.merge(new String[]{first, first}, layouts, root+"/summary_duplicate", 6);}
		catch(RuntimeException | AssertionError expected){failed=true;}
		check(failed && !new File(root+"/summary_duplicate.comparison.tsv").exists(), "Duplicate result input cannot silently count a fit twice");
		final String malformed=root+"/summary_bad.tsv";
		final String[] bad=row("failed_fit", 2048, 1000, 900, 120, -20).split("\t", -1);
		bad[23]="900";
		try(Output out=new Output(malformed)){
			out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(String.join("\t", bad));
		}
		failed=false;
		try{SummarizeNcrnaLayoutFits.merge(new String[]{malformed}, layouts, root+"/summary_rejected", 6);}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(root+"/summary_rejected.comparison.tsv").exists(),
			"A row that substitutes an unapproved fitted cost must fail before publication");
		final String replacement=root+"/summary_wrong_membership.tsv";
		try(Output out=new Output(replacement)){
			out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(row("substituted_id", 1024, 1000, 1000, 0, 0));
			out.writer.println(row("one_percent", 512, 1010, 1010, 0, 0));
			out.writer.println(row("three_percent", 256, 1030, 1030, 0, 0));
		}
		reject(new String[]{replacement, second}, layouts, root+"/summary_wrong_ids", "undeclared layout ID");
		for(final int field:new int[]{2, 37}){
			final String changed=root+"/summary_changed_"+field+".tsv";
			final String[] altered=row("best", 1024, 1000, 1000, 0, 0).split("\t", -1);
			altered[field]=field==2 ? "8" : "RAW_KEYS";
			try(Output out=new Output(changed)){
				out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
				out.writer.println(String.join("\t", altered));
				out.writer.println(row("one_percent", 512, 1010, 1010, 0, 0));
				out.writer.println(row("three_percent", 256, 1030, 1030, 0, 0));
			}
			reject(new String[]{changed, second}, layouts, root+"/summary_wrong_shape_"+field,
				field==2 ? "canonical layout metadata" : "query form");
		}
		final String badCounts=root+"/summary_bad_top_k.tsv";
		final String[] changedCounts=row("best", 1024, 1000, 1000, 0, 0).split("\t", -1);
		changedCounts[39]="3";
		try(Output out=new Output(badCounts)){
			out.writer.println(SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(String.join("\t", changedCounts));
		}
		reject(new String[]{badCounts, second}, layouts, root+"/summary_top_k_rejected", "Top-k counts");
		final String historical=root+"/summary_missing_top_k.tsv";
		try(Output out=new Output(historical)){
			out.writer.println(String.join("\t", Arrays.copyOf(SummarizeNcrnaLayoutFits.HEADER.split("\t", -1), 38)));
			out.writer.println(String.join("\t", Arrays.copyOf(changedCounts, 38)));
		}
		failed=false;
		try{SummarizeNcrnaLayoutFits.merge(new String[]{historical}, layouts, root+"/summary_historical_rejected", 6);}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(root+"/summary_historical_rejected.comparison.tsv").exists(),
			"Historical reports without top-k and available-population fields must fail instead of inventing zero counts");
		matchedSeeds(root, layouts, merged+".comparison.tsv");
		System.err.println("NCRNA_DDL_LAYOUT_PARETO_TEST_PASS");
	}

	/** The first seed's best recipe loses on the paired mean;the paired gain
	 * range must retain its negative side rather than hiding it in the mean. */
	private static void matchedSeeds(final String root, final String layouts, final String first){
		final String roster=root+"/recipe_seed_roster.tsv";
		try(Output out=new Output(roster)){
			out.writer.println("ddl_seed\tcomparison_table");
			for(final int seed:new int[]{1,7}){
				final String dir=root+"/recipe_seed_"+seed;
				check(new File(dir).mkdir(), "Create only the named saved-table seed fixture");
				try(Output sidecar=new Output(dir+"/ddl_seed.txt")){sidecar.writer.println(""+seed);}
				final String path=dir+"/comparison.tsv";
				if(seed==1){
					try(Input in=new Input(first); Output copy=new Output(path)){
						for(byte[] text; (text=in.reader.nextLine())!=null;){copy.writer.println(text);}
					}
				}else{
					try(Output copy=new Output(path)){
						copy.writer.println(SummarizeNcrnaLayoutFits.HEADER);
						copy.writer.println(row("best", 1024, 1200, 1200, 0, 0));
						copy.writer.println(row("one_percent", 512, 1010, 1010, 0, 0));
						copy.writer.println(row("three_percent", 256, 1030, 1030, 0, 0));
						copy.writer.println(row("dominated", 512, 1100, 1100, 0, 0));
						copy.writer.println(row("smallest", 128, 1100, 1100, 0, 0));
						copy.writer.println(row("failed_fit", 2048, 1400, 1300, 120, -20));
					}
				}
				out.writer.println(seed+"\t"+new File(path).getAbsolutePath());
			}
		}
		final String result=root+"/matched_recipes";
		SummarizeNcrnaLayoutFits.main(new String[]{"seedroster="+roster, "layouts="+layouts,
			"reference=best", "out="+result, "expectedlayouts=6"});
		try(Input in=new Input(result+".seed_choices.tsv")){
			check(in.reader.nextLine()!=null, "Matched choices need an explicit schema");
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).startsWith("available_rank_sum\tone_percent\t2\t1010\t"),
				"The mean across both seeds must replace a lucky first-seed winner");
		}
		boolean found=false;
		final parse.LineParser1 p=new parse.LineParser1('\t');
		try(Input in=new Input(result+".seed_summary.tsv")){
			in.reader.nextLine();
			for(byte[] text; (text=in.reader.nextLine())!=null;){
				p.set(text);
				if(!p.termEquals("one_percent", 0) || !p.termEquals("available_rank_sum", 1)){continue;}
				check(p.terms()==13 && p.parseInt(2)==2 && p.parseDouble(3)==1010
					&& p.parseDouble(6)==1100 && p.parseDouble(7)==90
					&& p.parseLong(8)==-10 && p.parseLong(9)==190 && p.parseInt(10)==1 && p.parseInt(12)==1,
					"The recipe comparison must report paired improvement and its variation on the identical seeds");
				found=true;
			}
		}
		check(found, "The complete paired summary must retain the expected recipe/metric row");
		final String mixed=root+"/recipe_seed_7/mixed_population.tsv", mixedRoster=root+"/mixed_seed_roster.tsv";
		try(Input in=new Input(root+"/recipe_seed_7/comparison.tsv"); Output out=new Output(mixed)){
			out.writer.println(in.reader.nextLine());
			for(byte[] text; (text=in.reader.nextLine())!=null;){
				final String[] fields=new String(text, StandardCharsets.UTF_8).split("\t", -1);
				// Add one unavailable,non-low parent consistently to this whole seed.
				for(final int field:new int[]{18,19,23}){fields[field]=""+(Long.parseLong(fields[field])+384);}
				fields[32]="11";
				out.writer.println(String.join("\t", fields));
			}
		}
		try(Output out=new Output(mixedRoster)){
			out.writer.println("ddl_seed\tcomparison_table");
			out.writer.println("1\t"+new File(root+"/recipe_seed_1/comparison.tsv").getAbsolutePath());
			out.writer.println("7\t"+new File(mixed).getAbsolutePath());
		}
		boolean rejected=false;
		try{SummarizeNcrnaLayoutFits.mergeSeeds(mixedRoster, layouts, "best", root+"/mixed_seed_rejected", 6);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("identical frozen population"), "The mismatched-seed witness must reach its intended population guard");
			rejected=true;
		}
		check(rejected && !new File(root+"/mixed_seed_rejected.seed_summary.tsv").exists()
			&& !new File(root+"/mixed_seed_rejected.seed_choices.tsv").exists(),
			"Internally valid seeds with different judging populations cannot be averaged or partially published");
	}

	private static void reject(final String[] inputs, final String layouts, final String output, final String reason){
		boolean failed=false;
		try{SummarizeNcrnaLayoutFits.merge(inputs, layouts, output, 6);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "A malformed roster must fail for its intended mismatch: "+expected.getMessage());
			failed=true;
		}
		check(failed, "Equal row counts cannot conceal changed canonical layout membership or geometry");
		for(final String suffix:SummarizeNcrnaLayoutFits.SUFFIXES){
			check(!new File(output+suffix).exists(), "Rejected merges must not publish any result: "+output+suffix);
		}
	}

	private static String metadata(final String id, final int bytes){
		check(bytes>0 && bytes%2==0, "Synthetic plain-key payloads contain whole two-byte bins");
		return "2\t"+id+"\tRF_TEST\t7\t10\tMIDDLE_RIGHT\tMIXED\t"+(bytes/2)
			+"\t0\t0\tBOTH\t256\t0\t4\t12\tWRAP\tmember,consensus\t0\t512\t"+(512L*bytes)
			+"\t262144\t0\t0\t0\t0\t0\t0\t0\t0\tNATIVE_DDL\tQUERY_DDL";
	}

	private static String row(final String id, final int bytes, final int unit, final int fitted, final int gainA, final int gainB){
		check(unit-fitted==gainA+gainB, "Literal costs and gains must describe the same judging halves");
		final boolean pass=gainA>0 && gainB>0;
		return id+"\tRF_TEST\t7\t10\tMIDDLE_RIGHT\tMIXED\t"+(bytes/2)+"\t0\t0\tBOTH\t256\t0\tmember,consensus\t2\tNATIVE_DDL\tWRAP\t"
			+(512L*bytes)+"\t"+bytes+"\t"+unit+"\t"+fitted+"\t"+gainA+"\t"+gainB+"\t"+pass+"\t"+(pass ? fitted : unit)
			+"\t10\t11\t"+(pass ? 11 : 10)+"\t10\t11\t"+(pass ? 11 : 10)+"\ttrue\t512\t10\t2\t5000\t20261113\tSYNTHETIC_COHORT\tQUERY_DDL"
			+"\t4\t5\t2\t2\t10\t"+(pass ? fitted : unit);
	}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
}
