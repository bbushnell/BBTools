package test.tracker;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import barcode.Barcode;
import tracker.ReadStats;

/** Literal ownership/repeated-merge and report-output checks for JT008.
 * @author Brian Bushnell, Jean */
public final class BarcodeMergeProbe {

	public static void main(String[] args) throws IOException{
		String out=null;
		boolean verify=false;
		for(String arg : args){
			if(arg.equals("mode=verify")){verify=true;}
			else if(arg.equals("mode=observe")){verify=false;}
			else if(arg.startsWith("out=")){out=arg.substring(4);}
			else{throw new IllegalArgumentException("Unknown JT008 probe argument: "+arg);}
		}
		System.out.println("case\tnondestructive_expected_or_control\tobserved\tagrees");
		try{
			reset(true);
			ReadStats a=worker("ACGTAC", 2), b=worker("ACGTAC", 3);
			ReadStats first=ReadStats.mergeAll();
			show("collision_first_sum", "5", count(first, "ACGTAC"));
			show("sources_after_first", "2,3", count(a, "ACGTAC")+","+count(b, "ACGTAC"));
			show("multi_worker_result_alias", "false", ""+(first.barcodeMap.get("ACGTAC")==a.barcodeMap.get("ACGTAC")));
			ReadStats second=ReadStats.mergeAll();
			show("collision_second_sum", "5", count(second, "ACGTAC"));
			show("prior_result_after_second", "5", count(first, "ACGTAC"));
			show("collision_third_sum", "5", count(ReadStats.mergeAll(), "ACGTAC"));

			reset(true);
			a=worker("ACGTAC", 2);
			b=worker("TGCAAA", 3);
			first=ReadStats.mergeAll();
			second=ReadStats.mergeAll();
			show("disjoint_repeated_control", "2,3,2,3", count(first, "ACGTAC")+","+count(first, "TGCAAA")+","+
				count(second, "ACGTAC")+","+count(second, "TGCAAA"));
			first.barcodeMap.get("ACGTAC").increment(7);
			show("result_mutation_source_count", "2", count(a, "ACGTAC"));
			show("result_mutation_other_result", "2", count(second, "ACGTAC"));

			reset(true);
			a=worker("ACGTAC", 2);
			show("singleton_borrowed_identity_control", "true", ""+(ReadStats.mergeAll()==a));
			show("singleton_count_control", "2", count(ReadStats.mergeAll(), "ACGTAC"));
			reset(true);
			show("empty_registry_control", "true", ""+(ReadStats.mergeAll()==null));

			reset(false);
			new ReadStats();
			new ReadStats();
			show("barcode_disabled_control", "true", ""+(ReadStats.mergeAll().barcodeMap==null));

			reset(true);
			a=worker("ACGTAC", 2);
			b=worker("ACGTAC", 3);
			Barcode decorated=new Barcode("ACGTAC", 2, 0, 17);
			decorated.frequency=0.25f;
			a.barcodeMap.put("ACGTAC", decorated);
			Barcode merged=ReadStats.mergeAll().barcodeMap.get("ACGTAC");
			show("first_metadata_control", "0,17,0.25", merged.expected+","+merged.tile+","+merged.frequency);

			if(out!=null){
				reset(true);
				worker("ACGTAC", 2);
				worker("ACGTAC", 3);
				first=ReadStats.mergeAll();
				first.writeBarcodesToFile(out+"/first-report.tsv");
				assert(!first.errorState) : "Barcode report writer failed; output cannot validate merge totals.";
				ReadStats.BARCODE_STATS_FILE=out+"/write-all-report.tsv";
				final boolean error=ReadStats.writeAll();
				assert(!error) : "writeAll failed; output cannot validate repeated merge totals.";
				final String expected="#Reads\t5\n#Barcodes\t1\nACGTAC\t5\n";
				show("first_report_control", escaped(expected), report(out+"/first-report.tsv"));
				show("write_all_after_merge_report", escaped(expected), report(out+"/write-all-report.tsv"));
			}
		}finally{ReadStats.clear();}
		if(verify && failures>0){throw new AssertionError("JT008 ownership/report mismatches: "+failures);}
	}

	private static String report(String path) throws IOException{
		return escaped(new String(Files.readAllBytes(Paths.get(path)), StandardCharsets.US_ASCII));
	}
	private static String escaped(String text){return text.replace("\t", "\\t").replace("\n", "\\n");}

	/** New workers register through the actual ReadStats constructor; fixture counts are literal. */
	private static ReadStats worker(String key, long count){
		assert(count>=0 && ReadStats.COLLECT_BARCODE_STATS) : "This fixture requires nonnegative counts and barcode collection.";
		ReadStats stats=new ReadStats();
		stats.barcodeMap.put(key, new Barcode(key, count));
		return stats;
	}

	private static String count(ReadStats stats, String key){return ""+stats.barcodeMap.get(key).count();}
	private static void reset(boolean barcodes){ReadStats.clear(); ReadStats.COLLECT_BARCODE_STATS=barcodes;}

	/** Agreed literal expectations preserve existing fast paths and first-entry metadata. */
	private static void show(String name, String expected, String observed){
		assert(name!=null && expected!=null && observed!=null) : "Every observation needs a named literal comparison.";
		System.out.println(name+"\t"+expected+"\t"+observed+"\t"+expected.equals(observed));
		if(!expected.equals(observed)){failures++;}
	}
	private static int failures=0;
}
