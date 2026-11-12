package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;

import parse.LineParser1;
import prok.BuildNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Hand-calculated timing and byte-accounting fixtures for the saved-pilot reader.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceCostTest{
	NcrnaReferenceCostTest(final File root_){root=root_;}
	void run(){
		check(root.mkdir(), "Cost fixture root must be new");
		final File good=prepare("good", 64, 20, "0");
		invoke(good);
		final File models=new File(good, "result.model_costs.tsv");
		row(models, "RF00230", new int[]{2, 3, 4, 5, 6, 7}, new double[]{64, 100, 20, 15, 200000000, .055556});
		row(models, "RF00498", new int[]{2, 3, 4, 5, 6, 7}, new double[]{256, 400, 40, 30, 100000000, .027778});
		final File summary=new File(good, "result.cost_summary.tsv");
		final String[] keys={"pilot_elapsed", "pilot_CPU", "preparation_and_extraction_CPU", "projected_all_new_models",
			"projected_without_two_viral_candidates", "projected_T_box", "final_sample_native_reports", "projected_native_report_bytes"};
		final double[] expected={120, 15, 6, .083333, .055556, .055556, 70, 400};
		for(int i=0; i<keys.length; i++){row(summary, keys[i], new int[]{1}, new double[]{expected[i]});}
		row(summary, "previous_report_model_rows", new int[]{1}, new double[]{2});
		reject(prepare("unknown_stage", 3, 20, "0"), "Unknown timing sample");
		reject(prepare("short", 64, 9, "0"), "short or incomplete");
		reject(prepare("failed", 64, 20, "1"), "must have exited successfully");
		reject(prepare("changed_projection", 64, 20, "0"), "numeric values disagree");
		expanded();
		System.err.println("NCRNA_REFERENCE_COST_TEST_PASS");
	}
	/** The same 1000-base population gives 200+100 projected CPU seconds and
	 * 300+100 output bytes. Partial cases must never fill the other model costs.
	 */
	private void expanded(){
		final File good=prepareExpanded("expanded_good", false);
		invokeExpanded(good);
		final File summary=new File(good, "result.cost_summary.tsv");
		row(summary, "pilot_CPU", new int[]{1}, new double[]{9});
		row(summary, "measured_models", new int[]{1}, new double[]{2});
		row(summary, "unpriced_models", new int[]{1}, new double[]{0});
		row(summary, "projected_all_models", new int[]{1}, new double[]{.083333});
		row(summary, "projected_native_report_bytes", new int[]{1}, new double[]{400});
		final File partial=prepareExpanded("expanded_partial", true);
		invokeExpanded(partial);
		final File partialSummary=new File(partial, "result.cost_summary.tsv");
		row(partialSummary, "measured_models", new int[]{1}, new double[]{1});
		row(partialSummary, "unpriced_models", new int[]{1}, new double[]{3});
		row(partialSummary, "projected_measured_subset", new int[]{1}, new double[]{.055556});
		term(partialSummary, "pilot_CPU", 1, "NA");
		term(partialSummary, "projected_measured_subset", 3, "MEASURED_SUBSET_ONLY_NOT_WHOLE_ROSTER_COST");
		final File partialModels=new File(partial, "result.model_costs.tsv");
		term(partialModels, "RF00014", 9, "NA");
		term(partialModels, "RF00014", 11, "UNPRICED_TOO_SHORT");
		term(partialModels, "RF00018", 7, "5");
		term(partialModels, "RF00018", 8, "NA");
		term(partialModels, "RF00018", 12, "PARTIAL_NATIVE_OUTPUTS");
		term(partialModels, "RF00019", 7, "NA");
		term(partialModels, "RF00019", 12, "NOT_STARTED");
		reject(prepareExpanded("expanded_short_measured", false), "qualification contradicts", true);
		reject(prepareExpanded("expanded_missing_bytes", false), "both measured output sizes", true);
		reject(prepareExpanded("expanded_false_complete", true), "cannot be labeled complete", true);
		reject(prepareExpanded("expanded_missing_row", true), "lacks its selected timing row", true);
	}
	private File prepareExpanded(final String name, final boolean partial){
		final File dir=new File(root, name);
		check(dir.mkdir(), "Expanded cost case must be new");
		write(dir, "exit.txt", partial ? "1\n" : "0\n");
		write(dir, "report.exit", "0\n");
		write(dir, "result.txt", "BROAD_TIMING_TASKS_FINISHED_CHECK_UNPRICED_ROWS\n");
		write(dir, "start.txt", "2026-11-08T00:00:00Z\n");
		write(dir, "end.txt", "2026-11-08T00:02:00Z\n");
		write(dir, "census.tsv", "genomes\tinput_file_bytes\tnucleotide_bases\tfasta_records\tboth_strands_search_space_Mb\n1000\t1000\t1000\t1000\t0.002\n");
		write(dir, "population.tsv", "stage\tgenomes\tnucleotide_bases\tfasta_records\n64\t64\t100\t64\n256\t256\t400\t256\n");
		final String assetSuffix="\t15.1\t10\t20\t00000000000000000000\tfixture.cm\n";
		write(dir, "model_assets.tsv", "accession\trelease\tCLEN\tGA_bits\tsha80\tpath\nRF00230"+assetSuffix+"RF00014"+assetSuffix
			+(partial ? "RF00018"+assetSuffix+"RF00019"+assetSuffix : ""));
		final String firstCpu=name.equals("expanded_short_measured") ? "9" : "20";
		final String hits=name.equals("expanded_missing_bytes") ? "NA" : "10";
		final String partialStatus=name.equals("expanded_false_complete") ? "COMPLETE_NATIVE_OUTPUTS" : "PARTIAL_NATIVE_OUTPUTS";
		write(dir, "costs.tsv", "model\tstage\tuser_CPU_seconds\tsystem_CPU_seconds\twall_seconds\tmax_RSS_KiB\ttiming_status\tresult_directory\thits_bytes\tsearch_report_bytes\toutput_status\n"
			+"RF00230\t64\t"+firstCpu+"\t0\t15\t10\tMEASURED_REUSED\tprior/startup\t"+hits+"\t20\tCOMPLETE_NATIVE_OUTPUTS\n"
			+(partial ? "RF00014\t256\t9\t0\t8\t20\tUNPRICED_TOO_SHORT\tsearch_256/RF00014\t20\t20\tCOMPLETE_NATIVE_OUTPUTS\n"
				+"RF00018\t64\tNA\tNA\tNA\tNA\tUNPRICED_NOT_COMPLETED\tsearch_64/RF00018\t5\tNA\t"+partialStatus+"\n"
				+(name.equals("expanded_missing_row") ? "" : "RF00019\t64\tNA\tNA\tNA\tNA\tUNPRICED_NOT_COMPLETED\tsearch_64/RF00019\tNA\tNA\tNOT_STARTED\n")
				: "RF00014\t256\t39\t1\t30\t20\tMEASURED\tsearch_256/RF00014\t20\t20\tCOMPLETE_NATIVE_OUTPUTS\n"));
		for(final String phase:new String[]{"prepare", "startup", "searches"}){write(dir, phase+".time", "1\t2\t3\t10\n");}
		return dir;
	}
	private static void invokeExpanded(final File dir){
		BuildNcrnaReference.main(new String[]{"mode=cost", "in="+dir.getPath(),
			"census="+new File(dir, "census.tsv").getPath(), "out="+new File(dir, "result").getPath(), "t=1"});
	}
	private File prepare(final String name, final int stage, final int cpu, final String pointExit){
		final File dir=new File(root, name);
		check(dir.mkdir(), "Cost case must be new");
		write(dir, "exit.txt", "0\n");
		write(dir, "result.txt", "ADDITIONAL_MODEL_COST_PILOT_COMPLETE\n");
		write(dir, "start.txt", "2026-11-08T00:00:00Z\n");
		write(dir, "end.txt", "2026-11-08T00:02:00Z\n");
		write(dir, "census.tsv", "genomes\tinput_file_bytes\tnucleotide_bases\tfasta_records\tboth_strands_search_space_Mb\n1000\t1000\t1000\t1000\t0.002\n");
		write(dir, "population.tsv", "stage\tgenomes\tnucleotide_bases\tfasta_records\n64\t64\t100\t64\n256\t256\t400\t256\n");
		write(dir, "candidates.tsv", "RF00230\tT_box\tfixture\t20\tleader\t\t10\nRF00498\tviral\tfixture\t20\tleader\t\t10\n");
		write(dir, "costs.tsv", "model\tstage\tuser_CPU_seconds\tsystem_CPU_seconds\twall_seconds\tmax_RSS_KiB\ttiming_status\n"
			+"RF00230\t"+stage+'\t'+cpu+"\t0\t15\t10\tMEASURED\nRF00498\t256\t39\t1\t30\t10\tMEASURED\n");
		write(dir, "previous.tsv", "accession\tname\tsample_genomes\tsample_bases\tCPU_seconds\twall_seconds\tCPU_seconds_per_Gbp\tprojected_full_CPU_hours\tstatus\n"
			+"RF00230\tT_box\t64\t100\t20.00\t15.00\t200000000.000000\t"
			+(name.equals("changed_projection") ? "0.065556" : "0.055556")+"\tLINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN\n"
			+"RF00498\tviral\t256\t400\t40.00\t30.00\t100000000.000000\t0.027778\tLINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN\n");
		for(final String phase:new String[]{"prepare", "extract", "startup", "models64", "models256"}){
			write(dir, phase+".time", "1\t2\t3\t10\n");
		}
		check(new File(dir, "search_64/RF00230").mkdirs() && new File(dir, "search_256/RF00498").mkdirs(), "Native point fixture directories must be new");
		write(dir, "search_64/RF00230/exit.txt", pointExit+'\n');
		write(dir, "search_256/RF00498/exit.txt", "0\n");
		write(dir, "search_64/RF00230/hits.tbl", "0123456789");
		write(dir, "search_64/RF00230/search.out", "01234567890123456789");
		write(dir, "search_256/RF00498/hits.tbl", "01234567890123456789");
		write(dir, "search_256/RF00498/search.out", "01234567890123456789");
		return dir;
	}
	private static void invoke(final File dir){
		BuildNcrnaReference.main(new String[]{"mode=cost", "in="+dir.getPath(),
			"census="+new File(dir, "census.tsv").getPath(), "compare="+new File(dir, "previous.tsv").getPath(),
			"out="+new File(dir, "result").getPath(), "t=1"});
	}
	private static void reject(final File dir, final String fragment){
		reject(dir, fragment, false);
	}
	private static void reject(final File dir, final String fragment, final boolean expanded){
		try{
			if(expanded){invokeExpanded(dir);}else{invoke(dir);}
		}catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Wrong cost failure: "+failure);
			check(!new File(dir, "result.cost_summary.tsv").exists(), "Failed input cannot publish a completed cost summary");
			return;
		}
		throw new AssertionError("Expected cost failure: "+fragment);
	}
	/** Literal field expectations distinguish NA, zero, partial and complete values. */
	private static void term(final File file, final String key, final int column, final String expected){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(file.getPath())){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(key, 0)){continue;}
				found++;
				check(p.termEquals(expected, column), "Expected "+expected+" for "+key+" column "+column);
			}
		}
		check(found==1, "Expected one literal-field row for "+key);
	}
	private static void write(final File dir, final String name, final String text){
		try(Output out=new Output(new File(dir, name).getPath())){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	/** Compare numeric values, not the native formatter's optional decimal places. */
	private static void row(final File file, final String key, final int[] columns, final double[] expected){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(file.getPath())){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(key, 0)){continue;}
				found++;
				for(int i=0; i<columns.length; i++){
					final double actual=p.parseDouble(columns[i]);
					check(Math.abs(actual-expected[i])<1e-9, "Hand-calculated value for "+key+" column "+columns[i]+": expected "+expected[i]+", got "+actual);
				}
			}
		}
		check(found==1, "Expected one output row for "+key);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
}
