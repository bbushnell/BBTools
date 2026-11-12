package test.prok;

import java.io.File;

import parse.LineParser1;

import prok.CompareCmSearchCosts;
import prok.CompareCmSearchCosts.Counts;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Finite report comparisons with literal multiplicity and cost oracles.
 * Synthetic scores deliberately include opposite flags at the same printed GA.
 * No native search, biological admission or performance claim is made.
 * @author Brian Bushnell, Raiden
 */
public final class CompareCmSearchCostsTest{

	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){
			if(arg.startsWith("out=")){out=arg.substring(4);}
		}
		check(out!=null && new File(out).mkdirs(), "Use out=fresh_fixture_directory");
		run(out);
		System.err.println("COMPARE_CM_SEARCH_COSTS_TEST_PASS");
	}

	private static void run(final String out){
		final String control=out+"/control", pilot=out+"/pilot";
		final String point=pilot+"/search_2/RF00050";
		check(new File(control).mkdirs() && new File(point).mkdirs(), "Fresh point directories");
		final String old=hit("A", "+", "20.0", "!")+hit("A", "+", "20.0", "!")+hit("B", "+", "30.0", "!");
		final String next=hit("A", "+", "20.0", "!")+hit("A", "+", "20.0", "?")
			+hit("B", "-", "30.0", "!")+hit("C", "+", "40.0", "!");
		write(control+"/hits.tbl", old+"# [ok]\n");
		write(point+"/hits.tbl", next+"# [ok]\n");
		final Counts counts=CompareCmSearchCosts.compareHits(control+"/hits.tbl", point+"/hits.tbl", "RF00050", 20, null, null);
		check(counts.oldHits==3 && counts.newHits==4, "Raw rows must retain '?' as well as '!'");
		check(counts.oldIncluded==3 && counts.newIncluded==3 && counts.retained==1 && counts.lost==2 && counts.gained==2,
			"Duplicate multiplicity and opposite-strand coordinates must not collapse");
		check(counts.roundedDisagreements==1, "Printed-at-GA '?' must remain excluded");
		write(out+"/empty.tbl", "# [ok]\n");
		final Counts empty=CompareCmSearchCosts.compareHits(out+"/empty.tbl", out+"/empty.tbl", "RF00050", 20, null, null);
		check(empty.oldHits==0 && empty.newHits==0 && empty.retained==0, "Completed empty reports are valid");
		write(out+"/incomplete.tbl", old);
		expectFailure(new Runnable(){public void run(){
			CompareCmSearchCosts.compareHits(out+"/incomplete.tbl", out+"/empty.tbl", "RF00050", 20, null, null);
		}}, "Incomplete table must fail");

		final String oldCommand="cmsearch --cpu 1 --cut_ga --nohmmonly --noali --notextw -Z 0.004 --tformat fasta --tblout old.tbl -o old.out model.cm sample.fa";
		final String newCommand="cmsearch --cpu 1 --default --nohmmonly --noali --notextw -Z 0.004 -T 15 --incT 20 --tformat fasta --tblout new.tbl -o new.out model.cm sample.fa";
		check(CompareCmSearchCosts.checkCommand(oldCommand, "model.cm", 20, 15, false)
			.equals(CompareCmSearchCosts.checkCommand(newCommand, "model.cm", 20, 15, true)), "Declared threshold changes preserve other profile tokens");
		expectFailure(new Runnable(){public void run(){
			CompareCmSearchCosts.checkCommand(newCommand.replace("--incT 20", "--incT 19"), "model.cm", 20, 15, true);
		}}, "Wrong inclusion cutoff must fail");
		expectFailure(new Runnable(){public void run(){
			CompareCmSearchCosts.checkCommand(newCommand.replace("--incT 20 ", ""), "model.cm", 20, 15, true);
		}}, "Missing inclusion cutoff must fail");
		expectFailure(new Runnable(){public void run(){
			CompareCmSearchCosts.checkCommand(newCommand.replace("--default", "--max"), "model.cm", 20, 15, true);
		}}, "Unmatched acceleration profile must fail");

		write(control+"/command.txt", oldCommand+"\n");
		write(point+"/command.txt", newCommand+"\n");
		write(control+"/time.tsv", "19\t1\t21\t1000\n");
		write(point+"/time.tsv", "39\t1\t41\t1000\n");
		write(control+"/exit.txt", "0\n");
		write(point+"/exit.txt", "0\n");
		write(control+"/search.out", "[ok]\n");
		write(point+"/search.out", "[ok]\n");
		write(pilot+"/exit.txt", "0\n");
		write(pilot+"/result.txt", "SUBGA_NATIVE_TIMING_COMPLETE_NOT_REFERENCE_OR_ADMISSION\n");
		write(pilot+"/bindings.tsv", "accession\tGA_bits\treport_bits\tmodel_sha80\tmodel_path\tcontrol_command\tcontrol_hits\tcontrol_report\n"
			+"RF00050\t20\t15\t00000000000000000000\tmodel.cm\t"+control+"/command.txt\t"+control+"/hits.tbl\t"+control+"/search.out\n");
		write(out+"/census.tsv", "genomes\tinput_file_bytes\tnucleotide_bases\tfasta_records\tboth_strands_search_space_Mb\n4\t2200\t2000\t4\t0.004\n");
		write(out+"/population.tsv", "stage\tgenomes\tnucleotide_bases\tfasta_records\n2\t2\t1000\t2\n");
		CompareCmSearchCosts.main(new String[]{"in="+pilot, "out="+out+"/report", "census="+out+"/census.tsv", "population="+out+"/population.tsv", "stage=2"});
		metric(out+"/report.cost_summary.tsv", "projected_new_pass_CPU", .022222, "CPU_hours", "LINEAR_BASE_SCALING_NOT_MEASURED_FULL_SCAN");
		metric(out+"/report.cost_summary.tsv", "candidate_to_control_CPU_ratio", 2, "ratio", "NONCONCURRENT_TIMINGS");
	}

	private static String hit(final String target, final String strand, final String score, final String included){
		return target+" - FMN RF00050 cm 1 100 1 100 "+strand+" no 1 0.50 0.0 "+score+" 1e-5 "+included+" synthetic\n";
	}
	private static void write(final String path, final String text){
		try(Output out=new Output(path)){out.writer.print(text);}
	}
	/** Compares numeric values at the declared precision, allowing unpadded integers. */
	private static void metric(final String path, final String name, final double expected,
		final String units, final String interpretation){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(!p.termEquals(name, 0)){continue;}
				found++;
				check(p.terms()==4 && p.termEquals(units, 2) && p.termEquals(interpretation, 3), "Metric units/meaning changed: "+name);
				check(Math.abs(p.parseDouble(1)-expected)<0.00000051, "Hand-calculated metric differs: "+name);
			}
		}
		check(found==1, "Expected one metric row: "+name);
	}
	private static void expectFailure(final Runnable action, final String message){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, message);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
