package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Summarizes matched five-arm QUAST reports without selecting a model or cutoff.
 * @author Fischl
 */
public final class SameKAssemblyReport {

	/** Standard BBTools argument/config parsing; input is a collected matched assembly panel. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKAssemblyReport.class, false);
		try{new SameKAssemblyReport(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires a new output directory, protecting both native inputs and previous results. */
	private SameKAssemblyReport(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String key=split[0].toLowerCase(Locale.ROOT);
			while(key.startsWith("-")){key=key.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(key.equals("expected")){expected=Integer.parseInt(value);}
			else if(key.equals("arms")){arms=value.split(",", -1);}
			else if(key.equals("fold")){fold=value;}
			else if(key.equals("timings")){timings=Parse.parseBoolean(value);}
			else if(key.equals("timingexclude")){timingExclusions.addAll(Arrays.asList(value.split(",")));}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown report argument: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(!"validation".equals(fold) && !"test".equals(fold)){
			throw new IllegalArgumentException("fold must be validation or test; held-out reports cannot silently include training genomes.");
		}
		final HashSet<String> names=new HashSet<String>();
		if(arms.length!=5 || !arms[0].equals("baseline") || !arms[4].equals("ordinary")){
			throw new IllegalArgumentException("Five arms required, with baseline first and ordinary last.");
		}
		for(String name : arms){
			if(!name.matches("[a-z][a-z0-9_]*") || !names.add(name)){throw new IllegalArgumentException("Invalid/duplicate arm: "+name);}
		}
		if(input==null || output==null || expected<1 || !new File(input).isDirectory() || new File(output).exists()){
			throw new IllegalArgumentException("Required in=collected_directory out=NEW_directory expected=genomes.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads every frozen genome before opening outputs; incomplete panels cannot silently shrink. */
	private void run(){
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<double[][]> panel=new ArrayList<double[][]>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile manifest=ByteFile.makeByteFile(input+"/manifest.tsv", true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			lp.set(manifest.nextLine());
			if(lp.terms()!=6 || !lp.termEquals("genome", 0) || !lp.termEquals("fold", 1)){
				throw new IllegalArgumentException("Wrong assembly manifest header.");
			}
			for(byte[] line=manifest.nextLine(); line!=null; line=manifest.nextLine()){
				lp.set(line);
				if(lp.terms()!=6 || !lp.termEquals(fold, 1)){
					throw new IllegalArgumentException("Expected six-column "+fold+" manifest row; mixed or mislabeled folds are forbidden.");
				}
				final String id=lp.parseString(0);
				if(!id.matches("g[0-9]{4}") || !seen.add(id)){throw new IllegalArgumentException("Unsafe/duplicate genome: "+id);}
				ids.add(id);
				panel.add(readReport(input+"/runs/"+id+"/quast/report.tsv", arms));
			}
		}finally{
			if(manifest.close()){throw new IllegalStateException("Manifest read failed.");}
		}
		if(ids.size()!=expected){throw new IllegalArgumentException("Expected "+expected+" genomes, found "+ids.size());}
		if(!seen.containsAll(timingExclusions) || (!timings && !timingExclusions.isEmpty())){
			throw new IllegalArgumentException("Timing exclusions must name panel genomes and require timings=t.");
		}
		final TimingData time=timings ? readTimings(ids) : null;
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create output: "+output);}
		writeRows(ids, panel);
		writeSummary(panel);
		if(time!=null){writeTimings(ids, time);}
		outstream.println("SAME_K_ASSEMBLY_REPORT_PASS genomes="+ids.size()+" arms="+arms.length);
	}

	/** Retains per-genome metrics, including a separate flag for undefined NGA50. */
	private void writeRows(final ArrayList<String> ids, final ArrayList<double[][]> panel){
		assert(ids.size()==panel.size()) : "Each frozen genome must have one matched five-arm report.";
		final ByteStreamWriter out=open(output+"/metrics.tsv");
		final ByteBuilder bb=new ByteBuilder();
		out.print("genome\tarm\tcontigs\tnga50_zero\tglobal\tlocal\tgenome_fraction\tbases\tnga50_defined\n");
		try{
			for(int g=0; g<ids.size(); g++){
				for(int a=0; a<arms.length; a++){
					bb.clear().append(ids.get(g)).tab().append(arms[a]);
					for(double value : panel.get(g)[a]){bb.tab().append(value, 6);}
					out.print(bb.nl());
				}
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Metrics write failed.");}}
	}

	/** Reports unweighted genome means/medians and paired wins, without optimizing thresholds. */
	private void writeSummary(final ArrayList<double[][]> panel){
		assert(!panel.isEmpty()) : "Summary denominators require a nonempty matched genome panel.";
		final ByteStreamWriter out=open(output+"/summary.tsv");
		final ByteBuilder bb=new ByteBuilder();
		out.print("arm\tgenomes\tundefined_nga50\tglobal_total\tlocal_total\tcontigs_total\tcontigs_median\tnga50_zero_mean\tnga50_zero_median\tgenome_fraction_mean\tnga50_better_raw\tnga50_equal_raw\tnga50_worse_raw\tnga50_better_ordinary\tnga50_equal_ordinary\tnga50_worse_ordinary\tglobal_more_raw\tglobal_less_raw\tglobal_more_ordinary\tglobal_less_ordinary\n");
		try{
			for(int a=0; a<arms.length; a++){
				final double[] nga=new double[panel.size()], counts=new double[panel.size()];
				long global=0, local=0, contigs=0;
				int undefined=0;
				double sum=0, gf=0;
				final int[] comparisons=new int[10];
				for(int g=0; g<panel.size(); g++){
					final double[][] group=panel.get(g);
					final double[] row=group[a];
					counts[g]=row[0]; nga[g]=row[1];
					contigs+=(long)row[0]; global+=(long)row[2]; local+=(long)row[3];
					sum+=row[1]; gf+=row[4]; undefined+=row[6]==0 ? 1 : 0;
					comparisons[compareSlot(row[1], group[0][1])]++;
					comparisons[3+compareSlot(row[1], group[4][1])]++;
					if(row[2]>group[0][2]){comparisons[6]++;}
					if(row[2]<group[0][2]){comparisons[7]++;}
					if(row[2]>group[4][2]){comparisons[8]++;}
					if(row[2]<group[4][2]){comparisons[9]++;}
				}
				bb.clear().append(arms[a]).tab().append(panel.size()).tab().append(undefined).tab()
					.append(global).tab().append(local).tab().append(contigs).tab().append(median(counts), 6)
					.tab().append(sum/panel.size(), 6).tab().append(median(nga), 6).tab().append(gf/panel.size(), 6);
				for(int value : comparisons){bb.tab().append(value);}
				out.print(bb.nl());
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
	}

	/** Reads native GNU-time triples before opening outputs; excludes incomparable pairs explicitly. */
	private TimingData readTimings(final ArrayList<String> ids){
		assert(!ids.isEmpty()) : "Timing rows must follow the already validated nonempty genome manifest.";
		final TimingData data=new TimingData(ids.size());
		for(int g=0; g<ids.size(); g++){
			final String base=input+"/runs/"+ids.get(g)+"/";
			for(int a=0; a<2; a++){
				final String arm=a==0 ? "live" : "ordinary";
				final byte[] line=singleLine(base+arm+"/elapsed.tsv");
				final LineParser1 lp=new LineParser1(' ');
				lp.set(line);
				if(lp.terms()!=3){throw new IllegalArgumentException("Expected elapsed/user/system seconds: "+base+arm);}
				for(int m=0; m<3; m++){
					final double value=lp.parseDouble(m);
					if(!Double.isFinite(value) || value<0 || (m==0 && value==0)){
						throw new IllegalArgumentException("Invalid native timing: "+base+arm+"="+value);
					}
					data.seconds[g][a*3+m]=value;
				}
				data.hosts[g][a]=new String(singleLine(base+arm+"/host.txt"), java.nio.charset.StandardCharsets.US_ASCII);
				if(!data.hosts[g][a].matches("[A-Za-z0-9.-]+")){throw new IllegalArgumentException("Invalid timing host: "+base+arm);}
			}
			data.reason[g]=timingExclusions.contains(ids.get(g)) ? "excluded_by_design" :
				!data.hosts[g][0].equals(data.hosts[g][1]) ? "different_hosts" :
				Math.min(data.seconds[g][0], data.seconds[g][3])<10 ? "under_10_seconds" : "included";
			if(data.reason[g].equals("included")){data.included++;}
		}
		if(data.included==0){throw new IllegalArgumentException("No comparable timing pairs remain; refusing an empty speed estimate.");}
		return data;
	}

	/** Writes every pair and the explicitly filtered summary; each genome has equal weight. */
	private void writeTimings(final ArrayList<String> ids, final TimingData data){
		assert(ids.size()==data.seconds.length && data.included>0) : "Timing summary requires one row per genome and a nonempty comparable subset.";
		final ByteStreamWriter rows=open(output+"/timings.tsv");
		final ByteBuilder bb=new ByteBuilder();
		final double[] live=new double[data.included], ordinary=new double[data.included], ratios=new double[data.included];
		double sumLive=0, sumOrdinary=0, sumRatio=0;
		int index=0;
		rows.print("genome\tlive_host\tordinary_host\tlive_wall\tlive_user\tlive_system\tordinary_wall\tordinary_user\tordinary_system\tlive_over_ordinary\tstatus\n");
		try{
			for(int g=0; g<ids.size(); g++){
				final double[] seconds=data.seconds[g];
				final double ratio=seconds[0]/seconds[3];
				bb.clear().append(ids.get(g)).tab().append(data.hosts[g][0]).tab().append(data.hosts[g][1]);
				for(double value : seconds){bb.tab().append(value, 6);}
				rows.print(bb.tab().append(ratio, 6).tab().append(data.reason[g]).nl());
				if(!data.reason[g].equals("included")){continue;}
				live[index]=seconds[0]; ordinary[index]=seconds[3]; ratios[index]=ratio; index++;
				sumLive+=seconds[0]; sumOrdinary+=seconds[3]; sumRatio+=ratio;
			}
		}finally{if(rows.poisonAndWait()){throw new IllegalStateException("Timing rows write failed.");}}
		assert(index==data.included) : "All included timing pairs must contribute exactly once to summary denominators.";
		final ByteStreamWriter summary=open(output+"/timing_summary.tsv");
		try{
			summary.print("pairs\texcluded\tlive_mean_seconds\tordinary_mean_seconds\tlive_median_seconds\tordinary_median_seconds\tmean_paired_ratio\tmedian_paired_ratio\tratio_of_totals\n");
			bb.clear().append(index).tab().append(ids.size()-index).tab().append(sumLive/index, 6)
				.tab().append(sumOrdinary/index, 6).tab().append(median(live), 6).tab().append(median(ordinary), 6)
				.tab().append(sumRatio/index, 6).tab().append(median(ratios), 6).tab().append(sumLive/sumOrdinary, 6);
			summary.print(bb.nl());
		}finally{if(summary.poisonAndWait()){throw new IllegalStateException("Timing summary write failed.");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Native timing/host files contain exactly one record; extra lines indicate a failed command or wrong artifact. */
	private static byte[] singleLine(final String path){
		assert(path!=null) : "Native artifact paths are constructed from validated manifest IDs.";
		final ByteFile file=ByteFile.makeByteFile(path, true);
		try{
			final byte[] line=file.nextLine();
			if(line==null || line.length==0 || file.nextLine()!=null){throw new IllegalArgumentException("Expected one nonempty line: "+path);}
			return line;
		}finally{if(file.close()){throw new IllegalStateException("Native timing read failed: "+path);}}
	}

	/** Requires fixed arm order and every named metric exactly once; only NGA50 may be undefined. */
	static double[][] readReport(final String path, final String[] arms){
		assert(arms.length==5) : "The matched reporter compares five ordered arms, with raw first and ordinary last.";
		final double[][] data=new double[arms.length][7];
		final boolean[] found=new boolean[METRICS.length];
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			lp.set(file.nextLine());
			if(lp.terms()!=6 || !lp.termEquals("Assembly", 0)){throw new IllegalArgumentException("Wrong QUAST header: "+path);}
			for(int a=0; a<arms.length; a++){
				if(!lp.termEquals(arms[a], a+1)){throw new IllegalArgumentException("Unexpected QUAST arm order: "+path);}
			}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				if(lp.terms()!=6){throw new IllegalArgumentException("Wrong QUAST row width: "+path);}
				for(int m=0; m<METRICS.length; m++){
					if(!lp.termEquals(METRICS[m], 0)){continue;}
					if(found[m]){throw new IllegalArgumentException("Duplicate QUAST metric: "+METRICS[m]);}
					found[m]=true;
					for(int a=0; a<arms.length; a++){
						final boolean missing=lp.termEquals("-", a+1);
						if(missing && m!=1){throw new IllegalArgumentException("Undefined required metric: "+METRICS[m]);}
						final double value=missing ? 0 : lp.parseDouble(a+1);
						if(!Double.isFinite(value) || value<0 || (m!=4 && value!=Math.rint(value)) || (m==4 && value>100)){
							throw new IllegalArgumentException("Invalid QUAST metric: "+METRICS[m]+"="+value);
						}
						data[a][m]=value;
						if(m==1){data[a][6]=missing ? 0 : 1;}
					}
				}
			}
		}finally{if(file.close()){throw new IllegalStateException("QUAST read failed: "+path);}}
		for(int m=0; m<found.length; m++){
			if(!found[m]){throw new IllegalArgumentException("Missing QUAST metric: "+METRICS[m]);}
		}
		return data;
	}

	/** Maps greater/equal/less values into the explicit three-column summary order. */
	private static int compareSlot(final double value, final double baseline){return value>baseline ? 0 : value==baseline ? 1 : 2;}

	/** Sorts a disposable copy of one metric and handles even-sized panels without rounding. */
	static double median(final double[] values){
		assert(values.length>0) : "Median requires at least one genome; missing NGA50 values are explicitly zero, not dropped.";
		Arrays.sort(values);
		final int mid=values.length/2;
		return values.length%2==0 ? (values[mid-1]+values[mid])/2 : values[mid];
	}

	/** Opens a new output and requires callers to observe terminal writer errors. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Small retained per-genome timing records, separate from assembly quality. */
	private static final class TimingData {
		/** Allocates fixed-size records only after the manifest has been validated. */
		TimingData(final int size){
			assert(size>0) : "Timing data corresponds to the validated nonempty panel.";
			seconds=new double[size][6]; hosts=new String[size][2]; reason=new String[size];
		}
		final double[][] seconds;
		final String[][] hosts;
		final String[] reason;
		int included;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final String input, output;
	private int expected=153;
	private String fold="validation";
	private boolean timings=false;
	private final HashSet<String> timingExclusions=new HashSet<String>();
	private String[] arms={"baseline", "oracle", "nn_three", "nn_four", "ordinary"};
	private static final String[] METRICS={"# contigs", "NGA50", "# misassemblies", "# local misassemblies", "Genome fraction (%)", "Total length"};
}
