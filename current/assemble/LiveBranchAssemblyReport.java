package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import assemble.SameKAssemblyReport.ReportRow;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Summarizes two to four named QUAST arms, retaining empty and undefined outcomes.
 * The first arm is the comparison baseline. Empty-marker arms contribute zero metrics;
 * undefined NGA50 also contributes zero, with a separate defined indicator. No genome
 * is dropped from quality summaries because of either condition or a timing exclusion.
 * Validates collected files rather than rerunning QUAST or verifying assembly inputs,
 * run settings or biological fold separation. NGA50, base and timing metrics at or above
 * 2^63 use round-tripping Java double text; smaller values retain the six-decimal
 * formatter. Nonfinite timing ratios reject before output, even for excluded pairs.
 * Count metrics retain exact longs; their medians retain exact integer/half-integer text.
 * Quality and timing fields guard decimal-accumulator overflow before domain checks;
 * their stored doubles still have ordinary binary rounding.
 * Output files are written separately.
 * @author Fischl
 */
public final class LiveBranchAssemblyReport{

	/** Uses standard config expansion and rejects unknown flags. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchAssemblyReport.class, false);
		try{new LiveBranchAssemblyReport(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Requires an existing input directory, a new output directory and an explicit positive
	 * expected panel size. Two to four distinct arm labels are allowed, with the first used
	 * as baseline; fold defaults to validation and may instead be confirmation.
	 * @param args Native directory options plus expected, arms, fold and timing controls
	 * @param outstream_ Borrowed diagnostic destination, closed by main when owned
	 */
	private LiveBranchAssemblyReport(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String key=split[0].toLowerCase(Locale.ROOT);
			while(key.startsWith("-")){key=key.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(key.equals("expected")){expected=Integer.parseInt(value);}
			else if(key.equals("arms")){arms=value==null ? new String[0] : value.split(",", -1);}
			else if(key.equals("fold")){fold=value;}
			else if(key.equals("timings")){timings=Parse.parseBoolean(value);}
			else if(key.equals("timingexclude")){
				if(value==null || value.isEmpty()){throw new IllegalArgumentException("timingexclude requires genome IDs.");}
				timingExclusions.addAll(Arrays.asList(value.split(",", -1)));
			}else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown report argument: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(!"validation".equals(fold) && !"confirmation".equals(fold)){
			throw new IllegalArgumentException("fold must be validation or confirmation; training data are not an evaluation panel.");
		}
		if(arms.length<2 || arms.length>4){
			throw new IllegalArgumentException("arms requires two to four distinct lowercase labels; first is the baseline.");
		}
		final HashSet<String> labels=new HashSet<String>();
		for(String arm : arms){
			if(!arm.matches("[a-z][a-z0-9_]*") || !labels.add(arm)){
				throw new IllegalArgumentException("Unsafe or duplicate arm label: "+arm);
			}
		}
		if(input==null || output==null || expected<1 || !new File(input).isDirectory() || new File(output).exists()){
			throw new IllegalArgumentException("Required in=collected_directory out=NEW_directory expected=genomes.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads exactly expected manifest genomes before creating output. Checks seven columns,
	 * the genome/fold/run header fields, unique gNNNN IDs, matching fold, K in [24,127] and
	 * positive depth. Other row fields are not interpreted; reports are read from runs/ID
	 * regardless of the run-column value. Optional timing inputs are also read before output.
	 * Separate writes can leave an incomplete result directory if a later write fails.
	 */
	private void run(){
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<ReportRow[]> panel=new ArrayList<ReportRow[]>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile manifest=ByteFile.makeByteFile(input+"/manifest.tsv", true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			lp.set(manifest.nextLine());
			if(lp.terms()!=7 || !lp.termEquals("genome", 0) || !lp.termEquals("fold", 1) || !lp.termEquals("run", 6)){
				throw new IllegalArgumentException("Expected the seven-column assembly-evaluation manifest.");
			}
			for(byte[] line=manifest.nextLine(); line!=null; line=manifest.nextLine()){
				lp.set(line);
				if(lp.terms()!=7 || !lp.termEquals(fold, 1) || lp.parseInt(2)<24 || lp.parseInt(2)>127 || lp.parseInt(3)<1){
					throw new IllegalArgumentException("Wrong fold, K, depth or manifest width.");
				}
				final String id=lp.parseString(0);
				if(!id.matches("g[0-9]{4}") || !seen.add(id)){throw new IllegalArgumentException("Unsafe/duplicate genome: "+id);}
				ids.add(id); panel.add(readReport(input+"/runs/"+id, arms));
			}
		}finally{if(manifest.close()){throw new IllegalStateException("Manifest read failed.");}}
		if(ids.size()!=expected){throw new IllegalArgumentException("Expected "+expected+" genomes, found "+ids.size());}
		if(!seen.containsAll(timingExclusions) || (!timings && !timingExclusions.isEmpty())){
			throw new IllegalArgumentException("Timing exclusions must name panel genomes and require timings=t.");
		}
		final TimingData timing=timings ? readTimings(ids) : null;
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create output: "+output);}
		writeRows(ids, panel); writeSummary(panel);
		if(timing!=null){writeTimings(ids, timing);}
		outstream.println("LIVE_BRANCH_ASSEMBLY_REPORT_PASS genomes="+ids.size()+" arms="+arms.length);
	}

	/** Emits each genome/arm, including explicit undefined-NGA50 flags. */
	private void writeRows(final ArrayList<String> ids, final ArrayList<ReportRow[]> panel){
		assert(ids.size()==panel.size()) : "Each validated manifest genome supplies exactly one paired report.";
		final ByteStreamWriter out=open(output+"/metrics.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			out.print("genome\tarm\tcontigs\tnga50_zero\tglobal\tlocal\tgenome_fraction\tbases\tnga50_defined\n");
			for(int g=0; g<ids.size(); g++){
				for(int a=0; a<arms.length; a++){
					row.clear().append(ids.get(g)).tab().append(arms[a]);
					final ReportRow values=panel.get(g)[a];
					row.tab().append(values.contigs).tab();
					SameKAssemblyReport.appendMetric(row, values.nga50).tab().append(values.global).tab().append(values.local)
						.tab().append(values.genomeFraction, 6).tab();
					SameKAssemblyReport.appendMetric(row, values.bases).tab().append(values.nga50Defined ? 1 : 0);
					out.print(row.nl());
				}
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Metrics write failed.");}}
	}

	/**
	 * Unweighted genome means/medians and count totals retain missing NGA50 as zero.
	 * Comparisons use the first arm; genome fraction is a percentage. Global-error equality
	 * is not separately counted. Empty-marker arms retain their zero rows in all summaries.
	 * Count totals reject long overflow; earlier output files may already have been written.
	 */
	private void writeSummary(final ArrayList<ReportRow[]> panel){
		assert(!panel.isEmpty()) : "A complete nonempty manifest is required for means and medians.";
		final ByteStreamWriter out=open(output+"/summary.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			out.print("arm\tgenomes\tundefined_nga50\tglobal_total\tlocal_total\tcontigs_total\tcontigs_median\tnga50_zero_mean\tnga50_zero_median\tgenome_fraction_mean\tnga50_better_"+arms[0]+
				"\tnga50_equal_"+arms[0]+"\tnga50_worse_"+arms[0]+"\tglobal_more_"+arms[0]+"\tglobal_less_"+arms[0]+"\n");
			for(int a=0; a<arms.length; a++){
				final double[] nga=new double[panel.size()];
				final long[] counts=new long[panel.size()];
				long global=0, local=0, contigs=0;
				int undefined=0, better=0, equal=0, worse=0, moreErrors=0, fewerErrors=0;
				double sum=0, gf=0;
				for(int g=0; g<panel.size(); g++){
					final ReportRow values=panel.get(g)[a], baseline=panel.get(g)[0];
					counts[g]=values.contigs; nga[g]=values.nga50;
					contigs=Math.addExact(contigs, values.contigs);
					global=Math.addExact(global, values.global);
					local=Math.addExact(local, values.local);
					sum+=values.nga50; gf+=values.genomeFraction; undefined+=values.nga50Defined ? 0 : 1;
					if(values.nga50>baseline.nga50){better++;}else if(values.nga50<baseline.nga50){worse++;}else{equal++;}
					if(values.global>baseline.global){moreErrors++;}else if(values.global<baseline.global){fewerErrors++;}
				}
				// ASM033: retain the finite fast path and scale before median sorts NGA50 values.
				final double mean=Double.isFinite(sum) ? sum/panel.size() : SameKAssemblyReport.timingMean(nga, sum);
				row.clear().append(arms[a]).tab().append(panel.size()).tab().append(undefined).tab().append(global).tab()
					.append(local).tab().append(contigs).tab();
				SameKAssemblyReport.appendCountMedian(row, counts).tab();
				SameKAssemblyReport.appendMetric(row, mean).tab();
				SameKAssemblyReport.appendMetric(row, SameKAssemblyReport.median(nga)).tab().append(gf/panel.size(), 6)
					.tab().append(better).tab().append(equal).tab().append(worse).tab().append(moreErrors).tab().append(fewerErrors);
				out.print(row.nl());
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
	}

	/** Validates wall/user/system triples and classifies each treatment/baseline pair.
	 * The ten-second threshold uses wall time, separately for each pair; matching
	 * host and duration alone do not establish equal inputs, threads or JVM settings, nor
	 * ten seconds of compute. Exclusion precedence is requested genome, different host,
	 * short wall time, then included. All arms' files are required even for excluded genomes.
	 * Each nonbaseline arm must retain at least one pair; their retained subsets may differ.
	 * Every arm/baseline ratio must be finite, including excluded pairs, before output. */
	private TimingData readTimings(final ArrayList<String> ids){
		assert(!ids.isEmpty()) : "Timing dimensions follow the validated nonempty genome manifest.";
		final TimingData data=new TimingData(ids.size(), arms.length);
		final LineParser1 lp=new LineParser1(' ');
		for(int g=0; g<ids.size(); g++){
			for(int a=0; a<arms.length; a++){
				final String base=input+"/runs/"+ids.get(g)+"/runs/"+arms[a]+"/";
				lp.set(singleLine(base+"elapsed.tsv"));
				if(lp.terms()!=3){throw new IllegalArgumentException("Expected elapsed/user/system seconds: "+base);}
				for(int m=0; m<3; m++){
					final double value=SameKAssemblyReport.parseReportDouble(lp, m);
					if(!Double.isFinite(value) || value<0 || (m==0 && value==0)){
						throw new IllegalArgumentException("Invalid native timing: "+base+"="+value);
					}
					data.seconds[g][a*3+m]=value;
				}
				data.hosts[g][a]=new String(singleLine(base+"host.txt"), java.nio.charset.StandardCharsets.US_ASCII);
				if(!data.hosts[g][a].matches("[A-Za-z0-9.-]+")){throw new IllegalArgumentException("Invalid host: "+base);}
			}
			for(int a=1; a<arms.length; a++){
				final double armSeconds=data.seconds[g][a*3], baseSeconds=data.seconds[g][0];
				if(!Double.isFinite(armSeconds/baseSeconds)){
					throw new IllegalArgumentException("Nonfinite timing ratio: genome="+ids.get(g)+
						" arm="+arms[a]+" baseline="+arms[0]+" wall="+armSeconds+" baseline_wall="+baseSeconds);
				}
				data.reasons[g][a]=timingExclusions.contains(ids.get(g)) ? "excluded_by_design" :
					!data.hosts[g][0].equals(data.hosts[g][a]) ? "different_hosts" :
					Math.min(data.seconds[g][0], data.seconds[g][a*3])<10 ? "under_10_seconds" : "included";
				if(data.reasons[g][a].equals("included")){data.included[a]++;}
			}
		}
		for(int a=1; a<arms.length; a++){
			if(data.included[a]==0){throw new IllegalArgumentException("No comparable timing pairs for "+arms[a]);}
		}
		return data;
	}

	/**
	 * Retains excluded rows and uses a separately matched denominator for each nonbaseline arm.
	 * Ratios are arm/baseline wall seconds. Mean/median paired ratios treat genomes equally;
	 * ratio_of_totals divides summed arm seconds by summed baseline seconds. User/system
	 * seconds are retained in rows but do not define these wall-time summaries.
	 * Finite sums retain their original arithmetic; overflow uses scaled means.
	 * Large finite timing values use the same formatter as assembly-quality metrics.
	 */
	private void writeTimings(final ArrayList<String> ids, final TimingData data){
		assert(ids.size()==data.seconds.length) : "Every manifest genome must retain its timing/exclusion rows.";
		final ByteStreamWriter rows=open(output+"/timings.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			rows.print("genome\tarm\tbaseline_host\tarm_host\tbaseline_wall\tbaseline_user\tbaseline_system\tarm_wall\tarm_user\tarm_system\tarm_over_baseline\tstatus\n");
			for(int g=0; g<ids.size(); g++){
				for(int a=1; a<arms.length; a++){
					final double[] seconds=data.seconds[g];
					row.clear().append(ids.get(g)).tab().append(arms[a]).tab().append(data.hosts[g][0]).tab().append(data.hosts[g][a]);
					for(int m=0; m<3; m++){SameKAssemblyReport.appendMetric(row.tab(), seconds[m]);}
					for(int m=0; m<3; m++){SameKAssemblyReport.appendMetric(row.tab(), seconds[a*3+m]);}
					rows.print(SameKAssemblyReport.appendMetric(row.tab(), seconds[a*3]/seconds[0]).tab().append(data.reasons[g][a]).nl());
				}
			}
		}finally{if(rows.poisonAndWait()){throw new IllegalStateException("Timing rows write failed.");}}
		final ByteStreamWriter summary=open(output+"/timing_summary.tsv");
		try{
			summary.print("arm\tpairs\texcluded\tbaseline_mean_seconds\tarm_mean_seconds\tbaseline_median_seconds\tarm_median_seconds\tmean_paired_ratio\tmedian_paired_ratio\tratio_of_totals\n");
			for(int a=1; a<arms.length; a++){
				final int n=data.included[a];
				final double[] baseline=new double[n], treatment=new double[n], ratios=new double[n];
				double sumBase=0, sumArm=0, sumRatio=0;
				int index=0;
				for(int g=0; g<ids.size(); g++){
					if(!data.reasons[g][a].equals("included")){continue;}
					baseline[index]=data.seconds[g][0]; treatment[index]=data.seconds[g][a*3];
					ratios[index]=treatment[index]/baseline[index];
					// Keep original addition order; timingMean handles positive overflow below.
					sumBase+=baseline[index]; sumArm+=treatment[index]; sumRatio+=ratios[index]; index++;
				}
				assert(index==n && n>0) : "Each included genome contributes once to the paired timing denominator.";
				// Sums match the retained arrays in their original order, before median sorts.
				final double meanBase=SameKAssemblyReport.timingMean(baseline, sumBase);
				final double meanArm=SameKAssemblyReport.timingMean(treatment, sumArm);
				final double meanRatio=SameKAssemblyReport.timingMean(ratios, sumRatio);
				final double totalRatio=Double.isFinite(sumArm) && Double.isFinite(sumBase) ?
						sumArm/sumBase : meanArm/meanBase;
				if(!Double.isFinite(totalRatio)){
					throw new IllegalStateException("Timing total ratio must be finite for matched wall times >=10: arm="+arms[a]+" ratio="+totalRatio);
				}
				row.clear().append(arms[a]).tab().append(n).tab().append(ids.size()-n).tab();
				SameKAssemblyReport.appendMetric(row, meanBase).tab();
				SameKAssemblyReport.appendMetric(row, meanArm).tab();
				SameKAssemblyReport.appendMetric(row, SameKAssemblyReport.median(baseline)).tab();
				SameKAssemblyReport.appendMetric(row, SameKAssemblyReport.median(treatment)).tab();
				SameKAssemblyReport.appendMetric(row, meanRatio).tab();
				SameKAssemblyReport.appendMetric(row, SameKAssemblyReport.median(ratios)).tab();
				SameKAssemblyReport.appendMetric(row, totalRatio);
				summary.print(row.nl());
			}
		}finally{if(summary.poisonAndWait()){throw new IllegalStateException("Timing summary write failed.");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra or absent timing records indicate a malformed native artifact, not a missing observation. */
	private static byte[] singleLine(final String path){
		assert(path!=null) : "Timing paths are constructed from validated genome and arm labels.";
		final ByteFile file=ByteFile.makeByteFile(path, true);
		try{
			final byte[] line=file.nextLine();
			if(line==null || line.length==0 || file.nextLine()!=null){throw new IllegalArgumentException("Expected one nonempty line: "+path);}
			return line;
		}finally{if(file.close()){throw new IllegalStateException("Timing read failed: "+path);}}
	}

	/**
	 * Uses regular ARM.empty marker files to identify omitted arms; marker contents and
	 * underlying assemblies are not inspected. Marked arms retain seven zero values.
	 * If every arm is marked empty, the QUAST report must be absent; otherwise columns contain
	 * exactly the unmarked arms in configured order. Each required metric must occur once;
	 * other metric names are ignored after checking row width. Only NGA50 may be '-', stored
	 * as zero with an unset defined flag. Values are finite/nonnegative, count/length metrics
	 * integral, and genome fraction at most 100 percent. Contig/global/local counts must
	 * be below 2^63 after parsing to double; near-limit integer tokens that round to 2^63
	 * are also rejected. Counts additionally require exact integral tokens and retain longs;
	 * the parsed-double cap remains, so not every nonnegative signed-long token is admitted.
	 * @param directory Collected per-genome directory, containing markers and quast/report.tsv
	 * @param arms Ordered expected arm labels, including any explicitly empty arms
	 * @return Per-arm exact counts and double quality metrics, including empty/undefined flags
	 */
	private static ReportRow[] readReport(final String directory, final String[] arms){
		assert(directory!=null) : "The report directory is constructed from a validated genome ID.";
		assert(arms.length>=2 && arms.length<=4) : "Validated arm labels map QUAST columns and explicit empty markers.";
		final boolean[] empty=new boolean[arms.length];
		final ReportRow[] data=new ReportRow[arms.length];
		for(int a=0; a<arms.length; a++){data[a]=new ReportRow();}
		int n=0;
		for(int a=0; a<arms.length; a++){
			empty[a]=new File(directory+"/"+arms[a]+".empty").isFile();
			if(!empty[a]){n++;}
		}
		final String path=directory+"/quast/report.tsv";
		if(n==0){
			if(new File(path).exists()){throw new IllegalArgumentException("QUAST report contradicts all-empty markers: "+directory);}
			return data;
		}
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		final boolean[] found=new boolean[METRICS.length];
		try{
			lp.set(file.nextLine());
			if(lp.terms()!=n+1 || !lp.termEquals("Assembly", 0)){throw new IllegalArgumentException("Wrong QUAST header: "+path);}
			for(int a=0, column=1; a<arms.length; a++){
				if(!empty[a] && !lp.termEquals(arms[a], column++)){throw new IllegalArgumentException("Wrong QUAST arm order: "+path);}
			}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				if(lp.terms()!=n+1){throw new IllegalArgumentException("Wrong QUAST row width: "+path);}
				for(int m=0; m<METRICS.length; m++){
					if(!lp.termEquals(METRICS[m], 0)){continue;}
					if(found[m]){throw new IllegalArgumentException("Duplicate metric: "+METRICS[m]);}
					found[m]=true;
					for(int a=0, column=1; a<arms.length; a++){
						if(empty[a]){continue;}
						final boolean missing=lp.termEquals("-", column);
						if(missing && m!=1){throw new IllegalArgumentException("Undefined required metric: "+METRICS[m]);}
						final double value=missing ? 0 : (m==1 || m==4 || m==5) ?
							SameKAssemblyReport.parseReportDouble(lp, column) : lp.parseDouble(column);
						if(!Double.isFinite(value) || value<0 || (m!=4 && value!=Math.rint(value)) || (m==4 && value>100)){
							throw new IllegalArgumentException("Invalid QUAST metric: "+METRICS[m]+"="+value);
						}
						if((m==0 || m==2 || m==3) && value>=0x1p63){
							throw new IllegalArgumentException("QUAST count must be below 2^63 after double parsing: "+METRICS[m]+"="+value);
						}
						SameKAssemblyReport.setMetric(data[a], m, value, missing, lp, column, path);
						column++;
					}
				}
			}
		}finally{if(file.close()){throw new IllegalStateException("QUAST read failed.");}}
		for(int m=0; m<found.length; m++){if(!found[m]){throw new IllegalArgumentException("Missing metric: "+METRICS[m]);}}
		return data;
	}

	/** Output files are new and every caller checks the writer's terminal error state. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter out=new ByteStreamWriter(path, false, false, false);
		out.start(); return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Small retained native records; all filtering is explicit and independent of assembly quality. */
	private static final class TimingData{
		/** Uses the validated manifest and arm count, never a guessed number of files. */
		TimingData(final int genomes, final int arms){
			assert(genomes>0 && arms>=2 && arms<=4) : "Timing rows require a nonempty panel and two to four validated arms.";
			seconds=new double[genomes][arms*3]; hosts=new String[genomes][arms];
			reasons=new String[genomes][arms]; included=new int[arms];
		}
		/** Per-genome wall/user/system triples in configured arm order. */
		final double[][] seconds;
		/** Literal host labels and first applicable pair exclusions; baseline reasons remain null. */
		final String[][] hosts, reasons;
		/** Included pair count per treatment arm; baseline slot remains zero. */
		final int[] included;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Collected panel root and separately created result directory. */
	private final String input, output;
	/** Borrowed destination for completion diagnostics. */
	private final PrintStream outstream;
	/** Exact required manifest row count; caller must replace the invalid default. */
	private int expected=-1;
	/** Required label on every manifest row, without biological partition verification. */
	private String fold="validation";
	/** Enables timing collection for every configured arm. */
	private boolean timings=false;
	/** Genome IDs excluded from all timing comparisons, but retained in quality summaries. */
	private final HashSet<String> timingExclusions=new HashSet<String>();
	/** Ordered quality/timing labels; the first arm anchors comparisons. */
	private String[] arms={"native", "live"};

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Required QUAST row labels in the shared setMetric order; only NGA50 may be undefined. */
	private static final String[] METRICS={"# contigs", "NGA50", "# misassemblies", "# local misassemblies", "Genome fraction (%)", "Total length"};
}
