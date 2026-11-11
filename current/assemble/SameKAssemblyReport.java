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
import shared.Tools;
import structures.ByteBuilder;

/**
 * Summarizes collected five-arm QUAST reports without selecting a model or cutoff.
 * Validates manifest IDs/fold labels and report columns; it does not rerun QUAST or
 * verify the underlying assembly inputs, run settings or biological fold separation.
 * Undefined NGA50 contributes zero to means, medians and comparisons, with a separate
 * defined flag/count. Timing exclusions affect only timing summaries, not assembly quality.
 * NGA50, base and timing metrics at or above 2^63 use round-tripping Java double text;
 * smaller values retain the six-decimal formatter. All inputs are loaded before
 * output creation, including finite timing-ratio checks even for excluded pairs.
 * Count metrics retain exact longs; their medians retain exact integer/half-integer text.
 * Quality and timing fields guard decimal-accumulator overflow before domain checks;
 * their stored doubles still have ordinary binary rounding.
 * Output files are written separately.
 * @author Fischl
 */
public final class SameKAssemblyReport{

	/** Standard BBTools argument/config parsing; input is a collected matched assembly panel. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKAssemblyReport.class, false);
		try{new SameKAssemblyReport(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Requires an existing collected-input directory and a new output directory, five distinct
	 * arm names with baseline first/ordinary last, and a validation/test fold label.
	 * The default fold is validation; timing exclusions are optional panel genome IDs.
	 * @param args Native directory options plus expected, arms, fold and timing controls
	 * @param stream Borrowed diagnostic destination, closed by main when owned
	 */
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

	/**
	 * Loads every manifest genome before opening outputs, requiring exactly expected rows.
	 * Manifest validation checks six columns, the genome/fold header fields, unique gNNNN IDs
	 * and matching fold labels; the other four manifest columns are not interpreted here.
	 * Reads reports beneath runs/ID/quast, then optional timing inputs. Separate output writes
	 * can leave an incomplete directory if a later write fails.
	 */
	private void run(){
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<ReportRow[]> panel=new ArrayList<ReportRow[]>();
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

	/** Emits all genome/arm metrics in manifest order, retaining zero-filled NGA50 and its defined flag. */
	private void writeRows(final ArrayList<String> ids, final ArrayList<ReportRow[]> panel){
		assert(ids.size()==panel.size()) : "Each frozen genome must have one matched five-arm report.";
		final ByteStreamWriter out=open(output+"/metrics.tsv");
		final ByteBuilder bb=new ByteBuilder();
		out.print("genome\tarm\tcontigs\tnga50_zero\tglobal\tlocal\tgenome_fraction\tbases\tnga50_defined\n");
		try{
			for(int g=0; g<ids.size(); g++){
				for(int a=0; a<arms.length; a++){
					bb.clear().append(ids.get(g)).tab().append(arms[a]);
					final ReportRow values=panel.get(g)[a];
					bb.tab().append(values.contigs).tab();
					appendMetric(bb, values.nga50).tab().append(values.global).tab().append(values.local)
						.tab().append(values.genomeFraction, 6).tab();
					appendMetric(bb, values.bases).tab().append(values.nga50Defined ? 1 : 0);
					out.print(bb.nl());
				}
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Metrics write failed.");}}
	}

	/**
	 * Reports unweighted genome means/medians, count totals and paired comparisons against
	 * arms 0 and 4. Zero-filled undefined NGA50 stays in every denominator and comparison;
	 * genome fraction remains a percentage. Global-error equality is not a separate summary column.
	 * Count totals reject long overflow; earlier output files may already have been written.
	 */
	private void writeSummary(final ArrayList<ReportRow[]> panel){
		assert(!panel.isEmpty()) : "Summary denominators require a nonempty matched genome panel.";
		final ByteStreamWriter out=open(output+"/summary.tsv");
		final ByteBuilder bb=new ByteBuilder();
		out.print("arm\tgenomes\tundefined_nga50\tglobal_total\tlocal_total\tcontigs_total\tcontigs_median\tnga50_zero_mean\tnga50_zero_median\tgenome_fraction_mean\tnga50_better_raw\tnga50_equal_raw\tnga50_worse_raw\tnga50_better_ordinary\tnga50_equal_ordinary\tnga50_worse_ordinary\tglobal_more_raw\tglobal_less_raw\tglobal_more_ordinary\tglobal_less_ordinary\n");
		try{
			for(int a=0; a<arms.length; a++){
				final double[] nga=new double[panel.size()];
				final long[] counts=new long[panel.size()];
				long global=0, local=0, contigs=0;
				int undefined=0;
				double sum=0, gf=0;
				final int[] comparisons=new int[10];
				for(int g=0; g<panel.size(); g++){
					final ReportRow[] group=panel.get(g);
					final ReportRow row=group[a];
					counts[g]=row.contigs; nga[g]=row.nga50;
					contigs=Math.addExact(contigs, row.contigs);
					global=Math.addExact(global, row.global);
					local=Math.addExact(local, row.local);
					sum+=row.nga50; gf+=row.genomeFraction; undefined+=row.nga50Defined ? 0 : 1;
					comparisons[compareSlot(row.nga50, group[0].nga50)]++;
					comparisons[3+compareSlot(row.nga50, group[4].nga50)]++;
					if(row.global>group[0].global){comparisons[6]++;}
					if(row.global<group[0].global){comparisons[7]++;}
					if(row.global>group[4].global){comparisons[8]++;}
					if(row.global<group[4].global){comparisons[9]++;}
				}
				// ASM033: scale only on overflow, before median sorts the original-order NGA50 array.
				final double mean=Double.isFinite(sum) ? sum/panel.size() : timingMean(nga, sum);
				bb.clear().append(arms[a]).tab().append(panel.size()).tab().append(undefined).tab()
					.append(global).tab().append(local).tab().append(contigs).tab();
				appendCountMedian(bb, counts).tab();
				appendMetric(bb, mean).tab();
				appendMetric(bb, median(nga)).tab().append(gf/panel.size(), 6);
				for(int value : comparisons){bb.tab().append(value);}
				out.print(bb.nl());
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
	}

	/** Reads wall/user/system seconds from the fixed live/ordinary directories before
	 * opening outputs. Excludes requested genomes, differing hosts and pairs below
	 * ten wall seconds, in that precedence order; this filter does not establish equal
	 * compute configurations or ten seconds of compute. Paths remain live/ordinary even
	 * with custom quality-arm names. Excluded pairs still require valid input files
	 * and finite live/ordinary ratios; unrepresentable ratios reject before output. */
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
					final double value=parseReportDouble(lp, m);
					if(!Double.isFinite(value) || value<0 || (m==0 && value==0)){
						throw new IllegalArgumentException("Invalid native timing: "+base+arm+"="+value);
					}
					data.seconds[g][a*3+m]=value;
				}
				data.hosts[g][a]=new String(singleLine(base+arm+"/host.txt"), java.nio.charset.StandardCharsets.US_ASCII);
				if(!data.hosts[g][a].matches("[A-Za-z0-9.-]+")){throw new IllegalArgumentException("Invalid timing host: "+base+arm);}
			}
			final double live=data.seconds[g][0], ordinary=data.seconds[g][3];
			if(!Double.isFinite(live/ordinary)){
				throw new IllegalArgumentException("Nonfinite timing ratio: genome="+ids.get(g)+
					" arm=live baseline=ordinary wall="+live+" baseline_wall="+ordinary);
			}
			data.reason[g]=timingExclusions.contains(ids.get(g)) ? "excluded_by_design" :
				!data.hosts[g][0].equals(data.hosts[g][1]) ? "different_hosts" :
				Math.min(data.seconds[g][0], data.seconds[g][3])<10 ? "under_10_seconds" : "included";
			if(data.reason[g].equals("included")){data.included++;}
		}
		if(data.included==0){throw new IllegalArgumentException("No comparable timing pairs remain; refusing an empty speed estimate.");}
		return data;
	}

	/**
	 * Writes every timing pair with its exclusion reason, then summarizes included pairs.
	 * Means/medians of paired ratios weight genomes equally; ratio_of_totals instead divides
	 * summed live wall seconds by summed ordinary wall seconds. User/system seconds are
	 * retained in the row table but do not define these wall-time ratios.
	 * Finite sums retain their original arithmetic; overflow uses scaled means.
	 * Large finite timing values use the same formatter as assembly-quality metrics.
	 */
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
				for(double value : seconds){appendMetric(bb.tab(), value);}
				rows.print(appendMetric(bb.tab(), ratio).tab().append(data.reason[g]).nl());
				if(!data.reason[g].equals("included")){continue;}
				live[index]=seconds[0]; ordinary[index]=seconds[3]; ratios[index]=ratio; index++;
				// Keep the original addition order; timingMean handles positive overflow below.
				sumLive+=seconds[0]; sumOrdinary+=seconds[3]; sumRatio+=ratio;
			}
		}finally{if(rows.poisonAndWait()){throw new IllegalStateException("Timing rows write failed.");}}
		assert(index==data.included) : "All included timing pairs must contribute exactly once to summary denominators.";
		// Each sum belongs to its retained array in the same order, before median sorts it.
		final double meanLive=timingMean(live, sumLive), meanOrdinary=timingMean(ordinary, sumOrdinary);
		final double meanRatio=timingMean(ratios, sumRatio);
		final double totalRatio=Double.isFinite(sumLive) && Double.isFinite(sumOrdinary) ?
				sumLive/sumOrdinary : meanLive/meanOrdinary;
		if(!Double.isFinite(totalRatio)){
			throw new IllegalStateException("Timing total ratio must be finite for matched wall times >=10: "+totalRatio);
		}
		final ByteStreamWriter summary=open(output+"/timing_summary.tsv");
		try{
			summary.print("pairs\texcluded\tlive_mean_seconds\tordinary_mean_seconds\tlive_median_seconds\tordinary_median_seconds\tmean_paired_ratio\tmedian_paired_ratio\tratio_of_totals\n");
			bb.clear().append(index).tab().append(ids.size()-index).tab();
			appendMetric(bb, meanLive).tab();
			appendMetric(bb, meanOrdinary).tab();
			appendMetric(bb, median(live)).tab();
			appendMetric(bb, median(ordinary)).tab();
			appendMetric(bb, meanRatio).tab();
			appendMetric(bb, median(ratios)).tab();
			appendMetric(bb, totalRatio);
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

	/**
	 * Requires fixed arm order and each required metric exactly once. All rows need six columns;
	 * other metric names are ignored. Only NGA50 may be '-', represented by zero plus an unset
	 * defined flag. Values must be finite/nonnegative; count/length metrics integral and genome
	 * fraction at most 100 percent. Contig/global/local counts must be below 2^63 after
	 * parsing to double; near-limit integer tokens that round to 2^63 are also rejected.
	 * Counts also require exact integral tokens and are retained as longs. The parsed-double
	 * cap is preserved, so this does not admit every nonnegative signed-long token.
	 * @param path Collected QUAST report.tsv
	 * @param arms Five expected assembly-column names in order
	 * @return Per-arm exact counts and double quality metrics, with undefined NGA50 flagged
	 */
	static ReportRow[] readReport(final String path, final String[] arms){
		assert(arms.length==5) : "The matched reporter compares five ordered arms, with raw first and ordinary last.";
		final ReportRow[] data=new ReportRow[arms.length];
		for(int a=0; a<arms.length; a++){data[a]=new ReportRow();}
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
						final double value=missing ? 0 : (m==1 || m==4 || m==5) ? parseReportDouble(lp, a+1) : lp.parseDouble(a+1);
						if(!Double.isFinite(value) || value<0 || (m!=4 && value!=Math.rint(value)) || (m==4 && value>100)){
							throw new IllegalArgumentException("Invalid QUAST metric: "+METRICS[m]+"="+value);
						}
						if((m==0 || m==2 || m==3) && value>=0x1p63){
							throw new IllegalArgumentException("QUAST count must be below 2^63 after double parsing: "+METRICS[m]+"="+value);
						}
						setMetric(data[a], m, value, missing, lp, a+1, path);
					}
				}
			}
		}finally{if(file.close()){throw new IllegalStateException("QUAST read failed: "+path);}}
		for(int m=0; m<found.length; m++){
			if(!found[m]){throw new IllegalArgumentException("Missing QUAST metric: "+METRICS[m]);}
		}
		return data;
	}

	/**
	 * Parses a quality or timing field without wrapping either default decimal accumulator.
	 * Java mode and fields of at most 18 bytes retain the original parser directly.
	 * Longer default fields use Java conversion only for a decimal token whose whole
	 * or fractional accumulator would overflow, or whose fraction exceeds the shared
	 * multiplier table. Other syntax stays on the original parser route; this does not
	 * enable Java-only forms in default mode. Safe decimal fields retain their original
	 * rounding. The caller validates finiteness, sign and its field-specific domain:
	 * integral NGA50/bases, genome fraction at most 100, or strictly positive wall time.
	 * User/system times may be zero. This does not provide arbitrary-integer precision.
	 * @param lp Borrowed parser positioned on the current quality or timing row
	 * @param term Field index; any missing-NGA50 marker has already been handled
	 * @return Parsed double, subject to the caller's existing domain checks
	 */
	static double parseReportDouble(final LineParser1 lp, final int term){
		assert(lp!=null) : "readReport and readTimings supply a field before validating its numeric domain.";
		lp.setBounds(term);
		final int a=lp.a(), b=lp.b();
		if(Tools.FORCE_JAVA_PARSE_DOUBLE || b-a<=18){return lp.parseDouble(term);}
		final byte[] line=lp.line();
		int i=a, places=0;
		if(line[i]=='-'){i++;}
		long value=0;
		boolean point=false, digits=false, overflow=false;
		for(; i<b; i++){
			final byte c=line[i];
			if(c=='.' && !point){point=true; value=0; continue;}
			if(c<'0' || c>'9'){return lp.parseDouble(term);}
			digits=true;
			if(point){places++;}
			final int digit=c-'0';
			// Never perform an overflowing operation; the flag survives the decimal point.
			if(value>(Long.MAX_VALUE-digit)/10){overflow=true;}
			else{value=value*10+digit;}
		}
		if(!digits){return lp.parseDouble(term);}
		return overflow || places>=ByteBuilder.decimalInvMult.length ? Parse.parseDoubleSlow(line, a, b) : lp.parseDouble(term);
	}

	/** Stores one already validated metric; count tokens additionally require exact integral parsing. */
	static void setMetric(final ReportRow row, final int metric, final double value, final boolean missing,
			final LineParser1 lp, final int term, final String path){
		assert(row!=null && metric>=0 && metric<METRICS.length) : "readReport maps each known QUAST metric to one retained arm.";
		switch(metric){
			case 0: row.contigs=parseCount(lp, term, METRICS[metric], path); break;
			case 1: row.nga50=value; row.nga50Defined=!missing; break;
			case 2: row.global=parseCount(lp, term, METRICS[metric], path); break;
			case 3: row.local=parseCount(lp, term, METRICS[metric], path); break;
			case 4: row.genomeFraction=value; break;
			case 5: row.bases=value; break;
			default: throw new IllegalArgumentException("Unknown QUAST metric index: "+metric);
		}
	}

	/**
	 * Reads an exact nonnegative count without floating arithmetic or unchecked accumulation.
	 * The caller first applies its existing parseDouble mode, finite/integral and below-2^63
	 * checks; this preserves the conservative admission cap and parser-mode syntax. Here
	 * decimal/hex significands are normalized by trailing zero digits and their exponent.
	 * Nonzero fractional residues and values outside long range reject before output.
	 * No strings or big-number objects are allocated on a successful parse.
	 * @param lp Borrowed line parser; its field bounds are set to term
	 * @param term Count field index
	 * @param metric Diagnostic metric name
	 * @param path Diagnostic report path
	 * @return Exact count; negative zero is zero
	 */
	private static long parseCount(final LineParser1 lp, final int term, final String metric, final String path){
		assert(lp!=null) : "readReport supplies the current field after native parser-mode validation.";
		lp.setBounds(term);
		final byte[] line=lp.line();
		int a=lp.a(), b=lp.b();
		while(a<b && (line[a]&255)<=32){a++;}
		while(b>a && (line[b-1]&255)<=32){b--;}
		if(a==b){throw invalidCount(lp, metric, path);}
		final byte suffix=line[b-1];
		if(suffix=='f' || suffix=='F' || suffix=='d' || suffix=='D'){b--;}
		final boolean negative=line[a]=='-';
		if(line[a]=='-' || line[a]=='+'){a++;}
		final boolean hex=b-a>=2 && line[a]=='0' && (line[a+1]=='x' || line[a+1]=='X');
		if(hex){a+=2;}
		int end=a;
		while(end<b && !(hex ? line[end]=='p' || line[end]=='P' : line[end]=='e' || line[end]=='E')){end++;}
		if(a>=end || (hex && end==b)){throw invalidCount(lp, metric, path);}
		int digits=0, fractional=0, trailing=0, first=-1, last=-1, firstDigit=0, lastDigit=0;
		boolean point=false;
		for(int i=a; i<end; i++){
			if(line[i]=='.'){
				if(point){throw invalidCount(lp, metric, path);}
				point=true; continue;
			}
			final int digit=countDigit(line[i], hex);
			if(digit<0){throw invalidCount(lp, metric, path);}
			digits++; if(point){fractional++;}
			if(digit!=0){
				if(first<0){first=i; firstDigit=digit;}
				last=i; lastDigit=digit; trailing=0;
			}else if(first>=0){trailing++;}
		}
		if(digits==0){throw invalidCount(lp, metric, path);}
		long exponent=0;
		if(end<b){
			int e=end+1;
			final boolean minus=e<b && line[e]=='-';
			if(e<b && (line[e]=='-' || line[e]=='+')){e++;}
			if(e==b){throw invalidCount(lp, metric, path);}
			// Beyond this magnitude no nonzero token can cancel enough mantissa digits to fit.
			final long limit=4L*(lp.b()-lp.a())+64;
			for(; e<b; e++){
				final int digit=countDigit(line[e], false);
				if(digit<0){throw invalidCount(lp, metric, path);}
				exponent=Math.min(limit, exponent*10+digit);
			}
			if(minus){exponent=-exponent;}
		}
		if(first<0){return 0;}
		if(negative){throw invalidCount(lp, metric, path);}
		final long scale=exponent+(hex ? 4L : 1L)*(trailing-fractional);
		long value=0;
		if(hex){
			int significant=0;
			for(int i=first; i<=last; i++){if(line[i]!='.'){significant++;}}
			final long bits=4L*(significant-1)+32-Integer.numberOfLeadingZeros(firstDigit);
			if(bits+scale>63 || bits+scale<1 || (scale<0 && -scale>Integer.numberOfTrailingZeros(lastDigit))){
				throw invalidCount(lp, metric, path);
			}
			// At most three low bits remain after removing trailing zero hex digits. Divide
			// during accumulation, so a 65/66-bit significand can still yield an exact long.
			final int cut=scale<0 ? (int)-scale : 0, mask=(1<<cut)-1;
			int carry=0;
			for(int i=first; i<=last; i++){
				if(line[i]=='.'){continue;}
				final int digit=countDigit(line[i], true);
				value=(value<<4)|((long)carry<<(4-cut))|(digit>>>cut);
				carry=digit&mask;
			}
			if(scale>0){value<<=(int)scale;}
		}else{
			// The normalized decimal significand ends in a nonzero digit; a negative scale
			// is therefore fractional, even if the preceding double parse rounded it away.
			if(scale<0 || scale>18){throw invalidCount(lp, metric, path);}
			for(int i=first; i<=last; i++){
				if(line[i]=='.'){continue;}
				final int digit=countDigit(line[i], false);
				if(value>(Long.MAX_VALUE-digit)/10){throw invalidCount(lp, metric, path);}
				value=value*10+digit;
			}
			for(int i=0; i<scale; i++){
				if(value>Long.MAX_VALUE/10){throw invalidCount(lp, metric, path);}
				value*=10;
			}
		}
		return value;
	}

	/** ASCII digit value for the selected radix, or -1 for any other byte. */
	private static int countDigit(final byte c, final boolean hex){
		return c>='0' && c<='9' ? c-'0' : hex && c>='a' && c<='f' ? c-'a'+10 :
			hex && c>='A' && c<='F' ? c-'A'+10 : -1;
	}

	/** Allocates diagnostics only on an invalid count; the parser retains original field bounds. */
	private static IllegalArgumentException invalidCount(final LineParser1 lp, final String metric, final String path){
		return new IllegalArgumentException("QUAST count must be an exact nonnegative integer within long range: "+metric+"="+
			new String(lp.line(), lp.a(), lp.b()-lp.a(), java.nio.charset.StandardCharsets.US_ASCII)+" in "+path);
	}

	/** Sorts disposable nonnegative counts and emits their exact integer or half-integer median. */
	static ByteBuilder appendCountMedian(final ByteBuilder out, final long[] values){
		assert(out!=null && values.length>0) : "Count median writes one nonempty validated genome panel to its caller-owned row.";
		Arrays.sort(values);
		final long low=values[(values.length-1)/2], high=values[values.length/2];
		out.append(low+(high-low)/2);
		if(((low^high)&1)!=0){out.append(".500000");}
		return out;
	}

	/** Maps greater/equal/less values into the explicit three-column summary order. */
	private static int compareSlot(final double value, final double baseline){return value>baseline ? 0 : value==baseline ? 1 : 2;}

	/**
	 * Appends a finite nonnegative quality or timing value. The ordinary six-decimal
	 * path is unchanged below 2^63; at and above that boundary the fast formatter
	 * saturates its long casts, so Java double text preserves the parsed value.
	 * @param out Borrowed row builder
	 * @param value NGA50, bases, timing seconds or a derived mean, median or ratio
	 * @return The same builder for chaining
	 * @throws IllegalArgumentException If a nonfinite or negative value would be published
	 */
	static ByteBuilder appendMetric(final ByteBuilder out, final double value){
		assert(out!=null) : "Report metrics append to the caller-owned output row.";
		if(!Double.isFinite(value) || value<0){throw new IllegalArgumentException("Report metric must be finite and nonnegative: "+value);}
		return value<0x1p63 ? out.append(value, 6) : out.appendSlow(value);
	}

	/**
	 * Returns a finite mean without changing the borrowed array. Finite sums retain
	 * the original sum/length calculation; positive overflow uses max-scaled values.
	 * Zero-only inputs use the finite path. Scaling has ordinary double rounding.
	 * @param values Nonempty finite nonnegative timing values, included paired ratios or NGA50 values
	 * @param sum Caller's left-to-right sum from +0, in this array's original order;
	 * consistency with values is the caller's obligation, not recomputed here
	 * @return Finite nonnegative mean; does not certify decimal formatting
	 * @throws IllegalArgumentException If explicit input preconditions are violated
	 * @throws IllegalStateException If the resulting mean is not finite/nonnegative
	 */
	static double timingMean(final double[] values, final double sum){
		if(values==null || values.length==0){throw new IllegalArgumentException("Timing mean requires a nonempty array.");}
		if(Double.isNaN(sum) || sum<0){throw new IllegalArgumentException("Timing mean requires a nonnegative sum or positive overflow: "+sum);}
		double max=0;
		for(int i=0; i<values.length; i++){
			final double value=values[i];
			if(!Double.isFinite(value) || value<0){
				throw new IllegalArgumentException("Timing mean requires finite nonnegative values: index="+i+" value="+value);
			}
			max=Math.max(max, value);
		}
		final double mean;
		if(Double.isFinite(sum)){mean=sum/values.length;}
		else{
			if(max==0){throw new IllegalArgumentException("Timing mean overflow requires a positive maximum.");}
			double scaledSum=0;
			for(double value : values){scaledSum+=value/max;}
			mean=max*(scaledSum/values.length);
		}
		if(!Double.isFinite(mean) || mean<0){throw new IllegalStateException("Timing mean must be finite and nonnegative: "+mean);}
		return mean;
	}

	/** Sorts a nonempty array of finite nonnegative metrics in place. Even panels
	 * retain sum-then-halving when the central sum is finite; on overflow, halving
	 * the operands first keeps the median finite. Callers supply disposable arrays. */
	static double median(final double[] values){
		assert(values.length>0) : "Median requires at least one genome; missing NGA50 values are explicitly zero, not dropped.";
		Arrays.sort(values);
		final int mid=values.length/2;
		if(values.length%2!=0){return values[mid];}
		final double a=values[mid-1], b=values[mid], sum=a+b;
		return Double.isFinite(sum) ? sum/2 : a*0.5+b*0.5;
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

	/** One arm's metrics, filled during input loading and then retained without mutation. */
	static final class ReportRow{
		/** Exact bounded counts; never rounded through the quality-value doubles. */
		long contigs, global, local;
		/** Existing finite-double quality metrics; undefined NGA50 is represented by zero. */
		double nga50, genomeFraction, bases;
		/** Whether the NGA50 input was defined; false also represents an empty arm. */
		boolean nga50Defined;
	}

	/** Small retained per-genome timing records, separate from assembly quality. */
	private static final class TimingData{
		/** Allocates fixed-size records only after the manifest has been validated. */
		TimingData(final int size){
			assert(size>0) : "Timing data corresponds to the validated nonempty panel.";
			seconds=new double[size][6]; hosts=new String[size][2]; reason=new String[size];
		}
		/** Six seconds values per genome: live wall/user/system, then ordinary wall/user/system. */
		final double[][] seconds;
		/** Literal live/ordinary host names; equality is not proof of identical execution settings. */
		final String[][] hosts;
		/** First applicable exclusion reason, or included. */
		final String[] reason;
		/** Number of pairs retained in timing summary denominators. */
		int included;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed destination for completion diagnostics. */
	private final PrintStream outstream;
	/** Collected panel root and separately created result directory. */
	private final String input, output;
	/** Exact required manifest row count, not a truncation limit. */
	private int expected=153;
	/** Required label on every manifest row; no biological partition verification. */
	private String fold="validation";
	/** Enables fixed-directory live/ordinary timing collection. */
	private boolean timings=false;
	/** Explicit timing-only exclusions, which must name manifest genomes. */
	private final HashSet<String> timingExclusions=new HashSet<String>();
	/** Expected report columns; positions 0 and 4 anchor quality comparisons. */
	private String[] arms={"baseline", "oracle", "nn_three", "nn_four", "ordinary"};

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Required QUAST row labels in setMetric order; only NGA50 may be undefined. */
	private static final String[] METRICS={"# contigs", "NGA50", "# misassemblies", "# local misassemblies", "Genome fraction (%)", "Total length"};
}
