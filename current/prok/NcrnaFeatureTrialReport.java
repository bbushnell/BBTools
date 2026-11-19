package prok;

import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Paired three-seed feature comparisons from existing complete Java summaries.
 * Reads no models, matrices or raw inference scores. Identity has only its one
 * retained start canary: its count-one spread is explicit, never pooled as three.
 * Ratios and RMS differences are descriptive trial results, not caller accuracy.
 * RMS differences retain the precision of the saved summaries; formatting extra
 * decimal places does not recover precision lost by their original formatting.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaFeatureTrialReport{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFeatureTrialReport.class, false);
		try{
			new NcrnaFeatureTrialReport(pp.args).process();
			pp.outstream.println("NCRNA_FEATURE_TRIAL_REPORT_COMPLETE_SAVED_SUMMARIES_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}
	public NcrnaFeatureTrialReport(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("expected")){
				final String[] fields=b.split(",", -1);
				require(fields.length==2, "expected= binds start and stop trial counts");
				expected=new long[]{Long.parseLong(fields[0]), Long.parseLong(fields[1])};
			}else{require(parser.parse(arg, a, b), "Unknown feature-report option: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(input!=null && output!=null && expected!=null && expected[0]>0 && expected[1]>0
			&& parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Require complete in=manifest out=fresh_prefix expected=start_count,stop_count");
		require(Tools.testInputFiles(false, true, input) && Tools.testOutputFiles(false, false, false,
			output+".all.tsv", output+".paired.tsv", output+".spread.tsv"), "Unreadable manifest or existing report outputs");
	}
	public void process(){
		final Cell[][][] cells=new Cell[4][3][2];
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		try(Input manifest=new Input(input)){
			RrnaResourceIO.header(manifest.reader, "variant\tseed\tend\tsummary", input);
			for(byte[] line=manifest.reader.nextLine(); line!=null; line=manifest.reader.nextLine()){
				p.set(line);
				require(p.terms()==4, "Every manifest row binds variant, seed, end and saved summary");
				final int variant=index(p, 0, VARIANTS);
				final int seed=p.parseInt(1)-58;
				final int end=index(p, 2, ENDS);
				require(seed>=0 && seed<3 && (variant<3 || seed==0 && end==0) && cells[variant][seed][end]==null,
					"The complete design is six baseline/length/GC cells each, plus only identity seed58/start");
				final String path=p.parseString(3);
				require(!path.isEmpty() && paths.add(path), "Each scientific cell must bind its own complete saved summary");
				cells[variant][seed][end]=read(path, end);
				count++;
			}
		}
		require(count==19, "Require exactly19saved cells: six baselines, twelve length/GC fits and one identity canary");
		for(int variant=0; variant<4; variant++){
			for(int seed=0; seed<3; seed++){
				for(int end=0; end<2; end++){
					final Cell cell=cells[variant][seed][end];
					if(variant==3 && (seed!=0 || end!=0)){require(cell==null, "Omitted identity fits cannot be invented");continue;}
					require(cell!=null && Arrays.equals(cell.raw[0], cells[0][0][end].raw[0]),
						"Every cell must retain the same end-specific do-nothing population and distances");
				}
			}
		}
		write(cells);
	}
	private Cell read(final String path, final int end){
		final Cell cell=new Cell();
		final LineParser1 p=new LineParser1('\t');
		try(Input file=new Input(path)){
			final byte[] header=file.reader.nextLine();
			require(header!=null, "Saved summary must have its complete column header");
			if(canonicalHeader==null){canonicalHeader=header.clone();}
			else{require(Arrays.equals(header, canonicalHeader), "All compared summaries require the same extended schema");}
			p.set(header);
			final int width=p.terms();
			final ObjectIntMap<String> columns=new ObjectIntMap<String>(String.class);
			for(int i=0; i<width; i++){
				final String key=p.parseString(i);
				require(!key.isEmpty() && !columns.contains(key), "Summary column names must be explicit and unique");
				columns.put(key, i);
			}
			require(columns.get("policy")==0 && columns.get("end")==1 && columns.get("usable_jitter_trials")==2,
				"Summary identity columns must retain the qualified evaluator order");
			final boolean legacyFraction=columns.contains("better_fraction");
			require(legacyFraction!=columns.contains("better_fraction_if_changed"),
				"Summary must declare exactly one legacy or conditional movement-fraction column");
			final int[] positions=new int[METRICS.length];
			for(int i=0; i<positions.length; i++){
				positions[i]=columns.get(i==FRACTION && !legacyFraction ? "better_fraction_if_changed" : METRICS[i]);
				require(positions[i]>=0, "Extended summary is missing metric "+METRICS[i]);
			}
			for(byte[] line=file.reader.nextLine(); line!=null; line=file.reader.nextLine()){
				p.set(line);
				require(p.terms()==width && p.termEquals(ENDS[end], 1) && p.parseLong(2)==expected[end],
					"Saved summary end, complete trial count or column count differs: "+path);
				final int policy=index(p, 0, POLICIES);
				require(cell.raw[policy]==null, "Duplicate policy row in saved summary");
				boolean fractionPresent=false;
				for(int i=0; i<positions.length; i++){
					final int col=positions[i];
					if(i==FRACTION && (p.length(col)==0 || legacyFraction && p.termEquals("NA", col))){continue;}
					require(p.length(col)>0, "A measured saved metric cannot be empty");
					p.setBounds(col);
					final double value=Parse.parseDoubleSlow(line, p.a(), p.b());
					require(Double.isFinite(value), "A measured saved metric must be finite");
					cell.values[policy][i]=value;
					if(i==FRACTION){fractionPresent=true;}
				}
				final double[] v=cell.values[policy];
				for(final int i:new int[]{EXACT, NEAR, BETTER, WORSE, TIED}){
					require(v[i]>=0 && v[i]<=expected[end] && v[i]==Math.rint(v[i]), "Saved count is outside its complete integer cohort");
				}
				require(v[EXACT]<=v[NEAR] && v[BETTER]+v[WORSE]+v[TIED]==expected[end]
					&& Math.abs(v[EXACT_PCT]-100*v[EXACT]/expected[end])<=0.000001
					&& Math.abs(v[NEAR_PCT]-100*v[NEAR]/expected[end])<=0.000001
					&& v[RMS]>=0, "Saved percentages and movement outcomes must conserve their declared trials");
				final double changed=v[BETTER]+v[WORSE];
				require(changed==0 ? !fractionPresent : fractionPresent && Math.abs(v[FRACTION]-v[BETTER]/changed)<=0.000001,
					"A better fraction exists exactly when better+worse is positive, and must match those paired counts");
				if(fractionPresent){cell.raw[policy]=line.clone();}
				else{
					p.setBounds(positions[FRACTION]);
					cell.raw[policy]=new ByteBuilder().append(line, 0, p.a()).append(line, p.b(), line.length-p.b()).toBytes();
				}
			}
		}
		for(final byte[] row:cell.raw){require(row!=null, "All three saved policies are required, including do-nothing");}
		return cell;
	}
	private void write(final Cell[][][] cells){
		final ByteBuilder b=new ByteBuilder();
		try(Output all=new Output(output+".all.tsv"); Output paired=new Output(output+".paired.tsv");
				Output spread=new Output(output+".spread.tsv")){
			final LineParser1 header=new LineParser1('\t').set(canonicalHeader);
			b.clear().append("variant\tseed");
			for(int i=0; i<header.terms(); i++){
				b.tab().append(header.termEquals("better_fraction", i) ? "better_fraction_if_changed" : header.parseString(i));
			}
			all.writer.print(b.nl());
			paired.writer.println("variant\tseed\tend\tpolicy\ttrials\texact\tbaseline_exact\texact_delta\texact_pct\texact_delta_pp\twithin3\tbaseline_within3\twithin3_delta\twithin3_delta_pp\trms_after_bp\tbaseline_rms_after_bp\trms_delta_bp\tbetter\tbaseline_better\tworse\tbaseline_worse\ttied\tbaseline_tied\tbetter_fraction_if_changed\tbaseline_better_fraction_if_changed");
			spread.writer.println("variant\tend\tpolicy\tseeds\texact_pct_mean\texact_delta_pp_mean\texact_delta_pp_min\texact_delta_pp_max\trms_delta_bp_mean\trms_delta_bp_min\trms_delta_bp_max\tqualification");
			for(int variant=0; variant<4; variant++){
				for(int end=0; end<2; end++){
					for(int seed=0; seed<3; seed++){
						final Cell cell=cells[variant][seed][end];
						if(cell==null){continue;}
						for(int policy=0; policy<3; policy++){
							all.writer.print(b.clear().append(VARIANTS[variant]).tab().append(seed+58).tab().append(cell.raw[policy]).nl());
							if(variant==0 || policy==0){continue;}
							final double[] a=cell.values[policy];
							final double[] base=cells[0][seed][end].values[policy];
							b.clear().append(VARIANTS[variant]).tab().append(seed+58).tab().append(ENDS[end]).tab().append(POLICIES[policy]).tab().append(expected[end]);
							b.tab().append((long)a[EXACT]).tab().append((long)base[EXACT]).tab().append((long)(a[EXACT]-base[EXACT]));
							b.tab().append(a[EXACT_PCT], 8).tab().append(100*(a[EXACT]-base[EXACT])/expected[end], 8);
							b.tab().append((long)a[NEAR]).tab().append((long)base[NEAR]).tab().append((long)(a[NEAR]-base[NEAR]));
							b.tab().append(100*(a[NEAR]-base[NEAR])/expected[end], 8);
							b.tab().append(a[RMS], 8).tab().append(base[RMS], 8).tab().append(a[RMS]-base[RMS], 8);
							for(final int i:new int[]{BETTER, WORSE, TIED}){b.tab().append((long)a[i]).tab().append((long)base[i]);}
							fraction(b.tab(), a);
							fraction(b.tab(), base);
							paired.writer.print(b.nl());
						}
					}
					if(variant==0 || variant==3 && end==1){continue;}
					for(int policy=1; policy<3; policy++){
						int n=0;
						double exact=0, deltaExact=0, deltaRms=0;
						double minExact=0, maxExact=0, minRms=0, maxRms=0;
						for(int seed=0; seed<3; seed++){
							final Cell cell=cells[variant][seed][end];
							if(cell==null){continue;}
							final double[] a=cell.values[policy];
							final double[] base=cells[0][seed][end].values[policy];
							final double de=100*(a[EXACT]-base[EXACT])/expected[end];
							final double dr=a[RMS]-base[RMS];
							if(n==0){
								minExact=maxExact=de;
								minRms=maxRms=dr;
							}
							n++;
							exact+=100*a[EXACT]/expected[end];
							deltaExact+=de;
							deltaRms+=dr;
							minExact=Math.min(minExact, de);
							maxExact=Math.max(maxExact, de);
							minRms=Math.min(minRms, dr);
							maxRms=Math.max(maxRms, dr);
						}
						assert(n==(variant==3 ? 1 : 3)) : "The validated design fixes the number of paired seeds, never inferred missingness";
						require(Double.isFinite(exact) && Double.isFinite(deltaExact) && Double.isFinite(deltaRms),
							"Computed paired-seed totals must be finite before their means are written");
						spread.writer.print(b.clear().append(VARIANTS[variant]).tab().append(ENDS[end]).tab().append(POLICIES[policy]).tab().append(n)
							.tab().append(exact/n, 8).tab().append(deltaExact/n, 8).tab().append(minExact, 8).tab().append(maxExact, 8)
							.tab().append(deltaRms/n, 8).tab().append(minRms, 8).tab().append(maxRms, 8)
							.tab().append(variant==3 ? "ONE_CANARY_AFFINE_STANDARDIZATION" : "PAIRED_SEEDS_CONDITIONAL_TRIALS_NO_CALLER").nl());
					}
				}
			}
		}
	}
	/** No-change cohorts have no movement fraction; their counts retain the presence condition. */
	private static void fraction(final ByteBuilder b, final double[] values){
		if(values[BETTER]+values[WORSE]>0){b.append(values[FRACTION], 8);}
	}
	private static int index(final LineParser1 p, final int field, final String[] choices){
		for(int i=0; i<choices.length; i++){if(p.termEquals(choices[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown manifest/policy identity in field "+field);
	}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private static final class Cell{
		final byte[][] raw=new byte[3][];
		final double[][] values=new double[3][METRICS.length];
	}
	private static final String[] VARIANTS={"baseline","length","gc","identity"};
	private static final String[] ENDS={"5prime","3prime"};
	private static final String[] POLICIES={"do_nothing","margin0","margin0.1"};
	private static final String[] METRICS={"exact","exact_pct","within3","within3_pct","rms_after_bp","better","worse","tied","better_fraction"};
	private static final int EXACT=0, EXACT_PCT=1, NEAR=2, NEAR_PCT=3, RMS=4, BETTER=5, WORSE=6, TIED=7, FRACTION=8;
	private String input, output;
	private long[] expected;
	private byte[] canonicalHeader;
}
