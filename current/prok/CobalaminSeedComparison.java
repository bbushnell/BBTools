package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Collects saved native grades and actual caller workload for paired seed sets.
 * Ratios use scheduled family-window events and all enrolled genomes, never
 * emitted-call counts. Union workload sums the three family counters; a shared
 * physical interval can contribute once per family and is not deduplicated.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSeedComparison {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminSeedComparison.class, false);
		String in=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("out")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}
		}
		require(in!=null && out!=null && expected>0, "Bind complete run manifest, fresh output and expected run count");Shared.setThreads(1);
		try{new CobalaminSeedComparison().run(in, out, expected);}finally{Shared.closeStream(pp.outstream);}
	}
	private void run(String path, String out, int expected){
		final ArrayList<Run> runs=new ArrayList<Run>();final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "run\tseed_minimum\tseed_keys\tmap_k\ttop_n\tid_pass\tgrade\tstderr\ttime", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==9, "Run manifest requires nine explicit fields");final Run r=new Run(p);
				require(names.add(r.name), "Repeated run would duplicate comparison rows");r.work(p.parseString(7));r.timing(p.parseString(8));r.grades(p.parseString(6));runs.add(r);
			}
		}finally{require(!in.close(), "Run manifest read failed");}
		require(runs.size()==expected, "Comparison must include every declared run");RrnaResourceIO.fresh(out, "");
		final ByteStreamWriter writer=RrnaResourceIO.writer(out);final ByteBuilder b=new ByteBuilder();
		try{
			writer.println("run\tseed_minimum\tseed_keys\tmap_k\ttop_n\tid_pass\t"+GRADE_HEADER+"\tscheduled_family_window_events\talignments\talignments_per_family_window\talignments_per_genome\tdiagnostic_wall_seconds\tdiagnostic_cpu_seconds");
			for(Run r:runs){for(Grade g:r.grades){
				final long windows=r.windows[g.index], alignments=r.alignments[g.index];
				b.clear().append(r.name).tab().append(r.minimum).tab().append(r.keys).tab().append(r.mapK).tab().append(r.top).tab().append(r.identity, 4).tab().append(g.raw)
					.tab().append(windows).tab().append(alignments).tab();
				if(windows==0){b.append("NA");}else{b.append((double)alignments/windows, 6);}
				writer.print(b.tab().append((double)alignments/g.genomes, 6).tab().append(r.wall, 3).tab().append(r.cpu, 3).nl());
			}}
		}finally{RrnaResourceIO.finish(writer);}
		System.err.println("COBALAMIN_SEED_COMPARISON_PASS runs="+runs.size()+"; diagnostic wall time is not production marginal cost");
	}
	private static final class Run {
		Run(LineParser1 p){
			name=p.parseString(0);minimum=p.parseInt(1);keys=p.parseInt(2);mapK=p.parseInt(3);top=p.parseInt(4);identity=p.parseDouble(5);
			require(!name.isEmpty() && minimum>=1 && keys>=minimum && mapK>0 && top>0 && identity>0 && identity<=1, "Every comparison setting must be explicit and valid");
		}
		void work(String path){
			final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final boolean[][] seen=new boolean[2][3];
			try{for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);if(p.terms()!=2){continue;}
				for(int i=0; i<3; i++){for(int metric=0; metric<2; metric++){
					if(!p.termEquals(metric==0?WINDOW_LABELS[i]:ALIGNMENT_LABELS[i], 0)){continue;}
					require(!seen[metric][i], "Caller must emit one workload total per B12 subtype");seen[metric][i]=true;
					final long value=Long.parseLong(p.parseString(1).trim());require(value>=0, "Caller workload cannot be negative");
					final long[] values=metric==0?windows:alignments;values[i]=value;values[3]=Math.addExact(values[3], value);
				}}
			}}finally{require(!in.close(), "Caller counter read failed");}
			for(int i=0; i<3; i++){
				require(seen[0][i] && seen[1][i], "Missing enabled-family counters are not zeros");
				require(windows[i]>0 || alignments[i]==0, "An alignment requires a scheduled window");
			}
		}
		void timing(String path){
			final ByteFile in=RrnaResourceIO.open(path);
			try{
				final byte[] row=in.nextLine();require(row!=null && in.nextLine()==null, "One exact GNU time row is required");
				final String[] fields=new String(row, java.nio.charset.StandardCharsets.US_ASCII).split(" ");
				require(fields.length==3, "Timing schema must retain wall,user,system seconds");
				wall=number(fields[0], "wall_seconds=");cpu=number(fields[1], "user_seconds=")+number(fields[2], "system_seconds=");
				require(Double.isFinite(cpu), "CPU sum overflow");
			}finally{require(!in.close(), "Timing read failed");}
		}
		void grades(String path){
			final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
			try{RrnaResourceIO.header(in, GRADE_HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==16 && (p.termEquals("SOURCE", 0) || p.termEquals("FULL_CM", 0)), "Use the two unchanged native truth views");
				int index=-1;for(int i=0; i<SCOPES.length; i++){if(p.termEquals(SCOPES[i], 1)){index=i;}}
				require(index>=0 && seen.add(p.parseString(0)+'\t'+p.parseString(1)), "One grade row per known view/subtype or union");
				final long truth=p.parseLong(3), calls=p.parseLong(4), tp=p.parseLong(5), fp=p.parseLong(6), fn=p.parseLong(7), unscored=p.parseLong(8);
				final int genomes=p.parseInt(2);
				require(genomes==300 && truth>=0 && tp>=0 && fp>=0 && fn>=0 && unscored>=0 && tp+fn==truth && tp+fp+unscored==calls,
					"Native one-to-one counts must conserve the full300-genome denominator");
				grades.add(new Grade(row.clone(), index, genomes));
			}}finally{require(!in.close(), "Native grade read failed");}
			require(grades.size()==8, "Both truth views must contain all three subtypes and the union");
		}
		final String name;final int minimum, keys, mapK, top;final double identity;double wall, cpu;
		final long[] windows=new long[4], alignments=new long[4];final ArrayList<Grade> grades=new ArrayList<Grade>();
	}
	private static double number(String value, String prefix){
		require(value.startsWith(prefix), "Timing field order/units changed");final double n=Double.parseDouble(value.substring(prefix.length()));
		require(Double.isFinite(n) && n>=0, "Timing observations must be finite and nonnegative");return n;
	}
	private static final class Grade {
		Grade(byte[] r, int i, int n){raw=r;index=i;genomes=n;}final byte[] raw;final int index, genomes;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final String[] SCOPES={"RF00174", "RF01482", "RF01689", "union"};
	private static final String[] WINDOW_LABELS={"ncRNA Windows (cobalamin_rf00174):", "ncRNA Windows (cobalamin_rf01482):", "ncRNA Windows (cobalamin_rf01689):"};
	private static final String[] ALIGNMENT_LABELS={"ncRNA Alignments (cobalamin_rf00174):", "ncRNA Alignments (cobalamin_rf01482):", "ncRNA Alignments (cobalamin_rf01689):"};
	private static final String GRADE_HEADER="view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact";
}
