package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.function.Consumer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Summarizes four saved CallGenes cost reports without rerunning callers.
 * Process CPU, worker CPU and selected startup spans remain distinct populations.
 * Signed per-input deltas are observations, not attribution of all process cost.
 * @author Brian Bushnell, Ganyu
 */
public final class CallGenesProductionCosts {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CallGenesProductionCosts.class, false);
		String in=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else{throw new IllegalArgumentException(key);}
		}
		require(in!=null && out!=null && expected>0, "Bind the four saved runs and their input denominator");Shared.setThreads(1);
		try{final CallGenesProductionCosts x=new CallGenesProductionCosts();x.load(in, expected);x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void load(String path, int expected){
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		read(path, RUN_HEADER, p->{
			final String name=p.parseString(0);final int threads=p.parseInt(1);final boolean enabled=bool(p, 2);
			require(!name.isEmpty() && names.add(name) && (threads==1 || threads==16), "Each off/on t1/t16 arm must be unique");
			final int index=(threads==1?0:2)+(enabled?1:0);require(runs[index]==null, "Duplicate thread/enable arm");
			final Run r=new Run(name, threads, enabled, expected);runs[index]=r;
			final String prefix=p.parseString(3);r.inputs(prefix+".inputs.tsv");r.families(prefix+".families.tsv");r.startup(prefix+".startup.tsv");r.time(p.parseString(4));
		});
		for(Run r:runs){require(r!=null, "All four arms are required before comparing process and per-input costs");}
		for(int i=1; i<4; i++){for(int j=0; j<expected; j++){
			final Input a=runs[0].inputs[j], b=runs[i].inputs[j];
			require(a.path.equals(b.path) && a.reads==b.reads && a.bases==b.bases,
				"Input ordinals, paths and processed bases/reads must agree across every arm");
		}}
		for(int i=0; i<4; i+=2){
			final Run off=runs[i], on=runs[i+1];
			for(String name:off.familyOrder){require(on.byFamily.contains(name), "Enabling B12 cannot silently drop another diagnostic family");}
			for(String name:on.familyOrder){require(off.byFamily.contains(name) || b12(name)>=0, "Only B12 may add diagnostic family slots");}
		}
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".runs.tsv", ".startup.tsv", ".inputs.tsv", ".deltas.tsv", ".families.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".runs.tsv"), startup=RrnaResourceIO.writer(out+".startup.tsv");
		final ByteStreamWriter inputs=RrnaResourceIO.writer(out+".inputs.tsv"), delta=RrnaResourceIO.writer(out+".deltas.tsv"), families=RrnaResourceIO.writer(out+".families.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			summary.println("run\tthreads\tb12\tinputs\twall_seconds\tuser_seconds\tsystem_seconds\tprocess_cpu_seconds\tworker_cpu_seconds_sum\tselected_startup_process_cpu_seconds\tselected_startup_elapsed_seconds\tlazy_init_cpu_seconds_sum\tB12_lazy_init_cpu_seconds\tB12_alignments\tB12_family_window_events");
			startup.println("run\tthreads\tb12\tstage\tstatus\telapsed_seconds\tprocess_cpu_seconds");
			inputs.println("run\tthreads\tb12\tinput_index\tinput\tworkers\treads\tbases\tworker_cpu_seconds\tworker_span_elapsed_seconds_including_waits\tprocess_read_inner_elapsed_seconds_sum\tlazy_init_cpu_seconds\tB12_lazy_init_cpu_seconds\tB12_alignments\tB12_family_window_events");
			delta.println("threads\tinput_index\tinput\treads\tbases\tworker_cpu_off_seconds\tworker_cpu_on_seconds\tworker_cpu_delta_seconds\tspan_off_seconds_including_waits\tspan_on_seconds_including_waits\tlazy_init_off_seconds\tlazy_init_on_seconds\tB12_lazy_init_on_seconds\tB12_alignments_on\tB12_family_window_events_on");
			families.println("run\tthreads\tb12\tfamily\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\tfamily_window_events\tlazy_init_cpu_seconds");
			for(Run r:runs){
				long worker=0, init=0, b12init=0, aligns=0, windows=0;
				for(Input x:r.inputs){worker=Math.addExact(worker, x.cpu);init=Math.addExact(init, x.init);b12init=Math.addExact(b12init, x.b12init);
					aligns=Math.addExact(aligns, x.alignments);windows=Math.addExact(windows, x.windows);}
				long startupCpu=0, startupWall=0;
				for(Stage s:r.stages){
					startupCpu=Math.addExact(startupCpu, s.cpu);startupWall=Math.addExact(startupWall, s.elapsed);
					prefix(b, r).tab().append(s.name).tab().append(s.status);seconds(b, s.elapsed);seconds(b, s.cpu);startup.print(b.nl());
				}
				prefix(b, r).tab().append(r.inputs.length).tab().append(r.wall, 6).tab().append(r.user, 6).tab().append(r.system, 6).tab().append(r.user+r.system, 6);
				seconds(b, worker);seconds(b, startupCpu);seconds(b, startupWall);seconds(b, init);seconds(b, b12init);summary.print(b.tab().append(aligns).tab().append(windows).nl());
				for(int i=0; i<r.inputs.length; i++){
					final Input x=r.inputs[i];prefix(b, r).tab().append(i).tab().append(x.path).tab().append(x.workers).tab().append(x.reads).tab().append(x.bases);
					seconds(b, x.cpu);seconds(b, x.span);seconds(b, x.calling);seconds(b, x.init);seconds(b, x.b12init);inputs.print(b.tab().append(x.alignments).tab().append(x.windows).nl());
				}
				Collections.sort(r.familyOrder);
				for(String name:r.familyOrder){
					final Family f=r.byFamily.get(name);prefix(b, r).tab().append(name);
					for(int i=0; i<6; i++){b.tab().append(f.total[i]);}seconds(b, f.total[6]);families.print(b.nl());
				}
			}
			for(int pair=0; pair<4; pair+=2){
				final Run off=runs[pair], on=runs[pair+1];
				for(int i=0; i<off.inputs.length; i++){
					final Input a=off.inputs[i], z=on.inputs[i];b.clear().append(off.threads).tab().append(i).tab().append(a.path).tab().append(a.reads).tab().append(a.bases);
					seconds(b, a.cpu);seconds(b, z.cpu);seconds(b, Math.subtractExact(z.cpu, a.cpu));seconds(b, a.span);seconds(b, z.span);
					seconds(b, a.init);seconds(b, z.init);seconds(b, z.b12init);delta.print(b.tab().append(z.alignments).tab().append(z.windows).nl());
				}
			}
		}finally{RrnaResourceIO.finish(summary, startup, inputs, delta, families);}
		System.err.println("CALLGENES_PRODUCTION_COSTS_PASS runs=4 inputs="+runs[0].inputs.length+"; process, worker and selected startup clocks remain distinct");
	}
	private static ByteBuilder prefix(ByteBuilder b, Run r){return b.clear().append(r.name).tab().append(r.threads).tab().append(r.enabled);}
	private static void seconds(ByteBuilder b, long ns){b.tab().append(ns/1e9, 9);}
	private static final class Run {
		Run(String n, int t, boolean e, int count){name=n;threads=t;enabled=e;inputs=new Input[count];}
		void inputs(String path){
			read(path, INPUT_HEADER, p->{
				final int index=p.parseInt(0);require(index>=0 && index<inputs.length && inputs[index]==null, "Each original input ordinal must appear once");
				final Input x=new Input(p.parseString(1));x.workers=p.parseInt(2);x.reads=p.parseLong(3);x.bases=p.parseLong(4);
				x.cpu=p.parseLong(5);x.span=p.parseLong(6);x.calling=p.parseLong(7);
				require(!x.path.isEmpty() && x.workers>0 && x.reads>=0 && x.bases>=0 && x.cpu>=0 && x.span>=0 && x.calling>=0, "Invalid completed worker measurement");inputs[index]=x;
			});
			for(Input x:inputs){require(x!=null, "Missing input cannot disappear from timing denominators");}
		}
		void families(String path){
			read(path, FAMILY_HEADER, p->{
				final int index=p.parseInt(0);require(index>=0 && index<inputs.length, "Family counters need an original input ordinal");
				final Input x=inputs[index];require(p.termEquals(x.path, 1), "Family counters must match the exact input path");
				final String name=p.parseString(2);require(!name.isEmpty() && x.seen.add(name), "Duplicate family row would inflate work");
				Family f=byFamily.get(name);if(f==null){f=new Family();byFamily.put(name, f);familyOrder.add(name);}
				final int member=b12(name);if(member>=0){require(enabled, "Disabled B12 cannot report B12 family work");x.b12Mask|=1<<member;}
				for(int i=0; i<7; i++){
					require(!p.termEquals("NA", i+3), "Required CPU diagnostic was disabled; missing measurements are not zeros");
					final long value=p.parseLong(i+3);require(value>=0, "Work counters cannot be negative");f.total[i]=Math.addExact(f.total[i], value);
					if(i==6){x.init=Math.addExact(x.init, value);if(member>=0){x.b12init=Math.addExact(x.b12init, value);}}
					else if(member>=0 && i==0){x.alignments=Math.addExact(x.alignments, value);}else if(member>=0 && i==5){x.windows=Math.addExact(x.windows, value);}
				}
			});
			require(!familyOrder.isEmpty(), "Production controls require default RNA family diagnostics");
			for(Input x:inputs){
				require(x.seen.size()==familyOrder.size() && x.b12Mask==(enabled?7:0), "Every input must report every enabled family, including zero work");
				for(String name:familyOrder){require(x.seen.contains(name), "Input family roster differs within a run");}
			}
		}
		void startup(String path){
			read(path, STARTUP_HEADER, p->{
				final int index=p.termEquals("B12_RESOURCES", 0)?0:p.termEquals("SHARED_SEED_INDEX", 0)?1:-1;
				require(index>=0 && stages[index]==null, "Each selected startup span must be measured once");
				final boolean skipped=index==0 && !enabled;final long elapsed=p.parseLong(2), cpu=p.parseLong(3);
				require(p.termEquals(skipped?"SKIPPED_DISABLED":"COMPLETE", 1) && elapsed>=0 && cpu>=0 && (!skipped || elapsed==0 && cpu==0), "Startup status must reflect the enabled configuration");
				stages[index]=new Stage(p.parseString(0), p.parseString(1), elapsed, cpu);
			});
			for(Stage s:stages){require(s!=null, "Missing startup span cannot be called zero-cost");}
		}
		void time(String path){
			final ByteFile in=RrnaResourceIO.open(path);
			try{
				final byte[] row=in.nextLine();require(row!=null && in.nextLine()==null, "One GNU timing row required");
				final String[] fields=new String(row, java.nio.charset.StandardCharsets.US_ASCII).split(" ");require(fields.length==3, "Preserve wall/user/system units");
				wall=number(fields[0], "wall_seconds=");user=number(fields[1], "user_seconds=");system=number(fields[2], "system_seconds=");
				require(Double.isFinite(user+system), "Process CPU total overflow");
			}finally{require(!in.close(), "Timing input failed");}
		}
		final String name;final int threads;final boolean enabled;final Input[] inputs;final Stage[] stages=new Stage[2];
		final ObjectMap<String,Family> byFamily=new ObjectMap<String,Family>(String.class, Family.class);final ArrayList<String> familyOrder=new ArrayList<String>();
		double wall, user, system;
	}
	private static final class Input {
		Input(String p){path=p;}final String path;int workers, b12Mask;long reads, bases, cpu, span, calling, init, b12init, alignments, windows;
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
	}
	private static final class Family {final long[] total=new long[7];}
	private static final class Stage {Stage(String n, String s, long e, long c){name=n;status=s;elapsed=e;cpu=c;}final String name, status;final long elapsed, cpu;}
	private static double number(String text, String prefix){require(text.startsWith(prefix), "GNU timing field changed");final double n=Double.parseDouble(text.substring(prefix.length()));require(Double.isFinite(n) && n>=0, "Timing must be finite and nonnegative");return n;}
	private static int b12(String name){for(int i=0; i<MEMBERS.length; i++){if(MEMBERS[i].equals(name)){return i;}}require(!name.startsWith("cobalamin"), "Unknown B12 diagnostic family");return -1;}
	private static boolean bool(LineParser1 p, int col){require(p.termEquals("true", col) || p.termEquals("false", col), "Use explicit true/false enable state");return p.termEquals("true", col);}
	private static void read(String path, String header, Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.UTF_8));final int terms=p.terms();boolean seen=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			if(row.length>0 && row[0]=='#'){continue;}
			if(!seen){require(header.equals(new String(row, java.nio.charset.StandardCharsets.UTF_8)), "Changed header: "+path);seen=true;continue;}
			p.set(row);require(p.terms()==terms, "Changed row schema: "+path);action.accept(p);
		}}finally{require(!in.close(), "Read failed: "+path);}require(seen, "Missing header: "+path);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final Run[] runs=new Run[4];
	private static final String[] MEMBERS={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	public static final String RUN_HEADER="run\tthreads\tb12\tcost_prefix\ttime";
	public static final String INPUT_HEADER="input_index\tinput\tworkers\treads\tbases\tworker_cpu_ns\tworker_span_elapsed_ns_including_waits\tprocess_read_inner_elapsed_ns_sum";
	public static final String FAMILY_HEADER="input_index\tinput\tfamily\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\twindows\tlazy_init_cpu_ns";
	public static final String STARTUP_HEADER="stage\tstatus\telapsed_ns\tprocess_cpu_ns";
}
