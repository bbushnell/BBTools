package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import parse.Parse;
import shared.Shared;

/**
 * In-memory scanner timing carrier for the AssemblyStats2 record-length guard.
 * Pilot mode runs baseline-only sizing and saves per-shape reps. Measure mode
 * reuses that exact plan in both phases so fixed never recalibrates on its own.
 * @author Jean
 */
public final class AssemblyStats2RecordLengthTiming{

	/** phase=baseline|fixed|control|candidate mode=pilot|measure outdir=NEW
	 * [plan=baseline timing-plan.tsv] [shape=fasta-long-a] [rounds=1]. */
	public static void main(final String[] args) throws Exception{
		String phase=null, mode=null, directory=null, plan=null, shape=null;
		int rounds=MEASURE_ROUNDS;
		boolean roundsSet=false;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("plan")){plan=pair[1];}
			else if(pair[0].equals("shape")){shape=pair[1];}
			else if(pair[0].equals("rounds")){
				rounds=Integer.parseInt(pair[1]);
				roundsSet=true;
			}
			else{throw new IllegalArgumentException("Unknown record-length timing option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase) ||
			"control".equals(phase) || "candidate".equals(phase),
			"Set phase=baseline|fixed|control|candidate.");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require(!"pilot".equals(mode) || "baseline".equals(phase), "Only baseline may run pilot mode.");
		require("pilot".equals(mode) == (plan==null), "Pilot writes a plan; measure requires plan=<baseline timing-plan.tsv>.");
		require(!"pilot".equals(mode) || shape==null, "shape selection is only supported in measure mode.");
		require(!"pilot".equals(mode) || !roundsSet, "rounds selection is only supported in measure mode.");
		require(rounds>0 && rounds<=MEASURE_ROUNDS, "rounds must be in 1.."+MEASURE_ROUNDS+".");
		require(shape==null || hasShape(shape), "Unknown shape: "+shape);
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() &&
			clock.getCurrentThreadCpuTime()>=0, "Current-thread CPU time must be supported and enabled.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats2recordlengthtiming.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyStats2RecordLengthTiming timing=
			new AssemblyStats2RecordLengthTiming(root, phase, mode, shape, rounds, clock);
		timing.run(plan==null ? null : Paths.get(plan).toAbsolutePath().normalize());
		System.out.print("ASSEMBLY_STATS2_RECORD_LENGTH_TIMING_"+(timing.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" mode="+mode+" shapes="+timing.selectedShapes.size()+
			" samples="+timing.samples+" accepted="+timing.accepted+
			" failures="+timing.failures+"\n");
		if(timing.failures>0){throw new AssertionError("Inspect "+root);}
	}

	/** Creates output files and materializes all four byte-array inputs once. */
	private AssemblyStats2RecordLengthTiming(final Path root_, final String phase_,
			final String mode_, final String selectedShapeName_,
			final int measureRounds_, final ThreadMXBean clock_) throws Exception{
		root=root_;
		phase=phase_;
		mode=mode_;
		selectedShapeName=selectedShapeName_;
		measureRounds=measureRounds_;
		clock=clock_;
		Files.write(root.resolve("timing-shapes.tsv"), SHAPES_HEADER.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("samples.tsv"), SAMPLE_HEADER.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		for(final ShapeSpec spec : SHAPES){shapes.put(spec.name, new Shape(spec));}
		if(selectedShapeName==null){
			selectedShapes.addAll(shapes.values());
		}else{
			final Shape selected=shapes.get(selectedShapeName);
			require(selected!=null, "Unknown shape: "+selectedShapeName);
			selectedShapes.add(selected);
		}
		Files.write(root.resolve("config.tsv"), config().getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Runs either baseline pilot sizing or one exact saved measurement plan. */
	private void run(final Path plan) throws Exception{
		writeShapes();
		if("pilot".equals(mode)){
			Files.write(root.resolve("timing-plan.tsv"), PLAN_HEADER.getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE_NEW);
			for(final Shape shape : selectedShapes){pilot(shape);}
		}else{
			final Map<String, PlanRow> rows=readPlan(plan);
			for(final Shape shape : selectedShapes){
				final PlanRow row=rows.get(shape.spec.name);
				require(row!=null, "Missing shape from saved plan: "+shape.spec.name);
				shape.checkPlan(row);
				for(int round=1; round<=measureRounds; round++){
					final Result result=sample(shape, row.reps, warmup(row.reps), round);
					if(result.cpuNanos>=MIN_ACCEPTED_CPU_NANOS){accepted++;}
					else{failures++;}
				}
			}
		}
	}

	/** Finds baseline reps for one shape, targeting twelve scanner CPU seconds. */
	private void pilot(final Shape shape) throws Exception{
		for(int reps=1; reps<=MAX_PILOT_REPS; reps*=2){
			final Result result=sample(shape, reps, warmup(reps), 0);
			if(result.cpuNanos>=MIN_PILOT_CPU_NANOS && result.wallNanos>=MIN_PILOT_WALL_NANOS){
				final int measuredReps=repsForTarget(reps, result.cpuNanos);
				Files.write(root.resolve("timing-plan.tsv"),
					(shape.spec.name+"\t"+measuredReps+"\t"+warmup(measuredReps)+"\t"+
					shape.spec.records+"\t"+shape.spec.bases+"\t"+shape.data.length+"\t"+
					shape.sha80+"\t"+result.cpuNanos+"\t"+result.wallNanos+"\n").getBytes(StandardCharsets.UTF_8),
					StandardOpenOption.APPEND);
				return;
			}
		}
		throw new AssertionError("Pilot never reached 250 ms CPU and wall for "+shape.spec.name);
	}

	/** Runs warmup reps, then scanner-only timed reps with setup outside each interval. */
	private Result sample(final Shape shape, final int reps, final int warmup, final int round) throws Exception{
		check(warmup<=MAX_MEASURE_WARMUP_REPS, "Warmup exceeded the per-sample cap: "+warmup);
		long checksum=0;
		for(int i=0; i<warmup; i++){checksum+=single(shape, false).checksum;}
		long cpu=0, wall=0;
		for(int i=0; i<reps; i++){
			final RunResult result=single(shape, true);
			checksum+=result.checksum;
			cpu+=result.cpuNanos;
			wall+=result.wallNanos;
		}
		SINK=checksum;
		final boolean acceptedSample=("measure".equals(mode) && cpu>=MIN_ACCEPTED_CPU_NANOS);
		final Result result=new Result(shape.spec.name, reps, warmup, round, cpu, wall, checksum, acceptedSample);
		Files.write(root.resolve("samples.tsv"), result.toString().getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
		samples++;
		return result;
	}

	/** Runs one count method call; scanner time starts after setup and stops before checks. */
	private RunResult single(final Shape shape, final boolean measured) throws Exception{
		final StaticState state=new StaticState();
		long cpu=0, wall=0;
		long[] counts=null;
		try{
			final AssemblyStats2 stats=new AssemblyStats2(SCANNER_ARGS);
			final InputStream input=new ByteArrayInputStream(shape.data);
			if(measured){
				final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
				counts=shape.spec.fastq ? stats.countFastq(input, null) : stats.countFasta(input, null);
				wall=System.nanoTime()-wallStart;
				cpu=clock.getCurrentThreadCpuTime()-cpuStart;
			}else{
				counts=shape.spec.fastq ? stats.countFastq(input, null) : stats.countFasta(input, null);
			}
		}finally{
			state.restore();
		}
		shape.checkCounts(counts);
		long checksum=0;
		for(final long value : counts){checksum+=value;}
		return new RunResult(cpu, wall, checksum);
	}

	/** Writes immutable input-shape oracles. */
	private void writeShapes() throws Exception{
		for(final Shape shape : shapes.values()){
			Files.write(root.resolve("timing-shapes.tsv"),
				(shape.spec.name+"\t"+(shape.spec.fastq ? "fastq" : "fasta")+"\t"+
				shape.spec.records+"\t"+shape.spec.bases+"\t"+shape.data.length+"\t"+
				shape.sha80+"\t"+Arrays.toString(shape.spec.expected)+"\n").getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.APPEND);
		}
	}

	/** Reads the successful baseline pilot's exact per-shape reps. */
	private Map<String, PlanRow> readPlan(final Path plan) throws Exception{
		final LinkedHashMap<String, PlanRow> rows=new LinkedHashMap<String, PlanRow>();
		final java.util.List<String> lines=Files.readAllLines(plan, StandardCharsets.UTF_8);
		require(!lines.isEmpty() && lines.get(0).equals(PLAN_HEADER.trim()), "Unrecognized timing plan header.");
		for(int i=1; i<lines.size(); i++){
			final PlanRow row=new PlanRow(lines.get(i).split("\t", -1));
			require(!rows.containsKey(row.shape), "Duplicate timing plan row: "+row.shape);
			rows.put(row.shape, row);
		}
		require(rows.size()==SHAPES.length, "Timing plan must have exactly "+SHAPES.length+" shapes.");
		return rows;
	}

	/** Returns an explicit configuration and bound record for controller artifacts. */
	private String config(){
		return "version\t1\n"+
			"selected_shape\t"+(selectedShapeName==null ? "all" : selectedShapeName)+"\n"+
			"selected_shapes\t"+selectedShapes.size()+"\n"+
			"target_cpu_nanos\t"+TARGET_CPU_NANOS+"\n"+
			"min_accepted_cpu_nanos\t"+MIN_ACCEPTED_CPU_NANOS+"\n"+
			"min_pilot_cpu_nanos\t"+MIN_PILOT_CPU_NANOS+"\n"+
			"min_pilot_wall_nanos\t"+MIN_PILOT_WALL_NANOS+"\n"+
			"max_pilot_measured_reps_per_shape\t"+MAX_PILOT_TOTAL_REPS_PER_SHAPE+"\n"+
			"max_pilot_warmup_reps_per_shape\t"+MAX_PILOT_WARMUP_REPS_PER_SHAPE+"\n"+
			"max_pilot_measured_reps_all_shapes\t"+(MAX_PILOT_TOTAL_REPS_PER_SHAPE*SHAPES.length)+"\n"+
			"max_pilot_warmup_reps_all_shapes\t"+(MAX_PILOT_WARMUP_REPS_PER_SHAPE*SHAPES.length)+"\n"+
			"measure_rounds\t"+measureRounds+"\n"+
			"max_measure_warmup_reps_per_phase\t"+
			(measureRounds*selectedShapes.size()*MAX_MEASURE_WARMUP_REPS)+"\n";
	}

	/** Scales a pilot to a twelve-second target and enforces the explicit cap. */
	private static int repsForTarget(final int reps, final long cpuNanos){
		require(cpuNanos>0, "Cannot scale a zero-CPU pilot.");
		final long scaled=(TARGET_CPU_NANOS*reps+cpuNanos-1)/cpuNanos;
		final long clamped=Math.max(MIN_MEASURE_REPS, scaled);
		require(clamped<=MAX_MEASURE_REPS, "Pilot selected "+clamped+" reps beyond cap "+MAX_MEASURE_REPS);
		return (int)clamped;
	}

	/** Chooses bounded untimed reps without changing measured work selection. */
	private static int warmup(final int reps){
		return Math.max(MIN_WARMUP_REPS, Math.min(MAX_MEASURE_WARMUP_REPS, (reps+15)/16));
	}

	/** Returns true when a shape selector names an existing immutable input. */
	private static boolean hasShape(final String name){
		for(final ShapeSpec spec : SHAPES){
			if(spec.name.equals(name)){return true;}
		}
		return false;
	}

	/** Records a mismatch without hiding subsequent shape evidence. */
	private void check(final boolean correct, final String message){
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Rejects setup errors before a misleading timing file can be used. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Saves and restores mutable statics touched by AssemblyStats2 construction and scanning. */
	private static final class StaticState{
		StaticState() throws Exception{
			for(final StaticRef ref : STATIC_FIELDS){
				final Field field=ref.clazz.getDeclaredField(ref.name);
				field.setAccessible(true);
				values.put(field, field.get(null));
			}
		}
		void restore() throws Exception{
			for(final Map.Entry<Field, Object> entry : values.entrySet()){
				entry.getKey().set(null, entry.getValue());
			}
		}
		private final LinkedHashMap<Field, Object> values=new LinkedHashMap<Field, Object>();
	}

	/** A mutable static field whose value must be restored after direct constructor calls. */
	private static final class StaticRef{
		StaticRef(final Class<?> clazz_, final String name_){clazz=clazz_; name=name_;}
		final Class<?> clazz;
		final String name;
	}

	/** One immutable timing input shape and its expected scanner counts. */
	private static final class ShapeSpec{
		ShapeSpec(final String name_, final boolean fastq_, final int records_, final int length_,
				final int a_, final int c_, final int n_, final long[] expected_){
			name=name_;
			fastq=fastq_;
			records=records_;
			length=length_;
			a=a_;
			c=c_;
			n=n_;
			expected=expected_;
			bases=(long)records*length;
		}
		final String name;
		final boolean fastq;
		final int records, length, a, c, n;
		final long bases;
		final long[] expected;
	}

	/** Materialized shape bytes plus their sha80. */
	private static final class Shape{
		Shape(final ShapeSpec spec_) throws Exception{
			spec=spec_;
			data=make(spec);
			sha80=sha80(data);
			checkBytes();
		}
		void checkCounts(final long[] counts){
			require(Arrays.equals(counts, spec.expected), spec.name+" counts differed: "+Arrays.toString(counts));
		}
		void checkPlan(final PlanRow row){
			require(row.records==spec.records && row.bases==spec.bases && row.bytes==data.length &&
				row.sha80.equals(sha80), "Saved timing plan shape oracle differs for "+spec.name);
		}
		void checkBytes(){
			require(data.length==expectedBytes(spec), spec.name+" byte oracle differed: "+data.length);
		}
		final ShapeSpec spec;
		final byte[] data;
		final String sha80;
	}

	/** One saved baseline pilot plan row. */
	private static final class PlanRow{
		PlanRow(final String[] fields){
			require(fields.length==9, "Bad timing plan row: "+Arrays.toString(fields));
			shape=fields[0]; reps=Integer.parseInt(fields[1]); warmup=Integer.parseInt(fields[2]);
			records=Integer.parseInt(fields[3]); bases=Long.parseLong(fields[4]);
			bytes=Integer.parseInt(fields[5]); sha80=fields[6];
			require(reps>=MIN_MEASURE_REPS && reps<=MAX_MEASURE_REPS,
				"Saved reps for "+shape+" must be in "+MIN_MEASURE_REPS+".."+MAX_MEASURE_REPS+": "+reps);
			require(warmup==warmup(reps), "Saved warmup does not match saved reps for "+shape);
		}
		final String shape, sha80;
		final int reps, warmup, records, bytes;
		final long bases;
	}

	/** One timed sample row. */
	private final class Result{
		Result(final String shape_, final int reps_, final int warmup_, final int round_,
				final long cpuNanos_, final long wallNanos_, final long checksum_, final boolean accepted_){
			shape=shape_; reps=reps_; warmup=warmup_; round=round_;
			cpuNanos=cpuNanos_; wallNanos=wallNanos_; checksum=checksum_; accepted=accepted_;
		}
		@Override
		public String toString(){
			return shape+"\t"+phase+"\t"+mode+"\t"+round+"\t"+reps+"\t"+warmup+
				"\t"+cpuNanos+"\t"+wallNanos+"\t"+checksum+"\t"+accepted+"\n";
		}
		final String shape;
		final int reps, warmup, round;
		final long cpuNanos, wallNanos, checksum;
		final boolean accepted;
	}

	/** One scanner call result. */
	private static final class RunResult{
		RunResult(final long cpuNanos_, final long wallNanos_, final long checksum_){
			cpuNanos=cpuNanos_; wallNanos=wallNanos_; checksum=checksum_;
		}
		final long cpuNanos, wallNanos, checksum;
	}

	/** Generates exactly one of the four fixed ASCII shapes. */
	private static byte[] make(final ShapeSpec spec) throws Exception{
		final java.io.ByteArrayOutputStream out=new java.io.ByteArrayOutputStream(expectedBytes(spec));
		for(int record=0; record<spec.records; record++){
			out.write(spec.fastq ? '@' : '>');
			out.write('r');
			out.write('\n');
			writeBases(out, spec.a, 'A');
			writeBases(out, spec.n, 'N');
			writeBases(out, spec.c, 'C');
			out.write('\n');
			if(spec.fastq){
				out.write('+');
				out.write('\n');
				writeBases(out, spec.length, 'I');
				out.write('\n');
			}
		}
		return out.toByteArray();
	}

	/** Returns the exact byte count for a shape. */
	private static int expectedBytes(final ShapeSpec spec){
		final int perRecord=(spec.fastq ? 3+spec.length+1+2+spec.length+1 : 3+spec.length+1);
		return perRecord*spec.records;
	}

	/** Writes one repeated ASCII byte without Java 11 APIs. */
	private static void writeBases(final java.io.ByteArrayOutputStream out,
			final int count, final int base){
		for(int i=0; i<count; i++){out.write(base);}
	}

	/** Returns the final eighty bits of a byte-array SHA-256 as lowercase hex. */
	private static String sha80(final byte[] data) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(data);
		final StringBuilder sb=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){sb.append(String.format("%02x", digest[i]&0xFF));}
		return sb.toString();
	}

	private final Path root;
	private final String phase, mode;
	private final String selectedShapeName;
	private final ThreadMXBean clock;
	private final LinkedHashMap<String, Shape> shapes=new LinkedHashMap<String, Shape>();
	private final ArrayList<Shape> selectedShapes=new ArrayList<Shape>();
	private final int measureRounds;
	private int samples, accepted, failures;
	private static volatile long SINK;

	private static final String[] SCANNER_ARGS={"in=stdin", "cutoff=2048", "n=2", "format=1",
		"showbbmap=f", "showspeed=f", "n_=f"};
	private static final long BASES=4194304L;
	private static final long TARGET_CPU_NANOS=12000000000L;
	private static final long MIN_ACCEPTED_CPU_NANOS=10000000000L;
	private static final long MIN_PILOT_CPU_NANOS=250000000L;
	private static final long MIN_PILOT_WALL_NANOS=250000000L;
	private static final int MIN_MEASURE_REPS=16;
	private static final int MAX_MEASURE_REPS=8192;
	private static final int MAX_PILOT_REPS=256;
	private static final int MIN_WARMUP_REPS=3;
	private static final int MAX_MEASURE_WARMUP_REPS=64;
	private static final int MAX_PILOT_TOTAL_REPS_PER_SHAPE=511;
	private static final int MAX_PILOT_WARMUP_REPS_PER_SHAPE=46;
	private static final int MEASURE_ROUNDS=3;
	private static final ShapeSpec[] SHAPES={
		new ShapeSpec("fasta-long-a", false, 1, (int)BASES, (int)BASES, 0, 0,
			new long[]{BASES, 0, 0, 0, 0, 0, 0, 0}),
		new ShapeSpec("fasta-split-n", false, 4096, 1024, 512, 510, 2,
			new long[]{2097152, 2088960, 0, 0, 0, 8192, 0, 0}),
		new ShapeSpec("fastq-long-a", true, 1024, 4096, 4096, 0, 0,
			new long[]{BASES, 0, 0, 0, 0, 0, 0, 0}),
		new ShapeSpec("fastq-split-n", true, 4096, 1024, 512, 510, 2,
			new long[]{2097152, 2088960, 0, 0, 0, 8192, 0, 0}),
	};
	private static final StaticRef[] STATIC_FIELDS={
		new StaticRef(AssemblyStats2.class, "GCFORMAT"),
		new StaticRef(AssemblyStats2.class, "GC_FRACTION"),
		new StaticRef(AssemblyStats2.class, "FORMAT"),
		new StaticRef(AssemblyStats2.class, "cutoff"),
		new StaticRef(AssemblyStats2.class, "LIMSUM"),
		new StaticRef(AssemblyStats2.class, "HEADERLENSUM"),
		new StaticRef(AssemblyStats2.class, "bbmapkmer"),
		new StaticRef(AssemblyStats2.class, "overwrite"),
		new StaticRef(AssemblyStats2.class, "append"),
		new StaticRef(AssemblyStats2.class, "useheader"),
		new StaticRef(AssemblyStats2.class, "addfilename"),
		new StaticRef(AssemblyStats2.class, "printL90"),
		new StaticRef(AssemblyStats2.class, "swapNL"),
		new StaticRef(AssemblyStats2.class, "printExtended"),
		new StaticRef(AssemblyStats2.class, "logSumOffset"),
		new StaticRef(AssemblyStats2.class, "powSumPower"),
		new StaticRef(AssemblyStats2.class, "logSumBase"),
		new StaticRef(AssemblyStats2.class, "squareLog"),
		new StaticRef(AssemblyStats2.class, "logPower"),
		new StaticRef(AssemblyStats2.class, "showspeed"),
		new StaticRef(AssemblyStats2.class, "printheadersize"),
		new StaticRef(AssemblyStats2.class, "skipDuplicateLines"),
		new StaticRef(AssemblyStats2.class, "N_UNDERSCORE"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMinLen"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMaxLen"),
		new StaticRef(AssemblyStats2.class, "alignedFraction"),
		new StaticRef(AssemblyStats2.class, "printAssemblyScore"),
		new StaticRef(AssemblyStats2.class, "lastLogSum"),
		new StaticRef(AssemblyStats2.class, "lastL90"),
		new StaticRef(AssemblyStats2.class, "lastL50"),
		new StaticRef(AssemblyStats2.class, "lastSize"),
		new StaticRef(AssemblyStats2.class, "lastContigs"),
		new StaticRef(AssemblyStats2.class, "lastMaxContig"),
		new StaticRef(AssemblyStats2.class, "outstream"),
		new StaticRef(ReadWrite.class, "USE_PIGZ"),
		new StaticRef(ReadWrite.class, "USE_UNPIGZ"),
		new StaticRef(Shared.class, "mainClass"),
		new StaticRef(Shared.class, "COMMAND_LINE")
	};
	private static final String SHAPES_HEADER=
		"shape\tformat\trecords\tsequenceBases\tbytes\tsha80\texpectedCounts\n";
	private static final String SAMPLE_HEADER=
		"shape\tphase\tmode\tround\treps\twarmup\tcpuNanos\twallNanos\tchecksum\taccepted\n";
	private static final String PLAN_HEADER=
		"shape\treps\twarmup\trecords\tsequenceBases\tbytes\tsha80\tpilotCpuNanos\tpilotWallNanos\n";
}
