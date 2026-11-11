package test.assemble;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import fileIO.FileFormat;
import kmer.AbstractKmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;

/** Matched reset-plus-reassembly kernel timing on six existing substitution/clean reads.
 * Loads real count tables and calls the unmodified shared pass and concrete inner methods.
 * Input loading, full validation and reporting are outside measured batches. This is not
 * whole-CLI throughput or biological accuracy evidence. Each invocation needs a fresh JVM.
 * @author Brian Bushnell, Sucrose
 */
public final class TadpoleReassemblyBench{

	/** Parses named fixture inputs and selects a sizing pilot or a fixed cycle count. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TadpoleReassemblyBench.class, false);
		try{
			String queries=null, reference=null, truth=null, phase="pilot";
			int k=0;
			long cycles=0;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2, "Expected flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("in")){queries=pair[1];}
				else if(pair[0].equalsIgnoreCase("ref")){reference=pair[1];}
				else if(pair[0].equalsIgnoreCase("truth")){truth=pair[1];}
				else if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("phase")){phase=pair[1];}
				else if(pair[0].equalsIgnoreCase("cycles")){cycles=Long.parseLong(pair[1]);}
				else{throw new IllegalArgumentException("Unknown benchmark option: "+arg);}
			}
			check(queries!=null && reference!=null && truth!=null && (k==31 || k==63),
					"Required in=queries.fq ref=counts.fa truth=expected.fq k=31|63.");
			final boolean pilot="pilot".equals(phase);
			check(pilot || "measure".equals(phase), "Use phase=pilot|measure.");
			check(pilot ? cycles==0 : cycles>0 && cycles<=1000000000L && cycles%BATCH==0,
					"Pilot requires cycles=0; measure requires a positive multiple of64, at most1e9.");
			final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
			check(clock.isCurrentThreadCpuTimeSupported(), "A real main-thread CPU clock is required.");
			if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
			check(clock.getCurrentThreadCpuTime()>=0, "Main-thread CPU clock is unavailable.");
			final long setupCpu=clock.getCurrentThreadCpuTime(), setupWall=System.nanoTime();
			final TadpoleReassemblyBench bench=new TadpoleReassemblyBench(queries, reference, truth, k);
			try{
				System.out.println("CONFIG k="+k+" phase="+phase+" cases=6 threads="+Shared.threads()+
						" packed="+Kmer.PACKED+" maskcore="+Kmer.MASK_CORE+" simd="+Shared.SIMD+
						" simd128="+simd.Vector.simd128+" simd256="+simd.Vector.simd256+
						" assertions="+TadpoleReassemblyBench.class.desiredAssertionStatus()+
						" qualities=uniform30 batch_cycles=64 requested_cycles="+cycles+
						" expected_corrections_per_cycle=6 unit=reset_plus_reassembly");
				bench.preflight(false);
				bench.preflight(true);
				System.out.println("SETUP wall_ns="+(System.nanoTime()-setupWall)+
						" main_cpu_ns="+(clock.getCurrentThreadCpuTime()-setupCpu));
				bench.measure(clock, pilot ? SECOND : 3*SECOND, 0, "WARMUP");
				bench.validate(true);
				final long[] measured=bench.measure(clock, pilot ? SECOND : 0, cycles, "SAMPLE");
				bench.validate(true);
				check(pilot || measured[2]>=10*SECOND,
						"Measured work was below ten main-thread CPU seconds; no performance conclusion.");
				System.out.println("REASSEMBLY_BENCH_PASS k="+k+" phase="+phase+
						" cases=6 corrections_per_cycle=6 output_checks=true");
			}finally{bench.tadpole.clearData();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Loads the same immutable support table once and materializes the six literal cases. */
	private TadpoleReassemblyBench(final String queries, final String reference, final String truth, final int k){
		final String[] options={"in="+reference, "out=null", "mode=correct", "k="+k, "packed=f",
				"t=1", "buildthreads=1", "ways=1", "initialsize=32768", "prealloc=f", "prefilter=0",
				"minprob=0", "minprobmain=f", "rcomp=t", "maskmiddle=f", "maskcore=f", "qin=33",
				"shave=f", "rinse=f", "interleaved=f", "pigz=f", "unpigz=f", "showstats=f",
				"aggressive=f", "conservative=f", "errorpath=1", "errormult1=16", "errormult2=2.6",
				"errormultqfactor=0.002", "errorlowerconst=4", "mincountcorrect=3", "deadzone=0",
				"window=12", "windowcount=6", "qualsum=80", "ecc=f", "ecco=f", "merge=f",
				"fixindels=f", "hpindel=f", "markerrors=f", "hollow=f", "hollowrescue=f"};
		tadpole=Tadpole.makeTadpole(options, false);
		try{
			check(Shared.threads()==1 && !Kmer.PACKED && !Kmer.MASK_CORE && !AbstractKmerTableSet.MASK_CORE,
					"Benchmark requires one thread and the declared nonpacked, unmasked geometry.");
			check(tadpole.tables().kbig()==k && (k==31 ? tadpole instanceof Tadpole1 : tadpole instanceof Tadpole2),
					"Factory did not select the required concrete representation and exact K.");
			tadpole.loadKmers(new Timer());
			final ArrayList<Read> inputs=ReadInputStream.toReads(queries, FileFormat.FASTQ, -1);
			final ArrayList<Read> truths=ReadInputStream.toReads(truth, FileFormat.FASTQ, -1);
			check(inputs.size()==18 && truths.size()==18, "Pinned panel requires eighteen query/truth records.");
			for(int i=0; i<NAMES.length; i++){
				cases[i]=new Case(find(inputs, NAMES[i]), find(truths, NAMES[i]), k, POSITIONS[i], tadpole);
			}
		}catch(RuntimeException e){tadpole.clearData(); throw e;}
		catch(Error e){tadpole.clearData(); throw e;}
	}

	/** Requires exactly one unpaired record with the selected literal ID. */
	private static Read find(final ArrayList<Read> reads, final String name){
		Read found=null;
		for(Read r : reads){if(name.equals(r.id)){
			check(found==null && r.mate==null, "Duplicate or paired fixture record: "+name);
			found=r;
		}}
		check(found!=null, "Missing fixture record: "+name);
		return found;
	}

	/** Executes and fully validates an untimed null or uniform-quality cycle. */
	private void preflight(final boolean quality){
		check(cycle(quality)==expectedChecksum(quality), "Initial cycle checksum violates the literal six-read workload.");
		validate(quality);
		System.out.println("PREFLIGHT qualities="+(quality ? "uniform30" : "null")+
				" cases=6 corrections=6 bases=true qualities=true counts=true tracker=true");
	}

	/** Measures batches with coarse clock reads; resets are intentionally inside both clocks. */
	private long[] measure(final ThreadMXBean clock, final long targetCpu, final long limit, final String phase){
		assert(targetCpu>0 || limit>0) : "Timing must terminate at a CPU target or fixed cycle count.";
		final long expected=expectedChecksum(true);
		long cycles=0, checksum=0;
		final long cpu0=clock.getCurrentThreadCpuTime(), wall0=System.nanoTime();
		long cpu;
		do{
			long batchSum=0;
			for(int i=0; i<BATCH; i++){batchSum+=cycle(true);}
			check(batchSum==BATCH*expected, "Timed batch checksum changed; reject the sample.");
			checksum+=batchSum;
			cycles+=BATCH;
			cpu=clock.getCurrentThreadCpuTime()-cpu0;
		}while(limit>0 ? cycles<limit : cpu<targetCpu);
		final long wall=System.nanoTime()-wall0;
		check(cpu>0 && wall>0 && checksum==Math.multiplyExact(cycles, expected), "Invalid timing clock or aggregate checksum.");
		System.out.println(phase+" cycles="+cycles+" reads="+(6*cycles)+" corrections="+(6*cycles)+
				" checksum="+checksum+" expected_per_cycle="+expected+" wall_ns="+wall+" main_cpu_ns="+cpu);
		return new long[]{cycles, wall, cpu};
	}

	/** Runs one real pass per case, including reset and constant-sized work checks. */
	private long cycle(final boolean quality){
		long checksum=0;
		for(Case c : cases){
			c.reset(quality);
			final int corrected=tadpole.reassemble_pass(c.bases, quality ? c.quals : null, c.left, c.right,
					c.rightCounts, c.counts, c.counts2, c.tracker, 5, c.kmer, c.regenKmer, c.changed);
			if(corrected!=c.positions.length || c.tracker.corrected()!=corrected || c.tracker.detected()!=corrected){
				throw new IllegalStateException("Correction/detection counts differ from the literal mutation workload: "+c.name);
			}
			checksum+=corrected;
			for(int position : c.positions){checksum+=c.bases[position]+(quality ? c.quals[position] : 0);}
		}
		return checksum;
	}

	/** Independent checksum from literal mutation positions and expected truth bases/qualities. */
	private long expectedChecksum(final boolean quality){
		long checksum=6;
		for(Case c : cases){for(int p : c.positions){checksum+=c.truth[p]+(quality ? 32 : 0);}}
		return checksum;
	}

	/** Checks every output base, quality, count, changed bit and tracker field outside timing. */
	private void validate(final boolean quality){
		for(Case c : cases){c.validate(quality);}
	}

	/** Always-on correctness guard; an unexpected legacy result stops sizing and timing. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Owns one immutable input/truth pair and all reused mutable pass state. */
	private static final class Case{
		Case(final Read input, final Read expected, final int k, final int[] positions_, final Tadpole tad){
			name=input.id;
			positions=positions_.clone();
			original=input.bases.clone(); truth=expected.bases.clone();
			final int length=name.startsWith("case0_") ? 241 : 244;
			check(original.length==length && truth.length==length, "Unexpected selected record length: "+name);
			int differences=0;
			for(int i=0; i<length; i++){
				check(isBase(original[i]) && isBase(truth[i]), "Undefined fixture base: "+name);
				if(original[i]!=truth[i]){
					check(differences<positions.length && positions[differences]==i, "Unexpected mutation coordinate: "+name);
					differences++;
				}
			}
			check(differences==positions.length, "Missing literal mutation: "+name);
			bases=new byte[length]; quals=new byte[length];
			left=new ByteBuilder(length); right=new ByteBuilder(length);
			counts=new IntList(length); counts2=new IntList(length);
			kmer=new Kmer(k); regenKmer=new Kmer(k);
			check(kmer.kbig==k, "Scratch Kmer geometry differs from table K.");
			tad.tables().fillCounts(original, counts, kmer);
			check(counts.size==length-k+1, "Input profile lacks a count for every valid window.");
			originalCounts=Arrays.copyOf(counts.array, counts.size);
			tad.tables().fillCounts(truth, counts, kmer);
			check(counts.size==length-k+1, "Truth profile lacks a count for every valid window.");
			truthCounts=Arrays.copyOf(counts.array, counts.size);
			for(int count : truthCounts){check(count>0, "Truth lacks real counted support: "+name);}
		}

		/** Restores the original state without new per-read arrays or table reconstruction. */
		void reset(final boolean quality){
			System.arraycopy(original, 0, bases, 0, bases.length);
			if(quality){Arrays.fill(quals, (byte)30);}
			System.arraycopy(originalCounts, 0, counts.array, 0, originalCounts.length);
			counts.size=originalCounts.length;
			tracker.clear(); changed.clear();
		}

		/** Literal bases/positions are independent of the correction implementation. */
		void validate(final boolean quality){
			check(Arrays.equals(bases, truth), "Corrected bases differ from independent truth: "+name);
			check(counts.size==truthCounts.length && changed.cardinality()==positions.length,
					"Incorrect count-list size or changed-position count: "+name);
			for(int i=0; i<counts.size; i++){check(counts.get(i)==truthCounts[i], "Regenerated count differs from loaded truth: "+name);}
			int p=0;
			for(int i=0; i<bases.length; i++){
				final boolean edited=p<positions.length && positions[p]==i;
				check(changed.get(i)==edited, "Unexpected changed-position bit: "+name);
				if(quality){check(quals[i]==(edited ? 32 : 30), "Legacy quality differs from bidirectional oracle: "+name);}
				if(edited){p++;}
			}
			check(tracker.correctedReassembleInner==positions.length && tracker.correctedReassembleOuter==0 &&
					tracker.corrected()==positions.length && tracker.detected()==positions.length &&
					tracker.detectedReassemble==positions.length && tracker.detectedPincer==0 &&
					tracker.detectedTail==0 && tracker.detectedBrute==0 && tracker.detectedSubstitute==0 &&
					tracker.correctedPincer==0 && tracker.correctedTail==0 && tracker.correctedBrute==0 && tracker.correctedSubstitute==0 &&
					!tracker.rollback && !tracker.rollbackSubstitute && tracker.rollbackTrigger==0 &&
					tracker.rollbackCorrectedSnapshot==0 && tracker.marked==0 && tracker.suspected==0,
					"Unexpected correction tracker state: "+name);
		}

		private static boolean isBase(final byte b){return b=='A' || b=='C' || b=='G' || b=='T';}
		final String name;
		final int[] positions, originalCounts, truthCounts;
		final byte[] original, truth, bases, quals;
		final ByteBuilder left, right;
		final IntList counts, counts2;
		final Kmer kmer, regenKmer;
		final int[] rightCounts=new int[4];
		final ErrorTracker tracker=new ErrorTracker();
		final BitSet changed=new BitSet();
	}

	final Tadpole tadpole;
	final Case[] cases=new Case[6];
	private static final long SECOND=1000000000L;
	private static final int BATCH=64;
	private static final String[] NAMES={"case0_F", "case0_RC", "case5_F", "case5_RC", "case8_F", "case8_RC"};
	private static final int[][] POSITIONS={{120}, {120}, {}, {}, {70, 170}, {73, 173}};
}
