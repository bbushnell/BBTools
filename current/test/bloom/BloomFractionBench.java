package test.bloom;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Arrays;
import java.util.Locale;

import bloom.BloomFilter;
import bloom.KCountArray7MTA;
import shared.Shared;
import stream.Read;
import ukmer.Kmer;

/** Matched fraction and average-depth queries over a small cache-warm read pool.
 * Tests final accounting overhead, not whole-program filtering throughput.
 * Run each k/profile in a fresh JVM; construction and input checks are untimed.
 * @author Brian Bushnell, Shinobu */
public final class BloomFractionBench{

	/** Warms up for two main-thread CPU seconds, then measures at least ten. */
	public static void main(final String[] args){
		int k=31;
		boolean ambiguous=false;
		boolean average=false;
		for(String arg : args){
			if(arg.startsWith("k=")){k=Integer.parseInt(arg.substring(2));}
			else if(arg.equals("profile=clean")){ambiguous=false;}
			else if(arg.equals("profile=ambiguous")){ambiguous=true;}
			else if(arg.equals("mode=average")){average=true;}
			else if(arg.equals("mode=fraction")){average=false;}
			else{throw new IllegalArgumentException("Expected k=31|63 profile=clean|ambiguous mode=fraction|average: "+arg);}
		}
		check(k==31 || k==63, "Only the two explicitly specified kernel widths are supported.");
		check(!average || k==31, "averageCount supports primitive k-mers only.");
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		check(clock.isCurrentThreadCpuTimeSupported(), "An actual main-thread CPU clock is required.");
		final boolean oldClock=clock.isThreadCpuTimeEnabled();
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		final boolean oldMask=Kmer.MASK_CORE, oldPacked=Kmer.PACKED, oldMix=Kmer.FULL_MIX;
		final boolean oldPrint=BloomFilter.printMem;
		try{
			check(Shared.threads()==1, "Run with -XX:ActiveProcessorCount=1 for the matched one-worker configuration.");
			BloomFilter.OVERRIDE_CELLS=100003;
			BloomFilter.printMem=false;
			Kmer.MASK_CORE=false;
			Kmer.PACKED=Kmer.FULL_MIX=true;
			KCountArray7MTA.setSeed(917);
			final BloomFractionBench bench=new BloomFractionBench(k, ambiguous, average);
			System.out.println("CONFIG mode="+(average ? "average" : "fraction")+" k="+k+" profile="+(ambiguous ? "ambiguous" : "clean")+
					" pool=32 lengths=128..159 canonical=true hashes=1 bits=8 threshold=2 fraction=0.8"+
					" loadedDepth=3 cells="+bench.filter.filter.cells+" validWindowsPerCycle="+bench.validWindows+
					" inputFingerprint="+bench.fingerprint+" seed=917 maskCore=false batch=4096"+
					" assertions="+BloomFractionBench.class.desiredAssertionStatus());
			bench.measure(clock, 2, "warmup");
			bench.measure(clock, 10, "sample");
		}finally{
			BloomFilter.OVERRIDE_CELLS=oldCells;
			BloomFilter.printMem=oldPrint;
			Kmer.MASK_CORE=oldMask;
			Kmer.PACKED=oldPacked;
			Kmer.FULL_MIX=oldMix;
			if(!oldClock){clock.setThreadCpuTimeEnabled(false);}
		}
	}

	/** Loads one repeated C/G key; all defined windows have exactly count three. */
	private BloomFractionBench(final int k, final boolean ambiguous, final boolean average_){
		average=average_;
		filter=new BloomFilter(k, k, 8, 1, 2, true, false, false, 0.01f);
		final Kmer key=new Kmer(k);
		for(int j=0; j<k; j++){key.addRight((byte)'C');}
		long raw=0;
		if(k<=31){for(int j=0; j<k; j++){raw=(raw<<2)|2;}}
		else{raw=key.xor();}
		filter.filter.increment(raw, 3);
		check(filter.getCount(raw)==3, "Loaded key must retain its literal depth.");
		int windows=0, hash=1;
		for(int i=0; i<reads.length; i++){
			final byte[] bases=new byte[128+i];
			Arrays.fill(bases, (byte)'C');
			if(ambiguous){bases[64]='N';}
			reads[i]=new Read(bases, null, "fraction-"+i, i, false);
			originals[i]=bases.clone();
			//Independent geometry: scan every literal window; no production fraction helper.
			int valid=0;
			for(int start=0; start+k<=bases.length; start++){
				boolean defined=true;
				for(int pos=start; pos<start+k; pos++){defined&=bases[pos]=='C';}
				if(defined){valid++;}
			}
			check(valid>0, "Each benchmark read must contain at least one complete C-only window.");
			windows+=valid;
			hash=31*hash+Arrays.hashCode(bases);
		}
		validWindows=windows;
		fingerprint=hash;
		validate();
	}

	/** Checks every result and input outside clocks, including long-key encoding. */
	private void validate(){
		for(int i=0; i<reads.length; i++){
			check(Arrays.equals(reads[i].bases, originals[i]), "Fraction query modified input "+i);
			check(reads[i].quality==null && reads[i].mate==null, "Query changed read ownership.");
			check(filter.hasHighCountFraction(reads[i], 2, 0.8f),
					"Every defined window is a loaded C/G key at depth three; query must pass.");
			if(average){check(filter.averageCount(reads[i].bases)==3f, "Uniform depths must average to exactly three.");}
		}
	}

	/** Consumes one result per read; no harness allocations or clocks within a batch. */
	private long batch(){
		long count=0;
		if(average){
			for(int i=0; i<BATCH; i++){count+=(long)filter.averageCount(reads[i&31].bases);}
		}else{
			for(int i=0; i<BATCH; i++){count+=filter.hasHighCountFraction(reads[i&31], 2, 0.8f) ? 1 : 0;}
		}
		return count;
	}

	/** Reports normalized cost because the fixed CPU duration permits different work counts. */
	private void measure(final ThreadMXBean clock, final int seconds, final String phase){
		validate();
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		check(cpuStart>=0, "CPU clock must return a usable timestamp.");
		long cpuEnd, work=0, checksum=0;
		do{
			final long result=batch();
			check(result==(average ? 3L*BATCH : BATCH), "Timed results must retain their literal fraction/depth checksum.");
			work+=BATCH;
			checksum+=result;
			sink=checksum;
			cpuEnd=clock.getCurrentThreadCpuTime();
		}while(cpuEnd-cpuStart<seconds*1000000000L);
		final long wall=System.nanoTime()-wallStart, cpu=cpuEnd-cpuStart;
		validate();
		System.out.printf(Locale.ROOT, "RESULT phase=%s work=%d checksum=%d wall_seconds=%.9f thread_cpu_seconds=%.9f wall_ns_per_read=%.3f cpu_ns_per_read=%.3f%n",
				phase, work, checksum, wall*1e-9, cpu*1e-9, wall/(double)work, cpu/(double)work);
	}

	/** Keeps fixture failures active under either assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private final BloomFilter filter;
	private final boolean average;
	private final Read[] reads=new Read[32];
	private final byte[][] originals=new byte[32][];
	private final int validWindows, fingerprint;
	private static final int BATCH=4096;
	private static volatile long sink;
}
