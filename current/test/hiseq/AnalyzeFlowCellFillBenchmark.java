package test.hiseq;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.ConcurrentHashMap;

import bloom.BloomFilter;
import hiseq.AnalyzeFlowCell;
import hiseq.FlowCell;
import hiseq.MicroTile;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import structures.AtomicStringNum;
import structures.ListNum;

/** Bounded-memory CPU workload for the real single-worker flow-cell filling stage.
 * Excludes file decoding, Bloom construction from reads, alignment and final filtering.
 * Reuses eight valid reads; this is a stage benchmark, not an end-to-end throughput claim.
 * @author Brian Bushnell, Shinobu */
public class AnalyzeFlowCellFillBenchmark{

	/** Reads and warmup are decimal multiples of eight; short default runs are calibration only. */
	public static void main(String[] args) throws Exception{
		long reads=100000, warmup=10000;
		boolean fixed=false;
		for(String arg : args){
			if(arg.startsWith("reads=")){reads=Long.parseLong(arg.substring(6));}
			else if(arg.startsWith("warmup=")){warmup=Long.parseLong(arg.substring(7));}
			else if(arg.equals("phase=fixed")){fixed=true;}
			else if(!arg.equals("phase=baseline")){throw new IllegalArgumentException(arg);}
		}
		require(reads>0 && reads%8==0 && reads<=Long.MAX_VALUE/150,
			"Reads must be a positive multiple of8 with representable base totals");
		require(warmup>=0 && warmup%8==0 && warmup<=Long.MAX_VALUE/150,
			"Warmup must be a nonnegative multiple of8 with representable base totals");
		if(warmup>0){run(warmup, fixed, false);}
		System.out.println("reads\tbases\tsampled_misses\tbarcode_count\tfill_wall_seconds");
		run(reads, fixed, true);
	}

	/** Runs real worker/model code and checks exact counters outside the timed region. */
	private static void run(long reads, boolean fixed, boolean report) throws Exception{
		final AnalyzeFlowCell tool=new AnalyzeFlowCell(new String[] {"in=/dev/null", "extin=fq",
			"threads=1", "fillthreads=1", "loadkmers=t", "align=f", "simd=f", "qin=33", "idmaskread=7"});
		final FlowCell flowcell=new FlowCell(31);
		set(tool, "flowcell", flowcell);
		final ConcurrentHashMap<String,AtomicStringNum> counts=new ConcurrentHashMap<String,AtomicStringNum>();
		set(tool, "barcodeMap", counts);
		final Method loadBarcodes=AnalyzeFlowCell.class.getDeclaredMethod("loadBarcodes", String.class);
		loadBarcodes.setAccessible(true);
		set(tool, "barcodeStats", loadBarcodes.invoke(tool, (Object)null));
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		final BloomFilter filter;
		try{
			BloomFilter.OVERRIDE_CELLS=4096;
			filter=new BloomFilter(31, 31, 2, 2, 1, true, false, false, 0.001f);
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
		set(tool, "bloomFilter", filter);
		final RepeatedReads input=new RepeatedReads(reads);
		long elapsed;
		try{
			input.start();
			final long start=System.nanoTime();
			tool.fillTilesInner(input, null);
			elapsed=System.nanoTime()-start;
		}finally{input.close(); filter.filter.shutdown();}
		require(!tool.errorState && tool.readsProcessed==reads && tool.basesProcessed==reads*150,
			"Worker counters must cover every generated read and base");
		long modelReads=0, modelBases=0, hits=0, misses=0, barcodes=0;
		for(MicroTile tile : flowcell.toList()){
			modelReads+=tile.readCount; modelBases+=tile.baseCount;
			hits+=tile.hits; misses+=tile.misses;
		}
		for(AtomicStringNum count : counts.values()){barcodes+=count.n.get();}
		require(modelReads==reads && modelBases==reads*150 && hits==0 && barcodes==reads,
			"Model and barcode totals must match the literal workload");
		//Across IDs0..7, each eligible end position is selected once: old31..149, fixed30..149.
		require(misses==(reads/8)*(fixed ? 120 : 119), "Unexpected sampled-window total: "+misses);
		if(report){System.out.printf(Locale.ROOT, "%d\t%d\t%d\t%d\t%.6f%n",
			reads, modelBases, misses, barcodes, elapsed*1e-9);}
	}

	/** Synthetic single-consumer source; reuses records only after their batch has been consumed. */
	private static class RepeatedReads implements Streamer{
		RepeatedReads(long count){
			remaining=count;
			final Read[] pool=new Read[8];
			for(int i=0; i<pool.length; i++){
				final byte[] bases=new byte[150], quality=new byte[150];
				Arrays.fill(quality, (byte)30);
				for(int j=0; j<bases.length; j++){bases[j]=(byte)"ACGT".charAt(j&3);}
				pool[i]=new Read(bases, quality, "INST:1:FC:1:1101:"+(i*500)+":0 1:N:0:ACGT+TGCA", i);
			}
			for(int i=0; i<1024; i++){batch.add(pool[i&7]);}
		}
		@Override
		public String fname(){return "generated-fill-workload";}
		@Override
		public void start(){}
		@Override
		public void close(){closed=true;}
		@Override
		public boolean paired(){return false;}
		@Override
		public int pairnum(){return 0;}
		@Override
		public long readsProcessed(){return delivered;}
		@Override
		public long basesProcessed(){return delivered*150;}
		@Override
		public void setSampleRate(float rate, long seed){throw new UnsupportedOperationException();}
		@Override
		public ListNum<Read> nextList(){
			if(!hasMore()){return null;}
			final int size=(int)Math.min(remaining, batch.size());
			final ArrayList<Read> list=size==batch.size() ? batch : new ArrayList<Read>(batch.subList(0, size));
			remaining-=size; delivered+=size;
			return new ListNum<Read>(list, nextID++);
		}
		@Override
		public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException();}
		@Override
		public boolean hasMore(){return !closed && remaining>0;}
		@Override
		public boolean errorState(){return false;}
		private final ArrayList<Read> batch=new ArrayList<Read>(1024);
		private long remaining, delivered, nextID;
		private boolean closed;
	}

	/** Installs setup state without exposing production-only testing APIs. */
	private static void set(Object target, String name, Object value) throws Exception{
		final Field field=target.getClass().getDeclaredField(name);
		field.setAccessible(true); field.set(target, value);
	}

	/** Keeps count/output checks active regardless of JVM assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
