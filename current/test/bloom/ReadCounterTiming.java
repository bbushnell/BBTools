package test.bloom;

import java.util.Locale;

import bloom.KCountArray;
import bloom.KCountArray7MTA;
import bloom.ReadCounter;
import parse.Parse;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import ukmer.Kmer;

/** One real unpaired input pass through ReadCounter for matched process-level timing.
 * Run in a fresh JVM and time the whole process externally. Pins quality, geometry,
 * worker and hash settings, prints exact observed counters, then the existing report.
 * This wrapper neither predicts eligibility nor hides changed work across versions.
 * @author Brian Bushnell, Sucrose */
public final class ReadCounterTiming{

	/** Requires all five workload settings; does not accept a mate input or perform a second pass. */
	public static void main(final String[] args) throws Exception{
		String input=null;
		long reads=0, cells=0;
		int k=0, bits=0, seen=0;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			check(equals>0 && equals<arg.length()-1, "Expected nonempty flag=value: "+arg);
			final String name=arg.substring(0, equals).toLowerCase(Locale.ROOT), value=arg.substring(equals+1);
			final int flag;
			if(name.equals("in")){
				input=value; flag=1;
			}else if(name.equals("reads")){
				reads=Parse.parseKMG(value); flag=2;
			}else if(name.equals("k")){
				k=Integer.parseInt(value); flag=4;
			}else if(name.equals("cells")){
				cells=Parse.parseKMG(value); flag=8;
			}else if(name.equals("bits")){
				bits=Integer.parseInt(value); flag=16;
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
			check((seen&flag)==0, "Duplicate workload parameter: "+name);
			seen|=flag;
		}
		check(seen==31, "Specify in=<FASTQ> reads=<positive limit> k=<long span> cells=<positive target> bits=<counter width>.");
		check(input!=null && !input.equalsIgnoreCase("null") && reads>0 && k>31 && cells>0,
				"This workload requires a real input, positive read/cell limits and long-k counting.");
		check(bits>0 && bits<=32 && Integer.bitCount(bits)==1, "Packed counter bits must be 1,2,4,8,16 or32.");
		configure();
		check(Kmer.getKbig(k)==k, "Packed geometry must retain the requested long-window span.");
		ReadCounter.maxReads=reads;
		//This is the only table constructed in this fresh JVM. Its mask counter starts at0.
		KCountArray7MTA.setSeed(0);
		final KCountArray table=KCountArray.makeNew(cells, bits, 1);
		System.out.println("CONFIG k="+k+" represented_k="+Kmer.getKbig(k)+" words="+Kmer.getMult(k)+
				" reads_limit="+reads+" requested_cells="+cells+" actual_cells="+table.cells+
				" bits="+table.cellBits+" hashes=1 seed_mask=0 table_ordinal=0 threads="+Shared.threads()+
				" qin="+FASTQ.ASCII_OFFSET+" min_called_quality="+Read.MIN_CALLED_QUALITY()+
				" max_called_quality="+Read.MAX_CALLED_QUALITY()+" min_prob="+ReadCounter.minProb+
				" packed="+Kmer.PACKED+" full_mix="+Kmer.FULL_MIX+" mask_core="+Kmer.MASK_CORE+
				" buffered="+ReadCounter.BUFFERED+" sketch="+ReadCounter.SKETCH_MODE+
				" keep_duplicates="+ReadCounter.KEEP_DUPLICATE_KMERS+" sampled_windows="+ReadCounter.KMERS_PER_READ+
				" idmask="+ReadCounter.IDMASK+" interleaved=false quality_detection=false"+
				" validation=true quality_capping=true overlap=false trusted=false dual_hash=false");
		final ReadCounter counter=new ReadCounter(k, true, false, false, false, false);
		try{
			check(counter.count(input, null, table)==table, "ReadCounter must return its borrowed destination.");
		}finally{
			//Table finalization is ours; the real count method owns and closes its input.
			table.shutdown();
			System.out.println("OBSERVED reads="+ReadCounter.readsProcessed+" keys="+ReadCounter.keysCounted+
					" increments="+ReadCounter.increments);
		}
		check(ReadCounter.readsProcessed>0 && ReadCounter.readsProcessed<=reads,
				"Unpaired observed reads must lie within the positive requested input limit.");
		check(ReadCounter.keysCounted>=0 && ReadCounter.increments==ReadCounter.keysCounted,
				"Unbuffered nonsketch single-key counting issues one increment per accepted occurrence.");
		counter.printStatistics(table);
	}

	/** Pins globals before table/reader creation; fresh-process use prevents cross-call state reuse. */
	private static void configure(){
		Locale.setDefault(Locale.ROOT);
		Shared.setThreads(1);
		FASTQ.ASCII_OFFSET=33;
		FASTQ.DETECT_QUALITY=FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		Read.VALIDATE_IN_CONSTRUCTOR=Read.CHANGE_QUALITY=true;
		Read.SKIP_SLOW_VALIDATION=false;
		Read.setMinCalledQuality(2);
		Read.setMaxCalledQuality(50);
		Kmer.PACKED=Kmer.FULL_MIX=true;
		Kmer.MASK_CORE=false;
		ReadCounter.maxShortKmerLength=31;
		ReadCounter.MAX_COUNT_THREADS=1;
		ReadCounter.minQuality=0;
		ReadCounter.minProb=0.5f;
		ReadCounter.KEEP_DUPLICATE_KMERS=ReadCounter.CANONICAL=true;
		ReadCounter.KMERS_PER_READ=ReadCounter.IDMASK=0;
		ReadCounter.BUFFERED=ReadCounter.SKETCH_MODE=ReadCounter.STORE_HASHED=ReadCounter.PREJOIN=false;
		check(ReadCounter.keysCounted==0 && ReadCounter.increments==0 && ReadCounter.readsProcessed==0,
				"Timing requires a fresh JVM with zero observed counters, not counters reset after unknown prior work.");
	}

	/** Rejects invalid setup regardless of JVM assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
}
